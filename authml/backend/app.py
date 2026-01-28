from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from werkzeug.security import generate_password_hash, check_password_hash
from db import users_collection
from otp_helper import generate_otp, send_otp_email
from datetime import datetime, timedelta
import os
from dotenv import load_dotenv
import redis
import json 
import jwt
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from joblib import load
import uuid
import subprocess
from dock_tasks import run_docking_task
from pymongo import MongoClient

# NEW IMPORTS FOR CACHING
from cache_utils import cache_get, cache_set
from chem_db import get_chem_result, save_chem_result

load_dotenv()

app = Flask(__name__)

MONGO_URI = os.getenv("MONGO_URI")
SECRET_KEY = os.getenv("SECRET_KEY")

app.config["SECRET_KEY"] = SECRET_KEY
client = MongoClient(MONGO_URI)
CORS(app)

# ===============================
# JWT & Redis Configuration
# ===============================
JWT_SECRET_KEY = os.getenv('JWT_SECRET_KEY', 'your-secret-key-change-this')
JWT_ALGORITHM = 'HS256'
JWT_EXPIRATION_HOURS = 24
r = redis.Redis(host='localhost', port=6379, decode_responses=True)

# ===============================
# Directories
# ===============================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RECEPTOR_DIR = os.path.join(BASE_DIR, "receptor")
OUTPUT_DIR = os.path.join(BASE_DIR, "results")
DOCKED_DIR = os.path.join(BASE_DIR, "docked_outputs")
LIGAND_GEN_DIR = os.path.join(BASE_DIR, "generated_ligands")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(DOCKED_DIR, exist_ok=True)
os.makedirs(LIGAND_GEN_DIR, exist_ok=True)

# ===============================
# ML Models
# ===============================
logp_model = load(os.path.join(BASE_DIR, "logp_model.pkl"))
pic50_model = load(os.path.join(BASE_DIR, "pic50_model.pkl"))

# ===============================
# JWT Helper Functions
# ===============================
def generate_jwt_token(email, name):
    payload = {
        "email": email,
        "name": name,
        "exp": datetime.utcnow() + timedelta(hours=JWT_EXPIRATION_HOURS),
        "iat": datetime.utcnow()
    }
    return jwt.encode(payload, JWT_SECRET_KEY, algorithm=JWT_ALGORITHM)

def verify_jwt_token(token):
    try:
        return jwt.decode(token, JWT_SECRET_KEY, algorithms=[JWT_ALGORITHM])
    except:
        return None

def token_required(f):
    def decorated(*args, **kwargs):
        token = None
        if "Authorization" in request.headers:
            bearer = request.headers["Authorization"]
            if bearer.startswith("Bearer "):
                token = bearer.split(" ")[1]

        if not token:
            return jsonify({"error": "Token missing"}), 401

        payload = verify_jwt_token(token)
        if not payload:
            return jsonify({"error": "Token invalid or expired"}), 401

        request.current_user = payload
        return f(*args, **kwargs)

    decorated.__name__ = f.__name__
    return decorated

# ===============================
# AUTH ROUTES
# ===============================
@app.route("/signup", methods=["POST"])
def signup():
    data = request.json
    name = data.get("name")
    email = data.get("email")
    password = data.get("password")

    if not (name and email and password):
        return jsonify({"error": "All fields are required"}), 400

    if users_collection.find_one({"email": email}):
        return jsonify({"error": "User already exists"}), 400

    otp = generate_otp()
    temp_user = {"name": name, "email": email, "password": generate_password_hash(password)}

    r.set(f"otp:{email}", otp, ex=300)
    r.set(f"tempuser:{email}", json.dumps(temp_user), ex=300)
    send_otp_email(email, otp)

    return jsonify({"message": "OTP sent to email"}), 200


@app.route("/verify-otp", methods=["POST"])
def verify_otp():
    data = request.json
    email = data.get("email")
    otp = data.get("otp")

    stored_otp = r.get(f"otp:{email}")
    if not stored_otp or str(stored_otp) != str(otp):
        return jsonify({"error": "Invalid or expired OTP"}), 400

    temp_user = json.loads(r.get(f"tempuser:{email}"))
    users_collection.insert_one({
        "name": temp_user["name"],
        "email": temp_user["email"],
        "password": temp_user["password"],
        "isVerified": True
    })

    r.delete(f"otp:{email}")
    r.delete(f"tempuser:{email}")

    token = generate_jwt_token(temp_user["email"], temp_user["name"])
    return jsonify({"message": "Email verified", "token": token, "user": temp_user})


@app.route("/login", methods=["POST"])
def login():
    data = request.json
    email = data.get("email")
    password = data.get("password")

    # 1. CHECK REDIS CACHE
    cached_user = cache_get(f"user:{email}")
    if cached_user:
        if check_password_hash(cached_user["password"], password):
            token = generate_jwt_token(email, cached_user["name"])
            return jsonify({
                "message": "Login successful (cached)",
                "token": token,
                "user": {"name": cached_user["name"], "email": email}
            })

    # 2. ORIGINAL LOGIN FLOW
    user = users_collection.find_one({"email": email})
    if not user:
        return jsonify({"error": "User not found"}), 404

    if not user.get("isVerified"):
        return jsonify({"error": "Email not verified"}), 403

    if not check_password_hash(user["password"], password):
        return jsonify({"error": "Invalid password"}), 401

    token = generate_jwt_token(user["email"], user["name"])

    # 3. STORE IN REDIS CACHE
    cache_set(f"user:{email}", {
        "name": user["name"],
        "email": user["email"],
        "password": user["password"]
    })

    return jsonify({
        "message": "Login successful",
        "token": token,
        "user": {"name": user["name"], "email": user["email"]}
    })

# ===============================
# SMILES → PDB LIGAND
# ===============================
@app.route("/convert_smiles", methods=["POST"])
@token_required
def convert_smiles():
    smiles = request.json.get("smiles")
    if not smiles:
        return jsonify({"error": "SMILES required"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES"}), 400

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    file_id = uuid.uuid4().hex
    pdb_path = os.path.join(LIGAND_GEN_DIR, f"{file_id}.pdb")
    Chem.MolToPDBFile(mol, pdb_path)

    return send_file(pdb_path, as_attachment=True)

# ===============================
# LOGP ML
# ===============================
@app.route("/predict_logp", methods=["POST"])
@token_required
def predict_logp():
    smiles = request.json.get("smiles")

    # 1. CHECK REDIS
    cached = cache_get(f"logp:{smiles}")
    if cached:
        return jsonify({"logp": cached["value"], "summary": "cached"})

    # 2. CHECK MONGO
    db_entry = get_chem_result(smiles)
    if db_entry and "logp" in db_entry:
        cache_set(f"logp:{smiles}", {"value": db_entry["logp"]})
        return jsonify({"logp": db_entry["logp"], "summary": "cached-db"})

    # 3. ORIGINAL PREDICTION
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES"}), 400
    mol = Chem.AddHs(mol)

    features = [
        Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol)
    ]

    value = float(logp_model.predict([features])[0])
    summary = (
        "Highly soluble" if value < 1 else
        "Moderately soluble" if value < 3 else
        "Low solubility"
    )

    # 4. STORE IN REDIS + MONGO
    cache_set(f"logp:{smiles}", {"value": value})
    save_chem_result(smiles, logp=value)

    return jsonify({"logp": value, "summary": summary})

# ===============================
# PIC50 ML
# ===============================
@app.route("/predict_pic50", methods=["POST"])
@token_required
def predict_pic50():
    smiles = request.json.get("smiles")

    # 1. REDIS
    cached = cache_get(f"pic50:{smiles}")
    if cached:
        return jsonify({"pic50": cached["value"], "remark": "cached"})

    # 2. MONGO
    db_entry = get_chem_result(smiles)
    if db_entry and "pic50" in db_entry:
        cache_set(f"pic50:{smiles}", {"value": db_entry["pic50"]})
        return jsonify({"pic50": db_entry["pic50"], "remark": "cached-db"})

    # 3. ORIGINAL PREDICTION
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES"}), 400
    mol = Chem.AddHs(mol)

    features = [
        Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
        Descriptors.NumRotatableBonds(mol), Descriptors.TPSA(mol)
    ]

    value = float(pic50_model.predict([features])[0])
    remark = (
        "Excellent inhibitor" if value > 9 else
        "Good inhibitor" if value > 7 else
        "Weak inhibitor"
    )

    # 4. STORE IN REDIS + MONGO
    cache_set(f"pic50:{smiles}", {"value": value})
    save_chem_result(smiles, pic50=value)

    return jsonify({"pic50": value, "remark": remark})

# ===============================
# Convert PDB → PDBQT
# ===============================
@app.route("/convert_pdb_to_pdbqt", methods=["POST"])
@token_required
def convert_pdb_to_pdbqt():
    file = request.files.get("file")
    if not file:
        return jsonify({"error": "File missing"}), 400

    pdb_path = os.path.join(OUTPUT_DIR, file.filename)
    file.save(pdb_path)

    pdbqt_name = file.filename.replace(".pdb", ".pdbqt")
    pdbqt_path = os.path.join(OUTPUT_DIR, pdbqt_name)

    OBABEL = r"C:\Users\Bellamkonda Sowmya\.conda\envs\docking\Library\bin\obabel.exe"

    subprocess.run([
        OBABEL,
        "-ipdb", pdb_path,
        "-opdbqt",
        "-O", pdbqt_path,
        "--gen3d"
    ], check=True)

    return jsonify({"ligand_filename": pdbqt_name})

# ===============================
# DOCKING
# ===============================

@app.route("/dock", methods=["POST"])
@token_required
def dock():
    protein = request.json.get("protein")
    ligand = request.json.get("ligand")

    task = run_docking_task.delay(protein, ligand)
    print("Flask received dock request!")
    result = task.get(timeout=600)

    return jsonify(result)

@app.route("/dock-status/<task_id>")
@token_required
def dock_status(task_id):
    task = run_docking_task.AsyncResult(task_id)

    if task.state == "PENDING":
        return {"status": "pending"}

    if task.state == "SUCCESS":
        return {
            "status": "completed",
            "result": task.result
        }

    return {"status": task.state}


# ===============================
# DOWNLOAD ROUTES
# ===============================
@app.route("/download/ligand/<name>")
@token_required
def download_ligand(name):
    path = os.path.join(OUTPUT_DIR, name)
    if not os.path.exists(path):
        return jsonify({"error": "Ligand not found"}), 404
    return send_file(path, as_attachment=True)

@app.route("/download/pose/<name>")
def download_pose(name):
    path = os.path.join(DOCKED_DIR, name)
    if not os.path.exists(path):
        return jsonify({"error": "Pose not found"}), 404
    return send_file(path, as_attachment=False)

@app.route("/download_protein", methods=["GET"])
def download_protein():
    protein = request.args.get("protein")
    filename = f"{protein}.pdbqt"
    file_path = os.path.join(RECEPTOR_DIR, filename)
    if not os.path.exists(file_path):
        return jsonify({"error": "Protein file not found"}), 404
    return send_file(file_path, as_attachment=False)

# ===============================
# Run Server
# ===============================
if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)
