import os
import uuid
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_ligand(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise Exception("Invalid SMILES")

    mol = Chem.AddHs(mol)

    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise Exception("3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol)

    os.makedirs("ligands", exist_ok=True)

    uid = str(uuid.uuid4())
    pdb_path = f"ligands/{uid}.pdb"
    pdbqt_path = f"ligands/{uid}.pdbqt"

    Chem.MolToPDBFile(mol, pdb_path)

    # 👇 Correct OpenBabel path (update if different)
    OBABEL = r"C:\Program Files\OpenBabel-3.1.1\obabel.exe"

    result = subprocess.run(
        [OBABEL, pdb_path, "-O", pdbqt_path], 
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        raise Exception("OpenBabel conversion failed: " + result.stderr)

    return {"pdb": pdb_path, "pdbqt": pdbqt_path}
