import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors

# Load models
logp_model = joblib.load("logp_model.pkl")
pic50_model = joblib.load("pic50_model.pkl")

# Descriptors (same as training & app.py)
def desc_logp(mol):
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
    ]

def desc_pic50(mol):
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.TPSA(mol),
    ]


def predict_values(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    
    logp_desc = desc_logp(mol)
    pic50_desc = desc_pic50(mol)

    logp_val = float(logp_model.predict([logp_desc])[0])
    pic50_val = float(pic50_model.predict([pic50_desc])[0])

    return logp_val, pic50_val


if __name__ == "__main__":
    smi = input("enter smiles: ")

    logp, pic50 = predict_values(smi)

    if logp is None:
        print("\nInvalid SMILES")
    else:
        print("\nvalues")
        print(f"logP   = {round(logp, 4)}")
        print(f"pIC50  = {round(pic50, 4)}")