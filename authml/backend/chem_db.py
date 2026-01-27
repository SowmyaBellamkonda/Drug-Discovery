from db import chem_results_collection

def get_chem_result(smiles):
    return chem_results_collection.find_one({"smiles": smiles})

def save_chem_result(smiles, logp=None, pic50=None):
    update_data = {}
    if logp is not None:
        update_data["logp"] = logp
    if pic50 is not None:
        update_data["pic50"] = pic50

    chem_results_collection.update_one(
        {"smiles": smiles},
        {"$set": update_data},
        upsert=True
    )
