import axios from "axios";
const API = "http://127.0.0.1:5000";

export const predictLogP = (smiles) => axios.post(`${API}/predict_logp`, { smiles });
export const predictPIC50 = (smiles) => axios.post(`${API}/predict_pic50`, { smiles });
export const getProtein = (protein) => axios.post(`${API}/get_protein`, { protein });
export const generateLigand = (smiles) => axios.post(`${API}/generate_ligand`, { smiles });
export const dock = (protein) => axios.post(`${API}/dock`, { protein });
