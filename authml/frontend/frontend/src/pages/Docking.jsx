import React, { useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import {
  Box,
  Button,
  Card,
  Typography,
  TextField,
  Select,
  MenuItem,
  InputLabel,
  FormControl,
  CircularProgress,
} from "@mui/material";
import axios from "axios";

/* ------------------------- NAVBAR COMPONENT ------------------------- */
const Navbar = ({ onLogout }) => (
  <Box
    sx={{
      display: "flex",
      alignItems: "center",
      justifyContent: "space-between",
      p: 2,
      px: "5%",
      background: "rgba(255,255,255,0.12)",
      backdropFilter: "blur(10px)",
      position: "sticky",
      top: 0,
      zIndex: 1000,
    }}
  >
    <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
      <img src="/LB_logo.png" alt="LB Logo" style={{ width: 60, height: 60 }} />
      <Typography variant="h4" sx={{ color: "#fff", fontWeight: "bold" }}>
        LigBind
      </Typography>
    </Box>

    <Button
      variant="outlined"
      color="inherit"
      onClick={onLogout}
      sx={{
        borderRadius: "25px",
        borderWidth: 2,
        color: "#fff",
        px: 3,
        py: 1,
        fontWeight: 600,
        "&:hover": { backgroundColor: "rgba(255,255,255,0.2)" },
      }}
    >
      Logout
    </Button>
  </Box>
);

export default function Docking() {
  const { state } = useLocation();
  const navigate = useNavigate();
  const token = localStorage.getItem("token");

  const [protein, setProtein] = useState("Thrombin");
  const [smiles, setSmiles] = useState(state?.smiles || "");
  const [ligand, setLigand] = useState(null);
  const [loading, setLoading] = useState(false);

  const handleLogout = () => {
    localStorage.removeItem("token");
    navigate("/");
  };

  /* --------------------- DOWNLOAD PROTEIN --------------------- */
  const downloadProtein = () => {
    const url = `http://127.0.0.1:5000/download_protein?protein=${protein}`;
    fetch(url, {
      method: "GET",
      headers: { Authorization: `Bearer ${token}` },
    })
      .then((res) => {
        if (!res.ok) throw new Error();
        return res.blob();
      })
      .then((blob) => {
        const fileURL = window.URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = fileURL;
        a.download = `${protein}.pdbqt`;
        a.click();
      })
      .catch(() => alert("Protein download failed"));
  };

  /* --------------------- GENERATE LIGAND --------------------- */
  const generateLigand = async () => {
    if (!smiles) return alert("Enter SMILES first");
    setLoading(true);

    try {
      const res = await fetch("http://127.0.0.1:5000/convert_smiles", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({ smiles }),
      });

      if (!res.ok) {
        setLoading(false);
        return alert("Ligand generation failed");
      }

      const blob = await res.blob();
      const ligandName = `ligand_${Date.now()}.pdb`;

      // Download .pdb
      setLigand(ligandName);
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = ligandName;
      a.click();

      // Convert .pdb → .pdbqt
      const formData = new FormData();
      formData.append("file", blob, ligandName);

      const convertRes = await fetch("http://127.0.0.1:5000/convert_pdb_to_pdbqt", {
        method: "POST",
        headers: { Authorization: `Bearer ${token}` },
        body: formData,
      });

      const json = await convertRes.json();
      setLigand(json.ligand_filename);
      alert("Ligand (.pdbqt) generated and ready for docking!");
    } catch {
      alert("Ligand generation failed");
    }
    setLoading(false);
  };

  /* --------------------- RUN DOCKING --------------------- */
  const runDocking = async () => {
    if (!ligand) return alert("Generate ligand first");
    setLoading(true);

    try {
      const res = await axios.post(
        "http://127.0.0.1:5000/dock",
        { protein, ligand },
        { headers: { Authorization: `Bearer ${token}` } }
      );

      navigate("/docking-result", {
        state: {
          affinity: res.data.affinity,
          summary: res.data.summary,
          pose_file: res.data.pose_file,
          protein,
        },
      });
    } catch {
      alert("Docking failed");
    }
    setLoading(false);
  };

  return (
    <Box sx={{ position: "relative", minHeight: "100vh", backgroundColor: "#0c0a21" }}>
      <Box
        sx={{
          position: "absolute",
          inset: 0,
          backgroundImage: 'url("/LB.png")',
          backgroundSize: "cover",
          backgroundPosition: "center",
          filter: "blur(6px) brightness(0.6)",
          zIndex: 0,
        }}
      />

      <Navbar onLogout={handleLogout} />

      <Box
        sx={{
          position: "relative",
          zIndex: 1,
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          minHeight: "calc(100vh - 90px)",
          px: 2,
        }}
      >
        <Card
          sx={{
            width: "100%",
            maxWidth: 600,
            p: 4,
            borderRadius: 3,
            backgroundColor: "rgba(26,26,46,0.85)",
            boxShadow: "0 8px 32px rgba(0,0,0,0.6)",
            display: "flex",
            flexDirection: "column",
            gap: 2,
            color: "white",
          }}
        >
          <Typography variant="h4" fontWeight="bold" textAlign="center" mb={2}>
            Protein–Ligand Docking
          </Typography>

          <FormControl fullWidth>
            <InputLabel sx={{ color: "white" }}>Select Protein</InputLabel>
            <Select
              value={protein}
              onChange={(e) => setProtein(e.target.value)}
              sx={{ color: "white", "& .MuiSvgIcon-root": { color: "white" } }}
            >
              <MenuItem value="Thrombin">Thrombin</MenuItem>
              <MenuItem value="Protease">Protease</MenuItem>
              <MenuItem value="Kinase">Kinase</MenuItem>
            </Select>
          </FormControl>

          <Button
            variant="contained"
            sx={{ backgroundColor: "#764ba2", "&:hover": { backgroundColor: "#667eea" } }}
            onClick={downloadProtein}
          >
            ⬇ Download Protein (.pdbqt)
          </Button>

          <TextField
            label="Enter SMILES"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            fullWidth
            InputLabelProps={{ style: { color: "white" } }}
            InputProps={{ style: { color: "white" } }}
          />

          <Button
            variant="contained"
            sx={{ backgroundColor: "#2196f3", "&:hover": { backgroundColor: "#1976d2" } }}
            onClick={generateLigand}
            disabled={loading}
          >
            {loading ? <CircularProgress size={24} color="inherit" /> : "Convert SMILES → Ligand (.pdb)"}
          </Button>

          {ligand && (
            <Button
              variant="contained"
              sx={{ backgroundColor: "#00e676", "&:hover": { backgroundColor: "#00c853" } }}
              onClick={runDocking}
              disabled={loading}
            >
              {loading ? <CircularProgress size={24} color="inherit" /> : "Run Docking"}
            </Button>
          )}
        </Card>
      </Box>
    </Box>
  );
}
