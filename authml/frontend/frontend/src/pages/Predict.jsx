import React, { useState, useEffect } from "react";
import axios from "axios";
import { useNavigate } from "react-router-dom";
import {
  Box,
  Paper,
  Typography,
  TextField,
  Button,
  Alert,
  CircularProgress,
  Card,
} from "@mui/material";

/* ------------------------- NAVBAR COMPONENT ------------------------- */
const Navbar = ({ onBack, onLogout }) => (
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
    {/* Logo + Title */}
    <Box
      sx={{ display: "flex", alignItems: "center", gap: 1, cursor: "pointer" }}
      onClick={onBack}
    >
      <img src="/LB_logo.png" alt="LB Logo" style={{ width: 60, height: 60 }} />
      <Typography variant="h4" sx={{ color: "#fff", fontWeight: "bold" }}>
        Predict
      </Typography>
    </Box>

    {/* Buttons */}
    <Box sx={{ display: "flex", gap: 2 }}>
      <Button
        variant="outlined"
        color="inherit"
        onClick={onBack}
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
        ⬅ Back
      </Button>

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
  </Box>
);

/* ------------------------- PREDICTION PAGE ------------------------- */
const Predict = () => {
  const navigate = useNavigate();
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");

  useEffect(() => {
    const token = localStorage.getItem("token");
    if (!token) navigate("/login");
  }, []);

  const handleLogout = () => {
    localStorage.removeItem("token");
    localStorage.removeItem("user");
    navigate("/");
  };

  const predictLogP = async () => {
    if (!smiles.trim()) {
      setError("Please enter a valid SMILES string");
      return;
    }

    const token = localStorage.getItem("token");
    if (!token) {
      setError("Authentication token missing. Please login again.");
      return;
    }

    setLoading(true);
    setError("");

    try {
      const response = await axios.post(
        "http://127.0.0.1:5000/predict_logp",
        { smiles: smiles.trim() },
        {
          headers: { Authorization: `Bearer ${token}` },
        }
      );

      navigate("/logp", {
        state: {
          logp: response.data.logp,
          remark: response.data.remark,
          smiles: smiles.trim(),
        },
      });
    } catch (err) {
      setError(err.response?.data?.error || "LogP prediction failed.");
    } finally {
      setLoading(false);
    }
  };

  const predictPIC50 = async () => {
    if (!smiles.trim()) {
      setError("Please enter a valid SMILES string");
      return;
    }

    const token = localStorage.getItem("token");
    if (!token) {
      setError("Authentication token missing. Please login again.");
      return;
    }

    setLoading(true);
    setError("");

    try {
      const response = await axios.post(
        "http://127.0.0.1:5000/predict_pic50",
        { smiles: smiles.trim() },
        {
          headers: { Authorization: `Bearer ${token}` },
        }
      );

      navigate("/pic50", {
        state: {
          pic50: response.data.pic50,
          remark: response.data.remark,
          smiles: smiles.trim(),
        },
      });
    } catch (err) {
      setError(err.response?.data?.error || "pIC50 prediction failed.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <Box sx={{ position: "relative", minHeight: "100vh", backgroundColor: "#0c0a21" }}>
      {/* Blurred Background */}
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

      {/* Main Content */}
      <Box sx={{ position: "relative", zIndex: 1 }}>
        <Navbar onBack={() => navigate("/home")} onLogout={handleLogout} />

        <Box
          sx={{
            display: "flex",
            justifyContent: "center",
            alignItems: "center",
            minHeight: "calc(100vh - 90px)",
            px: 3,
          }}
        >
          <Card
            sx={{
              width: "100%",
              maxWidth: 550,
              p: 4,
              borderRadius: 3,
              textAlign: "center",
              backgroundColor: "rgba(26, 26, 46, 0.85)",
              boxShadow: "0 8px 32px rgba(0,0,0,0.6)",
              color: "white",
            }}
          >
            <Typography variant="h4" fontWeight="bold" mb={2}>
              Predict LogP & pIC₅₀
            </Typography>

            <Typography variant="body1" color="gray" mb={3}>
              Enter a SMILES string to compute molecular properties
            </Typography>

            <TextField
              label="Enter SMILES String"
              fullWidth
              value={smiles}
              onChange={(e) => {
                setSmiles(e.target.value);
                setError("");
              }}
              sx={{
                mb: 3,
                input: { color: "white" },
                label: { color: "white" },
                "& .MuiOutlinedInput-root": {
                  "& fieldset": { borderColor: "white" },
                  "&:hover fieldset": { borderColor: "#764ba2" },
                  "&.Mui-focused fieldset": { borderColor: "#764ba2" },
                },
              }}
            />

            {error && (
              <Alert severity="error" sx={{ mb: 2 }}>
                {error}
              </Alert>
            )}

            {smiles && !error && (
              <Alert severity="info" sx={{ mb: 2 }}>
                🔍 <strong>SMILES Summary:</strong> Used to compute solubility,
                inhibition strength and docking affinity.
              </Alert>
            )}

            <Button
              variant="contained"
              fullWidth
              disabled={loading}
              onClick={predictLogP}
              sx={{
                mb: 2,
                py: 1.2,
                backgroundColor: "#764ba2",
                "&:hover": { backgroundColor: "#667eea" },
              }}
            >
              {loading ? <CircularProgress size={24} /> : "Predict LogP"}
            </Button>

            <Button
              variant="contained"
              color="secondary"
              fullWidth
              disabled={loading}
              onClick={predictPIC50}
              sx={{ py: 1.2 }}
            >
              {loading ? <CircularProgress size={24} /> : "Predict pIC₅₀"}
            </Button>
          </Card>
        </Box>
      </Box>
    </Box>
  );
};

export default Predict;
