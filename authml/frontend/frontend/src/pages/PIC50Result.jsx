import React, { useEffect, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import axios from "axios";
import {
  Box,
  Paper,
  Typography,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Button,
  Card,
} from "@mui/material";

/* ------------------------- NAVBAR COMPONENT ------------------------- */
const Navbar = ({ onBack }) => (
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
    <Box
      sx={{ display: "flex", alignItems: "center", gap: 1, cursor: "pointer" }}
      onClick={onBack}
    >
      <img src="/LB_logo.png" alt="LB Logo" style={{ width: 60, height: 60 }} />
      <Typography variant="h4" sx={{ color: "#fff", fontWeight: "bold" }}>
        pIC₅₀ Result
      </Typography>
    </Box>

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
  </Box>
);

/* ------------------------- MAIN RESULT PAGE ------------------------- */
const PIC50Result = () => {
  const { state } = useLocation();
  const navigate = useNavigate();

  const token = localStorage.getItem("token");
  const smiles = state?.smiles;

  const [pic50, setPic50] = useState(null);
  const [remark, setRemark] = useState("");

  useEffect(() => {
    if (!smiles) {
      navigate("/predict");
      return;
    }

    const getPrediction = async () => {
      try {
        const res = await axios.post(
          "http://127.0.0.1:5000/predict_pic50",
          { smiles },
          { headers: { Authorization: `Bearer ${token}` } }
        );

        setPic50(res.data.pic50);
        setRemark(res.data.remark);
      } catch {
        alert("pIC50 prediction failed");
      }
    };

    getPrediction();
  }, [smiles, navigate, token]);

  const finalRemark =
    remark ||
    (pic50 < 7
      ? "Weak inhibitor"
      : pic50 <= 9
      ? "Good inhibitor"
      : "Excellent inhibitor");

  return (
    <Box
      sx={{
        position: "relative",
        minHeight: "100vh",
        backgroundColor: "#0c0a21",
      }}
    >
      {/* Background Image Blur */}
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

      {/* MAIN UI */}
      <Box sx={{ position: "relative", zIndex: 1 }}>
        <Navbar onBack={() => navigate("/predict")} />

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
              maxWidth: 600,
              p: 4,
              borderRadius: 3,
              textAlign: "center",
              backgroundColor: "rgba(26, 26, 46, 0.85)",
              boxShadow: "0 8px 32px rgba(0,0,0,0.6)",
              color: "white",
            }}
          >
            <Typography variant="h4" fontWeight="bold" mb={2}>
              pIC₅₀ Prediction Result
            </Typography>

            <Typography variant="h5" fontWeight={600} color="#76b5ff" mb={1}>
              pIC₅₀ Value: {pic50 !== null ? pic50.toFixed(2) : "Predicting..."}
            </Typography>

            <Typography sx={{ fontSize: "18px", fontWeight: 500, mb: 3 }}>
              {pic50 !== null ? finalRemark : "Calculating interpretation..."}
            </Typography>

            <Typography variant="h6" sx={{ fontWeight: 600, textAlign: "left", mb: 1 }}>
              pIC₅₀ Range Table
            </Typography>

            {/* Table */}
            <TableContainer
              component={Paper}
              sx={{
                borderRadius: 2,
                mb: 3,
                backgroundColor: "rgba(20,20,45,0.9)",
                border: "2px solid #111",
              }}
            >
              <Table
                sx={{
                  "& th, & td": {
                    border: "none",
                    color: "white",
                    fontSize: "15px",
                  },
                  "& th": {
                    backgroundColor: "rgba(30, 30, 55, 0.95)",
                    fontWeight: 700,
                  },
                  "& td": {
                    backgroundColor: "rgba(15,15,35,0.9)",
                  },
                  "& tr:hover td": {
                    backgroundColor: "rgba(45,45,70,0.9)",
                  },
                }}
              >
                <TableHead>
                  <TableRow>
                    <TableCell>pIC₅₀ Range</TableCell>
                    <TableCell align="right">Interpretation</TableCell>
                  </TableRow>
                </TableHead>

                <TableBody>
                  <TableRow>
                    <TableCell>pIC₅₀ &lt; 7</TableCell>
                    <TableCell align="right">Weak inhibitor</TableCell>
                  </TableRow>
                  <TableRow>
                    <TableCell>7 – 9</TableCell>
                    <TableCell align="right">Good inhibitor</TableCell>
                  </TableRow>
                  <TableRow>
                    <TableCell>&gt; 9</TableCell>
                    <TableCell align="right">Excellent inhibitor</TableCell>
                  </TableRow>
                </TableBody>
              </Table>
            </TableContainer>

            {/* ACTION BUTTONS */}
            <Box sx={{ mt: 2, display: "flex", gap: 2 }}>
              <Button
                variant="contained"
                fullWidth
                sx={{
                  py: 1.3,
                  backgroundColor: "#764ba2",
                  "&:hover": { backgroundColor: "#667eea" },
                }}
                onClick={() => navigate("/predict")}
              >
                ← Back
              </Button>

              <Button
                variant="contained"
                fullWidth
                sx={{
                  py: 1.3,
                  backgroundColor: "#2196f3",
                  "&:hover": { backgroundColor: "#1976d2" },
                }}
                onClick={() => navigate("/docking", { state: { smiles } })}
              >
                Proceed to Docking →
              </Button>
            </Box>
          </Card>
        </Box>
      </Box>
    </Box>
  );
};

export default PIC50Result;
