

import React from "react";
import { useLocation, useNavigate } from "react-router-dom";
import {
  Box,
  Card,
  Typography,
  Button,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Paper,
  TableContainer,
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
        LogP Result
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

/* ------------------------- RESULT PAGE ------------------------- */
const LogPResult = () => {
  const { state } = useLocation();
  const navigate = useNavigate();

  const logp = state?.logp;
  const remark = state?.remark;

  const finalRemark = remark
    ? remark
    : logp < 2
    ? "Poor absorption"
    : logp <= 5
    ? "Optimal absorption (Drug-like)"
    : "High toxicity risk";

  return (
    <Box sx={{ position: "relative", minHeight: "100vh", backgroundColor: "#0c0a21" }}>
      
      {/* Background Blur Image */}
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

      {/* Main UI */}
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
              LogP Prediction Result
            </Typography>

            <Typography variant="h5" fontWeight={600} color="#76b5ff" mb={1}>
              LogP Value: {logp ? logp.toFixed(2) : "—"}
            </Typography>

            <Typography sx={{ fontSize: "18px", fontWeight: 500, mb: 3 }}>
              {finalRemark}
            </Typography>

            {/* --------------------------------------- */}
            {/*           TABLE INTEGRATED HERE         */}
            {/* --------------------------------------- */}
            <TableContainer
              component={Paper}
              sx={{
                borderRadius: 2,
                mb: 3,
                backgroundColor: "rgba(20,20,45,0.9)",
                border: "2px solid #111", // outer border only
              }}
            >
              <Table
                sx={{
                  "& th, & td": {
                    border: "none", // remove inner borders
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
                    <TableCell>LogP Range</TableCell>
                    <TableCell align="right">Interpretation</TableCell>
                  </TableRow>
                </TableHead>

                <TableBody>
                  <TableRow>
                    <TableCell>LogP &lt; 2</TableCell>
                    <TableCell align="right">Poor absorption</TableCell>
                  </TableRow>

                  <TableRow>
                    <TableCell>2 – 5</TableCell>
                    <TableCell align="right">Optimal absorption (Drug-like)</TableCell>
                  </TableRow>

                  <TableRow>
                    <TableCell>&gt; 5</TableCell>
                    <TableCell align="right">High toxicity risk</TableCell>
                  </TableRow>
                </TableBody>
              </Table>
            </TableContainer>

            {/* Buttons */}
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
                onClick={() => navigate("/pic50", { state })}
              >
                Proceed to pIC₅₀ →
              </Button>
            </Box>
          </Card>
        </Box>
      </Box>
    </Box>
  );
};

export default LogPResult;