import React, { useEffect, useRef } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import * as NGL from "ngl";
import {
  Box,
  Button,
  Card,
  CardContent,
  Typography,
  Table,
  TableBody,
  TableRow,
  TableCell,
} from "@mui/material";

/* ------------------------- NAVBAR COMPONENT ------------------------- */
const Navbar = ({ onLogout }) => (
  <Box
    sx={{
      display: "flex",
      justifyContent: "space-between",
      alignItems: "center",
      p: 2,
      px: "5%",
      background: "rgba(255,255,255,0.1)",
      backdropFilter: "blur(10px)",
      position: "sticky",
      top: 0,
      zIndex: 1000,
    }}
  >
    <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
      <img src="/LB_logo.png" alt="logo" style={{ width: 70, height: 70 }} />
      <Typography variant="h4" sx={{ fontWeight: "bold", color: "#fff" }}>
        LigBind
      </Typography>
    </Box>

    <Button
      variant="contained"
      onClick={onLogout}
      sx={{
        borderRadius: "25px",
        px: 3,
        py: 1,
        fontWeight: 600,
        color: "#8e44ad",
        backgroundColor: "#fdfdfd",
        boxShadow: "0 4px 15px rgba(0,0,0,0.2)",
        "&:hover": {
          transform: "translateY(-2px)",
          boxShadow: "0 6px 20px rgba(0,0,0,0.3)",
        },
      }}
    >
      Logout
    </Button>
  </Box>
);

export default function DockingResult() {
  const { state } = useLocation();
  const navigate = useNavigate();
  const viewerRef = useRef(null);

  const affinity = state?.affinity;
  const poseFile = state?.pose_file;
  const protein = state?.protein;
  const summary = state?.summary;
  const token = localStorage.getItem("token");

  useEffect(() => {
    if (!poseFile || !protein) return;

    const viewerDiv = document.getElementById("viewer");
    viewerDiv.innerHTML = "";

    const stage = new NGL.Stage("viewer", { backgroundColor: "black" });
    const timestamp = Date.now();

    const headers = { Authorization: `Bearer ${token}` };

    stage
      .loadFile(
        `http://127.0.0.1:5000/download_protein?protein=${protein}&t=${timestamp}`,
        { ext: "pdbqt", fetchParams: { headers } }
      )
      .then((o) => {
        o.addRepresentation("cartoon", { color: "skyblue", opacity: 0.8 });
        stage.autoView();
      });

    stage
      .loadFile(
        `http://127.0.0.1:5000/download/pose/${poseFile}?t=${timestamp}`,
        { ext: "pdbqt", fetchParams: { headers } }
      )
      .then((o) => {
        o.addRepresentation("ball+stick", {
          colorScheme: "element",
          scale: 2.0,
          bondScale: 1.0,
        });
        stage.autoView();
      });

    stage.setSpin(true);
    viewerRef.current = stage;

    const handleResize = () => stage.handleResize();
    window.addEventListener("resize", handleResize);

    return () => {
      window.removeEventListener("resize", handleResize);
      stage.dispose();
    };
  }, [poseFile, protein, token]);

  const handleLogout = () => {
    localStorage.removeItem("token");
    navigate("/");
  };

  const captureImage = () => {
    if (!viewerRef.current) return;
    viewerRef.current
      .makeImage({
        factor: 2,
        antialias: true,
        trim: false,
        transparent: false,
      })
      .then((blob) => {
        const link = document.createElement("a");
        link.href = URL.createObjectURL(blob);
        link.download = "docking_view.png";
        link.click();
      });
  };

  return (
    <Box sx={{ position: "relative", minHeight: "100vh", backgroundColor: "#0c0a21" }}>
      {/* Background */}
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
      <Box sx={{ position: "relative", zIndex: 1, px: 2 }}>
        <Navbar onLogout={handleLogout} />

        <Box sx={{ maxWidth: 1200, mx: "auto", py: 6 }}>
          <Card
            sx={{
              p: 4,
              borderRadius: 3,
              backgroundColor: "rgba(26,26,46,0.85)",
              color: "#fff",
              boxShadow: "0 10px 30px rgba(0,0,0,0.5)",
            }}
          >
            <Typography variant="h4" sx={{ textAlign: "center", mb: 4, fontWeight: "bold" }}>
              🔬 Docking Result
            </Typography>

            <Card sx={{ mb: 4, borderRadius: 2, backgroundColor: "rgba(255,255,255,0.1)" }}>
              <CardContent>
                <Typography variant="h5">Binding Affinity</Typography>
                <Typography variant="h4" sx={{ mt: 1, fontWeight: "bold", color: "#b998ff" }}>
                  {affinity} kcal/mol
                </Typography>
                <Typography sx={{ mt: 1, opacity: 0.8 }}>{summary}</Typography>
              </CardContent>
            </Card>

            {/* Strength Scale Table */}
            <Card
              sx={{
                mb: 4,
                borderRadius: 2,
                border: "1px solid rgba(200,160,255,0.3)",
                backgroundColor: "rgba(255,255,255,0.05)",
              }}
            >
              <CardContent>
                <Typography variant="h6" sx={{ mb: 2 }}>
                  💡 Binding Affinity Strength Scale
                </Typography>
                <Table>
                  <TableBody>
                    <TableRow>
                      <TableCell sx={{ color: "white" }}>&lt; -9</TableCell>
                      <TableCell sx={{ color: "#b998ff" }}>Excellent</TableCell>
                    </TableRow>
                    <TableRow>
                      <TableCell sx={{ color: "white" }}>-7 to -9</TableCell>
                      <TableCell sx={{ color: "#b998ff" }}>Good</TableCell>
                    </TableRow>
                    <TableRow>
                      <TableCell sx={{ color: "white" }}>-5 to -7</TableCell>
                      <TableCell sx={{ color: "#b998ff" }}>Moderate</TableCell>
                    </TableRow>
                    <TableRow>
                      <TableCell sx={{ color: "white" }}>&gt; -5</TableCell>
                      <TableCell sx={{ color: "#b998ff" }}>Weak</TableCell>
                    </TableRow>
                  </TableBody>
                </Table>
              </CardContent>
            </Card>

            {/* NGL Viewer */}
            <Box
              sx={{
                height: "500px",
                width: "100%",
                borderRadius: 3,
                overflow: "hidden",
                mb: 4,
                border: "1px solid rgba(255,255,255,0.2)",
              }}
            >
              <Box id="viewer" sx={{ width: "100%", height: "100%" }} />
            </Box>

            {/* Buttons */}
            <Box sx={{ display: "flex", justifyContent: "center", gap: 3, flexWrap: "wrap" }}>
              <Button
                variant="contained"
                sx={{
                  borderRadius: "25px",
                  px: 4,
                  py: 1.5,
                  fontWeight: "bold",
                  bgcolor: "#764ba2",
                  "&:hover": { transform: "scale(1.05)", bgcolor: "#6b3fa3" },
                }}
                onClick={captureImage}
              >
                📸 Save Image
              </Button>

              <Button
                variant="contained"
                sx={{
                  borderRadius: "25px",
                  px: 4,
                  py: 1.5,
                  fontWeight: "bold",
                  bgcolor: "#667eea",
                  "&:hover": { transform: "scale(1.05)", bgcolor: "#5b75e0" },
                }}
                onClick={() =>
                  window.open(
                    `http://127.0.0.1:5000/download/pose/${poseFile}`,
                    "_blank"
                  )
                }
              >
                ⬇ Download Docked Pose
              </Button>

              <Button
                variant="contained"
                sx={{
                  borderRadius: "25px",
                  px: 4,
                  py: 1.5,
                  fontWeight: "bold",
                  bgcolor: "#555",
                  "&:hover": { transform: "scale(1.05)", bgcolor: "#444" },
                }}
                onClick={() => navigate("/docking")}
              >
                ↩ Back
              </Button>
            </Box>
          </Card>
        </Box>
      </Box>
    </Box>
  );
}
