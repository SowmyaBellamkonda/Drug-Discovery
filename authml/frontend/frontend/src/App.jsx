import React from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate } from 'react-router-dom';
import { ThemeProvider, createTheme } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import Signup from './pages/Signup';
import OtpVerify from './pages/OtpVerify';
import Login from './pages/Login';
import LigbindHome from './pages/LigbindHome';
import LigBindLanding from './pages/LigBindLanding';
import DockingResult from './pages/DockingResult';
import Predict from './pages/Predict';
import LogPResult from './pages/LogPResult';
import PIC50Result from './pages/PIC50Result';
import Docking from './pages/Docking';

const theme = createTheme({
  palette: {
    primary: {
      main: '#667eea',
    },
    secondary: {
      main: '#764ba2',
    },
  },
  typography: {
    fontFamily: '"Roboto", "Helvetica", "Arial", sans-serif',
  },
});

function App() {
  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <Router>
        <Routes>
          {/* <Route path="/" element={<Navigate to="/LigBindLanding" replace />} /> */}
          <Route path="/" element={<LigBindLanding/>} />
          <Route path="/signup" element={<Signup />} />
          <Route path="/otp-verify" element={<OtpVerify />} />
          <Route path="/login" element={<Login />} />
          <Route path="/home" element={<LigbindHome />} />
          <Route path="/predict" element={<Predict />} />
          <Route path="/logp" element={<LogPResult />} />
          <Route path="/pic50" element={<PIC50Result />} />
          <Route path="/docking" element={<Docking />} />
          <Route path="/docking-result" element={<DockingResult />} />
        </Routes>
      </Router>
    </ThemeProvider>
  );
}

export default App;