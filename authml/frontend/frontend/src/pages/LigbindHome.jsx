

import React, { useEffect, useState } from 'react';
import { Box, Card, Typography, Button, Avatar } from '@mui/material';
import { useNavigate } from 'react-router-dom';

/* ------------------------- NAVBAR COMPONENT ------------------------- */
const Navbar = ({ onLogout }) => (
  <Box
    sx={{
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'space-between',
      p: 2,
      px: '5%',
      background: 'rgba(255,255,255,0.12)',
      backdropFilter: 'blur(10px)',
      position: 'sticky',
      top: 0,
      zIndex: 1000,
    }}
  >
    {/* Logo + Title */}
    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
      <img src="/LB_logo.png" alt="LB Logo" style={{ width: 60, height: 60 }} />
      <Typography variant="h4" sx={{ color: '#fff', fontWeight: 'bold' }}>
        LigBind
      </Typography>
    </Box>

    {/* Logout */}
    <Button
      variant="outlined"
      color="inherit"
      onClick={onLogout}
      sx={{
        borderRadius: '25px',
        borderWidth: 2,
        color: '#fff',
        px: 3,
        py: 1,
        fontWeight: 600,
        '&:hover': { backgroundColor: 'rgba(255,255,255,0.2)' },
      }}
    >
      Logout
    </Button>
  </Box>
);

/* ------------------------- HOME PAGE ------------------------- */
const LigbindHome = () => {
  const navigate = useNavigate();
  const [user, setUser] = useState(null);

  useEffect(() => {
    const userData = localStorage.getItem('user');
    if (!userData) {
      navigate('/login');
    } else {
      setUser(JSON.parse(userData));
    }
  }, [navigate]);

  const handleLogout = () => {
    localStorage.removeItem('user');
    navigate('/');
  };

  if (!user) return null;

  return (
    <Box sx={{ position: 'relative', minHeight: '100vh', backgroundColor: '#0c0a21' }}>
      
      {/* 🔥 Background Blur Image */}
      <Box
        sx={{
          position: 'absolute',
          inset: 0,
          backgroundImage: 'url("/LB.png")',
          backgroundSize: 'cover',
          backgroundPosition: 'center',
          filter: 'blur(6px) brightness(0.6)',
          zIndex: 0,
        }}
      />

      {/* Main Content */}
      <Box sx={{ position: 'relative', zIndex: 1 }}>
        <Navbar onLogout={handleLogout} />

        <Box
          sx={{
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
            minHeight: 'calc(100vh - 90px)',
            px: 3,
          }}
        >
          <Card
            sx={{
              width: '100%',
              maxWidth: 600,
              p: 4,
              borderRadius: 3,
              textAlign: 'center',
              backgroundColor: 'rgba(26, 26, 46, 0.85)',
              boxShadow: '0 8px 32px rgba(0,0,0,0.6)',
              color: 'white',
            }}
          >
            <Avatar
              sx={{
                width: 100,
                height: 100,
                margin: '0 auto 20px',
                bgcolor: '#764ba2',
                fontSize: '40px',
              }}
            >
              {user.name.charAt(0).toUpperCase()}
            </Avatar>

            <Typography variant="h4" fontWeight="bold" mb={1}>
              Welcome to LigBind!
            </Typography>

            <Typography variant="h6" color="gray" mb={1}>
              Hello, {user.name}
            </Typography>

            <Typography variant="body2" color="gray" mb={4}>
              {user.email}
            </Typography>

            {/* Cards Section */}
            <Box sx={{ display: 'flex', gap: 2, justifyContent: 'center', flexWrap: 'wrap' }}>
              <Card
                sx={{
                  padding: 3,
                  minWidth: 250,
                  backgroundColor: 'rgba(255,255,255,0.08)',
                  backdropFilter: 'blur(10px)',
                  color: '#fff',
                  boxShadow: 3,
                }}
              >
                <Typography variant="h6" fontWeight="bold" color="#764ba2" mb={1}>
                  Account Status
                </Typography>
                <Typography variant="body2">✓ Email Verified</Typography>
                <Typography variant="body2">✓ Account Active</Typography>
              </Card>

              <Card
                sx={{
                  padding: 3,
                  minWidth: 250,
                  backgroundColor: 'rgba(255,255,255,0.08)',
                  backdropFilter: 'blur(10px)',
                  color: '#fff',
                  boxShadow: 3,
                }}
              >
                <Typography variant="h6" fontWeight="bold" color="#764ba2" mb={1}>
                  Quick Actions
                </Typography>
                <Typography variant="body2" mb={1}>• View Dashboard</Typography>
                <Typography variant="body2" mb={1}>• Manage Projects</Typography>
                <Typography variant="body2">• Settings</Typography>
              </Card>
            </Box>

            {/* 🔬 Start Prediction Button */}
            <Button
              variant="contained"
              onClick={() => navigate('/predict')}
              sx={{
                mt: 4,
                px: 4,
                py: 1.5,
                textTransform: 'none',
                fontSize: '17px',
                fontWeight: 600,
                borderRadius: 2,
                backgroundColor: '#764ba2',
                '&:hover': { backgroundColor: '#667eea' },
              }}
            >
              🔬 Start Prediction
            </Button>
          </Card>
        </Box>
      </Box>
    </Box>
  );
};

export default LigbindHome;