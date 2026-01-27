

import React, { useState, useEffect } from 'react';
import { Box, Button, Card, Typography, Alert, TextField } from '@mui/material';
import { useNavigate, useLocation } from 'react-router-dom';
import api from '../api/axiosConfig';

/* ------------------------- NAVBAR COMPONENT ------------------------- */
const Navbar = ({ onLoginClick, onSignupClick }) => (
  <Box
    sx={{
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'space-between',   // Logo left - buttons right
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
      <img src="/LB_logo.png" alt="LB Logo" style={{ width: 70, height: 70 }} />
      <Typography variant="h4" sx={{ color: '#fff', fontWeight: 'bold' }}>
        LigBind
      </Typography>
    </Box>

    {/* Buttons Right */}
    <Box sx={{ display: 'flex', gap: 2 }}>
      <Button
        variant="outlined"
        color="inherit"
        onClick={onLoginClick}
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
        Login
      </Button>

      <Button
        variant="contained"
        onClick={onSignupClick}
        sx={{
          borderRadius: '25px',
          px: 3,
          py: 1,
          fontWeight: 600,
          color: '#8e44ad',
          backgroundColor: '#ffffff',
          boxShadow: '0 4px 15px rgba(0,0,0,0.2)',
          '&:hover': {
            transform: 'translateY(-2px)',
            boxShadow: '0 6px 20px rgba(0,0,0,0.3)',
          },
        }}
      >
        Signup
      </Button>
    </Box>
  </Box>
);

/* ------------------------- LOGIN PAGE ------------------------- */
const Login = () => {
  const navigate = useNavigate();
  const location = useLocation();

  const [formData, setFormData] = useState({ email: '', password: '' });
  const [error, setError] = useState('');
  const [success, setSuccess] = useState('');
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    if (location.state?.message) {
      setSuccess(location.state.message);
    }
  }, [location]);

  const handleChange = (e) => {
    setFormData({ ...formData, [e.target.name]: e.target.value });
    setError('');
  };

  const inputStyle = {
    input: { color: 'white' },
    label: { color: 'white' },
    '& .MuiOutlinedInput-root': {
      '& fieldset': { borderColor: 'white' },
      '&:hover fieldset': { borderColor: '#764ba2' },
      '&.Mui-focused fieldset': { borderColor: '#764ba2' },
    },
    mb: 2,
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    setError('');
    setSuccess('');
    setLoading(true);

    try {
      const response = await api.post('/login', formData);

      localStorage.setItem('token', response.data.token);
      localStorage.setItem('user', JSON.stringify(response.data.user));

      navigate('/home');
    } catch (err) {
      setError(err.response?.data?.error || 'Login failed. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <Box sx={{ position: 'relative', minHeight: '100vh', backgroundColor: '#0c0a21' }}>
      {/* Blurred Background Image */}
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

      {/* Content Layer */}
      <Box sx={{ position: 'relative', zIndex: 1 }}>
        <Navbar
          onLoginClick={() => navigate('/login')}
          onSignupClick={() => navigate('/signup')}
        />

        {/* Center Card */}
        <Box
          sx={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            minHeight: 'calc(100vh - 90px)',
            px: 2,
          }}
        >
          <Card
            sx={{
              maxWidth: 450,
              width: '100%',
              p: 4,
              borderRadius: 2,
              boxShadow: '0 8px 32px rgba(0,0,0,0.6)',
              backgroundColor: 'rgba(26,26,46,0.85)',
              color: 'white',
            }}
          >
            <Typography variant="h4" textAlign="center" fontWeight="bold" mb={1}>
              Welcome Back
            </Typography>
            <Typography variant="body2" textAlign="center" color="gray" mb={3}>
              Login to continue to LigBind
            </Typography>

            {success && <Alert severity="success" sx={{ mb: 2 }}>{success}</Alert>}
            {error && <Alert severity="error" sx={{ mb: 2 }}>{error}</Alert>}

            <form onSubmit={handleSubmit}>
              <TextField
                label="Email Address"
                name="email"
                type="email"
                value={formData.email}
                onChange={handleChange}
                required
                fullWidth
                sx={inputStyle}
              />

              <TextField
                label="Password"
                name="password"
                type="password"
                value={formData.password}
                onChange={handleChange}
                required
                fullWidth
                sx={inputStyle}
              />

              <Button
                type="submit"
                variant="contained"
                fullWidth
                size="large"
                disabled={loading}
                sx={{
                  mt: 2,
                  py: 1.5,
                  textTransform: 'none',
                  fontSize: '16px',
                  fontWeight: 600,
                  backgroundColor: '#764ba2',
                  '&:hover': { backgroundColor: '#667eea' },
                }}
              >
                {loading ? 'Logging in..a.' : 'Login'}
              </Button>
            </form>

            <Typography variant="body2" textAlign="center" mt={3} color="gray">
              Don't have an account?{' '}
              <Typography
                component="span"
                color="#764ba2"
                sx={{ cursor: 'pointer', fontWeight: 600 }}
                onClick={() => navigate('/signup')}
              >
                Sign Up
              </Typography>
            </Typography>
          </Card>
        </Box>
      </Box>
    </Box>
  );
};

export default Login;