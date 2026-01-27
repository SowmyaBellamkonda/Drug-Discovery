

import React from 'react';
import { Box, Button, Card, Typography } from '@mui/material';
import { useNavigate } from 'react-router-dom';

// ✅ Navbar Updated (Login/Signup moved to RIGHT)
const Navbar = ({ onLoginClick, onSignupClick }) => (
  <Box
    sx={{
      display: 'flex',
      justifyContent: 'space-between', // LEFT + RIGHT SEPARATION
      alignItems: 'center',
      p: 2,
      px: '5%',
      background: 'rgba(255,255,255,0.1)',
      backdropFilter: 'blur(10px)',
      position: 'sticky',
      top: 0,
      zIndex: 1000,
    }}
  >
    {/* LEFT: LOGO + TITLE */}
    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
      <img src="/LB_logo.png" alt="logo" style={{ width: 70, height: 70 }} />
      <Typography
        variant="h4"
        sx={{ fontWeight: 'bold', color: '#fff', letterSpacing: '-1px' }}
      >
        LigBind
      </Typography>
    </Box>

    {/* RIGHT: LOGIN + SIGNUP */}
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
          backgroundColor: '#fdfdfd',
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

const LigBindLanding = () => {
  const navigate = useNavigate();

  const handleLoginClick = () => navigate('/login');
  const handleSignupClick = () => navigate('/signup');

  const diseases = [
    { title: 'Blood Clotting', description: 'Stop dangerous blood clots', details: 'Find medications that thin your blood safely.', image: '/blood-clot.png' },
    { title: 'Tissue Damage', description: 'Speed up healing', details: 'Discover drugs that help your body repair injuries.', image: '/tissue-repair.png' },
    { title: 'Cancer', description: 'Fight cancer', details: 'Access the latest cancer-fighting drugs safely.', image: '/cancer-cell.png' },
  ];

  const howItWorks = [
    { step: '1', title: 'Pick your condition', description: 'Choose blood clots, injuries, or cancer.', image: '/step1.png' },
    { step: '2', title: 'We find matches', description: 'Smart system searches thousands of drugs.', image: '/step2.png' },
    { step: '3', title: 'See results', description: 'Clear results on which drugs work best.', image: '/step3.png' },
  ];

  return (
    <Box sx={{ position: 'relative', minHeight: '100vh', backgroundColor: '#0c0a21' }}>
      
      {/* Background Image */}
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
        <Navbar onLoginClick={handleLoginClick} onSignupClick={handleSignupClick} />

        {/* Hero Section */}
        <Box sx={{ textAlign: 'center', py: 12, px: 2 }}>
          <Typography variant="h2" sx={{ color: '#fff', fontWeight: 'bold', mb: 2 }}>
            Find the Right Medicine
          </Typography>
          <Typography variant="h5" sx={{ color: '#fff', mb: 3 }}>
            Simple. Fast. Safe.
          </Typography>
          <Button
            onClick={handleSignupClick}
            sx={{
              px: 6,
              py: 2,
              borderRadius: '30px',
              fontSize: '1.3rem',
              fontWeight: 'bold',
              color: '#8e44ad',
              backgroundColor: '#fff',
              boxShadow: '0 8px 25px rgba(0,0,0,0.3)',
              '&:hover': { transform: 'scale(1.05)', boxShadow: '0 12px 35px rgba(0,0,0,0.4)' },
            }}
          >
            Start Finding Your Medicine →
          </Button>
        </Box>

        {/* Disease Cards */}
        <Box
          sx={{
            maxWidth: 1200,
            mx: 'auto',
            px: 2,
            py: 6,
            display: 'grid',
            gap: 4,
            gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))',
          }}
        >
          {diseases.map((disease, index) => (
            <Card
              key={index}
              sx={{
                p: 4,
                borderRadius: 3,
                textAlign: 'center',
                backgroundColor: 'rgba(26,26,46,0.85)',
                color: '#fff',
                cursor: 'pointer',
                transition: 'all 0.3s ease',
                '&:hover': {
                  transform: 'translateY(-10px) scale(1.03)',
                  boxShadow: '0 20px 50px rgba(0,0,0,0.4)',
                },
              }}
              onClick={handleSignupClick}
            >
              <img
                src={disease.image}
                alt={disease.title}
                style={{
                  width: '80%',
                  height: '200px',
                  objectFit: 'cover',
                  borderRadius: '12px',
                  marginBottom: '1.5rem',
                }}
              />
              <Typography variant="h5" sx={{ mb: 1, fontWeight: 'bold' }}>
                {disease.title}
              </Typography>
              <Typography sx={{ mb: 1 }}>{disease.description}</Typography>
              <Typography sx={{ fontSize: '0.9rem', opacity: 0.8 }}>
                {disease.details}
              </Typography>
            </Card>
          ))}
        </Box>

        {/* How It Works Section */}
        <Box sx={{ maxWidth: 1200, mx: 'auto', px: 2, py: 6 }}>
          <Card sx={{ p: 6, borderRadius: 3, backgroundColor: 'rgba(255,255,255,0.95)' }}>
            <Typography
              variant="h4"
              sx={{ color: '#764ba2', textAlign: 'center', mb: 4, fontWeight: 'bold' }}
            >
              How It Works
            </Typography>

            <Box
              sx={{
                display: 'grid',
                gap: 4,
                gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))',
              }}
            >
              {howItWorks.map((step, index) => (
                <Card
                  key={index}
                  sx={{
                    p: 4,
                    borderRadius: 3,
                    textAlign: 'center',
                    position: 'relative',
                    backgroundColor: '#f8f9fa',
                  }}
                >
                  <Box
                    sx={{
                      position: 'absolute',
                      top: '-15px',
                      left: '20px',
                      width: 50,
                      height: 50,
                      borderRadius: '50%',
                      background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      color: '#fff',
                      fontWeight: 'bold',
                      boxShadow: '0 4px 15px rgba(102, 126, 234,0.4)',
                    }}
                  >
                    {step.step}
                  </Box>

                  <img
                    src={step.image}
                    alt={step.title}
                    style={{
                      width: 80,
                      height: 80,
                      margin: '0 auto 1.5rem',
                      display: 'block',
                    }}
                  />

                  <Typography sx={{ fontWeight: 'bold', mb: 1 }}>{step.title}</Typography>
                  <Typography sx={{ fontSize: '0.95rem', color: '#555' }}>{step.description}</Typography>
                </Card>
              ))}
            </Box>
          </Card>
        </Box>

        {/* Footer */}
        <Box sx={{ backgroundColor: 'rgba(0,0,0,0.2)', py: 6 }}>
          <Typography sx={{ color: '#fff', textAlign: 'center', mb: 3 }}>
            LigBind - Making drug discovery simple and accessible
          </Typography>
          <Box
            sx={{
              display: 'flex',
              justifyContent: 'center',
              gap: 4,
              flexWrap: 'wrap',
            }}
          >
            {['About Us', 'How It Works', 'Privacy', 'Contact', 'Help'].map((link, idx) => (
              <Typography
                key={idx}
                sx={{
                  color: '#fff',
                  cursor: 'pointer',
                  opacity: 0.9,
                  '&:hover': { opacity: 1 },
                }}
              >
                {link}
              </Typography>
            ))}
          </Box>
        </Box>

      </Box>
    </Box>
  );
};

export default LigBindLanding;