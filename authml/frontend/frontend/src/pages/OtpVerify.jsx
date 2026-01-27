import React, { useState, useEffect } from 'react';
import { Box, Button, Card, Typography, TextField, Alert } from '@mui/material';
import { useNavigate, useLocation } from 'react-router-dom';
import api from '../api/axiosConfig';

const OtpVerify = () => {
  const navigate = useNavigate();
  const location = useLocation();
  const email = location.state?.email;

  const [otp, setOtp] = useState('');
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [resendLoading, setResendLoading] = useState(false);
  const [timeLeft, setTimeLeft] = useState(300); // 5 minutes in seconds

  useEffect(() => {
    if (!email) {
      navigate('/signup');
      return;
    }

    const timer = setInterval(() => {
      setTimeLeft((prev) => {
        if (prev <= 1) {
          clearInterval(timer);
          return 0;
        }
        return prev - 1;
      });
    }, 1000);

    return () => clearInterval(timer);
  }, [email, navigate]);

  const formatTime = (seconds) => {
    const mins = Math.floor(seconds / 60);
    const secs = seconds % 60;
    return `${mins}:${secs.toString().padStart(2, '0')}`;
  };

  const handleVerify = async (e) => {
    e.preventDefault();
    setError('');
    setLoading(true);

    try {
      await api.post('/verify-otp', { email, otp });
      navigate('/login', {
        state: { message: 'Email verified successfully! Please login.' },
      });
    } catch (err) {
      setError(err.response?.data?.error || 'Verification failed. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  const handleResend = async () => {
    setError('');
    setResendLoading(true);

    try {
      await api.post('/resend-otp', { email });
      setTimeLeft(300);
    } catch (err) {
      setError(err.response?.data?.error || 'Failed to resend OTP.');
    } finally {
      setResendLoading(false);
    }
  };

  return (
    <Box
      sx={{
        minHeight: '100vh',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
        padding: 2,
      }}
    >
      <Card
        sx={{
          maxWidth: 450,
          width: '100%',
          padding: 4,
          borderRadius: 2,
          boxShadow: '0 8px 32px rgba(0,0,0,0.1)',
        }}
      >
        <Typography
          variant="h4"
          component="h1"
          textAlign="center"
          fontWeight="bold"
          color="primary"
          mb={1}
        >
          Verify Your Email
        </Typography>
        <Typography variant="body2" textAlign="center" color="text.secondary" mb={3}>
          Enter the 6-digit OTP sent to <br />
          <strong>{email}</strong>
        </Typography>

        {error && <Alert severity="error" sx={{ mb: 2 }}>{error}</Alert>}

        <form onSubmit={handleVerify}>
          <TextField
            fullWidth
            label="Enter OTP"
            value={otp}
            onChange={(e) => {
              setOtp(e.target.value.replace(/\D/g, '').slice(0, 6));
              setError('');
            }}
            inputProps={{
              maxLength: 6,
              style: { textAlign: 'center', fontSize: '24px', letterSpacing: '8px' },
            }}
            sx={{ mb: 2 }}
            required
          />

          <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
            <Typography variant="body2" color={timeLeft === 0 ? 'error' : 'text.secondary'}>
              {timeLeft > 0 ? `Time left: ${formatTime(timeLeft)}` : 'OTP Expired'}
            </Typography>

            <Button
              variant="text"
              size="small"
              onClick={handleResend}
              disabled={resendLoading || timeLeft > 0}
              sx={{ textTransform: 'none' }}
            >
              {resendLoading ? 'Sending...' : 'Resend OTP'}
            </Button>
          </Box>

          <Button
            type="submit"
            variant="contained"
            fullWidth
            size="large"
            disabled={loading || otp.length !== 6}
            sx={{
              py: 1.5,
              textTransform: 'none',
              fontSize: '16px',
              fontWeight: 600,
            }}
          >
            {loading ? 'Verifying...' : 'Verify OTP'}
          </Button>
        </form>

        <Typography variant="body2" textAlign="center" mt={3}>
          Wrong email?{' '}
          <Typography
            component="span"
            color="primary"
            sx={{ cursor: 'pointer', fontWeight: 600 }}
            onClick={() => navigate('/signup')}
          >
            Go Back
          </Typography>
        </Typography>
      </Card>
    </Box>
  );
};

export default OtpVerify;
