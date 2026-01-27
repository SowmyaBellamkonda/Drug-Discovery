import React from 'react';
import { TextField, InputAdornment } from '@mui/material';

const InputFieldWithIcon = ({ 
  icon, 
  label, 
  type = 'text', 
  value, 
  onChange, 
  name,
  error,
  helperText,
  ...props 
}) => {
  return (
    <TextField
      fullWidth
      label={label}
      type={type}
      name={name}
      value={value}
      onChange={onChange}
      error={error}
      helperText={helperText}
      InputProps={{
        startAdornment: (
          <InputAdornment position="start">
            <img 
              src={icon} 
              alt={label} 
              style={{ width: 20, height: 20, opacity: 0.6 }} 
            />
          </InputAdornment>
        ),
      }}
      sx={{ mb: 2 }}
      {...props}
    />
  );
};

export default InputFieldWithIcon;