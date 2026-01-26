LigBind Authentication System
A full-stack authentication system with OTP email verification built with React, Flask, and MongoDB.

Features
✅ User Signup with email verification
✅ 6-digit OTP sent via email
✅ OTP expires after 5 minutes
✅ Secure password hashing
✅ Login with verification check
✅ Material-UI components
✅ Responsive design
✅ Error handling
Project Structure
LIGBIND_PROJECT/
├── backend/
│   ├── app.py                 # Flask API routes
│   ├── db.py                  # MongoDB connection
│   ├── otp_helper.py          # OTP generation & email
│   ├── requirements.txt       # Python dependencies
│   └── .env                   # Environment variables
│
└── frontend/
    ├── index.html
    ├── package.json
    ├── vite.config.js
    │
    ├── public/
    │   ├── email.png          # Add your icon
    │   ├── user.png           # Add your icon
    │   └── lock.png           # Add your icon
    │
    └── src/
        ├── main.jsx
        ├── App.jsx
        │
        ├── api/
        │   └── axiosConfig.js
        │
        ├── pages/
        │   ├── Signup.jsx
        │   ├── OtpVerify.jsx
        │   ├── Login.jsx
        │   └── LigbindHome.jsx
        │
        └── components/
            └── InputFieldWithIcon.jsx
Setup Instructions
Backend Setup
Navigate to backend directory:
bash
   cd backend
Create virtual environment:
bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
Install dependencies:
bash
   pip install -r requirements.txt
Configure environment variables: Create .env file:
env
   EMAIL=your-email@gmail.com
   PASSWORD=your-app-password
   MONGO_URI=mongodb://localhost:27017/
Important for Gmail:

Enable 2-Factor Authentication
Go to Google Account → Security → App Passwords
Generate app password for "Mail"
Use that password in .env
Start MongoDB:
bash
   # If using local MongoDB
   mongod
   
   # Or use MongoDB Atlas cloud database
Run Flask server:
bash
   python app.py
Server runs on http://localhost:5000

Frontend Setup
Navigate to frontend directory:
bash
   cd frontend
Install dependencies:
bash
   npm install
Add PNG icons: Place these icons in public/ folder:
email.png
user.png
lock.png
Start development server:
bash
   npm run dev
App runs on http://localhost:3000

User Flow
1. Signup
User enters name, email, password
Backend generates 6-digit OTP
OTP sent to user's email
OTP expires in 5 minutes
Redirect to OTP verification page
2. OTP Verification
User enters 6-digit OTP
Backend validates OTP
If correct: isVerified = true, redirect to Login
If incorrect: Show error message
Option to resend OTP
3. Login
User enters email and password
Backend checks:
User exists?
Email verified?
Password correct?
If all checks pass: Redirect to Home page
Store user data in localStorage
4. Home Page
Display user information
Protected route (requires authentication)
Logout functionality
API Endpoints
POST /api/signup
Request:

json
{
  "name": "John Doe",
  "email": "john@example.com",
  "password": "password123"
}
Response:

json
{
  "message": "Signup successful. OTP sent to your email.",
  "email": "john@example.com"
}
POST /api/verify-otp
Request:

json
{
  "email": "john@example.com",
  "otp": "123456"
}
Response:

json
{
  "message": "Email verified successfully"
}
POST /api/login
Request:

json
{
  "email": "john@example.com",
  "password": "password123"
}
Response:

json
{
  "message": "Login successful",
  "user": {
    "name": "John Doe",
    "email": "john@example.com"
  }
}
POST /api/resend-otp
Request:

json
{
  "email": "john@example.com"
}
Response:

json
{
  "message": "OTP resent successfully"
}
Database Schema
Users Collection
javascript
{
  _id: ObjectId,
  name: String,
  email: String (unique),
  password: String (hashed),
  otp: Number (6 digits),
  otp_expiry: Date,
  isVerified: Boolean,
  created_at: Date
}

Security Features
✅ Password hashing with Werkzeug
✅ Unique email validation
✅ OTP expiration (5 minutes)
✅ Email verification required
✅ CORS enabled
✅ Input validation

Technologies Used
Frontend
React 18
Material-UI (MUI) v5
React Router v6
Axios
Vite
Backend
Python 3.x
Flask
Flask-CORS
PyMongo
SMTPLib
Werkzeug
Database
MongoDB
Common Issues & Solutions
Email Not Sending
Check Gmail app password is correct
Ensure 2FA is enabled
Check spam/junk folder
Verify email credentials in .env
MongoDB Connection Failed
Ensure MongoDB is running
Check MONGO_URI in .env
For Atlas: whitelist your IP
CORS Errors
Ensure Flask-CORS is installed
Backend should run on port 5000
Frontend should run on port 3000
OTP Not Expiring
Check system time is correct
Verify datetime is UTC
Development
To run in development mode:

Terminal 1 (Backend):

bash
cd backend
python app.py
Terminal 2 (Frontend):

bash
cd frontend
npm run dev
Production Build
Frontend:

bash
cd frontend
npm run build
The build files will be in frontend/dist/

License
MIT

Support
For issues or questions, please create an issue in the repository.

Author
Sowmya Bellamkonda
GitHub: SowmyaBellamkonda
