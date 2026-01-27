import random
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import os
from dotenv import load_dotenv

load_dotenv()

def generate_otp():
    """Generate a unique 6-digit OTP"""
    return random.randint(100000, 999999)

def send_otp_email(recipient_email, otp):
    """Send OTP via email using SMTP"""
    try:
        sender_email = os.getenv('EMAIL')
        sender_password = os.getenv('PASSWORD')
        
        if not sender_email or not sender_password:
            raise Exception("Email credentials not configured")
        
        # Create message
        message = MIMEMultipart('alternative')
        message['Subject'] = 'LigBind - Email Verification OTP'
        message['From'] = sender_email
        message['To'] = recipient_email
        
        # HTML email body
        html = f"""
        <html>
          <body>
            <div style="font-family: Arial, sans-serif; max-width: 600px; margin: 0 auto;">
              <h2 style="color: #1976d2;">LigBind Email Verification</h2>
              <p>Thank you for signing up with LigBind!</p>
              <p>Your One-Time Password (OTP) for email verification is:</p>
              <div style="background-color: #f5f5f5; padding: 20px; text-align: center; font-size: 32px; font-weight: bold; letter-spacing: 5px; color: #1976d2; margin: 20px 0;">
                {otp}
              </div>
              <p>This OTP will expire in <strong>5 minutes</strong>.</p>
              <p>If you didn't request this verification, please ignore this email.</p>
              <hr style="margin: 30px 0; border: none; border-top: 1px solid #ddd;">
              <p style="color: #666; font-size: 12px;">This is an automated email. Please do not reply.</p>
            </div>
          </body>
        </html>
        """
        
        # Attach HTML content
        html_part = MIMEText(html, 'html')
        message.attach(html_part)
        
        # Send email via Gmail SMTP
        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as server:
            server.login(sender_email, sender_password)
            server.sendmail(sender_email, recipient_email, message.as_string())
        
        print(f"OTP sent successfully to {recipient_email}")
        return True
        
    except Exception as e:
        print(f"Failed to send OTP email: {str(e)}")
        raise Exception(f"Failed to send OTP email: {str(e)}")