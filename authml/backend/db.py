from pymongo import MongoClient
import os
from dotenv import load_dotenv

load_dotenv()

# MongoDB connection
MONGO_URI = os.getenv('MONGO_URI', 'mongodb://127.0.0.1:27017/')
client = MongoClient(MONGO_URI)

# Database
db = client['ligbind_db']

# Existing users collection
users_collection = db['users']
users_collection.create_index('email', unique=True)

# NEW: Chemical predictions collection (added safely)
chem_results_collection = db['chem_results']
