from flask import Flask, jsonify
from flask_cors import CORS
import os
import json
from dotenv import load_dotenv
from extensions import db

# Function to load configuration from JSON file
def load_config(config_file='config.json'):
    try:
        with open(config_file, 'r') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading config: {e}")
        return {}

# Initialize Flask app
app = Flask(__name__)

# Load configuration
config = load_config()

# Configure SQLAlchemy with settings from config file
app.config['SQLALCHEMY_DATABASE_URI'] = config.get('SQLALCHEMY_DATABASE_URI', 'sqlite:///users.db')
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = config.get('SQLALCHEMY_TRACK_MODIFICATIONS', False)

# Initialize extensions
db.init_app(app)
load_dotenv()

# Configure CORS with settings from config
cors_config = config.get('CORS', {})
CORS(app, 
     origins=cors_config.get('origins'), 
     methods=cors_config.get('methods'), 
     allow_headers=cors_config.get('allow_headers'))

# Import blueprints after app creation
from user import user_bp
from opti import opti_bp

# Register blueprints
app.register_blueprint(user_bp)
app.register_blueprint(opti_bp)

@app.route('/')
def index():
    return "Optimol API"

@app.route('/test', methods=['POST'])
def optimize_test():
    return "POST test request received successfully!"

if __name__ == '__main__':
    with app.app_context():
        db.create_all()  # Create tables within the app context
    app.run(host="0.0.0.0", port=5000)