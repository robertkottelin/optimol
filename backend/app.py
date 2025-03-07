from flask import Flask, jsonify
from flask_cors import CORS
import os
from dotenv import load_dotenv
from flask_sqlalchemy import SQLAlchemy

# Initialize Flask app
app = Flask(__name__)

# Configure SQLite database
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# Initialize SQLAlchemy
db = SQLAlchemy(app)

# Initialize extensions
load_dotenv()
CORS(app, origins=["http://localhost:3000"], methods=["GET", "POST", "OPTIONS"], allow_headers=["Content-Type"])

# Import blueprints after db initialization to avoid circular imports
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