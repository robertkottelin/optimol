from flask import Flask, jsonify, request
from sqlalchemy import text
from flask_cors import CORS
import os
import json
from dotenv import load_dotenv
from extensions import db
from flask_jwt_extended import JWTManager

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

# JWT Configuration
app.config['JWT_SECRET_KEY'] = os.environ.get('JWT_SECRET_KEY', 'dev-key-change-in-production')
app.config['JWT_TOKEN_LOCATION'] = ['cookies']
app.config['JWT_COOKIE_SECURE'] = True
app.config['JWT_COOKIE_CSRF_PROTECT'] = False
app.config['JWT_COOKIE_SAMESITE'] = 'None'
app.config['JWT_ACCESS_TOKEN_EXPIRES'] = 30 * 24 * 60 * 60  # 30 days
app.config['JWT_COOKIE_DOMAIN'] = 'optimizemolecule.com'
app.config['JWT_COOKIE_PATH'] = '/'
app.config['JWT_ACCESS_COOKIE_NAME'] = 'access_token_cookie'
app.config['JWT_REFRESH_COOKIE_NAME'] = 'refresh_token_cookie'
jwt = JWTManager(app)

# Initialize extensions
db.init_app(app)
load_dotenv()

# CORS Configuration
CORS_CONFIG = {
    'origins': config.get('CORS', {}).get('origins', ["https://robertkottelin.github.io", "https://optimizemolecule.com"]),
    'methods': config.get('CORS', {}).get('methods', ["GET", "POST", "OPTIONS"]),
    'allow_headers': config.get('CORS', {}).get('allow_headers', ["Content-Type", "Authorization", "X-Requested-With"]),
    'expose_headers': ['Content-Type', 'Authorization'],
    'max_age': 600
}

# Initialize CORS with resource pattern matching
CORS(app, 
     resources={r"/*": {
         "origins": CORS_CONFIG['origins'],
         "methods": CORS_CONFIG['methods'], 
         "allow_headers": CORS_CONFIG['allow_headers'],
         "supports_credentials": True,
         "expose_headers": CORS_CONFIG['expose_headers']
     }})

# Global after_request handler for CORS
@app.after_request
def after_request(response):
    origin = request.headers.get('Origin')
    if origin and origin in CORS_CONFIG['origins']:
        response.headers.set('Access-Control-Allow-Origin', origin)
        response.headers.set('Access-Control-Allow-Headers', ', '.join(CORS_CONFIG['allow_headers']))
        response.headers.set('Access-Control-Allow-Methods', ', '.join(CORS_CONFIG['methods']))
        response.headers.set('Access-Control-Allow-Credentials', 'true')
        response.headers.set('Access-Control-Max-Age', str(CORS_CONFIG['max_age']))
    return response

# Route to handle specific preflight OPTIONS requests
@app.route('/me', methods=['OPTIONS'])
@app.route('/login', methods=['OPTIONS'])
@app.route('/register', methods=['OPTIONS'])
@app.route('/logout', methods=['OPTIONS']) 
@app.route('/subscribe', methods=['OPTIONS'])
@app.route('/check-subscription', methods=['OPTIONS'])
@app.route('/cancel-subscription', methods=['OPTIONS'])
@app.route('/optimize-molecule', methods=['OPTIONS'])
def handle_cors_preflight():
    response = jsonify({})
    origin = request.headers.get('Origin')
    if origin and origin in CORS_CONFIG['origins']:
        response.headers.set('Access-Control-Allow-Origin', origin)
        response.headers.set('Access-Control-Allow-Headers', ', '.join(CORS_CONFIG['allow_headers']))
        response.headers.set('Access-Control-Allow-Methods', ', '.join(CORS_CONFIG['methods']))
        response.headers.set('Access-Control-Allow-Credentials', 'true')
        response.headers.set('Access-Control-Max-Age', str(CORS_CONFIG['max_age']))
    return response, 200

# Global OPTIONS handler for all routes
@app.route('/<path:path>', methods=['OPTIONS'])
def handle_all_options(path):
    response = jsonify({})
    origin = request.headers.get('Origin')
    if origin and origin in CORS_CONFIG['origins']:
        response.headers.set('Access-Control-Allow-Origin', origin)
        response.headers.set('Access-Control-Allow-Headers', ', '.join(CORS_CONFIG['allow_headers']))
        response.headers.set('Access-Control-Allow-Methods', ', '.join(CORS_CONFIG['methods']))
        response.headers.set('Access-Control-Allow-Credentials', 'true')
        response.headers.set('Access-Control-Max-Age', str(CORS_CONFIG['max_age']))
    return response, 200

# Test endpoint for cookie verification
@app.route('/test-cookie')
def test_cookie():
    resp = jsonify({"success": True})
    resp.set_cookie(
        'test_cookie', 
        value='test_value',
        secure=True,
        httponly=False,
        samesite='None',
        path='/',
        max_age=3600
    )
    return resp

# Import blueprints after app creation
from user import user_bp
from opti import opti_bp

# Register blueprints
app.register_blueprint(user_bp)
app.register_blueprint(opti_bp)

# Initialize database tables during application startup
with app.app_context():
    db.create_all()
    print("Database tables initialized")

@app.route('/')
def index():
    return "Optimol API"

@app.route('/test', methods=['POST'])
def optimize_test():
    return "POST test request received successfully!"

@app.route('/health')
def health_check():
    try:
        db.session.execute(text('SELECT 1'))
        return jsonify({"status": "healthy"}), 200
    except Exception as e:
        return jsonify({"status": "unhealthy", "error": str(e)}), 500

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)