from flask import Blueprint, request, jsonify
import os
from datetime import datetime
from extensions import db
from models import User
from werkzeug.security import generate_password_hash, check_password_hash
from flask_jwt_extended import (
    create_access_token, set_access_cookies,
    unset_jwt_cookies, jwt_required, get_jwt_identity
)

# Create blueprint
user_bp = Blueprint('user', __name__)

@user_bp.route('/register', methods=['POST'])
def register():
    data = request.get_json()
    email = data.get('email')
    password = data.get('password')

    if not email or not password:
        return jsonify({"error": "Email and password required"}), 400

    existing_user = User.query.filter_by(email=email).first()
    if existing_user:
        return jsonify({"error": "Email already registered"}), 409

    # Create user
    user = User(email=email)
    user.password = generate_password_hash(password)
    db.session.add(user)
    db.session.commit()

    # Create token without setting cookie
    access_token = create_access_token(identity=str(user.id))

    # Return token in response body
    return jsonify({
        "success": True,
        "token": access_token,
        "user": {
            "email": user.email
        }
    }), 201

@user_bp.route('/login', methods=['POST'])
def login():
    data = request.get_json()
    email = data.get('email')
    password = data.get('password')

    user = User.query.filter_by(email=email).first()
    if not user or not check_password_hash(user.password, password):
        return jsonify({"error": "Invalid email or password"}), 401

    # Create token without setting cookie
    access_token = create_access_token(identity=str(user.id))

    # Return token in response body
    return jsonify({
        "success": True,
        "token": access_token,
        "user": {
            "email": user.email
        }
    })

@user_bp.route('/logout', methods=['POST'])
def logout():
    # Simple success response - token handling done client-side
    return jsonify({"success": True})

@user_bp.route('/me', methods=['GET'])
@jwt_required()
def get_current_user():
    # get_jwt_identity() returns the identity from the JWT which is now a string
    user_id = get_jwt_identity()
    # Convert string user_id back to integer for database query
    user = User.query.get(int(user_id))

    if not user:
        return jsonify({"error": "User not found"}), 404

    return jsonify({
        "email": user.email
    })

