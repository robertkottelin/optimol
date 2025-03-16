from flask import Blueprint, request, jsonify
import os
import stripe
from datetime import datetime
from extensions import db
from werkzeug.security import generate_password_hash, check_password_hash
from flask_jwt_extended import (
    create_access_token, set_access_cookies, 
    unset_jwt_cookies, jwt_required, get_jwt_identity
)

# Create blueprint
user_bp = Blueprint('user', __name__)

class User(db.Model):
    """User model for subscription management."""
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(120), unique=True, nullable=False)
    password_hash = db.Column(db.String(256), nullable=True)
    customer_id = db.Column(db.String(120), nullable=False)
    subscription_id = db.Column(db.String(120), nullable=False)
    subscription_status = db.Column(db.String(50), nullable=False)
    
    def set_password(self, password):
        self.password_hash = generate_password_hash(password)
        
    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    def __repr__(self):
        return f'<User {self.email}>'

# Initialize Stripe
stripe.api_key = os.getenv("STRIPE_SECRET")

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
    
    # Create unsubscribed user
    user = User(
        email=email,
        customer_id="unsubscribed",
        subscription_id="none",
        subscription_status="inactive"
    )
    user.set_password(password)
    db.session.add(user)
    db.session.commit()
    
    # Create token without setting cookie
    # FIXED: Convert user.id to string for JWT subject
    access_token = create_access_token(identity=str(user.id))
    
    # Return token in response body
    return jsonify({
        "success": True, 
        "token": access_token,
        "user": {
            "email": user.email, 
            "isSubscribed": False
        }
    }), 201

@user_bp.route('/login', methods=['POST'])
def login():
    data = request.get_json()
    email = data.get('email')
    password = data.get('password')
    
    user = User.query.filter_by(email=email).first()
    if not user or not user.check_password(password):
        return jsonify({"error": "Invalid email or password"}), 401
    
    # Create token without setting cookie
    # FIXED: Convert user.id to string for JWT subject
    access_token = create_access_token(identity=str(user.id))
    
    # Return token in response body
    return jsonify({
        "success": True,
        "token": access_token,
        "user": {
            "email": user.email, 
            "isSubscribed": user.subscription_status == "active"
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
        "email": user.email,
        "isSubscribed": user.subscription_status == "active"
    })

@user_bp.route('/register-and-subscribe', methods=['POST'])
def register_and_subscribe():
    """Combined endpoint to register a new user and create subscription in one step."""
    try:
        data = request.get_json()
        email = data.get('email')
        password = data.get('password')
        payment_method_id = data.get('paymentMethodId')
        
        if not email or not password or not payment_method_id:
            return jsonify({"error": "Email, password, and payment method required"}), 400
            
        existing_user = User.query.filter_by(email=email).first()
        if existing_user:
            return jsonify({"error": "Email already registered"}), 409
        
        # Create Stripe customer
        customer = stripe.Customer.create(email=email)
        
        # Attach payment method to customer
        stripe.PaymentMethod.attach(payment_method_id, customer=customer.id)
        
        # Set as default payment method
        stripe.Customer.modify(
            customer.id,
            invoice_settings={"default_payment_method": payment_method_id}
        )
        
        # Create subscription
        price_id = "price_1QYNn9JQZaUHxA2Ld9rV2MPd"
        subscription = stripe.Subscription.create(
            customer=customer.id,
            items=[{"price": price_id}],
            expand=["latest_invoice.payment_intent"],
        )
        
        # Create user with subscription already active
        user = User(
            email=email,
            customer_id=customer.id,
            subscription_id=subscription.id,
            subscription_status="active"
        )
        user.set_password(password)
        
        db.session.add(user)
        db.session.commit()
        
        # Create token
        access_token = create_access_token(identity=str(user.id))
        
        return jsonify({
            "success": True,
            "token": access_token,
            "subscriptionId": subscription.id,
            "clientSecret": subscription.latest_invoice.payment_intent.client_secret,
            "user": {
                "email": user.email,
                "isSubscribed": True
            }
        }), 201
        
    except stripe.error.StripeError as e:
        return jsonify({"error": str(e)}), 500
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@user_bp.route('/subscribe', methods=['POST'])
@jwt_required()
def subscribe_user():
    try:
        # Convert string user_id back to integer for database query
        user_id = get_jwt_identity()
        user = User.query.get(int(user_id))
        
        if not user:
            return jsonify({"error": "User not found"}), 404
        
        data = request.get_json()
        payment_method_id = data.get("paymentMethodId")
        
        price_id = "price_1QYNn9JQZaUHxA2Ld9rV2MPd"

        # Check for existing customers using authenticated user's email
        customers = stripe.Customer.list(email=user.email).data
        
        if customers:
            customer = customers[0]  # Use the existing customer
        else:
            customer = stripe.Customer.create(email=user.email)

        # Attach payment method
        stripe.PaymentMethod.attach(payment_method_id, customer=customer.id)

        stripe.Customer.modify(
            customer.id,
            invoice_settings={"default_payment_method": payment_method_id}
        )

        # Create subscription
        subscription = stripe.Subscription.create(
            customer=customer.id,
            items=[{"price": price_id}],
            expand=["latest_invoice.payment_intent"],
        )

        # Update authenticated user with subscription data
        user.customer_id = customer.id
        user.subscription_id = subscription.id
        user.subscription_status = "active"
        db.session.commit()

        return jsonify({
            "success": True,
            "subscriptionId": subscription.id,
            "clientSecret": subscription.latest_invoice.payment_intent.client_secret,
        })
    except stripe.error.StripeError as e:
        return jsonify({"error": str(e)}), 500
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@user_bp.route('/check-subscription', methods=['GET'])
@jwt_required()
def check_subscription():
    """Check if a user has an active subscription."""
    # Convert string user_id back to integer for database query
    user_id = get_jwt_identity()
    user = User.query.get(int(user_id))
    
    if not user:
        return jsonify({"error": "User not found"}), 404
        
    return jsonify({"isSubscribed": user.subscription_status == "active"}), 200

@user_bp.route('/cancel-subscription', methods=['POST'])
@jwt_required()
def cancel_subscription():
    """Cancel a user's subscription."""
    try:
        # Convert string user_id back to integer for database query
        user_id = get_jwt_identity()
        user = User.query.get(int(user_id))
        
        if not user:
            return jsonify({"error": "User not found"}), 404

        # Cancel the subscription in Stripe
        stripe.Subscription.delete(user.subscription_id)

        # Update the user's subscription status in the database
        user.subscription_status = "canceled"
        db.session.commit()

        return jsonify({"success": True, "message": "Subscription canceled successfully."}), 200
    except stripe.error.StripeError as e:
        return jsonify({"error": str(e)}), 500
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@user_bp.route('/webhook', methods=['POST'])
def stripe_webhook():
    """Handle Stripe webhook events for subscription lifecycle management."""
    payload = request.get_data(as_text=True)
    sig_header = request.headers.get('Stripe-Signature')

    endpoint_secret = os.getenv("STRIPE_WEBHOOK_SECRET")
    event = None

    try:
        event = stripe.Webhook.construct_event(
            payload, sig_header, endpoint_secret
        )
    except ValueError:
        return jsonify({"error": "Invalid payload"}), 400
    except stripe.error.SignatureVerificationError:
        return jsonify({"error": "Invalid signature"}), 400

    # Handle the event
    try:
        if event['type'] == 'invoice.payment_succeeded':
            invoice = event['data']['object']
            customer_id = invoice['customer']

            # Update subscription status to active
            user = User.query.filter_by(customer_id=customer_id).first()
            if user:
                user.subscription_status = "active"
                db.session.commit()

        elif event['type'] == 'customer.subscription.deleted':
            subscription = event['data']['object']
            customer_id = subscription['customer']

            # Update subscription status to canceled
            user = User.query.filter_by(customer_id=customer_id).first()
            if user:
                user.subscription_status = "canceled"
                db.session.commit()

        elif event['type'] == 'invoice.payment_failed':
            invoice = event['data']['object']
            customer_id = invoice['customer']

            # Update subscription status to payment_failed
            user = User.query.filter_by(customer_id=customer_id).first()
            if user:
                user.subscription_status = "payment_failed"
                db.session.commit()

        elif event['type'] == 'customer.subscription.updated':
            subscription = event['data']['object']
            customer_id = subscription['customer']

            # Update subscription status based on the new status
            user = User.query.filter_by(customer_id=customer_id).first()
            if user:
                user.subscription_status = subscription['status']
                db.session.commit()
                
        return jsonify({"status": "success"}), 200

    except Exception as e:
        return jsonify({"error": f"Webhook handling error: {str(e)}"}), 500