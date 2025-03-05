from flask import Flask, request, jsonify
from flask_cors import CORS
import numpy as np
import os
import stripe
from dotenv import load_dotenv
from flask_sqlalchemy import SQLAlchemy


app = Flask(__name__)

# Configure SQLite database
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# Initialize SQLAlchemy
db = SQLAlchemy(app)

class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(120), unique=True, nullable=False)
    customer_id = db.Column(db.String(120), nullable=False)
    subscription_id = db.Column(db.String(120), nullable=False)
    subscription_status = db.Column(db.String(50), nullable=False)

    def __repr__(self):
        return f'<User {self.email}>'

class Optimization(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    optimization_type = db.Column(db.String(50), nullable=False)  # 'classical' or 'quantum'
    parameters = db.Column(db.Text, nullable=False)  # JSON string of parameters
    result = db.Column(db.Text, nullable=True)  # JSON string of the result
    created_at = db.Column(db.DateTime, default=db.func.current_timestamp())

    def __repr__(self):
        return f'<Optimization {self.id} for User {self.user_id}>'

load_dotenv()
STRIPE_SECRET = os.getenv("STRIPE_SECRET")

stripe.api_key = STRIPE_SECRET
CORS(app, origins=["http://localhost:3000"], methods=["GET", "POST", "OPTIONS"], allow_headers=["Content-Type"])


# Include the existing routes for Stripe integration
@app.route('/subscribe', methods=['POST'])
def subscribe_user():
    try:
        data = request.get_json()
        
        if not data or "email" not in data or "paymentMethodId" not in data:
            return jsonify({"error": "Invalid input: Missing 'email' or 'paymentMethodId'"}), 400

        email = data["email"]
        payment_method_id = data["paymentMethodId"]
        price_id = "price_1QYNn9JQZaUHxA2Ld9rV2MPd"

        # Check for existing customers
        customers = stripe.Customer.list(email=email).data
        
        if customers:
            customer = customers[0]  # Use the existing customer
        else:
            customer = stripe.Customer.create(email=email)

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

        # Save user and subscription data in the database
        existing_user = User.query.filter_by(email=email).first()
        if existing_user:
            existing_user.customer_id = customer.id
            existing_user.subscription_id = subscription.id
            existing_user.subscription_status = "active"
        else:
            new_user = User(
                email=email,
                customer_id=customer.id,
                subscription_id=subscription.id,
                subscription_status="active"
            )
            db.session.add(new_user)
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

@app.route('/check-subscription', methods=['POST'])
def check_subscription():
    data = request.get_json()
    email = data.get("email")

    if not email:
        return jsonify({"error": "Email is required"}), 400

    user = User.query.filter_by(email=email).first()
    if user and user.subscription_status == "active":
        return jsonify({"isSubscribed": True}), 200
    return jsonify({"isSubscribed": False}), 200

@app.route('/cancel-subscription', methods=['POST'])
def cancel_subscription():
    try:
        data = request.get_json()
        email = data.get("email")

        if not email:
            return jsonify({"error": "Email is required"}), 400

        user = User.query.filter_by(email=email).first()
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

@app.route('/webhook', methods=['POST'])
def stripe_webhook():
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

@app.route('/test', methods=['POST'])
def optimize_test():
    return "POST test request received successfully!"

@app.route('/')
def index():
    return "Optimol API"

if __name__ == '__main__':
    with app.app_context():
        db.create_all()  # Create tables within the app context
    app.run(host="0.0.0.0", port=5000)