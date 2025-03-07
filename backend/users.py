from flask import Blueprint, request, jsonify
import os
import stripe
from datetime import datetime
from app import db

# Create blueprint
user_bp = Blueprint('user', __name__)

class User(db.Model):
    """User model for subscription management."""
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(120), unique=True, nullable=False)
    customer_id = db.Column(db.String(120), nullable=False)
    subscription_id = db.Column(db.String(120), nullable=False)
    subscription_status = db.Column(db.String(50), nullable=False)

    def __repr__(self):
        return f'<User {self.email}>'

# Initialize Stripe
stripe.api_key = os.getenv("STRIPE_SECRET")

@user_bp.route('/subscribe', methods=['POST'])
def subscribe_user():
    """Create or update a subscription for a user."""
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

@user_bp.route('/check-subscription', methods=['POST'])
def check_subscription():
    """Check if a user has an active subscription."""
    data = request.get_json()
    email = data.get("email")

    if not email:
        return jsonify({"error": "Email is required"}), 400

    user = User.query.filter_by(email=email).first()
    if user and user.subscription_status == "active":
        return jsonify({"isSubscribed": True}), 200
    return jsonify({"isSubscribed": False}), 200

@user_bp.route('/cancel-subscription', methods=['POST'])
def cancel_subscription():
    """Cancel a user's subscription."""
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