from flask import Flask, request, jsonify
from flask_cors import CORS
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from qiskit_aer import Aer
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit_algorithms.optimizers import COBYLA, SPSA
from qiskit_algorithms.minimum_eigensolvers import QAOA
from qiskit.primitives import Estimator, StatevectorSampler
import numpy as np
import os
import stripe
from dotenv import load_dotenv
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)

# Configure SQLite database
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'  # SQLite file will be `users.db`
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False  # Suppress warnings

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

# CORS(app, origins=["https://orca-app-pmrz6.ondigitalocean.app"], methods=["GET", "POST", "OPTIONS"], allow_headers=["Content-Type"])
# same but for localhost
# CORS(app, origins=["*"], methods=["GET", "POST", "OPTIONS"], allow_headers=["Content-Type"])
CORS(app, origins=["http://localhost:3000"], methods=["GET", "POST", "OPTIONS"], allow_headers=["Content-Type"])


def optimize_molecule(data):
    """
    Optimize the molecular geometry based on classical energy calculations.
    """
    atoms = []
    for atom in data["atoms"]:
        atoms.append((atom["element"], (atom["x"], atom["y"], atom["z"])))

    molecule = Atoms(symbols=[a[0] for a in atoms], positions=[a[1] for a in atoms])
    molecule.calc = EMT()

    # Validate initial geometry
    if any(np.linalg.norm(atom1 - atom2) < 0.5 for atom1, atom2 in zip(molecule.positions[:-1], molecule.positions[1:])):
        raise ValueError("Initial geometry contains overlapping atoms or invalid distances.")

    # Retrieve user-specified parameters or use defaults
    fmax = data.get("fmax", 0.005)
    steps = data.get("steps", 500)

    optimizer = BFGS(molecule)
    optimizer.run(fmax=fmax, steps=steps)

    # Validate positions after optimization
    if not np.all(np.isfinite(molecule.positions)):
        raise ValueError("Optimization failed: atomic positions contain NaN.")

    optimized_atoms = [
        {"id": i + 1, "element": atom.symbol, "x": atom.position[0], "y": atom.position[1], "z": atom.position[2]}
        for i, atom in enumerate(molecule)
    ]

    # Final validation
    if any(np.isnan(coord) for atom in optimized_atoms for coord in (atom["x"], atom["y"], atom["z"])):
        raise ValueError("Optimization result contains NaN in atomic positions.")

    return {"atoms": optimized_atoms}

from qiskit.primitives import Sampler
from qiskit.circuit.library import QAOAAnsatz
from qiskit.circuit import Parameter

from scipy.optimize import minimize
from scipy.optimize import minimize
def optimize_molecule_qaoa(data, optimizer='COBYLA', p=10):
    """
    Optimize a molecular geometry using QAOA and return adjusted positions.
    """
    try:
        maxiter = data.get("maxiter", 1000)  # Retrieve maxiter from request

        # Validate incoming data structure
        if "atoms" not in data:
            raise ValueError("Missing 'atoms' in input data.")
        
        # Debug incoming data
        print("Data received for quantum optimization:", data)

        # Define a diagonal problem Hamiltonian (example: Z + Z)
        problem = SparsePauliOp(Pauli('Z')) + SparsePauliOp(Pauli('Z'))
        print("Diagonal Problem Hamiltonian defined:", problem)

        # Initialize QAOA ansatz
        mixer = SparsePauliOp(Pauli('X'))  # Default mixer is a single Pauli-X term
        ansatz = QAOAAnsatz(cost_operator=problem, reps=p, mixer_operator=mixer)
        print("QAOA ansatz initialized.")

        # Initialize Estimator to compute expectation values
        estimator = Estimator()
        print("Estimator initialized.")

        # Define a parameter vector for optimization
        parameters = ansatz.parameters
        initial_point = np.random.uniform(-np.pi, np.pi, len(parameters))

        def objective_function(params):
            """Objective function for optimization."""
            # Assign parameters to the ansatz
            ansatz_with_params = ansatz.assign_parameters({param: val for param, val in zip(parameters, params)})

            # Compute the expectation value of the problem Hamiltonian
            job = estimator.run([ansatz_with_params], [problem])
            result = job.result()
            return result.values[0]

        # Optimize parameters using scipy.optimize.minimize
        optimization_result = minimize(
            fun=objective_function,
            x0=initial_point,
            method=optimizer.lower(),
            options={"maxiter": maxiter}
        )

        optimal_params = optimization_result.x
        min_energy = optimization_result.fun
        print("Optimization completed.")
        print("Optimal parameters:", optimal_params)
        print("Minimum energy:", min_energy)

        # Adjust atomic positions based on optimized parameters
        scale_factor = np.mean(optimal_params)
        adjusted_atoms = [
            {"id": i + 1, "element": atom["element"],
             "x": atom["x"] * (1 + scale_factor),
             "y": atom["y"] * (1 + scale_factor),
             "z": atom["z"] * (1 + scale_factor)}
            for i, atom in enumerate(data["atoms"])
        ]

        return {"atoms": adjusted_atoms, "optimal_params": optimal_params.tolist(), "min_energy": min_energy}
    except KeyError as ke:
        print("KeyError in optimize_molecule_qaoa:", str(ke))
        raise ValueError(f"Invalid input: Missing key {str(ke)}.")
    except Exception as e:
        print("Error in optimize_molecule_qaoa:", str(e))
        raise

import logging
logging.basicConfig(level=logging.DEBUG)

@app.route('/optimize', methods=['POST'])
def classical_optimize():
    logging.debug("Request received at /optimize")
    logging.debug("Headers: %s", request.headers)
    logging.debug("Data: %s", request.get_json())

    data = request.get_json()
    if not data or "file1" not in data:
        return jsonify({"error": "Invalid input"}), 400

    try:
        optimized_data = optimize_molecule(data["file1"])
        return jsonify({"optimized_file1": optimized_data})
    except Exception as e:
        logging.error("Error in /optimize: %s", str(e))
        return jsonify({"error": str(e)}), 500

@app.route('/quantum-optimize', methods=['POST'])
def quantum_optimize():
    """
    Optimize molecular energy using QAOA.
    """
    data = request.get_json()
    if not data or "file1" not in data:
        return jsonify({"error": "Invalid input: Missing 'file1' in request body"}), 400

    try:
        optimizer = data.get("optimizer", "COBYLA")
        p = data.get("p", 2)
        quantum_result = optimize_molecule_qaoa(data["file1"], optimizer, p)
        return jsonify({"optimized_file1": quantum_result})
    except Exception as e:
        # Log detailed error for debugging
        print("Error in quantum_optimize:", str(e))
        return jsonify({"error": str(e)}), 500

@app.route('/test', methods=['POST'])
def optimize_test():
    return "POST test request received successfully!"

@app.route('/subscribe', methods=['POST'])
def subscribe_user():
    """
    Create a subscription for a user using Stripe.
    """
    try:
        data = request.get_json()
        if not data or "email" not in data or "paymentMethodId" not in data:
            return jsonify({"error": "Invalid input: Missing 'email' or 'paymentMethodId'"}), 400

        email = data["email"]
        payment_method_id = data["paymentMethodId"]
        price_id = "price_1QYNn9JQZaUHxA2Ld9rV2MPd"  # Replace with your Stripe Price ID

        # Create or retrieve customer in Stripe
        customers = stripe.Customer.list(email=email).data
        if customers:
            customer = customers[0]  # Use the existing customer
        else:
            customer = stripe.Customer.create(email=email)

        # Attach the payment method to the customer
        stripe.PaymentMethod.attach(
            payment_method_id,
            customer=customer.id
        )

        # Set the payment method as the default for the customer
        stripe.Customer.modify(
            customer.id,
            invoice_settings={"default_payment_method": payment_method_id}
        )

        # Create the subscription
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


@app.route('/webhook', methods=['POST'])
def stripe_webhook():
    payload = request.get_data(as_text=True)
    sig_header = request.headers.get('Stripe-Signature')

    endpoint_secret = "your-webhook-signing-secret"  # Replace with your webhook secret
    event = None

    try:
        event = stripe.Webhook.construct_event(
            payload, sig_header, endpoint_secret
        )
    except ValueError as e:
        # Invalid payload
        return jsonify({"error": "Invalid payload"}), 400
    except stripe.error.SignatureVerificationError as e:
        # Invalid signature
        return jsonify({"error": "Invalid signature"}), 400

    # Handle the event
    if event['type'] == 'invoice.payment_succeeded':
        print("Payment succeeded:", event['data']['object'])
    elif event['type'] == 'customer.subscription.deleted':
        print("Subscription cancelled:", event['data']['object'])
    elif event['type'] == 'invoice.payment_failed':
        print("Payment failed:", event['data']['object'])

    return jsonify({"status": "success"}), 200

@app.route('/')
def index():
    return "Optimol API is running! Available endpoints: /test, /optimize, /quantum-optimize, /subscribe"

if __name__ == '__main__':
    with app.app_context():
        db.create_all()  # Create tables within the app context
    app.run(host="0.0.0.0", port=5000)
