import random
import traceback
from flask import Flask, request, jsonify
from flask_cors import CORS
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.cp2k import CP2K
import numpy as np
import os
import stripe
from dotenv import load_dotenv
from flask_sqlalchemy import SQLAlchemy
import tempfile
import json

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
CP2K_EXECUTABLE = os.getenv("CP2K_EXECUTABLE", "/usr/bin/cp2k.popt")

stripe.api_key = STRIPE_SECRET
CORS(app, origins=["http://localhost:3000"], methods=["GET", "POST", "OPTIONS"], allow_headers=["Content-Type"])

def setup_cp2k_calculator(basis_set="DZVP-MOLOPT-SR-GTH", 
                         functional="PBE", 
                         charge=0, 
                         uks=False):
    """
    Set up a CP2K calculator with specified parameters.
    """
    # Create input template
    inp = """&GLOBAL
  PROJECT cp2k-calc
  RUN_TYPE ENERGY_FORCE
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    &QS
      METHOD {}
    &END QS
    &POISSON
      PERIODIC NONE
      PSOLVER WAVELET
    &END POISSON
    &SCF
      MAX_SCF 200
      EPS_SCF 1.0E-6
    &END SCF
  &END DFT
&END FORCE_EVAL
""".format(functional)
    
    # Create calculator with input template
    calc = CP2K(
        label="cp2k-calc",
        command=f"{CP2K_EXECUTABLE} -i PREFIX.inp -o PREFIX.out",
        inp=inp,
        charge=charge,
        uks=uks
    )
    
    return calc

def optimize_molecule(data):
    """
    Optimize the molecular geometry using CP2K's DFT implementation.
    """
    # Extract atoms from input data
    atoms_data = data["atoms"]
    
    # Validate input structure
    if not atoms_data or len(atoms_data) < 2:
        raise ValueError("Invalid input: at least two atoms required for optimization")
    
    # Convert to ASE Atoms object
    symbols = [atom["element"] for atom in atoms_data]
    positions = np.array([[atom["x"], atom["y"], atom["z"]] for atom in atoms_data])
    
    molecule = Atoms(symbols=symbols, positions=positions)
    
    # Check for valid geometry
    if any(np.linalg.norm(atom1 - atom2) < 0.7 for atom1, atom2 in zip(molecule.positions[:-1], molecule.positions[1:])):
        raise ValueError("Initial geometry contains atoms that are too close")
    
    # Retrieve user-specified parameters or use defaults
    parameters = {
        "fmax": data.get("fmax", 0.001),  # Force convergence criterion in eV/Å
        "steps": data.get("steps", 200),   # Maximum optimization steps
        "basis_set": data.get("basis_set", "DZVP-MOLOPT-SR-GTH"),
        "functional": data.get("functional", "PBE"),
        "charge": data.get("charge", 0),
        "spin_polarized": data.get("spin_polarized", False)
    }
    
    # Set up CP2K calculator
    molecule.calc = setup_cp2k_calculator(
        basis_set=parameters["basis_set"],
        functional=parameters["functional"],
        charge=parameters["charge"],
        uks=parameters["spin_polarized"]
    )
    
    # Run geometry optimization
    optimizer = BFGS(molecule)
    optimizer.run(fmax=parameters["fmax"], steps=parameters["steps"])
    
    # Extract final energies
    total_energy = molecule.get_potential_energy()
    
    # Format optimized structure
    optimized_atoms = [
        {
            "id": i + 1, 
            "element": atom.symbol, 
            "x": float(atom.position[0]), 
            "y": float(atom.position[1]), 
            "z": float(atom.position[2])
        }
        for i, atom in enumerate(molecule)
    ]
    
    # Return optimized structure and energy
    return {
        "atoms": optimized_atoms,
        "energy": float(total_energy),
        "converged": optimizer.converged(),
        "parameters": parameters
    }

def setup_qmmm_calculation(protein_data, ligand_data, qm_indices=None):
    """
    Set up QM/MM calculation with CP2K for a protein-ligand system.
    """
    # Convert protein and ligand data to ASE Atoms objects
    protein_atoms = []
    for atom in protein_data["atoms"]:
        protein_atoms.append((atom["element"], (atom["x"], atom["y"], atom["z"])))
    
    ligand_atoms = []
    for atom in ligand_data["atoms"]:
        ligand_atoms.append((atom["element"], (atom["x"], atom["y"], atom["z"])))
    
    # Combine into one system
    all_symbols = [a[0] for a in protein_atoms + ligand_atoms]
    all_positions = [a[1] for a in protein_atoms + ligand_atoms]
    
    # Create ASE Atoms object
    system = Atoms(symbols=all_symbols, positions=all_positions)
    
    # Determine QM region (default to ligand atoms)
    if qm_indices is None:
        qm_indices = list(range(len(protein_atoms), len(protein_atoms) + len(ligand_atoms)))
    
    # Generate QM indices string for CP2K input
    qm_kinds = {}
    for idx in qm_indices:
        elem = all_symbols[idx]
        if elem not in qm_kinds:
            qm_kinds[elem] = []
        qm_kinds[elem].append(idx + 1)  # CP2K indices are 1-based
    
    # Build QM_KIND blocks for the input
    qm_kind_blocks = ""
    for elem, indices in qm_kinds.items():
        indices_str = " ".join(str(i) for i in indices)
        qm_kind_blocks += f"""    &QM_KIND {elem}
      MM_INDEX {indices_str}
    &END QM_KIND
"""
    
    # Create input template for QM/MM calculation
    inp = f"""&GLOBAL
  PROJECT cp2k-qmmm
  RUN_TYPE ENERGY_FORCE
&END GLOBAL
&FORCE_EVAL
  METHOD QMMM
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    &QS
      METHOD PBE
    &END QS
    &POISSON
      PERIODIC NONE
      PSOLVER WAVELET
    &END POISSON
    &SCF
      MAX_SCF 200
      EPS_SCF 1.0E-6
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &QMMM
    &CELL
      ABC 40.0 40.0 40.0
    &END CELL
{qm_kind_blocks}
  &END QMMM
&END FORCE_EVAL
"""
    
    # Setup CP2K calculator with only validated parameters
    calc = CP2K(
        label="cp2k-qmmm",
        command=f"{CP2K_EXECUTABLE} -i PREFIX.inp -o PREFIX.out",
        inp=inp
    )
    
    system.calc = calc
    return system

def optimize_binding(data):
    """
    Optimize protein-ligand binding using QM/MM approach.
    """
    try:
        if "protein" not in data or "ligand" not in data:
            raise ValueError("Both 'protein' and 'ligand' data required for binding optimization")
        
        # Extract QM region indices if specified
        qm_indices = data.get("qm_indices", None)
        
        # Set up QM/MM calculation
        system = setup_qmmm_calculation(data["protein"], data["ligand"], qm_indices)
        
        # Perform energy minimization
        fmax = data.get("fmax", 0.05)  # Higher tolerance for protein-ligand systems
        steps = data.get("steps", 100)  # Fewer steps for large systems
        
        optimizer = BFGS(system)
        optimizer.run(fmax=fmax, steps=steps)
        
        # Extract binding energy
        # In a real implementation, this would require multiple calculations
        # to properly account for all energy terms
        binding_energy = system.get_potential_energy()
        
        # Extract optimized coordinates
        num_protein_atoms = len(data["protein"]["atoms"])
        
        optimized_protein = [
            {
                "id": i + 1, 
                "element": atom.symbol, 
                "x": float(atom.position[0]), 
                "y": float(atom.position[1]), 
                "z": float(atom.position[2])
            }
            for i, atom in enumerate(system[:num_protein_atoms])
        ]
        
        optimized_ligand = [
            {
                "id": i + 1, 
                "element": atom.symbol, 
                "x": float(atom.position[0]), 
                "y": float(atom.position[1]), 
                "z": float(atom.position[2])
            }
            for i, atom in enumerate(system[num_protein_atoms:])
        ]
        
        return {
            "protein": {"atoms": optimized_protein},
            "ligand": {"atoms": optimized_ligand},
            "binding_energy": float(binding_energy),
            "converged": optimizer.converged()
        }
    
    except Exception as e:
        raise ValueError(f"Error in binding optimization: {str(e)}")

@app.route('/optimize', methods=['POST'])
def classical_optimize():
    """
    Endpoint for classical molecular optimization using CP2K.
    """
    try:
        data = request.get_json()
        if not data or "file1" not in data:
            return jsonify({"error": "Invalid input: Missing 'file1' in request body"}), 400

        # Run optimization
        optimized_data = optimize_molecule(data["file1"])
        
        # Store optimization in database for subscribed users
        if "email" in data:
            user = User.query.filter_by(email=data["email"]).first()
            if user and user.subscription_status == "active":
                optimization = Optimization(
                    user_id=user.id,
                    optimization_type="classical_cp2k",
                    parameters=json.dumps(data["file1"]),
                    result=json.dumps(optimized_data)
                )
                db.session.add(optimization)
                db.session.commit()
        
        return jsonify({"optimized_file1": optimized_data})
    
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/binding-optimize', methods=['POST'])
def binding_optimize_endpoint():
    """
    Endpoint for protein-ligand binding optimization using CP2K QM/MM.
    """
    try:
        data = request.get_json()
        if not data or "protein" not in data or "ligand" not in data:
            return jsonify({"error": "Invalid input: Both 'protein' and 'ligand' required"}), 400

        # Check subscription for premium features
        if "email" in data:
            user = User.query.filter_by(email=data["email"]).first()
            if not user or user.subscription_status != "active":
                return jsonify({"error": "Premium subscription required for binding optimization"}), 403
        else:
            return jsonify({"error": "User email required for authentication"}), 401

        # Use real implementation for binding optimization
        result = optimize_binding(data)
        
        # Store in database
        if "email" in data:
            user = User.query.filter_by(email=data["email"]).first()
            optimization = Optimization(
                user_id=user.id,
                optimization_type="binding_qmmm",
                parameters=json.dumps(data),
                result=json.dumps(result)
            )
            db.session.add(optimization)
            db.session.commit()
        
        return jsonify({"binding_result": result})
    
    except Exception as e:
        app.logger.error(f"Endpoint error: {str(e)}")
        app.logger.error(f"Traceback: {traceback.format_exc()}")
        return jsonify({"error": str(e)}), 500
                
@app.route('/quantum-optimize', methods=['POST'])
def quantum_optimize():
    """
    Enhanced quantum-optimize endpoint using CP2K for QM energy calculations.
    """
    try:
        data = request.get_json()
        if not data or "file1" not in data:
            return jsonify({"error": "Invalid input: Missing 'file1' in request body"}), 400

        # Check subscription for premium features
        if "email" in data:
            user = User.query.filter_by(email=data["email"]).first()
            if not user or user.subscription_status != "active":
                return jsonify({"error": "Premium subscription required for quantum optimization"}), 403
        
        # Use CP2K's DFT for quantum chemical calculations instead of QAOA
        # This is more appropriate for molecular energy calculations
        optimized_data = optimize_molecule(data["file1"])
        
        # Perform additional wavefunction analysis for quantum properties
        # In a real implementation, you would run additional CP2K calculations here
        
        return jsonify({"optimized_file1": optimized_data})
    
    except Exception as e:
        return jsonify({"error": str(e)}), 500

# Include the existing routes for Stripe integration
@app.route('/subscribe', methods=['POST'])
def subscribe_user():
    try:
        data = request.get_json()
        
        if not data or "email" not in data or "paymentMethodId" not in data:
            return jsonify({"error": "Invalid input: Missing 'email' or 'paymentMethodId'"}), 400

        email = data["email"]
        payment_method_id = data["paymentMethodId"]
        price_id = "price_1QYNn9JQZaUHxA2Ld9rV2MPd"  # Replace with your Stripe Price ID

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
    return "Optimol API with CP2K integration is running! Available endpoints: /optimize, /binding-optimize, /quantum-optimize, /subscribe"

if __name__ == '__main__':
    with app.app_context():
        db.create_all()  # Create tables within the app context
    app.run(host="0.0.0.0", port=5000)