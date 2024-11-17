from flask import Flask, request, jsonify
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT

app = Flask(__name__)
from flask_cors import CORS
CORS(app, origins=["*"], methods=["GET", "POST", "OPTIONS"], allow_headers=["Content-Type"])

def optimize_molecule(data):
    """
    Optimize the molecular geometry based on quantum energy calculations.
    This example uses ASE and EMT for demonstration.
    """
    # Convert input data to ASE Atoms object
    atoms = []
    for atom in data["atoms"]:
        atoms.append((atom["element"], (atom["x"], atom["y"], atom["z"])))

    molecule = Atoms(symbols=[a[0] for a in atoms], positions=[a[1] for a in atoms])

    # Assign a basic calculator (EMT) for optimization
    molecule.calc = EMT()

    # Optimize the geometry using BFGS algorithm
    optimizer = BFGS(molecule)
    optimizer.run(fmax=0.02)  # Stop when forces are below 0.02 eV/Å

    # Extract optimized positions
    optimized_atoms = []
    for i, atom in enumerate(molecule):
        optimized_atoms.append({
            "id": i + 1,
            "element": atom.symbol,
            "x": atom.position[0],
            "y": atom.position[1],
            "z": atom.position[2],
        })

    return {"atoms": optimized_atoms}

@app.route('/optimize', methods=['POST'])
def optimize():
    data = request.get_json()
    if not data or "file1" not in data:
        return jsonify({"error": "Invalid input"}), 400

    try:
        optimized_data = optimize_molecule(data["file1"])
        return jsonify({"optimized_file1": optimized_data})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('https://molcompute-bcbmhphmb7h5e0cy.northeurope-01.azurewebsites.net/')
def home():
    return "Flask app is running!", 200


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)
