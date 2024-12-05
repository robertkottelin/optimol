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

app = Flask(__name__)
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
    optimizer = BFGS(molecule)
    optimizer.run(fmax=0.02)

    optimized_atoms = [
        {"id": i + 1, "element": atom.symbol, "x": atom.position[0], "y": atom.position[1], "z": atom.position[2]}
        for i, atom in enumerate(molecule)
    ]
    return {"atoms": optimized_atoms}

from qiskit.primitives import Sampler
from qiskit.circuit.library import QAOAAnsatz
from qiskit.circuit import Parameter

from scipy.optimize import minimize
from scipy.optimize import minimize
def optimize_molecule_qaoa(data, optimizer='COBYLA', p=2):
    """
    Optimize a molecular geometry using QAOA and return adjusted positions.
    """
    try:
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
            options={"maxiter": 500}
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


@app.route('/optimize', methods=['POST'])
def classical_optimize():
    """
    Optimize molecular geometry using classical methods.
    """
    data = request.get_json()
    if not data or "file1" not in data:
        return jsonify({"error": "Invalid input"}), 400

    try:
        optimized_data = optimize_molecule(data["file1"])
        return jsonify({"optimized_file1": optimized_data})
    except Exception as e:
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

@app.route('/')
def index():
    return "Optimol API is running!"


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)
