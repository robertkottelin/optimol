from flask import Flask, request, jsonify

app = Flask(__name__)

def calculate_energy_and_optimize(molecule1):
    """
    Placeholder for energy calculation and coordinate optimization.
    Modify the coordinates slightly for demonstration.
    """
    for atom in molecule1["atoms"]:
        atom["x"] += 0.1
        atom["y"] += 0.1
        atom["z"] += 0.1

    return molecule1

@app.route('/optimize', methods=['POST'])
def optimize():
    data = request.get_json()

    if "file1" not in data:
        return jsonify({"error": "'file1' must be provided as JSON."}), 400

    molecule1 = data["file1"]

    # Perform optimization
    optimized_molecule1 = calculate_energy_and_optimize(molecule1)

    return jsonify({
        "optimized_file1": optimized_molecule1
    })

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)
