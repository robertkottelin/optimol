from flask import Flask, request, jsonify, send_file
from tempfile import NamedTemporaryFile
import os

app = Flask(__name__)

def calculate_energy_and_optimize(mol2_data1, mol2_data2):
    """
    Placeholder for energy calculation and coordinate optimization.
    Returns updated MOL2 content for both molecules.
    """
    # Mock optimization: shift the coordinates slightly
    optimized_mol2_1 = mol2_data1.replace("0.0000", "0.1000")
    optimized_mol2_2 = mol2_data2.replace("1.5000", "1.6000")
    return optimized_mol2_1, optimized_mol2_2

@app.route('/optimize', methods=['POST'])
def optimize():
    if 'file1' not in request.files or 'file2' not in request.files:
        return jsonify({"error": "Both 'file1' and 'file2' must be provided."}), 400

    file1 = request.files['file1']
    file2 = request.files['file2']

    mol2_data1 = file1.read().decode('utf-8')
    mol2_data2 = file2.read().decode('utf-8')

    # Perform optimization
    optimized_mol2_1, optimized_mol2_2 = calculate_energy_and_optimize(mol2_data1, mol2_data2)

    # Write the optimized MOL2 files to temporary files
    temp_file1 = NamedTemporaryFile(delete=False, suffix=".mol2")
    temp_file2 = NamedTemporaryFile(delete=False, suffix=".mol2")
    try:
        with open(temp_file1.name, 'w') as f:
            f.write(optimized_mol2_1)
        with open(temp_file2.name, 'w') as f:
            f.write(optimized_mol2_2)

        # Return the files as downloadable attachments
        return jsonify({
            "optimized_file1": send_file(temp_file1.name, as_attachment=True, download_name="optimized1.mol2"),
            "optimized_file2": send_file(temp_file2.name, as_attachment=True, download_name="optimized2.mol2")
        })
    finally:
        # Clean up the temporary files after the response is sent
        os.unlink(temp_file1.name)
        os.unlink(temp_file2.name)

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)
