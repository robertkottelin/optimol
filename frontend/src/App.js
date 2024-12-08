import React, { useState } from "react";
import axios from "axios";
import apiBaseUrl from "./config";


const App = () => {
  const [optimizedMolecule, setOptimizedMolecule] = useState(null); // State to store optimized molecule
  const [moleculeData, setMoleculeData] = useState(null); // State to store uploaded molecule data
  const [showInstructions, setShowInstructions] = useState(false); // State to toggle instructions

  const [fmax, setFmax] = useState(0.005);
  const [steps, setSteps] = useState(500);
  const [maxiter, setMaxiter] = useState(1000);
  const [qaoaLayers, setQaoaLayers] = useState(2);

  const handleFileUpload = (event) => {
    const file = event.target.files[0];

    if (!file) {
      alert("Please select a file.");
      return;
    }

    const reader = new FileReader();

    reader.onload = (e) => {
      try {
        const fileContent = e.target.result;
        const parsedData = JSON.parse(fileContent); // Parse the JSON file content

        if (!validateMoleculeJSON(parsedData)) {
          alert("Invalid molecule JSON format.");
          return;
        }

        setMoleculeData(parsedData); // Store parsed data for later use
        alert("File uploaded successfully. Choose an optimization method.");
      } catch (error) {
        console.error("Error processing the file:", error);
        alert("Error processing the file. Check the console for details.");
      }
    };

    reader.readAsText(file); // Read the file as text
  };

  const handleOptimize = async () => {
    if (!moleculeData) {
      alert("Please upload a file first.");
      return;
    }

    try {
      const payload = {
        file1: moleculeData.file1,
        fmax: fmax,
        steps: steps,
      };

      const response = await axios.post(`${apiBaseUrl}/optimize`, payload);
      setOptimizedMolecule(response.data.optimized_file1);
    } catch (error) {
      console.error("Error optimizing molecule:", error);
      alert("Error optimizing molecule. Check the console for details.");
    }
  };

  const handleQuantumOptimize = async () => {
    if (!moleculeData) {
      alert("Please upload a file first.");
      return;
    }

    try {
      const payload = {
        file1: moleculeData.file1,
        optimizer: "COBYLA",
        p: qaoaLayers,
        maxiter: maxiter,
      };

      const response = await axios.post(`${apiBaseUrl}/quantum-optimize`, payload);
      setOptimizedMolecule(response.data.optimized_file1);
    } catch (error) {
      console.error("Error optimizing molecule with quantum method:", error);
      alert("Error optimizing molecule with quantum method. Check the console for details.");
    }
  };

  const handleDownload = () => {
    if (!optimizedMolecule) {
      alert("No optimized molecule available to download.");
      return;
    }

    const jsonData = {
      file1: optimizedMolecule, // Wrap in the same format as input
    };
    const blob = new Blob([JSON.stringify(jsonData, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);

    const link = document.createElement("a");
    link.href = url;
    link.download = "optimized_molecule.json";
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  const validateMoleculeJSON = (data) => {
    // Basic validation to ensure the JSON structure matches the expected format
    return (
      data &&
      data.file1 &&
      Array.isArray(data.file1.atoms) &&
      data.file1.atoms.every(
        (atom) =>
          atom.id &&
          atom.element &&
          typeof atom.x === "number" &&
          typeof atom.y === "number" &&
          typeof atom.z === "number"
      )
    );
  };

  return (
    <div style={{ padding: "20px", textAlign: "center" }}>
      <h1>Optimize Molecule</h1>
      <h2>Based on Classical and Quantum Energy Calculations</h2>
      <input
        type="file"
        onChange={handleFileUpload}
        style={{
          padding: "10px",
          border: "2px dashed gray",
          cursor: "pointer",
          marginBottom: "20px",
        }}
      />
      <div style={{ display: "flex", gap: "20px" }}>
        <p>Classical optimization parameters:</p>
        <div>
          <label>fmax:</label>
          <input
            type="number"
            step="0.001"
            value={fmax}
            onChange={(e) => setFmax(Number(e.target.value))}
            placeholder="fmax"
          />
        </div>
        <div>
          <label>steps:</label>
          <input
            type="number"
            value={steps}
            onChange={(e) => setSteps(Number(e.target.value))}
            placeholder="steps"
          />
        </div>
        <p>Quantum optimization parameters:</p>
        <div>
          <label>maxiter:</label>
          <input
            type="number"
            value={maxiter}
            onChange={(e) => setMaxiter(Number(e.target.value))}
            placeholder="maxiter"
          />
        </div>
        <div>
          <label>QAOA Layers (p):</label>
          <input
            type="number"
            value={qaoaLayers}
            onChange={(e) => setQaoaLayers(Number(e.target.value))}
            placeholder="p"
          />
        </div>
      </div>
      <div
        style={{
          marginTop: "20px",
          display: "flex",
          justifyContent: "center",
          gap: "10px",
        }}
      >
        <button
          onClick={handleOptimize}
          style={{
            padding: "10px 20px",
            backgroundColor: "#007bff",
            color: "white",
            border: "none",
            borderRadius: "5px",
            cursor: "pointer",
          }}
        >
          Classical Optimize
        </button>
        <button
          onClick={handleQuantumOptimize}
          style={{
            padding: "10px 20px",
            backgroundColor: "#28a745",
            color: "white",
            border: "none",
            borderRadius: "5px",
            cursor: "pointer",
          }}
        >
          Quantum Optimize
        </button>
        <button
          onClick={() => setShowInstructions(!showInstructions)}
          style={{
            padding: "10px 20px",
            backgroundColor: "#ffc107",
            color: "black",
            border: "none",
            borderRadius: "5px",
            cursor: "pointer",
          }}
        >
          How To Use
        </button>
        <button
          onClick={handleDownload}
          style={{
            padding: "10px 20px",
            backgroundColor: "#6c757d",
            color: "white",
            border: "none",
            borderRadius: "5px",
            cursor: "pointer",
          }}
        >
          Download Optimized Molecule
        </button>
      </div>

      {showInstructions && (
        <div
          style={{
            marginTop: "20px",
            textAlign: "left",
            padding: "10px",
            border: "1px solid gray",
            borderRadius: "5px",
            backgroundColor: "#f9f9f9",
            maxWidth: "600px",
            margin: "auto",
          }}
        >
          <h2>How To Use This App</h2>
          <p>
            This app optimizes molecular geometry using either classical or quantum energy optimization techniques.
          </p>
          <h3>Steps to Use:</h3>
          <ol>
            <li>Prepare a JSON file containing your molecule data in the following format:</li>
            <pre>
              {`{
  "file1": {
    "atoms": [
      { "id": 1, "element": "H", "x": 0.0, "y": 0.0, "z": 0.0 },
      { "id": 2, "element": "H", "x": 0.0, "y": 0.0, "z": 0.74 },
      { "id": 3, "element": "O", "x": 0.0, "y": 0.74, "z": 0.0 }
    ]
  }
}`}
            </pre>
            <li>Upload your JSON file using the file upload button.</li>
            <li>Click on "Classical Optimize" or "Quantum Optimize" based on your preference.</li>
            <li>
              The results will be displayed below, including the optimized molecular geometry and, for quantum optimization, additional details like energy and parameters.
            </li>
          </ol>
          <h3>Optimization Methods:</h3>
          <ul>
            <li>
              <strong>Classical Optimize:</strong> Uses classical computational chemistry methods to minimize the energy of the molecular geometry.
              <ul>
                <li>
                  <strong>fmax:</strong> The force convergence criterion. Lower values require finer adjustments (default: 0.005).
                </li>
                <li>
                  <strong>steps:</strong> The maximum number of iterations for the optimizer (default: 500).
                </li>
              </ul>
            </li>
            <li>
              <strong>Quantum Optimize:</strong> Employs a quantum optimization algorithm (QAOA) to explore energy minimization with quantum methods.
              <ul>
                <li>
                  <strong>maxiter:</strong> The maximum number of iterations for the optimizer during parameter tuning (default: 1000).
                </li>
                <li>
                  <strong>p:</strong> The number of QAOA layers, determining the complexity of the quantum ansatz (default: 10).
                </li>
              </ul>
            </li>
          </ul>
          <h3>Download Optimized Molecule:</h3>
          <p>
            After optimization, you can download the optimized molecule data as a JSON file for further analysis or visualization.
          </p>
          <h3>Upcoming features:</h3>
          <ul>
            <li>Visualize the molecule geometry before and after optimization.</li>
            <li>
              <s>Customize optimization parameters (fmax, steps, maxiter, QAOA layers).</s>
            </li>            
            <li>Support for additional optimization methods (e.g., VQE, Nelder-Mead).</li>
            <li>Support for additional file formats (e.g., XYZ, PDB).</li>
            <li>Compare classical and quantum optimization results.</li>
            <li>Optimize multiple molecules simultaneously.</li>
          </ul>
        </div>
      )}

      {optimizedMolecule && (
        <div
          style={{
            marginTop: "20px",
            textAlign: "left",
            padding: "10px",
            border: "1px solid gray",
            borderRadius: "5px",
            backgroundColor: "#f9f9f9",
            maxWidth: "600px",
            margin: "auto",
          }}
        >
          <h2>Optimized Molecule</h2>
          <pre>{JSON.stringify(optimizedMolecule.atoms, null, 2)}</pre>

          {optimizedMolecule.min_energy !== undefined && (
            <>
              <h2>Optimization Details</h2>
              <p>
                <strong>Minimum Energy:</strong> {optimizedMolecule.min_energy.toFixed(6)}
              </p>
              <h3>Optimization Parameters:</h3>
              <ul>
                <li><strong>fmax:</strong> 0.02 (Classical)</li>
                <li><strong>steps:</strong> 100 (Classical)</li>
                <li><strong>maxiter:</strong> 500 (Quantum)</li>
                <li><strong>p:</strong> 10 (Quantum)</li>
              </ul>
              <h3>Optimal QAOA Parameters:</h3>
              <ul>
                {optimizedMolecule.optimal_params.map((param, index) => (
                  <li key={index}>
                    <strong>{index % 2 === 0 ? "γ" : "β"}:</strong> {param.toFixed(6)}
                  </li>
                ))}
              </ul>
              <p>
                The parameters represent the angles used in the QAOA ansatz:
                <ul>
                  <li><strong>γ (Cost Operator Angles):</strong> Guides energy minimization.</li>
                  <li><strong>β (Mixer Operator Angles):</strong> Guides exploration of feasible states.</li>
                </ul>
              </p>
            </>
          )}
        </div>
      )}
    </div>
  );
}
export default App;
