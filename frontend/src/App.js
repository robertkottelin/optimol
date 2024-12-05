import React, { useState } from "react";
import axios from "axios";
import apiBaseUrl from "./config";


const App = () => {
  const [optimizedMolecule, setOptimizedMolecule] = useState(null); // State to store optimized molecule
  const [moleculeData, setMoleculeData] = useState(null); // State to store uploaded molecule data
  const [showInstructions, setShowInstructions] = useState(false); // State to toggle instructions

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
      const response = await axios.post(`${apiBaseUrl}/optimize`, moleculeData);
      setOptimizedMolecule(response.data.optimized_file1); // Set optimized molecule in state
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
        file1: {
          atoms: moleculeData.file1.atoms, // Explicitly pass atoms array
        },
        optimizer: "COBYLA", // Default optimizer
        p: 2, // Default number of QAOA layers
      };
  
      const response = await axios.post(`${apiBaseUrl}/quantum-optimize`, payload);
      setOptimizedMolecule(response.data.optimized_file1); // Set optimized molecule in state
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
      <h1>Molecule Optimizer</h1>
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
          How To
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
            This app optimizes molecular geometry using either classical or quantum optimization techniques.
          </p>
          <h3>Steps to Use:</h3>
          <ol>
            <li>Prepare a JSON file containing your molecule data in the following format:</li>
            <pre>
              {`{
  "file1": {
    "atoms": [
      { "id": 1, "element": "H", "x": 0.0, "y": 0.0, "z": 0.0 },
      { "id": 2, "element": "H", "x": 0.0, "y": 0.0, "z": 0.74 }
    ]
  }
}`}
            </pre>
            <li>Upload your JSON file using the file upload button.</li>
            <li>Click on "Classical Optimize" or "Quantum Optimize" based on your preference.</li>
            <li>
              The results will be displayed below, including the optimized molecular geometry and, for quantum
              optimization, additional details like energy and parameters.
            </li>
          </ol>
          <h3>Optimization Methods:</h3>
          <ul>
            <li>
              <strong>Classical Optimize:</strong> Uses classical computational chemistry methods to minimize the energy
              of the molecular geometry.
            </li>
            <li>
              <strong>Quantum Optimize:</strong> Employs a quantum optimization algorithm (QAOA) to explore energy
              minimization with quantum methods.
            </li>
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
};

export default App;
