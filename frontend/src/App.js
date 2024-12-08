import React, { useState } from "react";
import axios from "axios";
import apiBaseUrl from "./config";

const App = () => {
  const [optimizedMolecule, setOptimizedMolecule] = useState(null); // State to store optimized molecule
  const [moleculeData, setMoleculeData] = useState(null); // State to store uploaded molecule data
  const [showInstructions, setShowInstructions] = useState(false); // State to toggle instructions
  const [selectedMethod, setSelectedMethod] = useState(null); // Track selected optimization method

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
      setSelectedMethod("classical");
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
      setSelectedMethod("quantum");
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

          <h2>Optimization Details</h2>

          {/* Display Minimum Energy */}
          {optimizedMolecule.min_energy !== undefined && (
            <p>
              <strong>Minimum Energy:</strong>{" "}
              <span style={{ color: "#28a745" }}>
                {optimizedMolecule.min_energy.toFixed(6)}
              </span>
            </p>
          )}

          {/* Display Optimization Parameters */}
          <h3>Optimization Parameters:</h3>
          <ul>
            {selectedMethod === "classical" && (
              <>
                <li>
                  <strong>fmax:</strong> {fmax} (Classical)
                </li>
                <li>
                  <strong>steps:</strong> {steps} (Classical)
                </li>
              </>
            )}
            {selectedMethod === "quantum" && (
              <>
                <li>
                  <strong>maxiter:</strong> {maxiter} (Quantum)
                </li>
                <li>
                  <strong>p:</strong> {qaoaLayers} (Quantum)
                </li>
              </>
            )}
          </ul>

          {/* Display Optimal QAOA Parameters */}
          {selectedMethod === "quantum" && optimizedMolecule.optimal_params && (
            <>
              <h3>Optimal QAOA Parameters:</h3>
              <ul>
                {optimizedMolecule.optimal_params.map((param, index) => (
                  <li key={index}>
                    <strong>{index % 2 === 0 ? "γ" : "β"}:</strong>{" "}
                    {param.toFixed(6)}
                  </li>
                ))}
              </ul>
              <div>
                <p>The parameters represent the angles used in the QAOA ansatz:</p>
                <ul>
                  <li>
                    <strong>γ (Cost Operator Angles):</strong> Guides energy
                    minimization.
                  </li>
                  <li>
                    <strong>β (Mixer Operator Angles):</strong> Guides exploration of
                    feasible states.
                  </li>
                </ul>
              </div>
            </>
          )}
        </div>
      )}
    </div>
  );
};

export default App;
