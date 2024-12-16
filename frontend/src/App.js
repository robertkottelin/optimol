import React, { useState, useRef, useEffect } from "react";
import axios from "axios";
import * as $3Dmol from '3dmol';
console.log("$3Dmol loaded:", $3Dmol);

const App = () => {
  const [optimizedMolecule, setOptimizedMolecule] = useState(null); // State to store optimized molecule
  const [moleculeData, setMoleculeData] = useState(null); // State to store uploaded molecule data
  const [showInstructions, setShowInstructions] = useState(false); // State to toggle instructions

  const [fmax, setFmax] = useState(0.005);
  const [steps, setSteps] = useState(500);
  const [maxiter, setMaxiter] = useState(1000);
  const [qaoaLayers, setQaoaLayers] = useState(2);

  const apiBaseUrl = "http://localhost:5000/";

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
        const parsedData = JSON.parse(fileContent);

        if (!validateMoleculeJSON(parsedData)) {
          alert("Invalid molecule JSON format.");
          return;
        }

        setMoleculeData(parsedData.file1);
        setOptimizedMolecule(null); // Reset optimized molecule when a new file is uploaded
      } catch (error) {
        console.error("Error processing the file:", error);
        alert("Error processing the file. Check the console for details.");
      }
    };

    reader.readAsText(file);
  };

  const handleOptimize = async () => {
    if (!moleculeData) {
      alert("Please upload a file first.");
      return;
    }

    try {
      const payload = {
        file1: moleculeData,
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
        file1: moleculeData,
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
      file1: optimizedMolecule,
    };
    const blob = new Blob([JSON.stringify(jsonData, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);

    const link = document.createElement("a");
    link.href = url;
    link.download = "optimized_molecule.json";
    link.click();
    URL.revokeObjectURL(url);
  };

  const validateMoleculeJSON = (data) => {
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

  const MoleculeViewer = ({ moleculeData }) => {
    const viewerRef = useRef();

    useEffect(() => {
      if (!moleculeData || !moleculeData.atoms) {
        console.warn("No molecule data found for visualization.");
        return;
      }

      // Initialize the viewer
      const viewer = $3Dmol.createViewer(viewerRef.current, {
        backgroundColor: "white",
      });
      viewer.clear();

      try {
        // Map molecule data to atoms for the viewer
        const atoms = moleculeData.atoms.map((atom) => ({
          elem: atom.element,
          x: atom.x,
          y: atom.y,
          z: atom.z,
        }));

        const model = viewer.addModel(); // Add a new model to the viewer
        model.addAtoms(atoms); // Add atoms to the model

        // Set styles for atoms and sticks (bonds)
        viewer.setStyle({}, {
          sphere: { radius: 0.2 }, // Smaller radius for atoms
          stick: { radius: 0.2 }, // Add sticks (bonds) between atoms with thinner lines
          line: { linewidth: 2.0 }, // Optional: add thin lines for bonds
        });
        viewer.zoomTo(); // Automatically adjust zoom to fit molecule
        viewer.render();
      } catch (error) {
        console.error("Error rendering molecule:", error);
      }
    }, [moleculeData]);

    return (
      <div
        ref={viewerRef}
        style={{
          width: "100%",
          height: "300px",
          position: "relative", // Ensure correct alignment
          border: "1px solid #ccc",
          margin: "10px auto", // Add margin for positioning
        }}
      ></div>
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

      {moleculeData && (
        <>
          <h3 style={{ marginBottom: "10px" }}>Molecule Visualization:</h3>
          <MoleculeViewer moleculeData={optimizedMolecule || moleculeData} />
        </>
      )}

      <div style={{ display: "flex", gap: "20px", justifyContent: "center", flexWrap: "wrap" }}>
        <p>Classical parameters:</p>
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
        <p>Quantum parameters:</p>
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
      <div style={{ marginTop: "20px", display: "flex", justifyContent: "center", gap: "10px" }}>
        <button onClick={handleOptimize}>Classical Optimize</button>
        <button onClick={handleQuantumOptimize}>Quantum Optimize</button>
        <button onClick={() => setShowInstructions(!showInstructions)}>How To Use</button>
        <button onClick={handleDownload}>Download Optimized Molecule</button>
      </div>
    </div>
  );
};

export default App;
