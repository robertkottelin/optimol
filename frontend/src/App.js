import React, { useState, useRef, useEffect } from "react";
import axios from "axios";
import * as $3Dmol from '3dmol';
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import SubscriptionForm from './SubscriptionForm';
console.log("$3Dmol loaded:", $3Dmol);

const stripePromise = loadStripe('pk_test_kbl0ETzPsoiTwU4ZJMvhsYJw006XVnV4Aq');

const App = () => {
  const [optimizedMolecule, setOptimizedMolecule] = useState(null); // State to store optimized molecule
  const [moleculeData, setMoleculeData] = useState(null); // State to store uploaded molecule data
  const [showInstructions, setShowInstructions] = useState(false); // State to toggle instructions
  const [isSubscribed, setIsSubscribed] = useState(false); // Track subscription status
  const [userEmail, setUserEmail] = useState(""); // Track user email

  const [fmax, setFmax] = useState(0.005);
  const [steps, setSteps] = useState(500);
  const [maxiter, setMaxiter] = useState(1000);
  const [qaoaLayers, setQaoaLayers] = useState(2);

  const handleSubscriptionSuccess = (email) => {
    setIsSubscribed(true);
    setUserEmail(email);
    localStorage.setItem("userEmail", email); // Save email in localStorage
  };

  const apiBaseUrl = "http://localhost:5000/";

  // Check subscription status on page load
  useEffect(() => {
    const email = localStorage.getItem("userEmail");
    if (email) {
      checkSubscriptionStatus(email);
    }
  }, []);

  const checkSubscriptionStatus = async (email) => {
    try {
      const response = await axios.post(`${apiBaseUrl}/check-subscription`, { email });
      if (response.data.isSubscribed) {
        setIsSubscribed(true);
        setUserEmail(email);
      } else {
        setIsSubscribed(false);
      }
    } catch (error) {
      console.error("Error checking subscription status:", error);
    }
  };

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
        backgroundColor: "grey",
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

        // Add labels for each atom
        moleculeData.atoms.forEach((atom) => {
          viewer.addLabel(atom.element, {
            position: { x: atom.x, y: atom.y, z: atom.z },
            fontSize: 12,
            fontColor: "white",
            backgroundColor: "black",
            borderRadius: 2,
            padding: 2,
            inFront: true, // Ensure labels are always visible in front
          });
        });

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
      {!isSubscribed ? (
        <Elements stripe={stripePromise}>
          <SubscriptionForm onSuccess={handleSubscriptionSuccess} />
        </Elements>
      ) : (
        <div>
          <p>Welcome, {userEmail}! You are subscribed.</p>
          {/* Render the main app features here */}
        </div>
      )}
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

      {showInstructions && (
        <div className="box">
          <h2>How To Use:</h2>
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
            <li>Modify optimization parameters (see description below).</li>
            <li>Click "Classical Optimize" or "Quantum Optimize".</li>
            <li>Download the optimized molecule data if needed.</li>
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
        </div>
      )}

      {optimizedMolecule && (
        <div className="box">
          <h2>Optimized Molecule</h2>
          <pre>{JSON.stringify(optimizedMolecule.atoms, null, 2)}</pre>
          <h3>Optimization Details</h3>
          <p>
            <strong>Minimum Energy:</strong> {optimizedMolecule.min_energy?.toFixed(6)}
          </p>
          <h3>Optimization Parameters:</h3>
          <ul>
            <li>
              <strong>fmax:</strong> {fmax} (Classical)
            </li>
            <li>
              <strong>steps:</strong> {steps} (Classical)
            </li>
            <li>
              <strong>maxiter:</strong> {maxiter} (Quantum)
            </li>
            <li>
              <strong>QAOA Layers (p):</strong> {qaoaLayers} (Quantum)
            </li>
          </ul>
        </div>
      )}

    </div>
  );
};


export default App;
