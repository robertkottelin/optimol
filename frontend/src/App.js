import React, { useState, useRef, useEffect } from "react";
import axios from "axios";
import * as $3Dmol from '3dmol';
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import ReactMarkdown from "react-markdown";
import SubscriptionForm from './SubscriptionForm';
console.log("$3Dmol loaded:", $3Dmol);

const stripePromise = loadStripe('pk_test_kbl0ETzPsoiTwU4ZJMvhsYJw006XVnV4Aq');

const DEFAULT_COMPUTE_LIMITS = {
  fmax: { min: 0.05, max: 0.05 }, // Fixed for non-subscribers
  steps: { min: 10, max: 100 },
  maxiter: { min: 10, max: 100 },
  qaoaLayers: { min: 1, max: 1 },
};

const SUBSCRIBED_COMPUTE_LIMITS = {
  fmax: { min: 0.001, max: 0.1 },
  steps: { min: 10, max: 1000000 },
  maxiter: { min: 10, max: 1000000 },
  qaoaLayers: { min: 1, max: 10000 },
};

const App = () => {
  const [optimizedMolecule, setOptimizedMolecule] = useState(null); // State to store optimized molecule
  const [moleculeData, setMoleculeData] = useState(null); // State to store uploaded molecule data
  const [isSubscribed, setIsSubscribed] = useState(false); // Track subscription status
  const [userEmail, setUserEmail] = useState(""); // Track user email
  const [limits, setLimits] = useState(DEFAULT_COMPUTE_LIMITS);
  const [isHowToUseVisible, setIsHowToUseVisible] = useState(false); // State to toggle the popup
  const [howToUseContent, setHowToUseContent] = useState(""); // Store the content of the readme


  // Loading states for each button
  const [isSubscribeLoading, setIsSubscribeLoading] = useState(false);
  const [isCancelLoading, setIsCancelLoading] = useState(false);
  const [isOptimizeLoading, setIsOptimizeLoading] = useState(false);
  const [isQuantumOptimizeLoading, setIsQuantumOptimizeLoading] = useState(false);

  // Define the state variables for compute limits
  const [fmax, setFmax] = useState(limits.fmax.min);
  const [steps, setSteps] = useState(limits.steps.min);
  const [maxiter, setMaxiter] = useState(limits.maxiter.min);
  const [qaoaLayers, setQaoaLayers] = useState(limits.qaoaLayers.min);

  // Fetch "how-to-use.md" content on demand
  const fetchHowToUse = async () => {
    try {
      const response = await axios.get('/how-to-use.md');
      setHowToUseContent(response.data);
      setIsHowToUseVisible(true); // Show the popup
    } catch (error) {
      console.error("Error fetching 'how-to-use.md':", error);
      alert("Failed to load instructions. Check console for details.");
    }
  };

  const handleClosePopup = () => {
    setIsHowToUseVisible(false);
  };

  const handleParameterChange = (setter, value, min, max) => {
    const parsedValue = parseFloat(value);
    if (parsedValue >= min && parsedValue <= max) {
      setter(parsedValue);
    }
  };

  const handleSubscriptionSuccess = (email) => {
    setIsSubscribed(true);
    setUserEmail(email);
    localStorage.setItem("userEmail", email);
    setLimits(SUBSCRIBED_COMPUTE_LIMITS);
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
        setLimits(SUBSCRIBED_COMPUTE_LIMITS);
      } else {
        setIsSubscribed(false);
        setLimits(DEFAULT_COMPUTE_LIMITS);
      }
    } catch (error) {
      console.error("Error checking subscription status:", error);
    }
  };

  const applyDefaultComputeLimits = () => {
    setFmax(DEFAULT_COMPUTE_LIMITS.fmax);
    setSteps(DEFAULT_COMPUTE_LIMITS.steps);
    setMaxiter(DEFAULT_COMPUTE_LIMITS.maxiter);
    setQaoaLayers(DEFAULT_COMPUTE_LIMITS.qaoaLayers);
  };

  const applyFullComputeLimits = () => {
    setFmax(0.005); // Higher precision
    setSteps(500);  // More steps
    setMaxiter(1000); // Higher iterations
    setQaoaLayers(2); // More QAOA layers
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
    setIsOptimizeLoading(true); // Start loading
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
    } finally {
      setIsOptimizeLoading(false); // Stop loading
    }
  };

  const handleQuantumOptimize = async () => {
    if (!moleculeData) {
      alert("Please upload a file first.");
      return;
    }
    setIsQuantumOptimizeLoading(true); // Start loading
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
    } finally {
      setIsQuantumOptimizeLoading(false); // Stop loading
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

  const handleCancelSubscription = async () => {
    const confirmCancel = window.confirm(
      "Are you sure you want to cancel your subscription? This action cannot be undone."
    );

    if (!confirmCancel) return;

    setIsCancelLoading(true); // Start loading

    try {
      const response = await axios.post(`${apiBaseUrl}/cancel-subscription`, {
        email: userEmail,
      });

      if (response.data.success) {
        alert("Your subscription has been canceled.");
        setIsSubscribed(false);
        setUserEmail(""); // Clear user email
        localStorage.removeItem("userEmail"); // Remove email from local storage
        setLimits(DEFAULT_COMPUTE_LIMITS);
      } else {
        alert("Failed to cancel subscription. Please try again.");
      }
    } catch (error) {
      console.error("Error canceling subscription:", error);
      alert("Error canceling subscription. Check the console for details.");
    } finally {
      setIsCancelLoading(false); // Stop loading
    }
  };

  const MoleculeViewer = ({ moleculeData }) => {
    const viewerRef = useRef();

    useEffect(() => {
      if (!moleculeData || !moleculeData.atoms) {
        console.warn("No molecule data found for visualization.");
        return;
      }

      const viewer = $3Dmol.createViewer(viewerRef.current, {
        backgroundColor: "grey",
      });
      viewer.clear();

      try {
        const atoms = moleculeData.atoms.map((atom) => ({
          elem: atom.element,
          x: atom.x,
          y: atom.y,
          z: atom.z,
        }));

        const model = viewer.addModel();
        model.addAtoms(atoms);

        moleculeData.atoms.forEach((atom) => {
          viewer.addLabel(atom.element, {
            position: { x: atom.x, y: atom.y, z: atom.z },
            fontSize: 12,
            fontColor: "white",
            backgroundColor: "black",
            borderRadius: 2,
            padding: 2,
            inFront: true,
          });
        });

        viewer.setStyle({}, {
          sphere: { radius: 0.2 },
          stick: { radius: 0.2 },
        });
        viewer.zoomTo();
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
          position: "relative",
          border: "1px solid #ccc",
          margin: "10px auto",
        }}
      ></div>
    );
  };

  return (
    <div style={{ padding: "20px", textAlign: "center" }}>
      <h1>Optimize Molecule</h1>
      <h2>Based on Classical and Quantum Energy Calculations</h2>
      {isSubscribed && (
        <button
          onClick={handleCancelSubscription}
          disabled={isCancelLoading}
          style={{
            position: "absolute", // Position it in the corner
            top: "10px", // Distance from the top
            right: "10px", // Distance from the right
            backgroundColor: isCancelLoading ? "#ccc" : "grey",
            color: "white",
            border: "none",
            borderRadius: "5px",
            padding: "10px 10px",
            fontSize: "10px",
            cursor: isCancelLoading ? "not-allowed" : "pointer",
            zIndex: 1000, // Ensure it's above other elements
          }}
        >
          {isCancelLoading ? "Cancelling..." : "Cancel Subscription"}
        </button>
      )}
      {!isSubscribed ? (
        <Elements stripe={stripePromise}>
          <SubscriptionForm
            onSuccess={(email) => {
              setIsSubscribeLoading(true); // Start loading
              handleSubscriptionSuccess(email);
              setIsSubscribeLoading(false); // Stop loading

            }}
          />
        </Elements>
      ) : (
        <div>
          <p>Welcome, {userEmail}! You are subscribed and have full compute power.</p>
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
            onChange={(e) =>
              handleParameterChange(setFmax, e.target.value, limits.fmax.min, limits.fmax.max)
            }
          />
        </div>
        <div>
          <label>steps:</label>
          <input
            type="number"
            value={steps}
            onChange={(e) =>
              handleParameterChange(setSteps, e.target.value, limits.steps.min, limits.steps.max)
            }
          />
        </div>
        <p>Quantum parameters:</p>
        <div>
          <label>maxiter:</label>
          <input
            type="number"
            value={maxiter}
            onChange={(e) =>
              handleParameterChange(setMaxiter, e.target.value, limits.maxiter.min, limits.maxiter.max)
            }
          />
        </div>
        <div>
          <label>QAOA Layers (p):</label>
          <input
            type="number"
            value={qaoaLayers}
            onChange={(e) =>
              handleParameterChange(setQaoaLayers, e.target.value, limits.qaoaLayers.min, limits.qaoaLayers.max)
            }
          />
        </div>
      </div>
      <div style={{ marginTop: "20px", display: "flex", justifyContent: "center", gap: "10px" }}>
        <button
          onClick={handleOptimize}
          disabled={isOptimizeLoading} // Disable button during loading
          style={{
            backgroundColor: isOptimizeLoading ? "#ccc" : "#007bff", // Grey out when loading
            color: "white",
            border: "none",
            borderRadius: "5px",
            padding: "10px 20px",
            cursor: isOptimizeLoading ? "not-allowed" : "pointer",
          }}
        >
          {isOptimizeLoading ? "Optimizing..." : "Classical Optimize"} {/* Dynamic text */}
        </button>
        <button
          onClick={handleQuantumOptimize}
          disabled={isQuantumOptimizeLoading} // Disable button during loading
          style={{
            backgroundColor: isQuantumOptimizeLoading ? "#ccc" : "#007bff", // Grey out when loading
            color: "white",
            border: "none",
            borderRadius: "5px",
            padding: "10px 20px",
            cursor: isQuantumOptimizeLoading ? "not-allowed" : "pointer",
          }}
        >
          {isQuantumOptimizeLoading ? "Quantum Optimizing..." : "Quantum Optimize"} {/* Dynamic text */}
        </button>
        <button
          onClick={fetchHowToUse}
          style={{
            position: "absolute",
            top: "10px",
            left: "10px",
            backgroundColor: "darkgreen",
            color: "white",
            border: "none",
            borderRadius: "5px",
            padding: "10px",
            fontSize: "12px",
            cursor: "pointer",
            zIndex: 1000,
          }}
        >
          How To Use & Theory
        </button>

        {isHowToUseVisible && (
          <div
            style={{
              position: "fixed",
              top: "50%",
              left: "50%",
              transform: "translate(-50%, -50%)",
              backgroundColor: "#333", // Dark grey background
              color: "white",
              border: "1px solid #444",
              boxShadow: "0 4px 8px rgba(0, 0, 0, 0.5)",
              zIndex: 1000,
              padding: "20px",
              width: "80%",
              height: "80%",
              overflowY: "auto",
              textAlign: "left", // Align text to the left
            }}
          >
            <button
              onClick={handleClosePopup}
              style={{
                position: "absolute",
                top: "10px",
                right: "10px",
                background: "red",
                color: "white",
                border: "none",
                borderRadius: "5px",
                padding: "5px 10px",
                cursor: "pointer",
              }}
            >
              Close
            </button>
            <ReactMarkdown style={{ textAlign: "left", color: "white" }}>
              {howToUseContent}
            </ReactMarkdown>
          </div>
        )}
        <button onClick={handleDownload}>Download Optimized Molecule</button>
      </div>

      {optimizedMolecule && (
        <div className="box">
          <h2>Optimized Molecule</h2>
          <pre>{JSON.stringify(optimizedMolecule.atoms, null, 2)}</pre>
          <h3>Optimization Details</h3>
          <p>
            <strong>Minimum Energy:</strong> {optimizedMolecule.min_energy?.toFixed(6)} eV
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
