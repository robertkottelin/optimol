import React, { useState, useRef, useEffect } from "react";
import axios from "axios";
import * as $3Dmol from '3dmol';
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import ReactMarkdown from "react-markdown";
import SubscriptionForm from './SubscriptionForm';

const stripePromise = loadStripe('pk_test_kbl0ETzPsoiTwU4ZJMvhsYJw006XVnV4Aq');

const App = () => {
  // State for molecule data
  const [moleculeData, setMoleculeData] = useState(null);
  const [optimizationResult, setOptimizationResult] = useState(null);
  const [activeView, setActiveView] = useState("original"); // original, classical, quantum
  
  // User and subscription state
  const [isSubscribed, setIsSubscribed] = useState(false);
  const [userEmail, setUserEmail] = useState("");
  
  // UI state
  const [isHowToUseVisible, setIsHowToUseVisible] = useState(false);
  const [howToUseContent, setHowToUseContent] = useState("");
  
  // Loading states
  const [isSubscribeLoading, setIsSubscribeLoading] = useState(false);
  const [isCancelLoading, setIsCancelLoading] = useState(false);
  const [isOptimizeLoading, setIsOptimizeLoading] = useState(false);

  const apiBaseUrl = "http://localhost:5000";

  // Check subscription status on page load
  useEffect(() => {
    const email = localStorage.getItem("userEmail");
    if (email) {
      checkSubscriptionStatus(email);
    }
  }, []);

  // Fetch "how-to-use.md" content on demand
  const fetchHowToUse = async () => {
    try {
      const response = await axios.get('/how-to-use.md');
      setHowToUseContent(response.data);
      setIsHowToUseVisible(true);
    } catch (error) {
      console.error("Error fetching 'how-to-use.md':", error);
      alert("Failed to load instructions. Check console for details.");
    }
  };

  const handleClosePopup = () => {
    setIsHowToUseVisible(false);
  };

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

  const handleSubscriptionSuccess = (email) => {
    setIsSubscribed(true);
    setUserEmail(email);
    localStorage.setItem("userEmail", email);
  };

  const handleCancelSubscription = async () => {
    const confirmCancel = window.confirm(
      "Are you sure you want to cancel your subscription? This action cannot be undone."
    );

    if (!confirmCancel) return;

    setIsCancelLoading(true);

    try {
      const response = await axios.post(`${apiBaseUrl}/cancel-subscription`, {
        email: userEmail,
      });

      if (response.data.success) {
        alert("Your subscription has been canceled.");
        setIsSubscribed(false);
        setUserEmail("");
        localStorage.removeItem("userEmail");
      } else {
        alert("Failed to cancel subscription. Please try again.");
      }
    } catch (error) {
      console.error("Error canceling subscription:", error);
      alert("Error canceling subscription. Check the console for details.");
    } finally {
      setIsCancelLoading(false);
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
        
        setMoleculeData(parsedData);
        setOptimizationResult(null);
        setActiveView("original");
      } catch (error) {
        console.error("Error processing the file:", error);
        alert("Error processing the file. Check the console for details.");
      }
    };

    reader.readAsText(file);
  };

  const validateMoleculeJSON = (data) => {
    // Check for either format: file1 structure or direct atoms array
    if (data.file1 && Array.isArray(data.file1.atoms)) {
      return data.file1.atoms.every(
        (atom) =>
          atom.id &&
          atom.element &&
          typeof atom.x === "number" &&
          typeof atom.y === "number" &&
          typeof atom.z === "number"
      );
    } else if (Array.isArray(data.atoms)) {
      return data.atoms.every(
        (atom) =>
          atom.id &&
          atom.element &&
          typeof atom.x === "number" &&
          typeof atom.y === "number" &&
          typeof atom.z === "number"
      );
    }
    return false;
  };

  const handleOptimize = async () => {
    if (!moleculeData) {
      alert("Please upload a molecule file first.");
      return;
    }
    
    if (!isSubscribed) {
      alert("Please subscribe to use optimization features.");
      return;
    }
    
    setIsOptimizeLoading(true);
    
    try {
      const payload = {
        email: userEmail,
        molecule: moleculeData,
      };
      
      const response = await axios.post(`${apiBaseUrl}/optimize-molecule`, payload);
      
      if (response.data.success) {
        setOptimizationResult(response.data);
        setActiveView("classical"); // Default to showing classical results first
      } else {
        alert("Optimization failed. " + (response.data.error || ""));
      }
    } catch (error) {
      console.error("Error optimizing molecule:", error);
      alert("Error optimizing molecule: " + (error.response?.data?.error || error.message));
    } finally {
      setIsOptimizeLoading(false);
    }
  };

  const handleDownload = (dataType) => {
    if (!optimizationResult) {
      alert("No optimization results available to download.");
      return;
    }

    let data;
    let filename;
    
    if (dataType === "classical") {
      data = {
        file1: {
          atoms: optimizationResult.classical.optimized_atoms,
          metadata: optimizationResult.classical.metadata
        }
      };
      filename = "classical_optimized_molecule.json";
    } else if (dataType === "quantum") {
      data = {
        file1: {
          atoms: optimizationResult.quantum.optimized_atoms,
          metadata: optimizationResult.quantum.metadata
        }
      };
      filename = "quantum_optimized_molecule.json";
    } else {
      alert("Invalid download type");
      return;
    }
    
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);

    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    link.click();
    
    URL.revokeObjectURL(url);
  };

  // Get atoms array based on active view
  const getAtoms = () => {
    if (activeView === "original") {
      if (moleculeData?.file1?.atoms) {
        return moleculeData.file1.atoms;
      } else if (moleculeData?.atoms) {
        return moleculeData.atoms;
      }
      return null;
    } else if (activeView === "classical" && optimizationResult?.classical?.optimized_atoms) {
      return optimizationResult.classical.optimized_atoms;
    } else if (activeView === "quantum" && optimizationResult?.quantum?.optimized_atoms) {
      return optimizationResult.quantum.optimized_atoms;
    }
    return null;
  };

  // Molecule visualization component
  const MoleculeViewer = () => {
    const viewerRef = useRef();
    const atoms = getAtoms();

    useEffect(() => {
      if (!atoms) {
        console.warn("No molecule data found for visualization.");
        return;
      }

      const viewer = $3Dmol.createViewer(viewerRef.current, {
        backgroundColor: "grey",
      });
      viewer.clear();

      try {
        // Convert to 3Dmol atom format
        const mol3dAtoms = atoms.map((atom) => ({
          elem: atom.element,
          x: atom.x,
          y: atom.y,
          z: atom.z,
        }));

        const model = viewer.addModel();
        model.addAtoms(mol3dAtoms);

        // Add element labels
        atoms.forEach((atom) => {
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

        // Style atoms and bonds
        viewer.setStyle({}, {
          sphere: { radius: 0.2 },
          stick: { radius: 0.2 },
        });
        
        viewer.zoomTo();
        viewer.render();
      } catch (error) {
        console.error("Error rendering molecule:", error);
      }
    }, [atoms]);

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
      <h1>Molecular Optimization System</h1>
      
      {/* Cancel Subscription Button */}
      {isSubscribed && (
        <button
          onClick={handleCancelSubscription}
          disabled={isCancelLoading}
          style={{
            position: "absolute",
            top: "10px",
            right: "10px",
            backgroundColor: isCancelLoading ? "#ccc" : "grey",
            color: "white",
            border: "none",
            borderRadius: "5px",
            padding: "10px 10px",
            fontSize: "10px",
            cursor: isCancelLoading ? "not-allowed" : "pointer",
            zIndex: 1000,
          }}
        >
          {isCancelLoading ? "Cancelling..." : "Cancel Subscription"}
        </button>
      )}
      
      {/* How To Use Button */}
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
      
      {/* Subscription Form or Welcome Message */}
      {!isSubscribed ? (
        <Elements stripe={stripePromise}>
          <SubscriptionForm
            onSuccess={(email) => {
              setIsSubscribeLoading(true);
              handleSubscriptionSuccess(email);
              setIsSubscribeLoading(false);
            }}
          />
        </Elements>
      ) : (
        <div>
          <p>Welcome, {userEmail}! You are subscribed and have full computational capabilities.</p>
        </div>
      )}
      
      <div>
        <h3>Molecule Optimization</h3>
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
            {/* Visualization tabs for switching between original and optimized structures */}
            {optimizationResult && (
              <div style={{ marginBottom: "10px" }}>
                <button 
                  onClick={() => setActiveView("original")}
                  style={{
                    backgroundColor: activeView === "original" ? "#007bff" : "#ccc",
                    color: "white",
                    border: "none",
                    borderRadius: "5px 0 0 5px",
                    padding: "5px 10px",
                    cursor: "pointer",
                  }}
                >
                  Original
                </button>
                <button 
                  onClick={() => setActiveView("classical")}
                  style={{
                    backgroundColor: activeView === "classical" ? "#28a745" : "#ccc",
                    color: "white",
                    border: "none",
                    padding: "5px 10px",
                    cursor: "pointer",
                  }}
                >
                  Classical
                </button>
                <button 
                  onClick={() => setActiveView("quantum")}
                  style={{
                    backgroundColor: activeView === "quantum" ? "#17a2b8" : "#ccc",
                    color: "white",
                    border: "none",
                    borderRadius: "0 5px 5px 0",
                    padding: "5px 10px",
                    cursor: "pointer",
                  }}
                >
                  Quantum
                </button>
              </div>
            )}
            
            <h3>{activeView.charAt(0).toUpperCase() + activeView.slice(1)} Structure:</h3>
            <MoleculeViewer />
          </>
        )}
        
        {/* Action Buttons */}
        <div style={{ marginTop: "20px", display: "flex", justifyContent: "center", gap: "10px" }}>
          <button
            onClick={handleOptimize}
            disabled={isOptimizeLoading || !moleculeData || !isSubscribed}
            style={{
              backgroundColor: isOptimizeLoading || !moleculeData || !isSubscribed ? "#ccc" : "#007bff",
              color: "white",
              border: "none",
              borderRadius: "5px",
              padding: "10px 20px",
              cursor: isOptimizeLoading || !moleculeData || !isSubscribed ? "not-allowed" : "pointer",
            }}
          >
            {isOptimizeLoading ? "Optimizing..." : "Optimize Molecule"}
          </button>
        </div>
        
        {/* Optimization Results */}
        {optimizationResult && (
          <div style={{ margin: "20px" }}>
            <h2>Optimization Results</h2>
            
            {/* Classical Results */}
            <div style={{ padding: "15px", border: "1px solid #ddd", borderRadius: "5px", textAlign: "left", marginBottom: "20px" }}>
              <h3>Classical Optimization</h3>
              {optimizationResult.classical.error ? (
                <p style={{ color: "red" }}>Error: {optimizationResult.classical.error}</p>
              ) : (
                <div>
                  <h4>Metadata:</h4>
                  <p><strong>Method:</strong> {optimizationResult.classical.metadata.method}</p>
                  <p><strong>Library:</strong> {optimizationResult.classical.metadata.library}</p>
                  <p><strong>Final Energy:</strong> {optimizationResult.classical.metadata.final_energy_kj_mol} kJ/mol</p>
                  <p><strong>Duration:</strong> {optimizationResult.classical.metadata.duration_seconds} seconds</p>
                  <p><strong>Status:</strong> {optimizationResult.classical.metadata.convergence || "Unknown"}</p>
                  
                  <button
                    onClick={() => handleDownload("classical")}
                    style={{
                      backgroundColor: "#28a745",
                      color: "white",
                      border: "none",
                      borderRadius: "5px",
                      padding: "5px 10px",
                      cursor: "pointer",
                    }}
                  >
                    Download Classical Result
                  </button>
                </div>
              )}
            </div>
            
            {/* Quantum Results */}
            <div style={{ padding: "15px", border: "1px solid #ddd", borderRadius: "5px", textAlign: "left" }}>
              <h3>Quantum Optimization</h3>
              {optimizationResult.quantum.error ? (
                <p style={{ color: "red" }}>Error: {optimizationResult.quantum.error}</p>
              ) : (
                <div>
                  <h4>Metadata:</h4>
                  <p><strong>Method:</strong> {optimizationResult.quantum.metadata.method}</p>
                  <p><strong>Library:</strong> {optimizationResult.quantum.metadata.library}</p>
                  <p><strong>Theory Level:</strong> {optimizationResult.quantum.metadata.theory_level}</p>
                  <p><strong>Final Energy:</strong> {optimizationResult.quantum.metadata.final_energy_hartree} Hartree</p>
                  <p><strong>Iterations:</strong> {optimizationResult.quantum.metadata.iterations}</p>
                  <p><strong>Converged:</strong> {optimizationResult.quantum.metadata.converged ? "Yes" : "No"}</p>
                  <p><strong>Duration:</strong> {optimizationResult.quantum.metadata.duration_seconds} seconds</p>
                  
                  <button
                    onClick={() => handleDownload("quantum")}
                    style={{
                      backgroundColor: "#17a2b8",
                      color: "white",
                      border: "none",
                      borderRadius: "5px",
                      padding: "5px 10px",
                      cursor: "pointer",
                    }}
                  >
                    Download Quantum Result
                  </button>
                </div>
              )}
            </div>
          </div>
        )}
      </div>
            
      {/* How To Use Popup */}
      {isHowToUseVisible && (
        <div
          style={{
            position: "fixed",
            top: "50%",
            left: "50%",
            transform: "translate(-50%, -50%)",
            backgroundColor: "#333",
            color: "white",
            border: "1px solid #444",
            boxShadow: "0 4px 8px rgba(0, 0, 0, 0.5)",
            zIndex: 1000,
            padding: "20px",
            width: "80%",
            height: "80%",
            overflowY: "auto",
            textAlign: "left",
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
    </div>
  );
};

export default App;