import React, { useState, useRef, useEffect } from "react";
import axios from "axios";
import * as $3Dmol from '3dmol';
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import ReactMarkdown from "react-markdown";
import SubscriptionForm from './SubscriptionForm';

const stripePromise = loadStripe('pk_test_kbl0ETzPsoiTwU4ZJMvhsYJw006XVnV4Aq');

// Computation limits based on subscription status
const DEFAULT_COMPUTE_LIMITS = {
};

const SUBSCRIBED_COMPUTE_LIMITS = {
};

const App = () => {
  // State for molecule data
  const [moleculeData, setMoleculeData] = useState(null);
  const [optimizedMolecule, setOptimizedMolecule] = useState(null);
  const [proteinData, setProteinData] = useState(null);
  const [ligandData, setLigandData] = useState(null);
  const [bindingResult, setBindingResult] = useState(null);
  
  // User and subscription state
  const [isSubscribed, setIsSubscribed] = useState(false);
  const [userEmail, setUserEmail] = useState("");
  const [limits, setLimits] = useState(DEFAULT_COMPUTE_LIMITS);
  
  // UI state
  const [isHowToUseVisible, setIsHowToUseVisible] = useState(false);
  const [howToUseContent, setHowToUseContent] = useState("");
  const [activeTab, setActiveTab] = useState("molecule");
  
  // Loading states
  const [isSubscribeLoading, setIsSubscribeLoading] = useState(false);
  const [isCancelLoading, setIsCancelLoading] = useState(false);
  const [isOptimizeLoading, setIsOptimizeLoading] = useState(false);
  const [isBindingOptimizeLoading, setIsBindingOptimizeLoading] = useState(false);

  const apiBaseUrl = "http://localhost:5000/";

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
        setLimits(SUBSCRIBED_COMPUTE_LIMITS);
        
        // Update parameters to use full capabilities

      } else {
        setIsSubscribed(false);
        setLimits(DEFAULT_COMPUTE_LIMITS);
        
        // Reset to basic parameters

      }
    } catch (error) {
      console.error("Error checking subscription status:", error);
    }
  };

  const handleSubscriptionSuccess = (email) => {
    setIsSubscribed(true);
    setUserEmail(email);
    localStorage.setItem("userEmail", email);
    setLimits(SUBSCRIBED_COMPUTE_LIMITS);
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
        setLimits(DEFAULT_COMPUTE_LIMITS);
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

  const handleParameterChange = (setter, value, min, max) => {
    const parsedValue = parseFloat(value);
    if (parsedValue >= min && parsedValue <= max) {
      setter(parsedValue);
    }
  };

  const handleFileUpload = (event, fileType) => {
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

        if (fileType === "molecule") {
          if (!validateMoleculeJSON(parsedData)) {
            alert("Invalid molecule JSON format.");
            return;
          }
          setMoleculeData(parsedData.file1);
          setOptimizedMolecule(null);
        } else if (fileType === "protein") {
          if (!validateMoleculeJSON(parsedData)) {
            alert("Invalid protein JSON format.");
            return;
          }
          setProteinData(parsedData.file1);
          setBindingResult(null);
        } else if (fileType === "ligand") {
          if (!validateMoleculeJSON(parsedData)) {
            alert("Invalid ligand JSON format.");
            return;
          }
          setLigandData(parsedData.file1);
          setBindingResult(null);
        }
      } catch (error) {
        console.error("Error processing the file:", error);
        alert("Error processing the file. Check the console for details.");
      }
    };

    reader.readAsText(file);
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

  const handleOptimize = async () => {
    if (!moleculeData) {
      alert("Please upload a molecule file first.");
      return;
    }
    
    setIsOptimizeLoading(true);
    
    try {
      const payload = {
        file1: moleculeData,
      };
      
      // Add authentication for subscribed users
      if (isSubscribed) {
        payload.email = userEmail;
      }
      
      const response = await axios.post(`${apiBaseUrl}/optimize`, payload);
      setOptimizedMolecule(response.data.optimized_file1);
    } catch (error) {
      console.error("Error optimizing molecule:", error);
      alert("Error optimizing molecule. Check the console for details.");
    } finally {
      setIsOptimizeLoading(false);
    }
  };

  const handleBindingOptimize = async () => {
    if (!proteinData || !ligandData) {
      alert("Please upload both protein and ligand files first.");
      return;
    }
    
    if (!isSubscribed) {
      alert("Binding optimization requires a subscription.");
      return;
    }
    
    setIsBindingOptimizeLoading(true);
    
    try {
      const payload = {
        protein: proteinData,
        ligand: ligandData,
        email: userEmail,
      };
      
      const response = await axios.post(`${apiBaseUrl}/binding-optimize`, payload);
      setBindingResult(response.data.binding_result);
    } catch (error) {
      console.error("Error optimizing binding:", error);
      alert("Error optimizing binding. Check the console for details.");
    } finally {
      setIsBindingOptimizeLoading(false);
    }
  };

  const handleDownload = (data, filename) => {
    if (!data) {
      alert("No data available to download.");
      return;
    }

    const jsonData = {
      file1: data,
    };
    
    const blob = new Blob([JSON.stringify(jsonData, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);

    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    link.click();
    
    URL.revokeObjectURL(url);
  };

  // Molecule visualization component
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

        // Add element labels
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

  // Binding system visualization component
  const BindingViewer = ({ proteinData, ligandData }) => {
    const viewerRef = useRef();

    useEffect(() => {
      if (!proteinData || !proteinData.atoms || !ligandData || !ligandData.atoms) {
        console.warn("Missing protein or ligand data for visualization.");
        return;
      }

      const viewer = $3Dmol.createViewer(viewerRef.current, {
        backgroundColor: "grey",
      });
      viewer.clear();

      try {
        // Add protein model
        const proteinModel = viewer.addModel();
        const proteinAtoms = proteinData.atoms.map((atom) => ({
          elem: atom.element,
          x: atom.x,
          y: atom.y,
          z: atom.z,
          properties: { protein: true }
        }));
        proteinModel.addAtoms(proteinAtoms);
        
        // Add ligand model
        const ligandModel = viewer.addModel();
        const ligandAtoms = ligandData.atoms.map((atom) => ({
          elem: atom.element,
          x: atom.x,
          y: atom.y,
          z: atom.z,
          properties: { ligand: true }
        }));
        ligandModel.addAtoms(ligandAtoms);

        // Style protein
        viewer.setStyle({properties: {protein: true}}, {
          cartoon: { color: 'blue' },
          stick: { radius: 0.1, color: 'lightblue' }
        });
        
        // Style ligand
        viewer.setStyle({properties: {ligand: true}}, {
          sphere: { radius: 0.4 },
          stick: { radius: 0.2, color: 'green' }
        });
        
        viewer.zoomTo();
        viewer.render();
      } catch (error) {
        console.error("Error rendering binding system:", error);
      }
    }, [proteinData, ligandData]);

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
      <h1>Molecular Simulation System</h1>
      
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
      
      {/* Tab Navigation */}
      <div style={{ margin: "20px 0", borderBottom: "1px solid #ccc" }}>
        <button 
          style={{
            padding: "10px 20px",
            backgroundColor: activeTab === "molecule" ? "#007bff" : "#f0f0f0",
            color: activeTab === "molecule" ? "white" : "black",
            border: "none",
            borderRadius: "5px 5px 0 0",
            marginRight: "5px"
          }}
          onClick={() => setActiveTab("molecule")}
        >
          Single Molecule
        </button>
        <button 
          style={{
            padding: "10px 20px",
            backgroundColor: activeTab === "binding" ? "#007bff" : "#f0f0f0",
            color: activeTab === "binding" ? "white" : "black",
            border: "none",
            borderRadius: "5px 5px 0 0",
            marginRight: "5px"
          }}
          onClick={() => setActiveTab("binding")}
        >
          Protein-Ligand Binding
        </button>
      </div>
      
      {/* Single Molecule Tab */}
      {activeTab === "molecule" && (
        <div>
          <h3>Single Molecule Optimization</h3>
          <input
            type="file"
            onChange={(e) => handleFileUpload(e, "molecule")}
            style={{
              padding: "10px",
              border: "2px dashed gray",
              cursor: "pointer",
              marginBottom: "20px",
            }}
          />
          
          {moleculeData && (
            <>
              <h3>Molecule Visualization:</h3>
              <MoleculeViewer moleculeData={optimizedMolecule || moleculeData} />
            </>
          )}
          
          
          {/* Action Buttons */}
          <div style={{ marginTop: "20px", display: "flex", justifyContent: "center", gap: "10px" }}>
            <button
              onClick={handleOptimize}
              disabled={isOptimizeLoading || !moleculeData}
              style={{
                backgroundColor: isOptimizeLoading || !moleculeData ? "#ccc" : "#007bff",
                color: "white",
                border: "none",
                borderRadius: "5px",
                padding: "10px 20px",
                cursor: isOptimizeLoading || !moleculeData ? "not-allowed" : "pointer",
              }}
            >
              {isOptimizeLoading ? "Optimizing..." : "Classical Optimize"}
            </button>
            <button
              onClick={() => handleDownload(optimizedMolecule, "optimized_molecule.json")}
              disabled={!optimizedMolecule}
              style={{
                backgroundColor: !optimizedMolecule ? "#ccc" : "#28a745",
                color: "white",
                border: "none",
                borderRadius: "5px",
                padding: "10px 20px",
                cursor: !optimizedMolecule ? "not-allowed" : "pointer",
              }}
            >
              Download Result
            </button>
          </div>
          
          {/* Optimization Results */}
          {optimizedMolecule && (
            <div style={{ margin: "20px", padding: "15px", border: "1px solid #ddd", borderRadius: "5px", textAlign: "left" }}>
              <h2>Optimization Results</h2>
              <p><strong>Energy:</strong> {optimizedMolecule.energy?.toFixed(6)} eV</p>
              <p><strong>Converged:</strong> {optimizedMolecule.converged ? "Yes" : "No"}</p>
              <h3>Parameters Used:</h3>
              <ul>

              </ul>
            </div>
          )}
        </div>
      )}
      
      {/* Protein-Ligand Binding Tab */}
      {activeTab === "binding" && (
        <div>
          <h3>Protein-Ligand Binding Optimization</h3>
          
          {!isSubscribed && (
            <div style={{ color: "red", margin: "10px 0" }}>
              Protein-ligand binding optimization requires subscription.
            </div>
          )}
          
          <div style={{ display: "flex", justifyContent: "center", gap: "20px", marginBottom: "20px" }}>
            <div>
              <h4>Protein File</h4>
              <input
                type="file"
                onChange={(e) => handleFileUpload(e, "protein")}
                style={{
                  padding: "10px",
                  border: "2px dashed gray",
                  cursor: "pointer",
                }}
                disabled={!isSubscribed}
              />
            </div>
            <div>
              <h4>Ligand File</h4>
              <input
                type="file"
                onChange={(e) => handleFileUpload(e, "ligand")}
                style={{
                  padding: "10px",
                  border: "2px dashed gray",
                  cursor: "pointer",
                }}
                disabled={!isSubscribed}
              />
            </div>
          </div>
          
          {/* Binding System Visualization */}
          {proteinData && ligandData && (
            <>
              <h3>Binding System Visualization:</h3>
              <BindingViewer 
                proteinData={bindingResult?.protein || proteinData} 
                ligandData={bindingResult?.ligand || ligandData} 
              />
            </>
          )}
          
          {/* Binding Parameters */}
          <div style={{ display: "flex", gap: "20px", justifyContent: "center", flexWrap: "wrap", marginTop: "20px" }}>
            <div>
              <label>Force Convergence (fmax):</label>
              <input
                type="number"
                step="0.001"
                value={fmax}
                onChange={(e) =>
                  handleParameterChange(setFmax, e.target.value, limits.fmax.min, limits.fmax.max)
                }
                style={{ marginLeft: "5px" }}
                disabled={!isSubscribed}
              />
            </div>
            <div>
              <label>Optimization Steps:</label>
              <input
                type="number"
                value={steps}
                onChange={(e) =>
                  handleParameterChange(setSteps, e.target.value, limits.steps.min, limits.steps.max)
                }
                style={{ marginLeft: "5px" }}
                disabled={!isSubscribed}
              />
            </div>
          </div>
          
          {/* Binding Action Buttons */}
          <div style={{ marginTop: "20px", display: "flex", justifyContent: "center", gap: "10px" }}>
            <button
              onClick={handleBindingOptimize}
              disabled={isBindingOptimizeLoading || !proteinData || !ligandData || !isSubscribed}
              style={{
                backgroundColor: isBindingOptimizeLoading || !proteinData || !ligandData || !isSubscribed ? "#ccc" : "#007bff",
                color: "white",
                border: "none",
                borderRadius: "5px",
                padding: "10px 20px",
                cursor: isBindingOptimizeLoading || !proteinData || !ligandData || !isSubscribed ? "not-allowed" : "pointer",
              }}
            >
              {isBindingOptimizeLoading ? "Optimizing Binding..." : "Optimize Binding"}
            </button>
            <button
              onClick={() => {
                if (bindingResult) {
                  handleDownload({
                    protein: bindingResult.protein,
                    ligand: bindingResult.ligand,
                    binding_energy: bindingResult.binding_energy,
                    converged: bindingResult.converged
                  }, "binding_result.json");
                }
              }}
              disabled={!bindingResult}
              style={{
                backgroundColor: !bindingResult ? "#ccc" : "#28a745",
                color: "white",
                border: "none",
                borderRadius: "5px",
                padding: "10px 20px",
                cursor: !bindingResult ? "not-allowed" : "pointer",
              }}
            >
              Download Binding Result
            </button>
          </div>
          
          {/* Binding Results */}
          {bindingResult && (
            <div style={{ margin: "20px", padding: "15px", border: "1px solid #ddd", borderRadius: "5px", textAlign: "left" }}>
              <h2>Binding Optimization Results</h2>
              <p><strong>Binding Energy:</strong> {bindingResult.binding_energy?.toFixed(6)} eV</p>
              <p><strong>Converged:</strong> {bindingResult.converged ? "Yes" : "No"}</p>
            </div>
          )}
        </div>
      )}
      
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