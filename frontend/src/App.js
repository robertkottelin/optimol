import React, { useState, useRef, useEffect } from "react";
import axios from "axios";
import * as $3Dmol from '3dmol';
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import ReactMarkdown from "react-markdown";
import SubscriptionForm from './SubscriptionForm';

const stripePromise = loadStripe('pk_test_kbl0ETzPsoiTwU4ZJMvhsYJw006XVnV4Aq');

// Constants for iteration limits
const ITERATION_LIMITS = {
  subscribed: {
    classical: 100000,
    quantum: 1000
  },
  unsubscribed: {
    classical: 100, // Reduced limit for non-subscribers
    quantum: 5      // Reduced limit for non-subscribers
  }
};

// Default optimization parameters
const defaultClassicalParams = {
  temperature: 300,
  max_iterations: 1000,  // Will be capped for non-subscribers
  bond_threshold: 0.2,
  bond_force_constant: 1000.0,
  angle_force_constant: 500.0
};

const defaultQuantumParams = {
  basis: "6-31g",
  max_iterations: 10,   // Will be capped for non-subscribers
  convergence_threshold: 0.00001,
  step_size: 0.1
};

const App = () => {
  // State for molecule data
  const [moleculeData, setMoleculeData] = useState(null);
  const [optimizationResult, setOptimizationResult] = useState(null);
  const [activeView, setActiveView] = useState("original"); // original, optimized
  const [optimizationType, setOptimizationType] = useState("classical"); // classical or quantum
  
  // Optimization parameters state
  const [classicalParams, setClassicalParams] = useState({...defaultClassicalParams});
  const [quantumParams, setQuantumParams] = useState({...defaultQuantumParams});
  const [showAdvancedParams, setShowAdvancedParams] = useState(false);
  
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

  // Helper function to consistently apply iteration limits
  const applyIterationLimits = (isUserSubscribed) => {
    // Apply limits to classical parameters
    const limitedClassicalParams = {...defaultClassicalParams};
    const classicalMaxIterations = isUserSubscribed 
      ? ITERATION_LIMITS.subscribed.classical
      : ITERATION_LIMITS.unsubscribed.classical;
      
    limitedClassicalParams.max_iterations = Math.min(
      limitedClassicalParams.max_iterations, 
      classicalMaxIterations
    );
    setClassicalParams(limitedClassicalParams);
    
    // Apply limits to quantum parameters
    const limitedQuantumParams = {...defaultQuantumParams};
    const quantumMaxIterations = isUserSubscribed 
      ? ITERATION_LIMITS.subscribed.quantum
      : ITERATION_LIMITS.unsubscribed.quantum;
      
    limitedQuantumParams.max_iterations = Math.min(
      limitedQuantumParams.max_iterations, 
      quantumMaxIterations
    );
    
    // Apply basis set restrictions for non-subscribers
    if (!isUserSubscribed && (limitedQuantumParams.basis === "6-311g" || limitedQuantumParams.basis === "cc-pvdz")) {
      limitedQuantumParams.basis = "6-31g";
    }
    
    setQuantumParams(limitedQuantumParams);
  };

  // Check subscription status on page load and apply limits
  useEffect(() => {
    const email = localStorage.getItem("userEmail");
    if (email) {
      checkSubscriptionStatus(email);
    } else {
      // For non-subscribed users, apply the default limits
      applyIterationLimits(false);
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
      const userIsSubscribed = response.data.isSubscribed;
      
      setIsSubscribed(userIsSubscribed);
      setUserEmail(email);
      
      // Apply appropriate limits based on subscription status
      applyIterationLimits(userIsSubscribed);
    } catch (error) {
      console.error("Error checking subscription status:", error);
      // If there's an error, assume user is not subscribed
      setIsSubscribed(false);
      applyIterationLimits(false);
    }
  };

  const handleSubscriptionSuccess = (email) => {
    setIsSubscribed(true);
    setUserEmail(email);
    localStorage.setItem("userEmail", email);
    
    // Apply subscriber limits
    applyIterationLimits(true);
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
        
        // Apply non-subscriber limits using the helper function
        applyIterationLimits(false);
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
    
    setIsOptimizeLoading(true);
    
    try {
      // Get the correct parameters based on selected optimization type
      const optimizationParams = 
        optimizationType === "classical" ? {...classicalParams} : {...quantumParams};
        
      // Apply iteration limits for all users
      if (optimizationType === "classical") {
        // Apply appropriate limits based on subscription status
        const maxIterations = isSubscribed 
          ? ITERATION_LIMITS.subscribed.classical
          : ITERATION_LIMITS.unsubscribed.classical;
          
        optimizationParams.max_iterations = Math.min(
          optimizationParams.max_iterations, 
          maxIterations
        );
      } else {
        // Apply appropriate limits based on subscription status
        const maxIterations = isSubscribed 
          ? ITERATION_LIMITS.subscribed.quantum
          : ITERATION_LIMITS.unsubscribed.quantum;
          
        optimizationParams.max_iterations = Math.min(
          optimizationParams.max_iterations, 
          maxIterations
        );
        
        // Enforce basis set restrictions for non-subscribers only
        if (!isSubscribed && (optimizationParams.basis === "6-311g" || optimizationParams.basis === "cc-pvdz")) {
          optimizationParams.basis = "6-31g";
        }
      }
      
      const payload = {
        email: userEmail || "guest@example.com",
        molecule: moleculeData,
        optimization_type: optimizationType,
        optimization_params: optimizationParams
      };
      
      const response = await axios.post(`${apiBaseUrl}/optimize-molecule`, payload);
      
      if (response.data.success) {
        setOptimizationResult(response.data);
        setActiveView("optimized"); // Show optimized results after successful optimization
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

  const handleParamChange = (type, paramName, value) => {
    if (type === "classical") {
      setClassicalParams(prev => ({
        ...prev,
        [paramName]: value
      }));
    } else {
      setQuantumParams(prev => ({
        ...prev,
        [paramName]: value
      }));
    }
  };

  const handleResetParams = (type) => {
    if (type === "classical") {
      const params = {...defaultClassicalParams};
      
      // Always apply appropriate limits based on subscription status
      const maxIterations = isSubscribed 
        ? ITERATION_LIMITS.subscribed.classical
        : ITERATION_LIMITS.unsubscribed.classical;
        
      params.max_iterations = Math.min(
        params.max_iterations, 
        maxIterations
      );
      
      setClassicalParams(params);
    } else {
      const params = {...defaultQuantumParams};
      
      // Always apply appropriate limits based on subscription status
      const maxIterations = isSubscribed 
        ? ITERATION_LIMITS.subscribed.quantum
        : ITERATION_LIMITS.unsubscribed.quantum;
        
      params.max_iterations = Math.min(
        params.max_iterations, 
        maxIterations
      );
      
      // Restrict to simpler basis sets for free users only
      if (!isSubscribed && (params.basis === "6-311g" || params.basis === "cc-pvdz")) {
        params.basis = "6-31g";
      }
      
      setQuantumParams(params);
    }
  };

  const handleDownload = () => {
    if (!optimizationResult) {
      alert("No optimization results available to download.");
      return;
    }

    const data = {
      file1: {
        atoms: optimizationResult.result.optimized_atoms,
        metadata: optimizationResult.result.metadata
      }
    };
    
    const filename = `${optimizationType}_optimized_molecule.json`;
    
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
    } else if (activeView === "optimized" && optimizationResult?.result?.optimized_atoms) {
      return optimizationResult.result.optimized_atoms;
    }
    return null;
  };

  // Subscription Limit Notice Component
  const SubscriptionLimitNotice = ({ isSubscribed, optimizationType }) => {
    if (isSubscribed) return null;
    
    const limit = optimizationType === "classical" 
      ? ITERATION_LIMITS.unsubscribed.classical 
      : ITERATION_LIMITS.unsubscribed.quantum;
    
    const fullLimit = optimizationType === "classical" 
      ? ITERATION_LIMITS.subscribed.classical 
      : ITERATION_LIMITS.subscribed.quantum;
    
    return (
      <div style={{
        backgroundColor: '#772200',
        color: 'white',
        padding: '8px 12px',
        borderRadius: '4px',
        marginBottom: '15px',
        fontSize: '14px',
        textAlign: 'left'
      }}>
        <strong>Free Account Limitation:</strong> Iterations capped at {limit} (vs. {fullLimit} for subscribers). 
        <a 
          href="#" 
          onClick={(e) => { e.preventDefault(); document.querySelector('.subscription-form').scrollIntoView(); }}
          style={{ color: '#ffcc00', marginLeft: '5px' }}
        >
          Subscribe for full capabilities
        </a>
      </div>
    );
  };

  // OptimizeButton Component
  const OptimizeButton = () => {
    // Determine the button text
    let buttonText = isOptimizeLoading 
      ? "Optimizing..." 
      : `Run ${optimizationType === "classical" ? "Classical" : "Quantum"} Optimization`;
      
    // Add subscription indicator for free users
    if (!isSubscribed && !isOptimizeLoading) {
      buttonText += " (Limited)";
    }
    
    return (
      <div style={{ marginTop: "20px", display: "flex", flexDirection: "column", alignItems: "center", gap: "10px" }}>
        <button
          onClick={handleOptimize}
          disabled={isOptimizeLoading || !moleculeData}
          style={{
            backgroundColor: isOptimizeLoading || !moleculeData ? "#ccc" : 
                            (optimizationType === "classical" ? "#28a745" : "#17a2b8"),
            color: "white",
            border: "none",
            borderRadius: "5px",
            padding: "10px 20px",
            cursor: isOptimizeLoading || !moleculeData ? "not-allowed" : "pointer",
          }}
        >
          {buttonText}
        </button>
        
        {!isSubscribed && (
          <div style={{ 
            fontSize: "12px", 
            color: "#aaa", 
            maxWidth: "400px",
            textAlign: "center"
          }}>
            Free users are limited to {ITERATION_LIMITS.unsubscribed.classical} iterations for classical and {ITERATION_LIMITS.unsubscribed.quantum} for quantum optimizations.
          </div>
        )}
      </div>
    );
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

  // Parameter Configuration Components
  const ClassicalParametersConfig = () => (
    <div style={{ 
      border: '1px solid #28a745', 
      padding: '15px', 
      borderRadius: '5px', 
      marginTop: '10px', 
      backgroundColor: '#242424',
      color: '#ffffff'
    }}>
      <h4 style={{ color: '#28a745', marginBottom: '15px' }}>Classical Optimization Parameters</h4>
      
      <SubscriptionLimitNotice isSubscribed={isSubscribed} optimizationType="classical" />
      
      <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
        <div>
          <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '180px' }}>Temperature (K): </label>
          <input 
            type="number" 
            value={classicalParams.temperature}
            min="1"
            max="1000"
            step="10"
            onChange={(e) => handleParamChange('classical', 'temperature', Number(e.target.value))}
            style={{ 
              width: '80px', 
              marginLeft: '5px',
              backgroundColor: '#444',
              color: '#ffffff',
              padding: '6px',
              borderRadius: '4px',
              border: '1px solid #28a745'
            }}
          />
        </div>
        
        <div>
          <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '180px' }}>Max Iterations: </label>
          <input 
            type="number" 
            value={classicalParams.max_iterations}
            min="100"
            max={isSubscribed ? ITERATION_LIMITS.subscribed.classical : ITERATION_LIMITS.unsubscribed.classical}
            step="100"
            onChange={(e) => handleParamChange('classical', 'max_iterations', Number(e.target.value))}
            style={{ 
              width: '80px', 
              marginLeft: '5px',
              backgroundColor: '#444',
              color: '#ffffff',
              padding: '6px',
              borderRadius: '4px',
              border: '1px solid #28a745'
            }}
          />
          {!isSubscribed && classicalParams.max_iterations > ITERATION_LIMITS.unsubscribed.classical && (
            <span style={{ color: '#ffcc00', marginLeft: '10px', fontSize: '12px' }}>
              Will be capped at {ITERATION_LIMITS.unsubscribed.classical}
            </span>
          )}
          {isSubscribed && classicalParams.max_iterations > ITERATION_LIMITS.subscribed.classical && (
            <span style={{ color: '#ffcc00', marginLeft: '10px', fontSize: '12px' }}>
              Will be capped at {ITERATION_LIMITS.subscribed.classical}
            </span>
          )}
        </div>
        
        {showAdvancedParams && (
          <>
            <div>
              <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '180px' }}>Bond Threshold (nm): </label>
              <input 
                type="number" 
                value={classicalParams.bond_threshold}
                min="0.1"
                max="0.5"
                step="0.01"
                onChange={(e) => handleParamChange('classical', 'bond_threshold', Number(e.target.value))}
                style={{ 
                  width: '80px', 
                  marginLeft: '5px',
                  backgroundColor: '#444',
                  color: '#ffffff',
                  padding: '6px',
                  borderRadius: '4px',
                  border: '1px solid #28a745'
                }}
              />
            </div>
            
            <div>
              <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '220px' }}>Bond Force Constant (kJ/mol/nm²): </label>
              <input 
                type="number" 
                value={classicalParams.bond_force_constant}
                min="100"
                max="10000"
                step="100"
                onChange={(e) => handleParamChange('classical', 'bond_force_constant', Number(e.target.value))}
                style={{ 
                  width: '80px', 
                  marginLeft: '5px',
                  backgroundColor: '#444',
                  color: '#ffffff',
                  padding: '6px',
                  borderRadius: '4px',
                  border: '1px solid #28a745'
                }}
              />
            </div>
            
            <div>
              <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '220px' }}>Angle Force Constant (kJ/mol/rad²): </label>
              <input 
                type="number" 
                value={classicalParams.angle_force_constant}
                min="50"
                max="5000"
                step="50"
                onChange={(e) => handleParamChange('classical', 'angle_force_constant', Number(e.target.value))}
                style={{ 
                  width: '80px', 
                  marginLeft: '5px',
                  backgroundColor: '#444',
                  color: '#ffffff',
                  padding: '6px',
                  borderRadius: '4px',
                  border: '1px solid #28a745'
                }}
              />
            </div>
          </>
        )}
      </div>
      
      <div style={{ marginTop: '15px' }}>
        <button
          onClick={() => handleResetParams('classical')}
          style={{
            backgroundColor: '#6c757d',
            color: 'white',
            border: 'none',
            borderRadius: '3px',
            padding: '8px 12px',
            marginRight: '10px',
            cursor: 'pointer',
          }}
        >
          Reset to Defaults
        </button>
        
        <button
          onClick={() => setShowAdvancedParams(!showAdvancedParams)}
          style={{
            backgroundColor: '#17a2b8',
            color: 'white',
            border: 'none',
            borderRadius: '3px',
            padding: '8px 12px',
            cursor: 'pointer',
          }}
        >
          {showAdvancedParams ? 'Hide Advanced Parameters' : 'Show Advanced Parameters'}
        </button>
      </div>
    </div>
  );

  const QuantumParametersConfig = () => (
    <div style={{ 
      border: '1px solid #33b5e5', 
      padding: '15px', 
      borderRadius: '5px', 
      marginTop: '10px', 
      backgroundColor: '#242424',
      color: '#ffffff'
    }}>
      <h4 style={{ color: '#33b5e5', marginBottom: '15px' }}>Quantum Optimization Parameters</h4>
      
      <SubscriptionLimitNotice isSubscribed={isSubscribed} optimizationType="quantum" />
      
      <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
        <div>
          <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '180px' }}>Basis Set: </label>
          <select 
            value={quantumParams.basis}
            onChange={(e) => handleParamChange('quantum', 'basis', e.target.value)}
            style={{ 
              marginLeft: '5px', 
              backgroundColor: '#444',
              color: '#ffffff',
              padding: '6px',
              borderRadius: '4px',
              border: '1px solid #33b5e5'
            }}
          >
            <option value="sto-3g">STO-3G (Minimal)</option>
            <option value="6-31g">6-31G (Standard)</option>
            {isSubscribed && <option value="6-311g">6-311G (Extended)</option>}
            {isSubscribed && <option value="cc-pvdz">cc-pVDZ (Double Zeta)</option>}
          </select>
          {!isSubscribed && (quantumParams.basis === "6-311g" || quantumParams.basis === "cc-pvdz") && (
            <span style={{ color: '#ffcc00', marginLeft: '10px', fontSize: '12px' }}>
              Extended basis sets require subscription
            </span>
          )}
        </div>
        
        <div>
          <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '180px' }}>Max Iterations: </label>
          <input 
            type="number" 
            value={quantumParams.max_iterations}
            min="1"
            max={isSubscribed ? ITERATION_LIMITS.subscribed.quantum : ITERATION_LIMITS.unsubscribed.quantum}
            onChange={(e) => handleParamChange('quantum', 'max_iterations', Number(e.target.value))}
            style={{ 
              width: '80px', 
              marginLeft: '5px',
              backgroundColor: '#444',
              color: '#ffffff',
              padding: '6px',
              borderRadius: '4px',
              border: '1px solid #33b5e5'
            }}
          />
          {!isSubscribed && quantumParams.max_iterations > ITERATION_LIMITS.unsubscribed.quantum && (
            <span style={{ color: '#ffcc00', marginLeft: '10px', fontSize: '12px' }}>
              Will be capped at {ITERATION_LIMITS.unsubscribed.quantum}
            </span>
          )}
          {isSubscribed && quantumParams.max_iterations > ITERATION_LIMITS.subscribed.quantum && (
            <span style={{ color: '#ffcc00', marginLeft: '10px', fontSize: '12px' }}>
              Will be capped at {ITERATION_LIMITS.subscribed.quantum}
            </span>
          )}
        </div>
        
        {showAdvancedParams && (
          <>
            <div>
              <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '180px' }}>Convergence Threshold: </label>
              <input 
                type="number" 
                value={quantumParams.convergence_threshold}
                min="0.000001"
                max="0.01"
                step="0.000001"
                onChange={(e) => handleParamChange('quantum', 'convergence_threshold', Number(e.target.value))}
                style={{ 
                  width: '120px', 
                  marginLeft: '5px',
                  backgroundColor: '#444',
                  color: '#ffffff',
                  padding: '6px',
                  borderRadius: '4px',
                  border: '1px solid #33b5e5'
                }}
              />
            </div>
            
            <div>
              <label style={{ color: '#ffffff', fontWeight: 'bold', display: 'inline-block', width: '180px' }}>Step Size: </label>
              <input 
                type="number" 
                value={quantumParams.step_size}
                min="0.01"
                max="1.0"
                step="0.01"
                onChange={(e) => handleParamChange('quantum', 'step_size', Number(e.target.value))}
                style={{ 
                  width: '80px', 
                  marginLeft: '5px',
                  backgroundColor: '#444',
                  color: '#ffffff',
                  padding: '6px',
                  borderRadius: '4px',
                  border: '1px solid #33b5e5'
                }}
              />
            </div>
          </>
        )}
      </div>
      
      <div style={{ marginTop: '15px' }}>
        <button
          onClick={() => handleResetParams('quantum')}
          style={{
            backgroundColor: '#6c757d',
            color: 'white',
            border: 'none',
            borderRadius: '3px',
            padding: '8px 12px',
            marginRight: '10px',
            cursor: 'pointer',
          }}
        >
          Reset to Defaults
        </button>
        
        <button
          onClick={() => setShowAdvancedParams(!showAdvancedParams)}
          style={{
            backgroundColor: '#17a2b8',
            color: 'white',
            border: 'none',
            borderRadius: '3px',
            padding: '8px 12px',
            cursor: 'pointer',
          }}
        >
          {showAdvancedParams ? 'Hide Advanced Parameters' : 'Show Advanced Parameters'}
        </button>
      </div>
    </div>
  );

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
        <div className="subscription-form">
          <Elements stripe={stripePromise}>
            <SubscriptionForm
              onSuccess={(email) => {
                setIsSubscribeLoading(true);
                handleSubscriptionSuccess(email);
                setIsSubscribeLoading(false);
              }}
            />
          </Elements>
        </div>
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
            {/* Optimization Type Selection */}
            <div style={{ marginBottom: "20px" }}>
              <h4>Select Optimization Method:</h4>
              <div style={{ display: "flex", justifyContent: "center", gap: "10px" }}>
                <button
                  onClick={() => setOptimizationType("classical")}
                  style={{
                    backgroundColor: optimizationType === "classical" ? "#28a745" : "#ccc",
                    color: "white",
                    border: "none",
                    borderRadius: "5px",
                    padding: "10px 20px",
                    cursor: "pointer",
                  }}
                >
                  Classical Optimization
                </button>
                <button
                  onClick={() => setOptimizationType("quantum")}
                  style={{
                    backgroundColor: optimizationType === "quantum" ? "#17a2b8" : "#ccc",
                    color: "white",
                    border: "none",
                    borderRadius: "5px",
                    padding: "10px 20px",
                    cursor: "pointer",
                  }}
                >
                  Quantum Optimization
                </button>
              </div>
            </div>
            
            {/* Parameter Configuration */}
            {optimizationType === "classical" ? 
              <ClassicalParametersConfig /> : 
              <QuantumParametersConfig />
            }
            
            {/* Visualization tabs for switching between original and optimized structures */}
            {optimizationResult && (
              <div style={{ marginBottom: "10px", marginTop: "20px" }}>
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
                  onClick={() => setActiveView("optimized")}
                  style={{
                    backgroundColor: activeView === "optimized" ? 
                      (optimizationType === "classical" ? "#28a745" : "#17a2b8") : "#ccc",
                    color: "white",
                    border: "none",
                    borderRadius: "0 5px 5px 0",
                    padding: "5px 10px",
                    cursor: "pointer",
                  }}
                >
                  {optimizationType === "classical" ? "Classical" : "Quantum"} Optimized
                </button>
              </div>
            )}
            
            <h3>{activeView === "original" ? "Original" : 
                 (optimizationType === "classical" ? "Classical" : "Quantum") + " Optimized"} Structure:</h3>
            <MoleculeViewer />
            
            {/* Action Buttons */}
            <OptimizeButton />
          </>
        )}
        
        {/* Optimization Results */}
        {optimizationResult && (
          <div style={{ margin: "20px" }}>
            <h2>{optimizationType === "classical" ? "Classical" : "Quantum"} Optimization Results</h2>
            
            <div style={{ padding: "15px", border: "1px solid #ddd", borderRadius: "5px", textAlign: "left" }}>
              {optimizationResult.result.error ? (
                <p style={{ color: "red" }}>Error: {optimizationResult.result.error}</p>
              ) : (
                <div>
                  <h4>Metadata:</h4>
                  <p><strong>Method:</strong> {optimizationResult.result.metadata.method}</p>
                  <p><strong>Library:</strong> {optimizationResult.result.metadata.library}</p>
                  
                  {optimizationType === "classical" ? (
                    <>
                      <p><strong>Temperature:</strong> {optimizationResult.result.metadata.parameters.temperature} K</p>
                      <p><strong>Max Iterations:</strong> {optimizationResult.result.metadata.parameters.max_iterations}</p>
                      <p><strong>Final Energy:</strong> {optimizationResult.result.metadata.final_energy_kj_mol} kJ/mol</p>
                      <p><strong>Convergence:</strong> {optimizationResult.result.metadata.convergence || "Unknown"}</p>
                      <p><strong>Bonds Detected:</strong> {optimizationResult.result.metadata.bonds_detected}</p>
                      <p><strong>Angles Detected:</strong> {optimizationResult.result.metadata.angles_detected}</p>
                    </>
                  ) : (
                    <>
                      <p><strong>Basis Set:</strong> {optimizationResult.result.metadata.parameters.basis}</p>
                      <p><strong>Max Iterations:</strong> {optimizationResult.result.metadata.parameters.max_iterations}</p>
                      <p><strong>Theory Level:</strong> {optimizationResult.result.metadata.theory_level}</p>
                      <p><strong>Final Energy:</strong> {optimizationResult.result.metadata.final_energy_hartree} Hartree</p>
                      <p><strong>Iterations:</strong> {optimizationResult.result.metadata.iterations}</p>
                      <p><strong>Converged:</strong> {optimizationResult.result.metadata.converged ? "Yes" : "No"}</p>
                    </>
                  )}
                  
                  <p><strong>Duration:</strong> {optimizationResult.result.metadata.duration_seconds} seconds</p>
                  
                  <button
                    onClick={handleDownload}
                    style={{
                      backgroundColor: optimizationType === "classical" ? "#28a745" : "#17a2b8",
                      color: "white",
                      border: "none",
                      borderRadius: "5px",
                      padding: "5px 10px",
                      cursor: "pointer",
                    }}
                  >
                    Download Result
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