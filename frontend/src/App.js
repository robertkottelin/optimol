import React, { useState, useRef, useEffect } from "react";
import axios from "axios";
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import ReactMarkdown from "react-markdown";
import howToUseContent from '!raw-loader!../how-to-use.md';

// Import style constants
import { styles } from './styles/components';
import { 
  ITERATION_LIMITS, 
  defaultClassicalParams, 
  defaultQuantumParams 
} from './styles/constants';

// Import components
import { Icons } from './components/Icons';
import SubscriptionForm from './components/SubscriptionForm';
import MoleculeViewer from './components/MoleculeViewer';
import { 
  ClassicalParametersConfig, 
  QuantumParametersConfig 
} from './components/ParameterConfig';
import OptimizationResults from './components/OptimizationResults';

// Initialize Stripe
const stripePromise = loadStripe('pk_test_kbl0ETzPsoiTwU4ZJMvhsYJw006XVnV4Aq');

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
  const [isDragActive, setIsDragActive] = useState(false);
  
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

  const handleShowHowToUse = () => {
    setIsHowToUseVisible(true);
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

  const handleDragOver = (e) => {
    e.preventDefault();
    setIsDragActive(true);
  };

  const handleDragLeave = () => {
    setIsDragActive(false);
  };

  const handleFileDrop = (e) => {
    e.preventDefault();
    setIsDragActive(false);
    
    if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
      const file = e.dataTransfer.files[0];
      
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
    }
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

  const handleOptimizationTypeChange = (newType) => {
    // If optimization type is changing, reset the active view to original
    if (newType !== optimizationType) {
      setActiveView("original");
    }
    
    setOptimizationType(newType);
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
        // Clear any previous result for the other optimization type
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

  // Main App render
  return (
    <div style={styles.app}>
      <div style={styles.decorativeBg}></div>
      
      {/* Decorative animated lines for cyberpunk effect */}
      <div style={{ ...styles.decorativeLine, top: "15%", animationDelay: "0s" }}></div>
      <div style={{ ...styles.decorativeLine, top: "35%", animationDelay: "0.5s" }}></div>
      <div style={{ ...styles.decorativeLine, top: "65%", animationDelay: "1s" }}></div>
      <div style={{ ...styles.decorativeLine, top: "85%", animationDelay: "1.5s" }}></div>
      
      <div style={styles.container}>
        {/* App Header */}
        <header style={styles.header} className="app-header">
          <h1 style={styles.headerTitle} className="app-title">Molecular Optimization System</h1>
          <p style={styles.headerSubtitle} className="app-subtitle">
            Advanced computational chemistry tools for structure optimization
          </p>
        </header>
        
        {/* Top Action Buttons */}
        <button
          onClick={handleShowHowToUse}
          style={styles.howToUseButton}
          className="float"
        >
          <span style={styles.howToUseIcon}><Icons.book /></span>
          Documentation & Theory
        </button>
        
        {isSubscribed && (
          <button
            onClick={handleCancelSubscription}
            disabled={isCancelLoading}
            style={{
              ...styles.cancelSubscriptionButton,
              opacity: isCancelLoading ? 0.7 : 1,
              cursor: isCancelLoading ? "not-allowed" : "pointer",
            }}
          >
            {isCancelLoading ? (
              <>
                <span className="spin" style={{ display: "inline-block", marginRight: "5px" }}>
                  <Icons.spinner />
                </span>
                Cancelling...
              </>
            ) : (
              <>
                <span style={styles.cancelSubscriptionIcon}><Icons.cancel /></span>
                Cancel Subscription
              </>
            )}
          </button>
        )}
        
        {/* Subscription Form or Welcome Message */}
        <div className="fade-in glass card" style={isSubscribed ? styles.cardWithGlow : styles.card}>
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
            <div className="welcome-message" style={styles.welcomeMessage}>
              <div style={styles.welcomeIcon}>
                <Icons.verified />
              </div>
              <p style={{ margin: 0 }}>Welcome, <strong>{userEmail}</strong>! You have full access to all computational capabilities with your premium subscription.</p>
            </div>
          )}
        </div>
        
        {/* Main Optimization Section */}
        <section style={styles.section}>
          <h2 style={styles.sectionTitle}>
            Molecule Optimization
            <span style={styles.sectionTitleUnderline}></span>
          </h2>
          
          {/* File Upload Area */}
          <div 
            className={`file-upload-area ${isDragActive ? 'file-upload-active' : ''}`}
            style={{
              ...styles.fileUpload,
              ...(isDragActive ? styles.fileUploadActive : {})
            }}
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            onDrop={handleFileDrop}
          >
            <div style={isDragActive ? { ...styles.fileUploadIcon, ...styles.fileUploadActiveIcon } : styles.fileUploadIcon}>
              <Icons.upload />
            </div>
            <p style={styles.fileUploadText}>Drag & drop a molecule file here, or click to select a file</p>
            <p style={styles.fileUploadSubtext}>
              Accepted format: JSON structured molecule data
            </p>
            <input
              type="file"
              onChange={handleFileUpload}
              style={{ 
                position: "absolute", 
                top: 0, 
                left: 0, 
                width: "100%", 
                height: "100%", 
                opacity: 0,
                cursor: "pointer" 
              }}
            />
          </div>
          
          {moleculeData && (
            <div className="slide-up" style={styles.contentWrapper}>
              {/* Optimization Method Selection */}
              <div style={styles.methodSelectionContainer}>
                <div 
                  className={`method-card ${optimizationType === "classical" ? 'classical-active' : ''}`}
                  style={styles.methodSelectionButton(optimizationType === "classical", "classical")}
                  onClick={() => handleOptimizationTypeChange("classical")}
                >
                  <span style={{
                    ...styles.methodIcon,
                    backgroundColor: optimizationType === "classical" ? "rgba(16, 185, 129, 0.1)" : "rgba(30, 41, 59, 0.5)",
                    border: optimizationType === "classical" ? "1px solid rgba(16, 185, 129, 0.3)" : "transparent",
                    padding: "8px",
                    borderRadius: "50%",
                    marginBottom: "16px"
                  }}>
                    <Icons.classical />
                  </span>
                  <div style={styles.methodTitle}>Classical Optimization</div>
                  <div style={styles.methodDescription}>
                    Molecular mechanics optimization using empirical force fields. Faster calculations suitable for larger molecular systems.
                  </div>
                </div>
                
                <div 
                  className={`method-card ${optimizationType === "quantum" ? 'quantum-active' : ''}`}
                  style={styles.methodSelectionButton(optimizationType === "quantum", "quantum")}
                  onClick={() => handleOptimizationTypeChange("quantum")}
                >
                  <span style={{
                    ...styles.methodIcon,
                    backgroundColor: optimizationType === "quantum" ? "rgba(56, 189, 248, 0.1)" : "rgba(30, 41, 59, 0.5)",
                    border: optimizationType === "quantum" ? "1px solid rgba(56, 189, 248, 0.3)" : "transparent",
                    padding: "8px",
                    borderRadius: "50%",
                    marginBottom: "16px"
                  }}>
                    <Icons.quantum />
                  </span>
                  <div style={styles.methodTitle}>Quantum Optimization</div>
                  <div style={styles.methodDescription}>
                    Ab initio quantum chemistry methods for precise electronic structure optimization with quantum mechanical accuracy.
                  </div>
                </div>
              </div>
              
              {/* Parameter Configuration */}
              {optimizationType === "classical" ? (
                <ClassicalParametersConfig 
                  isSubscribed={isSubscribed}
                  classicalParams={classicalParams}
                  showAdvancedParams={showAdvancedParams}
                  handleParamChange={handleParamChange}
                  handleResetParams={handleResetParams}
                  setShowAdvancedParams={setShowAdvancedParams}
                />
              ) : (
                <QuantumParametersConfig 
                  isSubscribed={isSubscribed}
                  quantumParams={quantumParams}
                  showAdvancedParams={showAdvancedParams}
                  handleParamChange={handleParamChange}
                  handleResetParams={handleResetParams}
                  setShowAdvancedParams={setShowAdvancedParams}
                />
              )}
              
              {/* Visualization */}
              <div style={styles.visualizationContainer} className="glass">
                <div style={styles.visualizationHeader}>
                  <div style={styles.visualizationTitle}>
                    <span style={styles.visualizationIcon}><Icons.molecule /></span>
                    {activeView === "original" ? "Original" : 
                     (optimizationType === "classical" ? "Classical" : "Quantum") + " Optimized"} Structure
                  </div>
                  
                  {optimizationResult && (
                    <div style={styles.tabs}>
                      <div 
                        style={styles.tab(activeView === "original", "#38bdf8")}
                        onClick={() => setActiveView("original")}
                      >
                        <Icons.molecule /> Original
                      </div>
                      <div 
                        style={styles.tab(
                          activeView === "optimized", 
                          optimizationType === "classical" ? "#10b981" : "#38bdf8"
                        )}
                        onClick={() => setActiveView("optimized")}
                      >
                        {optimizationType === "classical" ? <Icons.classical /> : <Icons.quantum />}
                        {optimizationType === "classical" ? "Classical" : "Quantum"} Optimized
                      </div>
                    </div>
                  )}
                </div>
                
                <MoleculeViewer atoms={getAtoms()} />
              </div>
              
              {/* Optimize Button */}
              <div style={styles.optimizeButtonContainer}>
                <button
                  onClick={handleOptimize}
                  disabled={isOptimizeLoading}
                  style={styles.optimizeButton(optimizationType, isOptimizeLoading)}
                  className={optimizationType === "classical" ? "classical-button" : "quantum-button"}
                >
                  {isOptimizeLoading ? (
                    <>
                      <span className="spin" style={{ display: "inline-block", marginRight: "12px" }}>
                        <Icons.spinner />
                      </span>
                      Optimizing...
                    </>
                  ) : (
                    <>
                      {`Run ${optimizationType === "classical" ? "Classical" : "Quantum"} Optimization${!isSubscribed ? " (Limited)" : ""}`}
                      <div style={styles.optimizeButtonShine}></div>
                    </>
                  )}
                </button>
                
                {!isSubscribed && (
                  <div style={styles.freeUserNotice}>
                    Free users are limited to {ITERATION_LIMITS.unsubscribed.classical.toLocaleString()} iterations for classical and {ITERATION_LIMITS.unsubscribed.quantum} for quantum optimizations.
                  </div>
                )}
              </div>
              
              {/* Optimization Results */}
              {optimizationResult && (
                <OptimizationResults 
                  optimizationResult={optimizationResult}
                  optimizationType={optimizationType}
                  handleDownload={handleDownload}
                />
              )}
            </div>
          )}
        </section>
      </div>
      
      {/* How To Use Popup */}
      {isHowToUseVisible && (
        <div style={styles.popup} className="popup-overlay">
          <div style={styles.popupContent} className="popup-content glass">
            <button
              onClick={handleClosePopup}
              style={styles.popupClose}
            >
              <Icons.close />
            </button>
            
            <div style={styles.popupScroll}>
              <ReactMarkdown>
                {howToUseContent}
              </ReactMarkdown>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default App;