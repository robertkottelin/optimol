import React, { useState, useEffect, useContext } from "react";
import axios from "axios";
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import ReactMarkdown from "react-markdown";

// Import authentication context
import { AuthContext } from './AuthContext';
import LoginForm from './components/LoginForm';
import RegisterForm from './components/RegisterForm';

// Import style constants
import { styles } from './styles/components';
import {
  ITERATION_LIMITS,
  defaultClassicalParams,
  defaultQuantumParams,
  TEST_MOLECULES
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

// Configure axios to include credentials with all requests
axios.defaults.withCredentials = true;

const App = () => {
  // Authentication context
  const { currentUser, isAuthenticated, isSubscribed, isLoading, logout, token } = useContext(AuthContext);
  
  // FIXED: Moved all useState declarations to component top level
  const [showLoginForm, setShowLoginForm] = useState(true);
  const [moleculeData, setMoleculeData] = useState(null);
  const [optimizationResult, setOptimizationResult] = useState(null);
  const [activeView, setActiveView] = useState("original"); // original, optimized
  const [optimizationType, setOptimizationType] = useState("classical"); // classical or quantum
  const [classicalParams, setClassicalParams] = useState({ ...defaultClassicalParams });
  const [quantumParams, setQuantumParams] = useState({ ...defaultQuantumParams });
  const [showAdvancedParams, setShowAdvancedParams] = useState(false);
  const [isHowToUseVisible, setIsHowToUseVisible] = useState(false);
  const [isDragActive, setIsDragActive] = useState(false);
  const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false);
  const [isSubscribeLoading, setIsSubscribeLoading] = useState(false);
  const [isCancelLoading, setIsCancelLoading] = useState(false);
  const [isOptimizeLoading, setIsOptimizeLoading] = useState(false);
  const [isMobile, setIsMobile] = useState(false);
  const [serverHealth, setServerHealth] = useState(null);
  const [isCheckingHealth, setIsCheckingHealth] = useState(false);
  const [serverHealthDetails, setServerHealthDetails] = useState(null);
  // FIXED: Moved howToUseContent useState before conditional return
  const [howToUseContent, setHowToUseContent] = useState("");

  const apiBaseUrl = "https://optimizemolecule.com";
  
  // Group all useEffect hooks together
  useEffect(() => {
    const checkIfMobile = () => {
      setIsMobile(window.innerWidth <= 768);
    };

    // Initial check
    checkIfMobile();

    // Add event listener for window resize
    window.addEventListener('resize', checkIfMobile);

    // Cleanup
    return () => window.removeEventListener('resize', checkIfMobile);
  }, []);
  
  useEffect(() => {
    // Load documentation from public directory
    fetch(`${process.env.PUBLIC_URL}/how-to-use.md`)
      .then(response => {
        if (!response.ok) {
          throw new Error('Documentation file not found');
        }
        return response.text();
      })
      .then(text => setHowToUseContent(text))
      .catch(error => {
        console.error("Failed to load documentation:", error);
        // Fallback content
        setHowToUseContent("# Molecular Optimization System\n\nDocumentation is currently unavailable.");
      });
  }, []);

  // Apply limits on initial load and whenever subscription status changes
  useEffect(() => {
    if (!isLoading) {
      applyIterationLimits(isSubscribed);
    }
  }, [isLoading, isSubscribed]);

  // Helper function to consistently apply iteration limits
  const applyIterationLimits = (isUserSubscribed) => {
    // Apply limits to classical parameters - preserving existing values
    setClassicalParams(prevParams => {
      const classicalMaxIterations = isUserSubscribed
        ? ITERATION_LIMITS.subscribed.classical
        : ITERATION_LIMITS.unsubscribed.classical;

      return {
        ...prevParams, // Preserve existing params including force_iterations
        max_iterations: Math.min(
          prevParams.max_iterations,
          classicalMaxIterations
        )
      };
    });

    // Apply limits to quantum parameters - preserving existing values
    setQuantumParams(prevParams => {
      const quantumMaxIterations = isUserSubscribed
        ? ITERATION_LIMITS.subscribed.quantum
        : ITERATION_LIMITS.unsubscribed.quantum;

      // Apply basis set restrictions for non-subscribers
      let updatedBasis = prevParams.basis;
      if (!isUserSubscribed && (updatedBasis === "6-311g" || updatedBasis === "cc-pvdz")) {
        updatedBasis = "6-31g";
      }

      return {
        ...prevParams, // Preserve existing params
        max_iterations: Math.min(
          prevParams.max_iterations,
          quantumMaxIterations
        ),
        basis: updatedBasis
      };
    });
  };

  const handleShowHowToUse = () => {
    setIsHowToUseVisible(true);
    setIsMobileMenuOpen(false);
  };

  const handleClosePopup = () => {
    setIsHowToUseVisible(false);
  };

  const toggleMobileMenu = () => {
    setIsMobileMenuOpen(!isMobileMenuOpen);
  };

  const checkServerHealth = async () => {
    setIsCheckingHealth(true);
    setServerHealthDetails(null);

    try {
      const response = await axios.get(`${apiBaseUrl}/health`);
      setServerHealth(response.data.status === "healthy");
      setServerHealthDetails({
        status: response.status,
        statusText: response.statusText,
        data: response.data
      });
      setTimeout(() => {
        setServerHealth(null);
        setServerHealthDetails(null);
      }, 5000); // Extended display time to 5 seconds for better readability
    } catch (error) {
      console.error("Error checking server health:", error);
      setServerHealth(false);

      // Capture detailed error information
      setServerHealthDetails({
        status: error.response?.status || 'Network Error',
        statusText: error.response?.statusText || error.message,
        data: error.response?.data || {}
      });

      setTimeout(() => {
        setServerHealth(null);
        setServerHealthDetails(null);
      }, 5000);
    } finally {
      setIsCheckingHealth(false);
    }
  };

  const handleSubscriptionSuccess = () => {
    // Reload user data from context to update subscription status
    window.location.reload();
  };

  const handleCancelSubscription = async () => {
    const confirmCancel = window.confirm(
      "Are you sure you want to cancel your subscription? This action cannot be undone."
    );

    if (!confirmCancel) return;

    setIsCancelLoading(true);

    try {
      // Token is handled by axios interceptor in AuthContext
      const response = await axios.post(`${apiBaseUrl}/cancel-subscription`, {});

      if (response.data.success) {
        alert("Your subscription has been canceled.");
        // Refresh page to update auth context
        window.location.reload();
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
      alert("Please upload or select a molecule first.");
      return;
    }

    setIsOptimizeLoading(true);

    try {
      // Get the correct parameters based on selected optimization type
      const optimizationParams =
        optimizationType === "classical" ? { ...classicalParams } : { ...quantumParams };

      // Apply iteration limits for all users
      if (optimizationType === "classical") {
        const maxIterations = isSubscribed
          ? ITERATION_LIMITS.subscribed.classical
          : ITERATION_LIMITS.unsubscribed.classical;

        optimizationParams.max_iterations = Math.min(
          optimizationParams.max_iterations,
          maxIterations
        );
      } else {
        const maxIterations = isSubscribed
          ? ITERATION_LIMITS.subscribed.quantum
          : ITERATION_LIMITS.unsubscribed.quantum;

        optimizationParams.max_iterations = Math.min(
          optimizationParams.max_iterations,
          maxIterations
        );

        if (!isSubscribed && (optimizationParams.basis === "6-311g" || optimizationParams.basis === "cc-pvdz")) {
          optimizationParams.basis = "6-31g";
        }
      }

      const payload = {
        molecule: moleculeData,
        optimization_type: optimizationType,
        optimization_params: optimizationParams
      };

      console.log('Optimization payload:', JSON.stringify(payload, null, 2));

      console.log("Attempting request to:", `${apiBaseUrl}/optimize-molecule`);
      
      // Token is handled by axios interceptor in AuthContext
      const response = await axios({
        method: 'post',
        url: `${apiBaseUrl}/optimize-molecule`,
        data: payload,
        headers: {
          'Content-Type': 'application/json',
          'Accept': 'application/json'
        }
      });

      console.log("Response received:", response);

      if (response.data.success) {
        setOptimizationResult(response.data);
        setActiveView("optimized");
      } else {
        alert("Optimization failed. " + (response.data.error || ""));
      }
    } catch (error) {
      console.error("Error optimizing molecule:", error);
      console.error("Request error details:", {
        message: error.message,
        response: error.response,
        request: error.request
      });
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
      // Keep current force_iterations setting when resetting other parameters
      setClassicalParams(prevParams => {
        const params = { ...defaultClassicalParams };

        // Preserve force_iterations setting
        params.force_iterations = prevParams.force_iterations;

        // Always apply appropriate limits based on subscription status
        const maxIterations = isSubscribed
          ? ITERATION_LIMITS.subscribed.classical
          : ITERATION_LIMITS.unsubscribed.classical;

        params.max_iterations = Math.min(
          params.max_iterations,
          maxIterations
        );

        return params;
      });
    } else {
      // For quantum parameters - apply subscription-based limits
      setQuantumParams(prevParams => {
        const params = { ...defaultQuantumParams };

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

        return params;
      });
    }
  };

  // Handle test molecule selection
  const handleTestMoleculeSelect = (moleculeKey) => {
    if (TEST_MOLECULES[moleculeKey]) {
      setMoleculeData(TEST_MOLECULES[moleculeKey]);
      setOptimizationResult(null);
      setActiveView("original");
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

  // FIXED: Now safe to use conditional returns after all hooks are declared
  // Display loading indicator while auth state is determined
  if (isLoading) {
    return (
      <div style={{ 
        display: 'flex', 
        justifyContent: 'center', 
        alignItems: 'center', 
        height: '100vh',
        backgroundColor: '#0c1021',
        color: '#f0f4f8'
      }}>
        <div>
          <div className="spinner" style={{ 
            width: '40px', 
            height: '40px', 
            border: '3px solid rgba(56, 189, 248, 0.3)',
            borderTopColor: '#38bdf8',
            borderRadius: '50%',
            animation: 'spin 1s linear infinite',
            margin: '0 auto 20px auto'
          }} />
          <div>Loading authentication...</div>
        </div>
      </div>
    );
  }

  // Display login/register if not authenticated
  if (!isAuthenticated) {
    return (
      <div style={{ ...styles.app, padding: isMobile ? '16px' : styles.app.padding }}>
        <div style={styles.decorativeBg}></div>
        <div style={{ ...styles.decorativeLine, top: "15%", animationDelay: "0s" }}></div>
        <div style={{ ...styles.decorativeLine, top: "35%", animationDelay: "0.5s" }}></div>
        <div style={{ ...styles.decorativeLine, top: "65%", animationDelay: "1s" }}></div>
        <div style={{ ...styles.decorativeLine, top: "85%", animationDelay: "1.5s" }}></div>
        
        <div style={styles.container}>
          <header style={styles.header} className="app-header">
            <h1 style={styles.headerTitle} className="app-title">Molecular Optimization System</h1>
            <p style={styles.headerSubtitle} className="app-subtitle">
              Advanced computational chemistry tools for structure optimization
            </p>
          </header>
          
          <div className="fade-in glass card" style={styles.card}>
            {showLoginForm ? (
              <LoginForm toggleForm={() => setShowLoginForm(false)} />
            ) : (
              <RegisterForm toggleForm={() => setShowLoginForm(true)} />
            )}
          </div>
        </div>
      </div>
    );
  }

  // Main App render
  return (
    <div style={{ ...styles.app, padding: isMobile ? '16px' : styles.app.padding }}>
      <div style={styles.decorativeBg}></div>

      {/* Decorative animated lines for cyberpunk effect */}
      <div style={{ ...styles.decorativeLine, top: "15%", animationDelay: "0s" }}></div>
      <div style={{ ...styles.decorativeLine, top: "35%", animationDelay: "0.5s" }}></div>
      <div style={{ ...styles.decorativeLine, top: "65%", animationDelay: "1s" }}></div>
      <div style={{ ...styles.decorativeLine, top: "85%", animationDelay: "1.5s" }}></div>

      <div style={styles.container}>
        {/* Top Action Buttons - Move before header for mobile */}
        <div className="top-buttons-container" style={{
          display: 'flex',
          flexWrap: 'wrap',
          gap: '10px',
          justifyContent: isMobile ? 'center' : 'flex-start',
          position: 'relative',
          zIndex: 10,
          marginBottom: '20px'
        }}>
          <button
            onClick={handleShowHowToUse}
            style={{
              ...styles.howToUseButton,
              position: 'static',
              marginRight: '5px'
            }}
            className="float howToUseButton"
          >
            <span style={styles.howToUseIcon}><Icons.book /></span>
            Documentation & Theory
          </button>

          <button
            onClick={checkServerHealth}
            disabled={isCheckingHealth}
            style={{
              ...styles.button,
              background: serverHealth === null
                ? "linear-gradient(145deg, rgba(59, 130, 246, 0.9), rgba(37, 99, 235, 0.9))"
                : serverHealth
                  ? "linear-gradient(145deg, rgba(16, 185, 129, 0.9), rgba(5, 150, 105, 0.9))"
                  : "linear-gradient(145deg, rgba(244, 63, 94, 0.9), rgba(225, 29, 72, 0.9))",
              color: "white",
              padding: "4px 12px",
              borderRadius: "6px",
              fontSize: "0.75rem",
              fontWeight: "600",
              position: 'static',
              boxShadow: "0 4px 6px rgba(0, 0, 0, 0.1)",
              display: "flex",
              alignItems: "center",
              gap: "5px",
              flexDirection: serverHealthDetails ? "column" : "row",
              alignItems: serverHealthDetails ? "flex-start" : "center",
              minWidth: serverHealthDetails ? "180px" : "auto"
            }}
            className="health-check-button"
          >
            {isCheckingHealth ? (
              <>
                <span className="spin" style={{ display: "inline-block" }}>
                  <Icons.spinner />
                </span>
                Checking...
              </>
            ) : (
              <>
                <div style={{ display: "flex", alignItems: "center", gap: "5px" }}>
                  {serverHealth === true && <span><Icons.checkmark /></span>}
                  {serverHealth === false && <span><Icons.warning /></span>}
                  {serverHealth === null && <span><Icons.info /></span>}
                  Server Health
                </div>

                {serverHealthDetails && (
                  <div style={{
                    fontSize: "0.7rem",
                    marginTop: "4px",
                    backgroundColor: "rgba(0,0,0,0.2)",
                    padding: "3px 6px",
                    borderRadius: "4px",
                    width: "100%"
                  }}>
                    <div>Status: {serverHealthDetails.status}</div>
                    <div style={{ wordBreak: "break-word" }}>
                      {serverHealthDetails.statusText}
                    </div>
                  </div>
                )}
              </>
            )}
          </button>

          <button
            onClick={logout}
            style={{
              ...styles.button,
              background: "linear-gradient(145deg, rgba(100, 116, 139, 0.9), rgba(71, 85, 105, 0.9))",
              color: "white",
              padding: "4px 12px",
              borderRadius: "6px",
              fontSize: "0.75rem",
              fontWeight: "600",
              position: 'static',
              boxShadow: "0 4px 6px rgba(0, 0, 0, 0.1)",
              display: "flex",
              alignItems: "center",
              gap: "5px"
            }}
          >
            <Icons.close />
            Logout
          </button>

          {isSubscribed && (
            <button
              onClick={handleCancelSubscription}
              disabled={isCancelLoading}
              style={{
                ...styles.cancelSubscriptionButton,
                opacity: isCancelLoading ? 0.7 : 1,
                cursor: isCancelLoading ? "not-allowed" : "pointer",
                position: 'static'
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
        </div>

        {/* App Header - Now comes after buttons in the DOM */}
        <header style={styles.header} className="app-header">
          <h1 style={styles.headerTitle} className="app-title">Molecular Optimization System</h1>
          <p style={styles.headerSubtitle} className="app-subtitle">
            Advanced computational chemistry tools for structure optimization
          </p>
        </header>

        {/* Subscription Form or Welcome Message */}
        <div className="fade-in glass card" style={isSubscribed ? styles.cardWithGlow : styles.card}>
          {!isSubscribed ? (
            <div className={`subscription-form ${isMobile ? 'mobile-smaller-padding' : ''}`}>
              <Elements stripe={stripePromise}>
                <SubscriptionForm
                  onSuccess={handleSubscriptionSuccess}
                  isMobile={isMobile}
                />
              </Elements>
            </div>
          ) : (
            <div className={`welcome-message ${isMobile ? 'mobile-stack mobile-text-center' : ''}`} style={styles.welcomeMessage}>
              <div style={styles.welcomeIcon}>
                <Icons.verified />
              </div>
              <p style={{ margin: 0 }}>Welcome, <strong>{currentUser.email}</strong>! You have full access to all computational capabilities with your premium subscription.</p>
            </div>
          )}
        </div>

        {/* Main Optimization Section */}
        <section style={styles.section}>
          <h2 style={styles.sectionTitle}>
            Molecule Optimization
            <span style={styles.sectionTitleUnderline}></span>
          </h2>

          {/* Test Molecules Buttons */}
          <div style={styles.testMoleculesContainer} className={isMobile ? 'mobile-smaller-padding' : ''}>
            <h3 style={styles.testMoleculesTitle}>
              Test Molecules
              <span style={styles.testMoleculesTitleIcon}></span>
            </h3>
            <div style={styles.testMoleculesButtonContainer} className={isMobile ? 'mobile-stack' : ''}>
              <button
                onClick={() => handleTestMoleculeSelect('water')}
                style={styles.testMoleculeButton}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
              >
                <span style={styles.testMoleculeIcon}><Icons.molecule /></span>
                Water (H₂O)
              </button>
              <button
                onClick={() => handleTestMoleculeSelect('ibuprofen')}
                style={styles.testMoleculeButton}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width' : ''}`}
              >
                <span style={styles.testMoleculeIcon}><Icons.molecule /></span>
                Ibuprofen (C₁₃H₁₈O₂)
              </button>
            </div>
          </div>

          {/* File Upload Area */}
          <div
            className={`file-upload-area ${isDragActive ? 'file-upload-active' : ''} ${isMobile ? 'mobile-smaller-padding' : ''}`}
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
            <p style={styles.fileUploadText}>
              {isMobile ? 'Upload a molecule file' : 'Drag & drop a molecule file here, or click to select a file'}
            </p>
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
              <div
                style={styles.methodSelectionContainer}
                className={isMobile ? "mobile-stack" : ""}
              >
                <div
                  className={`method-card ${optimizationType === "classical" ? 'classical-active' : ''} ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
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
                  className={`method-card ${optimizationType === "quantum" ? 'quantum-active' : ''} ${isMobile ? 'mobile-full-width' : ''}`}
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
                  isMobile={isMobile}
                />
              ) : (
                <QuantumParametersConfig
                  isSubscribed={isSubscribed}
                  quantumParams={quantumParams}
                  showAdvancedParams={showAdvancedParams}
                  handleParamChange={handleParamChange}
                  handleResetParams={handleResetParams}
                  setShowAdvancedParams={setShowAdvancedParams}
                  isMobile={isMobile}
                />
              )}

              {/* Visualization */}
              <div style={styles.visualizationContainer} className={`glass ${isMobile ? 'mobile-smaller-padding' : ''}`}>
                <div style={styles.visualizationHeader} className={isMobile ? 'mobile-stack' : ''}>
                  <div style={styles.visualizationTitle}>
                    <span style={styles.visualizationIcon}><Icons.molecule /></span>
                    {activeView === "original" ? "Original" :
                      (optimizationType === "classical" ? "Classical" : "Quantum") + " Optimized"} Structure
                  </div>

                  {optimizationResult && (
                    <div style={styles.tabs} className={isMobile ? 'mobile-full-width' : ''}>
                      <div
                        style={styles.tab(activeView === "original", "#38bdf8")}
                        onClick={() => setActiveView("original")}
                        className={isMobile ? 'mobile-smaller-text' : ''}
                      >
                        <Icons.molecule /> Original
                      </div>
                      <div
                        style={styles.tab(
                          activeView === "optimized",
                          optimizationType === "classical" ? "#10b981" : "#38bdf8"
                        )}
                        onClick={() => setActiveView("optimized")}
                        className={isMobile ? 'mobile-smaller-text' : ''}
                      >
                        {optimizationType === "classical" ? <Icons.classical /> : <Icons.quantum />}
                        {optimizationType === "classical" ? "Classical" : "Quantum"} Optimized
                      </div>
                    </div>
                  )}
                </div>

                <MoleculeViewer atoms={getAtoms()} isMobile={isMobile} />
              </div>

              {/* Optimize Button */}
              <div style={styles.optimizeButtonContainer} className={isMobile ? 'mobile-full-width' : ''}>
                <button
                  onClick={handleOptimize}
                  disabled={isOptimizeLoading}
                  style={styles.optimizeButton(optimizationType, isOptimizeLoading)}
                  className={`${optimizationType === "classical" ? "classical-button" : "quantum-button"} ${isMobile ? 'mobile-full-width' : ''}`}
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
                  <div style={styles.freeUserNotice} className={isMobile ? 'mobile-smaller-text mobile-text-center' : ''}>
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
                  isMobile={isMobile}
                />
              )}
            </div>
          )}
        </section>
      </div>

      {/* Mobile Menu Button */}
      {isMobile && (
        <button
          className="mobile-menu-button"
          onClick={toggleMobileMenu}
        >
          <Icons.menu />
        </button>
      )}

      {/* Mobile Navigation Menu */}
      {isMobile && isMobileMenuOpen && (
        <div className="mobile-nav">
          <button
            className={`mobile-nav-button ${moleculeData ? '' : 'active'}`}
            onClick={() => {
              setIsMobileMenuOpen(false);
              document.querySelector('.file-upload-area').scrollIntoView({ behavior: 'smooth' });
            }}
          >
            <span className="mobile-nav-icon"><Icons.upload /></span>
            Upload
          </button>
          {moleculeData && (
            <>
              <button
                className={`mobile-nav-button ${optimizationType === 'classical' ? 'active' : ''}`}
                onClick={() => {
                  setIsMobileMenuOpen(false);
                  handleOptimizationTypeChange("classical");
                }}
              >
                <span className="mobile-nav-icon"><Icons.classical /></span>
                Classical
              </button>
              <button
                className={`mobile-nav-button ${optimizationType === 'quantum' ? 'active' : ''}`}
                onClick={() => {
                  setIsMobileMenuOpen(false);
                  handleOptimizationTypeChange("quantum");
                }}
              >
                <span className="mobile-nav-icon"><Icons.quantum /></span>
                Quantum
              </button>
              <button
                className="mobile-nav-button"
                onClick={() => {
                  setIsMobileMenuOpen(false);
                  document.querySelector('.optimize-button-container').scrollIntoView({ behavior: 'smooth' });
                }}
              >
                <span className="mobile-nav-icon"><Icons.play /></span>
                Optimize
              </button>
            </>
          )}
          <button
            className="mobile-nav-button"
            onClick={handleShowHowToUse}
          >
            <span className="mobile-nav-icon"><Icons.book /></span>
            Help
          </button>
          <button
            className="mobile-nav-button"
            onClick={() => {
              setIsMobileMenuOpen(false);
              checkServerHealth();
            }}
            disabled={isCheckingHealth}
            style={{
              position: 'relative'
            }}
          >
            <span className="mobile-nav-icon">
              {isCheckingHealth ? <Icons.spinner /> :
                (serverHealth === true ? <Icons.checkmark /> :
                  (serverHealth === false ? <Icons.warning /> : <Icons.info />))}
            </span>
            Health
            {serverHealthDetails && (
              <div style={{
                position: 'absolute',
                top: '-35px',
                left: '0',
                backgroundColor: 'rgba(0,0,0,0.8)',
                padding: '3px 6px',
                borderRadius: '4px',
                fontSize: '0.6rem',
                whiteSpace: 'nowrap',
                zIndex: 100
              }}>
                {serverHealthDetails.status}: {serverHealthDetails.statusText}
              </div>
            )}
          </button>
        </div>
      )}

      {/* How To Use Popup */}
      {isHowToUseVisible && (
        <div style={styles.popup} className="popup-overlay">
          <div
            style={styles.popupContent}
            className={`popup-content glass ${isMobile ? 'mobile-smaller-padding' : ''}`}
          >
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