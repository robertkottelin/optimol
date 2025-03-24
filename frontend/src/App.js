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
  QUANTUM_SIZE_LIMITS,
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

// Initialize Stripe with public key
const stripePromise = loadStripe('pk_test_kbl0ETzPsoiTwU4ZJMvhsYJw006XVnV4Aq');

// Configure axios to include credentials with all requests
axios.defaults.withCredentials = true;

const App = () => {
  // Authentication context
  const { currentUser, isAuthenticated, isSubscribed, isLoading, logout, token } = useContext(AuthContext);

  // FIXED: Moved all useState declarations to component top level
  const [showAuthModal, setShowAuthModal] = useState(false);
  const [showLoginForm, setShowLoginForm] = useState(true);

  // Updated state for dual molecule support
  const [molecule1Data, setMolecule1Data] = useState(null);
  const [molecule2Data, setMolecule2Data] = useState(null);
  const [activeMolecule, setActiveMolecule] = useState(1); // 1 or 2, indicating which molecule is active
  const [interactionMode, setInteractionMode] = useState(false); // Whether to optimize molecule interaction
  const [positioningMode, setPositioningMode] = useState(false); // Whether positioning mode is active
  const [molecule2Offset, setMolecule2Offset] = useState({ x: 0, y: 0, z: 0 }); // Offset for molecule 2
  const [molecule2Rotation, setMolecule2Rotation] = useState({ x: 0, y: 0, z: 0 }); // Euler angles in degrees

  const [optimizationResult, setOptimizationResult] = useState(null);
  const [activeView, setActiveView] = useState("original"); // original, optimized
  const [optimizationType, setOptimizationType] = useState("classical"); // classical or quantum
  const [classicalParams, setClassicalParams] = useState({ ...defaultClassicalParams });
  const [quantumParams, setQuantumParams] = useState({ ...defaultQuantumParams });
  const [showAdvancedParams, setShowAdvancedParams] = useState(false);
  const [isHowToUseVisible, setIsHowToUseVisible] = useState(false);
  const [isTheoryVisible, setIsTheoryVisible] = useState(false);
  const [isAboutUsVisible, setIsAboutUsVisible] = useState(false);
  const [isDragActive, setIsDragActive] = useState(false);
  const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false);
  const [isSubscribeLoading, setIsSubscribeLoading] = useState(false);
  const [isCancelLoading, setIsCancelLoading] = useState(false);
  const [isOptimizeLoading, setIsOptimizeLoading] = useState(false);
  const [isMobile, setIsMobile] = useState(false);
  const [serverHealth, setServerHealth] = useState(null);
  const [isCheckingHealth, setIsCheckingHealth] = useState(false);
  const [serverHealthDetails, setServerHealthDetails] = useState(null);
  // FIXED: Moved content states before conditional return
  const [howToUseContent, setHowToUseContent] = useState("");
  const [theoryContent, setTheoryContent] = useState("");
  const [aboutUsContent, setAboutUsContent] = useState("");
  const [applyDefaultOffset, setApplyDefaultOffset] = useState(true);

  const [taskStatus, setTaskStatus] = useState(null);
  const [statusPollInterval, setStatusPollInterval] = useState(null);

  const countAtoms = (molecule) => {
    if (!molecule) return 0;
    return molecule.length;
  };

  // Check if molecule size is within quantum limits
  const validateQuantumMoleculeSize = (molecule1, molecule2, isInteractionMode) => {
    if (isInteractionMode) {
      // Check total atoms for interaction mode
      const totalAtoms = countAtoms(molecule1) + countAtoms(molecule2);
      if (totalAtoms > QUANTUM_SIZE_LIMITS.interactionTotal) {
        return {
          valid: false,
          message: `Quantum optimization for molecular interactions is limited to ${QUANTUM_SIZE_LIMITS.interactionTotal} total atoms (current: ${totalAtoms}). Please use classical optimization for larger systems.`
        };
      }

      // Check individual molecule sizes
      if (countAtoms(molecule1) > QUANTUM_SIZE_LIMITS.interactionPerMolecule) {
        return {
          valid: false,
          message: `Molecule 1 exceeds the limit of ${QUANTUM_SIZE_LIMITS.interactionPerMolecule} atoms for quantum interaction optimization (current: ${countAtoms(molecule1)}).`
        };
      }

      if (countAtoms(molecule2) > QUANTUM_SIZE_LIMITS.interactionPerMolecule) {
        return {
          valid: false,
          message: `Molecule 2 exceeds the limit of ${QUANTUM_SIZE_LIMITS.interactionPerMolecule} atoms for quantum interaction optimization (current: ${countAtoms(molecule2)}).`
        };
      }
    } else {
      // Single molecule mode - check the active molecule
      const activeMolecule = molecule1 || molecule2;
      const atomCount = countAtoms(activeMolecule);

      if (atomCount > QUANTUM_SIZE_LIMITS.singleMolecule) {
        return {
          valid: false,
          message: `Quantum optimization is limited to ${QUANTUM_SIZE_LIMITS.singleMolecule} atoms per molecule (current: ${atomCount}). Please use classical optimization for larger molecules.`
        };
      }
    }

    return { valid: true };
  };

  const apiBaseUrl = "/api";

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

  useEffect(() => {
    // Load theory documentation
    fetch(`${process.env.PUBLIC_URL}/theory.md`)
      .then(response => {
        if (!response.ok) {
          throw new Error('Theory documentation not found');
        }
        return response.text();
      })
      .then(text => setTheoryContent(text))
      .catch(error => {
        console.error("Failed to load theory:", error);
        // Fallback content
        setTheoryContent("# Theory\n\nTheory documentation is currently unavailable.");
      });
  }, []);

  useEffect(() => {
    // Load about us content
    fetch(`${process.env.PUBLIC_URL}/about-us.md`)
      .then(response => {
        if (!response.ok) {
          throw new Error('About Us content not found');
        }
        return response.text();
      })
      .then(text => setAboutUsContent(text))
      .catch(error => {
        console.error("Failed to load about us:", error);
        // Fallback content
        setAboutUsContent("# About Us\n\nInformation is currently unavailable.");
      });
  }, []);

  // Apply limits on initial load and whenever subscription status changes
  useEffect(() => {
    if (!isLoading) {
      applyIterationLimits(isSubscribed);
    }
  }, [isLoading, isSubscribed]);

  useEffect(() => {
    // Cleanup function for component unmount
    return () => {
      if (statusPollInterval) {
        clearInterval(statusPollInterval);
      }
    };
  }, [statusPollInterval]);

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
    setIsTheoryVisible(false);
    setIsAboutUsVisible(false);
    setIsMobileMenuOpen(false);
  };

  const handleShowTheory = () => {
    setIsHowToUseVisible(false);
    setIsTheoryVisible(true);
    setIsAboutUsVisible(false);
    setIsMobileMenuOpen(false);
  };

  const handleShowAboutUs = () => {
    setIsHowToUseVisible(false);
    setIsTheoryVisible(false);
    setIsAboutUsVisible(true);
    setIsMobileMenuOpen(false);
  };

  const handleClosePopup = () => {
    setIsHowToUseVisible(false);
    setIsTheoryVisible(false);
    setIsAboutUsVisible(false);
  };

  const toggleMobileMenu = () => {
    setIsMobileMenuOpen(!isMobileMenuOpen);
  };

  const getMoleculeName = (moleculeData) => {
    if (!moleculeData) return null;

    // Check if the molecule has metadata with name
    if (moleculeData.file1?.metadata?.name) {
      return moleculeData.file1.metadata.name;
    }

    // If no metadata name, check for a filename property (for uploaded files)
    if (moleculeData.filename) {
      return moleculeData.filename;
    }

    // If no name is found, return a default name
    return "Unknown Molecule";
  };

  const onMoleculeRotate = (newRotation) => {
    setMolecule2Rotation(newRotation);
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

        // Store file name with the molecule data
        parsedData.filename = file.name;

        // Target the active molecule
        if (activeMolecule === 1) {
          setMolecule1Data(parsedData);
        } else {
          setMolecule2Data(parsedData);
          // Set a default offset for molecule 2 only if applyDefaultOffset is true
          if (applyDefaultOffset) {
            setMolecule2Offset({ x: 5, y: 0, z: 0 });
          } else {
            setMolecule2Offset({ x: 0, y: 0, z: 0 });
          }
        }

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

          // Store file name with the molecule data
          parsedData.filename = file.name;

          // Target the active molecule
          if (activeMolecule === 1) {
            setMolecule1Data(parsedData);
          } else {
            setMolecule2Data(parsedData);
            // Set a default offset for molecule 2 only if applyDefaultOffset is true
            if (applyDefaultOffset) {
              setMolecule2Offset({ x: 5, y: 0, z: 0 });
            } else {
              setMolecule2Offset({ x: 0, y: 0, z: 0 });
            }
          }

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
    // If optimization type is changing
    if (newType !== optimizationType) {
      // Check if there's an optimization result for the current view
      if (activeView === "optimized" && optimizationResult?.result) {
        // Check if the optimization result was generated by the new type
        const resultMatchesNewType = (
          (newType === "classical" && optimizationResult.result.metadata?.method === "classical_molecular_dynamics") ||
          (newType === "quantum" && optimizationResult.result.metadata?.method === "quantum_chemistry")
        );

        // Only reset view if the optimization result doesn't match the new type
        if (!resultMatchesNewType) {
          setActiveView("original");
        }
        // Otherwise, keep the current "optimized" view to preserve positions
      } else {
        // If not in "optimized" view or no optimization result, reset to "original"
        setActiveView("original");
      }
    }

    setOptimizationType(newType);
  };

  const handleOptimize = async () => {
    // For interaction mode, ensure both molecules are loaded
    if (interactionMode && (!molecule1Data || !molecule2Data)) {
      alert("Please upload or select both molecules for interaction optimization.");
      return;
    }

    // For single molecule mode, at least one molecule must be loaded
    if (!interactionMode && !molecule1Data && !molecule2Data) {
      alert("Please upload or select at least one molecule.");
      return;
    }

    // Get the current atoms based on the active view for subsequent optimizations
    const { molecule1, molecule2 } = getAtoms();

    // For quantum optimization, check molecule size limits
    if (optimizationType === "quantum") {
      const sizeValidation = validateQuantumMoleculeSize(molecule1, molecule2, interactionMode);
      if (!sizeValidation.valid) {
        alert(sizeValidation.message);
        return;
      }
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

      // Prepare molecule1Data for optimization request
      let requestMolecule1 = null;
      if (molecule1) {
        // Create a deep clone of the source data to avoid modifying it directly
        requestMolecule1 = {
          atoms: molecule1
        };
      }

      // Helper functions for coordinate transformations
      const calculateCenterOfMass = (atoms) => {
        if (!atoms || atoms.length === 0) return { x: 0, y: 0, z: 0 };
        const sum = atoms.reduce((acc, atom) => ({
          x: acc.x + atom.x,
          y: acc.y + atom.y,
          z: acc.z + atom.z
        }), { x: 0, y: 0, z: 0 });
        return {
          x: sum.x / atoms.length,
          y: sum.y / atoms.length,
          z: sum.z / atoms.length
        };
      };

      const applyRotation = (coords, rotation, centerOfMass) => {
        // Convert degrees to radians
        const radX = rotation.x * Math.PI / 180;
        const radY = rotation.y * Math.PI / 180;
        const radZ = rotation.z * Math.PI / 180;

        // Translate to origin (center of mass)
        const centered = [
          coords[0] - centerOfMass.x,
          coords[1] - centerOfMass.y,
          coords[2] - centerOfMass.z
        ];

        // Apply rotations (ZYX order)
        // Z-axis rotation
        let nx = centered[0] * Math.cos(radZ) - centered[1] * Math.sin(radZ);
        let ny = centered[0] * Math.sin(radZ) + centered[1] * Math.cos(radZ);
        let nz = centered[2];

        // Y-axis rotation
        let mx = nx * Math.cos(radY) + nz * Math.sin(radY);
        let my = ny;
        let mz = -nx * Math.sin(radY) + nz * Math.cos(radY);

        // X-axis rotation
        let rx = mx;
        let ry = my * Math.cos(radX) - mz * Math.sin(radX);
        let rz = my * Math.sin(radX) + mz * Math.cos(radX);

        // Translate back from origin
        return [
          rx + centerOfMass.x,
          ry + centerOfMass.y,
          rz + centerOfMass.z
        ];
      };

      // Prepare molecule2Data for optimization request
      let requestMolecule2 = null;
      let processedMolecule2 = null;

      if (molecule2) {
        // Create base molecule2 data
        requestMolecule2 = {
          atoms: molecule2
        };

        // Only apply transformations if viewing original coordinates
        if (interactionMode) {
          if (activeView === "original") {
            // Calculate center of mass for rotation
            const centerOfMass = calculateCenterOfMass(molecule2);

            // Clone and transform molecule2 atoms
            const transformedAtoms = molecule2.map(atom => {
              // Deep clone atom
              const newAtom = { ...atom };

              // Apply rotation if needed
              if (molecule2Rotation && (molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0)) {
                const coords = [atom.x, atom.y, atom.z];
                const rotated = applyRotation(coords, molecule2Rotation, centerOfMass);
                newAtom.x = rotated[0];
                newAtom.y = rotated[1];
                newAtom.z = rotated[2];
              }

              // Apply offset
              newAtom.x += molecule2Offset.x;
              newAtom.y += molecule2Offset.y;
              newAtom.z += molecule2Offset.z;

              return newAtom;
            });

            // Create transformed molecule
            processedMolecule2 = {
              atoms: transformedAtoms
            };
          } else {
            // For optimized view, use coordinates as-is without additional transformations
            processedMolecule2 = {
              atoms: molecule2.map(atom => ({ ...atom }))
            };
          }
        }
      }

      // Payload now uses conditional transformation based on active view
      const payload = {
        molecule1: requestMolecule1,
        molecule2: interactionMode ? processedMolecule2 : null,
        optimization_type: optimizationType,
        optimization_params: optimizationParams,
        interaction_mode: interactionMode,
        // Keep these for backward compatibility with backend
        // Pass zero values when using optimized coordinates to avoid backend re-applying them
        molecule2_offset: interactionMode ? (activeView === "original" ? molecule2Offset : { x: 0, y: 0, z: 0 }) : null,
        molecule2_rotation: interactionMode ? (activeView === "original" ? molecule2Rotation : { x: 0, y: 0, z: 0 }) : null
      };

      console.log('Optimization payload:', JSON.stringify(payload, null, 2));
      console.log("Attempting request to:", `${apiBaseUrl}/optimize-molecule`);

      // Token is handled by axios interceptor in AuthContext if user is authenticated
      const headers = {
        'Content-Type': 'application/json',
        'Accept': 'application/json'
      };

      // Add Authorization token only if authenticated
      if (isAuthenticated && token) {
        headers['Authorization'] = `Bearer ${token}`;
      }

      const response = await axios({
        method: 'post',
        url: `${apiBaseUrl}/optimize-molecule`,
        data: payload,
        headers: headers
      });

      console.log("Response received:", response);

      // Handle task-based response (new asynchronous processing flow)
      if (response.data.success) {
        const taskId = response.data.task_id;
        const estimatedTime = response.data.estimated_time || 60; // seconds

        // Show task status to user
        setTaskStatus({
          id: taskId,
          status: 'pending',
          message: `Optimization in queue... (Estimated completion: ${Math.ceil(estimatedTime / 60)} min)`,
          progress: 0
        });

        // Begin polling for task status
        pollTaskStatus(taskId);
      } else {
        alert("Optimization failed. " + (response.data.error || ""));
        setIsOptimizeLoading(false);
      }
    } catch (error) {
      console.error("Error optimizing molecule:", error);
      console.error("Request error details:", {
        message: error.message,
        response: error.response,
        request: error.request
      });
      alert("Error optimizing molecule: " + (error.response?.data?.error || error.message));
      setIsOptimizeLoading(false);
    }
  };

  // Helper function to poll task status
  const pollTaskStatus = async (taskId) => {
    try {
      // Using interval to poll status every 3 seconds
      const intervalId = setInterval(async () => {
        try {
          const statusResponse = await axios({
            method: 'get',
            url: `${apiBaseUrl}/optimization-status/${taskId}`,
            headers: {
              'Accept': 'application/json',
              'Authorization': isAuthenticated && token ? `Bearer ${token}` : undefined
            }
          });

          const taskData = statusResponse.data;
          console.log(`Task ${taskId} status update:`, taskData);

          // Update status display
          setTaskStatus(prevStatus => ({
            ...prevStatus,
            status: taskData.status,
            message: getStatusMessage(taskData.status),
            progress: taskData.progress || calculateProgressFromStatus(taskData.status)
          }));

          // Handle task completion or failure
          if (taskData.status === 'completed' && taskData.result) {
            clearInterval(intervalId);

            // Store the molecule2Offset and molecule2Rotation with the result
            const result = {
              ...taskData.result,
              molecule2Offset: molecule2Offset,
              molecule2Rotation: molecule2Rotation
            };

            setOptimizationResult({
              result: result,
              success: true
            });

            setActiveView("optimized");
            setIsOptimizeLoading(false);

            // Disable positioning mode after optimization
            if (positioningMode) {
              setPositioningMode(false);
            }
          } else if (taskData.status === 'failed' || taskData.status === 'cancelled') {
            clearInterval(intervalId);
            alert("Optimization failed: " + (taskData.error || "Unknown error"));
            setIsOptimizeLoading(false);
          }
        } catch (pollError) {
          console.error("Error polling task status:", pollError);
          // Don't stop polling on temporary errors - backend might still be processing
        }
      }, 3000);

      // Store interval ID to clear on component unmount
      setStatusPollInterval(intervalId);

    } catch (error) {
      console.error("Error setting up task polling:", error);
      alert("Failed to track optimization status: " + error.message);
      setIsOptimizeLoading(false);
    }
  };

  // Helper function to determine status message
  const getStatusMessage = (status) => {
    switch (status) {
      case 'pending': return 'Task queued - waiting for processing...';
      case 'running': return 'Optimization in progress...';
      case 'completed': return 'Optimization completed successfully!';
      case 'failed': return 'Optimization failed. Please try again.';
      case 'cancelled': return 'Optimization was cancelled.';
      default: return 'Checking optimization status...';
    }
  };

  // Helper function to estimate progress based on status
  const calculateProgressFromStatus = (status) => {
    switch (status) {
      case 'pending': return 0.1;
      case 'running': return 0.5;
      case 'completed': return 1.0;
      default: return 0;
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
      const testMoleculeData = { ...TEST_MOLECULES[moleculeKey] };

      // Add a molecule name reference - though this is redundant since the test molecules
      // already have metadata.name, this ensures we have a consistent approach
      if (!testMoleculeData.filename) {
        testMoleculeData.filename = testMoleculeData.file1?.metadata?.name || moleculeKey;
      }

      if (activeMolecule === 1) {
        setMolecule1Data(testMoleculeData);
      } else {
        setMolecule2Data(testMoleculeData);
        // Set a default offset for molecule 2 only if applyDefaultOffset is true
        if (applyDefaultOffset) {
          setMolecule2Offset({ x: 5, y: 0, z: 0 });
        } else {
          setMolecule2Offset({ x: 0, y: 0, z: 0 });
        }
      }
      setOptimizationResult(null);
      setActiveView("original");
    }
  };

  const handleDownload = () => {
    if (!optimizationResult) {
      alert("No optimization results available to download.");
      return;
    }

    let downloadData;

    if (interactionMode && optimizationResult.result.molecule1_optimized_atoms && optimizationResult.result.molecule2_optimized_atoms) {
      // For interaction mode, create a combined file with both molecules
      downloadData = {
        file1: {
          atoms: optimizationResult.result.molecule1_optimized_atoms,
          metadata: {
            ...optimizationResult.result.metadata,
            molecule: "Molecule 1"
          }
        },
        file2: {
          atoms: optimizationResult.result.molecule2_optimized_atoms,
          metadata: {
            ...optimizationResult.result.metadata,
            molecule: "Molecule 2"
          }
        }
      };
    } else {
      // For single molecule, use the original format
      downloadData = {
        file1: {
          atoms: optimizationResult.result.optimized_atoms,
          metadata: optimizationResult.result.metadata
        }
      };
    }

    const filename = interactionMode ?
      `${optimizationType}_optimized_interaction.json` :
      `${optimizationType}_optimized_molecule.json`;

    const blob = new Blob([JSON.stringify(downloadData, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);

    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    link.click();

    URL.revokeObjectURL(url);
  };

  // Get atoms arrays based on active view
  const getAtoms = () => {
    // For molecule1, determine the correct atoms array based on active view
    const atoms1 = (() => {
      if (activeView === "original") {
        // Use original atoms from molecule1Data
        if (molecule1Data?.file1?.atoms) {
          return molecule1Data.file1.atoms;
        } else if (molecule1Data?.atoms) {
          return molecule1Data.atoms;
        }
      } else if (activeView === "optimized" && optimizationResult?.result) {
        // Use optimized atoms from previous optimization
        if (interactionMode && optimizationResult.result.molecule1_optimized_atoms) {
          return optimizationResult.result.molecule1_optimized_atoms;
        } else if (!interactionMode && optimizationResult.result.optimized_atoms) {
          // For backward compatibility with single molecule optimization
          return optimizationResult.result.optimized_atoms;
        }
      }
      return null;
    })();

    // For molecule2, determine the correct atoms array based on active view
    const atoms2 = (() => {
      if (activeView === "original") {
        // For "original" view, use the source molecule2Data with any user-applied
        // offset and rotation from the interaction mode
        if (molecule2Data?.file1?.atoms) {
          return molecule2Data.file1.atoms;
        } else if (molecule2Data?.atoms) {
          return molecule2Data.atoms;
        }
      } else if (activeView === "optimized" && optimizationResult?.result && interactionMode) {
        // For "optimized" view in interaction mode, use the optimized atoms for molecule2
        return optimizationResult.result.molecule2_optimized_atoms;
      }
      return null;
    })();

    return { molecule1: atoms1, molecule2: atoms2 };
  };

  // Toggle auth modal visibility
  const toggleAuthModal = () => {
    setShowAuthModal(!showAuthModal);
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
          <div>Loading...</div>
        </div>
      </div>
    );
  }

  // Main App render - works for both authenticated and guest users
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
              marginRight: '5px',
              background: "linear-gradient(145deg, rgba(56, 189, 248, 0.9), rgba(37, 99, 235, 0.9))"
            }}
            className="float howToUseButton"
          >
            <span style={styles.howToUseIcon}><Icons.book /></span>
            How To Use
          </button>

          <button
            onClick={handleShowTheory}
            style={{
              ...styles.howToUseButton,
              position: 'static',
              marginRight: '5px',
              background: "linear-gradient(145deg, rgba(16, 185, 129, 0.9), rgba(5, 150, 105, 0.9))"
            }}
            className="float theoryButton"
          >
            <span style={styles.howToUseIcon}><Icons.quantum /></span>
            Theory
          </button>

          <button
            onClick={handleShowAboutUs}
            style={{
              ...styles.howToUseButton,
              position: 'static',
              marginRight: '5px',
              background: "linear-gradient(145deg, rgba(99, 102, 241, 0.9), rgba(79, 70, 229, 0.9))"
            }}
            className="float aboutUsButton"
          >
            <span style={styles.howToUseIcon}><Icons.info /></span>
            About Us
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
            className="health-check-button float"
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

          {/* Conditional login/logout buttons */}
          {isAuthenticated ? (
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
                gap: "5px",
                animation: "float 8s ease-in-out infinite"
              }}
            >
              <Icons.close />
              Logout
            </button>
          ) : (
            <button
              onClick={toggleAuthModal}
              style={{
                ...styles.button,
                background: "linear-gradient(145deg, rgba(56, 189, 248, 0.9), rgba(37, 99, 235, 0.9))",
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
                animation: "float 8s ease-in-out infinite"
              }}
            >
              <Icons.verified />
              Login / Register
            </button>
          )}

          {isSubscribed && (
            <button
              onClick={handleCancelSubscription}
              disabled={isCancelLoading}
              style={{
                ...styles.cancelSubscriptionButton,
                opacity: isCancelLoading ? 0.7 : 1,
                cursor: isCancelLoading ? "not-allowed" : "pointer",
                position: 'static',
                animation: "float 8s ease-in-out infinite"
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
            Advanced computational chemistry tools for structure and drug-target energy optimization
          </p>
          <p style={styles.headerSubtitle} className="app-subtitle">
            -
          </p>
          <p style={styles.headerSubtitle} className="app-subtitle">
            Quantum energy optimization currently capped at 30 atoms total due to computational complexity. Will add more powerful servers and control of only quantum optimizing a subset of the system soon.
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
                  isAuthenticated={isAuthenticated}
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

          {/* Molecule Selection UI */}
          <div style={{
            backgroundColor: "rgba(15, 23, 42, 0.7)",
            backdropFilter: "blur(12px)",
            WebkitBackdropFilter: "blur(12px)",
            borderRadius: "16px",
            padding: "16px",
            marginBottom: "24px",
            border: "1px solid rgba(255, 255, 255, 0.1)",
            boxShadow: "0 15px 25px rgba(0, 0, 0, 0.25)",
          }}>
            <h3 style={{
              fontSize: "1rem",
              fontWeight: "600",
              marginBottom: "16px",
              color: "#f0f4f8",
              display: "flex",
              alignItems: "center",
              justifyContent: "center",
              gap: "8px"
            }}>
              <span style={{ color: "#38bdf8" }}><Icons.molecule /></span>
              Molecule Selection
            </h3>

            <div style={{
              display: "flex",
              justifyContent: "center",
              gap: "16px",
              marginBottom: "8px",
              flexWrap: "wrap"
            }} className={isMobile ? 'mobile-stack' : ''}>
              <div
                style={{
                  display: "flex",
                  flexDirection: "column",
                  width: isMobile ? "100%" : "180px",
                  marginBottom: isMobile ? "8px" : 0
                }}
              >
                <button
                  onClick={() => setActiveMolecule(1)}
                  style={{
                    padding: "8px 16px",
                    backgroundColor: activeMolecule === 1 ? "rgba(56, 189, 248, 0.2)" : "rgba(15, 23, 42, 0.7)",
                    color: activeMolecule === 1 ? "#38bdf8" : "#94a3b8",
                    borderRadius: "8px 8px 0 0",
                    fontWeight: activeMolecule === 1 ? "600" : "400",
                    cursor: "pointer",
                    transition: "all 0.3s cubic-bezier(0.4, 0, 0.2, 1)",
                    border: `1px solid ${activeMolecule === 1 ? "rgba(56, 189, 248, 0.3)" : "transparent"}`,
                    borderBottom: "none",
                    boxShadow: activeMolecule === 1 ? "0 4px 6px rgba(0, 0, 0, 0.1)" : "none",
                    display: "flex",
                    alignItems: "center",
                    gap: "6px",
                    justifyContent: "center",
                  }}
                  className={`molecule-selector-button ${isMobile ? 'mobile-full-width' : ''}`}
                >
                  <Icons.molecule />
                  Molecule 1 {molecule1Data ? " (Loaded)" : " (Empty)"}
                </button>

                {molecule1Data && (
                  <div style={{
                    backgroundColor: "rgba(56, 189, 248, 0.1)",
                    padding: "8px",
                    borderRadius: "0 0 8px 8px",
                    fontSize: "0.85rem",
                    textAlign: "center",
                    border: "1px solid rgba(56, 189, 248, 0.3)",
                    borderTop: "none",
                    display: "flex",
                    flexDirection: "column",
                    gap: "6px"
                  }}>
                    <div style={{ fontWeight: "500", color: "#38bdf8", whiteSpace: "nowrap", overflow: "hidden", textOverflow: "ellipsis" }}>
                      {getMoleculeName(molecule1Data)}
                    </div>
                    <button
                      onClick={() => setMolecule1Data(null)}
                      style={{
                        backgroundColor: "rgba(244, 63, 94, 0.1)",
                        color: "#f43f5e",
                        border: "1px solid rgba(244, 63, 94, 0.3)",
                        borderRadius: "4px",
                        padding: "4px 8px",
                        fontSize: "0.75rem",
                        cursor: "pointer",
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        gap: "4px",
                        transition: "all 0.2s ease"
                      }}
                    >
                      <Icons.close /> Clear
                    </button>
                  </div>
                )}
              </div>

              <div
                style={{
                  display: "flex",
                  flexDirection: "column",
                  width: isMobile ? "100%" : "180px",
                }}
              >
                <button
                  onClick={() => setActiveMolecule(2)}
                  style={{
                    padding: "8px 16px",
                    backgroundColor: activeMolecule === 2 ? "rgba(16, 185, 129, 0.2)" : "rgba(15, 23, 42, 0.7)",
                    color: activeMolecule === 2 ? "#10b981" : "#94a3b8",
                    borderRadius: "8px 8px 0 0",
                    fontWeight: activeMolecule === 2 ? "600" : "400",
                    cursor: "pointer",
                    transition: "all 0.3s cubic-bezier(0.4, 0, 0.2, 1)",
                    border: `1px solid ${activeMolecule === 2 ? "rgba(16, 185, 129, 0.3)" : "transparent"}`,
                    borderBottom: "none",
                    boxShadow: activeMolecule === 2 ? "0 4px 6px rgba(0, 0, 0, 0.1)" : "none",
                    display: "flex",
                    alignItems: "center",
                    gap: "6px",
                    justifyContent: "center",
                  }}
                  className={`molecule-selector-button ${isMobile ? 'mobile-full-width' : ''}`}
                >
                  <Icons.molecule />
                  Molecule 2 {molecule2Data ? " (Loaded)" : " (Empty)"}
                </button>

                {molecule2Data && (
                  <div style={{
                    backgroundColor: "rgba(16, 185, 129, 0.1)",
                    padding: "8px",
                    borderRadius: "0 0 8px 8px",
                    fontSize: "0.85rem",
                    textAlign: "center",
                    border: "1px solid rgba(16, 185, 129, 0.3)",
                    borderTop: "none",
                    display: "flex",
                    flexDirection: "column",
                    gap: "6px"
                  }}>
                    <div style={{ fontWeight: "500", color: "#10b981", whiteSpace: "nowrap", overflow: "hidden", textOverflow: "ellipsis" }}>
                      {getMoleculeName(molecule2Data)}
                    </div>
                    <button
                      onClick={() => setMolecule2Data(null)}
                      style={{
                        backgroundColor: "rgba(244, 63, 94, 0.1)",
                        color: "#f43f5e",
                        border: "1px solid rgba(244, 63, 94, 0.3)",
                        borderRadius: "4px",
                        padding: "4px 8px",
                        fontSize: "0.75rem",
                        cursor: "pointer",
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        gap: "4px",
                        transition: "all 0.2s ease"
                      }}
                    >
                      <Icons.close /> Clear
                    </button>
                  </div>
                )}
              </div>
            </div>

            {/* Clear both molecules button */}
            {(molecule1Data || molecule2Data) && (
              <div style={{
                display: "flex",
                justifyContent: "center",
                marginTop: "12px"
              }}>
                <button
                  onClick={() => {
                    setMolecule1Data(null);
                    setMolecule2Data(null);
                    setOptimizationResult(null);
                    setActiveView("original");
                  }}
                  style={{
                    backgroundColor: "rgba(244, 63, 94, 0.1)",
                    color: "#f43f5e",
                    border: "1px solid rgba(244, 63, 94, 0.3)",
                    borderRadius: "6px",
                    padding: "6px 12px",
                    fontSize: "0.875rem",
                    cursor: "pointer",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    gap: "6px",
                    transition: "all 0.2s ease"
                  }}
                >
                  <Icons.close /> Clear All Molecules
                </button>
              </div>
            )}

            <div style={{
              display: "flex",
              justifyContent: "center",
              marginTop: "16px",
            }}>
              <label style={{
                display: "flex",
                alignItems: "center",
                cursor: "pointer"
              }}>
                <input
                  type="checkbox"
                  checked={interactionMode}
                  onChange={(e) => setInteractionMode(e.target.checked)}
                  style={{
                    marginRight: "4px",
                    cursor: "pointer",
                    accentColor: "#38bdf8"
                  }}
                />
                <span>Optimize Molecular Interaction</span>
              </label>
            </div>
            {/* Default Offset Checkbox */}
            <div style={{
              display: "flex",
              justifyContent: "center",
              marginTop: "8px",
            }}>
              <label style={{
                display: "flex",
                alignItems: "center",
                cursor: "pointer"
              }}>
                <input
                  type="checkbox"
                  checked={applyDefaultOffset}
                  onChange={(e) => setApplyDefaultOffset(e.target.checked)}
                  style={{
                    marginRight: "4px",
                    cursor: "pointer",
                    accentColor: "#38bdf8"
                  }}
                />
                <span>Apply default offset (5Å) when adding molecule 2</span>
              </label>
            </div>
          </div>

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
                onClick={() => handleTestMoleculeSelect('aceticAcid')}
                style={styles.testMoleculeButton}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
              >
                <span style={styles.testMoleculeIcon}><Icons.molecule /></span>
                Acetic Acid (CH₃COOH)
              </button>
              <button
                onClick={() => handleTestMoleculeSelect('methanol')}
                style={styles.testMoleculeButton}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
              >
                <span style={styles.testMoleculeIcon}><Icons.molecule /></span>
                Methanol (CH₃OH)
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

          {(molecule1Data || molecule2Data) && (
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

              {/* Positioning Mode Controls */}
              {interactionMode && molecule1Data && molecule2Data && (
                <div style={{
                  backgroundColor: "rgba(15, 23, 42, 0.7)",
                  backdropFilter: "blur(12px)",
                  WebkitBackdropFilter: "blur(12px)",
                  borderRadius: "16px",
                  padding: "16px",
                  marginBottom: "16px",
                  border: `1px solid ${positioningMode ? "#38bdf8" : "rgba(255, 255, 255, 0.1)"}`,
                  boxShadow: positioningMode ? `0 0 10px rgba(56, 189, 248, 0.2)` : "0 15px 25px rgba(0, 0, 0, 0.25)",
                }}>
                  <div style={{
                    display: "flex",
                    flexDirection: isMobile ? "column" : "row",
                    alignItems: "center",
                    justifyContent: "space-between",
                    gap: "16px",
                  }}>
                    <div style={{
                      display: "flex",
                      alignItems: "center",
                      gap: "8px",
                    }}>
                      <button
                        onClick={() => setPositioningMode(!positioningMode)}
                        style={{
                          backgroundColor: positioningMode ? "rgba(56, 189, 248, 0.2)" : "rgba(15, 23, 42, 0.7)",
                          color: positioningMode ? "#38bdf8" : "#94a3b8",
                          border: `1px solid ${positioningMode ? "rgba(56, 189, 248, 0.3)" : "transparent"}`,
                          padding: "4px 16px",
                          borderRadius: "8px",
                          display: "flex",
                          alignItems: "center",
                          gap: "6px",
                          cursor: "pointer",
                          transition: "all 0.3s cubic-bezier(0.4, 0, 0.2, 1)",
                          fontSize: "16px", // Increased font size from default
                          fontWeight: "500" // Added for better visibility
                        }}
                      >
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                          <path d="M5 9l4-4 4 4M9 5v14M19 15l-4 4-4-4M15 19V5" />
                        </svg>
                        {positioningMode ? "Disable Positioning Mode" : "Enable Positioning Mode"}
                      </button>
                    </div>

                    {positioningMode && (
                      <div style={{
                        display: "flex",
                        gap: "8px",
                        alignItems: "center",
                        flexWrap: isMobile ? "wrap" : "nowrap",
                        justifyContent: isMobile ? "center" : "flex-end",
                        width: isMobile ? "100%" : "auto",
                      }}>
                        <button
                          onClick={() => setMolecule2Offset({ x: 0, y: 0, z: 0 })}
                          style={{
                            backgroundColor: "rgba(15, 23, 42, 0.7)",
                            color: "#94a3b8",
                            border: "1px solid rgba(255, 255, 255, 0.1)",
                            padding: "4px 8px",
                            borderRadius: "8px",
                            cursor: "pointer",
                          }}
                        >
                          Reset Position
                        </button>

                        <div style={{
                          display: "flex",
                          alignItems: "center",
                          gap: "4px",
                        }}>
                          <span>X:</span>
                          <input
                            type="number"
                            value={molecule2Offset.x.toFixed(2)}
                            onChange={(e) => setMolecule2Offset({
                              ...molecule2Offset,
                              x: parseFloat(e.target.value)
                            })}
                            style={{
                              width: "70px",
                              backgroundColor: "rgba(15, 23, 42, 0.7)",
                              color: "#f0f4f8",
                              border: "1px solid rgba(255, 255, 255, 0.1)",
                              borderRadius: "4px",
                              padding: "4px 4px",
                            }}
                            step="0.5"
                          />
                        </div>

                        <div style={{
                          display: "flex",
                          alignItems: "center",
                          gap: "4px",
                        }}>
                          <span>Y:</span>
                          <input
                            type="number"
                            value={molecule2Offset.y.toFixed(2)}
                            onChange={(e) => setMolecule2Offset({
                              ...molecule2Offset,
                              y: parseFloat(e.target.value)
                            })}
                            style={{
                              width: "70px",
                              backgroundColor: "rgba(15, 23, 42, 0.7)",
                              color: "#f0f4f8",
                              border: "1px solid rgba(255, 255, 255, 0.1)",
                              borderRadius: "4px",
                              padding: "4px 4px",
                            }}
                            step="0.5"
                          />
                        </div>

                        <div style={{
                          display: "flex",
                          alignItems: "center",
                          gap: "4px",
                        }}>
                          <span>Z:</span>
                          <input
                            type="number"
                            value={molecule2Offset.z.toFixed(2)}
                            onChange={(e) => setMolecule2Offset({
                              ...molecule2Offset,
                              z: parseFloat(e.target.value)
                            })}
                            style={{
                              width: "70px",
                              backgroundColor: "rgba(15, 23, 42, 0.7)",
                              color: "#f0f4f8",
                              border: "1px solid rgba(255, 255, 255, 0.1)",
                              borderRadius: "4px",
                              padding: "4px 4px",
                            }}
                            step="0.5"
                          />
                        </div>
                      </div>
                    )}
                  </div>

                  {/* Rotation Controls */}
                  {positioningMode && (
                    <div style={{
                      display: "flex",
                      flexDirection: "column",
                      gap: "8px",
                      padding: "8px",
                      backgroundColor: "rgba(15, 23, 42, 0.7)",
                      borderRadius: "8px",
                      marginTop: "8px",
                    }}>
                      <div style={{ fontWeight: "600", marginBottom: "4px" }}>Rotation (degrees):</div>
                      <div style={{
                        display: "flex",
                        gap: "8px",
                        flexWrap: isMobile ? "wrap" : "nowrap",
                        justifyContent: isMobile ? "center" : "flex-start",
                      }}>
                        <div style={{
                          display: "flex",
                          alignItems: "center",
                          gap: "4px",
                        }}>
                          <span>RX:</span>
                          <input
                            type="number"
                            value={molecule2Rotation.x}
                            onChange={(e) => setMolecule2Rotation({
                              ...molecule2Rotation,
                              x: parseFloat(e.target.value) % 360
                            })}
                            style={{
                              width: "70px",
                              backgroundColor: "rgba(15, 23, 42, 0.7)",
                              color: "#f0f4f8",
                              border: "1px solid rgba(255, 255, 255, 0.1)",
                              borderRadius: "4px",
                              padding: "4px 4px",
                            }}
                            step="15"
                          />
                        </div>

                        <div style={{
                          display: "flex",
                          alignItems: "center",
                          gap: "4px",
                        }}>
                          <span>RY:</span>
                          <input
                            type="number"
                            value={molecule2Rotation.y}
                            onChange={(e) => setMolecule2Rotation({
                              ...molecule2Rotation,
                              y: parseFloat(e.target.value) % 360
                            })}
                            style={{
                              width: "70px",
                              backgroundColor: "rgba(15, 23, 42, 0.7)",
                              color: "#f0f4f8",
                              border: "1px solid rgba(255, 255, 255, 0.1)",
                              borderRadius: "4px",
                              padding: "4px 4px",
                            }}
                            step="15"
                          />
                        </div>

                        <div style={{
                          display: "flex",
                          alignItems: "center",
                          gap: "4px",
                        }}>
                          <span>RZ:</span>
                          <input
                            type="number"
                            value={molecule2Rotation.z}
                            onChange={(e) => setMolecule2Rotation({
                              ...molecule2Rotation,
                              z: parseFloat(e.target.value) % 360
                            })}
                            style={{
                              width: "70px",
                              backgroundColor: "rgba(15, 23, 42, 0.7)",
                              color: "#f0f4f8",
                              border: "1px solid rgba(255, 255, 255, 0.1)",
                              borderRadius: "4px",
                              padding: "4px 4px",
                            }}
                            step="15"
                          />
                        </div>

                        <button
                          onClick={() => setMolecule2Rotation({ x: 0, y: 0, z: 0 })}
                          style={{
                            backgroundColor: "rgba(15, 23, 42, 0.7)",
                            color: "#94a3b8",
                            border: "1px solid rgba(255, 255, 255, 0.1)",
                            padding: "4px 8px",
                            borderRadius: "8px",
                            cursor: "pointer",
                          }}
                        >
                          Reset Rotation
                        </button>
                      </div>
                    </div>
                  )}

                  {positioningMode && (
                    <div style={{
                      backgroundColor: "rgba(56, 189, 248, 0.1)",
                      marginTop: "8px",
                      padding: "8px",
                      borderRadius: "8px",
                      fontSize: "0.85rem",
                      color: "#94a3b8",
                      border: "1px solid rgba(56, 189, 248, 0.2)",
                    }}>

                      <p style={{ margin: 0 }}>
                        <strong>Controls:</strong> Use keyboard to position Molecule 2: <kbd style={{
                          backgroundColor: "rgba(255, 255, 255, 0.1)",
                          padding: "2px 5px",
                          borderRadius: "3px",
                          fontSize: "0.8rem"
                        }}>←→</kbd> (X-axis), <kbd style={{
                          backgroundColor: "rgba(255, 255, 255, 0.1)",
                          padding: "2px 5px",
                          borderRadius: "3px",
                          fontSize: "0.8rem"
                        }}>↑↓</kbd> (Y-axis), <kbd style={{
                          backgroundColor: "rgba(255, 255, 255, 0.1)",
                          padding: "2px 5px",
                          borderRadius: "3px",
                          fontSize: "0.8rem"
                        }}>PgUp/PgDn</kbd> (Z-axis).
                        Hold <kbd style={{
                          backgroundColor: "rgba(255, 255, 255, 0.1)",
                          padding: "2px 5px",
                          borderRadius: "3px",
                          fontSize: "0.8rem"
                        }}>Shift</kbd> + arrows for rotation.
                      </p>

                    </div>
                  )}
                </div>
              )}


              {/* Visualization */}
              <div style={styles.visualizationContainer} className={`glass ${isMobile ? 'mobile-smaller-padding' : ''}`}>
                <div style={styles.visualizationHeader} className={isMobile ? 'mobile-stack' : ''}>
                  <div style={styles.visualizationTitle}>
                    <span style={styles.visualizationIcon}><Icons.molecule /></span>
                    {activeView === "original" ?
                      (interactionMode ? "Original Molecules" : "Original Structure") :
                      (optimizationType === "classical" ? "Classical" : "Quantum") +
                      (interactionMode ? " Optimized Interaction" : " Optimized Structure")}
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

                <MoleculeViewer
                  atoms={getAtoms()}
                  isMobile={isMobile}
                  positioningMode={positioningMode && interactionMode}
                  onMoleculeMove={setMolecule2Offset}
                  molecule2Offset={activeView === "original" ? molecule2Offset : { x: 0, y: 0, z: 0 }}
                  molecule2Rotation={activeView === "original" ? molecule2Rotation : { x: 0, y: 0, z: 0 }}
                  onMoleculeRotate={onMoleculeRotate}
                  molecule1Name={getMoleculeName(molecule1Data)}
                  molecule2Name={getMoleculeName(molecule2Data)}
                />
                {interactionMode && molecule1Data && molecule2Data && (
                  <div style={{
                    display: "flex",
                    justifyContent: "center",
                    marginTop: "8px",
                    padding: "4px",
                    backgroundColor: "rgba(15, 23, 42, 0.7)",
                    borderRadius: "8px",
                    fontSize: "0.85rem",
                  }}>
                    <span>Molecule 1: {molecule1Data.file1?.atoms?.length || molecule1Data.atoms?.length || 0} atoms</span>
                    <span style={{ margin: "0 10px" }}>•</span>
                    <span>Molecule 2: {molecule2Data.file1?.atoms?.length || molecule2Data.atoms?.length || 0} atoms</span>
                    {(molecule2Offset.x !== 0 || molecule2Offset.y !== 0 || molecule2Offset.z !== 0 ||
                      molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0) ? (
                      <>
                        <span style={{ margin: "0 10px" }}>•</span>
                        <span>Offset: ({molecule2Offset.x.toFixed(1)}, {molecule2Offset.y.toFixed(1)}, {molecule2Offset.z.toFixed(1)})</span>
                        {(molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0) && (
                          <>
                            <span style={{ margin: "0 10px" }}>•</span>
                            <span>Rotation: ({molecule2Rotation.x}°, {molecule2Rotation.y}°, {molecule2Rotation.z}°)</span>
                          </>
                        )}
                      </>
                    ) : null}
                  </div>
                )}
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

              {taskStatus && (
                <div style={{
                  marginTop: "10px",
                  padding: "10px",
                  backgroundColor: taskStatus.status === 'failed' ? "rgba(244, 63, 94, 0.1)" :
                    taskStatus.status === 'completed' ? "rgba(16, 185, 129, 0.1)" : "rgba(56, 189, 248, 0.1)",
                  borderRadius: "8px",
                  border: `1px solid ${taskStatus.status === 'failed' ? "rgba(244, 63, 94, 0.3)" :
                    taskStatus.status === 'completed' ? "rgba(16, 185, 129, 0.3)" : "rgba(56, 189, 248, 0.3)"}`,
                }}>
                  <div style={{ marginBottom: "5px" }}>{taskStatus.message}</div>
                  <div style={{
                    height: "4px",
                    backgroundColor: "rgba(255, 255, 255, 0.2)",
                    borderRadius: "2px",
                    overflow: "hidden"
                  }}>
                    <div style={{
                      height: "100%",
                      width: `${taskStatus.progress * 100}%`,
                      backgroundColor: taskStatus.status === 'failed' ? "#f43f5e" :
                        taskStatus.status === 'completed' ? "#10b981" : "#38bdf8",
                      transition: "width 0.3s ease"
                    }}></div>
                  </div>
                </div>
              )}

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
            className={`mobile-nav-button ${(molecule1Data || molecule2Data) ? '' : 'active'}`}
            onClick={() => {
              setIsMobileMenuOpen(false);
              document.querySelector('.file-upload-area').scrollIntoView({ behavior: 'smooth' });
            }}
          >
            <span className="mobile-nav-icon"><Icons.upload /></span>
            Upload
          </button>
          {(molecule1Data || molecule2Data) && (
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
            How To Use
          </button>
          <button
            className="mobile-nav-button"
            onClick={handleShowTheory}
          >
            <span className="mobile-nav-icon"><Icons.quantum /></span>
            Theory
          </button>
          <button
            className="mobile-nav-button"
            onClick={handleShowAboutUs}
          >
            <span className="mobile-nav-icon"><Icons.info /></span>
            About Us
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

      {/* Authentication Modal for Login/Register */}
      {showAuthModal && (
        <div style={styles.popup} className="popup-overlay">
          <div
            style={styles.popupContent}
            className={`popup-content glass ${isMobile ? 'mobile-smaller-padding' : ''}`}
          >
            <button
              onClick={toggleAuthModal}
              style={styles.popupClose}
            >
              <Icons.close />
            </button>

            <div style={styles.popupScroll}>
              {showLoginForm ? (
                <LoginForm
                  toggleForm={() => setShowLoginForm(false)}
                  onAuthSuccess={toggleAuthModal} // Add this callback
                />
              ) : (
                <RegisterForm
                  toggleForm={() => setShowLoginForm(true)}
                  onAuthSuccess={toggleAuthModal} // Add this callback
                />
              )}
            </div>
          </div>
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

      {/* Theory Popup */}
      {isTheoryVisible && (
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
                {theoryContent}
              </ReactMarkdown>
            </div>
          </div>
        </div>
      )}

      {/* About Us Popup */}
      {isAboutUsVisible && (
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
                {aboutUsContent}
              </ReactMarkdown>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default App;