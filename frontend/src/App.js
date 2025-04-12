import React, { useState, useEffect, useContext } from "react";
import axios from "axios";
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import ReactMarkdown from "react-markdown";

// Import authentication context
import { AuthContext } from './AuthContext';
import LoginForm from './components/LoginForm';
import RegisterForm from './components/RegisterForm';

// Import constants
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
const stripePromise = loadStripe('pk_live_Tfh90MeFSg6jVjRCaMExaGug0078PAfanh');

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

  const [optimizeMolecule1, setOptimizeMolecule1] = useState(true);
  const [optimizeMolecule2, setOptimizeMolecule2] = useState(true);

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

    if (interactionMode && !optimizeMolecule1 && !optimizeMolecule2) {
      alert("At least one molecule must be selected for optimization.");
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
        molecule2_rotation: interactionMode ? (activeView === "original" ? molecule2Rotation : { x: 0, y: 0, z: 0 }) : null,
        optimize_molecule1: optimizeMolecule1,
        optimize_molecule2: optimizeMolecule2
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
      <div className="loading-container">
        <div>
          <div className="spinner" />
          <div>Loading...</div>
        </div>
      </div>
    );
  }

  // Main App render - works for both authenticated and guest users
  return (
    <div className={`app ${isMobile ? 'mobile-app' : ''}`}>
      {/* Attribution element with bouncing animation */}
      <div className="attribution">
        Built by <a
          href="https://x.com/robertkottelin"
          target="_blank"
          rel="noopener noreferrer"
          className="attribution-link"
        >
          @robertkottelin
        </a>
        <a
          href="https://x.com/robertkottelin"
          target="_blank"
          rel="noopener noreferrer"
          className="attribution-sublink"
        >
          Send me a message about improvements and bug fixes
        </a>
      </div>
      <div className="decorative-bg"></div>
      {/* Decorative animated lines for cyberpunk effect */}
      <div className="decorative-line" style={{ top: "15%", animationDelay: "0s" }}></div>
      <div className="decorative-line" style={{ top: "35%", animationDelay: "0.5s" }}></div>
      <div className="decorative-line" style={{ top: "65%", animationDelay: "1s" }}></div>
      <div className="decorative-line" style={{ top: "85%", animationDelay: "1.5s" }}></div>

      <div className="container">
        {/* Top Action Buttons - Move before header for mobile */}
        <div className={`top-buttons-container ${isMobile ? 'mobile-center' : ''}`}>
          <button
            onClick={handleShowHowToUse}
            className="how-to-use-button float"
          >
            <span className="how-to-use-icon"><Icons.book /></span>
            How To Use
          </button>

          <button
            onClick={handleShowTheory}
            className="how-to-use-button float theory-button"
          >
            <span className="how-to-use-icon"><Icons.quantum /></span>
            Theory
          </button>

          <button
            onClick={handleShowAboutUs}
            className="how-to-use-button float about-us-button"
          >
            <span className="how-to-use-icon"><Icons.info /></span>
            About Us
          </button>

          <button
            onClick={checkServerHealth}
            disabled={isCheckingHealth}
            className={`health-check-button float ${serverHealth === true ? 'health-success' :
              serverHealth === false ? 'health-error' : ''}`}
          >
            {isCheckingHealth ? (
              <>
                <span className="spin">
                  <Icons.spinner />
                </span>
                Checking...
              </>
            ) : (
              <>
                <div className="health-status-indicator">
                  {serverHealth === true && <span><Icons.checkmark /></span>}
                  {serverHealth === false && <span><Icons.warning /></span>}
                  {serverHealth === null && <span><Icons.info /></span>}
                  Server Health
                </div>

                {serverHealthDetails && (
                  <div className="health-details">
                    <div>Status: {serverHealthDetails.status}</div>
                    <div>
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
              className="button logout-button float"
            >
              <Icons.close />
              Logout
            </button>
          ) : (
            <button
              onClick={toggleAuthModal}
              className="button login-button float"
            >
              <Icons.verified />
              Login / Register
            </button>
          )}

          {isSubscribed && (
            <button
              onClick={handleCancelSubscription}
              disabled={isCancelLoading}
              className={`cancel-subscription-button float ${isCancelLoading ? 'disabled' : ''}`}
            >
              {isCancelLoading ? (
                <>
                  <span className="spin cancel-subscription-icon">
                    <Icons.spinner />
                  </span>
                  Cancelling...
                </>
              ) : (
                <>
                  <span className="cancel-subscription-icon"><Icons.cancel /></span>
                  Cancel Subscription
                </>
              )}
            </button>
          )}
        </div>

        {/* App Header - Now comes after buttons in the DOM */}
        <header className="header">
          <h1 className="header-title">Molecular Optimization System</h1>
          <p className="header-subtitle">
            Advanced computational chemistry tools for structure and drug-target energy optimization
          </p>
          <p className="header-subtitle">
            -
          </p>
          <p className="header-subtitle">
            Quantum energy optimization currently capped at 30 atoms total due to computational complexity. Will add more powerful servers and control of only quantum optimizing a subset of the system soon.
          </p>
        </header>

        {/* Subscription Form or Welcome Message */}
        <div className={`fade-in glass card ${isSubscribed ? 'card-with-glow' : ''} ${isMobile ? 'mobile-smaller-padding' : ''}`}>
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
            <div className={`welcome-message ${isMobile ? 'mobile-stack mobile-text-center' : ''}`}>
              <div className="welcome-icon">
                <Icons.verified />
              </div>
              <p>Welcome, <strong>{currentUser.email}</strong>! You have full access to all computational capabilities with your premium subscription.</p>
            </div>
          )}
        </div>

        {/* Main Optimization Section */}
        <section className="section">
          <h2 className="section-title">
            Molecule Optimization
            <span className="section-title-underline"></span>
          </h2>

          {/* Molecule Selection UI */}
          <div className="molecule-selection-container">
            <h3 className="molecule-selection-title">
              <span className="molecule-selection-icon"><Icons.molecule /></span>
              Molecule Selection
            </h3>

            <div className={`molecule-selector-group ${isMobile ? 'mobile-stack' : ''}`}>
              <div className="molecule-selector-column">
                <button
                  onClick={() => setActiveMolecule(1)}
                  className={`molecule-selector-button ${activeMolecule === 1 ? 'active-molecule-1' : ''} ${isMobile ? 'mobile-full-width' : ''}`}
                >
                  <Icons.molecule />
                  Molecule 1 {molecule1Data ? " (Loaded)" : " (Empty)"}
                </button>

                {molecule1Data && (
                  <div className="molecule-info molecule-1-info">
                    <div className="molecule-name">
                      {getMoleculeName(molecule1Data)}
                    </div>
                    <button
                      onClick={() => setMolecule1Data(null)}
                      className="clear-molecule-button"
                    >
                      <Icons.close /> Clear
                    </button>
                  </div>
                )}
              </div>

              <div className="molecule-selector-column">
                <button
                  onClick={() => setActiveMolecule(2)}
                  className={`molecule-selector-button ${activeMolecule === 2 ? 'active-molecule-2' : ''} ${isMobile ? 'mobile-full-width' : ''}`}
                >
                  <Icons.molecule />
                  Molecule 2 {molecule2Data ? " (Loaded)" : " (Empty)"}
                </button>

                {molecule2Data && (
                  <div className="molecule-info molecule-2-info">
                    <div className="molecule-name">
                      {getMoleculeName(molecule2Data)}
                    </div>
                    <button
                      onClick={() => setMolecule2Data(null)}
                      className="clear-molecule-button"
                    >
                      <Icons.close /> Clear
                    </button>
                  </div>
                )}
              </div>
            </div>

            {/* Clear both molecules button */}
            {(molecule1Data || molecule2Data) && (
              <div className="clear-all-container">
                <button
                  onClick={() => {
                    setMolecule1Data(null);
                    setMolecule2Data(null);
                    setOptimizationResult(null);
                    setActiveView("original");
                  }}
                  className="clear-all-button"
                >
                  <Icons.close /> Clear All Molecules
                </button>
              </div>
            )}

            <div className="interaction-checkbox-container">
              <label className="interaction-checkbox-label">
                <input
                  type="checkbox"
                  checked={interactionMode}
                  onChange={(e) => setInteractionMode(e.target.checked)}
                  className="interaction-checkbox"
                />
                <span>Optimize Molecular Interaction</span>
              </label>
            </div>

            {/* Add molecule optimization controls when in interaction mode */}
            {interactionMode && molecule1Data && molecule2Data && (
              <div className="molecule-optimization-options" style={{ marginTop: '10px', padding: '8px', borderTop: '1px solid rgba(255, 255, 255, 0.1)' }}>
                <div style={{ marginBottom: '5px', fontWeight: 'bold' }}>Optimization Control:</div>
                <label className="optimization-checkbox-label" style={{ display: 'block', margin: '5px 0' }}>
                  <input
                    type="checkbox"
                    checked={optimizeMolecule1}
                    onChange={(e) => setOptimizeMolecule1(e.target.checked)}
                    className="optimization-checkbox"
                    style={{ marginRight: '8px' }}
                  />
                  <span>Optimize Molecule 1 Structure</span>
                </label>
                <label className="optimization-checkbox-label" style={{ display: 'block', margin: '5px 0' }}>
                  <input
                    type="checkbox"
                    checked={optimizeMolecule2}
                    onChange={(e) => setOptimizeMolecule2(e.target.checked)}
                    className="optimization-checkbox"
                    style={{ marginRight: '8px' }}
                  />
                  <span>Optimize Molecule 2 Structure</span>
                </label>
                {!optimizeMolecule1 && !optimizeMolecule2 && (
                  <div style={{ color: '#f87171', fontSize: '0.85rem', marginTop: '5px' }}>
                    Warning: At least one molecule should be optimized
                  </div>
                )}
              </div>
            )}

            {/* Default Offset Checkbox */}
            <div className="default-offset-checkbox-container">
              <label className="default-offset-checkbox-label">
                <input
                  type="checkbox"
                  checked={applyDefaultOffset}
                  onChange={(e) => setApplyDefaultOffset(e.target.checked)}
                  className="default-offset-checkbox"
                />
                <span>Apply default offset (5Å) when adding molecule 2</span>
              </label>
            </div>
          </div>

          {/* Test Molecules Buttons */}
          <div className={`test-molecules-container ${isMobile ? 'mobile-smaller-padding' : ''}`}>
            <h3 className="test-molecules-title">
              Test Molecules
              <span className="test-molecules-title-icon"></span>
            </h3>
            <div className={`test-molecules-button-container ${isMobile ? 'mobile-stack' : ''}`}>
              <button
                onClick={() => handleTestMoleculeSelect('water')}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
              >
                <span className="test-molecule-icon"><Icons.molecule /></span>
                Water (H₂O)
              </button>
              <button
                onClick={() => handleTestMoleculeSelect('aceticAcid')}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
              >
                <span className="test-molecule-icon"><Icons.molecule /></span>
                Acetic Acid (CH₃COOH)
              </button>
              <button
                onClick={() => handleTestMoleculeSelect('methanol')}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
              >
                <span className="test-molecule-icon"><Icons.molecule /></span>
                Methanol (CH₃OH)
              </button>
              <button
                onClick={() => handleTestMoleculeSelect('ibuprofen')}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
              >
                <span className="test-molecule-icon"><Icons.molecule /></span>
                Ibuprofen (C₁₃H₁₈O₂)
              </button>
              <button
                onClick={() => handleTestMoleculeSelect('cox2BindingSite')}
                className={`test-molecule-button ${isMobile ? 'mobile-full-width' : ''}`}
              >
                <span className="test-molecule-icon"><Icons.molecule /></span>
                COX-2 Binding Site
              </button>
            </div>
          </div>

          {/* File Upload Area */}
          <div
            className={`file-upload ${isDragActive ? 'file-upload-active' : ''} ${isMobile ? 'mobile-smaller-padding' : ''}`}
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            onDrop={handleFileDrop}
          >
            <div className={`file-upload-icon ${isDragActive ? 'file-upload-active-icon' : ''}`}>
              <Icons.upload />
            </div>
            <p className="file-upload-text">
              {isMobile ? 'Upload a molecule file' : 'Drag & drop a molecule file here, or click to select a file'}
            </p>
            <p className="file-upload-subtext">
              Accepted format: JSON structured molecule data
            </p>
            <input
              type="file"
              onChange={handleFileUpload}
              className="file-upload-input"
            />
          </div>

          {(molecule1Data || molecule2Data) && (
            <div className="content-wrapper">
              {/* Optimization Method Selection */}
              <div
                className={`method-selection-container ${isMobile ? "mobile-stack" : ""}`}
              >
                <div
                  className={`method-selection-button ${optimizationType === "classical" ? 'classical' : ''} ${isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}`}
                  onClick={() => handleOptimizationTypeChange("classical")}
                >
                  <span className={`method-icon ${optimizationType === "classical" ? 'classical-icon' : ''}`}>
                    <Icons.classical />
                  </span>
                  <div className="method-title">Classical Optimization</div>
                  <div className="method-description">
                    Molecular mechanics optimization using empirical force fields. Faster calculations suitable for larger molecular systems.
                  </div>
                </div>

                <div
                  className={`method-selection-button ${optimizationType === "quantum" ? 'quantum' : ''} ${isMobile ? 'mobile-full-width' : ''}`}
                  onClick={() => handleOptimizationTypeChange("quantum")}
                >
                  <span className={`method-icon ${optimizationType === "quantum" ? 'quantum-icon' : ''}`}>
                    <Icons.quantum />
                  </span>
                  <div className="method-title">Quantum Optimization</div>
                  <div className="method-description">
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
                <div className={`positioning-container ${positioningMode ? 'positioning-active' : ''}`}>
                  <div className={`positioning-controls ${isMobile ? 'mobile-stack' : ''}`}>
                    <div className="positioning-toggle">
                      <button
                        onClick={() => setPositioningMode(!positioningMode)}
                        className={`positioning-toggle-button ${positioningMode ? 'active' : ''}`}
                      >
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                          <path d="M5 9l4-4 4 4M9 5v14M19 15l-4 4-4-4M15 19V5" />
                        </svg>
                        {positioningMode ? "Disable Positioning Mode" : "Enable Positioning Mode"}
                      </button>
                    </div>

                    {positioningMode && (
                      <div className={`position-controls ${isMobile ? 'mobile-stack mobile-center' : ''}`}>
                        <button
                          onClick={() => setMolecule2Offset({ x: 0, y: 0, z: 0 })}
                          className="position-reset-button"
                        >
                          Reset Position
                        </button>

                        <div className="position-input-group">
                          <span>X:</span>
                          <input
                            type="number"
                            value={molecule2Offset.x.toFixed(2)}
                            onChange={(e) => setMolecule2Offset({
                              ...molecule2Offset,
                              x: parseFloat(e.target.value)
                            })}
                            className="position-input"
                            step="0.5"
                          />
                        </div>

                        <div className="position-input-group">
                          <span>Y:</span>
                          <input
                            type="number"
                            value={molecule2Offset.y.toFixed(2)}
                            onChange={(e) => setMolecule2Offset({
                              ...molecule2Offset,
                              y: parseFloat(e.target.value)
                            })}
                            className="position-input"
                            step="0.5"
                          />
                        </div>

                        <div className="position-input-group">
                          <span>Z:</span>
                          <input
                            type="number"
                            value={molecule2Offset.z.toFixed(2)}
                            onChange={(e) => setMolecule2Offset({
                              ...molecule2Offset,
                              z: parseFloat(e.target.value)
                            })}
                            className="position-input"
                            step="0.5"
                          />
                        </div>
                      </div>
                    )}
                  </div>

                  {/* Rotation Controls */}
                  {positioningMode && (
                    <div className="rotation-controls">
                      <div className="rotation-title">Rotation (degrees):</div>
                      <div className={`rotation-inputs ${isMobile ? 'mobile-stack mobile-center' : ''}`}>
                        <div className="rotation-input-group">
                          <span>RX:</span>
                          <input
                            type="number"
                            value={molecule2Rotation.x}
                            onChange={(e) => setMolecule2Rotation({
                              ...molecule2Rotation,
                              x: parseFloat(e.target.value) % 360
                            })}
                            className="rotation-input"
                            step="15"
                          />
                        </div>

                        <div className="rotation-input-group">
                          <span>RY:</span>
                          <input
                            type="number"
                            value={molecule2Rotation.y}
                            onChange={(e) => setMolecule2Rotation({
                              ...molecule2Rotation,
                              y: parseFloat(e.target.value) % 360
                            })}
                            className="rotation-input"
                            step="15"
                          />
                        </div>

                        <div className="rotation-input-group">
                          <span>RZ:</span>
                          <input
                            type="number"
                            value={molecule2Rotation.z}
                            onChange={(e) => setMolecule2Rotation({
                              ...molecule2Rotation,
                              z: parseFloat(e.target.value) % 360
                            })}
                            className="rotation-input"
                            step="15"
                          />
                        </div>

                        <button
                          onClick={() => setMolecule2Rotation({ x: 0, y: 0, z: 0 })}
                          className="rotation-reset-button"
                        >
                          Reset Rotation
                        </button>
                      </div>
                    </div>
                  )}

                  {positioningMode && (
                    <div className="positioning-help">
                      <p>
                        <strong>Controls:</strong> Use keyboard to position Molecule 2: <kbd>←→</kbd> (X-axis), <kbd>↑↓</kbd> (Y-axis), <kbd>PgUp/PgDn</kbd> (Z-axis).
                        Hold <kbd>Shift</kbd> + arrows for rotation.
                      </p>
                    </div>
                  )}
                </div>
              )}

              {/* Visualization */}
              <div className={`visualization-container glass ${isMobile ? 'mobile-smaller-padding' : ''}`}>
                <div className={`visualization-header ${isMobile ? 'mobile-stack' : ''}`}>
                  <div className="visualization-title">
                    <span className="visualization-icon"><Icons.molecule /></span>
                    {activeView === "original" ?
                      (interactionMode ? "Original Molecules" : "Original Structure") :
                      (optimizationType === "classical" ? "Classical" : "Quantum") +
                      (interactionMode ? " Optimized Interaction" : " Optimized Structure")}
                  </div>

                  {optimizationResult && (
                    <div className={`tabs ${isMobile ? 'mobile-full-width' : ''}`}>
                      <div
                        className={`tab ${activeView === "original" ? 'active-blue' : ''} ${isMobile ? 'mobile-smaller-text' : ''}`}
                        onClick={() => setActiveView("original")}
                      >
                        <Icons.molecule /> Original
                      </div>
                      <div
                        className={`tab ${activeView === "optimized" ? (optimizationType === "classical" ? 'active-green' : 'active-blue') : ''} ${isMobile ? 'mobile-smaller-text' : ''}`}
                        onClick={() => setActiveView("optimized")}
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
                  bondParams={classicalParams}
                />

                {interactionMode && molecule1Data && molecule2Data && (
                  <div className="molecule-stats">
                    <span>Molecule 1: {molecule1Data.file1?.atoms?.length || molecule1Data.atoms?.length || 0} atoms</span>
                    <span className="molecule-stats-separator">•</span>
                    <span>Molecule 2: {molecule2Data.file1?.atoms?.length || molecule2Data.atoms?.length || 0} atoms</span>
                    {(molecule2Offset.x !== 0 || molecule2Offset.y !== 0 || molecule2Offset.z !== 0 ||
                      molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0) ? (
                      <>
                        <span className="molecule-stats-separator">•</span>
                        <span>Offset: ({molecule2Offset.x.toFixed(1)}, {molecule2Offset.y.toFixed(1)}, {molecule2Offset.z.toFixed(1)})</span>
                        {(molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0) && (
                          <>
                            <span className="molecule-stats-separator">•</span>
                            <span>Rotation: ({molecule2Rotation.x}°, {molecule2Rotation.y}°, {molecule2Rotation.z}°)</span>
                          </>
                        )}
                      </>
                    ) : null}
                  </div>
                )}
              </div>

              {/* Optimize Button */}
              <div className={`optimize-button-container ${isMobile ? 'mobile-full-width' : ''}`}>
                <button
                  onClick={handleOptimize}
                  disabled={isOptimizeLoading}
                  className={`optimize-button ${optimizationType === "classical" ? 'classical' : 'quantum'} ${isMobile ? 'mobile-full-width' : ''} ${isOptimizeLoading ? 'disabled' : ''}`}
                >
                  {isOptimizeLoading ? (
                    <>
                      <span className="spin optimize-spinner">
                        <Icons.spinner />
                      </span>
                      Optimizing...
                    </>
                  ) : (
                    <>
                      {`Run ${optimizationType === "classical" ? "Classical" : "Quantum"} Optimization${!isSubscribed ? " (Limited)" : ""}`}
                      <div className="optimize-button-shine"></div>
                    </>
                  )}
                </button>

                {!isSubscribed && (
                  <div className={`free-user-notice ${isMobile ? 'mobile-smaller-text mobile-text-center' : ''}`}>
                    Free users are limited to {ITERATION_LIMITS.unsubscribed.classical.toLocaleString()} iterations for classical and {ITERATION_LIMITS.unsubscribed.quantum} for quantum optimizations.
                  </div>
                )}
              </div>

              {taskStatus && (
                <div className={`task-status-container task-status-${taskStatus.status}`}>
                  <div className="task-status-message">{taskStatus.message}</div>
                  <div className="task-status-progress-bar">
                    <div
                      className="task-status-progress"
                      style={{ width: `${taskStatus.progress * 100}%` }}
                    ></div>
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
              document.querySelector('.file-upload').scrollIntoView({ behavior: 'smooth' });
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
          >
            <span className="mobile-nav-icon">
              {isCheckingHealth ? <Icons.spinner /> :
                (serverHealth === true ? <Icons.checkmark /> :
                  (serverHealth === false ? <Icons.warning /> : <Icons.info />))}
            </span>
            Health
            {serverHealthDetails && (
              <div className="health-tooltip">
                {serverHealthDetails.status}: {serverHealthDetails.statusText}
              </div>
            )}
          </button>
        </div>
      )}

      {/* Authentication Modal for Login/Register */}
      {showAuthModal && (
        <div className="popup">
          <div
            className={`popup-content glass ${isMobile ? 'mobile-smaller-padding' : ''}`}
          >
            <button
              onClick={toggleAuthModal}
              className="popup-close"
            >
              <Icons.close />
            </button>

            <div className="popup-scroll">
              {showLoginForm ? (
                <LoginForm
                  toggleForm={() => setShowLoginForm(false)}
                  onAuthSuccess={toggleAuthModal}
                />
              ) : (
                <RegisterForm
                  toggleForm={() => setShowLoginForm(true)}
                  onAuthSuccess={toggleAuthModal}
                />
              )}
            </div>
          </div>
        </div>
      )}

      {/* How To Use Popup */}
      {isHowToUseVisible && (
        <div className="popup">
          <div
            className={`popup-content glass ${isMobile ? 'mobile-smaller-padding' : ''}`}
          >
            <button
              onClick={handleClosePopup}
              className="popup-close"
            >
              <Icons.close />
            </button>

            <div className="popup-scroll">
              <ReactMarkdown>
                {howToUseContent}
              </ReactMarkdown>
            </div>
          </div>
        </div>
      )}

      {/* Theory Popup */}
      {isTheoryVisible && (
        <div className="popup">
          <div
            className={`popup-content glass ${isMobile ? 'mobile-smaller-padding' : ''}`}
          >
            <button
              onClick={handleClosePopup}
              className="popup-close"
            >
              <Icons.close />
            </button>

            <div className="popup-scroll">
              <ReactMarkdown>
                {theoryContent}
              </ReactMarkdown>
            </div>
          </div>
        </div>
      )}

      {/* About Us Popup */}
      {isAboutUsVisible && (
        <div className="popup">
          <div
            className={`popup-content glass ${isMobile ? 'mobile-smaller-padding' : ''}`}
          >
            <button
              onClick={handleClosePopup}
              className="popup-close"
            >
              <Icons.close />
            </button>

            <div className="popup-scroll">
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