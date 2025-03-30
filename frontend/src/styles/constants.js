// Design system constants - extracted from App.js

// Color palette
export const COLORS = {
  // Primary palette
  primary: "#38bdf8",      // Sky blue
  secondary: "#818cf8",    // Indigo
  tertiary: "#c084fc",     // Purple

  // State colors
  success: "#10b981",      // Emerald green
  danger: "#f43f5e",       // Rose red
  warning: "#fbbf24",      // Amber
  info: "#3b82f6",         // Blue

  // Background colors - Darker, more premium
  background: "#0c1021",   // Deep blue-black
  surface: "#141b2d",      // Dark blue
  card: "#1e293b",         // Slate
  cardDark: "#0f172a",     // Darker slate

  // Text colors
  text: "#f0f4f8",         // Nearly white for text
  textSecondary: "#94a3b8", // Light slate for secondary text
  textTertiary: "#64748b",  // Even lighter for tertiary text

  // Contextual colors
  border: "rgba(255, 255, 255, 0.12)",
  borderHighlight: "rgba(56, 189, 248, 0.5)",
  overlay: "rgba(17, 24, 39, 0.8)",

  // Method-specific colors
  classical: "#10b981",    // Emerald green
  quantum: "#38bdf8",      // Sky blue

  // Gradient colors
  gradientStart: "#0c1021",
  gradientEnd: "#141b2d",
  gradientBlue: "linear-gradient(135deg, #38bdf8, #2563eb)",
  gradientIndigo: "linear-gradient(135deg, #818cf8, #4f46e5)",
  gradientEmerald: "linear-gradient(135deg, #10b981, #059669)",
  gradientRose: "linear-gradient(135deg, #f43f5e, #e11d48)",
};

// Shadow system
export const SHADOWS = {
  sm: "0 1px 2px rgba(0, 0, 0, 0.2)",
  md: "0 4px 6px -1px rgba(0, 0, 0, 0.2), 0 2px 4px -1px rgba(0, 0, 0, 0.1)",
  lg: "0 10px 25px -5px rgba(0, 0, 0, 0.3), 0 8px 10px -6px rgba(0, 0, 0, 0.2)",
  xl: "0 20px 25px -5px rgba(0, 0, 0, 0.3), 0 10px 10px -5px rgba(0, 0, 0, 0.2)",
  inner: "inset 0 2px 4px rgba(0, 0, 0, 0.2)",
  button: "0 4px 10px rgba(0, 0, 0, 0.3)",
  card: "0 20px 25px -5px rgba(0, 0, 0, 0.3), 0 8px 10px -6px rgba(0, 0, 0, 0.2)",
  highlight: "0 0 15px rgba(56, 189, 248, 0.5)",
  popup: "0 25px 50px -12px rgba(0, 0, 0, 0.5)",
};

// Typography system
export const FONTS = {
  heading: "'Manrope', -apple-system, BlinkMacSystemFont, sans-serif",
  body: "'Inter', -apple-system, BlinkMacSystemFont, sans-serif",
  mono: "'Space Mono', monospace",
  weightLight: 300,
  weightRegular: 400,
  weightMedium: 500,
  weightSemiBold: 600,
  weightBold: 700,
};

// Spacing system
export const SPACING = {
  xs: "4px",
  sm: "8px",
  md: "16px",
  lg: "24px",
  xl: "32px",
  xxl: "48px",
  xxxl: "64px",
};

// Border radius system
export const BORDER_RADIUS = {
  xs: "4px",
  sm: "6px",
  md: "8px",
  lg: "12px",
  xl: "16px",
  xxl: "24px",
  pill: "9999px",
  circle: "50%",
};

// Transition system
export const TRANSITIONS = {
  fast: "all 0.2s ease",
  medium: "all 0.3s cubic-bezier(0.4, 0, 0.2, 1)",
  slow: "all 0.5s cubic-bezier(0.65, 0, 0.35, 1)",
  bounce: "all 0.5s cubic-bezier(0.34, 1.56, 0.64, 1)",
};

// User type-based iteration limits
export const ITERATION_LIMITS = {
  subscribed: {
    classical: 1000000,
    quantum: 1000
  },
  unsubscribed: {
    classical: 1000,
    quantum: 5
  }
};

// Molecule size limits for quantum calculations
export const QUANTUM_SIZE_LIMITS = {
  singleMolecule: 20,  // Maximum atoms for single molecule quantum optimization
  interactionTotal: 30,  // Maximum total atoms for interaction quantum optimization
  interactionPerMolecule: 15  // Maximum atoms per molecule for interaction quantum optimization
};

// Default optimization parameters
export const defaultClassicalParams = {
  temperature: 300,
  max_iterations: 1000,
  bond_threshold: 0.2,
  bond_force_constant: 1000.0,
  angle_force_constant: 500.0,
  tolerance: 10.0,
  force_iterations: false,
  covalent_display_threshold: 1.1,  // in Angstroms
  hydrogen_display_threshold: 3.2,  // in Angstroms 
  angle_display_threshold: 45 // in degrees
};

export const defaultQuantumParams = {
  basis: "6-31g",
  max_iterations: 10,
  convergence_threshold: 0.00001,
  step_size: 0.1
};

// Test molecule structures
export const TEST_MOLECULES = {
  // Water molecule structure
  water: {
    file1: {
      atoms: [
        { id: 1, element: "O", x: 0.0000, y: 0.0000, z: 0.0000 },
        { id: 2, element: "H", x: 0.7572, y: 0.5860, z: 0.0000 },
        { id: 3, element: "H", x: -0.7572, y: 0.5860, z: 0.0000 }
      ],
      metadata: {
        name: "Water",
        formula: "H2O",
        description: "Water molecule in standard configuration with O-H bond length of 0.9584 Å and H-O-H angle of 104.5 degrees"
      }
    }
  },
  
  // Acetic Acid molecule structure
  aceticAcid: {
    file1: {
      atoms: [
        { id: 1, element: "C", x: 0.0000, y: 0.0000, z: 0.0000 },
        { id: 2, element: "C", x: 1.5000, y: 0.0000, z: 0.0000 },
        { id: 3, element: "O", x: 2.1500, y: 1.0500, z: 0.0000 },
        { id: 4, element: "O", x: 2.1500, y: -1.1500, z: 0.0000 },
        { id: 5, element: "H", x: -0.4000, y: 1.0000, z: 0.0000 },
        { id: 6, element: "H", x: -0.4000, y: -0.5000, z: 0.8800 },
        { id: 7, element: "H", x: -0.4000, y: -0.5000, z: -0.8800 },
        { id: 8, element: "H", x: 3.1200, y: -1.1000, z: 0.0000 }
      ],
      metadata: {
        name: "Acetic Acid",
        formula: "CH₃COOH",
        description: "Acetic acid with carboxyl group positioned for hydrogen bonding"
      }
    }
  },
  
  // Methanol molecule structure
  methanol: {
    file1: {
      atoms: [
        { id: 1, element: "C", x: 1.0000, y: 0.0000, z: 0.0000 },
        { id: 2, element: "O", x: 2.4000, y: 0.0000, z: 0.0000 },
        { id: 3, element: "H", x: 0.4000, y: 1.0000, z: 0.0000 },
        { id: 4, element: "H", x: 0.4000, y: -0.5000, z: 0.8800 },
        { id: 5, element: "H", x: 0.4000, y: -0.5000, z: -0.8800 },
        { id: 6, element: "H", x: 2.7500, y: -0.9000, z: 0.0000 }
      ],
      metadata: {
        name: "Methanol",
        formula: "CH₃OH",
        description: "Methanol with hydroxyl group positioned for hydrogen bonding"
      }
    }
  },
  
  // Ibuprofen molecule structure - positioned for binding
  ibuprofen: {
    file1: {
      atoms: [
        // Coordinates precisely positioned for correct binding pose with COX-2
        { id: 1, element: "C", x: 3.8322, y: -0.5536, z: 1.8151 },
        { id: 2, element: "C", x: 2.4540, y: -0.3860, z: 1.8391 },
        { id: 3, element: "C", x: 1.5830, y: -1.4776, z: 1.8391 },
        { id: 4, element: "C", x: 0.2048, y: -1.3099, z: 1.8061 },
        { id: 5, element: "C", x: -0.3027, y: 0.0494, z: 1.7892 },
        { id: 6, element: "C", x: 0.5682, y: 1.1410, z: 1.7891 },
        { id: 7, element: "C", x: 1.9464, y: 0.9733, z: 1.8382 },
        { id: 8, element: "C", x: -1.7899, y: 0.2272, z: 1.7892 },
        { id: 9, element: "C", x: -2.2843, y: 1.6793, z: 1.7550 },
        { id: 10, element: "C", x: -2.6511, y: -0.8623, z: 1.8233 },
        { id: 11, element: "C", x: -4.1383, y: -0.6845, z: 1.8233 },
        { id: 12, element: "C", x: -4.7278, y: -0.5455, z: 3.2216 },
        { id: 13, element: "C", x: -4.7278, y: -0.5455, z: 0.4251 },
        { id: 14, element: "O", x: -2.2326, y: -2.1200, z: 1.8606 },   // Carboxylate oxygen 1
        { id: 15, element: "O", x: -1.5949, y: 2.6467, z: 1.7183 },    // Carboxylate oxygen 2
        { id: 16, element: "H", x: 4.3166, y: -0.0546, z: 2.6707 },
        { id: 17, element: "H", x: 4.2856, y: -0.1326, z: 0.9084 },
        { id: 18, element: "H", x: 4.0959, y: -1.6242, z: 1.8457 },
        { id: 19, element: "H", x: 1.9679, y: -2.4897, z: 1.8681 },
        { id: 20, element: "H", x: -0.4595, y: -2.1729, z: 1.8062 },
        { id: 21, element: "H", x: 0.1833, y: 2.1533, z: 1.7602 },
        { id: 22, element: "H", x: 2.6108, y: 1.8365, z: 1.8382 },
        { id: 23, element: "H", x: -3.3751, y: 1.7683, z: 1.7549 },
        { id: 24, element: "H", x: -5.8165, y: -0.4203, z: 3.2216 },
        { id: 25, element: "H", x: -4.2887, y: 0.3239, z: 3.7210 },
        { id: 26, element: "H", x: -4.4834, y: -1.4454, z: 3.8059 },
        { id: 27, element: "H", x: -5.8165, y: -0.4203, z: 0.4251 },
        { id: 28, element: "H", x: -4.4834, y: -1.4454, z: -0.1593 },
        { id: 29, element: "H", x: -4.2887, y: 0.3239, z: -0.0743 },
        { id: 30, element: "H", x: -1.2744, y: -2.2019, z: 1.8606 },   // Carboxylic hydrogen
        { id: 31, element: "H", x: -2.0345, y: 3.5220, z: 1.6963 }
      ],
      metadata: {
        name: "Ibuprofen",
        formula: "C13H18O2",
        description: "Ibuprofen molecule positioned for binding to COX-2. The carboxylate group is oriented for hydrogen bond interaction with Arg120, and the hydrophobic tail extends into the binding pocket."
      }
    }
  },
  
  // COX-2 binding site - precisely positioned for ibuprofen binding
  cox2BindingSite: {
    file1: {
      atoms: [
        // Arg120 - positioned for hydrogen bonding with ibuprofen carboxylate
        { id: 1, element: "N", x: -1.2000, y: 3.0000, z: 1.8000 },     // Terminal nitrogen
        { id: 2, element: "C", x: -0.7000, y: 1.9000, z: 2.2000 },     // Carbon
        { id: 3, element: "C", x: -1.5000, y: 1.0000, z: 3.0000 },     // Carbon
        { id: 4, element: "N", x: -1.1000, y: -0.1500, z: 3.3000 },    // Nitrogen with H for H-bond
        { id: 5, element: "N", x: -2.7000, y: 1.4000, z: 3.3000 },     // Nitrogen with H for H-bond
        { id: 6, element: "H", x: -0.8000, y: 3.8000, z: 1.4000 },     // Hydrogen on terminal N
        { id: 7, element: "H", x: -0.2000, y: -0.2500, z: 3.2000 },    // Hydrogen on N4 for H-bond
        { id: 8, element: "H", x: -1.8000, y: -0.8000, z: 3.5000 },    // Hydrogen on N4 for H-bond
        { id: 9, element: "H", x: -3.0000, y: 2.3000, z: 3.1000 },     // Hydrogen on N5 for H-bond
        { id: 10, element: "H", x: -3.2000, y: 0.8000, z: 3.8000 },    // Hydrogen on N5 for H-bond
        
        // Tyr355 - positioned for hydrogen bond with carboxylate
        { id: 11, element: "C", x: -3.0000, y: 0.3000, z: 0.4000 },
        { id: 12, element: "C", x: -3.8000, y: -0.8000, z: 0.0000 },
        { id: 13, element: "C", x: -3.4000, y: -1.5000, z: -1.2000 },
        { id: 14, element: "O", x: -4.2000, y: -2.5000, z: -1.6000 },
        { id: 15, element: "H", x: -3.9000, y: -3.0000, z: -2.3000 },  // Hydroxyl H for H-bond
        
        // Ser530 - positioned near the binding site
        { id: 16, element: "C", x: 0.3000, y: -0.5000, z: 0.5000 },
        { id: 17, element: "C", x: 1.5000, y: -1.0000, z: -0.2000 },
        { id: 18, element: "O", x: 2.7000, y: -0.5000, z: 0.3000 },
        { id: 19, element: "H", x: 3.4000, y: -0.8000, z: -0.2000 },   // Hydroxyl H for H-bond
        
        // Val349 - hydrophobic pocket positioned for isobutyl group
        { id: 20, element: "C", x: -4.5000, y: -3.0000, z: 1.5000 },
        { id: 21, element: "C", x: -5.8000, y: -3.5000, z: 1.0000 },
        { id: 22, element: "C", x: -4.3000, y: -3.2000, z: 3.0000 },
        { id: 23, element: "H", x: -6.2000, y: -4.2000, z: 1.7000 },
        { id: 24, element: "H", x: -4.9000, y: -4.0000, z: 3.4000 },
        
        // Leu352 - positioned for hydrophobic interactions
        { id: 25, element: "C", x: -2.0000, y: -3.0000, z: 0.8000 },
        { id: 26, element: "C", x: -2.2000, y: -4.0000, z: -0.3000 },
        { id: 27, element: "C", x: -1.0000, y: -4.8000, z: -0.7000 },
        { id: 28, element: "C", x: -1.8000, y: -3.6000, z: 2.1000 },
        { id: 29, element: "H", x: -1.2000, y: -5.6000, z: -1.4000 },
        
        // Phe518 - positioned for π-stacking with ibuprofen's aromatic ring
        { id: 30, element: "C", x: 1.8000, y: 1.5000, z: 2.0000 },
        { id: 31, element: "C", x: 2.7000, y: 2.2000, z: 1.2000 },
        { id: 32, element: "C", x: 3.4000, y: 3.3000, z: 1.7000 },
        { id: 33, element: "C", x: 3.2000, y: 3.7000, z: 3.0000 },
        { id: 34, element: "C", x: 2.3000, y: 3.0000, z: 3.8000 },
        { id: 35, element: "C", x: 1.6000, y: 1.9000, z: 3.3000 },
        { id: 36, element: "H", x: 4.1000, y: 3.8000, z: 1.1000 }
      ],
      metadata: {
        name: "COX-2 Binding Site",
        formula: "Protein fragment",
        description: "Cyclooxygenase-2 (COX-2) binding site specifically positioned for hydrogen bonding with ibuprofen. Includes explicit hydrogens on Arg120 guanidinium group (atoms 7-10) for hydrogen bonding detection. The binding site is positioned within the optimal H-bond distance range (1.5-3.2Å) from ibuprofen's carboxylate group."
      }
    }
  }
};
