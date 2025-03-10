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
  
  // Default optimization parameters
  export const defaultClassicalParams = {
    temperature: 300,
    max_iterations: 1000,
    bond_threshold: 0.2,
    bond_force_constant: 1000.0,
    angle_force_constant: 500.0
  };
  
  export const defaultQuantumParams = {
    basis: "6-31g",
    max_iterations: 10,
    convergence_threshold: 0.00001,
    step_size: 0.1
  };