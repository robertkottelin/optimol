import React, { useState, useRef, useEffect } from "react";
import axios from "axios";
import * as $3Dmol from '3dmol';
import { loadStripe } from '@stripe/stripe-js';
import { Elements } from '@stripe/react-stripe-js';
import ReactMarkdown from "react-markdown";
import SubscriptionForm from './SubscriptionForm';

const stripePromise = loadStripe('pk_test_kbl0ETzPsoiTwU4ZJMvhsYJw006XVnV4Aq');

// Enhanced design system constants for premium appearance
const COLORS = {
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

const SHADOWS = {
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

const FONTS = {
  heading: "'Manrope', -apple-system, BlinkMacSystemFont, sans-serif",
  body: "'Inter', -apple-system, BlinkMacSystemFont, sans-serif",
  mono: "'Space Mono', monospace",
  weightLight: 300,
  weightRegular: 400,
  weightMedium: 500,
  weightSemiBold: 600,
  weightBold: 700,
};

const SPACING = {
  xs: "4px",
  sm: "8px",
  md: "16px",
  lg: "24px",
  xl: "32px",
  xxl: "48px",
  xxxl: "64px",
};

const BORDER_RADIUS = {
  xs: "4px",
  sm: "6px",
  md: "8px", 
  lg: "12px",
  xl: "16px",
  xxl: "24px",
  pill: "9999px",
  circle: "50%",
};

const TRANSITIONS = {
  fast: "all 0.2s ease",
  medium: "all 0.3s cubic-bezier(0.4, 0, 0.2, 1)",
  slow: "all 0.5s cubic-bezier(0.65, 0, 0.35, 1)",
  bounce: "all 0.5s cubic-bezier(0.34, 1.56, 0.64, 1)",
};

// Constants for iteration limits
const ITERATION_LIMITS = {
  subscribed: {
    classical: 100000,
    quantum: 1000
  },
  unsubscribed: {
    classical: 100,
    quantum: 5
  }
};

// Default optimization parameters
const defaultClassicalParams = {
  temperature: 300,
  max_iterations: 1000,
  bond_threshold: 0.2,
  bond_force_constant: 1000.0,
  angle_force_constant: 500.0
};

const defaultQuantumParams = {
  basis: "6-31g",
  max_iterations: 10,
  convergence_threshold: 0.00001,
  step_size: 0.1
};

// Component style objects with enhanced design
const styles = {
  app: {
    minHeight: "100vh",
    padding: SPACING.xl,
    background: `linear-gradient(135deg, ${COLORS.gradientStart} 0%, ${COLORS.gradientEnd} 100%)`,
    color: COLORS.text,
    fontFamily: FONTS.body,
    position: "relative",
    overflow: "hidden",
  },
  container: {
    maxWidth: "1200px",
    margin: "0 auto",
    position: "relative",
    zIndex: 2,
  },
  header: {
    marginBottom: SPACING.xxl,
    textAlign: "center",
    position: "relative",
  },
  headerTitle: {
    fontFamily: FONTS.heading,
    fontSize: "3rem",
    fontWeight: FONTS.weightBold,
    marginBottom: SPACING.md,
    letterSpacing: "-0.025em",
    background: `linear-gradient(to right, ${COLORS.primary}, ${COLORS.secondary}, ${COLORS.tertiary})`,
    WebkitBackgroundClip: "text",
    WebkitTextFillColor: "transparent",
    display: "inline-block",
    position: "relative",
  },
  headerSubtitle: {
    fontSize: "1.125rem",
    fontWeight: FONTS.weightLight,
    color: COLORS.textSecondary,
    maxWidth: "700px",
    margin: "0 auto",
  },
  card: {
    backgroundColor: "rgba(30, 41, 59, 0.7)",
    backdropFilter: "blur(12px)",
    WebkitBackdropFilter: "blur(12px)",
    borderRadius: BORDER_RADIUS.xl,
    boxShadow: SHADOWS.card,
    padding: SPACING.xl,
    marginBottom: SPACING.xl,
    border: `1px solid ${COLORS.border}`,
    position: "relative",
    overflow: "hidden",
    transition: TRANSITIONS.medium,
  },
  cardWithGlow: {
    backgroundColor: "rgba(30, 41, 59, 0.7)",
    backdropFilter: "blur(12px)",
    WebkitBackdropFilter: "blur(12px)",
    borderRadius: BORDER_RADIUS.xl,
    boxShadow: `${SHADOWS.card}, 0 0 20px rgba(56, 189, 248, 0.3)`,
    padding: SPACING.xl,
    marginBottom: SPACING.xl,
    border: `1px solid rgba(56, 189, 248, 0.3)`,
    position: "relative",
    overflow: "hidden",
    transition: TRANSITIONS.medium,
  },
  section: {
    marginBottom: SPACING.xl,
    animation: "fadeIn 0.5s ease-out",
  },
  sectionTitle: {
    fontSize: "1.5rem",
    fontWeight: FONTS.weightSemiBold,
    marginBottom: SPACING.lg,
    color: COLORS.text,
    fontFamily: FONTS.heading,
    position: "relative",
    display: "inline-block",
    paddingBottom: SPACING.xs,
  },
  sectionTitleUnderline: {
    position: "absolute",
    bottom: 0,
    left: 0,
    width: "60px",
    height: "3px",
    background: `linear-gradient(to right, ${COLORS.primary}, ${COLORS.secondary})`,
    borderRadius: BORDER_RADIUS.pill,
  },
  button: {
    display: "inline-flex",
    alignItems: "center",
    justifyContent: "center",
    padding: `${SPACING.sm} ${SPACING.lg}`,
    borderRadius: BORDER_RADIUS.md,
    fontWeight: FONTS.weightMedium,
    fontSize: "0.875rem",
    cursor: "pointer",
    transition: TRANSITIONS.medium,
    boxShadow: SHADOWS.button,
    border: "none",
    outline: "none",
    textDecoration: "none",
    userSelect: "none",
    position: "relative",
    overflow: "hidden",
  },
  primaryButton: {
    background: COLORS.gradientBlue,
    color: COLORS.text,
    "&:hover": {
      transform: "translateY(-2px)",
      boxShadow: SHADOWS.lg,
    },
    "&:active": {
      transform: "translateY(0)",
      boxShadow: SHADOWS.sm,
    },
  },
  successButton: {
    background: COLORS.gradientEmerald,
    color: COLORS.text,
  },
  infoButton: {
    background: COLORS.gradientIndigo,
    color: COLORS.text,
  },
  tabButton: (isActive, color) => ({
    background: isActive 
      ? (color === COLORS.classical 
          ? COLORS.gradientEmerald 
          : COLORS.gradientBlue) 
      : "rgba(30, 41, 59, 0.7)",
    color: COLORS.text,
    padding: `${SPACING.xs} ${SPACING.md}`,
    boxShadow: isActive ? SHADOWS.md : "none",
    fontWeight: isActive ? FONTS.weightSemiBold : FONTS.weightRegular,
    opacity: isActive ? 1 : 0.8,
    borderRadius: BORDER_RADIUS.md,
    backdropFilter: "blur(8px)",
    WebkitBackdropFilter: "blur(8px)",
  }),
  fileUpload: {
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    justifyContent: "center",
    padding: SPACING.xxl,
    borderRadius: BORDER_RADIUS.xl,
    border: `2px dashed ${COLORS.border}`,
    backgroundColor: "rgba(20, 27, 45, 0.5)",
    cursor: "pointer",
    transition: TRANSITIONS.medium,
    marginBottom: SPACING.xl,
    position: "relative",
    overflow: "hidden",
  },
  fileUploadActive: {
    borderColor: COLORS.primary,
    backgroundColor: `rgba(56, 189, 248, 0.1)`,
  },
  fileUploadIcon: {
    color: COLORS.textSecondary,
    fontSize: "2.5rem",
    marginBottom: SPACING.md,
    transition: TRANSITIONS.medium,
  },
  fileUploadActiveIcon: {
    color: COLORS.primary,
  },
  fileUploadText: {
    fontWeight: FONTS.weightMedium,
    marginBottom: SPACING.xs,
  },
  fileUploadSubtext: {
    fontSize: "0.875rem",
    color: COLORS.textSecondary,
  },
  inputGroup: {
    marginBottom: SPACING.md,
  },
  label: {
    display: "block",
    marginBottom: SPACING.xs,
    fontSize: "0.875rem",
    fontWeight: FONTS.weightMedium,
    color: COLORS.textSecondary,
  },
  input: {
    width: "100%",
    padding: `${SPACING.sm} ${SPACING.md}`,
    backgroundColor: "rgba(15, 23, 42, 0.7)",
    border: `1px solid ${COLORS.border}`,
    borderRadius: BORDER_RADIUS.md,
    color: COLORS.text,
    fontSize: "0.875rem",
    transition: TRANSITIONS.fast,
    "&:focus": {
      borderColor: COLORS.borderHighlight,
      boxShadow: `0 0 0 2px ${COLORS.primary}33`,
      outline: "none",
    },
  },
  parametersContainer: {
    backgroundColor: "rgba(15, 23, 42, 0.7)",
    backdropFilter: "blur(12px)",
    WebkitBackdropFilter: "blur(12px)",
    borderRadius: BORDER_RADIUS.xl,
    padding: SPACING.xl,
    marginBottom: SPACING.xl,
    border: `1px solid ${COLORS.border}`,
    boxShadow: SHADOWS.lg,
  },
  parameterGroup: {
    display: "flex",
    alignItems: "center", 
    marginBottom: SPACING.md,
  },
  parameterLabel: {
    flexBasis: "220px",
    flexShrink: 0,
    fontWeight: FONTS.weightMedium,
    color: COLORS.text,
  },
  parameterInput: {
    width: "120px",
    padding: `${SPACING.xs} ${SPACING.sm}`,
    backgroundColor: "rgba(15, 23, 42, 0.8)",
    border: `1px solid ${COLORS.border}`,
    borderRadius: BORDER_RADIUS.md,
    color: COLORS.text,
    fontSize: "0.875rem",
    fontFamily: FONTS.mono,
  },
  parameterSelect: {
    padding: `${SPACING.xs} ${SPACING.md}`,
    backgroundColor: "rgba(15, 23, 42, 0.8)",
    border: `1px solid ${COLORS.border}`,
    borderRadius: BORDER_RADIUS.md,
    color: COLORS.text,
    fontSize: "0.875rem",
    minWidth: "180px",
  },
  warningText: {
    color: COLORS.warning,
    fontSize: "0.75rem",
    marginLeft: SPACING.sm,
    display: "flex",
    alignItems: "center",
    gap: "6px",
  },
  badge: {
    display: "inline-flex",
    alignItems: "center",
    padding: `${SPACING.xs} ${SPACING.sm}`,
    borderRadius: BORDER_RADIUS.pill,
    fontSize: "0.75rem",
    fontWeight: FONTS.weightSemiBold,
    marginLeft: SPACING.sm,
  },
  warningBadge: {
    backgroundColor: `${COLORS.warning}22`,
    color: COLORS.warning,
    border: `1px solid ${COLORS.warning}44`,
  },
  dangerBadge: {
    backgroundColor: `${COLORS.danger}22`,
    color: COLORS.danger,
    border: `1px solid ${COLORS.danger}44`,
  },
  visualizationContainer: {
    backgroundColor: "rgba(15, 23, 42, 0.7)",
    borderRadius: BORDER_RADIUS.xl,
    overflow: "hidden",
    marginTop: SPACING.xl,
    marginBottom: SPACING.xl,
    boxShadow: SHADOWS.card,
    border: `1px solid ${COLORS.border}`,
    backdropFilter: "blur(12px)",
    WebkitBackdropFilter: "blur(12px)",
    transition: TRANSITIONS.medium,
  },
  visualizationHeader: {
    display: "flex",
    alignItems: "center",
    justifyContent: "space-between",
    padding: SPACING.md,
    borderBottom: `1px solid ${COLORS.border}`,
    backgroundColor: "rgba(15, 23, 42, 0.9)",
  },
  visualizationTitle: {
    fontSize: "1.125rem",
    fontWeight: FONTS.weightSemiBold,
    color: COLORS.text,
    display: "flex",
    alignItems: "center",
    gap: "8px",
  },
  visualizationIcon: {
    color: COLORS.primary,
  },
  visualizationContent: {
    padding: SPACING.md,
  },
  tabs: {
    display: "flex",
    gap: SPACING.xs,
  },
  tab: (isActive, color) => ({
    padding: `${SPACING.sm} ${SPACING.md}`,
    backgroundColor: isActive 
      ? (color === COLORS.classical 
          ? "rgba(16, 185, 129, 0.2)" 
          : "rgba(56, 189, 248, 0.2)")
      : "rgba(15, 23, 42, 0.7)",
    color: isActive 
      ? (color === COLORS.classical 
          ? COLORS.classical 
          : COLORS.quantum)
      : COLORS.textSecondary,
    borderRadius: BORDER_RADIUS.md,
    fontWeight: isActive ? FONTS.weightSemiBold : FONTS.weightRegular,
    cursor: "pointer",
    transition: TRANSITIONS.medium,
    marginRight: SPACING.xs,
    border: `1px solid ${isActive 
      ? (color === COLORS.classical 
          ? "rgba(16, 185, 129, 0.3)" 
          : "rgba(56, 189, 248, 0.3)")
      : "transparent"}`,
    boxShadow: isActive ? SHADOWS.sm : "none",
    display: "flex",
    alignItems: "center",
    gap: "6px",
  }),
  resultsContainer: {
    backgroundColor: "rgba(15, 23, 42, 0.7)",
    borderRadius: BORDER_RADIUS.xl,
    padding: SPACING.xl,
    marginTop: SPACING.xl,
    border: `1px solid ${COLORS.border}`,
    boxShadow: SHADOWS.lg,
    backdropFilter: "blur(12px)",
    WebkitBackdropFilter: "blur(12px)",
    position: "relative",
    overflow: "hidden",
  },
  resultTitle: {
    fontSize: "1.25rem",
    fontWeight: FONTS.weightSemiBold,
    marginBottom: SPACING.lg,
    paddingBottom: SPACING.sm,
    borderBottom: `1px solid ${COLORS.border}`,
    color: COLORS.text,
    fontFamily: FONTS.heading,
    display: "flex",
    alignItems: "center",
    gap: "10px",
  },
  resultIcon: {
    color: COLORS.primary,
  },
  resultItem: {
    display: "flex",
    justifyContent: "space-between",
    alignItems: "center",
    padding: `${SPACING.sm} 0`,
    borderBottom: `1px solid rgba(255, 255, 255, 0.05)`,
  },
  resultLabel: {
    fontWeight: FONTS.weightMedium,
    color: COLORS.textSecondary,
  },
  resultValue: {
    color: COLORS.text,
    fontFamily: FONTS.mono,
    padding: `${SPACING.xs} ${SPACING.sm}`,
    backgroundColor: "rgba(15, 23, 42, 0.5)",
    borderRadius: BORDER_RADIUS.sm,
    fontSize: "0.875rem",
  },
  tooltipContainer: {
    position: "relative",
    display: "inline-block",
  },
  tooltip: {
    position: "absolute",
    bottom: "calc(100% + 10px)",
    left: "50%",
    transform: "translateX(-50%)",
    padding: SPACING.sm,
    backgroundColor: COLORS.cardDark,
    color: COLORS.text,
    borderRadius: BORDER_RADIUS.md,
    fontSize: "0.75rem",
    boxShadow: SHADOWS.md,
    zIndex: 10,
    whiteSpace: "nowrap",
    pointerEvents: "none",
    opacity: 0,
    transition: "opacity 0.2s ease",
    border: `1px solid ${COLORS.border}`,
  },
  tooltipActive: {
    opacity: 1,
  },
  tooltipArrow: {
    position: "absolute",
    top: "100%",
    left: "50%",
    transform: "translateX(-50%)",
    width: 0,
    height: 0,
    borderLeft: "6px solid transparent",
    borderRight: "6px solid transparent",
    borderTop: `6px solid ${COLORS.cardDark}`,
  },
  welcomeMessage: {
    display: "flex",
    alignItems: "center",
    padding: SPACING.lg,
    borderRadius: BORDER_RADIUS.lg,
    border: `1px solid rgba(56, 189, 248, 0.3)`,
    background: "linear-gradient(145deg, rgba(56, 189, 248, 0.1), rgba(37, 99, 235, 0.05))",
    marginBottom: SPACING.xl,
    position: "relative",
    overflow: "hidden",
    boxShadow: `0 8px 16px rgba(0, 0, 0, 0.1), 0 0 0 1px rgba(56, 189, 248, 0.2)`,
  },
  welcomeIcon: {
    color: COLORS.primary,
    fontSize: "1.5rem",
    marginRight: SPACING.md,
    borderRadius: BORDER_RADIUS.circle,
    padding: SPACING.sm,
    backgroundColor: "rgba(56, 189, 248, 0.1)",
    border: `1px solid rgba(56, 189, 248, 0.2)`,
  },
  buttonGroup: {
    display: "flex",
    gap: SPACING.sm,
    marginTop: SPACING.lg,
  },
  buttonGroupItem: {
    flex: 1,
  },
  subscriptionNotice: {
    background: "linear-gradient(145deg, rgba(251, 191, 36, 0.1), rgba(245, 158, 11, 0.05))",
    border: `1px solid rgba(251, 191, 36, 0.3)`,
    borderRadius: BORDER_RADIUS.lg,
    padding: SPACING.lg,
    marginBottom: SPACING.lg,
    display: "flex",
    alignItems: "center",
    boxShadow: `0 4px 6px rgba(0, 0, 0, 0.1)`,
  },
  subscriptionNoticeIcon: {
    color: COLORS.warning,
    marginRight: SPACING.md,
    fontSize: "1.25rem",
    backgroundColor: "rgba(251, 191, 36, 0.1)",
    padding: SPACING.xs,
    borderRadius: BORDER_RADIUS.circle,
  },
  viewerContainer: {
    border: `1px solid ${COLORS.border}`,
    borderRadius: BORDER_RADIUS.lg,
    overflow: "hidden",
    height: "450px",
    backgroundColor: "rgba(15, 23, 42, 0.5)",
    boxShadow: `inset 0 0 20px rgba(0, 0, 0, 0.2)`,
  },
  optimizeButtonContainer: {
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    marginTop: SPACING.xl,
  },
  optimizeButton: (type, isLoading) => ({
    padding: `${SPACING.md} ${SPACING.xxl}`,
    background: isLoading 
      ? "linear-gradient(145deg, #475569, #334155)" 
      : (type === "classical" 
          ? COLORS.gradientEmerald 
          : COLORS.gradientBlue),
    color: COLORS.text,
    borderRadius: BORDER_RADIUS.lg,
    fontWeight: FONTS.weightSemiBold,
    cursor: isLoading ? "not-allowed" : "pointer",
    opacity: isLoading ? 0.7 : 1,
    transition: TRANSITIONS.bounce,
    border: "none",
    boxShadow: isLoading ? SHADOWS.sm : SHADOWS.lg,
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    minWidth: "260px",
    position: "relative",
    overflow: "hidden",
    fontSize: "1rem",
    letterSpacing: "0.5px",
  }),
  optimizeButtonShine: {
    position: "absolute",
    top: 0,
    left: 0,
    width: "100%",
    height: "100%",
    background: "linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.1), transparent)",
    transform: "translateX(-100%)",
    animation: "shimmer 2s infinite",
    pointerEvents: "none",
  },
  freeUserNotice: {
    fontSize: "0.75rem",
    color: COLORS.textSecondary,
    marginTop: SPACING.sm,
    textAlign: "center",
    maxWidth: "400px",
  },
  spinner: {
    width: "18px",
    height: "18px",
    border: `2px solid rgba(255, 255, 255, 0.3)`,
    borderTopColor: COLORS.text,
    borderRadius: "50%",
    animation: "spin 1s linear infinite",
    marginRight: SPACING.sm,
  },
  howToUseButton: {
    position: "absolute",
    top: SPACING.md,
    left: SPACING.md,
    padding: `${SPACING.xs} ${SPACING.sm}`,
    background: "linear-gradient(145deg, rgba(16, 185, 129, 0.9), rgba(5, 150, 105, 0.9))",
    color: "white",
    borderRadius: BORDER_RADIUS.md,
    fontSize: "0.75rem",
    fontWeight: FONTS.weightSemiBold,
    cursor: "pointer",
    zIndex: 5,
    display: "flex",
    alignItems: "center",
    boxShadow: SHADOWS.lg,
    border: "1px solid rgba(16, 185, 129, 0.3)",
    backdropFilter: "blur(8px)",
    WebkitBackdropFilter: "blur(8px)",
  },
  howToUseIcon: {
    marginRight: SPACING.xs,
  },
  cancelSubscriptionButton: {
    position: "absolute",
    top: SPACING.md,
    right: SPACING.md,
    padding: `${SPACING.xs} ${SPACING.sm}`,
    background: "linear-gradient(145deg, rgba(244, 63, 94, 0.9), rgba(225, 29, 72, 0.9))",
    color: "white",
    borderRadius: BORDER_RADIUS.md,
    fontSize: "0.75rem",
    fontWeight: FONTS.weightMedium,
    cursor: "pointer",
    zIndex: 5,
    display: "flex",
    alignItems: "center",
    boxShadow: SHADOWS.lg,
    border: "1px solid rgba(244, 63, 94, 0.3)",
    backdropFilter: "blur(8px)",
    WebkitBackdropFilter: "blur(8px)",
  },
  cancelSubscriptionIcon: {
    marginRight: SPACING.xs,
  },
  popup: {
    position: "fixed",
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    backgroundColor: "rgba(0, 0, 0, 0.8)",
    backdropFilter: "blur(10px)",
    WebkitBackdropFilter: "blur(10px)",
    zIndex: 1000,
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
  },
  popupContent: {
    backgroundColor: "rgba(15, 23, 42, 0.95)",
    borderRadius: BORDER_RADIUS.xxl,
    boxShadow: SHADOWS.popup,
    width: "80%",
    height: "80%",
    maxWidth: "1200px",
    padding: SPACING.xxl,
    position: "relative",
    overflow: "hidden",
    border: `1px solid ${COLORS.border}`,
    backdropFilter: "blur(16px)",
    WebkitBackdropFilter: "blur(16px)",
  },
  popupClose: {
    position: "absolute",
    top: SPACING.md,
    right: SPACING.md,
    width: "36px",
    height: "36px",
    borderRadius: BORDER_RADIUS.circle,
    background: "linear-gradient(145deg, rgba(244, 63, 94, 0.9), rgba(225, 29, 72, 0.9))",
    color: "white",
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    cursor: "pointer",
    zIndex: 1,
    border: "none",
    boxShadow: SHADOWS.md,
  },
  popupScroll: {
    height: "calc(100% - 20px)",
    overflowY: "auto",
    padding: `${SPACING.md} 0`,
  },
  decorativeBg: {
    position: "absolute",
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    background: `
      radial-gradient(circle at 10% 10%, rgba(56, 189, 248, 0.08), transparent 40%),
      radial-gradient(circle at 90% 90%, rgba(16, 185, 129, 0.08), transparent 40%),
      radial-gradient(circle at 90% 10%, rgba(168, 85, 247, 0.08), transparent 30%),
      radial-gradient(circle at 10% 90%, rgba(236, 72, 153, 0.08), transparent 30%)
    `,
    zIndex: 1,
    pointerEvents: "none",
  },
  decorativeLine: {
    position: "absolute",
    height: "1px",
    width: "100%",
    background: "linear-gradient(90deg, transparent, rgba(56, 189, 248, 0.2), transparent)",
    zIndex: 1,
    pointerEvents: "none",
  },
  methodSelectionContainer: {
    display: "flex",
    justifyContent: "center",
    gap: SPACING.lg,
    marginBottom: SPACING.xl,
  },
  methodSelectionButton: (isActive, type) => ({
    flex: 1,
    maxWidth: "280px",
    padding: SPACING.xl,
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    justifyContent: "center",
    backgroundColor: "rgba(15, 23, 42, 0.7)",
    backdropFilter: "blur(12px)",
    WebkitBackdropFilter: "blur(12px)",
    borderRadius: BORDER_RADIUS.xl,
    cursor: "pointer",
    transition: TRANSITIONS.bounce,
    border: `2px solid ${isActive 
      ? (type === "classical" 
          ? COLORS.classical 
          : COLORS.quantum) 
      : "transparent"}`,
    boxShadow: isActive ? SHADOWS.lg : SHADOWS.md,
    position: "relative",
    overflow: "hidden",
  }),
  methodIcon: {
    fontSize: "2.5rem",
    marginBottom: SPACING.md,
    color: COLORS.text,
    padding: SPACING.sm,
    borderRadius: BORDER_RADIUS.circle,
    backgroundColor: isActive => isActive 
      ? (type => type === "classical" 
          ? "rgba(16, 185, 129, 0.1)" 
          : "rgba(56, 189, 248, 0.1)")(type) 
      : "rgba(30, 41, 59, 0.5)",
    border: `1px solid ${isActive => isActive 
      ? (type => type === "classical" 
          ? "rgba(16, 185, 129, 0.3)" 
          : "rgba(56, 189, 248, 0.3)")(type) 
      : "transparent"}`,
  },
  methodTitle: {
    fontWeight: FONTS.weightSemiBold,
    marginBottom: SPACING.xs,
    color: COLORS.text,
    fontSize: "1.125rem",
  },
  methodDescription: {
    fontSize: "0.875rem",
    color: COLORS.textSecondary,
    textAlign: "center",
    lineHeight: 1.5,
  },
  fadeIn: {
    animation: "fadeIn 0.5s ease-in-out",
  },
  contentWrapper: {
    animation: "fadeIn 0.5s ease-in-out",
  },
  glowEffect: {
    position: "absolute",
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    background: `radial-gradient(circle at 50% 50%, ${COLORS.primary}22, transparent 70%)`,
    opacity: 0,
    transition: "opacity 0.5s ease",
    "&:hover": {
      opacity: 1,
    },
  },
};

// Enhanced SVG Icons for UI elements
const Icons = {
  classical: () => (
    <svg width="32" height="32" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" style={{ filter: "drop-shadow(0 0 2px rgba(16, 185, 129, 0.3))" }}>
      <circle cx="12" cy="12" r="8" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <circle cx="12" cy="12" r="3" fill="currentColor"/>
      <path d="M12 4V2" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M12 22V20" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M4 12H2" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M22 12H20" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M6.34315 6.34315L4.92893 4.92893" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M19.0711 19.0711L17.6569 17.6569" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M6.34315 17.6569L4.92893 19.0711" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M19.0711 4.92893L17.6569 6.34315" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
    </svg>
  ),
  quantum: () => (
    <svg width="32" height="32" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" style={{ filter: "drop-shadow(0 0 2px rgba(56, 189, 248, 0.3))" }}>
      <circle cx="12" cy="12" r="9" stroke="currentColor" strokeWidth="1.5" opacity="0.6"/>
      <circle cx="12" cy="12" r="3" stroke="currentColor" strokeWidth="1.5"/>
      <ellipse cx="12" cy="12" rx="6" ry="1.5" transform="rotate(30 12 12)" stroke="currentColor" strokeWidth="1.5"/>
      <ellipse cx="12" cy="12" rx="6" ry="1.5" transform="rotate(-30 12 12)" stroke="currentColor" strokeWidth="1.5"/>
      <ellipse cx="12" cy="12" rx="6" ry="1.5" transform="rotate(90 12 12)" stroke="currentColor" strokeWidth="1.5"/>
    </svg>
  ),
  info: () => (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <circle cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="1.5"/>
      <path d="M12 7V12" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <circle cx="12" cy="16" r="1" fill="currentColor"/>
    </svg>
  ),
  warning: () => (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M12 3L2 21H22L12 3Z" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M12 10V14" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <circle cx="12" cy="18" r="1" fill="currentColor"/>
    </svg>
  ),
  close: () => (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M18 6L6 18" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M6 6L18 18" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  ),
  checkmark: () => (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M5 12L10 17L19 8" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  ),
  reset: () => (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M3 12C3 7.02944 7.02944 3 12 3C16.9706 3 21 7.02944 21 12C21 16.9706 16.9706 21 12 21C7.02944 21 3 16.9706 3 12Z" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M16 12L8 12" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M11 9L8 12L11 15" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  ),
  download: () => (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M12 3V16" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M7 12L12 17L17 12" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M3 21H21" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
    </svg>
  ),
  book: () => (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M4 4V20C4 20 8 18 12 18C16 18 20 20 20 20V4C20 4 16 6 12 6C8 6 4 4 4 4Z" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M12 6V18" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  ),
  cancel: () => (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M15 9L9 15" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M9 9L15 15" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <circle cx="12" cy="12" r="9" stroke="currentColor" strokeWidth="1.5"/>
    </svg>
  ),
  upload: () => (
    <svg width="40" height="40" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" style={{ filter: "drop-shadow(0 0 2px rgba(0, 0, 0, 0.2))" }}>
      <rect x="3" y="3" width="18" height="18" rx="3" stroke="currentColor" strokeWidth="1.5"/>
      <path d="M12 16V8" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M8 12L12 8L16 12" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  ),
  spinner: () => (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <circle cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="1.5" opacity="0.25"/>
      <path d="M12 2C6.47715 2 2 6.47715 2 12" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
    </svg>
  ),
  molecule: () => (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <circle cx="12" cy="8" r="2.5" stroke="currentColor" strokeWidth="1.5"/>
      <circle cx="7" cy="16" r="2.5" stroke="currentColor" strokeWidth="1.5"/>
      <circle cx="17" cy="16" r="2.5" stroke="currentColor" strokeWidth="1.5"/>
      <path d="M12 10.5L7 13.5" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
      <path d="M12 10.5L17 13.5" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round"/>
    </svg>
  ),
  verified: () => (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M9 12L11 14L15 10" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M21 12C21 16.9706 16.9706 21 12 21C7.02944 21 3 16.9706 3 12C3 7.02944 7.02944 3 12 3C16.9706 3 21 7.02944 21 12Z" stroke="currentColor" strokeWidth="1.5"/>
    </svg>
  ),
  settings: () => (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M12 15C13.6569 15 15 13.6569 15 12C15 10.3431 13.6569 9 12 9C10.3431 9 9 10.3431 9 12C9 13.6569 10.3431 15 12 15Z" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
      <path d="M19.4 15C19.2669 15.3016 19.2272 15.6362 19.286 15.9606C19.3448 16.285 19.4995 16.5843 19.73 16.82L19.79 16.88C19.976 17.0657 20.1235 17.2863 20.2241 17.5291C20.3248 17.7719 20.3766 18.0322 20.3766 18.295C20.3766 18.5578 20.3248 18.8181 20.2241 19.0609C20.1235 19.3037 19.976 19.5243 19.79 19.71C19.6043 19.896 19.3837 20.0435 19.1409 20.1441C18.8981 20.2448 18.6378 20.2966 18.375 20.2966C18.1122 20.2966 17.8519 20.2448 17.6091 20.1441C17.3663 20.0435 17.1457 19.896 16.96 19.71L16.9 19.65C16.6643 19.4195 16.365 19.2648 16.0406 19.206C15.7162 19.1472 15.3816 19.1869 15.08 19.32C14.7842 19.4468 14.532 19.6572 14.3543 19.9255C14.1766 20.1938 14.0813 20.5082 14.08 20.83V21C14.08 21.5304 13.8693 22.0391 13.4942 22.4142C13.1191 22.7893 12.6104 23 12.08 23C11.5496 23 11.0409 22.7893 10.6658 22.4142C10.2907 22.0391 10.08 21.5304 10.08 21V20.91C10.0723 20.579 9.96512 20.258 9.77251 19.9887C9.5799 19.7194 9.31074 19.5143 9 19.4C8.69838 19.2669 8.36381 19.2272 8.03941 19.286C7.71502 19.3448 7.41568 19.4995 7.18 19.73L7.12 19.79C6.93425 19.976 6.71368 20.1235 6.47088 20.2241C6.22808 20.3248 5.96783 20.3766 5.705 20.3766C5.44217 20.3766 5.18192 20.3248 4.93912 20.2241C4.69632 20.1235 4.47575 19.976 4.29 19.79C4.10405 19.6043 3.95653 19.3837 3.85588 19.1409C3.75523 18.8981 3.70343 18.6378 3.70343 18.375C3.70343 18.1122 3.75523 17.8519 3.85588 17.6091C3.95653 17.3663 4.10405 17.1457 4.29 16.96L4.35 16.9C4.58054 16.6643 4.73519 16.365 4.794 16.0406C4.85282 15.7162 4.81312 15.3816 4.68 15.08C4.55324 14.7842 4.34276 14.532 4.07447 14.3543C3.80618 14.1766 3.49179 14.0813 3.17 14.08H3C2.46957 14.08 1.96086 13.8693 1.58579 13.4942C1.21071 13.1191 1 12.6104 1 12.08C1 11.5496 1.21071 11.0409 1.58579 10.6658C1.96086 10.2907 2.46957 10.08 3 10.08H3.09C3.42099 10.0723 3.742 9.96512 4.0113 9.77251C4.28059 9.5799 4.48572 9.31074 4.6 9C4.73312 8.69838 4.77282 8.36381 4.714 8.03941C4.65519 7.71502 4.50054 7.41568 4.27 7.18L4.21 7.12C4.02405 6.93425 3.87653 6.71368 3.77588 6.47088C3.67523 6.22808 3.62343 5.96783 3.62343 5.705C3.62343 5.44217 3.67523 5.18192 3.77588 4.93912C3.87653 4.69632 4.02405 4.47575 4.21 4.29C4.39575 4.10405 4.61632 3.95653 4.85912 3.85588C5.10192 3.75523 5.36217 3.70343 5.625 3.70343C5.88783 3.70343 6.14808 3.75523 6.39088 3.85588C6.63368 3.95653 6.85425 4.10405 7.04 4.29L7.1 4.35C7.33568 4.58054 7.63502 4.73519 7.95941 4.794C8.28381 4.85282 8.61838 4.81312 8.92 4.68H9C9.29577 4.55324 9.54802 4.34276 9.72569 4.07447C9.90337 3.80618 9.99872 3.49179 10 3.17V3C10 2.46957 10.2107 1.96086 10.5858 1.58579C10.9609 1.21071 11.4696 1 12 1C12.5304 1 13.0391 1.21071 13.4142 1.58579C13.7893 1.96086 14 2.46957 14 3V3.09C14.0013 3.41179 14.0966 3.72618 14.2743 3.99447C14.452 4.26276 14.7042 4.47324 15 4.6C15.3016 4.73312 15.6362 4.77282 15.9606 4.714C16.285 4.65519 16.5843 4.50054 16.82 4.27L16.88 4.21C17.0657 4.02405 17.2863 3.87653 17.5291 3.77588C17.7719 3.67523 18.0322 3.62343 18.295 3.62343C18.5578 3.62343 18.8181 3.67523 19.0609 3.77588C19.3037 3.87653 19.5243 4.02405 19.71 4.21C19.896 4.39575 20.0435 4.61632 20.1441 4.85912C20.2448 5.10192 20.2966 5.36217 20.2966 5.625C20.2966 5.88783 20.2448 6.14808 20.1441 6.39088C20.0435 6.63368 19.896 6.85425 19.71 7.04L19.65 7.1C19.4195 7.33568 19.2648 7.63502 19.206 7.95941C19.1472 8.28381 19.1869 8.61838 19.32 8.92V9C19.4468 9.29577 19.6572 9.54802 19.9255 9.72569C20.1938 9.90337 20.5082 9.99872 20.83 10H21C21.5304 10 22.0391 10.2107 22.4142 10.5858C22.7893 10.9609 23 11.4696 23 12C23 12.5304 22.7893 13.0391 22.4142 13.4142C22.0391 13.7893 21.5304 14 21 14H20.91C20.5882 14.0013 20.2738 14.0966 20.0055 14.2743C19.7372 14.452 19.5268 14.7042 19.4 15Z" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  ),
  chevronDown: () => (
    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M6 9L12 15L18 9" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  ),
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
      <div style={styles.subscriptionNotice}>
        <span style={styles.subscriptionNoticeIcon}>
          <Icons.warning />
        </span>
        <div>
          <strong>Free Account Limitation:</strong> Iterations capped at {limit.toLocaleString()} (vs. {fullLimit.toLocaleString()} for subscribers).{" "}
          <a 
            href="#" 
            onClick={(e) => { e.preventDefault(); document.querySelector('.subscription-form').scrollIntoView({ behavior: 'smooth' }); }}
            style={{ color: COLORS.warning, textDecoration: "underline" }}
          >
            Subscribe for full capabilities
          </a>
        </div>
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
        backgroundColor: "rgb(15, 23, 42)",
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
            fontSize: 14,
            fontColor: "white",
            backgroundColor: "rgba(0, 0, 0, 0.5)",
            borderRadius: 10,
            padding: 3,
            inFront: true,
          });
        });

        // Style atoms and bonds with enhanced visuals
        viewer.setStyle({}, {
          sphere: { radius: 0.35, scale: 0.9, colorscheme: 'Jmol' },
          stick: { radius: 0.15, colorscheme: 'Jmol' },
        });
        
        // Add better lighting
        viewer.setLightingPreset('default');
        
        // Add depth perception with fog
        viewer.fog(0.1, 100, "rgb(15, 23, 42)");
        
        // Apply ambient occlusion for better depth perception
        viewer.enableAO();
        
        // Add subtle rotation for 3D effect
        viewer.spin('y', 0.5);
        
        viewer.zoomTo();
        viewer.render();
      } catch (error) {
        console.error("Error rendering molecule:", error);
      }
    }, [atoms]);

    return (
      <div className="viewer-container" style={styles.viewerContainer}>
        <div
          ref={viewerRef}
          style={{
            width: "100%",
            height: "100%",
            borderRadius: BORDER_RADIUS.lg,
          }}
        ></div>
      </div>
    );
  };

  // Parameter Configuration Components with enhanced UI
  const ClassicalParametersConfig = () => (
    <div style={styles.parametersContainer}>
      <SubscriptionLimitNotice isSubscribed={isSubscribed} optimizationType="classical" />
      
      <div>
        <div style={styles.parameterGroup}>
          <label style={styles.parameterLabel}>Temperature (K):</label>
          <input 
            type="number" 
            value={classicalParams.temperature}
            min="1"
            max="1000"
            step="10"
            onChange={(e) => handleParamChange('classical', 'temperature', Number(e.target.value))}
            style={styles.parameterInput}
          />
        </div>
        
        <div style={styles.parameterGroup}>
          <label style={styles.parameterLabel}>Max Iterations:</label>
          <input 
            type="number" 
            value={classicalParams.max_iterations}
            min="100"
            max={isSubscribed ? ITERATION_LIMITS.subscribed.classical : ITERATION_LIMITS.unsubscribed.classical}
            step="100"
            onChange={(e) => handleParamChange('classical', 'max_iterations', Number(e.target.value))}
            style={styles.parameterInput}
          />
          {!isSubscribed && classicalParams.max_iterations > ITERATION_LIMITS.unsubscribed.classical && (
            <span style={styles.warningText}>
              <Icons.warning />
              Will be capped at {ITERATION_LIMITS.unsubscribed.classical.toLocaleString()}
            </span>
          )}
          {isSubscribed && classicalParams.max_iterations > ITERATION_LIMITS.subscribed.classical && (
            <span style={styles.warningText}>
              <Icons.warning />
              Will be capped at {ITERATION_LIMITS.subscribed.classical.toLocaleString()}
            </span>
          )}
        </div>
        
        {showAdvancedParams && (
          <>
            <div style={styles.parameterGroup}>
              <label style={styles.parameterLabel}>Bond Threshold (nm):</label>
              <input 
                type="number" 
                value={classicalParams.bond_threshold}
                min="0.1"
                max="0.5"
                step="0.01"
                onChange={(e) => handleParamChange('classical', 'bond_threshold', Number(e.target.value))}
                style={styles.parameterInput}
              />
            </div>
            
            <div style={styles.parameterGroup}>
              <label style={styles.parameterLabel}>Bond Force Constant (kJ/mol/nm²):</label>
              <input 
                type="number" 
                value={classicalParams.bond_force_constant}
                min="100"
                max="10000"
                step="100"
                onChange={(e) => handleParamChange('classical', 'bond_force_constant', Number(e.target.value))}
                style={styles.parameterInput}
              />
            </div>
            
            <div style={styles.parameterGroup}>
              <label style={styles.parameterLabel}>Angle Force Constant (kJ/mol/rad²):</label>
              <input 
                type="number" 
                value={classicalParams.angle_force_constant}
                min="50"
                max="5000"
                step="50"
                onChange={(e) => handleParamChange('classical', 'angle_force_constant', Number(e.target.value))}
                style={styles.parameterInput}
              />
            </div>
          </>
        )}
      </div>
      
      <div style={styles.buttonGroup}>
        <button
          onClick={() => handleResetParams('classical')}
          style={{
            ...styles.button,
            background: "rgba(15, 23, 42, 0.7)",
            color: COLORS.text,
            display: "flex",
            alignItems: "center",
            gap: "5px",
            padding: "10px 16px",
            borderRadius: BORDER_RADIUS.md,
            border: `1px solid ${COLORS.border}`,
          }}
        >
          <span><Icons.reset /></span>
          Reset to Defaults
        </button>
        
        <button
          onClick={() => setShowAdvancedParams(!showAdvancedParams)}
          style={{
            ...styles.button,
            background: "linear-gradient(135deg, rgba(56, 189, 248, 0.1), rgba(59, 130, 246, 0.1))",
            color: COLORS.info,
            display: "flex",
            alignItems: "center",
            gap: "5px",
            padding: "10px 16px",
            borderRadius: BORDER_RADIUS.md,
            border: `1px solid rgba(56, 189, 248, 0.3)`,
          }}
        >
          {showAdvancedParams ? (
            <>Hide Advanced Parameters <Icons.chevronDown /></>
          ) : (
            <>Advanced Parameters <Icons.settings /></>
          )}
        </button>
      </div>
    </div>
  );

  const QuantumParametersConfig = () => (
    <div style={styles.parametersContainer}>
      <SubscriptionLimitNotice isSubscribed={isSubscribed} optimizationType="quantum" />
      
      <div>
        <div style={styles.parameterGroup}>
          <label style={styles.parameterLabel}>Basis Set:</label>
          <select 
            value={quantumParams.basis}
            onChange={(e) => handleParamChange('quantum', 'basis', e.target.value)}
            style={styles.parameterSelect}
          >
            <option value="sto-3g">STO-3G (Minimal)</option>
            <option value="6-31g">6-31G (Standard)</option>
            {isSubscribed && <option value="6-311g">6-311G (Extended)</option>}
            {isSubscribed && <option value="cc-pvdz">cc-pVDZ (Double Zeta)</option>}
          </select>
          {!isSubscribed && (quantumParams.basis === "6-311g" || quantumParams.basis === "cc-pvdz") && (
            <span style={styles.dangerBadge}>
              Extended basis sets require subscription
            </span>
          )}
        </div>
        
        <div style={styles.parameterGroup}>
          <label style={styles.parameterLabel}>Max Iterations:</label>
          <input 
            type="number" 
            value={quantumParams.max_iterations}
            min="1"
            max={isSubscribed ? ITERATION_LIMITS.subscribed.quantum : ITERATION_LIMITS.unsubscribed.quantum}
            onChange={(e) => handleParamChange('quantum', 'max_iterations', Number(e.target.value))}
            style={styles.parameterInput}
          />
          {!isSubscribed && quantumParams.max_iterations > ITERATION_LIMITS.unsubscribed.quantum && (
            <span style={styles.warningText}>
              <Icons.warning />
              Will be capped at {ITERATION_LIMITS.unsubscribed.quantum}
            </span>
          )}
          {isSubscribed && quantumParams.max_iterations > ITERATION_LIMITS.subscribed.quantum && (
            <span style={styles.warningText}>
              <Icons.warning />
              Will be capped at {ITERATION_LIMITS.subscribed.quantum}
            </span>
          )}
        </div>
        
        {showAdvancedParams && (
          <>
            <div style={styles.parameterGroup}>
              <label style={styles.parameterLabel}>Convergence Threshold:</label>
              <input 
                type="number" 
                value={quantumParams.convergence_threshold}
                min="0.000001"
                max="0.01"
                step="0.000001"
                onChange={(e) => handleParamChange('quantum', 'convergence_threshold', Number(e.target.value))}
                style={styles.parameterInput}
              />
            </div>
            
            <div style={styles.parameterGroup}>
              <label style={styles.parameterLabel}>Step Size:</label>
              <input 
                type="number" 
                value={quantumParams.step_size}
                min="0.01"
                max="1.0"
                step="0.01"
                onChange={(e) => handleParamChange('quantum', 'step_size', Number(e.target.value))}
                style={styles.parameterInput}
              />
            </div>
          </>
        )}
      </div>
      
      <div style={styles.buttonGroup}>
        <button
          onClick={() => handleResetParams('quantum')}
          style={{
            ...styles.button,
            background: "rgba(15, 23, 42, 0.7)",
            color: COLORS.text,
            display: "flex",
            alignItems: "center",
            gap: "5px",
            padding: "10px 16px",
            borderRadius: BORDER_RADIUS.md,
            border: `1px solid ${COLORS.border}`,
          }}
        >
          <span><Icons.reset /></span>
          Reset to Defaults
        </button>
        
        <button
          onClick={() => setShowAdvancedParams(!showAdvancedParams)}
          style={{
            ...styles.button,
            background: "linear-gradient(135deg, rgba(56, 189, 248, 0.1), rgba(59, 130, 246, 0.1))",
            color: COLORS.info,
            display: "flex",
            alignItems: "center",
            gap: "5px",
            padding: "10px 16px",
            borderRadius: BORDER_RADIUS.md,
            border: `1px solid rgba(56, 189, 248, 0.3)`,
          }}
        >
          {showAdvancedParams ? (
            <>Hide Advanced Parameters <Icons.chevronDown /></>
          ) : (
            <>Advanced Parameters <Icons.settings /></>
          )}
        </button>
      </div>
    </div>
  );

  // Replace the existing OptimizationResults component with this updated version
  const OptimizationResults = () => {
    if (!optimizationResult || !optimizationResult.result) return null;
    
    const result = optimizationResult.result;
    
    if (result.error) {
      return (
        <div style={styles.resultsContainer}>
          <h3 style={styles.resultTitle}>
            <span style={{ ...styles.resultIcon, color: COLORS.danger }}>
              <Icons.warning />
            </span>
            Optimization Error
          </h3>
          <div style={{ ...styles.resultItem, color: COLORS.danger }}>
            {result.error}
          </div>
        </div>
      );
    }
    
    // Check if the optimization type matches the result data structure
    const resultMatchesType = (
      (optimizationType === "classical" && result.metadata.method === "classical_molecular_dynamics") ||
      (optimizationType === "quantum" && result.metadata.method === "quantum_chemistry")
    );
    
    // If types don't match, show type mismatch message
    if (!resultMatchesType) {
      return (
        <div style={styles.resultsContainer}>
          <h3 style={styles.resultTitle}>
            <span style={{ ...styles.resultIcon, color: COLORS.warning }}>
              <Icons.warning />
            </span>
            Results Not Available
          </h3>
          <div style={{ ...styles.resultItem }}>
            <p>Please run a {optimizationType} optimization to see {optimizationType} results.</p>
            <p>Currently viewing results from a previous {optimizationType === "classical" ? "quantum" : "classical"} optimization.</p>
          </div>
        </div>
      );
    }
    
    return (
      <div style={styles.resultsContainer}>
        <h3 style={styles.resultTitle}>
          <span style={{ ...styles.resultIcon, color: optimizationType === "classical" ? COLORS.classical : COLORS.quantum }}>
            {optimizationType === "classical" ? <Icons.classical /> : <Icons.quantum />}
          </span>
          {optimizationType === "classical" ? "Classical" : "Quantum"} Optimization Results
        </h3>
        
        <div style={styles.resultItem}>
          <span style={styles.resultLabel}>Method:</span>
          <span style={styles.resultValue}>{result.metadata.method}</span>
        </div>
        
        <div style={styles.resultItem}>
          <span style={styles.resultLabel}>Library:</span>
          <span style={styles.resultValue}>{result.metadata.library}</span>
        </div>
        
        {optimizationType === "classical" ? (
          <>
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Temperature:</span>
              <span style={styles.resultValue}>
                {result.metadata.parameters?.temperature !== undefined ? 
                  `${result.metadata.parameters.temperature} K` : 
                  "N/A"}
              </span>
            </div>
            
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Final Energy:</span>
              <span style={styles.resultValue}>
                {result.metadata.final_energy_kj_mol !== undefined ? 
                  `${result.metadata.final_energy_kj_mol.toFixed(4)} kJ/mol` : 
                  "N/A"}
              </span>
            </div>
            
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Bonds Detected:</span>
              <span style={styles.resultValue}>
                {result.metadata.bonds_detected !== undefined ? 
                  result.metadata.bonds_detected : 
                  "N/A"}
              </span>
            </div>
            
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Angles Detected:</span>
              <span style={styles.resultValue}>
                {result.metadata.angles_detected !== undefined ? 
                  result.metadata.angles_detected : 
                  "N/A"}
              </span>
            </div>
          </>
        ) : (
          <>
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Basis Set:</span>
              <span style={styles.resultValue}>
                {result.metadata.parameters?.basis !== undefined ?
                  result.metadata.parameters.basis :
                  "N/A"}
              </span>
            </div>
            
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Theory Level:</span>
              <span style={styles.resultValue}>
                {result.metadata.theory_level !== undefined ?
                  result.metadata.theory_level :
                  "N/A"}
              </span>
            </div>
            
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Final Energy:</span>
              <span style={styles.resultValue}>
                {result.metadata.final_energy_hartree !== undefined ? 
                  `${result.metadata.final_energy_hartree.toFixed(6)} Hartree` : 
                  "N/A"}
              </span>
            </div>
            
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Iterations:</span>
              <span style={styles.resultValue}>
                {result.metadata.iterations !== undefined ?
                  result.metadata.iterations :
                  "N/A"}
              </span>
            </div>
            
            <div style={styles.resultItem}>
              <span style={styles.resultLabel}>Converged:</span>
              <span style={styles.resultValue}>
                {result.metadata.converged !== undefined ?
                  (result.metadata.converged ? "Yes" : "No") :
                  "N/A"}
              </span>
            </div>
          </>
        )}
        
        <div style={styles.resultItem}>
          <span style={styles.resultLabel}>Duration:</span>
          <span style={styles.resultValue}>
            {result.metadata.duration_seconds !== undefined ?
              `${result.metadata.duration_seconds.toFixed(2)} seconds` :
              "N/A"}
          </span>
        </div>
        
        <div style={{ marginTop: SPACING.lg }}>
          <button
            onClick={handleDownload}
            style={{
              ...styles.button,
              background: optimizationType === "classical" ? COLORS.gradientEmerald : COLORS.gradientBlue,
              color: COLORS.text,
              display: "flex",
              alignItems: "center",
              gap: "8px",
              padding: "10px 20px",
              borderRadius: BORDER_RADIUS.md,
              fontWeight: FONTS.weightSemiBold,
            }}
          >
            <span><Icons.download /></span>
            Download Result
          </button>
        </div>
      </div>
    );
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
          onClick={fetchHowToUse}
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
                    padding: SPACING.sm,
                    borderRadius: BORDER_RADIUS.circle,
                    marginBottom: SPACING.md
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
                    padding: SPACING.sm,
                    borderRadius: BORDER_RADIUS.circle,
                    marginBottom: SPACING.md
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
              {optimizationType === "classical" ? 
                <ClassicalParametersConfig /> : 
                <QuantumParametersConfig />
              }
              
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
                        style={styles.tab(activeView === "original", COLORS.primary)}
                        onClick={() => setActiveView("original")}
                      >
                        <Icons.molecule /> Original
                      </div>
                      <div 
                        style={styles.tab(
                          activeView === "optimized", 
                          optimizationType === "classical" ? COLORS.classical : COLORS.quantum
                        )}
                        onClick={() => setActiveView("optimized")}
                      >
                        {optimizationType === "classical" ? <Icons.classical /> : <Icons.quantum />}
                        {optimizationType === "classical" ? "Classical" : "Quantum"} Optimized
                      </div>
                    </div>
                  )}
                </div>
                
                <MoleculeViewer />
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
              {optimizationResult && <OptimizationResults />}
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