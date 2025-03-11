import React from 'react';
import { styles } from '../styles/components';
import { ITERATION_LIMITS } from '../styles/constants';
import { Icons } from './Icons';

// Subscription limit notice component
export const SubscriptionLimitNotice = ({ isSubscribed, optimizationType }) => {
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
          style={{ color: "#fbbf24", textDecoration: "underline" }}
        >
          Subscribe for full capabilities
        </a>
      </div>
    </div>
  );
};

// Classical parameters configuration component
export const ClassicalParametersConfig = ({ 
  isSubscribed, 
  classicalParams, 
  showAdvancedParams, 
  handleParamChange, 
  handleResetParams, 
  setShowAdvancedParams 
}) => (
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
      </div>
      
      <div style={styles.parameterGroup}>
        <label style={styles.parameterLabel}>Force All Iterations:</label>
        <input 
          type="checkbox" 
          checked={classicalParams.force_iterations}
          onChange={(e) => handleParamChange('classical', 'force_iterations', e.target.checked)}
          style={styles.parameterCheckbox}
        />
        <span style={styles.parameterHelp}>
          Forces execution of exactly the specified number of iterations regardless of convergence
        </span>
      </div>
      
      {showAdvancedParams && (
        <>
          <div style={styles.parameterGroup}>
            <label style={styles.parameterLabel}>Convergence Tolerance (kJ/mol):</label>
            <input 
              type="number" 
              value={classicalParams.tolerance}
              min="0.001"
              max="10.0"
              step="0.1"
              onChange={(e) => handleParamChange('classical', 'tolerance', Number(e.target.value))}
              style={styles.parameterInput}
            />
          </div>
          
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
          color: "#f0f4f8",
          display: "flex",
          alignItems: "center",
          gap: "5px",
          padding: "10px 16px",
          borderRadius: "8px",
          border: `1px solid rgba(255, 255, 255, 0.12)`,
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
          color: "#3b82f6",
          display: "flex",
          alignItems: "center",
          gap: "5px",
          padding: "10px 16px",
          borderRadius: "8px",
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

// Quantum parameters configuration component
export const QuantumParametersConfig = ({ 
  isSubscribed, 
  quantumParams, 
  showAdvancedParams, 
  handleParamChange, 
  handleResetParams, 
  setShowAdvancedParams 
}) => (
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
          color: "#f0f4f8",
          display: "flex",
          alignItems: "center",
          gap: "5px",
          padding: "10px 16px",
          borderRadius: "8px",
          border: `1px solid rgba(255, 255, 255, 0.12)`,
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
          color: "#3b82f6",
          display: "flex",
          alignItems: "center",
          gap: "5px",
          padding: "10px 16px",
          borderRadius: "8px",
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