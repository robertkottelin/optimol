import React from 'react';
import { styles } from '../styles/components';
import { ITERATION_LIMITS } from '../styles/constants';
import { Icons } from './Icons';

// Subscription limit notice component
export const SubscriptionLimitNotice = ({ isSubscribed, optimizationType, isMobile }) => {
  if (isSubscribed) return null;
  
  const limit = optimizationType === "classical" 
    ? ITERATION_LIMITS.unsubscribed.classical 
    : ITERATION_LIMITS.unsubscribed.quantum;
  
  const fullLimit = optimizationType === "classical" 
    ? ITERATION_LIMITS.subscribed.classical 
    : ITERATION_LIMITS.subscribed.quantum;
  
  return (
    <div style={styles.subscriptionNotice} className={isMobile ? 'mobile-stack mobile-smaller-padding' : ''}>
      <span style={styles.subscriptionNoticeIcon}>
        <Icons.warning />
      </span>
      <div className={isMobile ? 'mobile-smaller-text' : ''}>
        <strong>Free Account Limitation:</strong> Iterations capped at {limit.toLocaleString()} (vs. {fullLimit.toLocaleString()} for subscribers).{" "}
        {!isMobile ? (
          <a 
            href="#" 
            onClick={(e) => { e.preventDefault(); document.querySelector('.subscription-form').scrollIntoView({ behavior: 'smooth' }); }}
            style={{ color: "#fbbf24", textDecoration: "underline" }}
          >
            Subscribe for full capabilities
          </a>
        ) : (
          <div style={{ marginTop: '4px' }}>
            <a 
              href="#" 
              onClick={(e) => { e.preventDefault(); document.querySelector('.subscription-form').scrollIntoView({ behavior: 'smooth' }); }}
              style={{ color: "#fbbf24", textDecoration: "underline", padding: '4px 0', display: 'inline-block' }}
            >
              Subscribe for full capabilities
            </a>
          </div>
        )}
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
  setShowAdvancedParams,
  isMobile
}) => (
  <div style={styles.parametersContainer} className={`parameters-container ${isMobile ? 'mobile-smaller-padding' : ''}`}>
    <SubscriptionLimitNotice 
      isSubscribed={isSubscribed} 
      optimizationType="classical"
      isMobile={isMobile} 
    />
    
    <div className="parameter-groups-container">
      {/* Temperature - using direct structure for mobile */}
      <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
        <label 
          style={{ 
            ...styles.parameterLabel, 
            display: 'block', 
            margin: isMobile ? '0 0 2px 0' : undefined,
            padding: 0
          }}
        >
          Temperature (K):
        </label>
        <input 
          type="number" 
          value={classicalParams.temperature}
          min="1"
          max="1000"
          step="10"
          onChange={(e) => handleParamChange('classical', 'temperature', Number(e.target.value))}
          style={{ 
            ...styles.parameterInput, 
            width: isMobile ? '100%' : undefined,
            margin: 0,
            height: isMobile ? '36px' : undefined,
            padding: isMobile ? '8px' : undefined
          }}
        />
      </div>
      
      {/* Max Iterations - using direct structure for mobile */}
      <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
        <label 
          style={{ 
            ...styles.parameterLabel, 
            display: 'block', 
            margin: isMobile ? '0 0 2px 0' : undefined,
            padding: 0
          }}
        >
          Max Iterations:
        </label>
        <input 
          type="number" 
          value={classicalParams.max_iterations}
          min="100"
          max={isSubscribed ? ITERATION_LIMITS.subscribed.classical : ITERATION_LIMITS.unsubscribed.classical}
          step="100"
          onChange={(e) => handleParamChange('classical', 'max_iterations', Number(e.target.value))}
          style={{ 
            ...styles.parameterInput, 
            width: isMobile ? '100%' : undefined,
            margin: 0,
            height: isMobile ? '36px' : undefined,
            padding: isMobile ? '8px' : undefined
          }}
        />
        {!isSubscribed && classicalParams.max_iterations > ITERATION_LIMITS.unsubscribed.classical && (
          <span style={{ ...styles.warningText, marginTop: '4px', display: 'block' }}>
            <Icons.warning />
            Will be capped at {ITERATION_LIMITS.unsubscribed.classical.toLocaleString()}
          </span>
        )}
      </div>
      
      {/* Force All Iterations - using direct structure for mobile */}
      <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
        <div style={{ display: 'flex', alignItems: 'center' }}>
          <input 
            type="checkbox" 
            checked={classicalParams.force_iterations}
            onChange={(e) => handleParamChange('classical', 'force_iterations', e.target.checked)}
            style={{ 
              ...styles.parameterCheckbox,
              margin: 0,
              marginRight: '8px'
            }}
          />
          <label style={{ margin: 0 }}>Force All Iterations:</label>
        </div>
        <span style={{ 
          ...styles.parameterHelp, 
          marginLeft: isMobile ? '24px' : undefined,
          display: 'block',
          fontSize: '0.8rem'
        }}>
          Forces execution of iterations regardless of convergence
        </span>
      </div>
      
      {/* Solution Environment Selector - NEW */}
      <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
        <label 
          style={{ 
            ...styles.parameterLabel, 
            display: 'block', 
            margin: isMobile ? '0 0 2px 0' : undefined,
            padding: 0
          }}
        >
          Solution Environment:
        </label>
        <select 
          value={classicalParams.solution}
          style={{ 
            ...styles.parameterSelect, 
            width: isMobile ? '100%' : undefined,
            margin: 0,
            height: isMobile ? '36px' : undefined,
            padding: isMobile ? '8px' : undefined
          }}
        >
          <option value="water">Water</option>
          <option value="blood">Blood Plasma</option>
          <option value="intracellular">Intracellular</option>
        </select>
        <span style={{ 
          ...styles.parameterHelp, 
          marginLeft: isMobile ? '0' : undefined,
          display: 'block',
          marginTop: '4px',
          fontSize: '0.8rem'
        }}>
          Solution affects molecular interactions
        </span>
      </div>
      
      {showAdvancedParams && (
        <>
          {/* Convergence Tolerance - using direct structure for mobile */}
          <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
            <label 
              style={{ 
                ...styles.parameterLabel, 
                display: 'block', 
                margin: isMobile ? '0 0 2px 0' : undefined,
                padding: 0
              }}
            >
              Convergence Tolerance (kJ/mol):
            </label>
            <input 
              type="number" 
              value={classicalParams.tolerance}
              min="0.001"
              max="10.0"
              step="0.1"
              onChange={(e) => handleParamChange('classical', 'tolerance', Number(e.target.value))}
              style={{ 
                ...styles.parameterInput, 
                width: isMobile ? '100%' : undefined,
                margin: 0,
                height: isMobile ? '36px' : undefined,
                padding: isMobile ? '8px' : undefined
              }}
            />
          </div>
          
          {/* Bond Threshold - using direct structure for mobile */}
          <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
            <label 
              style={{ 
                ...styles.parameterLabel, 
                display: 'block', 
                margin: isMobile ? '0 0 2px 0' : undefined,
                padding: 0
              }}
            >
              Bond Threshold (nm):
            </label>
            <input 
              type="number" 
              value={classicalParams.bond_threshold}
              min="0.1"
              max="0.5"
              step="0.01"
              onChange={(e) => handleParamChange('classical', 'bond_threshold', Number(e.target.value))}
              style={{ 
                ...styles.parameterInput, 
                width: isMobile ? '100%' : undefined,
                margin: 0,
                height: isMobile ? '36px' : undefined,
                padding: isMobile ? '8px' : undefined
              }}
            />
          </div>
          
          {/* Bond Force Constant - using direct structure for mobile */}
          <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
            <label 
              style={{ 
                ...styles.parameterLabel, 
                display: 'block', 
                margin: isMobile ? '0 0 2px 0' : undefined,
                padding: 0
              }}
            >
              Bond Force Constant (kJ/mol/nm²):
            </label>
            <input 
              type="number" 
              value={classicalParams.bond_force_constant}
              min="100"
              max="10000"
              step="100"
              onChange={(e) => handleParamChange('classical', 'bond_force_constant', Number(e.target.value))}
              style={{ 
                ...styles.parameterInput, 
                width: isMobile ? '100%' : undefined,
                margin: 0,
                height: isMobile ? '36px' : undefined,
                padding: isMobile ? '8px' : undefined
              }}
            />
          </div>
          
          {/* Angle Force Constant - using direct structure for mobile */}
          <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
            <label 
              style={{ 
                ...styles.parameterLabel, 
                display: 'block', 
                margin: isMobile ? '0 0 2px 0' : undefined,
                padding: 0
              }}
            >
              Angle Force Constant (kJ/mol/rad²):
            </label>
            <input 
              type="number" 
              value={classicalParams.angle_force_constant}
              min="50"
              max="5000"
              step="50"
              onChange={(e) => handleParamChange('classical', 'angle_force_constant', Number(e.target.value))}
              style={{ 
                ...styles.parameterInput, 
                width: isMobile ? '100%' : undefined,
                margin: 0,
                height: isMobile ? '36px' : undefined,
                padding: isMobile ? '8px' : undefined
              }}
            />
          </div>

          {/* New Bond Visualization Settings Section */}
          <div style={{ borderTop: '1px solid rgba(255, 255, 255, 0.1)', marginTop: '12px', paddingTop: '12px' }}>
            <div style={{ marginBottom: '8px', fontWeight: 'bold' }}>Bond Visualization Settings:</div>
            
            {/* Covalent Bond Threshold */}
            <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
              <label 
                style={{ 
                  ...styles.parameterLabel, 
                  display: 'block', 
                  margin: isMobile ? '0 0 2px 0' : undefined,
                  padding: 0
                }}
              >
                Covalent Bond Display Length (Å):
              </label>
              <input 
                type="number" 
                value={classicalParams.covalent_display_threshold}
                min="0.1"
                max="3.0"
                step="0.1"
                onChange={(e) => handleParamChange('classical', 'covalent_display_threshold', Number(e.target.value))}
                style={{ 
                  ...styles.parameterInput, 
                  width: isMobile ? '100%' : undefined,
                  margin: 0,
                  height: isMobile ? '36px' : undefined,
                  padding: isMobile ? '8px' : undefined
                }}
              />
            </div>
            
            {/* Hydrogen Bond Threshold */}
            <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
              <label 
                style={{ 
                  ...styles.parameterLabel, 
                  display: 'block', 
                  margin: isMobile ? '0 0 2px 0' : undefined,
                  padding: 0
                }}
              >
                Hydrogen Bond Display Length (Å):
              </label>
              <input 
                type="number" 
                value={classicalParams.hydrogen_display_threshold}
                min="1.5"
                max="4.0"
                step="0.1"
                onChange={(e) => handleParamChange('classical', 'hydrogen_display_threshold', Number(e.target.value))}
                style={{ 
                  ...styles.parameterInput, 
                  width: isMobile ? '100%' : undefined,
                  margin: 0,
                  height: isMobile ? '36px' : undefined,
                  padding: isMobile ? '8px' : undefined
                }}
              />
            </div>
            
            {/* Angle Display Threshold */}
            <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
              <label 
                style={{ 
                  ...styles.parameterLabel, 
                  display: 'block', 
                  margin: isMobile ? '0 0 2px 0' : undefined,
                  padding: 0
                }}
              >
                Angle Display Threshold (degrees):
              </label>
              <input 
                type="number" 
                value={classicalParams.angle_display_threshold}
                min="5"
                max="180"
                step="5"
                onChange={(e) => handleParamChange('classical', 'angle_display_threshold', Number(e.target.value))}
                style={{ 
                  ...styles.parameterInput, 
                  width: isMobile ? '100%' : undefined,
                  margin: 0,
                  height: isMobile ? '36px' : undefined,
                  padding: isMobile ? '8px' : undefined
                }}
              />
            </div>
          </div>
        </>
      )}
    </div>
    
    <div style={styles.buttonGroup} className={isMobile ? 'mobile-stack' : ''}>
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
        className={isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}
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
        className={isMobile ? 'mobile-full-width' : ''}
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
  setShowAdvancedParams,
  isMobile
}) => (
  <div style={styles.parametersContainer} className={`parameters-container ${isMobile ? 'mobile-smaller-padding' : ''}`}>
    <SubscriptionLimitNotice 
      isSubscribed={isSubscribed} 
      optimizationType="quantum"
      isMobile={isMobile}
    />
    
    <div className="parameter-groups-container">
      {/* Basis Set - using direct structure for mobile */}
      <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
        <label 
          style={{ 
            ...styles.parameterLabel, 
            display: 'block', 
            margin: isMobile ? '0 0 2px 0' : undefined,
            padding: 0
          }}
        >
          Basis Set:
        </label>
        <select 
          value={quantumParams.basis}
          onChange={(e) => handleParamChange('quantum', 'basis', e.target.value)}
          style={{ 
            ...styles.parameterSelect, 
            width: isMobile ? '100%' : undefined,
            margin: 0,
            height: isMobile ? '36px' : undefined,
            padding: isMobile ? '8px' : undefined
          }}
        >
          <option value="sto-3g">STO-3G (Minimal)</option>
          <option value="6-31g">6-31G (Standard)</option>
          {isSubscribed && <option value="6-311g">6-311G (Extended)</option>}
          {isSubscribed && <option value="cc-pvdz">cc-pVDZ (Double Zeta)</option>}
        </select>
        {!isSubscribed && (quantumParams.basis === "6-311g" || quantumParams.basis === "cc-pvdz") && (
          <span style={{ ...styles.dangerBadge, marginTop: '4px', display: 'block' }}>
            Extended basis sets require subscription
          </span>
        )}
      </div>
      
      {/* Max Iterations - using direct structure for mobile */}
      <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
        <label 
          style={{ 
            ...styles.parameterLabel, 
            display: 'block', 
            margin: isMobile ? '0 0 2px 0' : undefined,
            padding: 0
          }}
        >
          Max Iterations:
        </label>
        <input 
          type="number" 
          value={quantumParams.max_iterations}
          min="1"
          max={isSubscribed ? ITERATION_LIMITS.subscribed.quantum : ITERATION_LIMITS.unsubscribed.quantum}
          onChange={(e) => handleParamChange('quantum', 'max_iterations', Number(e.target.value))}
          style={{ 
            ...styles.parameterInput, 
            width: isMobile ? '100%' : undefined,
            margin: 0,
            height: isMobile ? '36px' : undefined,
            padding: isMobile ? '8px' : undefined
          }}
        />
        {!isSubscribed && quantumParams.max_iterations > ITERATION_LIMITS.unsubscribed.quantum && (
          <span style={{ ...styles.warningText, marginTop: '4px', display: 'block' }}>
            <Icons.warning />
            Will be capped at {ITERATION_LIMITS.unsubscribed.quantum}
          </span>
        )}
        {isSubscribed && quantumParams.max_iterations > ITERATION_LIMITS.subscribed.quantum && (
          <span style={{ ...styles.warningText, marginTop: '4px', display: 'block' }}>
            <Icons.warning />
            Will be capped at {ITERATION_LIMITS.subscribed.quantum}
          </span>
        )}
      </div>
      
      {/* Solution Environment Selector - NEW */}
      <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
        <label 
          style={{ 
            ...styles.parameterLabel, 
            display: 'block', 
            margin: isMobile ? '0 0 2px 0' : undefined,
            padding: 0
          }}
        >
          Solution Environment:
        </label>
        <select 
          value={quantumParams.solution}
          style={{ 
            ...styles.parameterSelect, 
            width: isMobile ? '100%' : undefined,
            margin: 0,
            height: isMobile ? '36px' : undefined,
            padding: isMobile ? '8px' : undefined
          }}
        >
          <option value="water">Water (pH 7.0)</option>
          <option value="blood">Blood Plasma (pH 7.4)</option>
          <option value="intracellular">Intracellular (pH 7.2)</option>
        </select>
        <span style={{ 
          ...styles.parameterHelp, 
          marginLeft: isMobile ? '0' : undefined,
          display: 'block',
          marginTop: '4px',
          fontSize: '0.8rem'
        }}>
          Solution affects molecular interactions
        </span>
      </div>
      
      {showAdvancedParams && (
        <>
          {/* Convergence Threshold - using direct structure for mobile */}
          <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
            <label 
              style={{ 
                ...styles.parameterLabel, 
                display: 'block', 
                margin: isMobile ? '0 0 2px 0' : undefined,
                padding: 0
              }}
            >
              Convergence Threshold:
            </label>
            <input 
              type="number" 
              value={quantumParams.convergence_threshold}
              min="0.000001"
              max="0.01"
              step="0.000001"
              onChange={(e) => handleParamChange('quantum', 'convergence_threshold', Number(e.target.value))}
              style={{ 
                ...styles.parameterInput, 
                width: isMobile ? '100%' : undefined,
                margin: 0,
                height: isMobile ? '36px' : undefined,
                padding: isMobile ? '8px' : undefined
              }}
            />
          </div>
          
          {/* Step Size - using direct structure for mobile */}
          <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
            <label 
              style={{ 
                ...styles.parameterLabel, 
                display: 'block', 
                margin: isMobile ? '0 0 2px 0' : undefined,
                padding: 0
              }}
            >
              Step Size:
            </label>
            <input 
              type="number" 
              value={quantumParams.step_size}
              min="0.01"
              max="1.0"
              step="0.01"
              onChange={(e) => handleParamChange('quantum', 'step_size', Number(e.target.value))}
              style={{ 
                ...styles.parameterInput, 
                width: isMobile ? '100%' : undefined,
                margin: 0,
                height: isMobile ? '36px' : undefined,
                padding: isMobile ? '8px' : undefined
              }}
            />
          </div>

          {/* New Bond Visualization Settings Section */}
          <div style={{ borderTop: '1px solid rgba(255, 255, 255, 0.1)', marginTop: '12px', paddingTop: '12px' }}>
            <div style={{ marginBottom: '8px', fontWeight: 'bold' }}>Bond Visualization Settings:</div>
            
            {/* Covalent Bond Threshold */}
            <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
              <label 
                style={{ 
                  ...styles.parameterLabel, 
                  display: 'block', 
                  margin: isMobile ? '0 0 2px 0' : undefined,
                  padding: 0
                }}
              >
                Covalent Bond Display Length (Å):
              </label>
              <input 
                type="number" 
                value={quantumParams.covalent_display_threshold}
                min="0.1"
                max="3.0"
                step="0.1"
                onChange={(e) => handleParamChange('quantum', 'covalent_display_threshold', Number(e.target.value))}
                style={{ 
                  ...styles.parameterInput, 
                  width: isMobile ? '100%' : undefined,
                  margin: 0,
                  height: isMobile ? '36px' : undefined,
                  padding: isMobile ? '8px' : undefined
                }}
              />
            </div>
            
            {/* Hydrogen Bond Threshold */}
            <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
              <label 
                style={{ 
                  ...styles.parameterLabel, 
                  display: 'block', 
                  margin: isMobile ? '0 0 2px 0' : undefined,
                  padding: 0
                }}
              >
                Hydrogen Bond Display Length (Å):
              </label>
              <input 
                type="number" 
                value={quantumParams.hydrogen_display_threshold}
                min="1.5"
                max="4.0"
                step="0.1"
                onChange={(e) => handleParamChange('quantum', 'hydrogen_display_threshold', Number(e.target.value))}
                style={{ 
                  ...styles.parameterInput, 
                  width: isMobile ? '100%' : undefined,
                  margin: 0,
                  height: isMobile ? '36px' : undefined,
                  padding: isMobile ? '8px' : undefined
                }}
              />
            </div>
            
            {/* Angle Display Threshold */}
            <div className="parameter-group-wrapper" style={{ margin: isMobile ? '0 0 8px 0' : undefined }}>
              <label 
                style={{ 
                  ...styles.parameterLabel, 
                  display: 'block', 
                  margin: isMobile ? '0 0 2px 0' : undefined,
                  padding: 0
                }}
              >
                Angle Display Threshold (degrees):
              </label>
              <input 
                type="number" 
                value={quantumParams.angle_display_threshold}
                min="5"
                max="180"
                step="5"
                onChange={(e) => handleParamChange('quantum', 'angle_display_threshold', Number(e.target.value))}
                style={{ 
                  ...styles.parameterInput, 
                  width: isMobile ? '100%' : undefined,
                  margin: 0,
                  height: isMobile ? '36px' : undefined,
                  padding: isMobile ? '8px' : undefined
                }}
              />
            </div>
          </div>
        </>
      )}
    </div>
    
    <div style={styles.buttonGroup} className={isMobile ? 'mobile-stack' : ''}>
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
        className={isMobile ? 'mobile-full-width mobile-margin-bottom' : ''}
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
        className={isMobile ? 'mobile-full-width' : ''}
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