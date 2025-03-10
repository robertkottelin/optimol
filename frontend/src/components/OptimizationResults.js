import React from 'react';
import { styles } from '../styles/components';
import { Icons } from './Icons';
import { COLORS } from '../styles/constants';

const OptimizationResults = ({ optimizationResult, optimizationType, handleDownload }) => {
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
      
      <div style={{ marginTop: "24px" }}>
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
            borderRadius: "8px",
            fontWeight: "600",
          }}
        >
          <span><Icons.download /></span>
          Download Result
        </button>
      </div>
    </div>
  );
};

export default OptimizationResults;