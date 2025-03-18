import React from 'react';
import { styles } from '../styles/components';
import { Icons } from './Icons';
import { COLORS, FONTS, SPACING } from '../styles/constants';

const OptimizationResults = ({ optimizationResult, optimizationType, handleDownload, isMobile }) => {
  if (!optimizationResult || !optimizationResult.result) return null;
  
  const result = optimizationResult.result;
  
  if (result.error) {
    return (
      <div style={styles.resultsContainer} className={isMobile ? 'mobile-smaller-padding' : ''}>
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
      <div style={styles.resultsContainer} className={isMobile ? 'mobile-smaller-padding' : ''}>
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
  
  const ResultItem = ({ label, value }) => (
    <div 
      style={styles.resultItem} 
      className={isMobile ? 'mobile-stack' : ''}
    >
      <span style={styles.resultLabel} className={isMobile ? 'mobile-full-width' : ''}>
        {label}:
      </span>
      <span 
        style={styles.resultValue} 
        className={isMobile ? 'mobile-full-width' : ''}
      >
        {value}
      </span>
    </div>
  );
  
  // Determine if we're in interaction mode (multiple molecules)
  const isInteractionMode = result.metadata && 
    (result.molecule1_optimized_atoms || result.molecule2_optimized_atoms || 
     result.metadata.molecules > 1);
  
  return (
    <div style={styles.resultsContainer} className={isMobile ? 'mobile-smaller-padding' : ''}>
      <h3 style={styles.resultTitle}>
        <span style={{ ...styles.resultIcon, color: optimizationType === "classical" ? COLORS.classical : COLORS.quantum }}>
          {optimizationType === "classical" ? <Icons.classical /> : <Icons.quantum />}
        </span>
        {optimizationType === "classical" ? "Classical" : "Quantum"} Optimization Results
      </h3>
      
      <ResultItem 
        label="Method" 
        value={result.metadata.method} 
      />
      
      <ResultItem 
        label="Library" 
        value={result.metadata.library} 
      />
      
      {optimizationType === "classical" ? (
        <>
          <ResultItem 
            label="Temperature" 
            value={result.metadata.parameters?.temperature !== undefined ? 
                `${result.metadata.parameters.temperature} K` : 
                "N/A"} 
          />
          
          <ResultItem 
            label="Iterations" 
            value={result.metadata.iterations_performed !== undefined ? 
                (result.metadata.parameters.force_iterations ? 
                  `${result.metadata.iterations_performed} (forced)` : 
                  `${result.metadata.iterations_performed} (converged)`) : 
                "N/A"} 
          />
          
          <ResultItem 
            label="Initial Energy" 
            value={result.metadata.initial_energy_kj_mol !== undefined ? 
                `${result.metadata.initial_energy_kj_mol.toFixed(4)} kJ/mol` : 
                "N/A"} 
          />
          
          <ResultItem 
            label="Final Energy" 
            value={result.metadata.final_energy_kj_mol !== undefined ? 
                `${result.metadata.final_energy_kj_mol.toFixed(4)} kJ/mol` : 
                "N/A"} 
          />
          
          <ResultItem 
            label="Energy Change" 
            value={result.metadata.energy_change_kj_mol !== undefined ? 
                `${result.metadata.energy_change_kj_mol.toFixed(4)} kJ/mol` : 
                "N/A"} 
          />
          
          <ResultItem 
            label="Bonds Detected" 
            value={result.metadata.bonds_detected !== undefined ? 
                result.metadata.bonds_detected : 
                "N/A"} 
          />
          
          <ResultItem 
            label="Angles Detected" 
            value={result.metadata.angles_detected !== undefined ? 
                result.metadata.angles_detected : 
                "N/A"} 
          />
        </>
      ) : (
          <>
          <ResultItem 
            label="Basis Set" 
            value={result.metadata.parameters?.basis !== undefined ?
                result.metadata.parameters.basis :
                "N/A"} 
          />
          
          <ResultItem 
            label="Theory Level" 
            value={result.metadata.theory_level !== undefined ?
                result.metadata.theory_level :
                "N/A"} 
          />
          
          <ResultItem 
            label="Final Energy" 
            value={result.metadata.final_energy_hartree !== undefined ? 
                `${result.metadata.final_energy_hartree.toFixed(6)} Hartree` : 
                "N/A"} 
          />
          
          <ResultItem 
            label="Iterations" 
            value={result.metadata.iterations !== undefined ?
                result.metadata.iterations :
                "N/A"} 
          />
          
          <ResultItem 
            label="Converged" 
            value={result.metadata.converged !== undefined ?
                (result.metadata.converged ? "Yes" : "No") :
                "N/A"} 
          />
        </>
      )}
      
      <ResultItem 
        label="Duration" 
        value={result.metadata.duration_seconds !== undefined ?
            `${result.metadata.duration_seconds.toFixed(2)} seconds` :
            "N/A"} 
      />
      
      {/* Molecular Interaction Section */}
      {isInteractionMode && (
        <>
          <h4 style={{
            marginTop: SPACING.lg,
            marginBottom: SPACING.sm,
            fontWeight: FONTS.weightSemiBold,
            borderBottom: `1px solid ${COLORS.border}`,
            paddingBottom: SPACING.xs
          }}>
            Molecular Interaction
          </h4>
          
          <ResultItem 
            label="Molecules" 
            value={result.metadata.molecules || "2"} 
          />
          
          <ResultItem 
            label="Molecule 1 Atoms" 
            value={result.metadata.molecule1_atom_count || 
              (result.molecule1_optimized_atoms ? result.molecule1_optimized_atoms.length : "N/A")} 
          />
          
          <ResultItem 
            label="Molecule 2 Atoms" 
            value={result.metadata.molecule2_atom_count || 
              (result.molecule2_optimized_atoms ? result.molecule2_optimized_atoms.length : "N/A")} 
          />
          
          {optimizationType === "classical" && (
            <>
              <ResultItem 
                label="Intermolecular Bonds" 
                value={result.metadata.intermolecular_bonds || "N/A"} 
              />
              
              <ResultItem 
                label="Interaction Energy" 
                value={result.metadata.interaction_energy_kj_mol !== undefined ? 
                    `${result.metadata.interaction_energy_kj_mol.toFixed(4)} kJ/mol` : 
                    "N/A"} 
              />
            </>
          )}
          
          {optimizationType === "quantum" && (
            <ResultItem 
              label="Interaction Energy" 
              value={result.metadata.interaction_energy_hartree !== undefined ? 
                  `${result.metadata.interaction_energy_hartree.toFixed(6)} Hartree` : 
                  "N/A"} 
            />
          )}
        </>
      )}
      
      <div style={{ marginTop: "24px", textAlign: isMobile ? 'center' : 'left' }}>
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
            width: isMobile ? '100%' : 'auto',
            justifyContent: isMobile ? 'center' : 'flex-start',
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