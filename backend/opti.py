from flask import Blueprint, request, jsonify
import json
import numpy as np
import logging
from datetime import datetime
from extensions import db
from user import User

# Classical optimization imports
import openmm as mm
import openmm.app as app
from openmm import unit

# Quantum optimization imports
from pyscf import gto, scf, grad
from pyscf import df  # Import density fitting module

from flask_jwt_extended import jwt_required, get_jwt_identity, verify_jwt_in_request
from constants import ITERATION_LIMITS

# Configure logging
logger = logging.getLogger(__name__)

# Create blueprint
opti_bp = Blueprint('opti', __name__)

# Atomic parameter constants
ELEMENT_MASSES = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
    'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904
}

ELEMENT_NONBONDED_PARAMS = {
    # Element -> [charge, sigma (nm), epsilon (kJ/mol)]
    'H': [0.0, 0.106, 0.0656],
    'C': [0.0, 0.340, 0.4577],
    'N': [0.0, 0.325, 0.7113],
    'O': [0.0, 0.296, 0.8786],
    'F': [0.0, 0.312, 0.255],
    'P': [0.0, 0.374, 0.8368],
    'S': [0.0, 0.356, 1.046],
    'Cl': [0.0, 0.347, 1.1087],
    'Br': [0.0, 0.391, 0.8368],
    'I': [0.0, 0.425, 0.6699]
}

# Default nonbonded parameters for unknown elements
DEFAULT_NONBONDED_PARAMS = [0.0, 0.3, 0.5]

# System size thresholds for quantum calculations based on basis set complexity
QUANTUM_ATOM_THRESHOLDS = {
    "sto-3g": 100,    # Most minimal basis
    "3-21g": 70,      
    "6-31g": 50,      # Medium basis
    "6-311g": 30,     # Larger basis
    "cc-pvdz": 25,    # Even larger basis
    "cc-pvtz": 15     # Very large basis
}

class Optimization(db.Model):
    """Optimization record model for tracking molecular optimizations."""
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=True)  # Modified to allow null for non-subscribers
    optimization_type = db.Column(db.String(50), nullable=False)  # 'classical' or 'quantum'
    parameters = db.Column(db.Text, nullable=False)  # JSON string of parameters
    result = db.Column(db.Text, nullable=True)  # JSON string of the result
    created_at = db.Column(db.DateTime, default=db.func.current_timestamp())

    def __repr__(self):
        return f'<Optimization {self.id} for User {self.user_id}>'


def get_element_mass(element):
    """Safely retrieve element mass with fallback for unknown elements."""
    return ELEMENT_MASSES.get(element, 0.0)


def get_element_nonbonded_params(element):
    """Safely retrieve nonbonded parameters with fallback for unknown elements."""
    return ELEMENT_NONBONDED_PARAMS.get(element, DEFAULT_NONBONDED_PARAMS)


def calculate_angle(vec1, vec2):
    """
    Calculate angle between two vectors with improved numerical stability.
    
    Args:
        vec1: First vector
        vec2: Second vector
        
    Returns:
        Angle in radians
    """
    norm1 = np.sqrt(np.sum(vec1*vec1))
    norm2 = np.sqrt(np.sum(vec2*vec2))
    
    # Prevent division by zero
    if norm1 < 1e-10 or norm2 < 1e-10:
        return 0.0
    
    dot = np.sum(vec1*vec2)
    cosine = dot / (norm1 * norm2)
    
    # Ensure numerical stability
    cosine = max(-1.0 + 1e-10, min(1.0 - 1e-10, cosine))
    
    return np.arccos(cosine)


def are_molecules_similar(molecule1_atoms, molecule2_atoms, position_threshold=0.01, element_match=True):
    """
    Determine if two molecules are identical or very similar.
    
    Args:
        molecule1_atoms: List of atoms for first molecule
        molecule2_atoms: List of atoms for second molecule
        position_threshold: Threshold (in Å) for considering positions similar
        element_match: Whether to require elements to match
    
    Returns:
        tuple: (is_similar, rmsd) 
               is_similar: Boolean indicating if molecules are similar
               rmsd: Root mean square deviation between atom positions (None if different lengths)
    """
    # Quick check: different number of atoms means different molecules
    if len(molecule1_atoms) != len(molecule2_atoms):
        return False, None
    
    positions1 = np.array([[atom['x'], atom['y'], atom['z']] for atom in molecule1_atoms])
    positions2 = np.array([[atom['x'], atom['y'], atom['z']] for atom in molecule2_atoms])
    
    # Check if elements match
    if element_match:
        elements1 = [atom['element'] for atom in molecule1_atoms]
        elements2 = [atom['element'] for atom in molecule2_atoms]
        if elements1 != elements2:
            return False, None
    
    # Calculate RMSD between molecules
    squared_diff = np.sum((positions1 - positions2)**2, axis=1)
    rmsd = np.sqrt(np.mean(squared_diff))
    
    # Consider molecules similar if RMSD is below threshold
    is_similar = rmsd < position_threshold
    
    return is_similar, rmsd


def displace_molecule(atoms, displacement_vector=None, magnitude=1.0):
    """
    Apply displacement to molecule atoms.
    
    Args:
        atoms: List of atom dictionaries to displace
        displacement_vector: Optional specific displacement vector (default: random unit vector)
        magnitude: Magnitude of displacement in Å
    
    Returns:
        List of displaced atom dictionaries
    """
    if displacement_vector is None:
        # Generate random unit vector
        displacement_vector = np.random.randn(3)
        displacement_vector = displacement_vector / np.linalg.norm(displacement_vector)
    
    # Scale to desired magnitude
    displacement_vector = displacement_vector * magnitude
    
    # Apply displacement to all atoms
    displaced_atoms = []
    for atom in atoms:
        displaced_atom = atom.copy()
        displaced_atom['x'] = atom['x'] + displacement_vector[0]
        displaced_atom['y'] = atom['y'] + displacement_vector[1]
        displaced_atom['z'] = atom['z'] + displacement_vector[2]
        displaced_atoms.append(displaced_atom)
    
    return displaced_atoms


def optimize_classical(atoms, params=None):
    """
    Optimize molecule using classical molecular dynamics with OpenMM.
    Template-free approach for arbitrary molecules.
    """
    try:   
        start_time = datetime.now()
        
        # Set default parameters if not provided
        if params is None:
            params = {}
            
        # Extract parameters with defaults
        temperature = params.get("temperature", 300)  # Kelvin
        max_iterations = params.get("max_iterations", 1000)  # For energy minimization
        bond_threshold = params.get("bond_threshold", 0.2)  # nm
        bond_force_constant = params.get("bond_force_constant", 1000.0)  # kJ/mol/nm^2
        angle_force_constant = params.get("angle_force_constant", 500.0)  # kJ/mol/radian^2
        tolerance = params.get("tolerance", 10.0)  # kJ/mol/nm (force units)
        force_iterations = params.get("force_iterations", False)  # Whether to force all iterations
        
        # Create system from atoms
        positions = []
        elements = []
        
        for atom in atoms:
            positions.append([atom["x"], atom["y"], atom["z"]])
            elements.append(atom["element"])
            
        positions = np.array(positions) * unit.angstrom
        
        # Create topology
        topology = app.Topology()
        chain = topology.addChain()
        residue = topology.addResidue('MOL', chain)
        
        # Add atoms with proper masses
        atom_objects = []
        for i, element in enumerate(elements):
            atom_obj = topology.addAtom(element, app.Element.getBySymbol(element), residue)
            atom_objects.append(atom_obj)
        
        # Create system
        system = mm.System()
        
        # Add atom masses to the system
        for i, element in enumerate(elements):
            mass = get_element_mass(element)
            system.addParticle(mass)
        
        # Determine bonds based on distance
        positions_nm = positions.value_in_unit(unit.nanometer)
        bonds = []
        
        for i in range(len(atom_objects)):
            for j in range(i+1, len(atom_objects)):
                dist = np.sqrt(np.sum((positions_nm[i] - positions_nm[j])**2))
                if dist < bond_threshold:  # Bond threshold in nm (now parameterized)
                    topology.addBond(atom_objects[i], atom_objects[j])
                    bonds.append((i, j, dist))
        
        # Add HarmonicBondForce
        bond_force = mm.HarmonicBondForce()
        for i, j, dist in bonds:
            # Use parameterized force constant
            bond_force.addBond(i, j, dist, bond_force_constant)
        system.addForce(bond_force)
        
        # Add HarmonicAngleForce
        angle_force = mm.HarmonicAngleForce()
        # Find all angle triplets i-j-k where i-j and j-k are bonded
        bond_partners = {}
        for i, j, _ in bonds:
            if i not in bond_partners:
                bond_partners[i] = []
            if j not in bond_partners:
                bond_partners[j] = []
            bond_partners[i].append(j)
            bond_partners[j].append(i)
        
        angles = []
        for j in range(len(atom_objects)):
            if j in bond_partners:
                partners = bond_partners[j]
                for i in range(len(partners)):
                    for k in range(i+1, len(partners)):
                        atom_i = partners[i]
                        atom_k = partners[k]
                        # Calculate current angle with improved numerical stability
                        vec1 = positions_nm[atom_i] - positions_nm[j]
                        vec2 = positions_nm[atom_k] - positions_nm[j]
                        angle = calculate_angle(vec1, vec2)
                        
                        # Use parameterized force constant
                        angle_force.addAngle(atom_i, j, atom_k, angle, angle_force_constant)
                        angles.append((atom_i, j, atom_k))
        
        system.addForce(angle_force)
        
        # Add NonbondedForce for non-bonded interactions
        nonbonded_force = mm.NonbondedForce()
        
        for i, element in enumerate(elements):
            params = get_element_nonbonded_params(element)
            nonbonded_force.addParticle(params[0], params[1], params[2])
        
        # Track exception pairs to avoid duplicates
        exception_pairs = set()
        
        # Create 1-2 exclusions for bonded atoms
        for i, j, _ in bonds:
            pair = (min(i, j), max(i, j))
            if pair not in exception_pairs:
                nonbonded_force.addException(i, j, 0.0, 0.1, 0.0)
                exception_pairs.add(pair)
        
        # Add 1-3 angle exclusions
        for i, j, k in angles:
            pair = (min(i, k), max(i, k))
            if pair not in exception_pairs:
                nonbonded_force.addException(i, k, 0.0, 0.1, 0.0)
                exception_pairs.add(pair)
            
        system.addForce(nonbonded_force)
        
        # Create Integrator with parameterized temperature
        integrator = mm.LangevinIntegrator(temperature*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        
        # Create simulation
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        
        # Get initial energy
        state = simulation.context.getState(getEnergy=True)
        start_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Custom iteration tracking
        iterations_performed = 0
        
        if force_iterations:
            # Manual optimization to force all iterations
            energy_history = []
            
            # Initial minimization to stabilize structure
            simulation.minimizeEnergy(maxIterations=100)
            
            # Create optimization simulation with small time step
            integrator_opt = mm.VerletIntegrator(0.001*unit.picoseconds)
            simulation_opt = app.Simulation(topology, system, integrator_opt)
            simulation_opt.context.setPositions(simulation.context.getState(getPositions=True).getPositions())
            
            # Execute exact number of iterations requested
            for i in range(max_iterations):
                simulation_opt.step(1)
                
                state = simulation_opt.context.getState(getEnergy=True)
                energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                energy_history.append(energy)
                iterations_performed += 1
                
                # More robust safety check for complex molecules
                if (i > 20 and np.isfinite(energy) and energy > energy_history[i-1] * 5.0) or not np.isfinite(energy):
                    logger.info(f"Breaking at iteration {i} due to energy instability: {energy}")
                    break
                    
            state = simulation_opt.context.getState(getPositions=True, getEnergy=True)
            minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
            final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
        else:
            # Standard energy minimization with correct force units for tolerance
            # Convert tolerance from kJ/mol/nm to force units as expected by OpenMM
            force_tolerance = tolerance * unit.kilojoule_per_mole / unit.nanometer
            
            # Skip passing tolerance if using default value to avoid unit errors
            if abs(tolerance - 10.0) < 0.01:  # If close to default value
                simulation.minimizeEnergy(maxIterations=max_iterations)
            else:
                # Use custom tolerance
                simulation.minimizeEnergy(tolerance=force_tolerance, maxIterations=max_iterations)
            
            # Get final state
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            energy_change = abs(start_energy - final_energy)
            
            # Estimate iterations based on energy change and tolerance
            # This is an approximation since OpenMM doesn't expose actual iteration count
            iterations_performed = min(max_iterations, int((energy_change/max(tolerance, 1e-6)) * 10))
            
            minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        
        # Format result with optimized coordinates
        optimized_atoms = []
        for i, atom in enumerate(atoms):
            optimized_atoms.append({
                "id": atom.get("id", i+1),
                "element": atom["element"],
                "x": float(minimized_positions[i][0]),
                "y": float(minimized_positions[i][1]),
                "z": float(minimized_positions[i][2])
            })
            
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        # Ensure we clean up the simulation objects to release GPU resources
        del simulation
        if force_iterations:
            del simulation_opt
            
        return {
            "optimized_atoms": optimized_atoms,
            "metadata": {
                "method": "classical_molecular_dynamics",
                "library": "OpenMM",
                "parameters": {
                    "temperature": temperature,
                    "max_iterations": max_iterations,
                    "bond_threshold": bond_threshold,
                    "bond_force_constant": bond_force_constant,
                    "angle_force_constant": angle_force_constant,
                    "tolerance": tolerance,
                    "force_iterations": force_iterations
                },
                "iterations_performed": iterations_performed,
                "final_energy_kj_mol": final_energy,
                "initial_energy_kj_mol": start_energy,
                "energy_change_kj_mol": abs(start_energy - final_energy),
                "duration_seconds": duration,
                "convergence": "energy_minimized",
                "bonds_detected": len(bonds),
                "angles_detected": len(angles)
            }
        }
        
    except Exception as e:
        logger.error(f"Classical optimization error: {str(e)}", exc_info=True)
        return {
            "error": str(e),
            "metadata": {
                "method": "classical_molecular_dynamics",
                "library": "OpenMM",
                "status": "failed"
            }
        }   


def optimize_classical_combined(molecule1_atoms, molecule2_atoms, params=None):
    """
    Optimize combined molecule system using classical molecular dynamics with OpenMM.
    
    Args:
        molecule1_atoms: Atoms from the first molecule (can be None)
        molecule2_atoms: Atoms from the second molecule (can be None)
        params: Optimization parameters
        
    Returns:
        Dictionary with optimized atoms and metadata
    """
    try:
        start_time = datetime.now()
        
        # Set default parameters if not provided
        if params is None:
            params = {}
            
        # Extract parameters with defaults
        temperature = params.get("temperature", 300)  # Kelvin
        max_iterations = params.get("max_iterations", 1000)  # For energy minimization
        bond_threshold = params.get("bond_threshold", 0.2)  # nm
        bond_force_constant = params.get("bond_force_constant", 1000.0)  # kJ/mol/nm^2
        angle_force_constant = params.get("angle_force_constant", 500.0)  # kJ/mol/radian^2
        tolerance = params.get("tolerance", 10.0)  # kJ/mol/nm (force units)
        force_iterations = params.get("force_iterations", False)  # Whether to force all iterations
        
        # Check if molecules are too similar (causes numerical instability)
        is_similar = False
        rmsd = None
        
        if molecule1_atoms and molecule2_atoms:
            is_similar, rmsd = are_molecules_similar(molecule1_atoms, molecule2_atoms)
            
            if is_similar:
                logger.info(f"Detected identical or very similar molecules (RMSD={rmsd}Å), applying displacement")
                # Apply a 3Å displacement to second molecule to avoid instability
                molecule2_atoms = displace_molecule(molecule2_atoms, magnitude=3.0)
        
        # Combine atoms from both molecules, tracking their origins
        all_atoms = []
        atom_origins = []  # Track which molecule each atom came from (1 or 2)
        
        # Add atoms from molecule 1 if present
        if molecule1_atoms is not None:
            for atom in molecule1_atoms:
                all_atoms.append(atom)
                atom_origins.append(1)
                
        # Add atoms from molecule 2 if present
        if molecule2_atoms is not None:
            for atom in molecule2_atoms:
                all_atoms.append(atom)
                atom_origins.append(2)
        
        # Create system from atoms
        positions = []
        elements = []
        
        for atom in all_atoms:
            positions.append([atom["x"], atom["y"], atom["z"]])
            elements.append(atom["element"])
            
        positions = np.array(positions) * unit.angstrom
        
        # Create topology
        topology = app.Topology()
        chain = topology.addChain()
        residue = topology.addResidue('MOL', chain)
        
        # Add atoms with proper masses
        atom_objects = []
        for i, element in enumerate(elements):
            atom_obj = topology.addAtom(element, app.Element.getBySymbol(element), residue)
            atom_objects.append(atom_obj)
        
        # Create system
        system = mm.System()
        
        # Add atom masses to the system
        for i, element in enumerate(elements):
            mass = get_element_mass(element)
            system.addParticle(mass)
        
        # Determine bonds based on distance with additional safeguards
        positions_nm = positions.value_in_unit(unit.nanometer)
        bond_set = set()  # Track unique bond pairs
        bonds = []
        intermolecular_bonds = 0  # Counter for bonds between molecules
        
        for i in range(len(atom_objects)):
            for j in range(i+1, len(atom_objects)):  # Ensure i < j to avoid duplicates
                dist = np.sqrt(np.sum((positions_nm[i] - positions_nm[j])**2))
                if dist < bond_threshold:  # Bond threshold in nm
                    # Ensure we haven't already added this bond
                    if (i, j) not in bond_set:
                        topology.addBond(atom_objects[i], atom_objects[j])
                        bonds.append((i, j, dist))
                        bond_set.add((i, j))
                        
                        # Check if this bond connects atoms from different molecules
                        if len(atom_origins) > max(i, j) and atom_origins[i] != atom_origins[j]:
                            intermolecular_bonds += 1
        
        # Add HarmonicBondForce
        bond_force = mm.HarmonicBondForce()
        for i, j, dist in bonds:
            # Use parameterized force constant
            bond_force.addBond(i, j, dist, bond_force_constant)
        system.addForce(bond_force)
        
        # Add HarmonicAngleForce
        angle_force = mm.HarmonicAngleForce()
        # Find all angle triplets i-j-k where i-j and j-k are bonded
        bond_partners = {}
        for i, j, _ in bonds:
            if i not in bond_partners:
                bond_partners[i] = []
            if j not in bond_partners:
                bond_partners[j] = []
            bond_partners[i].append(j)
            bond_partners[j].append(i)
        
        angles = []
        angle_set = set()  # Track unique angle triplets
        for j in range(len(atom_objects)):
            if j in bond_partners:
                partners = bond_partners[j]
                for i_idx in range(len(partners)):
                    atom_i = partners[i_idx]
                    for k_idx in range(i_idx+1, len(partners)):
                        atom_k = partners[k_idx]
                        
                        # Ensure canonical ordering for angle triplets
                        ordered_i = min(atom_i, atom_k)
                        ordered_k = max(atom_i, atom_k)
                        
                        # Check if we've already processed this angle
                        angle_key = (ordered_i, j, ordered_k)
                        if angle_key in angle_set:
                            continue
                            
                        angle_set.add(angle_key)
                        
                        # Calculate current angle with improved numerical stability
                        vec1 = positions_nm[atom_i] - positions_nm[j]
                        vec2 = positions_nm[atom_k] - positions_nm[j]
                        angle = calculate_angle(vec1, vec2)
                        
                        # Use parameterized force constant with ordered indices
                        angle_force.addAngle(ordered_i, j, ordered_k, angle, angle_force_constant)
                        angles.append((ordered_i, j, ordered_k))
        
        system.addForce(angle_force)
        
        # Add NonbondedForce for non-bonded interactions
        nonbonded_force = mm.NonbondedForce()
        
        for i, element in enumerate(elements):
            params = get_element_nonbonded_params(element)
            nonbonded_force.addParticle(params[0], params[1], params[2])
        
        # Track exception pairs to avoid duplicates
        exception_pairs = set()
        
        # Create 1-2 exclusions for bonded atoms
        for i, j, _ in bonds:
            # No need for min/max since bonds are already ordered with i < j
            if (i, j) not in exception_pairs:
                nonbonded_force.addException(i, j, 0.0, 0.1, 0.0)
                exception_pairs.add((i, j))
        
        # Add 1-3 angle exclusions
        for i, j, k in angles:
            # Angles are already in canonical order
            if (i, k) not in exception_pairs:
                nonbonded_force.addException(i, k, 0.0, 0.1, 0.0)
                exception_pairs.add((i, k))
            
        system.addForce(nonbonded_force)
        
        # Create Integrator with parameterized temperature
        integrator = mm.LangevinIntegrator(temperature*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        
        # Create simulation
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        
        # Get initial energy
        state = simulation.context.getState(getEnergy=True)
        start_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Get molecule indices for energy calculations
        molecule1_indices = [i for i, origin in enumerate(atom_origins) if origin == 1]
        molecule2_indices = [i for i, origin in enumerate(atom_origins) if origin == 2]
        
        # Calculate initial interaction energy if both molecules present
        initial_interaction_energy = None
        
        if molecule1_atoms and molecule2_atoms and len(molecule1_indices) > 0 and len(molecule2_indices) > 0:
            try:
                # Create a new copy of the system and context for energy calculations
                system_energy = mm.System()
                for i, element in enumerate(elements):
                    mass = get_element_mass(element)
                    system_energy.addParticle(mass)
                    
                # Add the same forces
                system_energy.addForce(mm.HarmonicBondForce())
                system_energy.addForce(mm.HarmonicAngleForce())
                
                # Add NonbondedForce but with modified interaction groups
                nb_force = mm.NonbondedForce()
                for i, element in enumerate(elements):
                    params = get_element_nonbonded_params(element)
                    nb_force.addParticle(params[0], params[1], params[2])
                    
                # Set up interaction groups
                nb_force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
                nb_force.setCutoffDistance(1.0 * unit.nanometer)
                system_energy.addForce(nb_force)
                
                # Create context for energy calculation
                integrator_energy = mm.VerletIntegrator(0.001 * unit.picoseconds)
                platform = mm.Platform.getPlatformByName('Reference')  # Use CPU platform for consistent results
                context_energy = mm.Context(system_energy, integrator_energy, platform)
                context_energy.setPositions(positions)
                
                # Calculate energy with all interactions
                state_all = context_energy.getState(getEnergy=True)
                energy_all = state_all.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                
                # Calculate energy with only molecule1 interactions
                # Use a separate set to track energy calculation exceptions
                energy_exception_pairs = set()
                
                for i in range(len(all_atoms)):
                    for j in range(i+1, len(all_atoms)):  # Only process each pair once with i < j
                        # Skip intra-molecular interactions (keep them)
                        if (i in molecule1_indices and j in molecule1_indices) or \
                           (i in molecule2_indices and j in molecule2_indices):
                            continue
                        
                        # This is an inter-molecular interaction we want to disable
                        if (i, j) not in energy_exception_pairs:
                            nb_force.addException(i, j, 0.0, 0.1, 0.0)
                            energy_exception_pairs.add((i, j))
                        
                context_energy.reinitialize(True)
                context_energy.setPositions(positions)
                state_separate = context_energy.getState(getEnergy=True)
                energy_separate = state_separate.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                
                # Interaction energy is the difference
                initial_interaction_energy = energy_all - energy_separate
                
                # Clean up energy calculation objects
                del context_energy
                del system_energy
                
            except Exception as energy_error:
                logger.warning(f"Initial interaction energy calculation failed: {str(energy_error)}")
                initial_interaction_energy = None
        
        # Custom iteration tracking
        iterations_performed = 0
        
        try:
            if force_iterations:
                # Manual optimization to force all iterations
                energy_history = []
                
                # Initial minimization to stabilize structure
                simulation.minimizeEnergy(maxIterations=100)
                
                # Execute exact number of iterations requested
                for i in range(max_iterations):
                    simulation.step(1)
                    
                    state = simulation.context.getState(getEnergy=True)
                    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                    energy_history.append(energy)
                    iterations_performed += 1
                    
                    # Safety check for complex molecules
                    if (i > 20 and np.isfinite(energy) and 
                        i < len(energy_history) and 
                        energy > energy_history[i-1] * 5.0) or not np.isfinite(energy):
                        logger.info(f"Breaking at iteration {i} due to energy instability: {energy}")
                        break
                        
                state = simulation.context.getState(getPositions=True, getEnergy=True)
                minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
                final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                
            else:
                # Apply a preliminary short equilibration to stabilize the system first
                # This helps prevent NaN issues with identical molecules during minimization
                if is_similar:
                    logger.info("Performing preliminary equilibration for identical molecules")
                    simulation.step(50)  # 50 very small steps to slightly perturb the structure
                
                # Standard energy minimization with correct force units for tolerance
                force_tolerance = tolerance * unit.kilojoule_per_mole / unit.nanometer
                
                try:
                    # Skip passing tolerance if using default value to avoid unit errors
                    if abs(tolerance - 10.0) < 0.01:  # If close to default value
                        simulation.minimizeEnergy(maxIterations=max_iterations)
                    else:
                        # Use custom tolerance
                        simulation.minimizeEnergy(tolerance=force_tolerance, maxIterations=max_iterations)
                except Exception as min_error:
                    # If minimization fails, try a more robust approach with staged minimization
                    logger.warning(f"Initial minimization failed: {str(min_error)}. Trying staged approach.")
                    
                    # Reset positions with slight perturbation to break exact symmetry
                    perturbed_positions = np.array(positions.value_in_unit(unit.nanometer))
                    # Add small random noise (0.01 nm = 0.1 Å)
                    perturbed_positions += np.random.normal(0, 0.01, perturbed_positions.shape)
                    simulation.context.setPositions(perturbed_positions * unit.nanometer)
                    
                    # Try staged minimization with increasing iterations
                    for stage_iterations in [50, 100, max_iterations]:
                        try:
                            simulation.minimizeEnergy(maxIterations=stage_iterations)
                        except Exception as stage_error:
                            logger.warning(f"Stage minimization with {stage_iterations} iterations failed: {str(stage_error)}")
                            # Continue to next stage anyway
                
                # Get final state
                state = simulation.context.getState(getPositions=True, getEnergy=True)
                final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                energy_change = abs(start_energy - final_energy)
                
                # Estimate iterations based on energy change and tolerance
                iterations_performed = min(max_iterations, int((energy_change/max(tolerance, 1e-6)) * 10))
                
                minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        
        except Exception as opt_error:
            # Final fallback if all optimization attempts fail
            logger.error(f"All optimization attempts failed: {str(opt_error)}")
            
            # Return with original positions
            minimized_positions = np.array(positions.value_in_unit(unit.angstrom))
            final_energy = start_energy
            iterations_performed = 0
            
        # Calculate final interaction energy
        final_interaction_energy = None
        
        if molecule1_atoms and molecule2_atoms and len(molecule1_indices) > 0 and len(molecule2_indices) > 0:
            try:
                # Create context for energy calculation (same as before)
                system_energy = mm.System()
                for i, element in enumerate(elements):
                    mass = get_element_mass(element)
                    system_energy.addParticle(mass)
                    
                # Add forces
                system_energy.addForce(mm.HarmonicBondForce())
                system_energy.addForce(mm.HarmonicAngleForce())
                
                # Add NonbondedForce
                nb_force = mm.NonbondedForce()
                for i, element in enumerate(elements):
                    params = get_element_nonbonded_params(element)
                    nb_force.addParticle(params[0], params[1], params[2])
                    
                # Set up interaction groups
                nb_force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
                nb_force.setCutoffDistance(1.0 * unit.nanometer)
                system_energy.addForce(nb_force)
                
                # Create context for energy calculation
                integrator_energy = mm.VerletIntegrator(0.001 * unit.picoseconds)
                platform = mm.Platform.getPlatformByName('Reference')
                context_energy = mm.Context(system_energy, integrator_energy, platform)
                context_energy.setPositions(minimized_positions * unit.angstrom)
                
                # Calculate energy with all interactions
                state_all = context_energy.getState(getEnergy=True)
                energy_all = state_all.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                
                # Calculate energy with only intramolecular interactions
                energy_exception_pairs = set()
                
                for i in range(len(all_atoms)):
                    for j in range(i+1, len(all_atoms)):
                        # Skip intra-molecular interactions (keep them)
                        if (i in molecule1_indices and j in molecule1_indices) or \
                           (i in molecule2_indices and j in molecule2_indices):
                            continue
                        
                        # This is an inter-molecular interaction we want to disable
                        if (i, j) not in energy_exception_pairs:
                            nb_force.addException(i, j, 0.0, 0.1, 0.0)
                            energy_exception_pairs.add((i, j))
                        
                context_energy.reinitialize(True)
                context_energy.setPositions(minimized_positions * unit.angstrom)
                state_separate = context_energy.getState(getEnergy=True)
                energy_separate = state_separate.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                
                # Interaction energy is the difference
                final_interaction_energy = energy_all - energy_separate
                
                # Clean up energy calculation objects
                del context_energy
                del system_energy
                
            except Exception as energy_error:
                logger.warning(f"Final interaction energy calculation failed: {str(energy_error)}")
                final_interaction_energy = None
        
        # Format result with optimized coordinates, separated by molecule
        molecule1_optimized_atoms = []
        molecule2_optimized_atoms = []
        
        for i, atom in enumerate(all_atoms):
            optimized_atom = {
                "id": atom.get("id", i+1),  # Default to index+1 if no ID
                "element": atom["element"],
                "x": float(minimized_positions[i][0]),
                "y": float(minimized_positions[i][1]),
                "z": float(minimized_positions[i][2])
            }
            
            # Use safe indexing to prevent index errors
            if i < len(atom_origins):
                if atom_origins[i] == 1:
                    molecule1_optimized_atoms.append(optimized_atom)
                else:
                    molecule2_optimized_atoms.append(optimized_atom)
            
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        # Ensure we clean up the simulation object to release GPU resources
        del simulation
            
        return {
            "molecule1_optimized_atoms": molecule1_optimized_atoms,
            "molecule2_optimized_atoms": molecule2_optimized_atoms,
            "metadata": {
                "method": "classical_molecular_dynamics",
                "library": "OpenMM",
                "parameters": {
                    "temperature": temperature,
                    "max_iterations": max_iterations,
                    "bond_threshold": bond_threshold,
                    "bond_force_constant": bond_force_constant,
                    "angle_force_constant": angle_force_constant,
                    "tolerance": tolerance,
                    "force_iterations": force_iterations
                },
                "molecules": 2 if molecule1_atoms and molecule2_atoms else 1,
                "molecule1_atom_count": len(molecule1_atoms) if molecule1_atoms else 0,
                "molecule2_atom_count": len(molecule2_atoms) if molecule2_atoms else 0,
                "iterations_performed": iterations_performed,
                "final_energy_kj_mol": float(final_energy),
                "initial_energy_kj_mol": float(start_energy),
                "energy_change_kj_mol": float(abs(start_energy - final_energy)),
                "duration_seconds": duration,
                "convergence": "energy_minimized",
                "identical_molecules_detected": bool(is_similar),  # Convert NumPy boolean to Python boolean
                "molecule_rmsd": float(rmsd) if rmsd is not None else None,
                "bonds_detected": len(bonds),
                "intermolecular_bonds": intermolecular_bonds,
                "angles_detected": len(angles),
                "initial_interaction_energy_kj_mol": float(initial_interaction_energy) if initial_interaction_energy is not None else None,
                "interaction_energy_kj_mol": float(final_interaction_energy) if final_interaction_energy is not None else None
            }
        }
        
    except Exception as e:
        logger.error(f"Classical combined optimization error: {str(e)}", exc_info=True)
        return {
            "error": str(e),
            "metadata": {
                "method": "classical_molecular_dynamics",
                "library": "OpenMM",
                "status": "failed"
            }
        }      


def optimize_quantum(atoms, params=None):
    """
    Optimize molecule using quantum chemistry with PySCF.
    
    Args:
        atoms: List of atom dictionaries with element and coordinates
        params: Dictionary of optimization parameters
        
    Returns:
        Dictionary containing optimized atoms and metadata
    """
    try:       
        start_time = datetime.now()
        
        # Set default parameters if not provided
        if params is None:
            params = {}
        # Extract parameters with defaults
        basis = params.get("basis", "6-31g")
        max_iterations = params.get("max_iterations", 10)
        convergence_threshold = params.get("convergence_threshold", 1e-5)
        step_size = params.get("step_size", 0.1)
        
        # New memory optimization parameters
        use_direct_scf = params.get("direct_scf", True)  # Enable direct SCF by default
        use_density_fitting = params.get("density_fitting", False)  # Optional density fitting
        
        # Check system size against thresholds
        atom_count = len(atoms)
        threshold = QUANTUM_ATOM_THRESHOLDS.get(basis, 30)  # Default to 30 if basis not in dict
        
        if atom_count > threshold:
            return {
                "error": f"System size of {atom_count} atoms exceeds the recommended limit of {threshold} for basis {basis}. Consider using classical optimization instead.",
                "metadata": {
                    "method": "quantum_chemistry",
                    "library": "PySCF",
                    "status": "failed",
                    "reason": "system_too_large"
                }
            }
        
        # Convert atoms to PySCF format
        atom_list = []
        for atom in atoms:
            atom_list.append([atom["element"], (atom["x"], atom["y"], atom["z"])])
        
        # Create molecule with parameterized basis set
        mol = gto.M(
            atom=atom_list,
            basis=basis,
            verbose=0
        )
        
        # Run Hartree-Fock calculation with memory optimizations
        mf = scf.RHF(mol)
        
        # Enable direct SCF to avoid storing integrals in memory
        if use_direct_scf:
            mf.direct_scf = True
            logger.info("Using direct SCF for quantum optimization")
            
        # Apply density fitting if requested (for even larger systems)
        if use_density_fitting:
            mf = scf.density_fit(mf)
            logger.info("Using density fitting approximation for quantum optimization")
        
        energy = mf.kernel()
        
        # Check if SCF converged
        if not mf.converged:
            logger.warning("Initial SCF calculation did not converge")
        
        # Calculate gradient
        g = grad.RHF(mf)
        gradients = g.kernel()
        
        # Perform geometry optimization with parameterized values
        converged = False
        iteration = 0
        current_atoms = atom_list.copy()
        current_energy = energy
        optimization_history = []
        
        while not converged and iteration < max_iterations:
            # Store current state in history
            optimization_history.append({
                'iteration': iteration,
                'energy': float(current_energy),
                'gradient_norm': float(np.linalg.norm(gradients))
            })
            
            # Get coordinates
            coords = np.array([atom[1] for atom in current_atoms])
            
            # Update geometry based on gradient using parameterized step size
            new_coords = coords - step_size * gradients
            
            # Update atom list
            new_atom_list = []
            for i, atom in enumerate(current_atoms):
                new_atom_list.append([atom[0], tuple(new_coords[i])])
            
            # Create new molecule
            new_mol = gto.M(
                atom=new_atom_list,
                basis=basis,  # Use parameterized basis
                verbose=0
            )
            
            # Run Hartree-Fock on new geometry with same memory optimizations
            new_mf = scf.RHF(new_mol)
            
            # Apply same memory optimization settings
            if use_direct_scf:
                new_mf.direct_scf = True
                
            if use_density_fitting:
                new_mf = scf.density_fit(new_mf)
                
            new_energy = new_mf.kernel()
            
            # Check if SCF converged
            if not new_mf.converged:
                logger.warning(f"SCF calculation at iteration {iteration+1} did not converge")
            
            # Calculate new gradient
            new_g = grad.RHF(new_mf)
            new_gradients = new_g.kernel()
            
            # Check convergence using parameterized threshold
            energy_change = abs(new_energy - current_energy)
            gradient_norm = np.linalg.norm(new_gradients)
            
            # Use a more robust convergence criterion
            if (energy_change < convergence_threshold and gradient_norm < convergence_threshold):
                converged = True
                logger.info(f"Geometry optimization converged at iteration {iteration+1}")
            
            # Energy increase check - break if energy increases significantly
            if new_energy > current_energy + 0.1:  # 0.1 Hartree is a significant increase
                logger.warning(f"Breaking optimization at iteration {iteration+1} due to energy increase")
                # Don't update to the higher energy state
                break
            
            # Update for next iteration
            current_atoms = new_atom_list
            current_energy = new_energy
            gradients = new_gradients
            iteration += 1
        
        # Create result with optimized atoms
        optimized_atoms = []
        for i, atom in enumerate(atoms):
            optimized_atoms.append({
                "id": atom.get("id", i+1),
                "element": atom["element"],
                "x": float(current_atoms[i][1][0]),
                "y": float(current_atoms[i][1][1]),
                "z": float(current_atoms[i][1][2])
            })
            
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        return {
            "optimized_atoms": optimized_atoms,
            "metadata": {
                "method": "quantum_chemistry",
                "library": "PySCF",
                "parameters": {
                    "basis": basis,
                    "max_iterations": max_iterations,
                    "convergence_threshold": convergence_threshold,
                    "step_size": step_size,
                    "direct_scf": use_direct_scf,
                    "density_fitting": use_density_fitting
                },
                "theory_level": f"RHF/{basis}",
                "final_energy_hartree": float(current_energy),
                "iterations": iteration,
                "converged": bool(converged),  # Convert NumPy boolean to Python boolean
                "duration_seconds": duration,
                "optimization_history": optimization_history
            }
        }
        
    except Exception as e:
        logger.error(f"Quantum optimization error: {str(e)}", exc_info=True)
        return {
            "error": str(e),
            "metadata": {
                "method": "quantum_chemistry",
                "library": "PySCF",
                "status": "failed"
            }
        }
    

def optimize_quantum_combined(molecule1_atoms, molecule2_atoms, params=None):
    """
    Optimize combined molecule system using quantum chemistry with PySCF.
    
    Args:
        molecule1_atoms: Atoms from the first molecule (can be None)
        molecule2_atoms: Atoms from the second molecule (can be None)
        params: Optimization parameters
        
    Returns:
        Dictionary with optimized atoms and metadata
    """
    try:
        start_time = datetime.now()
        
        # Set default parameters if not provided
        if params is None:
            params = {}
            
        # Extract parameters with defaults
        basis = params.get("basis", "6-31g")
        max_iterations = params.get("max_iterations", 10)
        convergence_threshold = params.get("convergence_threshold", 1e-5)
        step_size = params.get("step_size", 0.1)
        
        # New memory optimization parameters
        use_direct_scf = params.get("direct_scf", True)  # Enable direct SCF by default
        use_density_fitting = params.get("density_fitting", False)  # Optional density fitting
        
        # Check if molecules are too similar (causes numerical instability)
        is_similar = False
        rmsd = None
        
        if molecule1_atoms and molecule2_atoms:
            is_similar, rmsd = are_molecules_similar(molecule1_atoms, molecule2_atoms)
            
            if is_similar:
                logger.info(f"Detected identical or very similar molecules (RMSD={rmsd}Å), applying displacement")
                # Apply a 3Å displacement to second molecule to avoid instability
                molecule2_atoms = displace_molecule(molecule2_atoms, magnitude=3.0)
        
        # Combine atoms from both molecules, tracking their origins
        all_atoms = []
        atom_origins = []  # Track which molecule each atom came from (1 or 2)
        
        # Add atoms from molecule 1 if present
        if molecule1_atoms is not None:
            for atom in molecule1_atoms:
                all_atoms.append(atom)
                atom_origins.append(1)
                
        # Add atoms from molecule 2 if present
        if molecule2_atoms is not None:
            for atom in molecule2_atoms:
                all_atoms.append(atom)
                atom_origins.append(2)
        
        # Check system size against thresholds
        total_atoms = len(all_atoms)
        threshold = QUANTUM_ATOM_THRESHOLDS.get(basis, 30)  # Default to 30 if basis not in dict
        
        if total_atoms > threshold:
            return {
                "error": f"System size of {total_atoms} atoms exceeds the recommended limit of {threshold} for basis {basis}. Consider using classical optimization instead.",
                "metadata": {
                    "method": "quantum_chemistry",
                    "library": "PySCF",
                    "status": "failed",
                    "reason": "system_too_large"
                }
            }
        
        # Get indices for atoms from each molecule
        molecule1_indices = [i for i, origin in enumerate(atom_origins) if origin == 1]
        molecule2_indices = [i for i, origin in enumerate(atom_origins) if origin == 2]
        
        # Convert atoms to PySCF format for all atoms
        all_atom_list = []
        for atom in all_atoms:
            all_atom_list.append([atom["element"], (atom["x"], atom["y"], atom["z"])])
        
        # Separate molecules for interaction energy calculation
        molecule1_atom_list = []
        if molecule1_atoms:
            for atom in molecule1_atoms:
                molecule1_atom_list.append([atom["element"], (atom["x"], atom["y"], atom["z"])])
                
        molecule2_atom_list = []
        if molecule2_atoms:
            for atom in molecule2_atoms:
                molecule2_atom_list.append([atom["element"], (atom["x"], atom["y"], atom["z"])])
        
        # Calculate initial interaction energy if both molecules present
        initial_interaction_energy = None
        initial_energy1 = None
        initial_energy2 = None
        
        if molecule1_atoms and molecule2_atoms and len(molecule1_indices) > 0 and len(molecule2_indices) > 0:
            try:
                # Create molecule for combined system
                mol_combined = gto.M(
                    atom=all_atom_list,
                    basis=basis,
                    verbose=0
                )
                
                # Run Hartree-Fock on combined system with memory optimizations
                mf_combined = scf.RHF(mol_combined)
                
                # Apply memory optimizations
                if use_direct_scf:
                    mf_combined.direct_scf = True
                    
                if use_density_fitting:
                    mf_combined = scf.density_fit(mf_combined)
                    
                energy_combined = mf_combined.kernel()
                
                # Create and calculate energy for molecule 1 with memory optimizations
                mol1 = gto.M(
                    atom=molecule1_atom_list,
                    basis=basis,
                    verbose=0
                )
                mf1 = scf.RHF(mol1)
                
                # Apply same memory optimizations
                if use_direct_scf:
                    mf1.direct_scf = True
                    
                if use_density_fitting:
                    mf1 = scf.density_fit(mf1)
                    
                energy1 = mf1.kernel()
                initial_energy1 = energy1
                
                # Create and calculate energy for molecule 2 with memory optimizations
                mol2 = gto.M(
                    atom=molecule2_atom_list,
                    basis=basis,
                    verbose=0
                )
                mf2 = scf.RHF(mol2)
                
                # Apply same memory optimizations
                if use_direct_scf:
                    mf2.direct_scf = True
                    
                if use_density_fitting:
                    mf2 = scf.density_fit(mf2)
                    
                energy2 = mf2.kernel()
                initial_energy2 = energy2
                
                # Interaction energy is the difference
                initial_interaction_energy = energy_combined - (energy1 + energy2)
            except Exception as energy_error:
                logger.warning(f"Initial interaction energy calculation failed: {str(energy_error)}")
                initial_interaction_energy = None
        
        # Create molecule with complete system
        mol = gto.M(
            atom=all_atom_list,
            basis=basis,
            verbose=0
        )
        
        # Run Hartree-Fock calculation with memory optimizations
        mf = scf.RHF(mol)
        
        # Apply memory optimizations
        if use_direct_scf:
            mf.direct_scf = True
            logger.info("Using direct SCF for quantum combined optimization")
            
        if use_density_fitting:
            mf = scf.density_fit(mf)
            logger.info("Using density fitting approximation for quantum combined optimization")
            
        energy = mf.kernel()
        
        # Check if SCF converged
        if not mf.converged:
            logger.warning("Initial SCF calculation did not converge")
        
        # Calculate gradient
        g = grad.RHF(mf)
        gradients = g.kernel()
        
        # Perform geometry optimization with parameterized values
        converged = False
        iteration = 0
        current_atoms = all_atom_list.copy()
        current_energy = energy
        optimization_history = []
        
        while not converged and iteration < max_iterations:
            # Store current state in history
            optimization_history.append({
                'iteration': iteration,
                'energy': float(current_energy),
                'gradient_norm': float(np.linalg.norm(gradients))
            })
            
            # Get coordinates
            coords = np.array([atom[1] for atom in current_atoms])
            
            # Update geometry based on gradient using parameterized step size
            new_coords = coords - step_size * gradients
            
            # Update atom list
            new_atom_list = []
            for i, atom in enumerate(current_atoms):
                new_atom_list.append([atom[0], tuple(new_coords[i])])
            
            # Create new molecule
            new_mol = gto.M(
                atom=new_atom_list,
                basis=basis,
                verbose=0
            )
            
            # Run Hartree-Fock on new geometry with same memory optimizations
            new_mf = scf.RHF(new_mol)
            
            # Apply same memory optimization settings
            if use_direct_scf:
                new_mf.direct_scf = True
                
            if use_density_fitting:
                new_mf = scf.density_fit(new_mf)
                
            new_energy = new_mf.kernel()
            
            # Check if SCF converged
            if not new_mf.converged:
                logger.warning(f"SCF calculation at iteration {iteration+1} did not converge")
            
            # Calculate new gradient
            new_g = grad.RHF(new_mf)
            new_gradients = new_g.kernel()
            
            # Check convergence using parameterized threshold
            energy_change = abs(new_energy - current_energy)
            gradient_norm = np.linalg.norm(new_gradients)
            
            # Use a more robust convergence criterion
            if (energy_change < convergence_threshold and gradient_norm < convergence_threshold):
                converged = True
                logger.info(f"Geometry optimization converged at iteration {iteration+1}")
            
            # Energy increase check - break if energy increases significantly
            if new_energy > current_energy + 0.1:  # 0.1 Hartree is a significant increase
                logger.warning(f"Breaking optimization at iteration {iteration+1} due to energy increase")
                # Don't update to the higher energy state
                break
            
            # Update for next iteration
            current_atoms = new_atom_list
            current_energy = new_energy
            gradients = new_gradients
            iteration += 1
        
        # Calculate final interaction energy
        final_interaction_energy = None
        final_energy1 = None
        final_energy2 = None
        
        if molecule1_atoms and molecule2_atoms and len(molecule1_indices) > 0 and len(molecule2_indices) > 0:
            try:
                # Extract optimized coordinates for each molecule
                optimized_coords = np.array([atom[1] for atom in current_atoms])
                
                # Create atom lists for optimized individual molecules
                mol1_opt_atoms = []
                for i in molecule1_indices:
                    if i < len(current_atoms):
                        mol1_opt_atoms.append([current_atoms[i][0], current_atoms[i][1]])
                    
                mol2_opt_atoms = []
                for i in molecule2_indices:
                    if i < len(current_atoms):
                        mol2_opt_atoms.append([current_atoms[i][0], current_atoms[i][1]])
                
                # Calculate energy for optimized combined system with memory optimizations
                mol_combined_opt = gto.M(
                    atom=current_atoms,
                    basis=basis,
                    verbose=0
                )
                mf_combined_opt = scf.RHF(mol_combined_opt)
                
                # Apply memory optimizations
                if use_direct_scf:
                    mf_combined_opt.direct_scf = True
                    
                if use_density_fitting:
                    mf_combined_opt = scf.density_fit(mf_combined_opt)
                    
                energy_combined_opt = mf_combined_opt.kernel()
                
                # Calculate energy for optimized molecule 1 with memory optimizations
                if mol1_opt_atoms:
                    mol1_opt = gto.M(
                        atom=mol1_opt_atoms,
                        basis=basis,
                        verbose=0
                    )
                    mf1_opt = scf.RHF(mol1_opt)
                    
                    # Apply memory optimizations
                    if use_direct_scf:
                        mf1_opt.direct_scf = True
                        
                    if use_density_fitting:
                        mf1_opt = scf.density_fit(mf1_opt)
                        
                    energy1_opt = mf1_opt.kernel()
                    final_energy1 = energy1_opt
                else:
                    energy1_opt = 0
                
                # Calculate energy for optimized molecule 2 with memory optimizations
                if mol2_opt_atoms:
                    mol2_opt = gto.M(
                        atom=mol2_opt_atoms,
                        basis=basis,
                        verbose=0
                    )
                    mf2_opt = scf.RHF(mol2_opt)
                    
                    # Apply memory optimizations
                    if use_direct_scf:
                        mf2_opt.direct_scf = True
                        
                    if use_density_fitting:
                        mf2_opt = scf.density_fit(mf2_opt)
                        
                    energy2_opt = mf2_opt.kernel()
                    final_energy2 = energy2_opt
                else:
                    energy2_opt = 0
                
                # Interaction energy is the difference
                final_interaction_energy = energy_combined_opt - (energy1_opt + energy2_opt)
            except Exception as energy_error:
                logger.warning(f"Final interaction energy calculation failed: {str(energy_error)}")
                final_interaction_energy = None
        
        # Create result with optimized atoms
        molecule1_optimized_atoms = []
        molecule2_optimized_atoms = []
        
        for i, atom in enumerate(all_atoms):
            # Index safety check
            if i < len(current_atoms):
                optimized_atom = {
                    "id": atom.get("id", i+1),
                    "element": atom["element"],
                    "x": float(current_atoms[i][1][0]),
                    "y": float(current_atoms[i][1][1]),
                    "z": float(current_atoms[i][1][2])
                }
                
                # Use safe indexing to prevent index errors
                if i < len(atom_origins):
                    if atom_origins[i] == 1:
                        molecule1_optimized_atoms.append(optimized_atom)
                    else:
                        molecule2_optimized_atoms.append(optimized_atom)
            
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        return {
            "molecule1_optimized_atoms": molecule1_optimized_atoms,
            "molecule2_optimized_atoms": molecule2_optimized_atoms,
            "metadata": {
                "method": "quantum_chemistry",
                "library": "PySCF",
                "parameters": {
                    "basis": basis,
                    "max_iterations": max_iterations,
                    "convergence_threshold": convergence_threshold,
                    "step_size": step_size,
                    "direct_scf": use_direct_scf,
                    "density_fitting": use_density_fitting
                },
                "molecules": 2 if molecule1_atoms and molecule2_atoms else 1,
                "molecule1_atom_count": len(molecule1_atoms) if molecule1_atoms else 0,
                "molecule2_atom_count": len(molecule2_atoms) if molecule2_atoms else 0,
                "identical_molecules_detected": bool(is_similar),  # Convert NumPy boolean to Python boolean
                "molecule_rmsd": float(rmsd) if rmsd is not None else None,
                "theory_level": f"RHF/{basis}",
                "final_energy_hartree": float(current_energy),
                "initial_molecule1_energy_hartree": float(initial_energy1) if initial_energy1 is not None else None,
                "initial_molecule2_energy_hartree": float(initial_energy2) if initial_energy2 is not None else None,
                "final_molecule1_energy_hartree": float(final_energy1) if final_energy1 is not None else None,
                "final_molecule2_energy_hartree": float(final_energy2) if final_energy2 is not None else None,
                "iterations": iteration,
                "converged": bool(converged),  # Convert NumPy boolean to Python boolean
                "duration_seconds": duration,
                "initial_interaction_energy_hartree": float(initial_interaction_energy) if initial_interaction_energy is not None else None,
                "interaction_energy_hartree": float(final_interaction_energy) if final_interaction_energy is not None else None,
                "optimization_history": optimization_history
            }
        }
        
    except Exception as e:
        logger.error(f"Quantum combined optimization error: {str(e)}", exc_info=True)
        return {
            "error": str(e),
            "metadata": {
                "method": "quantum_chemistry",
                "library": "PySCF",
                "status": "failed"
            }
        }


def extract_molecule_atoms(molecule_data):
    """
    Extract atoms from molecule data with robust error handling.
    
    Args:
        molecule_data: Dictionary containing molecule data
        
    Returns:
        List of atom dictionaries or None if not found/invalid
    """
    if not molecule_data:
        return None
        
    try:
        # Try different possible structures
        if "file1" in molecule_data and "atoms" in molecule_data["file1"]:
            return molecule_data["file1"]["atoms"]
        elif "atoms" in molecule_data:
            return molecule_data["atoms"]
        elif isinstance(molecule_data, list):
            # Direct atom list
            return molecule_data
        else:
            # Try to iterate through keys to find atoms
            for key in molecule_data:
                if isinstance(molecule_data[key], dict) and "atoms" in molecule_data[key]:
                    return molecule_data[key]["atoms"]
        
        # No atoms found
        return None
    except Exception as e:
        logger.error(f"Error extracting molecule atoms: {str(e)}")
        return None


@opti_bp.route('/optimize-molecule', methods=['POST'])
def optimize_molecule():
    """
    Endpoint to optimize molecular structures using either classical or quantum methods.
    Now supports both single molecule optimization and interaction between two molecules.
    
    Expected request format:
    {
        "optimization_type": "classical" or "quantum",
        "optimization_params": {
            // Parameters for the selected optimization type
        },
        "molecule1": {
            "file1": {
                "atoms": [
                    {"id": 1, "element": "C", "x": 1.1, "y": 0.0, "z": 0.0},
                    ...
                ]
            }
        },
        "molecule2": {
            "file1": {
                "atoms": [
                    {"id": 1, "element": "C", "x": 1.1, "y": 0.0, "z": 0.0},
                    ...
                ]
            }
        },
        "interaction_mode": true/false  // Whether to optimize molecular interaction
    }
    """
    try:
        data = request.get_json()
        if not data:
            return jsonify({"error": "Invalid input: No data provided"}), 400
            
        molecule1_data = data.get("molecule1")
        molecule2_data = data.get("molecule2")
        optimization_type = data.get("optimization_type")
        optimization_params = data.get("optimization_params", {})
        interaction_mode = data.get("interaction_mode", False)
        
        # Extract atoms with robust parsing
        molecule1_atoms = extract_molecule_atoms(molecule1_data)
        molecule2_atoms = extract_molecule_atoms(molecule2_data)
            
        # Validate inputs
        if not molecule1_atoms and not molecule2_atoms:
            return jsonify({"error": "Invalid input: At least one molecule with atoms must be provided"}), 400
            
        if interaction_mode and (not molecule1_atoms or not molecule2_atoms):
            return jsonify({"error": "Invalid input: Both molecules with atoms must be provided for interaction mode"}), 400
            
        if not optimization_type or optimization_type not in ["classical", "quantum"]:
            return jsonify({"error": "Invalid input: Missing or invalid 'optimization_type' (must be 'classical' or 'quantum')"}), 400
        
        # Determine user subscription status - default to unsubscribed
        is_subscribed = False
        user_id = None
        
        # Try to verify JWT token without raising an exception if missing
        try:
            verify_jwt_in_request(optional=True)
            # If we have a valid token, get the user
            user_identity = get_jwt_identity()
            if user_identity:
                try:
                    # Convert string user_id back to integer for database query with proper error handling
                    user_id_int = int(user_identity) if isinstance(user_identity, str) else user_identity
                    user = User.query.get(user_id_int)
                    if user and user.subscription_status == "active":
                        is_subscribed = True
                        user_id = user.id
                except (ValueError, TypeError) as e:
                    logger.warning(f"Invalid user ID format in JWT: {str(e)}")
        except Exception as e:
            # If token verification fails, proceed as unauthenticated
            logger.warning(f"JWT verification error (proceeding as unauthenticated): {str(e)}")
        
        # Apply subscription-based limits to parameters
        if optimization_type == "classical":
            max_iterations_limit = ITERATION_LIMITS['subscribed' if is_subscribed else 'unsubscribed']['classical']
            # Cap max_iterations to the appropriate limit
            optimization_params['max_iterations'] = min(
                optimization_params.get('max_iterations', 1000), 
                max_iterations_limit
            )
        
        elif optimization_type == "quantum":
            max_iterations_limit = ITERATION_LIMITS['subscribed' if is_subscribed else 'unsubscribed']['quantum']
            # Cap max_iterations to the appropriate limit
            optimization_params['max_iterations'] = min(
                optimization_params.get('max_iterations', 10), 
                max_iterations_limit
            )
            
            # Restrict advanced basis sets for non-subscribers
            if not is_subscribed and optimization_params.get('basis') in ["6-311g", "cc-pvdz"]:
                optimization_params['basis'] = "6-31g"
                logger.info("Downgraded basis set to 6-31g for non-subscriber")
                
            # Force direct SCF for quantum calculations to prevent memory issues
            optimization_params['direct_scf'] = True
            
            # Enable density fitting for larger systems (over 25 atoms) or when requested
            total_atoms = (len(molecule1_atoms) if molecule1_atoms else 0) + \
                         (len(molecule2_atoms) if molecule2_atoms else 0)
            if total_atoms > 25 or optimization_params.get('density_fitting', False):
                optimization_params['density_fitting'] = True
                logger.info(f"Enabling density fitting for system with {total_atoms} atoms")
        
        # Log the optimization request
        logger.info(f"Starting {optimization_type} optimization, interaction_mode={interaction_mode}, "
                  f"molecule1_atoms={len(molecule1_atoms) if molecule1_atoms else 0}, "
                  f"molecule2_atoms={len(molecule2_atoms) if molecule2_atoms else 0}")
        
        # Check for similar/identical molecules if in interaction mode
        if interaction_mode and molecule1_atoms and molecule2_atoms:
            is_similar, rmsd = are_molecules_similar(molecule1_atoms, molecule2_atoms)
            if is_similar:
                logger.info(f"Detected identical or very similar molecules (RMSD={rmsd}Å) in interaction mode")
        
        # Perform requested optimization with parameters
        result = None
        if interaction_mode:
            # For molecular interaction, use the combined optimization functions
            if optimization_type == "classical":
                result = optimize_classical_combined(molecule1_atoms, molecule2_atoms, optimization_params)
            elif optimization_type == "quantum":
                result = optimize_quantum_combined(molecule1_atoms, molecule2_atoms, optimization_params)
        else:
            # For single molecule optimization, use the original functions
            # Use molecule1 if available, otherwise use molecule2
            atoms = molecule1_atoms if molecule1_atoms else molecule2_atoms
            if optimization_type == "classical":
                result = optimize_classical(atoms, optimization_params)
            elif optimization_type == "quantum":
                result = optimize_quantum(atoms, optimization_params)
        
        # Check if optimization was successful
        if result and "error" in result:
            logger.error(f"Optimization failed: {result['error']}")
        
        # Store optimization results in database
        try:
            optimization = Optimization(
                user_id=user_id,
                optimization_type=optimization_type,
                parameters=json.dumps({
                    "optimization_params": optimization_params,
                    "interaction_mode": interaction_mode,
                    "molecule1_atom_count": len(molecule1_atoms) if molecule1_atoms else 0,
                    "molecule2_atom_count": len(molecule2_atoms) if molecule2_atoms else 0
                }),
                result=json.dumps(result)
            )
            db.session.add(optimization)
            db.session.commit()
            optimization_id = optimization.id
        except Exception as db_error:
            logger.error(f"Database error: {str(db_error)}")
            # Continue without database storage
            optimization_id = None
        
        # Return results
        return jsonify({
            "success": True if not (result and "error" in result) else False,
            "optimizationId": optimization_id,
            "result": result
        }), 200
        
    except Exception as e:
        logger.error(f"Unexpected error in optimize_molecule endpoint: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500