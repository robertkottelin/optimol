from flask import Blueprint, request, jsonify
import json
import numpy as np
from datetime import datetime
from extensions import db
from user import User

# Classical optimization imports
import openmm as mm
import openmm.app as app
from openmm import unit

# Quantum optimization imports
from pyscf import gto, scf, grad

from flask_jwt_extended import jwt_required, get_jwt_identity, verify_jwt_in_request
from constants import ITERATION_LIMITS

# Create blueprint
opti_bp = Blueprint('opti', __name__)

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
        
        # Element -> Mass mapping (g/mol)
        element_masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
            'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904
        }
        
        # Add atoms with proper masses
        atom_objects = []
        for i, element in enumerate(elements):
            mass = element_masses.get(element, 0.0)
            atom_obj = topology.addAtom(element, app.Element.getBySymbol(element), residue)
            atom_objects.append(atom_obj)
        
        # Create system
        system = mm.System()
        
        # Add atom masses to the system
        for i, element in enumerate(elements):
            mass = element_masses.get(element, 0.0)
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
                        # Calculate current angle
                        vec1 = positions_nm[atom_i] - positions_nm[j]
                        vec2 = positions_nm[atom_k] - positions_nm[j]
                        norm1 = np.sqrt(np.sum(vec1*vec1))
                        norm2 = np.sqrt(np.sum(vec2*vec2))
                        dot = np.sum(vec1*vec2)
                        cosine = dot / (norm1 * norm2)
                        cosine = max(-1.0, min(1.0, cosine))
                        angle = np.arccos(cosine)
                        # Use parameterized force constant
                        angle_force.addAngle(atom_i, j, atom_k, angle, angle_force_constant)
                        angles.append((atom_i, j, atom_k))
        
        system.addForce(angle_force)
        
        # Add NonbondedForce for non-bonded interactions
        nonbonded_force = mm.NonbondedForce()
        
        # Element -> [charge, sigma (nm), epsilon (kJ/mol)]
        element_params = {
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
        
        for i, element in enumerate(elements):
            params = element_params.get(element, [0.0, 0.3, 0.5])
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
                
                # MODIFIED: More robust safety check for complex molecules
                # Only break if energy increases dramatically (5x) and we're past initial fluctuations
                if (i > 20 and energy > energy_history[i-1] * 5.0) or not np.isfinite(energy):
                    print(f"Breaking at iteration {i} due to energy instability: {energy}")
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
            iterations_performed = min(max_iterations, int((energy_change/tolerance) * 10))
            
            minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        
        # Format result with optimized coordinates
        optimized_atoms = []
        for i, atom in enumerate(atoms):
            optimized_atoms.append({
                "id": atom["id"],
                "element": atom["element"],
                "x": float(minimized_positions[i][0]),
                "y": float(minimized_positions[i][1]),
                "z": float(minimized_positions[i][2])
            })
            
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
            
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
        
        # Element -> Mass mapping (g/mol)
        element_masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
            'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904
        }
        
        # Add atoms with proper masses
        atom_objects = []
        for i, element in enumerate(elements):
            mass = element_masses.get(element, 0.0)
            atom_obj = topology.addAtom(element, app.Element.getBySymbol(element), residue)
            atom_objects.append(atom_obj)
        
        # Create system
        system = mm.System()
        
        # Add atom masses to the system
        for i, element in enumerate(elements):
            mass = element_masses.get(element, 0.0)
            system.addParticle(mass)
        
        # Determine bonds based on distance
        positions_nm = positions.value_in_unit(unit.nanometer)
        bonds = []
        intermolecular_bonds = 0  # Counter for bonds between molecules
        
        for i in range(len(atom_objects)):
            for j in range(i+1, len(atom_objects)):
                dist = np.sqrt(np.sum((positions_nm[i] - positions_nm[j])**2))
                if dist < bond_threshold:  # Bond threshold in nm
                    topology.addBond(atom_objects[i], atom_objects[j])
                    bonds.append((i, j, dist))
                    
                    # Check if this bond connects atoms from different molecules
                    if atom_origins[i] != atom_origins[j]:
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
        for j in range(len(atom_objects)):
            if j in bond_partners:
                partners = bond_partners[j]
                for i in range(len(partners)):
                    for k in range(i+1, len(partners)):
                        atom_i = partners[i]
                        atom_k = partners[k]
                        # Calculate current angle
                        vec1 = positions_nm[atom_i] - positions_nm[j]
                        vec2 = positions_nm[atom_k] - positions_nm[j]
                        norm1 = np.sqrt(np.sum(vec1*vec1))
                        norm2 = np.sqrt(np.sum(vec2*vec2))
                        dot = np.sum(vec1*vec2)
                        cosine = dot / (norm1 * norm2)
                        cosine = max(-1.0, min(1.0, cosine))
                        angle = np.arccos(cosine)
                        # Use parameterized force constant
                        angle_force.addAngle(atom_i, j, atom_k, angle, angle_force_constant)
                        angles.append((atom_i, j, atom_k))
        
        system.addForce(angle_force)
        
        # Add NonbondedForce for non-bonded interactions
        nonbonded_force = mm.NonbondedForce()
        
        # Element -> [charge, sigma (nm), epsilon (kJ/mol)]
        element_params = {
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
        
        for i, element in enumerate(elements):
            params = element_params.get(element, [0.0, 0.3, 0.5])
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
        
        # Calculate initial energy for each molecule separately to get interaction energy
        molecule1_indices = [i for i, origin in enumerate(atom_origins) if origin == 1]
        molecule2_indices = [i for i, origin in enumerate(atom_origins) if origin == 2]
        
        # Create interactions groups for molecule1 and molecule2
        initial_interaction_energy = None
        
        if molecule1_atoms and molecule2_atoms:
            # Create a new copy of the system and context for energy calculations
            system_energy = mm.System()
            for i, element in enumerate(elements):
                mass = element_masses.get(element, 0.0)
                system_energy.addParticle(mass)
                
            # Add the same forces
            system_energy.addForce(mm.HarmonicBondForce())
            system_energy.addForce(mm.HarmonicAngleForce())
            
            # Add NonbondedForce but with modified interaction groups
            nb_force = mm.NonbondedForce()
            for i, element in enumerate(elements):
                params = element_params.get(element, [0.0, 0.3, 0.5])
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
            energy_exception_pairs = set()

            # Calculate energy with only molecule1 interactions
            for i in range(len(all_atoms)):
                for j in range(i+1, len(all_atoms)):  # Only process each pair once with i < j
                    # Skip intra-molecular interactions (keep them)
                    if (i in molecule1_indices and j in molecule1_indices) or \
                    (i in molecule2_indices and j in molecule2_indices):
                        continue
                        
                    # This is an inter-molecular interaction we want to disable
                    pair = (i, j)  # No need for min/max here since i < j is enforced
                    if pair not in energy_exception_pairs:
                        nb_force.addException(i, j, 0.0, 0.1, 0.0)
                        energy_exception_pairs.add(pair)
                    
            context_energy.reinitialize(True)
            context_energy.setPositions(positions)
            state_separate = context_energy.getState(getEnergy=True)
            energy_separate = state_separate.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            # Interaction energy is the difference
            initial_interaction_energy = energy_all - energy_separate
        
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
                
                # Safety check for complex molecules
                if (i > 20 and energy > energy_history[i-1] * 5.0) or not np.isfinite(energy):
                    print(f"Breaking at iteration {i} due to energy instability: {energy}")
                    break
                    
            state = simulation_opt.context.getState(getPositions=True, getEnergy=True)
            minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
            final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
        else:
            # Standard energy minimization with correct force units for tolerance
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
            iterations_performed = min(max_iterations, int((energy_change/tolerance) * 10))
            
            minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        
        # Calculate final interaction energy
        final_interaction_energy = None
        
        if molecule1_atoms and molecule2_atoms:
            # Create context for energy calculation (same as before)
            system_energy = mm.System()
            for i, element in enumerate(elements):
                mass = element_masses.get(element, 0.0)
                system_energy.addParticle(mass)
                
            # Add forces
            system_energy.addForce(mm.HarmonicBondForce())
            system_energy.addForce(mm.HarmonicAngleForce())
            
            # Add NonbondedForce
            nb_force = mm.NonbondedForce()
            for i, element in enumerate(elements):
                params = element_params.get(element, [0.0, 0.3, 0.5])
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
            for i in range(len(all_atoms)):
                for j in range(len(all_atoms)):
                    if i in molecule1_indices and j in molecule1_indices:
                        continue  # Keep molecule1-molecule1 interactions
                    if i in molecule2_indices and j in molecule2_indices:
                        continue  # Keep molecule2-molecule2 interactions
                    nb_force.addException(i, j, 0.0, 0.1, 0.0)  # Turn off cross-molecule interactions
                    
            context_energy.reinitialize(True)
            context_energy.setPositions(minimized_positions * unit.angstrom)
            state_separate = context_energy.getState(getEnergy=True)
            energy_separate = state_separate.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            # Interaction energy is the difference
            final_interaction_energy = energy_all - energy_separate
        
        # Format result with optimized coordinates, separated by molecule
        molecule1_optimized_atoms = []
        molecule2_optimized_atoms = []
        
        for i, atom in enumerate(all_atoms):
            optimized_atom = {
                "id": atom["id"],
                "element": atom["element"],
                "x": float(minimized_positions[i][0]),
                "y": float(minimized_positions[i][1]),
                "z": float(minimized_positions[i][2])
            }
            
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
                "final_energy_kj_mol": final_energy,
                "initial_energy_kj_mol": start_energy,
                "energy_change_kj_mol": abs(start_energy - final_energy),
                "duration_seconds": duration,
                "convergence": "energy_minimized",
                "bonds_detected": len(bonds),
                "intermolecular_bonds": intermolecular_bonds,
                "angles_detected": len(angles),
                "initial_interaction_energy_kj_mol": initial_interaction_energy,
                "interaction_energy_kj_mol": final_interaction_energy
            }
        }
        
    except Exception as e:
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
        
        # Run Hartree-Fock calculation
        mf = scf.RHF(mol)
        energy = mf.kernel()
        
        # Calculate gradient
        g = grad.RHF(mf)
        gradients = g.kernel()
        
        # Perform geometry optimization with parameterized values
        converged = False
        iteration = 0
        current_atoms = atom_list.copy()
        current_energy = energy
        
        while not converged and iteration < max_iterations:
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
            
            # Run Hartree-Fock on new geometry
            new_mf = scf.RHF(new_mol)
            new_energy = new_mf.kernel()
            
            # Calculate new gradient
            new_g = grad.RHF(new_mf)
            new_gradients = new_g.kernel()
            
            # Check convergence using parameterized threshold
            energy_change = abs(new_energy - current_energy)
            gradient_norm = np.linalg.norm(new_gradients)
            
            if energy_change < convergence_threshold and gradient_norm < convergence_threshold:
                converged = True
            
            # Update for next iteration
            current_atoms = new_atom_list
            current_energy = new_energy
            gradients = new_gradients
            iteration += 1
        
        # Create result with optimized atoms
        optimized_atoms = []
        for i, atom in enumerate(atoms):
            optimized_atoms.append({
                "id": atom["id"],
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
                    "step_size": step_size
                },
                "theory_level": f"RHF/{basis}",
                "final_energy_hartree": float(current_energy),
                "iterations": iteration,
                "converged": converged,
                "duration_seconds": duration
            }
        }
        
    except Exception as e:
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
        if molecule1_atoms and molecule2_atoms:
            # Create molecule for combined system
            mol_combined = gto.M(
                atom=all_atom_list,
                basis=basis,
                verbose=0
            )
            
            # Run Hartree-Fock on combined system
            mf_combined = scf.RHF(mol_combined)
            energy_combined = mf_combined.kernel()
            
            # Create and calculate energy for molecule 1
            mol1 = gto.M(
                atom=molecule1_atom_list,
                basis=basis,
                verbose=0
            )
            mf1 = scf.RHF(mol1)
            energy1 = mf1.kernel()
            
            # Create and calculate energy for molecule 2
            mol2 = gto.M(
                atom=molecule2_atom_list,
                basis=basis,
                verbose=0
            )
            mf2 = scf.RHF(mol2)
            energy2 = mf2.kernel()
            
            # Interaction energy is the difference
            initial_interaction_energy = energy_combined - (energy1 + energy2)
        
        # Create molecule with complete system
        mol = gto.M(
            atom=all_atom_list,
            basis=basis,
            verbose=0
        )
        
        # Run Hartree-Fock calculation
        mf = scf.RHF(mol)
        energy = mf.kernel()
        
        # Calculate gradient
        g = grad.RHF(mf)
        gradients = g.kernel()
        
        # Perform geometry optimization with parameterized values
        converged = False
        iteration = 0
        current_atoms = all_atom_list.copy()
        current_energy = energy
        
        while not converged and iteration < max_iterations:
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
            
            # Run Hartree-Fock on new geometry
            new_mf = scf.RHF(new_mol)
            new_energy = new_mf.kernel()
            
            # Calculate new gradient
            new_g = grad.RHF(new_mf)
            new_gradients = new_g.kernel()
            
            # Check convergence using parameterized threshold
            energy_change = abs(new_energy - current_energy)
            gradient_norm = np.linalg.norm(new_gradients)
            
            if energy_change < convergence_threshold and gradient_norm < convergence_threshold:
                converged = True
            
            # Update for next iteration
            current_atoms = new_atom_list
            current_energy = new_energy
            gradients = new_gradients
            iteration += 1
        
        # Calculate final interaction energy
        final_interaction_energy = None
        if molecule1_atoms and molecule2_atoms:
            # Extract optimized coordinates for each molecule
            optimized_coords = np.array([atom[1] for atom in current_atoms])
            
            # Create atom lists for optimized individual molecules
            mol1_opt_atoms = []
            for i in molecule1_indices:
                mol1_opt_atoms.append([current_atoms[i][0], current_atoms[i][1]])
                
            mol2_opt_atoms = []
            for i in molecule2_indices:
                mol2_opt_atoms.append([current_atoms[i][0], current_atoms[i][1]])
            
            # Calculate energy for optimized combined system
            mol_combined_opt = gto.M(
                atom=current_atoms,
                basis=basis,
                verbose=0
            )
            mf_combined_opt = scf.RHF(mol_combined_opt)
            energy_combined_opt = mf_combined_opt.kernel()
            
            # Calculate energy for optimized molecule 1
            mol1_opt = gto.M(
                atom=mol1_opt_atoms,
                basis=basis,
                verbose=0
            )
            mf1_opt = scf.RHF(mol1_opt)
            energy1_opt = mf1_opt.kernel()
            
            # Calculate energy for optimized molecule 2
            mol2_opt = gto.M(
                atom=mol2_opt_atoms,
                basis=basis,
                verbose=0
            )
            mf2_opt = scf.RHF(mol2_opt)
            energy2_opt = mf2_opt.kernel()
            
            # Interaction energy is the difference
            final_interaction_energy = energy_combined_opt - (energy1_opt + energy2_opt)
        
        # Create result with optimized atoms
        molecule1_optimized_atoms = []
        molecule2_optimized_atoms = []
        
        for i, atom in enumerate(all_atoms):
            optimized_atom = {
                "id": atom["id"],
                "element": atom["element"],
                "x": float(current_atoms[i][1][0]),
                "y": float(current_atoms[i][1][1]),
                "z": float(current_atoms[i][1][2])
            }
            
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
                    "step_size": step_size
                },
                "molecules": 2 if molecule1_atoms and molecule2_atoms else 1,
                "molecule1_atom_count": len(molecule1_atoms) if molecule1_atoms else 0,
                "molecule2_atom_count": len(molecule2_atoms) if molecule2_atoms else 0,
                "theory_level": f"RHF/{basis}",
                "final_energy_hartree": float(current_energy),
                "iterations": iteration,
                "converged": converged,
                "duration_seconds": duration,
                "initial_interaction_energy_hartree": float(initial_interaction_energy) if initial_interaction_energy is not None else None,
                "interaction_energy_hartree": float(final_interaction_energy) if final_interaction_energy is not None else None
            }
        }
        
    except Exception as e:
        return {
            "error": str(e),
            "metadata": {
                "method": "quantum_chemistry",
                "library": "PySCF",
                "status": "failed"
            }
        }

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
        
        # Validate inputs
        if not molecule1_data and not molecule2_data:
            return jsonify({"error": "Invalid input: At least one molecule must be provided"}), 400
            
        if interaction_mode and (not molecule1_data or not molecule2_data):
            return jsonify({"error": "Invalid input: Both molecules must be provided for interaction mode"}), 400
            
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
                # Convert string user_id back to integer for database query
                user = User.query.get(int(user_identity))
                if user and user.subscription_status == "active":
                    is_subscribed = True
                    user_id = user.id
        except Exception as e:
            # If token verification fails, proceed as unauthenticated
            print(f"JWT verification error (proceeding as unauthenticated): {str(e)}")
            pass
        
        # Extract atoms from molecule data
        try:
            molecule1_atoms = None
            if molecule1_data:
                if "file1" in molecule1_data:
                    molecule1_atoms = molecule1_data["file1"]["atoms"]
                elif "atoms" in molecule1_data:
                    molecule1_atoms = molecule1_data["atoms"]
            
            molecule2_atoms = None
            if molecule2_data:
                if "file1" in molecule2_data:
                    molecule2_atoms = molecule2_data["file1"]["atoms"]
                elif "atoms" in molecule2_data:
                    molecule2_atoms = molecule2_data["atoms"]
                
            # Ensure at least one molecule has atoms
            if not molecule1_atoms and not molecule2_atoms:
                return jsonify({"error": "Invalid molecule data: No atoms found in either molecule"}), 400
                
        except Exception as e:
            return jsonify({"error": f"Error parsing molecule data: {str(e)}"}), 400
        
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
        
        # Store optimization results in database
        optimization = Optimization(
            user_id=user_id,
            optimization_type=optimization_type,
            parameters=json.dumps({
                "molecule1": molecule1_data,
                "molecule2": molecule2_data,
                "optimization_params": optimization_params,
                "interaction_mode": interaction_mode
            }),
            result=json.dumps(result)
        )
        db.session.add(optimization)
        db.session.commit()
        
        # Return results
        return jsonify({
            "success": True,
            "optimizationId": optimization.id,
            "result": result,
            "original_molecule1": molecule1_data,
            "original_molecule2": molecule2_data
        }), 200
        
    except Exception as e:
        return jsonify({"error": str(e)}), 500