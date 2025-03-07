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

# Create blueprint
opti_bp = Blueprint('opti', __name__)

class Optimization(db.Model):
    """Optimization record model for tracking molecular optimizations."""
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    optimization_type = db.Column(db.String(50), nullable=False)  # 'classical' or 'quantum'
    parameters = db.Column(db.Text, nullable=False)  # JSON string of parameters
    result = db.Column(db.Text, nullable=True)  # JSON string of the result
    created_at = db.Column(db.DateTime, default=db.func.current_timestamp())

    def __repr__(self):
        return f'<Optimization {self.id} for User {self.user_id}>'

def optimize_classical(atoms):
    """
    Optimize molecule using classical molecular dynamics with OpenMM.
    Template-free approach for arbitrary molecules.
    """
    try:   
        start_time = datetime.now()
        
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
                if dist < 0.2:  # Bond threshold in nm
                    topology.addBond(atom_objects[i], atom_objects[j])
                    bonds.append((i, j, dist))
        
        # Add HarmonicBondForce
        bond_force = mm.HarmonicBondForce()
        for i, j, dist in bonds:
            # Default parameters: K = 1000 kJ/mol/nm^2
            bond_force.addBond(i, j, dist, 1000.0)
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
                        # Add angle with K = 500 kJ/mol/radian^2
                        angle_force.addAngle(atom_i, j, atom_k, angle, 500.0)
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
        
        # Create Integrator (300K, 1 ps⁻¹ friction, 2 fs timestep)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        
        # Create simulation
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        
        # Minimize energy
        simulation.minimizeEnergy(maxIterations=1000)
        
        # Get minimized positions and energy
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Format result
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
                "final_energy_kj_mol": final_energy,
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
    
def optimize_quantum(atoms):
    """
    Optimize molecule using quantum chemistry with PySCF.
    
    Args:
        atoms: List of atom dictionaries with element and coordinates
        
    Returns:
        Dictionary containing optimized atoms and metadata
    """
    try:       
        start_time = datetime.now()
        
        # Convert atoms to PySCF format
        atom_list = []
        for atom in atoms:
            atom_list.append([atom["element"], (atom["x"], atom["y"], atom["z"])])
        
        # Create molecule
        mol = gto.M(
            atom=atom_list,
            basis='6-31g',
            verbose=0
        )
        
        # Run Hartree-Fock calculation
        mf = scf.RHF(mol)
        energy = mf.kernel()
        
        # Calculate gradient
        g = grad.RHF(mf)
        gradients = g.kernel()
        
        # Initialize geometry optimizer
        max_iterations = 10
        convergence_threshold = 1e-5
        
        # Perform geometry optimization
        converged = False
        iteration = 0
        current_atoms = atom_list.copy()
        current_energy = energy
        
        while not converged and iteration < max_iterations:
            # Get coordinates
            coords = np.array([atom[1] for atom in current_atoms])
            
            # Update geometry based on gradient
            step_size = 0.1
            new_coords = coords - step_size * gradients
            
            # Update atom list
            new_atom_list = []
            for i, atom in enumerate(current_atoms):
                new_atom_list.append([atom[0], tuple(new_coords[i])])
            
            # Create new molecule
            new_mol = gto.M(
                atom=new_atom_list,
                basis='6-31g',
                verbose=0
            )
            
            # Run Hartree-Fock on new geometry
            new_mf = scf.RHF(new_mol)
            new_energy = new_mf.kernel()
            
            # Calculate new gradient
            new_g = grad.RHF(new_mf)
            new_gradients = new_g.kernel()
            
            # Check convergence
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
                "theory_level": "RHF/6-31G",
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

@opti_bp.route('/optimize-molecule', methods=['POST'])
def optimize_molecule():
    """
    Endpoint to optimize molecular structures using both classical and quantum methods.
    Requires an active subscription.
    
    Expected request format:
    {
        "email": "user@example.com",
        "molecule": {
            "file1": {
                "atoms": [
                    {"id": 1, "element": "C", "x": 1.1, "y": 0.0, "z": 0.0},
                    ...
                ]
            }
        }
    }
    """
    try:
        # Check if user is authenticated and subscribed
        data = request.get_json()
        if not data:
            return jsonify({"error": "Invalid input: No data provided"}), 400
            
        email = data.get("email")
        molecule_data = data.get("molecule")
        
        if not email or not molecule_data:
            return jsonify({"error": "Invalid input: Missing 'email' or 'molecule' data"}), 400
            
        # Verify subscription status
        user = User.query.filter_by(email=email).first()
        if not user or user.subscription_status != "active":
            return jsonify({"error": "Active subscription required for molecule optimization"}), 403
            
        # Extract atoms from molecule data
        try:
            # Handle structure from example file format
            if "file1" in molecule_data:
                atoms = molecule_data["file1"]["atoms"]
            # Handle direct atoms array
            elif "atoms" in molecule_data:
                atoms = molecule_data["atoms"]
            else:
                return jsonify({"error": "Invalid molecule data structure"}), 400
                
            if not atoms:
                return jsonify({"error": "Invalid molecule data: No atoms found"}), 400
        except Exception as e:
            return jsonify({"error": f"Error parsing molecule data: {str(e)}"}), 400
            
        # Perform classical optimization
        classical_result = optimize_classical(atoms)
        
        # Perform quantum optimization
        quantum_result = optimize_quantum(atoms)
        
        # Store optimization results in database
        optimization = Optimization(
            user_id=user.id,
            optimization_type="both",  # Both classical and quantum
            parameters=json.dumps(molecule_data),
            result=json.dumps({
                "classical": classical_result,
                "quantum": quantum_result
            })
        )
        db.session.add(optimization)
        db.session.commit()
        
        # Return results
        return jsonify({
            "success": True,
            "optimizationId": optimization.id,
            "classical": classical_result,
            "quantum": quantum_result,
            "original_molecule": molecule_data
        }), 200
        
    except Exception as e:
        return jsonify({"error": str(e)}), 500