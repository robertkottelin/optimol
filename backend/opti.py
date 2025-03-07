from flask import Blueprint, request, jsonify
import json
import numpy as np
from datetime import datetime
from app import db
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
    
    Args:
        atoms: List of atom dictionaries with element and coordinates
        
    Returns:
        Dictionary containing optimized atoms and metadata
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
        
        atom_objects = []
        for i, element in enumerate(elements):
            atom_obj = topology.addAtom(element, app.Element.getBySymbol(element), residue)
            atom_objects.append(atom_obj)
        
        # Add bonds (simplified approach - determining bonds based on distance)
        positions_nm = positions.value_in_unit(unit.nanometer)
        for i in range(len(atom_objects)):
            for j in range(i+1, len(atom_objects)):
                dist = np.sqrt(np.sum((positions_nm[i] - positions_nm[j])**2))
                # Typical bond length is ~0.15 nm
                if dist < 0.2:  # Bond threshold
                    topology.addBond(atom_objects[i], atom_objects[j])
        
        # Create system with force field
        forcefield = app.ForceField('amber14-all.xml')
        system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff)
        
        # Create integrator for minimization
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        
        # Create simulation
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        
        # Minimize energy
        simulation.minimizeEnergy(maxIterations=1000)
        
        # Get minimized positions
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        minimized_positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Create result with optimized atoms
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
                "convergence": "energy_minimized"
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
        max_iterations = 50
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