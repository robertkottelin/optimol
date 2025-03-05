## How To Use & Theory

This system performs molecular structure optimization and binding simulations using CP2K, a quantum chemistry and solid state physics software package implementing DFT-based electronic structure methods.

---

### Technical Overview

The application implements a multi-scale computational framework providing:

1. **DFT-Based Structural Optimization**: Geometry optimization using Kohn-Sham Density Functional Theory
2. **QM/MM Binding Simulations**: Hybrid quantum mechanics/molecular mechanics for protein-ligand interactions
3. **Tiered Computational Parameters**: Resource allocation based on subscription status

---

### Input File Specification

#### Single Molecule JSON Schema

```json
{
  "file1": {
    "atoms": [
      { "id": 1, "element": "H", "x": 0.0, "y": 0.0, "z": 0.0 },
      { "id": 2, "element": "H", "x": 0.0, "y": 0.0, "z": 0.74 },
      { "id": 3, "element": "O", "x": 0.0, "y": 0.65, "z": 0.0 }
    ]
  }
}
```

#### Protein-Ligand JSON Schema (Two Separate Files)

Protein file:
```json
{
  "file1": {
    "atoms": [
      { "id": 1, "element": "C", "x": 0.0, "y": 0.0, "z": 0.0 },
      { "id": 2, "element": "N", "x": 1.32, "y": 0.0, "z": 0.0 },
      ...
    ]
  }
}
```

Ligand file (separate upload):
```json
{
  "file1": {
    "atoms": [
      { "id": 1, "element": "C", "x": 10.0, "y": 5.0, "z": 3.0 },
      { "id": 2, "element": "O", "x": 11.2, "y": 5.2, "z": 3.0 },
      ...
    ]
  }
}
```

---

### Computational Parameters

#### 1. Structural Optimization Parameters

| Parameter | Description | Default | Standard Range | Premium Range |
|-----------|-------------|---------|---------------|---------------|
| `fmax` | Force convergence criterion (eV/Å) | 0.005 | Fixed: 0.05 | 0.001-0.1 |
| `steps` | Maximum geometry optimization steps | 100 | 10-100 | 10-1000 |
| `basis_set` | DFT basis set | SZV-MOLOPT-GTH | SZV only | Multiple options |
| `functional` | Exchange-correlation functional | PBE | PBE only | PBE, BLYP, B3LYP, PBE0 |
| `charge` | System charge | 0 | Fixed: 0 | Any integer |
| `spin_polarized` | Unrestricted Kohn-Sham | false | Fixed: false | true/false |

#### 2. Binding Simulation Parameters (Premium Only)

| Parameter | Description | Default | Range |
|-----------|-------------|---------|-------|
| `fmax` | Force convergence (QM/MM) | 0.05 | 0.01-0.1 |
| `steps` | Maximum optimization steps | 100 | 10-1000 |
| `qm_indices` | Optional QM atom indices | Auto (ligand) | Custom selection |

---

### Technical Methodology

#### 1. DFT-Based Structure Optimization

CP2K implements the **Gaussian and Plane Waves** (GPW) method for efficient DFT calculations:

- **Electron Density**: Represented on a real-space grid using plane waves
- **Orbitals**: Represented using Gaussian basis functions
- **Core-Electron Treatment**: GTH pseudopotentials
- **SCF Procedure**: Converges electron density to minimize total energy

The optimization process uses the **BFGS** algorithm to minimize forces on atoms:

1. Calculate energy and forces using CP2K's DFT implementation
2. Update atomic positions based on BFGS quasi-Newton algorithm
3. Iterate until forces fall below threshold (`fmax`) or maximum steps reached
4. Return optimized structure and energy

#### 2. QM/MM Binding Simulations (Premium Feature)

QM/MM simulations divide the system into regions treated at different levels of theory:

- **QM Region**: Typically the ligand and binding site residues
- **MM Region**: Remainder of protein and solvent environment

The technical implementation uses:

1. Electronic embedding with electrostatic coupling between regions
2. DFT treatment of QM region with selected functional and basis set
3. Classical force field treatment of MM region
4. Energy minimization of the combined system
5. Extraction of binding energy and optimized structures

---

### Computational Theory

#### Density Functional Theory

CP2K implements Kohn-Sham DFT, which solves for the ground state electron density ρ(r) by minimizing the energy functional:

E[ρ] = T[ρ] + Eext[ρ] + EH[ρ] + Exc[ρ]

Where:
- T[ρ]: Kinetic energy of non-interacting electrons
- Eext[ρ]: External potential energy (electron-nuclei interaction)
- EH[ρ]: Hartree energy (classical electron-electron repulsion)
- Exc[ρ]: Exchange-correlation energy

#### Exchange-Correlation Functionals

Available functionals represent different approximations to Exc[ρ]:

- **PBE**: Generalized gradient approximation (GGA) with efficient scaling
- **BLYP**: Combined Becke exchange with Lee-Yang-Parr correlation
- **B3LYP**: Hybrid functional with exact exchange mixing (Premium)
- **PBE0**: Hybrid functional with 25% exact exchange (Premium)

#### Basis Sets

CP2K uses molecularly optimized (MOLOPT) basis sets of increasing accuracy:

- **SZV**: Single-zeta valence (minimal basis, fastest)
- **DZVP**: Double-zeta valence with polarization functions
- **TZVP**: Triple-zeta valence with polarization functions
- **TZV2P**: Triple-zeta valence with double polarization (highest accuracy)

---

### Workflow Instructions

#### Single Molecule Optimization

1. Navigate to "Single Molecule" tab
2. Upload molecular structure JSON
3. Select computational parameters:
   - fmax: Force convergence criterion
   - steps: Maximum optimization steps
   - basis_set: Select appropriate basis functions
   - functional: Select exchange-correlation functional
   - charge: Set molecular charge
   - spin_polarized: Enable for open-shell systems
4. Execute "Classical Optimize" for standard DFT calculation
5. Premium users can select "Quantum Optimize" for enhanced accuracy
6. Download optimized structure and energy data

#### Protein-Ligand Binding (Premium Feature)

1. Navigate to "Protein-Ligand Binding" tab
2. Upload separate protein and ligand structure files
3. Configure binding parameters:
   - fmax: Force convergence criterion
   - steps: Maximum optimization steps
4. Execute "Optimize Binding" to run QM/MM simulation
5. Download binding results with optimized structures and binding energy

---

### Output Data Format

#### Single Molecule Optimization Output

```json
{
  "atoms": [
    { "id": 1, "element": "H", "x": 0.012, "y": 0.003, "z": -0.087 },
    { "id": 2, "element": "H", "x": -0.002, "y": 0.005, "z": 0.962 },
    { "id": 3, "element": "O", "x": -0.005, "y": 0.744, "z": 0.001 }
  ],
  "energy": -75.091253,
  "converged": true,
  "parameters": {
    "fmax": 0.005,
    "steps": 200,
    "basis_set": "DZVP-MOLOPT-GTH",
    "functional": "PBE",
    "charge": 0,
    "spin_polarized": false
  }
}
```

#### Binding Simulation Output

```json
{
  "protein": {
    "atoms": [
      { "id": 1, "element": "C", "x": 0.127, "y": 0.043, "z": 0.012 },
      ...
    ]
  },
  "ligand": {
    "atoms": [
      { "id": 1, "element": "C", "x": 9.872, "y": 5.126, "z": 3.041 },
      ...
    ]
  },
  "binding_energy": -12.347,
  "converged": true
}
```

---

### Technical Notes

1. **Resource Allocation**: Computational intensity scales with:
   - System size (number of atoms)
   - Basis set complexity (SZV < DZVP < TZVP < TZV2P)
   - Functional complexity (GGA < hybrid)
   - QM region size in QM/MM simulations

2. **Accuracy Considerations**:
   - DFT calculations approximate electron correlation
   - Larger basis sets reduce basis set superposition error
   - Hybrid functionals generally improve accuracy at higher computational cost
   - QM/MM boundaries may introduce artificial effects

3. **Performance Optimization**:
   - Use minimal basis sets for initial optimization
   - Progress to more complete basis sets for final refinement
   - Select appropriate functional for specific chemical properties
   - Optimize QM region selection for binding studies

---

### Development Roadmap

- Integration of additional CP2K capabilities:
  - Ab initio molecular dynamics
  - Excited state calculations (TDDFT)
  - Advanced sampling methods
  - Reaction path optimization
- Expanded visualization capabilities
- Custom force field implementation for MM region
- Machine learning acceleration for QM calculations