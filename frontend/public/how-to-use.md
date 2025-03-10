# Optimol Molecular Optimization System Documentation

## 1. System Overview

The Optimol Molecular Optimization System provides computational tools for molecular structure optimization using both classical molecular dynamics and quantum chemical methods. This document outlines operational procedures, input specifications, and theoretical foundations.

## 2. Input Specifications

### 2.1 File Format

Input files must be JSON formatted with the following structure:

```json
{
  "file1": {
    "atoms": [
      {
        "id": 1,
        "element": "C",
        "x": 0.0,
        "y": 0.0,
        "z": 0.0
      },
      {
        "id": 2,
        "element": "H",
        "x": 1.09,
        "y": 0.0,
        "z": 0.0
      }
      // Additional atoms
    ]
  }
}
```

Alternatively:

```json
{
  "atoms": [
    {
      "id": 1,
      "element": "C",
      "x": 0.0,
      "y": 0.0,
      "z": 0.0
    },
    // Additional atoms
  ]
}
```

### 2.2 Required Fields

| Field | Type | Description |
|-------|------|-------------|
| id | Integer | Unique atom identifier |
| element | String | Chemical element symbol |
| x | Number | X-coordinate in Angstroms |
| y | Number | Y-coordinate in Angstroms |
| z | Number | Z-coordinate in Angstroms |

### 2.3 Supported Elements

The system supports standard chemical elements including but not limited to: H, C, N, O, F, P, S, Cl, Br, I.

### 2.4 Example Input File

```json
{
  "file1": {
    "atoms": [
      {"id": 1, "element": "C", "x": 0.000, "y": 0.000, "z": 0.000},
      {"id": 2, "element": "H", "x": 0.000, "y": 0.000, "z": 1.089},
      {"id": 3, "element": "H", "x": 1.026, "y": 0.000, "z": -0.363},
      {"id": 4, "element": "H", "x": -0.513, "y": -0.889, "z": -0.363},
      {"id": 5, "element": "H", "x": -0.513, "y": 0.889, "z": -0.363}
    ]
  }
}
```

This example represents methane (CH₄).

## 3. Application Usage

### 3.1 Workflow Overview

1. Upload molecular structure file
2. Select optimization method
3. Configure optimization parameters
4. Execute optimization
5. Visualize and analyze results
6. Download optimized structure

### 3.2 Optimization Method Selection

The system provides two optimization methods:

1. **Classical Optimization**: Uses molecular mechanics force fields
2. **Quantum Optimization**: Employs quantum chemistry calculations

### 3.3 Parameter Configuration

#### 3.3.1 Classical Parameters

| Parameter | Description | Range | Default |
|-----------|-------------|-------|---------|
| Temperature | Simulation temperature (K) | 1-1000 | 300 |
| Max Iterations | Maximum energy minimization steps | 100-1000 | 1000/100* |
| Bond Threshold | Distance threshold for bond detection (nm) | 0.1-0.5 | 0.2 |
| Bond Force Constant | Bond strength parameter (kJ/mol/nm²) | 100-10000 | 1000.0 |
| Angle Force Constant | Angle stiffness parameter (kJ/mol/rad²) | 50-5000 | 500.0 |

#### 3.3.2 Quantum Parameters

| Parameter | Description | Range | Default |
|-----------|-------------|-------|---------|
| Basis Set | Atomic orbital basis functions | STO-3G, 6-31G, 6-311G**, cc-pVDZ** | 6-31G |
| Max Iterations | Maximum optimization steps | 1-10/3* | 10/3* |
| Convergence Threshold | Energy convergence criterion | 0.000001-0.01 | 0.00001 |
| Step Size | Geometry update magnitude | 0.01-1.0 | 0.1 |

*Subscription-dependent limits  
**Available only with subscription

### 3.4 Result Interpretation

#### 3.4.1 Classical Optimization Results

- Final energy (kJ/mol)
- Number of detected bonds and angles
- Convergence status
- Computation duration

#### 3.4.2 Quantum Optimization Results

- Final energy (Hartree)
- Iteration count
- Convergence status
- Computation duration

### 3.5 Output Format

The optimized structure is available for download in the same JSON format as the input, with updated atomic coordinates and additional metadata.

## 4. Theoretical Background

### 4.1 Classical Molecular Optimization

#### 4.1.1 Methodology

The classical optimization employs molecular mechanics principles using the OpenMM library. The system implements a template-free approach suitable for arbitrary molecules.

#### 4.1.2 Force Field Implementation

The molecular system is parameterized using:

1. **HarmonicBondForce**: Models covalent bonds as harmonic springs
   ```
   E_bond = k_b(r - r_0)²
   ```
   where k_b is the bond force constant, r is the bond length, and r_0 is the equilibrium bond length.

2. **HarmonicAngleForce**: Models bond angles as harmonic potentials
   ```
   E_angle = k_θ(θ - θ_0)²
   ```
   where k_θ is the angle force constant, θ is the bond angle, and θ_0 is the equilibrium angle.

3. **NonbondedForce**: Models van der Waals and electrostatic interactions
   ```
   E_nonbonded = 4ε[(σ/r)¹² - (σ/r)⁶] + q₁q₂/4πε₀r
   ```
   where ε and σ are Lennard-Jones parameters, q₁ and q₂ are charges, and r is the interatomic distance.

#### 4.1.3 Energy Minimization

The optimization uses the Langevin integrator in combination with an energy minimization algorithm to find the local energy minimum. The temperature parameter affects the stochastic dynamics during energy minimization.

### 4.2 Quantum Chemical Optimization

#### 4.2.1 Methodology

Quantum optimization employs ab initio quantum chemistry methods implemented in PySCF. The system uses Hartree-Fock (HF) theory with specified basis sets.

#### 4.2.2 Hartree-Fock Theory

The HF method approximates the many-electron wavefunction as a Slater determinant of single-electron orbitals. The electronic energy is:

```
E = Σᵢ hᵢᵢ + 1/2 Σᵢⱼ (2Jᵢⱼ - Kᵢⱼ)
```

where:
- hᵢᵢ represents one-electron integrals
- Jᵢⱼ represents Coulomb integrals
- Kᵢⱼ represents exchange integrals

#### 4.2.3 Basis Sets

The available basis sets provide different levels of accuracy:

- **STO-3G**: Minimal basis set, fastest but least accurate
- **6-31G**: Standard split-valence basis, good balance of speed and accuracy
- **6-311G**: Extended split-valence basis, improved accuracy (subscription only)
- **cc-pVDZ**: Correlation-consistent basis, high accuracy (subscription only)

#### 4.2.4 Geometry Optimization

The optimization employs a gradient-based method where:

1. The electronic energy is calculated using the HF method
2. Energy gradients (forces) with respect to nuclear coordinates are computed
3. Molecular geometry is updated using:
   ```
   R_new = R_old - step_size × gradient
   ```
4. Process repeats until convergence criteria are satisfied:
   - Energy change < convergence_threshold
   - Gradient norm < convergence_threshold

## 5. Subscription Features

The application implements a tiered feature model:

### 5.1 Free Account Limitations

- Classical optimization: Max 100 iterations
- Quantum optimization: Max 3 iterations
- Basis set selection limited to STO-3G and 6-31G

### 5.2 Subscribed Account Features

- Classical optimization: Max 1000 iterations
- Quantum optimization: Max 10 iterations
- Full basis set access including 6-311G and cc-pVDZ
- Enhanced accuracy and convergence capabilities

## 6. Technical Implementation

The system architecture consists of:

- Frontend: React.js with 3DMol.js for molecular visualization
- Backend: Flask with OpenMM and PySCF for computational operations
- Database: SQLite for user and optimization record storage
- Authentication: Stripe for subscription management

## 7. Troubleshooting

### 7.1 Common Issues

- **Invalid JSON format**: Ensure structure follows required schema
- **Optimization failure**: Consider adjusting parameters or simplifying molecular structure
- **Visualization issues**: Verify that atom coordinates are in proper units (Angstroms)

### 7.2 Known Limitations

- Maximum molecule size: ~100 atoms for quantum calculations
- Supported atom types: Main group elements
- No support for metal complexes or periodic systems

### 7.3 Notes
- Developers are pushing to prod constantly so please save your work!