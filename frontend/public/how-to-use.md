# Molecular Optimization System - User Guide

## System Overview

The Molecular Optimization System provides computational tools for structural optimization of molecular systems. The application implements both classical molecular mechanics and quantum chemistry methodologies to perform energy minimization calculations, enabling:
- Single-molecule optimization with selectable force field and quantum mechanical methods
- Dual-molecule interaction analysis with precise positional control
- Energy calculation with geometry convergence verification
- Export of optimized structural coordinates for external analysis workflows

This guide provides comprehensive instructions for effective system utilization across all functional components.

## Authentication System

### Access Tiers

The system implements a JWT-based authentication system with the following access tiers:

**Guest Access**: 
Limited computational resources, parameter constraints. No authentication required.

**Registered Account**: 
Standard computational access, persistent data storage. Requires email and password.

**Premium Subscription**: 
Advanced feature set, expanded computational limits. Requires email, password, and payment method.

### Registration Protocol

1. Select "Login / Register" in the application header
2. Choose "Register now" on the authentication modal
3. Input email address and secure password (8+ characters)
4. Submit credentials to create account

### Subscription Specifications

**Free Tier Limitations**:
- Classical optimization: Maximum 1,000 iterations
- Quantum optimization: Maximum 5 iterations
- Basis set access: Limited to STO-3G and 6-31G
- System size constraints based on optimization method

**Premium Tier ($10/month)**:
- Classical optimization: Maximum 1,000,000 iterations
- Quantum optimization: Maximum 1,000 iterations
- Full basis set access: STO-3G, 6-31G, 6-311G, cc-pVDZ
- Increased system size allowances
- Priority computational task scheduling

## Molecular Data Input

### Supported File Formats

The system accepts JSON-formatted molecular data with the following schema structures:

**Standard Format**:
```json
{
  "file1": {
    "atoms": [
      {"id": 1, "element": "C", "x": 1.1, "y": 0.0, "z": 0.0},
      {"id": 2, "element": "H", "x": 2.1, "y": 0.0, "z": 0.0}
    ],
    "metadata": {
      "name": "Molecule Name",
      "formula": "Chemical Formula"
    }
  }
}
```

**Alternative Format**:
```json
{
  "atoms": [
    {"id": 1, "element": "C", "x": 1.1, "y": 0.0, "z": 0.0},
    {"id": 2, "element": "H", "x": 2.1, "y": 0.0, "z": 0.0}
  ]
}
```

### Import Methods

1. **Direct File Upload**:
   - Click the upload area or drag-and-drop a JSON file
   - System validates JSON structure and molecular data integrity
   - File size limit: 5MB
   - Coordinate units: Ångströms

2. **Test Molecule Selection**:
   - Click on any predefined test molecule button:
     - Water (H₂O): 3-atom system
     - Acetic Acid (CH₃COOH): 8-atom system
     - Methanol (CH₃OH): 6-atom system
     - Ibuprofen (C₁₃H₁₈O₂): 31-atom system
   - The molecule will be automatically loaded into the active slot

3. **Molecule Selection Control**:
   - Toggle between Molecule 1 and Molecule 2 using the selector tabs
   - The active molecule slot is highlighted with its respective color
   - Clear individual molecules using the "Clear" button
   - Clear all molecules using the "Clear All Molecules" button

## Molecular Visualization Interface

### Viewer Controls

The 3D molecular viewer implements the following interaction modes:

**Rotation**: Click and drag with mouse to rotate the view
**Zoom**: Scroll mouse wheel to zoom in/out
**Pan**: Press and drag with middle button to pan
**Reset View**: Double-click to reset to default orientation

### Molecular Display Features

- Element-specific atomic labels
- Color-coded molecules (blue: molecule 1, green: molecule 2)
- Bond visualization with length measurements (Å)
- Hydrogen bond representation via dashed line visualization
- Intermolecular hydrogen bonds shown in gold color

### Bond Display Controls

- Toggle bond lengths display with the checkbox
- Toggle controls and bond type legend with the checkbox
- Bond type legend shows:
  - Covalent bonds (solid lines)
  - Hydrogen bonds (white dashed lines)
  - Intermolecular hydrogen bonds (gold dashed lines)

## Optimization Methods

### Method Selection

The system provides two primary optimization methodologies:

**Classical Molecular Mechanics**:
- Force field-based optimization using OpenMM engine
- Suitable for large molecular systems (100+ atoms)
- Faster computation with lower precision
- Optimizes molecular geometry using empirical potentials

**Quantum Chemistry Methods**:
- Ab initio quantum mechanical calculations using PySCF
- Suitable for smaller molecular systems (<30 atoms)
- Higher computational cost with improved accuracy
- Accounts for electronic structure effects

### Selection Guidelines

**Classical Optimization Recommended For**:
- Systems with >30 atoms
- Moderate precision requirements
- Higher priority on computational speed
- Initial structure refinement
- Large molecular complexes

**Quantum Optimization Recommended For**:
- Systems with <30 atoms
- High precision requirements
- Electronic property analysis
- Small molecule interactions
- Systems where electronic effects predominate

## Parameter Configuration

### Classical Parameters

**Temperature (K)**:
- Function: Simulation temperature
- Default: 300K
- Range: 1-1000K

**Max Iterations**:
- Function: Optimization cycle limit
- Default: 1000
- Range: 100-1000000 (subscription dependent)

**Force All Iterations**:
- Function: Override convergence criteria
- Default: False
- Options: True/False

**Advanced Parameters**:
- Convergence Tolerance (kJ/mol): Energy gradient threshold (default: 10.0)
- Bond Threshold (nm): Maximum bond detection distance (default: 0.2)
- Bond Force Constant (kJ/mol/nm²): Spring constant for bond terms (default: 1000.0)
- Angle Force Constant (kJ/mol/rad²): Spring constant for angle terms (default: 500.0)

### Quantum Parameters

**Basis Set**:
- Function: Atomic orbital basis functions
- Default: 6-31G
- Options: 
  - All users: STO-3G, 6-31G
  - Premium only: 6-311G, cc-pVDZ

**Max Iterations**:
- Function: Geometry optimization steps
- Default: 10
- Range: 
  - Free tier: 1-5
  - Premium tier: 1-1000

**Advanced Parameters**:
- Convergence Threshold: Energy and gradient convergence criteria (default: 0.00001)
- Step Size: Geometric step magnitude in optimization (default: 0.1)

## Molecular Interaction Analysis

### Enabling Interaction Mode

1. Load two molecular structures (one in each slot)
2. Enable "Optimize Molecular Interaction" checkbox
3. Note: Molecule 1 serves as reference frame; Molecule 2 is positionally adjustable
4. Option to apply default offset (5Å) when loading Molecule 2

### Positional Control System

To access positional controls:
1. Enable "Enable Positioning Mode" button
2. Use the following controls:

**Keyboard Navigation**:
- Arrow keys: XY-plane translation (±0.5Å per keystroke)
- Page Up/Down: Z-axis translation (±0.5Å per keystroke)
- Shift + arrows: Rotational control (±15° per keystroke)

**UI Controls**:
- XYZ position input fields: Direct coordinate input
- XYZ rotation input fields: Direct angle input
- Reset Position button: Reset to zero offset
- Reset Rotation button: Reset to zero rotation

The system displays current coordinates and rotational angles with decimal notation precision.

## Optimization Execution

### Execution Workflow

1. Verify molecular structure(s) are loaded properly
2. Select optimization method (Classical/Quantum)
3. Configure parameter settings as needed
4. Initiate optimization via "Run X Optimization" button

### Progress Tracking

During computation:
- Progress indicator displays "Optimizing..." state
- UI remains responsive but parameter changes are disabled
- Computational duration depends on:
  - Molecular complexity (atom count)
  - Optimization method selection
  - Parameter configuration (iterations)
  - Server resource availability/load

### Error Handling

Common optimization errors and resolutions:
- Invalid molecular structure: Verify JSON format
- Convergence failure: Increase iteration count
- Size limitation exceeded: Switch to classical method for large systems
- Subscription limits: Consider premium upgrade for extended capabilities

## Results Analysis

### Results Display

Upon successful optimization, results are presented in the Results panel:

**Classical Optimization Results**:
- Method identification and library
- Temperature and iteration metrics
- Energy values (kJ/mol): initial, final, delta
- Structural analysis: bond and angle counts
- Computational performance metrics
- For interaction mode: interaction energy, intermolecular bonds

**Quantum Optimization Results**:
- Theory level and basis set specification
- Final energy in Hartree units
- Iteration count with convergence status
- Computational duration
- For interaction mode: intermolecular interaction energy

### Visualization Toggle

Switch between original and optimized structures:
- Original: Pre-optimization molecular configuration
- Optimized: Post-optimization molecular geometry
- For interaction analyses: Both molecules displayed in optimized state

### Data Export

Export optimized structures via the "Download Result" button:
- File format: JSON with identical schema to input
- Contents:
  - Optimized atomic coordinates
  - Complete parameter and metadata
  - Energy and convergence metrics
  - Timestamp and optimization method details

## Troubleshooting Guide

### Common Issues

**File upload failures**:
- Verify JSON format matches required schema
- Check for malformed coordinates or missing element definitions
- Ensure file size is below 5MB limit

**Slow computation**:
- Reduce system size or iteration count
- Switch to classical method for large systems
- Check server health status

**Position controls unresponsive**:
- Ensure positioning mode is enabled
- Verify both molecules are loaded
- Try resetting browser cache

**Convergence failures**:
- Increase iterations or adjust convergence criteria
- Verify initial structure has reasonable geometry
- Try classical pre-optimization before quantum optimization

**Visualization errors**:
- Try clearing and reloading molecules
- Ensure WebGL is enabled in browser
- Update browser to latest version

### Server Health Check

The system provides a server health check feature:
1. Click "Server Health" button in top navigation
2. System will query backend services
3. Results appear as:
   - Green: All systems operational
   - Red: Service disruption detected
   - Details provided for diagnostic purposes

## Interface Reference

### Mobile Device Considerations

The system implements responsive design with mobile-specific optimizations:
- Collapsed menu available via menu button
- Simplified controls for touch interfaces
- Reduced visualization complexity
- Vertical layout optimization

### Keyboard Shortcuts

**Molecule 2 X-axis position**: Left/Right arrow keys
**Molecule 2 Y-axis position**: Up/Down arrow keys
**Molecule 2 Z-axis position**: Page Up/Down keys
**Molecule 2 X-axis rotation**: Shift + Up/Down arrow keys
**Molecule 2 Y-axis rotation**: Shift + Left/Right arrow keys
**Molecule 2 Z-axis rotation**: Shift + Page Up/Down keys
