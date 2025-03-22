import React, { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';
import { styles } from '../styles/components';

const MoleculeViewer = ({
  atoms,
  isMobile,
  positioningMode,
  onMoleculeMove,
  molecule2Offset,
  molecule2Rotation,
  onMoleculeRotate
}) => {
  const viewerRef = useRef();
  const containerRef = useRef();
  const [isDragging, setIsDragging] = useState(false);
  const [dragStartPos, setDragStartPos] = useState({ x: 0, y: 0 });
  const [viewerInstance, setViewerInstance] = useState(null);
  const [isInitialRender, setIsInitialRender] = useState(true);
  const viewerIdRef = useRef(`molecule-viewer-${Math.random().toString(36).substr(2, 9)}`);
  const cameraStateRef = useRef(null);

  // Control button styles
  const positionControlButtonStyle = {
    width: isMobile ? '48px' : '40px',
    height: isMobile ? '48px' : '40px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    backgroundColor: 'rgba(56, 189, 248, 0.2)',
    border: '1px solid rgba(56, 189, 248, 0.3)',
    borderRadius: '8px',
    color: '#f0f4f8',
    fontSize: isMobile ? '20px' : '18px',
    cursor: 'pointer',
    padding: '4px',
    touchAction: 'manipulation'
  };

  const rotationControlButtonStyle = {
    width: isMobile ? '46px' : '40px',
    height: isMobile ? '36px' : '30px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    backgroundColor: 'rgba(16, 185, 129, 0.2)',
    border: '1px solid rgba(16, 185, 129, 0.3)',
    borderRadius: '8px',
    color: '#f0f4f8',
    fontSize: isMobile ? '16px' : '14px',
    cursor: 'pointer',
    margin: '2px',
    padding: '4px',
    touchAction: 'manipulation'
  };

  const resetButtonStyle = {
    backgroundColor: 'rgba(15, 23, 42, 0.7)',
    border: '1px solid rgba(255, 255, 255, 0.1)',
    borderRadius: '8px',
    color: '#f0f4f8',
    padding: '8px 12px',
    cursor: 'pointer',
    fontSize: isMobile ? '14px' : '12px',
    touchAction: 'manipulation'
  };

  // Extract atoms for both molecules
  const { molecule1, molecule2 } = atoms || { molecule1: null, molecule2: null };

  // Function to apply rotation transformation to coordinates
  const applyRotation = (coords, rotation, centerOfMass) => {
    // Convert degrees to radians
    const radX = rotation.x * Math.PI / 180;
    const radY = rotation.y * Math.PI / 180;
    const radZ = rotation.z * Math.PI / 180;

    // Translate to origin (center of mass)
    const centered = [
      coords[0] - centerOfMass.x,
      coords[1] - centerOfMass.y,
      coords[2] - centerOfMass.z
    ];

    // Apply rotations (ZYX order)
    // Z-axis rotation
    let nx = centered[0] * Math.cos(radZ) - centered[1] * Math.sin(radZ);
    let ny = centered[0] * Math.sin(radZ) + centered[1] * Math.cos(radZ);
    let nz = centered[2];

    // Y-axis rotation
    let mx = nx * Math.cos(radY) + nz * Math.sin(radY);
    let my = ny;
    let mz = -nx * Math.sin(radY) + nz * Math.cos(radY);

    // X-axis rotation
    let rx = mx;
    let ry = my * Math.cos(radX) - mz * Math.sin(radX);
    let rz = my * Math.sin(radX) + mz * Math.cos(radX);

    // Translate back from origin
    const result = [
      rx + centerOfMass.x,
      ry + centerOfMass.y,
      rz + centerOfMass.z
    ];

    return result;
  };

  const enableViewerInteractions = (viewer) => {
    if (!viewer) return;
    
    try {
      // First ensure current_hover is properly nullified
      if (viewer.current_hover) {
        viewer.current_hover = null;
      }
      
      // Remove overlay
      const viewerElement = viewerRef.current;
      const overlay = viewerElement.querySelector('.event-blocker-overlay');
      if (overlay) {
        viewerElement.removeChild(overlay);
      }
      
      // Reinitialize viewer events AFTER hover state is reset
      viewer.render();
      
      // Enable hover with empty callback to prevent undefined function calls
      const emptyCallback = function() {};
      viewer.setHoverable({}, true, emptyCallback, emptyCallback);
    } catch (e) {
      console.error("Error enabling viewer interactions:", e);
    }
  };
  
  const disableViewerInteractions = (viewer) => {
    if (!viewer) return;
    
    try {
      // Properly disable hover by setting false WITHOUT callbacks
      viewer.setHoverable({}, false);
      
      // Explicitly null the hover state
      viewer.current_hover = null;
      
      // Create event-blocking overlay
      const viewerElement = viewerRef.current;
      let overlay = viewerElement.querySelector('.event-blocker-overlay');
      
      if (!overlay) {
        overlay = document.createElement('div');
        overlay.className = 'event-blocker-overlay';
        overlay.style.position = 'absolute';
        overlay.style.top = '0';
        overlay.style.left = '0';
        overlay.style.width = '100%';
        overlay.style.height = '100%';
        overlay.style.zIndex = '100';
        
        // Attach event handlers to overlay
        overlay.addEventListener('mousemove', e => e.stopPropagation(), true);
        overlay.addEventListener('mousedown', e => e.stopPropagation(), true);
        overlay.addEventListener('mouseup', e => e.stopPropagation(), true);
        
        viewerElement.appendChild(overlay);
      }
    } catch (e) {
      console.error("Error disabling viewer interactions:", e);
    }
  };
    
  // Calculate center of mass for a molecule
  const calculateCenterOfMass = (atoms) => {
    if (!atoms || atoms.length === 0) return { x: 0, y: 0, z: 0 };

    // Simple average of atom coordinates
    const sum = atoms.reduce((acc, atom) => {
      return {
        x: acc.x + atom.x,
        y: acc.y + atom.y,
        z: acc.z + atom.z
      };
    }, { x: 0, y: 0, z: 0 });

    return {
      x: sum.x / atoms.length,
      y: sum.y / atoms.length,
      z: sum.z / atoms.length
    };
  };

  // Calculate covalent bonds between atoms based on distance and element types
  const calculateCovalentBonds = (atoms) => {
    const bonds = [];
    
    // Define covalent radii for common elements (in Angstroms)
    const covalentRadii = {
      'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 
      'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20
    };
    
    // Maximum bond distance factor (multiplier for sum of covalent radii)
    const bondDistanceFactor = 1.3;
    
    // Check all possible atom pairs
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const atom1 = atoms[i];
        const atom2 = atoms[j];
        
        // Calculate distance between atoms
        const distance = Math.sqrt(
          Math.pow(atom1.x - atom2.x, 2) +
          Math.pow(atom1.y - atom2.y, 2) +
          Math.pow(atom1.z - atom2.z, 2)
        );
        
        // Get covalent radii (using default if element not in table)
        const radius1 = covalentRadii[atom1.element] || 0.75;
        const radius2 = covalentRadii[atom2.element] || 0.75;
        
        // Calculate bond threshold as sum of covalent radii with tolerance
        const bondThreshold = (radius1 + radius2) * bondDistanceFactor;
        
        // Check if distance is within bond threshold
        if (distance <= bondThreshold) {
          bonds.push({
            atom1Index: i,
            atom2Index: j,
            distance: distance
          });
        }
      }
    }
    
    return bonds;
  };

  // Calculate hydrogen bonds between donor-H and acceptor atoms
  const calculateHydrogenBonds = (atoms) => {
    const hBonds = [];
    const hBondLengthMax = 3.2; // Maximum H-bond length in Angstroms
    const hBondLengthMin = 1.5; // Minimum length to avoid counting covalent bonds
    
    // Identify potential hydrogen bond donors (H attached to O, N)
    const donors = [];
    const acceptors = [];
    
    // Find all H atoms that are bonded to O or N (potential donors)
    for (let i = 0; i < atoms.length; i++) {
      if (atoms[i].element === 'H') {
        // Check if H is bonded to O or N
        for (let j = 0; j < atoms.length; j++) {
          if (i !== j && (atoms[j].element === 'O' || atoms[j].element === 'N')) {
            const distance = Math.sqrt(
              Math.pow(atoms[i].x - atoms[j].x, 2) +
              Math.pow(atoms[i].y - atoms[j].y, 2) +
              Math.pow(atoms[i].z - atoms[j].z, 2)
            );
            
            // If H is close enough to O/N, it's a potential donor
            if (distance < 1.2) { // Typical H-O/H-N bond length
              donors.push({
                hIndex: i,
                donorHeavyIndex: j // The O or N atom index
              });
            }
          }
        }
      }
    }
    
    // Find all potential acceptors (O or N atoms with lone pairs)
    for (let i = 0; i < atoms.length; i++) {
      if (atoms[i].element === 'O' || atoms[i].element === 'N') {
        acceptors.push(i);
      }
    }
    
    // Check for hydrogen bonds between donors and acceptors
    for (const donor of donors) {
      for (const acceptorIndex of acceptors) {
        // Skip if acceptor is the donor atom itself
        if (acceptorIndex === donor.donorHeavyIndex) continue;
        
        const hAtom = atoms[donor.hIndex];
        const acceptorAtom = atoms[acceptorIndex];
        
        // Calculate H to acceptor distance
        const distance = Math.sqrt(
          Math.pow(hAtom.x - acceptorAtom.x, 2) +
          Math.pow(hAtom.y - acceptorAtom.y, 2) +
          Math.pow(hAtom.z - acceptorAtom.z, 2)
        );
        
        // Check if distance is within hydrogen bond range
        if (distance >= hBondLengthMin && distance <= hBondLengthMax) {
          // Calculate angles to check for linearity
          const donorAtom = atoms[donor.donorHeavyIndex];
          
          // Calculate vectors for angle determination
          const donorToH = {
            x: hAtom.x - donorAtom.x,
            y: hAtom.y - donorAtom.y,
            z: hAtom.z - donorAtom.z
          };
          
          const hToAcceptor = {
            x: acceptorAtom.x - hAtom.x,
            y: acceptorAtom.y - hAtom.y,
            z: acceptorAtom.z - hAtom.z
          };
          
          // Normalize vectors
          const donorToHMag = Math.sqrt(
            donorToH.x * donorToH.x + donorToH.y * donorToH.y + donorToH.z * donorToH.z
          );
          const hToAcceptorMag = Math.sqrt(
            hToAcceptor.x * hToAcceptor.x + hToAcceptor.y * hToAcceptor.y + hToAcceptor.z * hToAcceptor.z
          );
          
          const donorToHNorm = {
            x: donorToH.x / donorToHMag,
            y: donorToH.y / donorToHMag,
            z: donorToH.z / donorToHMag
          };
          
          const hToAcceptorNorm = {
            x: hToAcceptor.x / hToAcceptorMag,
            y: hToAcceptor.y / hToAcceptorMag,
            z: hToAcceptor.z / hToAcceptorMag
          };
          
          // Calculate dot product to get cosine of angle
          const dotProduct = 
            donorToHNorm.x * hToAcceptorNorm.x + 
            donorToHNorm.y * hToAcceptorNorm.y + 
            donorToHNorm.z * hToAcceptorNorm.z;
          
          // For hydrogen bond, angle should be reasonably linear (cos > 0.7 is < 45°)
          if (dotProduct > 0.7) {
            hBonds.push({
              hIndex: donor.hIndex,
              donorHeavyIndex: donor.donorHeavyIndex,
              acceptorIndex: acceptorIndex,
              distance: distance
            });
          }
        }
      }
    }
    
    return hBonds;
  };

  /**
   * Calculate hydrogen bonds between two different molecules
   * @param {Array} molecule1 - First molecule atoms
   * @param {Array} molecule2 - Second molecule atoms (with transformations applied)
   * @returns {Array} Array of intermolecular hydrogen bonds
   */
  const calculateIntermolecularHydrogenBonds = (molecule1, molecule2) => {
    const hBonds = [];
    const hBondLengthMax = 3.2; // Maximum H-bond length in Angstroms
    const hBondLengthMin = 1.5; // Minimum length to avoid counting covalent bonds
    
    // Identify donors and acceptors in molecule 1
    const donors1 = [];
    const acceptors1 = [];
    
    // Find all H atoms in molecule1 that are bonded to O or N (potential donors)
    for (let i = 0; i < molecule1.length; i++) {
      if (molecule1[i].element === 'H') {
        // Check if H is bonded to O or N
        for (let j = 0; j < molecule1.length; j++) {
          if (i !== j && (molecule1[j].element === 'O' || molecule1[j].element === 'N')) {
            const distance = Math.sqrt(
              Math.pow(molecule1[i].x - molecule1[j].x, 2) +
              Math.pow(molecule1[i].y - molecule1[j].y, 2) +
              Math.pow(molecule1[i].z - molecule1[j].z, 2)
            );
            
            // If H is close enough to O/N, it's a potential donor
            if (distance < 1.2) { // Typical H-O/H-N bond length
              donors1.push({
                hIndex: i,
                donorHeavyIndex: j,
                molecule: 1 // Mark as belonging to molecule 1
              });
            }
          }
        }
      }
    }
    
    // Find all potential acceptors in molecule 1 (O or N atoms)
    for (let i = 0; i < molecule1.length; i++) {
      if (molecule1[i].element === 'O' || molecule1[i].element === 'N') {
        acceptors1.push({
          index: i,
          molecule: 1 // Mark as belonging to molecule 1
        });
      }
    }
    
    // Identify donors and acceptors in molecule 2
    const donors2 = [];
    const acceptors2 = [];
    
    // Find all H atoms in molecule2 that are bonded to O or N (potential donors)
    for (let i = 0; i < molecule2.length; i++) {
      if (molecule2[i].element === 'H') {
        // Check if H is bonded to O or N
        for (let j = 0; j < molecule2.length; j++) {
          if (i !== j && (molecule2[j].element === 'O' || molecule2[j].element === 'N')) {
            const distance = Math.sqrt(
              Math.pow(molecule2[i].x - molecule2[j].x, 2) +
              Math.pow(molecule2[i].y - molecule2[j].y, 2) +
              Math.pow(molecule2[i].z - molecule2[j].z, 2)
            );
            
            // If H is close enough to O/N, it's a potential donor
            if (distance < 1.2) { // Typical H-O/H-N bond length
              donors2.push({
                hIndex: i,
                donorHeavyIndex: j,
                molecule: 2 // Mark as belonging to molecule 2
              });
            }
          }
        }
      }
    }
    
    // Find all potential acceptors in molecule 2 (O or N atoms)
    for (let i = 0; i < molecule2.length; i++) {
      if (molecule2[i].element === 'O' || molecule2[i].element === 'N') {
        acceptors2.push({
          index: i,
          molecule: 2 // Mark as belonging to molecule 2
        });
      }
    }
    
    // Check for hydrogen bonds from molecule 1 donors to molecule 2 acceptors
    for (const donor of donors1) {
      for (const acceptor of acceptors2) {
        const hAtom = molecule1[donor.hIndex];
        const acceptorAtom = molecule2[acceptor.index];
        
        // Calculate H to acceptor distance
        const distance = Math.sqrt(
          Math.pow(hAtom.x - acceptorAtom.x, 2) +
          Math.pow(hAtom.y - acceptorAtom.y, 2) +
          Math.pow(hAtom.z - acceptorAtom.z, 2)
        );
        
        // Check if distance is within hydrogen bond range
        if (distance >= hBondLengthMin && distance <= hBondLengthMax) {
          // Calculate angles to check for linearity
          const donorAtom = molecule1[donor.donorHeavyIndex];
          
          // Calculate vectors for angle determination
          const donorToH = {
            x: hAtom.x - donorAtom.x,
            y: hAtom.y - donorAtom.y,
            z: hAtom.z - donorAtom.z
          };
          
          const hToAcceptor = {
            x: acceptorAtom.x - hAtom.x,
            y: acceptorAtom.y - hAtom.y,
            z: acceptorAtom.z - hAtom.z
          };
          
          // Normalize vectors
          const donorToHMag = Math.sqrt(
            donorToH.x * donorToH.x + donorToH.y * donorToH.y + donorToH.z * donorToH.z
          );
          const hToAcceptorMag = Math.sqrt(
            hToAcceptor.x * hToAcceptor.x + hToAcceptor.y * hToAcceptor.y + hToAcceptor.z * hToAcceptor.z
          );
          
          const donorToHNorm = {
            x: donorToH.x / donorToHMag,
            y: donorToH.y / donorToHMag,
            z: donorToH.z / donorToHMag
          };
          
          const hToAcceptorNorm = {
            x: hToAcceptor.x / hToAcceptorMag,
            y: hToAcceptor.y / hToAcceptorMag,
            z: hToAcceptor.z / hToAcceptorMag
          };
          
          // Calculate dot product to get cosine of angle
          const dotProduct = 
            donorToHNorm.x * hToAcceptorNorm.x + 
            donorToHNorm.y * hToAcceptorNorm.y + 
            donorToHNorm.z * hToAcceptorNorm.z;
          
          // For hydrogen bond, angle should be reasonably linear (cos > 0.7 is < 45°)
          if (dotProduct > 0.7) {
            hBonds.push({
              donor: {
                hIndex: donor.hIndex,
                donorHeavyIndex: donor.donorHeavyIndex,
                molecule: 1
              },
              acceptor: {
                index: acceptor.index,
                molecule: 2
              },
              distance: distance
            });
          }
        }
      }
    }
    
    // Check for hydrogen bonds from molecule 2 donors to molecule 1 acceptors
    for (const donor of donors2) {
      for (const acceptor of acceptors1) {
        const hAtom = molecule2[donor.hIndex];
        const acceptorAtom = molecule1[acceptor.index];
        
        // Calculate H to acceptor distance
        const distance = Math.sqrt(
          Math.pow(hAtom.x - acceptorAtom.x, 2) +
          Math.pow(hAtom.y - acceptorAtom.y, 2) +
          Math.pow(hAtom.z - acceptorAtom.z, 2)
        );
        
        // Check if distance is within hydrogen bond range
        if (distance >= hBondLengthMin && distance <= hBondLengthMax) {
          // Calculate angles to check for linearity
          const donorAtom = molecule2[donor.donorHeavyIndex];
          
          // Calculate vectors for angle determination
          const donorToH = {
            x: hAtom.x - donorAtom.x,
            y: hAtom.y - donorAtom.y,
            z: hAtom.z - donorAtom.z
          };
          
          const hToAcceptor = {
            x: acceptorAtom.x - hAtom.x,
            y: acceptorAtom.y - hAtom.y,
            z: acceptorAtom.z - hAtom.z
          };
          
          // Normalize vectors
          const donorToHMag = Math.sqrt(
            donorToH.x * donorToH.x + donorToH.y * donorToH.y + donorToH.z * donorToH.z
          );
          const hToAcceptorMag = Math.sqrt(
            hToAcceptor.x * hToAcceptor.x + hToAcceptor.y * hToAcceptor.y + hToAcceptor.z * hToAcceptor.z
          );
          
          const donorToHNorm = {
            x: donorToH.x / donorToHMag,
            y: donorToH.y / donorToHMag,
            z: donorToH.z / donorToHMag
          };
          
          const hToAcceptorNorm = {
            x: hToAcceptor.x / hToAcceptorMag,
            y: hToAcceptor.y / hToAcceptorMag,
            z: hToAcceptor.z / hToAcceptorMag
          };
          
          // Calculate dot product to get cosine of angle
          const dotProduct = 
            donorToHNorm.x * hToAcceptorNorm.x + 
            donorToHNorm.y * hToAcceptorNorm.y + 
            donorToHNorm.z * hToAcceptorNorm.z;
          
          // For hydrogen bond, angle should be reasonably linear (cos > 0.7 is < 45°)
          if (dotProduct > 0.7) {
            hBonds.push({
              donor: {
                hIndex: donor.hIndex,
                donorHeavyIndex: donor.donorHeavyIndex,
                molecule: 2
              },
              acceptor: {
                index: acceptor.index,
                molecule: 1
              },
              distance: distance
            });
          }
        }
      }
    }
    
    return hBonds;
  };

  // Function to store current camera state
  const saveCameraState = (viewer) => {
    if (!viewer) return null;

    try {
      // Extract camera parameters from the viewer
      const camState = {
        position: viewer.getView(),
        model: viewer.getModel()
      };
      return camState;
    } catch (e) {
      console.error("Error saving camera state:", e);
      return null;
    }
  };

  // Function to restore camera state
  const restoreCameraState = (viewer, state) => {
    if (!viewer || !state) return;

    try {
      // Restore camera parameters
      viewer.setView(state.position);
      viewer.render();
    } catch (e) {
      console.error("Error restoring camera state:", e);
    }
  };

  // Initial viewer setup - only runs once
  useEffect(() => {
    // Skip if no molecules are provided
    if (!molecule1 && !molecule2) {
      console.warn("No molecule data found for visualization.");
      return;
    }

    // Set up viewer instance only once
    const viewer = $3Dmol.createViewer(viewerRef.current, {
      backgroundColor: "rgb(15, 23, 42)",
    });

    // Store viewer instance for later use
    setViewerInstance(viewer);
    setIsInitialRender(true);

    // Cleanup function
    return () => {
      if (viewer) {
        try {
          viewer.clear();
        } catch (e) {
          console.error("Error cleaning up viewer:", e);
        }
      }
    };
  }, []); // Empty dependency array means this runs once

  // Effect for updating molecule visualization based on changes
  useEffect(() => {
    // Skip if no molecules are provided or viewer not initialized
    if ((!molecule1 && !molecule2) || !viewerInstance) {
      return;
    }

    // Add a slight delay to ensure container is fully rendered
    const timer = setTimeout(() => {
      try {
        // Save current camera state before updating
        if (!isInitialRender) {
          cameraStateRef.current = saveCameraState(viewerInstance);
        }

        // Clear previous content
        viewerInstance.clear();
        
        // Initialize transformed molecule2 variable to null at the function scope level
        // This ensures it's available throughout the entire rendering logic
        let transformedMolecule2 = null;

        // Handle molecule1
        if (molecule1) {
          // Transform molecule1 data to 3Dmol format
          let m1Data = molecule1.map(atom => ({
            elem: atom.element,
            x: atom.x,
            y: atom.y,
            z: atom.z,
            properties: { molecule: 1 }
          }));

          // Add molecule1 with its own model
          let model1 = viewerInstance.addModel();
          model1.addAtoms(m1Data);
          
          // Calculate covalent bonds
          const covalentBonds1 = calculateCovalentBonds(molecule1);
          
          // Add covalent bonds as cylinders instead of using addBond
          covalentBonds1.forEach(bond => {
            const atom1 = molecule1[bond.atom1Index];
            const atom2 = molecule1[bond.atom2Index];
            
            viewerInstance.addCylinder({
              start: { x: atom1.x, y: atom1.y, z: atom1.z },
              end: { x: atom2.x, y: atom2.y, z: atom2.z },
              radius: isMobile ? 0.12 : 0.15,
              fromCap: 1,
              toCap: 1,
              color: "0x38bdf8"  // Blue for molecule 1
            });
          });
          
          // Calculate and visualize hydrogen bonds
          const hydrogenBonds1 = calculateHydrogenBonds(molecule1);
          
          // Add hydrogen bonds visualization with dashed lines
          hydrogenBonds1.forEach(bond => {
            viewerInstance.addCylinder({
              start: {
                x: molecule1[bond.hIndex].x,
                y: molecule1[bond.hIndex].y,
                z: molecule1[bond.hIndex].z
              },
              end: {
                x: molecule1[bond.acceptorIndex].x,
                y: molecule1[bond.acceptorIndex].y,
                z: molecule1[bond.acceptorIndex].z
              },
              radius: isMobile ? 0.05 : 0.07,  // Thinner than covalent bonds
              fromCap: 1,
              toCap: 1,
              color: "0xFFFFFF",  // White color for hydrogen bonds
              dashed: true,
              dashLength: 0.15,   // Length of dash segments
              gapLength: 0.15,    // Length of gaps
            });
          });

          // Add element labels for molecule 1
          molecule1.forEach((atom) => {
            viewerInstance.addLabel(atom.element, {
              position: { x: atom.x, y: atom.y, z: atom.z },
              fontSize: isMobile ? 10 : 12,
              fontColor: "white",
              backgroundColor: "rgba(56, 189, 248, 0.5)",  // Blue background for molecule 1
              borderRadius: 10,
              padding: isMobile ? 1 : 2,
              inFront: true,
            });
          });

          // Add molecule1 label if in interaction mode
          if (molecule1 && molecule2) {
            viewerInstance.addLabel("Molecule 1 (Fixed)", {
              position: { x: molecule1[0].x, y: molecule1[0].y, z: molecule1[0].z + 5 },
              fontSize: isMobile ? 14 : 16,
              fontColor: "white",
              backgroundColor: "rgba(56, 189, 248, 0.7)",
              borderRadius: 10,
              padding: isMobile ? 2 : 4,
              inFront: true,
              fixedPosition: true,
            });
          }
        }

        // Handle molecule2
        if (molecule2) {
          // Calculate center of mass for rotation
          const centerOfMass = calculateCenterOfMass(molecule2);

          // Process each atom with rotation and offset
          let m2Data = molecule2.map(atom => {
            // Apply rotation if needed
            let x = atom.x, y = atom.y, z = atom.z;

            if (molecule2Rotation && (molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0)) {
              const coords = [x, y, z];
              const rotated = applyRotation(coords, molecule2Rotation, centerOfMass);
              x = rotated[0];
              y = rotated[1];
              z = rotated[2];
            }

            // Apply offset and return atom
            return {
              elem: atom.element,
              x: x + molecule2Offset.x,
              y: y + molecule2Offset.y,
              z: z + molecule2Offset.z,
              properties: { molecule: 2 }
            };
          });

          // Add molecule2 with its own model
          let model2 = viewerInstance.addModel();
          model2.addAtoms(m2Data);

          // Get properly transformed coordinates with rotation and offset applied
          const transformedMolecule2 = molecule2.map(atom => {
            // Apply transformation as in the existing code
            let x = atom.x, y = atom.y, z = atom.z;
            
            if (molecule2Rotation && (molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0)) {
              const coords = [x, y, z];
              const rotated = applyRotation(coords, molecule2Rotation, centerOfMass);
              x = rotated[0];
              y = rotated[1];
              z = rotated[2];
            }
            
            return {
              element: atom.element,
              x: x + molecule2Offset.x,
              y: y + molecule2Offset.y,
              z: z + molecule2Offset.z
            };
          });
          
          // Calculate bonds for transformed coordinates
          const covalentBonds2 = calculateCovalentBonds(transformedMolecule2);
          
          // Add covalent bonds using cylinders instead of addBond
          covalentBonds2.forEach(bond => {
            const atom1 = transformedMolecule2[bond.atom1Index];
            const atom2 = transformedMolecule2[bond.atom2Index];
            
            viewerInstance.addCylinder({
              start: { x: atom1.x, y: atom1.y, z: atom1.z },
              end: { x: atom2.x, y: atom2.y, z: atom2.z },
              radius: isMobile ? 0.12 : 0.15,
              fromCap: 1,
              toCap: 1,
              color: "0x10b981"  // Green for molecule 2
            });
          });
          
          // Calculate hydrogen bonds
          const hydrogenBonds2 = calculateHydrogenBonds(transformedMolecule2);
          
          // Add hydrogen bond visualization
          hydrogenBonds2.forEach(bond => {
            viewerInstance.addCylinder({
              start: {
                x: transformedMolecule2[bond.hIndex].x,
                y: transformedMolecule2[bond.hIndex].y,
                z: transformedMolecule2[bond.hIndex].z
              },
              end: {
                x: transformedMolecule2[bond.acceptorIndex].x,
                y: transformedMolecule2[bond.acceptorIndex].y,
                z: transformedMolecule2[bond.acceptorIndex].z
              },
              radius: isMobile ? 0.05 : 0.07,  
              fromCap: 1,
              toCap: 1,
              color: "0xFFFFFF",  
              dashed: true,
              dashLength: 0.15,   
              gapLength: 0.15,   
            });
          });

          // Add element labels for molecule2
          molecule2.forEach((atom) => {
            // Apply same transformations as for the atoms
            let x = atom.x, y = atom.y, z = atom.z;

            if (molecule2Rotation && (molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0)) {
              const coords = [x, y, z];
              const rotated = applyRotation(coords, molecule2Rotation, centerOfMass);
              x = rotated[0];
              y = rotated[1];
              z = rotated[2];
            }

            viewerInstance.addLabel(atom.element, {
              position: {
                x: x + molecule2Offset.x,
                y: y + molecule2Offset.y,
                z: z + molecule2Offset.z
              },
              fontSize: isMobile ? 10 : 12,
              fontColor: "white",
              backgroundColor: "rgba(16, 185, 129, 0.5)",  // Green background
              borderRadius: 10,
              padding: isMobile ? 1 : 2,
              inFront: true,
            });
          });

          // Add molecule2 label in interaction mode
          if (molecule1 && molecule2) {
            // Calculate position for label with transformations
            let labelX = molecule2[0].x, labelY = molecule2[0].y, labelZ = molecule2[0].z;

            if (molecule2Rotation && (molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0)) {
              const coords = [labelX, labelY, labelZ];
              const rotated = applyRotation(coords, molecule2Rotation, centerOfMass);
              labelX = rotated[0];
              labelY = rotated[1];
              labelZ = rotated[2];
            }

            viewerInstance.addLabel("Molecule 2 (Use Arrow Keys)", {
              position: {
                x: labelX + molecule2Offset.x,
                y: labelY + molecule2Offset.y,
                z: labelZ + molecule2Offset.z + 5
              },
              fontSize: isMobile ? 14 : 16,
              fontColor: "white",
              backgroundColor: "rgba(16, 185, 129, 0.7)",
              borderRadius: 10,
              padding: isMobile ? 2 : 4,
              inFront: true,
              fixedPosition: true,
            });
          }
        }

        // Calculate and visualize intermolecular hydrogen bonds
        // Move this section AFTER both molecule1 and molecule2 have been fully processed
        // This ensures transformedMolecule2 is properly initialized before use
        if (molecule1 && molecule2 && transformedMolecule2 && Array.isArray(molecule1) && Array.isArray(transformedMolecule2)) {
          console.log("Calculating intermolecular hydrogen bonds with:", 
                     "molecule1.length:", molecule1.length, 
                     "transformedMolecule2.length:", transformedMolecule2.length);
          // Calculate intermolecular hydrogen bonds
          const intermolecularHBonds = calculateIntermolecularHydrogenBonds(molecule1, transformedMolecule2);
          
          // Add intermolecular hydrogen bond visualization with dashed lines
          intermolecularHBonds.forEach(bond => {
            // Get start coordinates (hydrogen atom)
            let startX, startY, startZ;
            if (bond.donor.molecule === 1) {
              startX = molecule1[bond.donor.hIndex].x;
              startY = molecule1[bond.donor.hIndex].y;
              startZ = molecule1[bond.donor.hIndex].z;
            } else {
              startX = transformedMolecule2[bond.donor.hIndex].x;
              startY = transformedMolecule2[bond.donor.hIndex].y;
              startZ = transformedMolecule2[bond.donor.hIndex].z;
            }
            
            // Get end coordinates (acceptor atom)
            let endX, endY, endZ;
            if (bond.acceptor.molecule === 1) {
              endX = molecule1[bond.acceptor.index].x;
              endY = molecule1[bond.acceptor.index].y;
              endZ = molecule1[bond.acceptor.index].z;
            } else {
              endX = transformedMolecule2[bond.acceptor.index].x;
              endY = transformedMolecule2[bond.acceptor.index].y;
              endZ = transformedMolecule2[bond.acceptor.index].z;
            }
            
            // Add dashed line for hydrogen bond
            viewerInstance.addCylinder({
              start: { x: startX, y: startY, z: startZ },
              end: { x: endX, y: endY, z: endZ },
              radius: isMobile ? 0.05 : 0.07,  // Thinner than covalent bonds
              fromCap: 1,
              toCap: 1,
              color: "0xFFD700",  // Gold color for intermolecular hydrogen bonds (for distinction)
              dashed: true,
              dashLength: 0.15,
              gapLength: 0.15,
            });
            
            // Add a label for the hydrogen bond if needed
            if (bond.distance <= 2.5) { // Only label stronger H-bonds
              viewerInstance.addLabel(`${bond.distance.toFixed(2)}Å`, {
                position: {
                  x: (startX + endX) / 2,
                  y: (startY + endY) / 2,
                  z: (startZ + endZ) / 2
                },
                fontSize: isMobile ? 8 : 10,
                fontColor: "white",
                backgroundColor: "rgba(255, 215, 0, 0.5)",
                borderRadius: 10,
                padding: 2,
                inFront: true,
              });
            }
          });
        }

        // No longer need 3D-based legend as we're using a static HTML overlay instead

        // Style both molecules - must be done after adding all models
        viewerInstance.setStyle({ properties: { molecule: 1 } }, {
          sphere: {
            radius: isMobile ? 0.30 : 0.35,
            scale: isMobile ? 0.85 : 0.9,
            color: "0x38bdf8"  // Blue for molecule 1
          },
          stick: {
            radius: isMobile ? 0.12 : 0.15,
            color: "0x38bdf8",
            smooth: true
          },
        });

        viewerInstance.setStyle({ properties: { molecule: 2 } }, {
          sphere: {
            radius: isMobile ? 0.30 : 0.35,
            scale: isMobile ? 0.85 : 0.9,
            color: "0x10b981"  // Green for molecule 2
          },
          stick: {
            radius: isMobile ? 0.12 : 0.15,
            color: "0x10b981",
            smooth: true
          },
        });

        // Handle mouse interaction behavior and positioning mode
        try {
          if (positioningMode) {
            // Use the dedicated function to disable all interactions properly
            disableViewerInteractions(viewerInstance);
            
            // Add positioning mode instructions if needed
            if (molecule2) {
              viewerInstance.addLabel("Use arrow keys to position Molecule 2", {
                position: { x: 0, y: 0, z: 10 },
                fontSize: 16,
                fontColor: "white",
                backgroundColor: "rgba(0, 0, 0, 0.7)",
                borderRadius: 10,
                padding: 5,
                inFront: true,
                fixedPosition: true,
              });
            }
          } else {
            // Reset to normal behavior
            enableViewerInteractions(viewerInstance);
          }
        } catch (e) {
          console.error("3DMol interaction handler error:", e);
        }

        // Initial view or restore camera
        if (isInitialRender) {
          viewerInstance.zoomTo();
          setIsInitialRender(false);
        } else if (cameraStateRef.current) {
          restoreCameraState(viewerInstance, cameraStateRef.current);
        }

        // Force render and resize
        viewerInstance.render();
        viewerInstance.resize();
      } catch (error) {
        console.error("Error rendering molecule:", error);
      }
    }, 50);

    // Set up responsive resize handler
    const handleResize = () => {
      if (viewerRef.current) {
        const viewer = $3Dmol.viewers[viewerRef.current.id];
        if (viewer) {
          const tempCameraState = saveCameraState(viewer);
          viewer.resize();
          if (tempCameraState) {
            restoreCameraState(viewer, tempCameraState);
          } else {
            viewer.render();
          }
        }
      }
    };

    window.addEventListener('resize', handleResize);

    // Cleanup
    return () => {
      window.removeEventListener('resize', handleResize);
      clearTimeout(timer);
    };
  }, [atoms, isMobile, molecule2Offset, molecule2Rotation, positioningMode, viewerInstance, isInitialRender]);

  // Effect to handle keyboard events for positioning mode
  useEffect(() => {
    if (!positioningMode || !molecule2) return;

    const container = containerRef.current;
    container.focus();

    const handleKeyDown = (e) => {
      const moveStep = 0.5; // Angstroms per keypress
      const rotationStep = 15; // Degrees per keypress

      if (e.shiftKey) {
        // Shift key held - handle rotation
        let newRotation = { ...molecule2Rotation };

        switch (e.key) {
          case 'ArrowLeft':
            newRotation.y = ((newRotation.y - rotationStep) % 360 + 360) % 360;
            break;
          case 'ArrowRight':
            newRotation.y = ((newRotation.y + rotationStep) % 360 + 360) % 360;
            break;
          case 'ArrowUp':
            newRotation.x = ((newRotation.x - rotationStep) % 360 + 360) % 360;
            break;
          case 'ArrowDown':
            newRotation.x = ((newRotation.x + rotationStep) % 360 + 360) % 360;
            break;
          case 'PageUp':
            newRotation.z = ((newRotation.z + rotationStep) % 360 + 360) % 360;
            break;
          case 'PageDown':
            newRotation.z = ((newRotation.z - rotationStep) % 360 + 360) % 360;
            break;
          default:
            return; // Exit if not a handled key
        }

        onMoleculeRotate(newRotation);
        e.preventDefault();
      } else {
        // No shift key - handle translation
        let newOffset = { ...molecule2Offset };

        switch (e.key) {
          case 'ArrowLeft':
            newOffset.x -= moveStep;
            break;
          case 'ArrowRight':
            newOffset.x += moveStep;
            break;
          case 'ArrowUp':
            newOffset.y += moveStep;
            break;
          case 'ArrowDown':
            newOffset.y -= moveStep;
            break;
          case 'PageUp':
            newOffset.z += moveStep;
            break;
          case 'PageDown':
            newOffset.z -= moveStep;
            break;
          default:
            return; // Exit if not a handled key
        }

        onMoleculeMove(newOffset);
        e.preventDefault();
      }
    };

    container.addEventListener('keydown', handleKeyDown, true);

    return () => {
      container.removeEventListener('keydown', handleKeyDown, true);
    };
  }, [positioningMode, molecule2, molecule2Offset, molecule2Rotation, onMoleculeMove, onMoleculeRotate]);

  return (
    <>
      {positioningMode && molecule2 && (
        <div style={{
          display: 'flex',
          flexDirection: 'column',
          gap: '10px',
          marginBottom: '16px',
          backgroundColor: 'rgba(15, 23, 42, 0.7)',
          padding: isMobile ? '12px' : '16px',
          borderRadius: '8px',
          border: '1px solid rgba(255, 255, 255, 0.1)'
        }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <h4 style={{ margin: '0 0 8px 0', fontSize: '16px' }}>Position Controls</h4>
          </div>
          
          {/* XY position control buttons */}
          <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: '4px' }}>
            <button
              onClick={() => onMoleculeMove({ ...molecule2Offset, y: molecule2Offset.y + 0.5 })}
              style={positionControlButtonStyle}>
              ↑
            </button>
            <div style={{ display: 'flex', gap: '4px' }}>
              <button
                onClick={() => onMoleculeMove({ ...molecule2Offset, x: molecule2Offset.x - 0.5 })}
                style={positionControlButtonStyle}>
                ←
              </button>
              <button
                onClick={() => onMoleculeMove({ ...molecule2Offset, y: molecule2Offset.y - 0.5 })}
                style={positionControlButtonStyle}>
                ↓
              </button>
              <button
                onClick={() => onMoleculeMove({ ...molecule2Offset, x: molecule2Offset.x + 0.5 })}
                style={positionControlButtonStyle}>
                →
              </button>
            </div>
          </div>
          
          {/* Z position control buttons */}
          <div style={{ display: 'flex', justifyContent: 'center', gap: '4px', marginTop: '4px' }}>
            <button
              onClick={() => onMoleculeMove({ ...molecule2Offset, z: molecule2Offset.z - 0.5 })}
              style={positionControlButtonStyle}>
              Z−
            </button>
            <button
              onClick={() => onMoleculeMove({ ...molecule2Offset, z: molecule2Offset.z + 0.5 })}
              style={positionControlButtonStyle}>
              Z+
            </button>
          </div>
          
          {/* Rotation control buttons */}
          <div style={{ marginTop: '8px', padding: '8px 0', borderTop: '1px solid rgba(255, 255, 255, 0.1)' }}>
            <h4 style={{ margin: '0 0 8px 0', fontSize: '16px' }}>Rotation Controls</h4>
            <div style={{ display: 'flex', justifyContent: 'center', gap: '8px' }}>
              <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
                <button
                  onClick={() => onMoleculeRotate({ ...molecule2Rotation, x: (molecule2Rotation.x - 15 + 360) % 360 })}
                  style={rotationControlButtonStyle}>
                  X−
                </button>
                <button
                  onClick={() => onMoleculeRotate({ ...molecule2Rotation, x: (molecule2Rotation.x + 15) % 360 })}
                  style={rotationControlButtonStyle}>
                  X+
                </button>
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
                <button
                  onClick={() => onMoleculeRotate({ ...molecule2Rotation, y: (molecule2Rotation.y - 15 + 360) % 360 })}
                  style={rotationControlButtonStyle}>
                  Y−
                </button>
                <button
                  onClick={() => onMoleculeRotate({ ...molecule2Rotation, y: (molecule2Rotation.y + 15) % 360 })}
                  style={rotationControlButtonStyle}>
                  Y+
                </button>
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
                <button
                  onClick={() => onMoleculeRotate({ ...molecule2Rotation, z: (molecule2Rotation.z - 15 + 360) % 360 })}
                  style={rotationControlButtonStyle}>
                  Z−
                </button>
                <button
                  onClick={() => onMoleculeRotate({ ...molecule2Rotation, z: (molecule2Rotation.z + 15) % 360 })}
                  style={rotationControlButtonStyle}>
                  Z+
                </button>
              </div>
            </div>
          </div>
          
          {/* Reset buttons */}
          <div style={{ 
            display: 'flex', 
            justifyContent: 'center', 
            gap: '8px',
            marginTop: '8px', 
            paddingTop: '8px', 
            borderTop: '1px solid rgba(255, 255, 255, 0.1)' 
          }}>
            <button
              onClick={() => onMoleculeMove({ x: 0, y: 0, z: 0 })}
              style={resetButtonStyle}>
              Reset Position
            </button>
            <button
              onClick={() => onMoleculeRotate({ x: 0, y: 0, z: 0 })}
              style={resetButtonStyle}>
              Reset Rotation
            </button>
          </div>
        </div>
      )}

      <div
        ref={containerRef}
        className="viewer-container"
        tabIndex="0"
        style={{
          position: 'relative',
          width: '100%',
          height: isMobile ? '350px' : '450px',
          margin: '0 auto',
          borderRadius: '12px',
          overflow: 'hidden',
          backgroundColor: "rgba(15, 23, 42, 0.5)",
          boxShadow: "0 10px 15px -3px rgba(0, 0, 0, 0.2)",
          border: "1px solid rgba(255, 255, 255, 0.1)",
          cursor: "auto",
          outline: "none"
        }}
      >
        <div
          ref={viewerRef}
          id={viewerIdRef.current}
          style={{
            width: "100%",
            height: "100%",
            position: "absolute",
            top: 0,
            left: 0,
            borderRadius: "12px",
            pointerEvents: "auto"
          }}
        ></div>

        {/* Static Bond Legend Overlay */}
        {(molecule1 || molecule2) && (
          <div style={{
            position: 'absolute',
            top: '10px',
            left: '10px',
            backgroundColor: 'rgba(0,0,0,0.7)',
            borderRadius: '5px',
            padding: '8px 10px',
            color: 'white',
            fontSize: '12px',
            pointerEvents: 'none',
            zIndex: 20,
            maxWidth: isMobile ? '160px' : '180px'
          }}>
            <div style={{ fontWeight: 'bold', marginBottom: '6px' }}>Bond Types:</div>
            
            {/* Covalent Bond */}
            <div style={{ display: 'flex', alignItems: 'center', marginBottom: '6px' }}>
              <div style={{ 
                width: '30px', 
                height: isMobile ? '4px' : '5px', 
                backgroundColor: molecule1 ? '#38bdf8' : '#10b981',
                borderRadius: '4px',
                marginRight: '8px'
              }}></div>
              <div>Covalent</div>
            </div>
            
            {/* Hydrogen Bond */}
            <div style={{ display: 'flex', alignItems: 'center', marginBottom: '6px' }}>
              <div style={{ 
                width: '30px', 
                height: '0', 
                borderTop: isMobile ? '2px dashed white' : '3px dashed white',
                marginRight: '8px'
              }}></div>
              <div>Hydrogen</div>
            </div>
            
            {/* Intermolecular Hydrogen Bond - only show if both molecules present */}
            {molecule1 && molecule2 && (
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <div style={{ 
                  width: '30px', 
                  height: '0', 
                  borderTop: isMobile ? '2px dashed #FFD700' : '3px dashed #FFD700',
                  marginRight: '8px'
                }}></div>
                <div>Intermolecular H-bond</div>
              </div>
            )}
          </div>
        )}

        {positioningMode && molecule2 && (
          <div style={{
            position: 'absolute',
            bottom: '10px',
            right: '10px',
            backgroundColor: 'rgba(0,0,0,0.7)',
            borderRadius: '5px',
            padding: '8px 10px',
            color: 'white',
            fontSize: '12px',
            pointerEvents: 'none',
            zIndex: 20
          }}>
            <div>Offset: X: {molecule2Offset.x.toFixed(2)}, Y: {molecule2Offset.y.toFixed(2)}, Z: {molecule2Offset.z.toFixed(2)}</div>
            <div>Rotation: X: {molecule2Rotation.x}°, Y: {molecule2Rotation.y}°, Z: {molecule2Rotation.z}°</div>
            <div style={{ marginTop: '4px', borderTop: '1px solid rgba(255,255,255,0.2)', paddingTop: '4px' }}>
              <strong>Keyboard Controls:</strong>
              <div>• Arrow keys (←→↑↓): X/Y position</div>
              <div>• PageUp/PageDown: Z position</div>
              <div>• Hold Shift + arrows: Rotation</div>
            </div>
          </div>
        )}
      </div>
    </>
  );
};

export default MoleculeViewer;