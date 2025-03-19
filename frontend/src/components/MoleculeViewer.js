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

        // Style both molecules - must be done after adding all models
        viewerInstance.setStyle({properties: {molecule: 1}}, {
          sphere: { 
            radius: isMobile ? 0.30 : 0.35, 
            scale: isMobile ? 0.85 : 0.9,
            color: "0x38bdf8"  // Blue for molecule 1
          },
          stick: { 
            radius: isMobile ? 0.12 : 0.15,
            color: "0x38bdf8"  
          },
        });
        
        viewerInstance.setStyle({properties: {molecule: 2}}, {
          sphere: { 
            radius: isMobile ? 0.30 : 0.35, 
            scale: isMobile ? 0.85 : 0.9,
            color: "0x10b981"  // Green for molecule 2
          },
          stick: { 
            radius: isMobile ? 0.12 : 0.15,
            color: "0x10b981"  
          },
        });
        
        // Add positioning mode instructions if active
        if (positioningMode && molecule2) {
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
        
        // Disable default mouse handling in positioning mode
        if (positioningMode) {
          viewerInstance.setClickable(false, true);
          viewerInstance.setHoverable(false, true);
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
        let newRotation = {...molecule2Rotation};
        
        switch(e.key) {
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
        let newOffset = {...molecule2Offset};
        
        switch(e.key) {
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
  );
};

export default MoleculeViewer;