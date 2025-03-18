import React, { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';
import { styles } from '../styles/components';

const MoleculeViewer = ({ 
  atoms, 
  isMobile, 
  positioningMode, 
  onMoleculeMove, 
  molecule2Offset 
}) => {
  const viewerRef = useRef();
  const [isDragging, setIsDragging] = useState(false);
  const [dragStartPos, setDragStartPos] = useState({ x: 0, y: 0 });
  const [viewerInstance, setViewerInstance] = useState(null);
  const [modelInstance, setModelInstance] = useState(null);
  
  // Extract atoms for both molecules
  const { molecule1, molecule2 } = atoms || { molecule1: null, molecule2: null };

  // Effect for initial molecule rendering
  useEffect(() => {
    // Skip if no molecules are provided
    if (!molecule1 && !molecule2) {
      console.warn("No molecule data found for visualization.");
      return;
    }

    // Add a slight delay to ensure container is fully rendered
    const timer = setTimeout(() => {
      const viewer = $3Dmol.createViewer(viewerRef.current, {
        backgroundColor: "rgb(15, 23, 42)",
      });
      
      // Store viewer instance for later use
      setViewerInstance(viewer);
      
      viewer.clear();

      try {
        const model = viewer.addModel();
        setModelInstance(model);
        
        // Add atoms from molecule 1 if present
        if (molecule1) {
          // Apply any offset to molecule1 (typically fixed at origin)
          const mol1Atoms = molecule1.map((atom) => ({
            elem: atom.element,
            x: atom.x,
            y: atom.y,
            z: atom.z,
            properties: { molecule: 1 }  // Add property to identify molecule
          }));
          
          model.addAtoms(mol1Atoms);
          
          // Add element labels for molecule 1
          molecule1.forEach((atom) => {
            viewer.addLabel(atom.element, {
              position: { x: atom.x, y: atom.y, z: atom.z },
              fontSize: isMobile ? 10 : 12,
              fontColor: "white",
              backgroundColor: "rgba(56, 189, 248, 0.5)",
              borderRadius: 10,
              padding: isMobile ? 1 : 2,
              inFront: true,
            });
          });
        }
        
        // Add atoms from molecule 2 if present
        if (molecule2) {
          // Apply stored offset to molecule2
          const mol2Atoms = molecule2.map((atom) => ({
            elem: atom.element,
            x: atom.x + molecule2Offset.x,
            y: atom.y + molecule2Offset.y,
            z: atom.z + molecule2Offset.z,
            properties: { molecule: 2 }  // Add property to identify molecule
          }));
          
          model.addAtoms(mol2Atoms);
          
          // Add element labels for molecule 2
          molecule2.forEach((atom) => {
            viewer.addLabel(atom.element, {
              position: { 
                x: atom.x + molecule2Offset.x, 
                y: atom.y + molecule2Offset.y, 
                z: atom.z + molecule2Offset.z 
              },
              fontSize: isMobile ? 10 : 12,
              fontColor: "white",
              backgroundColor: "rgba(16, 185, 129, 0.5)",
              borderRadius: 10,
              padding: isMobile ? 1 : 2,
              inFront: true,
            });
          });
        }

        // Style atoms with different colors based on molecule
        viewer.setStyle({properties: {molecule: 1}}, {
          sphere: { 
            radius: isMobile ? 0.30 : 0.35, 
            scale: isMobile ? 0.85 : 0.9,
            color: "0x38bdf8"  // Blue color for molecule 1
          },
          stick: { 
            radius: isMobile ? 0.12 : 0.15,
            color: "0x38bdf8"  // Blue color for molecule 1
          },
        });
        
        viewer.setStyle({properties: {molecule: 2}}, {
          sphere: { 
            radius: isMobile ? 0.30 : 0.35, 
            scale: isMobile ? 0.85 : 0.9,
            color: "0x10b981"  // Green color for molecule 2
          },
          stick: { 
            radius: isMobile ? 0.12 : 0.15,
            color: "0x10b981"  // Green color for molecule 2
          },
        });
        
        // Add molecule labels if both molecules are present
        if (molecule1 && molecule2) {
          viewer.addLabel("Molecule 1 (Fixed)", {
            position: { x: molecule1[0].x, y: molecule1[0].y, z: molecule1[0].z + 5 },
            fontSize: isMobile ? 14 : 16,
            fontColor: "white",
            backgroundColor: "rgba(56, 189, 248, 0.7)",
            borderRadius: 10,
            padding: isMobile ? 2 : 4,
            inFront: true,
            fixedPosition: true,
          });
          
          viewer.addLabel(positioningMode ? "Molecule 2 (Draggable)" : "Molecule 2", {
            position: { 
              x: molecule2[0].x + molecule2Offset.x, 
              y: molecule2[0].y + molecule2Offset.y, 
              z: molecule2[0].z + molecule2Offset.z + 5 
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
        
        // Add positioning mode instructions if active
        if (positioningMode && molecule2) {
          viewer.addLabel("Drag to position Molecule 2", {
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
        
        // Setup bonds between atoms
        model.calculateBonds();
        
        // Center and zoom to fit the molecule(s)
        viewer.zoomTo();
        viewer.render();
        
        // Force a resize to ensure proper fit
        viewer.resize();
        
        // Additional render after short delay to ensure everything is displayed correctly
        setTimeout(() => {
          viewer.resize();
          viewer.render();
        }, 100);
      } catch (error) {
        console.error("Error rendering molecule:", error);
      }
    }, 50);

    // Set up responsive resize handler
    const handleResize = () => {
      if (viewerRef.current) {
        const viewer = $3Dmol.viewers[viewerRef.current.id];
        if (viewer) {
          viewer.resize();
          viewer.render();
        }
      }
    };
    
    window.addEventListener('resize', handleResize);
    
    // Cleanup
    return () => {
      window.removeEventListener('resize', handleResize);
      clearTimeout(timer);
    };
  }, [atoms, isMobile, molecule2Offset, positioningMode]);

  // Effect to handle mouse/touch events for dragging in positioning mode
  useEffect(() => {
    if (!positioningMode || !viewerInstance || !molecule2) return;

    const viewer = viewerRef.current;
    
    // Get viewer dimensions for calculations
    const viewerRect = viewer.getBoundingClientRect();
    const viewerWidth = viewerRect.width;
    const viewerHeight = viewerRect.height;

    // Mouse event handlers
    const handleMouseDown = (e) => {
      if (positioningMode) {
        setIsDragging(true);
        setDragStartPos({
          x: e.clientX,
          y: e.clientY
        });
        e.preventDefault();
      }
    };
    
    const handleMouseMove = (e) => {
      if (isDragging && positioningMode) {
        // Calculate drag delta in screen coordinates
        const deltaX = e.clientX - dragStartPos.x;
        const deltaY = e.clientY - dragStartPos.y;
        
        // Scale factor for movement (adjust for sensitivity)
        const scaleFactor = 0.05;
        
        // Get viewer orientation to map screen coordinates to model coordinates
        // This is a simplified approach that maps screen X,Y to model X,Z
        const newOffset = {
          x: molecule2Offset.x + (deltaX * scaleFactor),
          y: molecule2Offset.y - (deltaY * scaleFactor), // Invert Y axis
          z: molecule2Offset.z
        };
        
        // Update the offset
        onMoleculeMove(newOffset);
        
        // Update drag start position
        setDragStartPos({
          x: e.clientX,
          y: e.clientY
        });
        
        e.preventDefault();
      }
    };
    
    const handleMouseUp = (e) => {
      if (isDragging) {
        setIsDragging(false);
        e.preventDefault();
      }
    };
    
    // Touch event handlers (for mobile)
    const handleTouchStart = (e) => {
      if (positioningMode && e.touches.length === 1) {
        setIsDragging(true);
        setDragStartPos({
          x: e.touches[0].clientX,
          y: e.touches[0].clientY
        });
        e.preventDefault();
      }
    };
    
    const handleTouchMove = (e) => {
      if (isDragging && positioningMode && e.touches.length === 1) {
        // Calculate drag delta
        const deltaX = e.touches[0].clientX - dragStartPos.x;
        const deltaY = e.touches[0].clientY - dragStartPos.y;
        
        // Scale factor for movement
        const scaleFactor = 0.05;
        
        // Update the offset
        const newOffset = {
          x: molecule2Offset.x + (deltaX * scaleFactor),
          y: molecule2Offset.y - (deltaY * scaleFactor), // Invert Y axis
          z: molecule2Offset.z
        };
        
        onMoleculeMove(newOffset);
        
        // Update drag start position
        setDragStartPos({
          x: e.touches[0].clientX,
          y: e.touches[0].clientY
        });
        
        e.preventDefault();
      }
    };
    
    const handleTouchEnd = (e) => {
      if (isDragging) {
        setIsDragging(false);
        e.preventDefault();
      }
    };
    
    // Keyboard controls for more precise positioning
    const handleKeyDown = (e) => {
      if (!positioningMode) return;
      
      const moveStep = 0.5; // Angstroms per keypress
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
    };

    // Add event listeners
    viewer.addEventListener('mousedown', handleMouseDown);
    document.addEventListener('mousemove', handleMouseMove);
    document.addEventListener('mouseup', handleMouseUp);
    
    viewer.addEventListener('touchstart', handleTouchStart);
    document.addEventListener('touchmove', handleTouchMove);
    document.addEventListener('touchend', handleTouchEnd);
    
    document.addEventListener('keydown', handleKeyDown);
    
    // Cleanup function
    return () => {
      viewer.removeEventListener('mousedown', handleMouseDown);
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
      
      viewer.removeEventListener('touchstart', handleTouchStart);
      document.removeEventListener('touchmove', handleTouchMove);
      document.removeEventListener('touchend', handleTouchEnd);
      
      document.removeEventListener('keydown', handleKeyDown);
    };
  }, [positioningMode, isDragging, dragStartPos, viewerInstance, molecule2, molecule2Offset, onMoleculeMove]);

  return (
    <div 
      className="viewer-container" 
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
        cursor: positioningMode && molecule2 ? "move" : "auto"
      }}
    >
      <div
        ref={viewerRef}
        id={`molecule-viewer-${Math.random().toString(36).substr(2, 9)}`}
        style={{
          width: "100%",
          height: "100%",
          position: "absolute",
          top: 0,
          left: 0,
          borderRadius: "12px",
        }}
      ></div>
      
      {positioningMode && molecule2 && (
        <div style={{
          position: 'absolute',
          bottom: '10px',
          right: '10px',
          backgroundColor: 'rgba(0,0,0,0.7)',
          borderRadius: '5px',
          padding: '5px 10px',
          color: 'white',
          fontSize: '12px'
        }}>
          Offset: X: {molecule2Offset.x.toFixed(2)}, Y: {molecule2Offset.y.toFixed(2)}, Z: {molecule2Offset.z.toFixed(2)}
        </div>
      )}
    </div>
  );
};

export default MoleculeViewer;