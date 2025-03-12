import React, { useEffect, useRef } from 'react';
import * as $3Dmol from '3dmol';
import { styles } from '../styles/components';

const MoleculeViewer = ({ atoms, isMobile }) => {
  const viewerRef = useRef();

  useEffect(() => {
    if (!atoms) {
      console.warn("No molecule data found for visualization.");
      return;
    }

    const viewer = $3Dmol.createViewer(viewerRef.current, {
      backgroundColor: "rgb(15, 23, 42)",
    });
    viewer.clear();

    try {
      // Convert to 3Dmol atom format
      const mol3dAtoms = atoms.map((atom) => ({
        elem: atom.element,
        x: atom.x,
        y: atom.y,
        z: atom.z,
      }));

      const model = viewer.addModel();
      model.addAtoms(mol3dAtoms);

      // Add element labels - scale down for mobile
      atoms.forEach((atom) => {
        viewer.addLabel(atom.element, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          fontSize: isMobile ? 12 : 14,
          fontColor: "white",
          backgroundColor: "rgba(0, 0, 0, 0.5)",
          borderRadius: 10,
          padding: isMobile ? 2 : 3,
          inFront: true,
        });
      });

      // Style atoms and bonds with enhanced visuals
      // Use smaller sphere radius on mobile for better rendering
      viewer.setStyle({}, {
        sphere: { 
          radius: isMobile ? 0.30 : 0.35, 
          scale: isMobile ? 0.85 : 0.9, 
          colorscheme: 'Jmol' 
        },
        stick: { 
          radius: isMobile ? 0.12 : 0.15, 
          colorscheme: 'Jmol' 
        },
      });
      
      // Add controls info for touch devices
      if (isMobile) {
        viewer.addLabel("Touch: Rotate | Pinch: Zoom", {
          position: { x: 0, y: 0, z: 0 },
          fontSize: 10,
          fontColor: "white",
          backgroundColor: "rgba(0, 0, 0, 0.7)",
          padding: 4,
          inFront: true,
          fixed: true, // stays in place during rotation
        });
        
        // Clear the help text after 5 seconds
        setTimeout(() => {
          viewer.removeAllLabels();
          // Re-add the atom labels
          atoms.forEach((atom) => {
            viewer.addLabel(atom.element, {
              position: { x: atom.x, y: atom.y, z: atom.z },
              fontSize: 12,
              fontColor: "white",
              backgroundColor: "rgba(0, 0, 0, 0.5)",
              borderRadius: 10,
              padding: 2,
              inFront: true,
            });
          });
          viewer.render();
        }, 5000);
      }
      
      viewer.zoomTo();
      viewer.render();
      
      // Set up responsive resize handler
      const handleResize = () => {
        viewer.resize();
        viewer.render();
      };
      
      window.addEventListener('resize', handleResize);
      
      // Cleanup
      return () => {
        window.removeEventListener('resize', handleResize);
      };
    } catch (error) {
      console.error("Error rendering molecule:", error);
    }
  }, [atoms, isMobile]);

  return (
    <div 
      className="viewer-container" 
      style={{
        ...styles.viewerContainer,
        height: isMobile ? '350px' : '450px'
      }}
    >
      <div
        ref={viewerRef}
        style={{
          width: "100%",
          height: "100%",
          borderRadius: "12px",
        }}
      ></div>
    </div>
  );
};

export default MoleculeViewer;