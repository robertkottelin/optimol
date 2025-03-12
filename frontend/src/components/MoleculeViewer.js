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

    // Add a slight delay to ensure container is fully rendered
    const timer = setTimeout(() => {
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
        
        // Center and zoom to fit the molecule
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
  }, [atoms, isMobile]);

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
        border: "1px solid rgba(255, 255, 255, 0.1)"
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
    </div>
  );
};

export default MoleculeViewer;