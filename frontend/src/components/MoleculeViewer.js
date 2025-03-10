import React, { useEffect, useRef } from 'react';
import * as $3Dmol from '3dmol';
import { styles } from '../styles/components';

const MoleculeViewer = ({ atoms }) => {
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

      // Add element labels
      atoms.forEach((atom) => {
        viewer.addLabel(atom.element, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          fontSize: 14,
          fontColor: "white",
          backgroundColor: "rgba(0, 0, 0, 0.5)",
          borderRadius: 10,
          padding: 3,
          inFront: true,
        });
      });

      // Style atoms and bonds with enhanced visuals
      viewer.setStyle({}, {
        sphere: { radius: 0.35, scale: 0.9, colorscheme: 'Jmol' },
        stick: { radius: 0.15, colorscheme: 'Jmol' },
      });
      
      // Add subtle rotation for 3D effect
    //   viewer.spin('y', 0.5);
      
      viewer.zoomTo();
      viewer.render();
    } catch (error) {
      console.error("Error rendering molecule:", error);
    }
  }, [atoms]);

  return (
    <div className="viewer-container" style={styles.viewerContainer}>
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