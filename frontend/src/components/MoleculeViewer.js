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
  
  // Measurement state
  const [measurementMode, setMeasurementMode] = useState(null); // null, 'distance', 'angle', 'dihedral'
  const [selectedAtoms, setSelectedAtoms] = useState([]);
  const [measurements, setMeasurements] = useState([]);
  const [measurementLabels, setMeasurementLabels] = useState({});
  const spheresRef = useRef([]);

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

  // Calculate the distance between two points in 3D space
  const calculateDistance = (point1, point2) => {
    const dx = point2.x - point1.x;
    const dy = point2.y - point1.y;
    const dz = point2.z - point1.z;
    return Math.sqrt(dx*dx + dy*dy + dz*dz);
  };

  // Calculate the angle between three points in 3D space
  const calculateAngle = (point1, point2, point3) => {
    // Create vectors from point2 to point1 and point2 to point3
    const v1 = {
      x: point1.x - point2.x,
      y: point1.y - point2.y,
      z: point1.z - point2.z
    };
    
    const v2 = {
      x: point3.x - point2.x,
      y: point3.y - point2.y,
      z: point3.z - point2.z
    };
    
    // Calculate magnitudes
    const v1Mag = Math.sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
    const v2Mag = Math.sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
    
    // Calculate dot product
    const dotProduct = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    
    // Calculate angle in radians and convert to degrees
    const cosTheta = dotProduct / (v1Mag * v2Mag);
    // Ensure cosTheta is within -1 to 1 range due to potential floating point errors
    const clampedCosTheta = Math.max(-1, Math.min(1, cosTheta));
    return Math.acos(clampedCosTheta) * (180 / Math.PI);
  };

  // Calculate the dihedral angle between four points in 3D space
  const calculateDihedral = (point1, point2, point3, point4) => {
    // Create vectors
    const b1 = {
      x: point2.x - point1.x,
      y: point2.y - point1.y,
      z: point2.z - point1.z
    };
    
    const b2 = {
      x: point3.x - point2.x,
      y: point3.y - point2.y,
      z: point3.z - point2.z
    };
    
    const b3 = {
      x: point4.x - point3.x,
      y: point4.y - point3.y,
      z: point4.z - point3.z
    };
    
    // Calculate normal vectors to the two planes
    const n1 = {
      x: b1.y*b2.z - b1.z*b2.y,
      y: b1.z*b2.x - b1.x*b2.z,
      z: b1.x*b2.y - b1.y*b2.x
    };
    
    const n2 = {
      x: b2.y*b3.z - b2.z*b3.y,
      y: b2.z*b3.x - b2.x*b3.z,
      z: b2.x*b3.y - b2.y*b3.x
    };
    
    // Calculate magnitudes
    const n1Mag = Math.sqrt(n1.x*n1.x + n1.y*n1.y + n1.z*n1.z);
    const n2Mag = Math.sqrt(n2.x*n2.x + n2.y*n2.y + n2.z*n2.z);
    
    // Direction of the middle bond
    const b2Norm = Math.sqrt(b2.x*b2.x + b2.y*b2.y + b2.z*b2.z);
    const b2UnitX = b2.x / b2Norm;
    const b2UnitY = b2.y / b2Norm;
    const b2UnitZ = b2.z / b2Norm;
    
    // Calculate the components
    const m1 = {
      x: n1.x / n1Mag,
      y: n1.y / n1Mag,
      z: n1.z / n1Mag
    };
    
    const m2 = {
      x: n2.x / n2Mag,
      y: n2.y / n2Mag,
      z: n2.z / n2Mag
    };
    
    // Calculate the cross product of m2 and b2 unit vector
    const cp = {
      x: m2.y*b2UnitZ - m2.z*b2UnitY,
      y: m2.z*b2UnitX - m2.x*b2UnitZ,
      z: m2.x*b2UnitY - m2.y*b2UnitX
    };
    
    // Calculate the angle
    const cosTheta = m1.x*m2.x + m1.y*m2.y + m1.z*m2.z;
    const sinTheta = m1.x*cp.x + m1.y*cp.y + m1.z*cp.z;
    
    const dihedral = Math.atan2(sinTheta, cosTheta) * (180 / Math.PI);
    
    // Convert to the range of 0 to 360 degrees
    return dihedral < 0 ? dihedral + 360 : dihedral;
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

  // Function to handle atom clicks for measurements
  const handleAtomClick = (atom) => {
    if (!measurementMode) return;
    
    // Create a copy of the selected atoms
    const updatedSelection = [...selectedAtoms];
    updatedSelection.push(atom);
    
    // Check if we've selected enough atoms for the current measurement type
    let calculation = null;
    
    if (measurementMode === 'distance' && updatedSelection.length === 2) {
      // Calculate distance between two atoms
      const distance = calculateDistance(
        updatedSelection[0].coords,
        updatedSelection[1].coords
      );
      
      calculation = {
        type: 'distance',
        value: distance.toFixed(3),
        atoms: updatedSelection.map(a => ({ 
          id: a.id, 
          element: a.element,
          coords: {...a.coords} 
        }))
      };
    } 
    else if (measurementMode === 'angle' && updatedSelection.length === 3) {
      // Calculate angle between three atoms
      const angle = calculateAngle(
        updatedSelection[0].coords,
        updatedSelection[1].coords,
        updatedSelection[2].coords
      );
      
      calculation = {
        type: 'angle',
        value: angle.toFixed(2),
        atoms: updatedSelection.map(a => ({ 
          id: a.id, 
          element: a.element,
          coords: {...a.coords} 
        }))
      };
    } 
    else if (measurementMode === 'dihedral' && updatedSelection.length === 4) {
      // Calculate dihedral angle between four atoms
      const dihedral = calculateDihedral(
        updatedSelection[0].coords,
        updatedSelection[1].coords,
        updatedSelection[2].coords,
        updatedSelection[3].coords
      );
      
      calculation = {
        type: 'dihedral',
        value: dihedral.toFixed(2),
        atoms: updatedSelection.map(a => ({ 
          id: a.id, 
          element: a.element,
          coords: {...a.coords} 
        }))
      };
    }
    
    if (calculation) {
      // Add the measurement
      setMeasurements(prev => [...prev, calculation]);
      
      // Reset selection if we've completed a measurement
      setSelectedAtoms([]);
      
      // Update 3D display to show the measurement
      displayMeasurement(calculation, viewerInstance);
    } else {
      // Still collecting atoms for this measurement
      setSelectedAtoms(updatedSelection);
      
      // Highlight the selected atom
      if (viewerInstance) {
        const sphere = viewerInstance.addSphere({
          center: {x: atom.coords.x, y: atom.coords.y, z: atom.coords.z},
          radius: 0.4,
          color: '#FFFF00',
          opacity: 0.7
        });
        
        spheresRef.current.push(sphere);
        viewerInstance.render();
      }
    }
  };

  // Function to display a measurement on the 3D viewer
  const displayMeasurement = (measurement, viewer) => {
    if (!viewer) return;
    
    try {
      // Remove selection indicator spheres
      spheresRef.current.forEach(sphere => {
        viewer.removeShape(sphere);
      });
      spheresRef.current = [];
      
      // Calculate label position (midpoint or center)
      let labelPosition;
      let color;
      let text;
      
      const atoms = measurement.atoms;
      
      if (measurement.type === 'distance') {
        // Distance: midpoint between the two atoms
        labelPosition = {
          x: (atoms[0].coords.x + atoms[1].coords.x) / 2,
          y: (atoms[0].coords.y + atoms[1].coords.y) / 2,
          z: (atoms[0].coords.z + atoms[1].coords.z) / 2
        };
        color = '#FFFF00'; // Yellow
        text = `${measurement.value} Å`;
        
        // Draw a line between the two atoms
        viewer.addCylinder({
          start: {
            x: atoms[0].coords.x,
            y: atoms[0].coords.y,
            z: atoms[0].coords.z
          },
          end: {
            x: atoms[1].coords.x,
            y: atoms[1].coords.y,
            z: atoms[1].coords.z
          },
          radius: 0.1,
          fromCap: 1,
          toCap: 1,
          color: color,
          opacity: 0.7
        });
      } 
      else if (measurement.type === 'angle') {
        // Angle: position slightly away from middle atom
        const v1 = {
          x: atoms[0].coords.x - atoms[1].coords.x,
          y: atoms[0].coords.y - atoms[1].coords.y,
          z: atoms[0].coords.z - atoms[1].coords.z
        };
        
        const v2 = {
          x: atoms[2].coords.x - atoms[1].coords.x,
          y: atoms[2].coords.y - atoms[1].coords.y,
          z: atoms[2].coords.z - atoms[1].coords.z
        };
        
        // Normalize and add the vectors, then scale to position label
        const sumMag = Math.sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z) + 
                        Math.sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
        
        labelPosition = {
          x: atoms[1].coords.x + (v1.x + v2.x) / sumMag * 2,
          y: atoms[1].coords.y + (v1.y + v2.y) / sumMag * 2,
          z: atoms[1].coords.z + (v1.z + v2.z) / sumMag * 2
        };
        
        color = '#00FFFF'; // Cyan
        text = `${measurement.value}°`;
        
        // Draw lines connecting the atoms
        viewer.addCylinder({
          start: {
            x: atoms[0].coords.x,
            y: atoms[0].coords.y,
            z: atoms[0].coords.z
          },
          end: {
            x: atoms[1].coords.x,
            y: atoms[1].coords.y,
            z: atoms[1].coords.z
          },
          radius: 0.1,
          fromCap: 1,
          toCap: 1,
          color: color,
          opacity: 0.7
        });
        
        viewer.addCylinder({
          start: {
            x: atoms[1].coords.x,
            y: atoms[1].coords.y,
            z: atoms[1].coords.z
          },
          end: {
            x: atoms[2].coords.x,
            y: atoms[2].coords.y,
            z: atoms[2].coords.z
          },
          radius: 0.1,
          fromCap: 1,
          toCap: 1,
          color: color,
          opacity: 0.7
        });
      } 
      else if (measurement.type === 'dihedral') {
        // Dihedral: center of the four atoms
        labelPosition = {
          x: (atoms[0].coords.x + atoms[1].coords.x + atoms[2].coords.x + atoms[3].coords.x) / 4,
          y: (atoms[0].coords.y + atoms[1].coords.y + atoms[2].coords.y + atoms[3].coords.y) / 4,
          z: (atoms[0].coords.z + atoms[1].coords.z + atoms[2].coords.z + atoms[3].coords.z) / 4 + 1
        };
        
        color = '#FF00FF'; // Magenta
        text = `${measurement.value}°`;
        
        // Draw lines connecting the atoms
        viewer.addCylinder({
          start: {
            x: atoms[0].coords.x,
            y: atoms[0].coords.y,
            z: atoms[0].coords.z
          },
          end: {
            x: atoms[1].coords.x,
            y: atoms[1].coords.y,
            z: atoms[1].coords.z
          },
          radius: 0.1,
          fromCap: 1,
          toCap: 1,
          color: color,
          opacity: 0.7
        });
        
        viewer.addCylinder({
          start: {
            x: atoms[1].coords.x,
            y: atoms[1].coords.y,
            z: atoms[1].coords.z
          },
          end: {
            x: atoms[2].coords.x,
            y: atoms[2].coords.y,
            z: atoms[2].coords.z
          },
          radius: 0.1,
          fromCap: 1,
          toCap: 1,
          color: color,
          opacity: 0.7
        });
        
        viewer.addCylinder({
          start: {
            x: atoms[2].coords.x,
            y: atoms[2].coords.y,
            z: atoms[2].coords.z
          },
          end: {
            x: atoms[3].coords.x,
            y: atoms[3].coords.y,
            z: atoms[3].coords.z
          },
          radius: 0.1,
          fromCap: 1,
          toCap: 1,
          color: color,
          opacity: 0.7
        });
      }
      
      // Add the measurement label
      const labelId = viewer.addLabel(text, {
        position: labelPosition,
        backgroundColor: color,
        fontColor: 'black',
        fontSize: 14,
        borderRadius: 10,
        padding: 5,
        inFront: true
      });
      
      // Store the label ID to allow removal later
      setMeasurementLabels(prev => ({
        ...prev,
        [measurement.type + '-' + measurements.length]: labelId
      }));
      
      viewer.render();
    } catch (e) {
      console.error("Error displaying measurement:", e);
    }
  };

  // Function to clear all measurements
  const clearMeasurements = () => {
    if (viewerInstance) {
      try {
        // Remove selection indicator spheres
        spheresRef.current.forEach(sphere => {
          viewerInstance.removeShape(sphere);
        });
        spheresRef.current = [];
        
        // Remove measurement labels and shapes
        Object.values(measurementLabels).forEach(label => {
          viewerInstance.removeLabel(label);
        });
        
        // Force re-render the scene by completely reloading the molecules
        updateMoleculeDisplay(viewerInstance);
        
        // Clear state
        setMeasurementLabels({});
        setMeasurements([]);
        setSelectedAtoms([]);
        
        viewerInstance.render();
      } catch (e) {
        console.error("Error clearing measurements:", e);
      }
    }
  };

  // Handle measurement mode change
  const handleMeasurementModeChange = (mode) => {
    // If selecting the same mode, toggle it off
    if (mode === measurementMode) {
      setMeasurementMode(null);
    } else {
      setMeasurementMode(mode);
    }
    
    // Clear current selection when changing modes
    setSelectedAtoms([]);
    
    // Remove selection indicator spheres
    if (viewerInstance) {
      spheresRef.current.forEach(sphere => {
        viewerInstance.removeShape(sphere);
      });
      spheresRef.current = [];
      viewerInstance.render();
    }
  };

  // Helper function to set up atom click handling
  const setupAtomClickHandling = (viewer) => {
    if (!viewer) return;
    
    // Enable atom picking
    viewer.setClickable({}, true, (atom) => {
      if (measurementMode) {
        const atomIndex = atom.serial;
        let atomData;
        
        // Need to map back to the original atom data from the props
        // First, look in molecule1
        if (molecule1 && atomIndex < molecule1.length) {
          atomData = {
            id: molecule1[atomIndex].id || atomIndex,
            element: molecule1[atomIndex].element,
            coords: {
              x: molecule1[atomIndex].x,
              y: molecule1[atomIndex].y, 
              z: molecule1[atomIndex].z
            },
            moleculeIndex: 1
          };
        } 
        // Then, look in molecule2
        else if (molecule2) {
          const mol2Index = atomIndex - (molecule1 ? molecule1.length : 0);
          if (mol2Index < molecule2.length) {
            // Apply any transformations
            let x = molecule2[mol2Index].x;
            let y = molecule2[mol2Index].y;
            let z = molecule2[mol2Index].z;
            
            // Apply rotation if needed
            if (molecule2Rotation && (molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0)) {
              const centerOfMass = calculateCenterOfMass(molecule2);
              const coords = [x, y, z];
              const rotated = applyRotation(coords, molecule2Rotation, centerOfMass);
              x = rotated[0];
              y = rotated[1];
              z = rotated[2];
            }
            
            // Apply offset
            atomData = {
              id: molecule2[mol2Index].id || mol2Index,
              element: molecule2[mol2Index].element,
              coords: {
                x: x + molecule2Offset.x,
                y: y + molecule2Offset.y,
                z: z + molecule2Offset.z
              },
              moleculeIndex: 2
            };
          }
        }
        
        if (atomData) {
          handleAtomClick(atomData);
        }
      }
    });
  };

  // Update molecule display with all models and styles
  const updateMoleculeDisplay = (viewer) => {
    if (!viewer) return;
    
    try {
      // Clear previous content
      viewer.clear();

      // Handle molecule1
      if (molecule1) {
        // Transform molecule1 data to 3Dmol format
        let m1Data = molecule1.map((atom, index) => ({
          elem: atom.element,
          x: atom.x,
          y: atom.y,
          z: atom.z,
          properties: { molecule: 1 },
          serial: index // Add serial for atom identification
        }));
        
        // Add molecule1 with its own model
        let model1 = viewer.addModel();
        model1.addAtoms(m1Data);
        
        // Add element labels for molecule 1
        molecule1.forEach((atom) => {
          viewer.addLabel(atom.element, {
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
        }
      }

      // Handle molecule2
      if (molecule2) {
        // Calculate center of mass for rotation
        const centerOfMass = calculateCenterOfMass(molecule2);
        
        // Process each atom with rotation and offset
        let m2Data = molecule2.map((atom, index) => {
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
            properties: { molecule: 2 },
            serial: (molecule1 ? molecule1.length : 0) + index // Add serial for atom identification
          };
        });
        
        // Add molecule2 with its own model
        let model2 = viewer.addModel();
        model2.addAtoms(m2Data);
        
        // Add element labels for molecule2
        molecule2.forEach((atom, index) => {
          // Apply same transformations as for the atoms
          let x = atom.x, y = atom.y, z = atom.z;
          
          if (molecule2Rotation && (molecule2Rotation.x !== 0 || molecule2Rotation.y !== 0 || molecule2Rotation.z !== 0)) {
            const coords = [x, y, z];
            const rotated = applyRotation(coords, molecule2Rotation, centerOfMass);
            x = rotated[0];
            y = rotated[1];
            z = rotated[2];
          }
          
          viewer.addLabel(atom.element, {
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
          
          viewer.addLabel("Molecule 2 (Use Arrow Keys)", {
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
      viewer.setStyle({properties: {molecule: 1}}, {
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
      
      viewer.setStyle({properties: {molecule: 2}}, {
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
        viewer.addLabel("Use arrow keys to position Molecule 2", {
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
      
      // Restore any existing measurements
      measurements.forEach(measurement => {
        displayMeasurement(measurement, viewer);
      });
      
      // Setup click handling for measurements
      if (measurementMode) {
        setupAtomClickHandling(viewer);
      }
      
      // Disable default mouse handling in positioning mode
      if (positioningMode) {
        viewer.setClickable(false, true);
        viewer.setHoverable(false, true);
      }
      
      // Force render and resize
      viewer.render();
      viewer.resize();
    } catch (error) {
      console.error("Error rendering molecule:", error);
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
        
        // Update the display
        updateMoleculeDisplay(viewerInstance);
        
        // Initial view or restore camera
        if (isInitialRender) {
          viewerInstance.zoomTo();
          setIsInitialRender(false);
        } else if (cameraStateRef.current) {
          restoreCameraState(viewerInstance, cameraStateRef.current);
        }
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

  // Effect to handle measurement mode changes
  useEffect(() => {
    if (!viewerInstance) return;
    
    // If measurement mode is active, set up click handling
    if (measurementMode) {
      setupAtomClickHandling(viewerInstance);
    } 
    // Otherwise disable atom clicking
    else {
      viewerInstance.setClickable({}, false);
      
      // Clear any partial selections
      if (selectedAtoms.length > 0) {
        setSelectedAtoms([]);
        
        // Clear selection indicator spheres
        spheresRef.current.forEach(sphere => {
          viewerInstance.removeShape(sphere);
        });
        spheresRef.current = [];
        viewerInstance.render();
      }
    }
  }, [measurementMode, viewerInstance, molecule1, molecule2]);

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
    <div>
      {/* Measurement Tool Controls */}
      {(!positioningMode && (molecule1 || molecule2)) && (
        <div style={{
          display: 'flex',
          flexWrap: 'wrap',
          gap: '8px',
          marginBottom: '12px',
          justifyContent: isMobile ? 'center' : 'flex-start',
          backgroundColor: 'rgba(16, 24, 39, 0.7)',
          padding: '8px',
          borderRadius: '8px',
          border: '1px solid rgba(255, 255, 255, 0.1)'
        }}>
          <div style={{ fontWeight: 500, marginRight: '8px', display: 'flex', alignItems: 'center' }}>
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" style={{ marginRight: '6px' }}>
              <path d="M22 12h-4l-3 9L9 3l-3 9H2" />
            </svg>
            Measurement:
          </div>
          
          <button
            onClick={() => handleMeasurementModeChange('distance')}
            style={{
              backgroundColor: measurementMode === 'distance' ? 'rgba(255, 255, 0, 0.2)' : 'rgba(16, 24, 39, 0.7)',
              color: measurementMode === 'distance' ? '#FFFF00' : '#f0f4f8',
              border: `1px solid ${measurementMode === 'distance' ? 'rgba(255, 255, 0, 0.5)' : 'rgba(255, 255, 255, 0.1)'}`,
              padding: '4px 10px',
              borderRadius: '6px',
              fontSize: '0.85rem',
              display: 'flex',
              alignItems: 'center',
              gap: '5px'
            }}
          >
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M2 12h20" />
              <path d="M2 18h20" />
              <path d="M2 6h20" />
            </svg>
            Distance
            {measurementMode === 'distance' && selectedAtoms.length > 0 && (
              <span style={{ marginLeft: '4px', fontSize: '0.7rem', opacity: 0.8 }}>
                ({selectedAtoms.length}/2)
              </span>
            )}
          </button>
          
          <button
            onClick={() => handleMeasurementModeChange('angle')}
            style={{
              backgroundColor: measurementMode === 'angle' ? 'rgba(0, 255, 255, 0.2)' : 'rgba(16, 24, 39, 0.7)',
              color: measurementMode === 'angle' ? '#00FFFF' : '#f0f4f8',
              border: `1px solid ${measurementMode === 'angle' ? 'rgba(0, 255, 255, 0.5)' : 'rgba(255, 255, 255, 0.1)'}`,
              padding: '4px 10px',
              borderRadius: '6px',
              fontSize: '0.85rem',
              display: 'flex',
              alignItems: 'center',
              gap: '5px'
            }}
          >
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M21 7L14 14M12.5 5L15 2M19 12l3 -2.5M5 21l7 -7" />
            </svg>
            Angle
            {measurementMode === 'angle' && selectedAtoms.length > 0 && (
              <span style={{ marginLeft: '4px', fontSize: '0.7rem', opacity: 0.8 }}>
                ({selectedAtoms.length}/3)
              </span>
            )}
          </button>
          
          <button
            onClick={() => handleMeasurementModeChange('dihedral')}
            style={{
              backgroundColor: measurementMode === 'dihedral' ? 'rgba(255, 0, 255, 0.2)' : 'rgba(16, 24, 39, 0.7)',
              color: measurementMode === 'dihedral' ? '#FF00FF' : '#f0f4f8',
              border: `1px solid ${measurementMode === 'dihedral' ? 'rgba(255, 0, 255, 0.5)' : 'rgba(255, 255, 255, 0.1)'}`,
              padding: '4px 10px',
              borderRadius: '6px',
              fontSize: '0.85rem',
              display: 'flex',
              alignItems: 'center',
              gap: '5px'
            }}
          >
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M4 14h6l4 -8l3 16h7" />
            </svg>
            Dihedral
            {measurementMode === 'dihedral' && selectedAtoms.length > 0 && (
              <span style={{ marginLeft: '4px', fontSize: '0.7rem', opacity: 0.8 }}>
                ({selectedAtoms.length}/4)
              </span>
            )}
          </button>
          
          {measurements.length > 0 && (
            <button
              onClick={clearMeasurements}
              style={{
                backgroundColor: 'rgba(244, 63, 94, 0.2)',
                color: '#f43f5e',
                border: '1px solid rgba(244, 63, 94, 0.5)',
                padding: '4px 10px',
                borderRadius: '6px',
                fontSize: '0.85rem',
                marginLeft: 'auto',
                display: 'flex',
                alignItems: 'center',
                gap: '5px'
              }}
            >
              <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M18 6L6 18" />
                <path d="M6 6l12 12" />
              </svg>
              Clear All
            </button>
          )}
        </div>
      )}
      
      {/* Measurement Instructions */}
      {measurementMode && (
        <div style={{
          backgroundColor: 'rgba(59, 130, 246, 0.1)',
          border: '1px solid rgba(59, 130, 246, 0.3)',
          borderRadius: '8px',
          padding: '8px 12px',
          marginBottom: '12px',
          fontSize: '0.85rem',
          color: '#93c5fd'
        }}>
          <div style={{ fontWeight: 600, marginBottom: '4px', display: 'flex', alignItems: 'center', gap: '6px' }}>
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="12" cy="12" r="10" />
              <path d="M12 16v-4" />
              <path d="M12 8h.01" />
            </svg>
            Click on atoms to measure {measurementMode}
          </div>
          {measurementMode === 'distance' && (
            <div>Select 2 atoms to measure the distance between them.</div>
          )}
          {measurementMode === 'angle' && (
            <div>Select 3 atoms to measure the angle formed. The second atom will be the vertex.</div>
          )}
          {measurementMode === 'dihedral' && (
            <div>Select 4 atoms to measure the dihedral angle between the planes formed by atoms 1-2-3 and 2-3-4.</div>
          )}
        </div>
      )}
      
      {/* Measurements List */}
      {measurements.length > 0 && (
        <div style={{
          backgroundColor: 'rgba(16, 24, 39, 0.7)',
          border: '1px solid rgba(255, 255, 255, 0.1)',
          borderRadius: '8px',
          padding: '8px',
          marginBottom: '12px',
          maxHeight: isMobile ? '120px' : '150px',
          overflowY: 'auto'
        }}>
          <div style={{ fontWeight: 600, marginBottom: '6px', fontSize: '0.9rem' }}>Measurements:</div>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
            {measurements.map((m, i) => (
              <div key={i} style={{
                backgroundColor: m.type === 'distance' ? 'rgba(255, 255, 0, 0.1)' : 
                                m.type === 'angle' ? 'rgba(0, 255, 255, 0.1)' : 'rgba(255, 0, 255, 0.1)',
                border: `1px solid ${m.type === 'distance' ? 'rgba(255, 255, 0, 0.3)' : 
                                    m.type === 'angle' ? 'rgba(0, 255, 255, 0.3)' : 'rgba(255, 0, 255, 0.3)'}`,
                borderRadius: '6px',
                padding: '4px 8px',
                fontSize: '0.85rem'
              }}>
                <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                  <span style={{ fontWeight: 500 }}>
                    {m.type === 'distance' ? 'Distance:' : 
                     m.type === 'angle' ? 'Angle:' : 'Dihedral:'}
                  </span>
                  <span style={{ fontWeight: 600 }}>
                    {m.value}{m.type === 'distance' ? ' Å' : '°'}
                  </span>
                </div>
                <div style={{ fontSize: '0.75rem', opacity: 0.7 }}>
                  {m.atoms.map(a => a.element).join('-')}
                </div>
              </div>
            ))}
          </div>
        </div>
      )}
      
      {/* 3D Viewer Container */}
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
          cursor: measurementMode ? "crosshair" : "auto",
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
        
        {measurementMode && (
          <div style={{
            position: 'absolute',
            bottom: '10px',
            left: '10px',
            backgroundColor: measurementMode === 'distance' ? 'rgba(255, 255, 0, 0.2)' : 
                            measurementMode === 'angle' ? 'rgba(0, 255, 255, 0.2)' : 'rgba(255, 0, 255, 0.2)',
            borderRadius: '5px',
            padding: '8px 10px',
            color: 'white',
            fontSize: '12px',
            border: `1px solid ${measurementMode === 'distance' ? 'rgba(255, 255, 0, 0.5)' : 
                                measurementMode === 'angle' ? 'rgba(0, 255, 255, 0.5)' : 'rgba(255, 0, 255, 0.5)'}`,
            zIndex: 20
          }}>
            <div style={{ fontWeight: 'bold' }}>
              {measurementMode === 'distance' ? 'Distance Measurement' : 
              measurementMode === 'angle' ? 'Angle Measurement' : 'Dihedral Measurement'}
            </div>
            <div>
              Click on {
                measurementMode === 'distance' ? '2' : 
                measurementMode === 'angle' ? '3' : '4'
              } atoms to measure
            </div>
            {selectedAtoms.length > 0 && (
              <div style={{ marginTop: '4px' }}>
                Selected: {selectedAtoms.length}/{
                  measurementMode === 'distance' ? '2' : 
                  measurementMode === 'angle' ? '3' : '4'
                }
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
};

export default MoleculeViewer;