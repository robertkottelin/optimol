import React, { useState } from "react";
import axios from "axios";

const App = () => {
  const [optimizedMolecule, setOptimizedMolecule] = useState(null); // State to store optimized molecule

  const handleFileUpload = (event) => {
    const file = event.target.files[0];

    if (!file) {
      alert("Please select a file.");
      return;
    }

    const reader = new FileReader();

    reader.onload = async (e) => {
      try {
        const fileContent = e.target.result;
        const moleculeData = JSON.parse(fileContent); // Parse the JSON file content

        if (!validateMoleculeJSON(moleculeData)) {
          alert("Invalid molecule JSON format.");
          return;
        }

        const response = await axios.post("http://localhost:5000/optimize", moleculeData);

        // Set optimized molecule in state to display it
        setOptimizedMolecule(response.data.optimized_file1);
      } catch (error) {
        console.error("Error processing the file or sending data to API:", error);
        alert("Error uploading the file. Check the console for details.");
      }
    };

    reader.readAsText(file); // Read the file as text
  };

  const validateMoleculeJSON = (data) => {
    // Basic validation to ensure the JSON structure matches the expected format
    return (
      data &&
      data.file1 &&
      Array.isArray(data.file1.atoms) &&
      data.file1.atoms.every(
        (atom) =>
          atom.id &&
          atom.element &&
          typeof atom.x === "number" &&
          typeof atom.y === "number" &&
          typeof atom.z === "number"
      )
    );
  };

  return (
    <div style={{ padding: "20px", textAlign: "center" }}>
      <h1>Molecule Optimizer</h1>
      <input
        type="file"
        onChange={handleFileUpload}
        style={{
          padding: "10px",
          border: "2px dashed gray",
          cursor: "pointer",
          marginBottom: "20px",
        }}
      />


      {optimizedMolecule && (
        <div
          style={{
            marginTop: "20px",
            textAlign: "left",
            padding: "10px",
            border: "1px solid gray",
            borderRadius: "5px",
            backgroundColor: "#f9f9f9",
            maxWidth: "600px",
            margin: "auto",
          }}
        >
          <h2>Optimized Molecule</h2>
          <pre>{JSON.stringify(optimizedMolecule, null, 2)}</pre>
        </div>
      )}
    </div>
  );
};

export default App;