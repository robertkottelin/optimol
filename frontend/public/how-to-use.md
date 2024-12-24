## How To Use & Theory:

This app optimizes molecular geometry using either classical or quantum energy optimization techniques.

---

### Steps to Use:

1. **Prepare a JSON file** containing your molecular geometry in the following format. The file should include atomic positions (x, y, z) and element types (H, O, etc.). This is only an example format; you can expand the file structure to represent bigger molecules as needed, or just add random atoms to random positions and see what happens!

   ```json
   {
     "file1": {
       "atoms": [
         { "id": 1, "element": "H", "x": 0.0, "y": 0.0, "z": 0.0 },
         { "id": 2, "element": "H", "x": 0.0, "y": 0.0, "z": 0.74 },
         { "id": 3, "element": "O", "x": 0.0, "y": 0.74, "z": 0.0 }
       ]
     }
   }
   ```

2. **Upload the JSON file** using the "Upload File" button in the app.

3. **Modify optimization parameters** (described below).

4. Click **"Classical Optimize"** or **"Quantum Optimize"** depending on the optimization method you want to use.

5. **Download the optimized molecule data** as a JSON file for further analysis or visualization.

---

### Classical Optimization (BFGS Optimizer)

The classical optimization process uses the **BFGS (Broyden–Fletcher–Goldfarb–Shanno)** algorithm to minimize the molecular geometry energy. 

#### Theory:
- Classical optimization aims to find a stable molecular structure by minimizing the potential energy using the **Effective Medium Theory (EMT)** calculator from ASE.
- The optimization process iteratively adjusts the atomic positions until the forces acting on atoms fall below a specified convergence criterion (**fmax**) or a maximum number of steps is reached.

#### Parameters:
1. **fmax** (Force Convergence Criterion):
   - Specifies the maximum force tolerance before considering the structure optimized.
   - **Default:** 0.005
   - **Subscribed Users Range:** 0.001 to 0.1
   - **Unsubscribed Users Range:** Fixed at 0.05.

2. **steps** (Maximum Iterations):
   - The maximum number of optimization steps allowed.
   - **Default:** 500
   - **Subscribed Users Range:** 10 to 1,000,000
   - **Unsubscribed Users Range:** 10 to 100.

#### Example Calculation:
- Uses **gradient descent-like techniques** to update atomic positions.
- Iteratively reduces the molecular potential energy calculated by the EMT model.

---

### Quantum Optimization (QAOA Algorithm)

The quantum optimization process employs the **Quantum Approximate Optimization Algorithm (QAOA)** to explore the molecular geometry energy minimization problem.

#### Theory:
- **QAOA** is a hybrid quantum-classical algorithm used to solve optimization problems by combining quantum circuits with classical parameter optimization.
- The algorithm constructs a **cost Hamiltonian** (representing the energy landscape of the molecular structure) and a **mixer Hamiltonian** (introducing state evolution).
- The optimization minimizes the expectation value of the cost Hamiltonian using variational parameters in the quantum circuit.

#### Parameters:
1. **maxiter** (Maximum Classical Iterations):
   - The number of iterations allowed for the classical optimization process.
   - **Default:** 1,000
   - **Subscribed Users Range:** 10 to 1,000,000
   - **Unsubscribed Users Range:** Fixed at 100.

2. **p** (QAOA Layers or Repetitions):
   - Determines the number of alternating layers of the cost and mixer Hamiltonians in the quantum circuit.
   - **Default:** 2
   - **Subscribed Users Range:** 1 to 1,000
   - **Unsubscribed Users Range:** Fixed at 1.

3. **optimizer** (Classical Optimizer):
   - Specifies the classical optimizer for parameter tuning.
   - Options include:
     - **COBYLA (Constrained Optimization BY Linear Approximations):** Suitable for constrained optimization.
     - **SPSA (Simultaneous Perturbation Stochastic Approximation):** Efficient for noisy gradients.
   - **Default:** COBYLA.

#### Example Calculations:
- **Cost Hamiltonian**: \( H_c = \sum Z_i + \sum Z_j \), where \( Z \) represents Pauli-Z operators.
- **Mixer Hamiltonian**: \( H_m = \sum X_i \), where \( X \) represents Pauli-X operators.
- **Objective Function**: The expectation value of \( H_c \) given the variational parameters.

---

### Optimization Process

1. **Classical Optimization**:
   - Uses the **BFGS optimizer** to iteratively refine atomic positions based on force and energy calculations.
   - Ensures the final structure is stable and energetically favorable.

2. **Quantum Optimization**:
   - Constructs a quantum circuit based on the **QAOA Ansatz**.
   - Minimizes the molecular energy using a hybrid quantum-classical approach.
   - Outputs variational parameters and the adjusted atomic positions.

---

### Output and Results

- The app returns a JSON file with the optimized atomic positions, molecular energy, and (for quantum optimization) the variational parameters.

#### Example Output:
```json
{
  "atoms": [
    { "id": 1, "element": "H", "x": 0.0, "y": 0.0, "z": 0.0 },
    { "id": 2, "element": "H", "x": 0.0, "y": 0.0, "z": 0.75 },
    { "id": 3, "element": "O", "x": 0.0, "y": 0.75, "z": 0.0 }
  ],
  "optimal_params": [1.2, -0.5, 0.8],
  "min_energy": -1.23456
}
```

---

### How to Use the App

1. Upload the molecular geometry JSON file.
2. Select the desired optimization method (Classical or Quantum).
3. Adjust parameters such as `fmax`, `steps`, `maxiter`, and `p` if needed.
4. Click the optimization button to start the calculation.
5. Download the optimized molecule data for further analysis or visualization.

By using the **Optimize Molecule** app, you can explore molecular structures with both classical and quantum computational techniques, providing flexibility and precision for your research needs.

### Developer Roadmap:
- Expand quantum calculations to include more advanced quantum algorithms.
- Include user input for more tunable parameters.
- Implement additional classical optimization methods for comparison.
- Enhance visualization tools for optimized molecular structures.
- Molecule file converter (SDF, XYZ, etc.) for broader compatibility.
- Crypto payments if some don't want to use traditional payment methods.