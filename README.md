# Simulating Mid-Latitude Cyclones using the Shallow Water Model

Welcome to the official repository for **Simulating Mid-Latitude Cyclones using the Shallow Water Model**!

This project models the dynamics of mid-latitude cyclones by solving the shallow water equations on a rotating sphere. The simulation explores cyclone evolution under varying physical parameters such as friction, solar forcing, and fluid depth.

For a full technical explanation, detailed methodology, validation, and analysis, **please refer to the accompanying research paper**:  
📄 **SMLRCSWM Documentation**.

---

## Repository Structure

- **`textures/`**  
  Contains `.jpg` images of the Earth used for rendering during the simulation.

- **`results/`**  
  Stores simulation output such as diagnostic data, figures, and screenshots used in the research paper.

- **`sources/`**
  - **`source_old/`**  
    Early prototypes and developmental versions of the simulation.
  - **`source_current/`**  
    Main project code.

---

## Important Files

- **`cyclone_sim_varnormalV2.m`**  
  The latest and best version of the simulation model.  
  ➔ This file includes both **validation mode** and **study mode** (see SMLRCSWM Documentation for details).

- **`*_nonnormalVX.m` Files**  
  Older versions where velocity values were not normalized. These are mainly for historical reference or exploratory purposes.

---

## Getting Started

1. Clone the repository:
    ```bash
    git clone https://github.com/nsmlord/Cyclone-Simulation.git
    ```
2. Navigate to the `sources/source_current/` folder.
3. Open and run **`cyclone_sim_varnormalV2.m`** in MATLAB.
4. Adjust parameters for validation or study mode as explained in the SMLRCSWM Documentation.

---

## How to Contribute

Feel free to fork the repo, explore different simulation parameters, improve visualizations, or add new model variations!  
Pull requests are welcome.

---

## Acknowledgments

Special thanks to the development of classical shallow water modeling techniques and to all resources referenced in the **SMLRCSWM Documentation**.

---

# 🚀 Happy simulating!

---
