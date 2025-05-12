# autoxidation-microkinetics

Python scripts to model the autoxidation of unsaturated lipids at representative atmospheric and intracellular concentrations.

## Overview

This repository contains two standalone scripts:

- **rate_solver_atmosphere.py**  
  Solves a coupled ODE system for lipid autoxidation under atmospheric‐like conditions (O₂ ≈ 0.2 atm).

- **rate_solver_cells.py**  
  Solves the same ODE system for intracellular‐like conditions (includes explicit O₂ consumption).

Each script defines a `rate_equations(…)` function, integrates the system with SciPy’s `solve_ivp` (Radau), and produces a time‐series plot of all species.

## Requirements

- **Python** 3.13.1  
- **NumPy**  
- **SciPy** (≥ 1.10)  
- **Matplotlib**

Install with:

```bash
pip install numpy scipy matplotlib
```

## Usage

1. Clone the repo
```
git clone https://github.com/asterlingchem/autoxidation-microkinetics.git
cd autoxidation-microkinetics
```

2. Atmospheric model
```
python rate_solver_atmosphere.py
```
→ produces coupled_rate_equations_atmosphere.pdf and .png

3. Intracellular model
 ```
python rate_solver_cells.py
```  
→ produces coupled_rate_equations_cells.pdf and .png

## Scripts

rate_solver_atmosphere.py

• Constant O₂ pressure = 0.2 atm

• 10 ODEs for species:

  R, OH, ROH, RO2, RO22, ALD, RO2_OH, POZ, VHP, VO

• Branching fractions for RO₂ + OH → RO₂_OH or POZ

• Integrates over 0–10000 s (1000 points)

• Saves a plot of species vs. time

rate_solver_cells.py

• O₂ as an explicit dynamic variable (intracellular concentration)

• 11 ODEs for species:

  R, OH, ROH, O2, RO2, RO22, ALD, RO2_OH, POZ, VHP, VO
  
• Same solver settings

• Saves a plot of species vs. time

## Parameters

All rate constants (k1…k7, kd) are defined at the top of each script.
Edit them in‐place to explore different barrier heights, diffusion limits, etc.

## Output

Each script outputs a .pdf and .png showing the time‐evolution of all reactive species (in atm).
Use these figures for analysis or publication.

## Citation

If you use this code in published work, please cite:

Gong, Z., Sterling, A.J., Wong, J., Head-Gordon, M., “Formation and consumption of Criegee intermediates in the autoxidation of unsaturated lipids” _Journal_ **2025**, DOI: TBC.

## License

This project is distributed under the MIT License.

## Contributing & Contact

Issues and pull requests are welcome!
For questions, open an issue or email alistair.sterling@utdallas.edu.

_This readme was written with assistance from ChatGPT_
