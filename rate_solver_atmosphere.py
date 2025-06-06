import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# rate equations for full model
def rate_equations(t, y, k1, k2, k3, k4, k5, k6, k7, kd):
    R, OH, ROH, RO2, RO22, ALD, RO2_OH, POZ, VHP, VO = y

    O2 = 0.2  # use a constant pressure of 0.2 atm for O2

    # create branching fractions for RO2 + OH consumption
    kfr = kd + k3  # combined rate constant for RO2 + OH reaction
    f1 = kd / kfr  # branching fraction for RO2 + OH --> RO2_OH
    f2 = k3 / kfr  # branching fraction for RO2 + OH --> POZ

    # Rate equations
    dR_dt = - kd * R * OH
    dOH_dt = - kd * R * OH + k1 * RO2 - kfr * RO2 * OH + k2 * RO2_OH + k2 * VHP + k6 * VO * O2
    dROH_dt = kd * R * OH - kd * ROH * O2
    dRO2_dt = kd * ROH * O2 - k1 * RO2 - kfr * RO2 * OH + k2 * RO2_OH - 2 * kd * RO2**2 + k7 * RO22
    dRO22_dt = kd * RO2**2 - k7 * RO22 - k4 * RO22
    dALD_dt = k1 * RO2 + k5 * POZ + k6 * VO * O2
    dRO2_OH_dt = f1 * kfr * RO2 * OH - k2 * RO2_OH
    dPOZ_dt = k4 * RO22 + f2 * kfr * RO2 * OH - k5 * POZ
    dVHP_dt = k5 * POZ - k2 * VHP
    dVO_dt = k2 * VHP - k6 * VO * O2

    return [dR_dt, dOH_dt, dROH_dt, dRO2_dt, dRO22_dt, dALD_dt, dRO2_OH_dt, dPOZ_dt, dVHP_dt, dVO_dt]

# Initial concentrations
R_0 = 4.05*10**-12  # Initial concentration of R, in atm
OH_0 = 4.05*10**-12   # Initial concentration of OH
ROH_0 = 0.0  # Initial concentration of ROH
RO2_0 = 0.0  # Initial concentration of RO2
RO22_0 = 0.0  # Initial concentration of RO2 dimer
ALD_0 = 0.0  # Initial concentration of aldehyde+other products
RO2_OH_0 = 0.0  # Initial concentration of RO2_OH
POZ_0 = 0.0  # Initial concentration of POZ
VHP_0 = 0.0  # Initial concentration of VHP
VO_0 = 0.0  # Initial concentration of VO

# Gas phase rate constants for complex reaction
k1 = 3.14*10**-4  # 1/s, RO2 --> P_uni + OH
k2 = 1.0*10**-5 # 1/s, RO2_OH + hv --> RO2 + OH and VHP + hv --> VO + OH
k3 = 2.16*10**5  # 1/(atm*s), RO2 + OH --> POZ
k4 = 6.27*10**1  # 1/(atm*s), est. 15 kcal/mol barrier, (RO2)2 --> POZ
k5 = 6.53*10**-2  # 1/s, POZ --> VHP  
k6 = 1.79*10**-3  # 1/(atm*s), VO + O2 --> aldehyde + OH
k7 = 10**8  # 1/*s, (RO2)2 --> RO2 + RO2
kd = 1.0*10**9  # 1/(atm*s), diffusion limited reaction

# Time span
t_span = (0, 10000)  # Simulate from t=0 to t=10000 seconds
t_eval = np.linspace(*t_span, 10000)  # Generate 10000 time points

# Solve ODE system
solution = solve_ivp(rate_equations, t_span, [R_0, OH_0, ROH_0, RO2_0, RO22_0, ALD_0, RO2_OH_0, POZ_0, VHP_0, VO_0], args=(k1, k2, k3, k4, k5, k6, k7, kd), t_eval=t_eval, method='Radau', rtol=1e-6, atol=1e-12)

# Concentrations
# Unpack species from solution.y
R       = solution.y[0]
OH      = solution.y[1]
ROH     = solution.y[2]
RO2     = solution.y[3]
RO22    = solution.y[4]
ALD     = solution.y[5]
RO2_OH  = solution.y[6]
POZ     = solution.y[7]
VHP     = solution.y[8]
VO      = solution.y[9]

# Calculate total reactive species (excluding O2, since it's in large excess)
total = R + OH + ROH + RO2 + RO22 + ALD + RO2_OH + POZ + VHP + VO
print(f"Total reactive species (start): {total[0]} atm")
print(f"Total reactive species (end): {total[-1]} atm")

# Calculate total carbon atoms
total_C_atoms = R + ROH + RO2 + RO2_OH + ALD + POZ + VHP + VO + 2 * RO22
print(f"Total carbon atoms (start): {total_C_atoms[0]} atm-equivalent")
print(f"Total carbon atoms (end): {total_C_atoms[-1]} atm-equivalent")

# Plot results
plt.figure(figsize=(8, 6))

plt.plot(solution.t, solution.y[0], label='R', color='blue', alpha=0.8)
plt.plot(solution.t, solution.y[1], label='OH', color='black', ls='--', alpha=0.8)
plt.plot(solution.t, solution.y[2], label='ROH', color='green', alpha=0.8)
# plt.plot(solution.t, solution.y[3], label='O2', color='pink', alpha=0.8)
plt.plot(solution.t, solution.y[3], label='RO2', color='purple', alpha=0.8)
plt.plot(solution.t, solution.y[4], label='RO2 dimer', color='orange', alpha=0.8)
plt.plot(solution.t, solution.y[5], label='ALD', color='red', alpha=0.8)
plt.plot(solution.t, solution.y[6], label='RO2_OH', color='blue', alpha=0.8, ls='--')
plt.plot(solution.t, solution.y[7], label='POZ', color='brown', alpha=0.8)
plt.plot(solution.t, solution.y[8], label='VHP', color='gray', alpha=0.8)
plt.plot(solution.t, solution.y[9], label='VO', color='black', alpha=0.8)

plt.xlabel('Time / s', size=16)
plt.ylabel('Concentration / atm', size=16)
plt.xticks(size=14)
plt.yticks(size=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('coupled_rate_equations_atmosphere.png')
plt.savefig('coupled_rate_equations_atmosphere.pdf')
plt.show()
