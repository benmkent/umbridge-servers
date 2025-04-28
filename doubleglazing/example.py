from doubleglazingpde import DoubleGlazingPDE
import numpy as np


# Define structured grid for output
n = 401  # Grid resolution
x_vals = np.linspace(-1, 1, n)
y_vals = np.linspace(-1, 1, n)

dt=1e-5

model = DoubleGlazingPDE(400, 1)
# Set up cookie problem
model.setupProblem('parameter', [0.01,1,-1,-1,-1,-1], 1, varcoeffs=0,bcrate=0.1,advection=1)
# Solve linear system with preconditioning pc and solver tolerance tol
model.solveTimeSimple(dt,10.0)
# Sample u(x, y) on the structured grid
data = []
for y in y_vals[::-1]:  # Reverse order for correct PGFPlots orientation
    for x in x_vals:
        point = np.array([x, y])
        u_val = model.u(point)  # Evaluate function
        data.append(f"{x} {y} {u_val}\n")
    data.append("\n")  # Blank line to separate rows
# Save as a .dat file
with open("surface-1.dat", "w") as f:
    f.writelines(data)

model = DoubleGlazingPDE(400, 1)
# Set up cookie problem
model.setupProblem('parameter', [0.1,1,0,0,0,0], 1, varcoeffs=0,bcrate=0.1,advection=1)
# Solve linear system with preconditioning pc and solver tolerance tol
model.solveTimeSimple(dt,10.0)
# Sample u(x, y) on the structured grid
data = []
for y in y_vals[::-1]:  # Reverse order for correct PGFPlots orientation
    for x in x_vals:
        point = np.array([x, y])
        u_val = model.u(point)  # Evaluate function
        data.append(f"{x} {y} {u_val}\n")
    data.append("\n")  # Blank line to separate rows
# Save as a .dat file
with open("surface0.dat", "w") as f:
    f.writelines(data)

model = DoubleGlazingPDE(400, 1)
# Set up cookie problem
model.setupProblem('parameter', [0.1*1.9,1,1,1,1,1], 1, varcoeffs=0,bcrate=0.1,advection=1)
# Solve linear system with preconditioning pc and solver tolerance tol
model.solveTimeSimple(dt,10.0)
# Sample u(x, y) on the structured grid
data = []
for y in y_vals[::-1]:  # Reverse order for correct PGFPlots orientation
    for x in x_vals:
        point = np.array([x, y])
        u_val = model.u(point)  # Evaluate function
        data.append(f"{x} {y} {u_val}\n")
    data.append("\n")  # Blank line to separate rows
# Save as a .dat file
with open("surface1.dat", "w") as f:
    f.writelines(data)