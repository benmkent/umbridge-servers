from AdvectionDiffusionPDE import AdvDiffPDE
from numpy.polynomial.chebyshev import chebfit, chebval
import matplotlib.pyplot as plt
import numpy as np

import pyvista
from dolfinx import plot

d = AdvDiffPDE()
d.setup_function_spaces(ndim=2)
ux = lambda x,y: (1-x)
uy = lambda x,y: (1-y)
d.setup_problem(ux,uy, 0.1)
d.solve()
d.write_solution()
u = d.get_u()
print(u)