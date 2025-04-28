from DarcyPDE import DarcyPDE
from numpy.polynomial.chebyshev import chebfit, chebval
import matplotlib.pyplot as plt
import numpy as np

d = DarcyPDE()
d.setup_function_spaces(ndim=1)
d.setup_permiability([0.1,10,1,0.5],"patch",[4])
d.setup_problem()
d.solve()
d.write_solution()
u = d.get_u()
p = d.get_p()

coords, values = p
# interpolant = interp1d(coords[:,0], values, kind='linear', fill_value="extrapolate")

# Fit a Chebyshev series to the data
degree = 10  # Degree of the Chebyshev polynomial
coeffs = chebfit(coords[:,0], values, degree)

# Plot the original data and the Chebyshev fit
x = np.linspace(-1, 1, 100)
y_cheb = chebval(x, coeffs)
plt.plot(x, y_cheb, label=f'Chebyshev fit (degree={degree})', linestyle='dashed')
plt.plot(coords, values, label=f'Data', linestyle='solid')
plt.legend()
plt.show()