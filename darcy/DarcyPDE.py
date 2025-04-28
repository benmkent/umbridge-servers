from mpi4py import MPI
from petsc4py import PETSc

import numpy as np

from basix.ufl import element, mixed_element
from dolfinx import fem, io, mesh, plot
from dolfinx.fem.petsc import LinearProblem
from ufl import (Measure, SpatialCoordinate, TestFunction, TrialFunction, TestFunctions, TrialFunctions,
                 div, exp, inner, grad, Constant, dot, FacetNormal)

class DarcyPDE:
    """
    DarcyPDE is a class for solving Darcy flow problems using dolfinx
    """
    def __init__(self):
        """
        Constructor for the DarcyPDE class
        """
        return
        
    def setup_function_spaces(self,ndim=2):
        """
        Sets up function spaces for solving a Darcy PDE problem.
        Supports 1D and 2D domains, creating appropriate meshes and spaces.
        Parameters:
            ndim (int): Spatial dimension (1 or 2). Defaults to 2.
        Attributes: V, Q, P, domain, ndim.
        """
        k = 1
        if ndim == 2:
            domain = mesh.create_rectangle(MPI.COMM_WORLD, [[-1.0, -1.0],[1.0,1.0]], [8, 8], cell_type=mesh.CellType.quadrilateral)
            Q_el = element("RT", domain.basix_cell(), k)
            P_el = element("DG", domain.basix_cell(), k - 1)
        elif ndim == 1:
            domain = mesh.create_interval(MPI.COMM_WORLD, 32, [-1.0, 1.0])
            Q_el = element("CG", domain.basix_cell(), k)
            P_el = element("DG", domain.basix_cell(), k - 1)
        else:
            raise ValueError("ndim must be 1 or 2")
        
        V_el = mixed_element([Q_el, P_el])
        V = fem.functionspace(domain, V_el)
        Q, _ = V.sub(0).collapse()
        P, _ = V.sub(1).collapse()

        self.V = V
        self.Q = Q
        self.P = P
        self.domain = domain
        self.ndim = ndim

    def setup_permeability(self, y, basis, n_basis):
        """
        Sets up the permeability function for the Darcy PDE problem.

        Args:
            y (array-like): Permeability values for each patch.
            basis (str): Basis type, only "patch" is supported.
            n_basis (tuple): Number of patches in each spatial dimension.

        Outputs:
            Interpolates and assigns permeability as `self.permeability`.
        """
        permeability = fem.Function(self.P)
                   
        if basis == "patch":
            def z(x):
                return np.array([patchwork_perm(xi, y, self.ndim, n_basis) for xi in x.T])
        else:
            ValueError("Unknown input basis")

        # Interpolate the permeability function
        permeability.interpolate(z)
        permeability.name = "permeability"

        # Assign the computed permeability function to the class
        self.permeability = permeability

    def setup_problem(self):
        """
        Sets up the finite element problem for solving Darcy's equation.
        """
        (u, p) = TrialFunctions(self.V)
        (v, q) = TestFunctions(self.V)

        x = SpatialCoordinate(self.domain)
        # f = 1.0 * exp(-((x[0] - 0.5) * (x[0] - 0.5)))
        f = fem.Function(self.P)
        f.interpolate(lambda x:  0.0+0.0* x[0])

        dx = Measure("dx", self.domain)

        # Define boundary conditions
        LEFT = 1
        fdim = self.domain.topology.dim - 1
        def left_boundary(x):
            return np.isclose(x[0], -1.0)
        def right_boundary(x):
            return np.isclose(x[0], 1.0)
        left_facets = mesh.locate_entities_boundary(self.domain, fdim, left_boundary)
        right_facets = mesh.locate_entities_boundary(self.domain, fdim, right_boundary)
        facet_values = np.full(left_facets.size, LEFT, dtype=np.int32)
        facet_tag = mesh.meshtags(self.domain, fdim, left_facets, facet_values)
        ds = Measure("ds", domain=self.domain, subdomain_data=facet_tag)

        # Define bilinear and linear forms
        x = SpatialCoordinate(self.domain)
        n = FacetNormal(self.domain)
        if self.ndim == 2:
            a = 1/self.permeability * inner(u, v) * dx - inner(p, div(v)) * dx - inner(q, div(u)) * dx
            L = inner(f, q) * dx  - inner((1-x[1]**4) * v, n) * ds(LEFT)
        elif self.ndim == 1:
            a = 1/self.permeability * inner(u, v) * dx - inner(p, grad(v)[0]) * dx - inner(q, grad(u)[0]) * dx
            L = inner(f, q) * dx  + v * ds(LEFT)
        else:
            raise ValueError("ndim must be 1 or 2")
        
        # Construct the linear problem
        self.problem = LinearProblem(a, L, bcs=[], petsc_options={"ksp_type": "gmres", "pc_type": "none"})

    def solve(self):
        """
        Solves the Darcy PDE problem and splits the solution into velocity and pressure components.
        """
        w_h = self.problem.solve()
        self.w_h = w_h
        self.u_h, self.p_h = w_h.split()

    def write_solution(self):
        """
        Write the solution to XDMF files for velocity and pressure.
        """
        u_h, p_h = self.w_h.split()
        u_h.name = "u_h"
        p_h.name = "p_h"

        # Export to XDMF
        if self.ndim == 2:
            el = element("CG", self.domain.basix_cell(), 1, shape=(2,))
            u_h_out = fem.Function(fem.functionspace(self.domain, el))
            u_h_out.interpolate(u_h)
        elif self.ndim == 1:
            u_h_out = u_h
        else:
            raise ValueError("ndim must be 1 or 2")
        u_h_out.name = "u"

        el = element("CG", self.domain.basix_cell(), 1)
        p_h_out = fem.Function(fem.functionspace(self.domain, el))
        p_h_out.interpolate(p_h)
        p_h_out.name = "p"

        with io.XDMFFile(self.domain.comm, "darcy/u.xdmf", "w") as file:
            file.write_mesh(self.domain)
            file.write_function(u_h_out)

        with io.XDMFFile(self.domain.comm, "darcy/p.xdmf", "w") as file:
            file.write_mesh(self.domain)
            file.write_function(p_h_out)

        with io.XDMFFile(self.domain.comm, "darcy/permeability.xdmf", "w") as file:
            file.write_mesh(self.domain)
            file.write_function(self.permeability)

    def get_u(self):
        """
        Retrieve the coordinates and values of the velocity field.

        Returns:
            - coords (numpy.ndarray): The coordinates of the velocity field.
            - values (numpy.ndarray): The velocity values at the corresponding coordinates.
        """
        f = self.u_h

        # Interpolate RT field into CG1 vector space
        # Does this make sense?
        if self.ndim == 2:
            V_cg = fem.functionspace(self.domain, element("CG", self.domain.basix_cell(), 1, shape=(2,)))

            u_interp = fem.Function(V_cg)
            u_interp.interpolate(f)

            # Get coordinates and values for velocity
            coords = V_cg.tabulate_dof_coordinates()
            values = u_interp.x.array[:].reshape(-1, 2)
        elif self.ndim == 1:
            fs, dofmap = f.function_space.collapse()
            values = f.x.array[dofmap]
            coords = fs.tabulate_dof_coordinates()
        else:
            raise ValueError("ndim must be 1 or 2")

        return coords, values
    
    def get_p(self):
        """
        Retrieve the coordinates and values of the pressure field.

        Returns:
            - coords (numpy.ndarray): The coordinates of the pressure field.
            - values (numpy.ndarray): The pressure values at the corresponding coordinates.
        """
        f = self.p_h

        fs, dofmap = f.function_space.collapse()
        values = f.x.array[dofmap]

        coords = fs.tabulate_dof_coordinates()

        return coords, values

# Define the patchwork grid
def patchwork_perm(x, y, ndim, n_basis):
    """
    Generate a permeability field based on a patchwork grid.
    Parameters:
    - x (ndarray): Coordinates of shape (2, N), where x[0] and x[1] are x and y coordinates.
    - y (ndarray): Permeability values for each grid cell.
    - ndim (int): Number of dimensions (1 or 2).
    - n_basis (tuple): Number of grid cells along each dimension.
    Returns:
    - ndarray: Permeability field corresponding to the input coordinates.
    """
    x0 = x[0]
    grid_size_x = n_basis[0]
    xmin, xmax = -1.0, 1.0
    patch_width = (xmax - xmin) / grid_size_x
    grid_size_y = 1
    patch_height = 0
    if ndim == 2:
        x1 = x[1]
        grid_size_y = n_basis[1]
        ymin, ymax = -1.0, 1.0
        patch_height = (ymax - ymin) / grid_size_y
    
    # Define permeability values for each grid cell (can be modified)
    perm_values = y
    
    # Initialize the permeability field
    z = np.zeros_like(x0)
    
    # Iterate over the grid to assign permeability
    for i in range(grid_size_x):
        for j in range(grid_size_y):
            # Determine which patch the current point belongs to
            mask_x = (x0 >= xmin + i * patch_width) & (x0 < xmin + (i + 1) * patch_width)
            if ndim == 2:
                mask_y = (x1 >= ymin + j * patch_height) & (x1 < ymin + (j + 1) * patch_height)
            else:
                mask_y = np.ones_like(x0, dtype=bool)
            
            # Set permeability value for this patch
            mask = mask_x & mask_y
            z[mask] = perm_values[i * grid_size_y + j]  
    return z