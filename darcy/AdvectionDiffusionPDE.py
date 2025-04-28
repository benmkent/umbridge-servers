from mpi4py import MPI
from petsc4py import PETSc

import numpy as np
from numpy.polynomial.chebyshev import chebval2d

from basix.ufl import element, mixed_element
from dolfinx import fem, io, mesh, plot
from dolfinx.fem.petsc import LinearProblem
from ufl import (Measure, SpatialCoordinate, TestFunction, TrialFunction, TestFunctions, TrialFunctions,
                 div, exp, inner, grad, Constant, dot, FacetNormal)

class AdvDiffPDE:
    """
    Define and solve an advection-diffusion PDE using FEniCSx.
    """

    def __init__(self):
        """
        Initialize the AdvDiffPDE class
        """
        return
    
    def setup_function_spaces(self, ndim=2):
        """
        Set up the function spaces for the PDE.

        Parameters:
        ndim (int): Dimension of the problem (1D or 2D).
        """
        k = 1  # Polynomial degree for finite element space
        if ndim == 2:
            # Create a 2D rectangular mesh
            domain = mesh.create_rectangle(MPI.COMM_WORLD, [[-1.0, -1.0], [1.0, 1.0]], [8, 8], cell_type=mesh.CellType.quadrilateral)
            Q_el = element("CG", domain.basix_cell(), k)
        elif ndim == 1:
            # Create a 1D interval mesh
            domain = mesh.create_interval(MPI.COMM_WORLD, 32, [-1.0, 1.0])
            Q_el = element("CG", domain.basix_cell(), k)
        else:
            raise ValueError("ndim must be 1 or 2")  # Error for unsupported dimensions
        
        V = fem.functionspace(domain, Q_el)
        self.V = V
        self.domain = domain

    def setup_problem(self, ux, uy, epsilon):
        """
        Set up the advection-diffusion problem.

        Parameters:
        ux (callable): x-component of the advection field.
        uy (callable): y-component of the advection field.
        epsilon (float): Diffusion coefficient.
        """
        u = TrialFunction(self.V)  # Trial function
        v = TestFunction(self.V)  # Test function

        x = SpatialCoordinate(self.domain)
        dx = Measure("dx", self.domain)

        # Set up boundary condition
        LEFT = 1
        fdim = self.domain.topology.dim - 1
        def left_boundary(x):
            return np.isclose(x[0], -1.0)
        def right_boundary(x):
            return np.isclose(x[0], 1.0)
        left_facets = mesh.locate_entities_boundary(self.domain, fdim, left_boundary)
        left_dofs = fem.locate_dofs_topological(self.V, fdim, left_facets)

        # Define Dirichlet boundary condition on the left wall
        leftwall = fem.Function(self.V)
        leftwall.interpolate(lambda x: 1 - x[1]**4)
        bc = fem.dirichletbc(leftwall, left_dofs)
        facet_values = np.full(left_facets.size, LEFT, dtype=np.int32)
        facet_tag = mesh.meshtags(self.domain, fdim, left_facets, facet_values)

        # Define source term
        f = fem.Function(self.V)
        f.interpolate(lambda x: 0.0 + 0.0 * x[0])

        # Create advection field from ux, uy
        w = create_advection_field(ux, uy, self.domain)
        self.advfn = w

        # Write advection field to file
        with io.XDMFFile(self.domain.comm, "darcy/wind_advdiff.xdmf", "w") as file:
            file.write_mesh(self.domain)
            file.write_function(self.advfn)

        # Define variational problem
        a = epsilon * inner(grad(u), grad(v)) * dx + inner(w, grad(u)) * v * dx
        L = f * v * dx

        # SUPG
        h = 2.0/np.sqrt(domain.topology.index_map(0).size) # Approximate
        tau = h / dot(w, w)
        a += tau * inner(inner(w, grad(u)) * grad(v), dx)
        
        # Create linear problem
        self.problem = LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "gmres", "pc_type": "none"})

    def solve(self):
        """
        Solve the weak formulation.
        """
        self.u_h = self.problem.solve()

    def write_solution(self):
        """
        Write the solution to an XDMF file.
        """
        self.u_h.name = "u_h_advdiff"  # Name the solution

        with io.XDMFFile(self.domain.comm, "darcy/u_advdiff.xdmf", "w") as file:
            file.write_mesh(self.domain)
            file.write_function(self.u_h)

    def get_u(self):
        """
        Get the solution values and coordinates.

        Returns:
        tuple: Coordinates and solution values.
        """
        u = self.u_h  # Solution function
        values = u.x.array[:]  # Extract solution values
        coords = self.V.tabulate_dof_coordinates()  # Get coordinates of degrees of freedom

        return coords, values
    
def create_advection_field(ux, uy, domain):
    """
    Create the advection field as a vector function.

    Parameters:
    ux (callable): x-component of the advection field.
    uy (callable): y-component of the advection field.
    domain (dolfinx.mesh.Mesh): Mesh domain.

    Returns:
    dolfinx.fem.Function: Advection field.
    """
    V_vec = fem.functionspace(domain, element("CG", domain.basix_cell(), 1, shape=(2,)))  # Vector function space

    advection = fem.Function(V_vec)  # Advection field function

    def adv_expression(x):
        return np.vstack((
            ux(x[0], x[1]),  # x-component
            uy(x[0], x[1]),  # y-component
        ))

    advection.interpolate(adv_expression)  # Interpolate advection field
    return advection