import umbridge
from ellipticpde import EllipticPDE
from dolfin import *


class CookieBenchmark(umbridge.Model):
    """
    Class representing the Cookie Benchmark model for elliptic PDEs.
    """

    def __init__(self):
        """
        Initialize the CookieBenchmark model.
        """
        super().__init__("cookiebenchmark")

    def get_input_sizes(self, config):
        """
        Get the sizes of the input parameters.

        Parameters:
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the size of the input parameter vector.
        """
        return [8]

    def get_output_sizes(self, config):
        """
        Get the sizes of the output parameters.

        Parameters:
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the size of the output parameter vector.
        """
        return [1]

    def __call__(self, parameters, config):
        """
        Evaluate the Cookie Benchmark model.

        Parameters:
        parameters (list): List of input parameters for the model.
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the quantity of interest.
        """
        config = verifyConfig(config)
        # Initialize PDE model
        model = EllipticPDE(config['N'])
        # Set up cookie problem
        model.setupProblem(
            'cookie', parameters[0], config['quad_degree'], config['coeffs'])
        # Solve linear system with preconditioning pc and solver tolerance tol
        model.solve(config['pc'], config['tol'])
        # Compute quantity of interest (QoI) on solution
        integral = model.computebenchmarkqoi()
        # Return QoI
        return [[integral]]

    def supports_evaluate(self):
        """
        Check if the model supports evaluation.

        Returns:
        bool: True, indicating that the model supports evaluation.
        """
        return True

    def supports_gradient(self):
        """
        Check if the model supports gradient computation.

        Returns:
        bool: False, indicating that the model does not support gradient computation.
        """
        return False


def verifyConfig(config):
    if config is None:
        config = {}

    if ~('N' in config):
        config['N'] = 500
    if ~('quad_degree' in config):
        config['quad_degree'] = 5
    if ~('coeffs' in config):
        config['coeffs'] = None
    if ~('pc' in config):
        config['pc'] = "none"
    if ~('tol' in config):
        config['tol'] = "LU"
    return config


cookiebenchmark = CookieBenchmark()

umbridge.serve_models([cookiebenchmark], 4242)
