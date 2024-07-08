import umbridge
from ellipticpde import EllipticPDE
from dolfin import *


class CookieForward(umbridge.Model):
    """
    Class representing the forward cookie model of the elliptic PDEs.
    """

    def __init__(self):
        """
        Initialize the CookieForward model.
        """
        super().__init__("forward")

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
        Evaluate the CookieForward model.

        Parameters:
        parameters (list): List of input parameters for the model.
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the quantity of interest.
        """
        # Fill missing entries in config with default values
        config = verifyConfig(config)
        # Initialize PDE model
        model = EllipticPDE(config['N'], config['BasisDegree'], config['quad_degree'])
        # Set up cookie problem
        model.setupProblem('cookie', parameters[0], config['quad_degree'], config['coeffs'])
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


class CookieBenchmark(umbridge.Model):
    """
    Class representing the Cookie Benchmark model for elliptic PDEs.
    """

    def __init__(self):
        """
        Initialize the CookieBenchmark model.
        """
        super().__init__("benchmark")

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
        # Use default values by ensuring dictionary is empty
        config={}
        config = verifyConfig(config)
        # Initialize PDE model
        model = EllipticPDE(config['N'], config['BasisDegree'], config['quad_degree'])
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

class CookieTime(umbridge.Model):
    """
    A model for Cookie Time.

    Inherits from umbridge.Model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initializes the CookieTime object.

        Args:
            None

        Returns:
            None
        """
        super().__init__("cookietime")

    def get_input_sizes(self, config):
        """
        Returns the input sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the input sizes.
        """
        return [8]

    def get_output_sizes(self, config):
        """
        Returns the output sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the output sizes.
        """
        return [1]
    
    def __call__(self, parameters, config):
        """
        Invokes the model to perform computations.

        Args:
            parameters (list): A list of parameters.
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the computed integral.
        """
        config = verifyConfig(config)

        model = EllipticPDE(config['N'], config['BasisDegree'], config['quad_degree'])
        model.setupProblem('cookie', parameters[0], advection=True)
        # u = model.solveTime(config['tol'],config['T'])
        u = model.solveTimeSimple(config['letol'],config['T'])
        integral = model.computebenchmarkqoi()
        # model.writeSln("outputFinal")
        return [[integral]]

    def supports_evaluate(self):
        """
        Checks if the model supports evaluation.

        Returns:
            bool: True if the model supports evaluation, False otherwise.
        """
        return True

    def supports_gradient(self):
        """
        Checks if the model supports gradient computation.

        Returns:
            bool: True if the model supports gradient computation, False otherwise.
        """
        return False

class CookieTimeBenchmark(umbridge.Model):
    """
    A benchmark model for Cookie Time.

    Inherits from umbridge.Model.

    Attributes:
        None
    """

    def __init__(self):
        """
        Initializes the CookieTimeBenchmark object.

        Args:
            None

        Returns:
            None
        """
        super().__init__("cookietimebenchmark")

    def get_input_sizes(self, config):
        """
        Returns the input sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the input sizes.
        """
        return [8]

    def get_output_sizes(self, config):
        """
        Returns the output sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the output sizes.
        """
        return [1]
    
    def __call__(self, parameters, config):
        """
        Invokes the model to perform computations.

        Args:
            parameters (list): A list of parameters.
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the computed integral.
        """
        # Use default values by ensuring dictionary is empty
        config={}
        config = verifyConfig(config)

        model = EllipticPDE(config['N'], config['BasisDegree'], config['quad_degree'])
        model.setupProblem('cookie', parameters[0], advection=True)
        # u = model.solveTime(config['tol'],config['T'])
        u = model.solveTimeSimple(config['letol'],config['T'])
        integral = model.computebenchmarkqoi()
        # model.writeSln("outputFinal")
        return [[integral]]

    def supports_evaluate(self):
        """
        Checks if the model supports evaluation.

        Returns:
            bool: True if the model supports evaluation, False otherwise.
        """
        return True

    def supports_gradient(self):
        """
        Checks if the model supports gradient computation.

        Returns:
            bool: True if the model supports gradient computation, False otherwise.
        """
        return False


def verifyConfig(config):
    if config is None:
        config = {}

    if 'Fidelity' not in config:
        config['N'] = 400
    else:
        config['N'] = 100 * config['Fidelity']

    if 'BasisDegree' not in config:
        config['BasisDegree'] = 1
    if 'quad_degree' not in config:
        config['quad_degree'] = 8
    if 'coeffs' not in config:
        config['coeffs'] = None
    if 'pc' not in config:
        config['pc'] = "none"
    if 'tol' not in config:
        config['tol'] = "LU"
    if 'letol' not in config:
        config['letol'] = 1e-4
    if 'T' not in config:
        config['T'] = 10.0

    print(config)
    return config

# Initialise UM-BRIDGE models
cookieforward = CookieForward()
cookiebenchmark= CookieBenchmark()
cookietimebenchmark = CookieTimeBenchmark()
cookietime = CookieTime()

# Start UM-BRIDGE server
umbridge.serve_models([cookieforward,cookiebenchmark,cookietime,cookietimebenchmark], 4242)
