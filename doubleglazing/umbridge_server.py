import umbridge
import json
import os
from doubleglazingpde import DoubleGlazingPDE
from dolfin import *


class DoubleGlazingForward(umbridge.Model):
    """
    Class representing the forward DG elliptic PDE model.
    """

    def __init__(self):
        """
        Initialize the DoubleGlazingForward model.
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
        config = verifyConfig(config)
        n = 1 + config['advection']
        return [n]

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
        Evaluate the DoubleGlazing model.

        Parameters:
        parameters (list): List of input parameters for the model.
        config (dict): Configuration dictionary.

        Returns:
        list: List containing the quantity of interest.
        """
        # Fill missing entries in config with default values
        config = verifyConfig(config)
        output_dir = './outputdata'
        configstring = json.dumps(config, sort_keys=True)+str(parameters[0])
        configstring=str(hash(configstring))
        filename = output_dir+'/'+configstring+'.csv'
        if os.path.exists(filename):
            print("Opening file "+configstring+"\n")
            with open(filename, 'r') as file:
                q = float(file.read())
        else:
            # Initialize PDE model
            model = DoubleGlazingPDE(config['N'], config['BasisDegree'])
            # Set up cookie problem
            model.setupProblem('parameter', parameters[0], config['quad_degree'], varcoeffs=config['diffzero'])
            # Solve linear system with preconditioning pc and solver tolerance tol
            model.solve(config['directsolver'], config['pc'], config['tol'])
            q = model.computepointqoi()
            with open(filename, 'w') as file:
                file.write(str(q))
                print("Writing file "+configstring+"\n")
        # Compute quantity of interest (QoI) on solution
        pointval = q

        # Return QoI
        return [[pointval]]

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

class DoubleGlazingTime(umbridge.Model):
    """
    A model for parabolic formulation of the DG model.

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
        super().__init__("forwardparabolic")

    def get_input_sizes(self, config):
        """
        Returns the input sizes for the model.

        Args:
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the input sizes.
        """
        config = verifyConfig(config)

        n = 1 + config['advection']
        return [n]

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
        Evaluates the parabolic DG model.

        Args:
            parameters (list): A list of parameters.
            config (dict): A dictionary containing configuration parameters.

        Returns:
            list: A list containing the computed integral.
        """
        # Verfiy config and fills in empty keys.
        config = verifyConfig(config)
        print(config)
        output_dir = './outputdata'
        configstring = json.dumps(config, sort_keys=True)+str(parameters[0])
        configstring=str(hash(configstring+'parabolic'))
        filename = output_dir+'/'+configstring+'.csv'
        if os.path.exists(filename):
            print("Opening file "+configstring+"\n")
            with open(filename, 'r') as file:
                q = float(file.read())
        else:
            # Initialize PDE model
            model = DoubleGlazingPDE(config['N'], config['BasisDegree'])
            # Set up cookie problem
            model.setupProblem('parameter', parameters[0], config['quad_degree'], varcoeffs=config['diffzero'], advection=config['advection'], bcrate=config['bcrate'], point=config['qoipoint'])
            # Use the custom TR-AB2 solver. Optional solveTime function uses built in PETSC TS solver.
            u = model.solveTimeSimple(config['letol'],config['T'])
            # Compute QoI at finalTime T
            pointval = model.computepointqoi()

            with open(filename, 'w') as file:
                file.write(str(q))
                print("Writing file "+configstring+"\n")

        # Return QoI
        return [[pointval]]

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

    # Use direct solver by default
    if 'directsolver' not in config:
        config['directsolver'] = 1

    # Use 400x400 mesh by default
    if 'Fidelity' not in config:
        config['N'] = 400
    else:
        config['N'] = int(100 * config['Fidelity'])

    # Use Q1 approximation by default
    if 'BasisDegree' not in config:
        config['BasisDegree'] = 1
    
    # Use degree 8 quadrature by default
    if 'quad_degree' not in config:
        config['quad_degree'] = 8
    
    # Use background diffusion 1.0 by default
    if 'diffzero' not in config:
        config['diffzero'] = [1.0]
    
    # Use no preconditioning by default
    if 'pc' not in config:
        config['pc'] = "none"

    # Use GM-RES tol 1e-4 by default
    if 'tol' not in config:
        config['tol'] = 1e-4
    
    # Use local timestepping error tolerance 1e-4 by default
    if 'letol' not in config:
        config['letol'] = 1e-4
    
    # Use a finalTime T = 10.0 by default
    if 'T' not in config:
        config['T'] = 10.0

    if 'advection' not in config:
        config['advection'] = 1 + 4

    if 'bcrate' not in config:
        config['bcrate'] = 0

    if 'qoipoint' not in config:
        config['qoipoint'] = [0.5,-0.5]

    print(config)
    return config

# Initialise UM-BRIDGE models
dgelliptic = DoubleGlazingForward()
dgparabolic = DoubleGlazingTime()

# Start UM-BRIDGE server
umbridge.serve_models([dgelliptic,dgparabolic], 4242)
