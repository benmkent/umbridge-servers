import umbridge
import os

class TestModel(umbridge.Model):

    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [1]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        # arguments = " ".join([str(x) for x in parameters[0]]);

        # num_threads = str(config.get("NumThreads",1));
        # basis_degree = str(config.get("BasisDegree",4));
        # fidelity = str(config.get("Fidelity",2));

        # arguments = arguments + " " + basis_degree + " " + fidelity

        # System call
        os.system('openfoam2312 simpleFoam -case ./2DWMH/pitzDailyModified/')

        # Read second line of output file
        #with open('/poisson_lorenzo_results.dat', 'r') as f:
        #    line = f.readline() # Read first line

        return [[float(0)]]


    def supports_evaluate(self):
        return True

testmodel = TestModel()


umbridge.serve_models([testmodel], 4242)