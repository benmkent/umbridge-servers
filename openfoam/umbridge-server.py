import umbridge
import os
from postprocess_openfoam import extract_reattachment_point

class TestModel(umbridge.Model):

    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        # Decide on fidelity
        if config['Fidelity'] == 1:
            casefile = "./NASA_hump_data_coarse3"
        elif config['Fidelity'] == 2:
            casefile = "./NASA_hump_data_coarse2"        
        elif config['Fidelity'] == 3:
            casefile = "./NASA_hump_data_coarse1"
        elif config['Fidelity'] == 4:
            casefile = "./NASA_hump_data_coarse1"
        else:
            AssertionError("Unknown config")

        if 'res_tol' not in config:
            config['res_tol'] = 1e-10;

        # Copy folder to use as realisation
        tempcasefile = "./caserealisation"
        os.system('cp -r ' + casefile + ' ' + tempcasefile)

        # For realisation assign parameters
        input_file= tempcasefile+"/system/controlDict"
        output_file=input_file
        replacement_value=str(parameters[0][0])

        # Use sed to replace $JET_MAG with the replacement value
        print("Writing jet mag "+replacement_value)
        sed_command = f"sed 's/JET_MAG/{replacement_value}/g' {input_file} > {output_file}"
        os.system(sed_command)

        input_file= tempcasefile+"/system/fvSolution"
        output_file=input_file
        replacement_value=str(config['res_tol'])

        # Use sed to replace $JET_MAG with the replacement value
        print("Writing residual tol "+replacement_value)
        sed_command = f"sed 's/RES_TOL/{replacement_value}/g' {input_file} > {output_file}"
        os.system(sed_command)

        # Set up boundary conditions
        print("Enforcing boundary conditions jetNasaHump")
        os.system('openfoam2406 jetNasaHump -case '+tempcasefile)

        # Run simple foam
        print("Run simplefoam")
        os.system('openfoam2406 simpleFoam -case '+tempcasefile)

        # Extract quantity of interest (reattachment point)
        print("Extract reattachment point")
        x = extract_reattachment_point(tempcasefile)
        display(x)

        # Clean up
        os.system('rm -r '+ tempcasefile)

        return [[x]]


    def supports_evaluate(self):
        return True

testmodel = TestModel()


umbridge.serve_models([testmodel], 4242)