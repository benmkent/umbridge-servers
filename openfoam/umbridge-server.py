import umbridge
import os
import sys
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
            config['res_tol'] = 1e-10

        # Copy folder to use as realisation
        tempcasefile = "./caserealisation"
        os.system('cp -r ' + casefile + ' ' + tempcasefile)

        # For realisation assign parameters
        input_file = tempcasefile+"/system/controlDict"
        output_file = input_file
        replacement_value = str(parameters[0][0])
        replace_jet_mag(input_file, output_file, replacement_value)

        input_file = tempcasefile+"/system/fvSolution"
        output_file = input_file
        replacement_value = str(config['res_tol'])
        replace_jet_mag(input_file, output_file, replacement_value)

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
        os.system('rm -r ' + tempcasefile)

        return [[x]]

    def supports_evaluate(self):
        return True


testmodel = TestModel()


def replace_jet_mag(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'JET_MAG' with 'replacement_value'
    modified_data = file_data.replace('JET_MAG', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'JET_MAG' with '{replacement_value}' in '{
          input_file}'. Output saved to '{output_file}'.")


def replace_res_tol(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'JET_MAG' with 'replacement_value'
    modified_data = file_data.replace('RES_TOL', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'RES_TOL' with '{replacement_value}' in '{
          input_file}'. Output saved to '{output_file}'.")


umbridge.serve_models([testmodel], 4242)
