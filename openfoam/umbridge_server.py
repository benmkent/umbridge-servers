import umbridge
import os
import sys
import csv
import numpy as np
from postprocess_openfoam import extract_reattachment_point, extract_reattachment_point_from_dataseries
from postprocess_openfoam import extract_cf, extract_cf_from_dataseries


class ReattachmentModel(umbridge.Model):

    def __init__(self):
        super().__init__("forwardopenfoam")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        output_dir = './outputdata'
        try:
            # Decide on fidelity
            if config['Fidelity'] == 0:
                casefile = "./NASA_hump_data_coarse4"
                print("Selecting fidelity 0")
            elif config['Fidelity'] == 1:
                casefile = "./NASA_hump_data_coarse3"
                print("Selecting fidelity 1")
            elif config['Fidelity'] == 2:
                casefile = "./NASA_hump_data_coarse2"
                print("Selecting fidelity 2")
            elif config['Fidelity'] == 3:
                casefile = "./NASA_hump_data_coarse1"
                print("Selecting fidelity 3")
            elif config['Fidelity'] == 4:
                casefile = "./NASA_hump_data_baseline"
                print("Selecting fidelity 4")
            else:
                AssertionError("Unknown config")

            if 'res_tol' not in config:
                config['res_tol'] = 1e-10
                res_tol_str = ''
            elif abs(config['res_tol'] - 1e-10) < 1e-14:
                res_tol_str = ''
            else:
                res_tol_str = '_restol' + str(config['res_tol'])
            print("Residual tolerance "+str(config['res_tol']))

            if 'abs_tol' not in config:
                config['abs_tol'] = 1e-10
                abs_tol_str = ''
            elif abs(config['abs_tol'] - 1e-10) < 1e-14:
                abs_tol_str = ''
            else:
                abs_tol_str = '_abstol' + str(config['abs_tol'])
            print("Iterative solver tolerance "+str(config['abs_tol']))

            # Copy folder to use as realisation
            tempcasefile = "./caserealisation"
            os.system('cp -r ' + casefile + ' ' + tempcasefile)
            print("Create temporary case files")

            # For realisation assign parameters
            input_file = tempcasefile+"/system/controlDict"
            output_file = input_file
            replacement_value_jet = str(parameters[0][0])
            replace_jet_mag(input_file, output_file, replacement_value_jet)
            replacement_value_inflow = str(parameters[0][1])
            replace_inflow_mag(input_file, output_file,
                                replacement_value_inflow)
            # replacement_value = str(config['final_time'])
            # replace_final_time(input_file, output_file, replacement_value)

            input_file = tempcasefile+"/system/fvSolution"
            output_file = input_file
            replacement_value = str(config['res_tol'])
            replace_res_tol(input_file, output_file, replacement_value)
            replacement_value = str(config['abs_tol'])
            replace_abs_tol(input_file, output_file, replacement_value)

            replacement_value_jet = f"{float(replacement_value_jet):.12g}"
            replacement_value_inflow = f"{float(replacement_value_inflow):.12g}"

            filename = output_dir+'/wallshearJet'+ str(replacement_value_jet) + 'Inflow' + str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + res_tol_str + abs_tol_str +'.csv'
            filename_console = output_dir+'/console_'+ str(replacement_value_jet) + 'Inflow' + str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + res_tol_str + abs_tol_str +'.log'

            print("Case configured")

            if os.path.exists(filename):
                X = []
                Tx = []

                print("Reuse precomputed wall shear stress datafile")
                with open(filename, mode='r') as file:
                    reader = csv.reader(file, delimiter=' ')
                    for row in reader:
                        if row:  # Check if row is not empty
                            X.append(row[0])  # First column
                            Tx.append(row[1])  # Second column
                x = extract_reattachment_point_from_dataseries(X,Tx)
            else:
                # Set up boundary conditions
                print("Enforcing boundary conditions jetNasaHump")
                os.system('openfoam2406 jetNasaHump -case '+tempcasefile + ' | tee -a ' + filename_console)

                # Run simple foam
                print("Run simplefoam")
                os.system('openfoam2406 simpleFoam -case '+tempcasefile + ' | tee -a ' + filename_console)

                # Extract quantity of interest (reattachment point)
                print("Extract reattachment point")
                (x, X, Tx) = extract_reattachment_point(tempcasefile, 5000)

                # Step 3: Stack the vectors as columns
                wall_shear = np.column_stack((X, Tx))
                np.savetxt(filename, wall_shear, fmt='%f')

            print("Reattachment point: " + str(x))

            # Clean up
            print("Clean up temporary case file")
            os.system('rm -r ' + tempcasefile)
        except Exception as e:
            # Code to handle any exception
            print(f"An error occurred: {e}")
            raise Exception("A generic error occurred.")
        
        return [[x]]

    def supports_evaluate(self):
        return True

class CfModel(umbridge.Model):

    def __init__(self):
        super().__init__("forwardcf")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        if 'ninterp' not in config:
            config['ninterp'] = 1000
        return [config['ninterp']]

    def __call__(self, parameters, config):
        output_dir = './outputdata'
        try:
            # Decide on fidelity
            if config['Fidelity'] == 0:
                casefile = "./NASA_hump_data_coarse4"
                print("Selecting fidelity 0")
            elif config['Fidelity'] == 1:
                casefile = "./NASA_hump_data_coarse3"
                print("Selecting fidelity 1")
            elif config['Fidelity'] == 2:
                casefile = "./NASA_hump_data_coarse2"
                print("Selecting fidelity 2")
            elif config['Fidelity'] == 3:
                casefile = "./NASA_hump_data_coarse1"
                print("Selecting fidelity 3")
            elif config['Fidelity'] == 4:
                casefile = "./NASA_hump_data_baseline"
                print("Selecting fidelity 4")
            else:
                AssertionError("Unknown config")

            if 'res_tol' not in config:
                config['res_tol'] = 1e-10
                res_tol_str = ''
            elif abs(config['res_tol'] - 1e-10) < 1e-14:
                res_tol_str = ''
            else:
                res_tol_str = '_restol' + str(config['res_tol'])
            print("Residual tolerance "+str(config['res_tol']))

            if 'abs_tol' not in config:
                config['abs_tol'] = 1e-10
                abs_tol_str = ''
            elif abs(config['abs_tol'] - 1e-10) < 1e-14:
                abs_tol_str = ''
            else:
                abs_tol_str = '_abstol' + str(config['abs_tol'])
            print("Iterative solver tolerance "+str(config['abs_tol']))

            # Copy folder to use as realisation
            tempcasefile = "./caserealisation"
            os.system('cp -r ' + casefile + ' ' + tempcasefile)
            print("Create temporary case files")

            # For realisation assign parameters
            input_file = tempcasefile+"/system/controlDict"
            output_file = input_file
            replacement_value_jet = str(parameters[0][0])
            replace_jet_mag(input_file, output_file, replacement_value_jet)
            replacement_value_inflow = str(parameters[0][1])
            replace_inflow_mag(input_file, output_file,
                                replacement_value_inflow)
            # replacement_value = str(config['final_time'])
            # replace_final_time(input_file, output_file, replacement_value)

            input_file = tempcasefile+"/system/fvSolution"
            output_file = input_file
            replacement_value = str(config['res_tol'])
            replace_res_tol(input_file, output_file, replacement_value)
            replacement_value = str(config['abs_tol'])
            replace_abs_tol(input_file, output_file, replacement_value)

            replacement_value_jet = f"{float(replacement_value_jet):.12g}"
            replacement_value_inflow = f"{float(replacement_value_inflow):.12g}"

            filename = output_dir+'/wallshearJet'+ str(replacement_value_jet) + 'Inflow' + str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + res_tol_str + abs_tol_str +'.csv'
            filename_console = output_dir+'/console_'+ str(replacement_value_jet) + 'Inflow' + str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + res_tol_str + abs_tol_str +'.log'

            print("Case configured")

            # Interpolation parameters
            if 'ninterp' not in config:
                config['ninterp'] = 1000
            if 'xmin' not in config:
                config['xmin'] = -6.4
            if 'xmax' not in config:
                config['xmax'] = 4.9

            xmin = config['xmin']
            xmax = config['xmax']
            ninterp = config['ninterp']

            if os.path.exists(filename):
                X = []
                Tx = []

                print("Reuse precomputed wall shear stress datafile")
                with open(filename, mode='r') as file:
                    reader = csv.reader(file, delimiter=' ')
                    for row in reader:
                        if row:  # Check if row is not empty
                            X.append(row[0])  # First column
                            Tx.append(row[1])  # Second column
                
                cf = extract_cf_from_dataseries(X,Tx,xmin,xmax,ninterp)
            else:
                # Set up boundary conditions
                print("Enforcing boundary conditions jetNasaHump")
                os.system('openfoam2406 jetNasaHump -case '+tempcasefile + ' | tee -a ' + filename_console)

                # Run simple foam
                print("Run simplefoam")
                os.system('openfoam2406 simpleFoam -case '+tempcasefile + ' | tee -a ' + filename_console)

                # Extract quantity of interest (reattachment point)
                print("Extract reattachment point")
                (cf, X, Tx) = extract_cf(tempcasefile, 5000, xmin, xmax, ninterp)

                # Step 3: Stack the vectors as columns
                wall_shear = np.column_stack((X, Tx))
                np.savetxt(filename, wall_shear, fmt='%f')

            print("Reattachment point: " + str(x))

            # Clean up
            print("Clean up temporary case file")
            os.system('rm -r ' + tempcasefile)
        except Exception as e:
            # Code to handle any exception
            print(f"An error occurred: {e}")
            raise Exception("A generic error occurred.")
        
        return [[cf]]

    def supports_evaluate(self):
        return True


reattachment_model = ReattachmentModel()
cf_model = CfModel()


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

    # Replace all occurrences of 'RES_TOL' with 'replacement_value'
    modified_data = file_data.replace('RES_TOL', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'RES_TOL' with '{replacement_value}' in '{
          input_file}'. Output saved to '{output_file}'.")


def replace_final_time(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'FINAL_TIME' with 'replacement_value'
    modified_data = file_data.replace('FINAL_TIME', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'FINAL_TIME' with '{replacement_value}' in '{
          input_file}'. Output saved to '{output_file}'.")


def replace_abs_tol(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'ABS_TOL' with 'replacement_value'
    modified_data = file_data.replace('ABS_TOL', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'ABS_TOL' with '{replacement_value}' in '{
          input_file}'. Output saved to '{output_file}'.")


def replace_inflow_mag(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'INFLOW_MAG' with 'replacement_value'
    modified_data = file_data.replace('INFLOW_MAG', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'INFLOW_MAG' with '{replacement_value}' in '{input_file}'. Output saved to '{output_file}'.")


umbridge.serve_models([reattachment_model,cf_model], 4242)
