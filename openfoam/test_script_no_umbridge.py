import os
import csv
import numpy as np
from postprocess_openfoam import extract_reattachment_point

def main():
    x_a = 0.0
    x_b = 1.0
    n = 2**1+1
    print(n)
    x = [np.cos(float(x) * np.pi/(n-1)) for x in range(0,n)]
    print(x)
    x = [(x_b-x_a)/2*x + (x_a+x_b)/2 for x in x]
    print(x)

    print("Define config dictionary")
    config = {'Fidelity':1}
    print(config)

    # Specify the filename
    filename = "reattachment_points.csv"
    reattachment_pt = np.zeros([n,1])
    # Open the file in write mode
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)

        for (ii,z) in enumerate(x):
            print(str(ii))
            u = 23.4*(0.1 + (1-0.1)*z)
            v = 34.6
            x = run_openfoam_simulation([[u,v]], config)
            print("Simulation for "+str(u)+" jet, " + str(v) + " inflow")
            # x = run_openfoam_simulation_empty([[u,v]], config)
            print("Reattachment point " + str(x[0][0]))

            reattachment_pt[ii] = x[0][0]

            # Write each pair of elements to the CSV file
            writer.writerow([u, v, x[0][0]])

def run_openfoam_simulation_empty(parameters, config):
    return [[0]]

def run_openfoam_simulation(parameters, config):
    output_dir = './outputdata'

    # Decide on fidelity
    if config['Fidelity'] == 0:
        casefile = "./NASA_hump_data_coarse4"
    elif config['Fidelity'] == 1:
        casefile = "./NASA_hump_data_coarse3"
    elif config['Fidelity'] == 2:
        casefile = "./NASA_hump_data_coarse2"
    elif config['Fidelity'] == 3:
        casefile = "./NASA_hump_data_coarse1"
    elif config['Fidelity'] == 4:
        casefile = "./NASA_hump_data_baseline"
    else:
        AssertionError("Unknown config")

    if 'res_tol' not in config:
        config['res_tol'] = 1e-10
    if 'abs_tol' not in config:
        config['abs_tol'] = 1e-10

    # Copy folder to use as realisation
    tempcasefile = "./caserealisation"
    os.system('cp -r ' + casefile + ' ' + tempcasefile)

    # For realisation assign parameters
    input_file = tempcasefile+"/system/controlDict"
    output_file = input_file
    replacement_value_jet = str(parameters[0][0])
    replace_jet_mag(input_file, output_file, replacement_value_jet)
    replacement_value_inflow = str(parameters[0][1])
    replace_inflow_mag(input_file, output_file, replacement_value_inflow)
    # replacement_value = str(config['final_time'])
    # replace_final_time(input_file, output_file, replacement_value)

    input_file = tempcasefile+"/system/fvSolution"
    output_file = input_file
    replacement_value = str(config['res_tol'])
    replace_res_tol(input_file, output_file, replacement_value)
    replacement_value = str(config['abs_tol'])
    replace_abs_tol(input_file, output_file, replacement_value)

    # Set up boundary conditions
    print("Enforcing boundary conditions jetNasaHump")
    os.system('openfoam2406 jetNasaHump -case '+tempcasefile)

    # Run simple foam
    print("Run simplefoam")
    os.system('openfoam2406 simpleFoam -case '+tempcasefile)

    # Extract quantity of interest (reattachment point)
    print("Extract reattachment point")
    (x, X, Tx) = extract_reattachment_point(tempcasefile, 5000)
    print("Reattachment point: " + str(x))

    # Step 3: Stack the vectors as columns
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    wall_shear = np.column_stack((X, Tx))
    np.savetxt(output_dir+'/wallshearJet' + str(replacement_value_jet) + 'Inflow' +
               str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + '.csv', wall_shear, fmt='%f')

    # Clean up
    print("Clean up temporary case file")
    #os.system('rm -r ' + tempcasefile)

    return [[x]]


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

    print(f"Replaced 'INFLOW_MAG' with '{replacement_value}' in '{
          input_file}'. Output saved to '{output_file}'.")


if __name__ == "__main__":
    main()