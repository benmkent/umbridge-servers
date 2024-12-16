import umbridge
import os
import sys
import csv
import numpy as np
from postprocess_openfoam import extract_reattachment_point, extract_reattachment_point_from_dataseries
from postprocess_openfoam import extract_cf, extract_cf_from_dataseries
from postprocess_openfoam import extract_cp, extract_cp_from_dataseries
from postprocess_openfoam import extract_yplus, extract_pWall

def configure_case(config):
    # Decide on fidelity
    if config['Fidelity'] == 0:
        casefile = "./NASA_hump_data_coarse4"
        print("Selecting fidelity 0", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 1:
        casefile = "./NASA_hump_data_coarse3"
        print("Selecting fidelity 1", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 2:
        casefile = "./NASA_hump_data_coarse2"
        print("Selecting fidelity 2", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 3:
        casefile = "./NASA_hump_data_coarse1"
        print("Selecting fidelity 3", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 4:
        casefile = "./NASA_hump_data_baseline"
        print("Selecting fidelity 4", file=sys.stdout, flush=True)
    elif config['Fidelity'] == -1:
        casefile = "./NASA_hump_data_coarse5"
        print("Selecting fidelity -1", file=sys.stdout, flush=True)
    else:
        AssertionError("Unknown config")

    if 'res_tol' not in config:
        config['res_tol'] = 1e-10
        res_tol_str = ''
    elif abs(config['res_tol'] - 1e-10) < 1e-14:
        res_tol_str = ''
    else:
        res_tol_str = '_restol' + str(config['res_tol'])
    print("Residual tolerance "+str(config['res_tol']), file=sys.stdout, flush=True)

    if 'abs_tol' not in config:
        config['abs_tol'] = 1e-10
        abs_tol_str = ''
    elif abs(config['abs_tol'] - 1e-10) < 1e-14:
        abs_tol_str = ''
    else:
        abs_tol_str = '_abstol' + str(config['abs_tol'])
    print("Iterative solver tolerance "+str(config['abs_tol']), file=sys.stdout, flush=True)

    # Copy folder to use as realisation
    tempcasefile = "./caserealisation"
    os.system('cp -r ' + casefile + ' ' + tempcasefile)
    print("Create temporary case files", file=sys.stdout, flush=True)

    # For realisation assign parameters
    input_file = tempcasefile+"/system/controlDict"
    output_file = input_file
    replacement_value_jet = str(parameters[0][0])
    replace_jet_mag(input_file, output_file, replacement_value_jet)
    replacement_value_inflow = str(parameters[0][1])
    replace_inflow_mag(input_file, output_file,
                        replacement_value_inflow)
    wpd = 0.0
    replacement_value = str(wpd)
    replace_wallpressuredist(input_file, output_file, replacement_value)
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

    filename = 'Jet'+str(replacement_value_jet) + 'Inflow' + str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + res_tol_str + abs_tol_str +'.csv'
    filename_console = 'console_'+ str(replacement_value_jet) + 'Inflow' + str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + res_tol_str + abs_tol_str +'.log'

    return filename, filename_console

def run_case(tempcasefile, filename_console, filename, parameters):
    output_dir = './outputdata'

    # Set up boundary conditions
    print("Enforcing boundary conditions jetNasaHump", file=sys.stdout, flush=True)
    os.system('openfoam2406 jetNasaHump -case '+tempcasefile + ' | tee -a ' + filename_console)

    # Run simple foam
    print("Run simplefoam", file=sys.stdout, flush=True)
    os.system('openfoam2406 simpleFoam -case '+tempcasefile + ' | tee -a ' + filename_console)

    # Save out wallshear data
    print("Extract reattachment point", file=sys.stdout, flush=True)
    (rp_value, X, Tx) = extract_reattachment_point(tempcasefile, 5000)
    wall_shear = np.column_stack((X, Tx))
    np.savetxt(output_dir+'/wallshear'+filename, wall_shear, fmt='%f')
    
    # Save out nearWallField data
    uinf = parameters[0][1]
    rhoinf = 1.185
    print("Extract nearWallFields Cp", file=sys.stdout, flush=True)
    (cp_value_nwf, X, Tx) = extract_cp(tempcasefile, 5000, rhoinf, uinf)
    cp = np.column_stack((X, Tx))
    np.savetxt(output_dir+'/nearWallField'+filename, cp, fmt='%f')

    print("Extract pressure field at wall", file=sys.stdout, flush=True)
    (X, Tx) = extract_pWall(tempcasefile, 5000)
    pwall= np.column_stack((X, Tx))
    np.savetxt(output_dir+'/wallPressure'+filename, pwall, fmt='%f')

    print("Extract yPlus at wall", file=sys.stdout, flush=True)
    (X,Tx) = extract_yplus(tempcasefile, 5000)
    yp= np.column_stack((X, Tx))
    np.savetxt(output_dir+'/yPlus'+filename, yp, fmt='%f')

    # Clean up
    print("Clean up temporary case file", file=sys.stdout, flush=True)
    os.system('rm -r ' + tempcasefile)

class ReattachmentModel(umbridge.Model):

    def __init__(self):
        super().__init__("forwardopenfoam")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        try:
            filename, filename_console = configure_case(config)
            print("Case configured", file=sys.stdout, flush=True)

            if not os.path.exists(output_dir+'/wallshear'+filename):
                run_case(tempcasefile, filename_console, filename, parameters)

            X = []
            Tx = []
            print("Reuse precomputed wall shear stress datafile", file=sys.stdout, flush=True)
            with open(filename, mode='r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    if row:  # Check if row is not empty
                        X.append(row[0])  # First column
                        Tx.append(row[1])  # Second column
            cp_value = extract_reattachment_point_from_dataseries(X,Tx)

            print("Reattachment point: " + str(cp_value), file=sys.stdout, flush=True)

        except Exception as e:
            # Code to handle any exception
            print(f"An error occurred: {e}", file=sys.stdout, flush=True)
            raise Exception("A generic error occurred.")
        
        return [[cp_value]]

    def supports_evaluate(self):
        return True

class CfModel(umbridge.Model):

    def __init__(self):
        super().__init__("forwardcf")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1000,1000]

    def __call__(self, parameters, config):
        try:
            filename, filename_console = configure_case(config)
            print("Case configured", file=sys.stdout, flush=True)

            if not os.path.exists(output_dir+'/wallshear'+filename):
                run_case(tempcasefile, filename_console, filename, parameters)

            X = []
            Tx = []
            print("Reuse precomputed wall shear stress datafile", file=sys.stdout, flush=True)
            with open(output_dir+'/wallshear'+filename, mode='r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    if row:  # Check if row is not empty
                        X.append(row[0])  # First column
                        Tx.append(row[1])  # Second column
            cf = extract_cf_from_dataseries(X,Tx,rhoinf,uinf)

            X_return = np.zeros(1000)
            X_return[0:len(X)] = X
            cf_return = np.zeros(1000)
            cf_return[0:len(X)] = cf
            X_out = X_return.tolist()
            cf_out = cf_return.tolist()

        except Exception as e:
            # Code to handle any exception
            print(f"An error occurred: {e}")
            raise Exception("A generic error occurred.", file=sys.stdout, flush=True)
               
        return [X_out, cf_out]
    
    def supports_evaluate(self):
        return True

class CpModel(umbridge.Model):

    def __init__(self):
        super().__init__("forwardcp")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1000,1000]

    def __call__(self, parameters, config):
        try:
            filename, filename_console = configure_case(config)
            print("Case configured", file=sys.stdout, flush=True)

            if not os.path.exists(output_dir+'/nearWallField'+filename):
                run_case(tempcasefile, filename_console, filename, parameters)

            X = []
            Tx = []

            print("Reuse precomputed wall pressure datafile", file=sys.stdout, flush=True)
            with open(output_dir+'/nearWallField'+filename, mode='r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    if row:  # Check if row is not empty
                        X.append(row[0])  # First column
                        Tx.append(row[1])  # Second column
            
            print("Extract Cp", file=sys.stdout, flush=True)
            uinf = parameters[0][1]
            rhoinf = 1.185
            cp = extract_cp_from_dataseries(X,Tx, rhoinf, uinf)


            X_return = np.zeros(1000)
            X_return[0:len(X)] = X
            cp_return = np.zeros(1000)
            cp_return[0:len(X)] = cp
            x_out = X_return.tolist()
            cp_out = cp_return.tolist()
        except Exception as e:
            # Code to handle any exception
            print(f"An error occurred: {e}", file=sys.stdout, flush=True)
            raise Exception("A generic error occurred.")

        return [x_out, cp_out]
    
    def supports_evaluate(self):
        return True
    
class PwallModel(umbridge.Model):

    def __init__(self):
        super().__init__("forwardpwall")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1000,1000]

    def __call__(self, parameters, config):
        try:
            filename, filename_console = configure_case(config)
            print("Case configured", file=sys.stdout, flush=True)

            if not os.path.exists(output_dir+'/wallPressure'+filename):
                run_case(tempcasefile, filename_console, filename, parameters)

            X = []
            Tx = []

            print("Reuse precomputed wall pressure datafile", file=sys.stdout, flush=True)
            with open(output_dir+'/wallPressure'+filename, mode='r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    if row:  # Check if row is not empty
                        X.append(row[0])  # First column
                        Tx.append(row[1])  # Second column
            
            print("Extract wall pressure", file=sys.stdout, flush=True)

            X_return = np.zeros(1000)
            X_return[0:len(X)] = X
            pwall_return = np.zeros(1000)
            pwall_return[0:len(X)] = Tx
            x_out = X_return.tolist()
            pwall_out = pwall_return.tolist()
        except Exception as e:
            # Code to handle any exception
            print(f"An error occurred: {e}", file=sys.stdout, flush=True)
            raise Exception("A generic error occurred.")

        return [x_out, pwall_out]
    
    def supports_evaluate(self):
        return True
    
class yPlusModel(umbridge.Model):

    def __init__(self):
        super().__init__("forwardyplus")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1000,1000]

    def __call__(self, parameters, config):
        try:
            filename, filename_console = configure_case(config)
            print("Case configured", file=sys.stdout, flush=True)

            if not os.path.exists(output_dir+'/yPlus'+filename):
                run_case(tempcasefile, filename_console, filename, parameters)

            X = []
            Tx = []

            print("Reuse precomputed yPlus datafile", file=sys.stdout, flush=True)
            with open(output_dir+'/yPlus'+filename, mode='r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    if row:  # Check if row is not empty
                        X.append(row[0])  # First column
                        Tx.append(row[1])  # Second column
            
            print("Extract yPlus", file=sys.stdout, flush=True)

            X_return = np.zeros(1000)
            X_return[0:len(X)] = X
            yplus_return = np.zeros(1000)
            yplus_return[0:len(X)] = Tx
            x_out = X_return.tolist()
            yplus_out = yplus_return.tolist()
        except Exception as e:
            # Code to handle any exception
            print(f"An error occurred: {e}", file=sys.stdout, flush=True)
            raise Exception("A generic error occurred.")

        return [x_out, yplus_out]
    
    def supports_evaluate(self):
        return True

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
          input_file}'. Output saved to '{output_file}'.", file=sys.stdout, flush=True)


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
          input_file}'. Output saved to '{output_file}'.", file=sys.stdout, flush=True)


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
          input_file}'. Output saved to '{output_file}'.", file=sys.stdout, flush=True)


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
          input_file}'. Output saved to '{output_file}'.", file=sys.stdout, flush=True)


def replace_inflow_mag(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'INFLOW_MAG' with 'replacement_value'
    modified_data = file_data.replace('INFLOW_MAG', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'INFLOW_MAG' with '{replacement_value}' in '{input_file}'. Output saved to '{output_file}'.", file=sys.stdout, flush=True)

def replace_wallpressuredist(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'WALLPRESSUREDIST' with 'replacement_value'
    modified_data = file_data.replace('WALLPRESSUREDIST', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'WALLPRESSUREDIST' with '{replacement_value}' in '{input_file}'. Output saved to '{output_file}'.", file=sys.stdout, flush=True)


# Define UM-BRIDGE Models and serve
reattachment_model = ReattachmentModel()
cf_model = CfModel()
cp_model = CpModel()
pWall_model = PwallModel()
yplus_model = yPlusModel()
umbridge.serve_models([reattachment_model,cf_model,cp_model,yplus_model,pWall_model], 4242)
