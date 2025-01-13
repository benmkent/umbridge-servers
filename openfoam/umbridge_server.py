import umbridge
import os
import sys
import csv
import numpy as np
import time
from postprocess_openfoam import extract_reattachment_point, extract_reattachment_point_from_dataseries
from postprocess_openfoam import extract_cf, extract_cf_from_dataseries
from postprocess_openfoam import extract_cp, extract_cp_from_dataseries
from postprocess_openfoam import extract_yplus, extract_pWall
from postprocess_openfoam import get_largest_number_subdirectory
from postprocess_openfoam import extract_integrals
from postprocess_openfoam import extract_forces
from postprocess_openfoam import extract_linesamples

def configure_case(config,parameters):
    print("==========START CONFIGURE CASE==========")

    # Decide on fidelity
    if config['Fidelity'] == 10:
        casefile = "./NASA_hump_data_fine"
        print("Selecting fidelity <<fine>>", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 11:
        casefile = "./NASA_hump_data_coarse"
        print("Selecting fidelity <<coarse>>", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 12:
        casefile = "./NASA_hump_data_coarse_2"
        print("Selecting fidelity <<coarse 2>>", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 20:
        casefile = "./NASA_hump_data_jan25_fine"
        print("Selecting fidelity <<jan25_fine>>", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 21:
        casefile = "./NASA_hump_data_jan25_coarse_1"
        print("Selecting fidelity <<jan25_coarse_1>>", file=sys.stdout, flush=True)
    elif config['Fidelity'] == 22:
        casefile = "./NASA_hump_data_jan25_coarse_2"
        print("Selecting fidelity <<jan25_coarse_2>>", file=sys.stdout, flush=True)
    else:
        AssertionError("Unknown config")

    # Copy folder to use as realisation
    tempcasefile = "./caserealisation"
    os.system('rm -rf ' + tempcasefile)
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

    # Ste up yslice
    if 'yslice_x' not in config:
        config['yslice_x'] = 0.5
    replacement_value = str(config['yslice_x'])
    input_file = tempcasefile+"/system/controlDict"
    output_file = input_file
    replace_xforyslice(input_file, output_file, replacement_value)

    # Set up residual tolerances (could be simplified)
    if 'res_tol' not in config:
        config['res_tol'] = 1e-10
        res_tol_str = ''
    elif abs(config['res_tol'] - 1e-10) < 1e-14:
        config['res_tol'] = 1e-10
        res_tol_str = ''
    else:
        res_tol_str = '_restol_' + str(config['res_tol'])
    if 'res_tol_u' not in config:
        config['res_tol_u'] = config['res_tol']
        res_tol_str = res_tol_str+''
    elif abs(config['res_tol_u'] - 1e-10) < 1e-14:
        config['res_tol_u'] = 1e-10
        res_tol_str = res_tol_str+''
    else:
        res_tol_str = res_tol_str+'_resU_' + str(config['res_tol_u'])

    if 'res_tol_p' not in config:
        config['res_tol_p'] = config['res_tol']
        res_tol_str = res_tol_str+''
    elif abs(config['res_tol_p'] - 1e-10) < 1e-14:
        config['res_tol_p'] = 1e-10
        res_tol_str = res_tol_str+''
    else:
        res_tol_str = res_tol_str+'_resP_' + str(config['res_tol_p'])

    print("Residual tolerances "+res_tol_str, file=sys.stdout, flush=True)

    if 'abs_tol' not in config:
        config['abs_tol'] = 1e-10
        abs_tol_str = ''
    elif abs(config['abs_tol'] - 1e-10) < 1e-14:
        abs_tol_str = ''
    else:
        abs_tol_str = '_abstol_' + str(config['abs_tol'])
    print("Iterative solver tolerance "+str(config['abs_tol']), file=sys.stdout, flush=True)

    input_file = tempcasefile+"/system/fvSolution"
    output_file = input_file
    replacement_value = str(config['res_tol_u'])
    replace_res_tol_u(input_file, output_file, replacement_value)
    replacement_value = str(config['res_tol_p'])
    replace_res_tol_p(input_file, output_file, replacement_value)
    replacement_value = str(config['res_tol'])
    replace_res_tol(input_file, output_file, replacement_value)
    replacement_value = str(config['abs_tol'])
    replace_abs_tol(input_file, output_file, replacement_value)

    replacement_value_jet = f"{float(replacement_value_jet):.12g}"
    replacement_value_inflow = f"{float(replacement_value_inflow):.12g}"

    if 'qoi' in config:
        # nasa2dwmh naming convention
        filename = 'Jet_'+str(replacement_value_jet) + '_Inflow_' + str(replacement_value_inflow) + '_Fidelity_' + str(config['Fidelity']) +'_' +res_tol_str + abs_tol_str +'.csv'
        filename_console = 'console_Jet_'+str(replacement_value_jet) + '_Inflow_' + str(replacement_value_inflow) + '_Fidelity_' + str(config['Fidelity']) +res_tol_str + abs_tol_str +'.log'
    else:
        # legacy naming convention
        filename = str(replacement_value_jet) + 'Inflow' + str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + res_tol_str + abs_tol_str +'.csv'
        filename_console = 'console_'+ str(replacement_value_jet) + 'Inflow' + str(replacement_value_inflow) + 'Fidelity' + str(config['Fidelity']) + res_tol_str + abs_tol_str +'.log'

    # BlockMesh for info
    # os.system('openfoam2406 blockMesh -case '+tempcasefile)
    os.system('openfoam2406 topoSet -case '+tempcasefile)


    print('======CASE CONFIGURED======')
    return filename, filename_console

def copy_case(foldername):
    tempcasefile = "./caserealisation"
    os.system('cp -r outputdata/'+ foldername +'/* '+tempcasefile)


def run_case(filename_console, filename, parameters,copywholecase=False):
    output_dir = './outputdata'
    tempcasefile = "./caserealisation"

    # BlockMesh for info
    os.system('cat '+ tempcasefile+'/system/blockMeshDict' + ' | tee -a ' + output_dir+'/'+filename_console)
    os.system('openfoam2406 blockMesh -case '+tempcasefile + ' | tee -a ' + output_dir+'/'+filename_console)
    os.system('openfoam2406 topoSet -case '+tempcasefile + ' | tee -a ' + output_dir+'/'+filename_console)

    # Set up boundary conditions
    print("Enforcing boundary conditions jetNasaHump", file=sys.stdout, flush=True)
    os.system('openfoam2406 jetNasaHump -case '+tempcasefile + ' | tee -a ' + output_dir+'/'+filename_console)

    # Run simple foam
    print("Run simplefoam", file=sys.stdout, flush=True)
    os.system('openfoam2406 simpleFoam -case '+tempcasefile + ' | tee -a ' + output_dir+'/'+filename_console)

    # Copy out console file (already prints directly to PV)
    # os.system('cp '+tempcasefile+'/'+filename_console+' '+output_dir+'/'+filename_console)

    # Save out wallshear data
    print("Extract reattachment point", file=sys.stdout, flush=True)
    (rp_value, X, Tx) = extract_reattachment_point(tempcasefile, 5000)
    wall_shear = np.column_stack((X, Tx))
    np.savetxt(output_dir+'/wallshearJet'+filename, wall_shear, fmt='%f')
    print("Reattachment point "+str(rp_value))

    # Save out nearWallField data
    uinf = parameters[0][1]
    rhoinf = 1.185
    print("Extract nearWallFields Cp", file=sys.stdout, flush=True)
    (cp_value_nwf, X, Tx) = extract_cp(tempcasefile, 5000, rhoinf, uinf)
    cp = np.column_stack((X, Tx))
    np.savetxt(output_dir+'/nearWallField'+filename, cp, fmt='%f')
    print("Cp max min "+str(max(Tx))+" "+str(min(Tx)))

    print("Extract pressure field at wall", file=sys.stdout, flush=True)
    (X, Tx) = extract_pWall(tempcasefile, 5000)
    pwall= np.column_stack((X, Tx))
    np.savetxt(output_dir+'/wallPressure'+filename, pwall, fmt='%f')
    print("p wall max min "+str(max(Tx))+" "+str(min(Tx)))

    print("Extract yPlus at wall", file=sys.stdout, flush=True)
    (X,Tx) = extract_yplus(tempcasefile, 5000)
    yp= np.column_stack((X, Tx))
    np.savetxt(output_dir+'/yPlus'+filename, yp, fmt='%f')
    print("yPlus max min "+str(max(Tx))+" "+str(min(Tx)))

    max_iteration_number = get_largest_number_subdirectory(tempcasefile)
    foldername = filename[0:-4]
    os.system('mkdir -p outputdata/'+foldername)
    if copywholecase == False:
        os.system('cp -r ' + tempcasefile +'/' + max_iteration_number +' outputdata/'+foldername+'/' + max_iteration_number)
        os.system('cp -r ' + tempcasefile +'/postProcessing' +' outputdata/'+foldername+'/postProcessing')
    else:
        os.system('cp -r ' + tempcasefile +'/* outputdata/'+foldername)

class Nasa2DWMHModel(umbridge.Model):
    def __init__(self):
        super().__init__("forward2dwmh")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        if config['qoi'] == 'reattachmentpoint' or config['qoi'] == 'exectime':
            return [1]
        elif config['qoi'] == 'cf' or config['qoi'] == 'cp' or config['qoi'] == 'yplus' or config['qoi'] == 'p':
            return [1000,1000]
        elif config['qoi'] == 'pressureintegrals':
            return [2]
        elif config['qoi'] == 'forces':
            return [3,3,3]
        elif config['qoi'] == 'yslice':
            return [200,200,200,200,200]
        else:
            raise Exception('unknown qoi')
        
    def __call__(self, parameters, config):
        # First configure case based on input parameters
        filename, filename_console = configure_case(config,parameters)
        print("Case configured: "+filename, file=sys.stdout, flush=True)
        foldername = filename[0:-4]

        # Now check if we have saved case data
        if os.path.exists('outputdata/'+foldername) and os.path.isdir('outputdata/'+foldername):
            # Copy the final saved case data
            copy_case(foldername)
        else:
            # We must run the case and save out the data
            t = time.time()
            if 'saveall' not in config:
                run_case(filename_console, filename, parameters)
            elif config['saveall'] == 1:
                run_case(filename_console, filename, parameters, copywholecase=True)
            else:
                run_case(filename_console, filename, parameters)

            elapsed = time.time() - t
            # Write the number to the file
            with open('outputdata/exectime_'+filename, mode="w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow([elapsed])  # Write the number as a single row

            if 'debug' not in config:
                # Run case usually deletes the temporary case so we recreate it and copy the data
                filename, filename_console = configure_case(config,parameters)
                copy_case(foldername)

        # Use post process to evaluate the OpenFOAM functions
        tempcasefile = './caserealisation'
        os.system('openfoam2406 simpleFoam -postProcess -case '+tempcasefile + ' | tee -a outputdata/'+filename_console)

        # Now extract the proper QoI
        uinf=parameters[0][1]
        rhoinf=1.0
        if config['qoi'] == 'reattachmentpoint':
            print("Extract reattachment point", file=sys.stdout, flush=True)
            (rp_value, X, Tx) = extract_reattachment_point(tempcasefile, 5000)
            print("Reattachment point "+str(rp_value))
            return [[rp_value]]
        elif config['qoi'] == 'cf':
            (cf, X, Tx) = extract_cf(tempcasefile, 999999, rhoinf, uinf)
            X_return = np.zeros(1000)
            X_return[0:len(X)] = X
            cf_return = np.zeros(1000)
            cf_return[0:len(X)] = cf
            X_out = X_return.tolist()
            cf_out = cf_return.tolist()
            return [X_out, cf_out]
        elif config['qoi'] == 'cp':
            (cp, X, Tx) = extract_cp(tempcasefile, 999999, rhoinf, uinf)
            X_return = np.zeros(1000)
            X_return[0:len(X)] = X
            cp_return = np.zeros(1000)
            cp_return[0:len(X)] = cp
            X_out = X_return.tolist()
            cp_out = cp_return.tolist()
            return [X_out, cp_out]            
        elif config['qoi'] == 'yplus':
            (X,Tx) = extract_yplus(tempcasefile, 5000)
            X_return = np.zeros(1000)
            X_return[0:len(X)] = X
            yplus_return = np.zeros(1000)
            yplus_return[0:len(X)] = Tx
            X_out = X_return.tolist()
            yplus_out = yplus_return.tolist()
            return [X_out,yplus_out]
        elif config['qoi'] == 'p':
            (X, Tx) = extract_pWall(tempcasefile, 5000)
            X_return = np.zeros(1000)
            X_return[0:len(X)] = X
            pwall_return = np.zeros(1000)
            pwall_return[0:len(X)] = Tx
            X_out = X_return.tolist()
            pwall_out = pwall_return.tolist()
            return [X_out,pwall_out]
        elif config['qoi'] == 'pressureintegrals':
            integrals = extract_integrals(tempcasefile)
            return [integrals]
        elif config['qoi'] == 'forces':
            [total_forces,pressure_forces,viscous_forces] = extract_forces(tempcasefile)
            return [total_forces,pressure_forces,viscous_forces]
        elif config['qoi'] == 'exectime':
            # Read the number from the file
            with open('outputdata/exectime_'+filename, mode="r") as file:
                reader = csv.reader(file)
                for row in reader:
                    if row:  # Skip empty rows
                        t = float(row[0])  # Convert the first value to a float
            return [[t]]
        elif config['qoi'] == 'yslice':
            (y, p, ux, uy, uz) = extract_linesamples(tempcasefile, 5000)
            y_return = np.zeros(200)
            y_return[0:len(y)] = y
            p_return = np.zeros(200)
            p_return[0:len(p)] = p
            ux_return = np.zeros(200)
            ux_return[0:len(ux)] = ux
            uy_return = np.zeros(200)
            uy_return[0:len(uy)] = uy
            uz_return = np.zeros(200)
            uz_return[0:len(uz)] = uz
            return [y_return.tolist(),p_return.tolist(),ux_return.tolist(),uy_return.tolist(),uz_return.tolist()]
        else:
            raise Exception('unknown qoi')
        
    def supports_evaluate(self):
        return True
    

class ReattachmentModel(umbridge.Model):

    def __init__(self):
        super().__init__("forwardopenfoam")

    def get_input_sizes(self, config):
        return [2]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        try:
            filename, filename_console = configure_case(config,parameters)
            print("Case configured: "+filename, file=sys.stdout, flush=True)
            
            output_dir = './outputdata'

            if not os.path.exists(output_dir+'/wallshearJet'+filename):
                run_case(filename_console, filename, parameters)

            X = []
            Tx = []
            print("Reuse precomputed wall shear stress datafile", file=sys.stdout, flush=True)
            with open(output_dir+'/wallshearJet'+filename, mode='r') as file:
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
            filename, filename_console = configure_case(config,parameters)
            print("Case configured: "+filename, file=sys.stdout, flush=True)

            output_dir = './outputdata'

            if not os.path.exists(output_dir+'/wallshearJet'+filename):
                run_case(filename_console, filename, parameters)

            X = []
            Tx = []
            print("Reuse precomputed wall shear stress datafile", file=sys.stdout, flush=True)
            with open(output_dir+'/wallshearJet'+filename, mode='r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    if row:  # Check if row is not empty
                        X.append(row[0])  # First column
                        Tx.append(row[1])  # Second column
            
            rhoinf = 1.0
            uinf = parameters[0][1]
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
            filename, filename_console = configure_case(config,parameters)
            print("Case configured: "+filename, file=sys.stdout, flush=True)

            output_dir = './outputdata'

            if not os.path.exists(output_dir+'/nearWallField'+filename):
                run_case(filename_console, filename, parameters)

            X = []
            Tx = []

            print("Reuse precomputed nearWallField datafile", file=sys.stdout, flush=True)
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
            filename, filename_console = configure_case(config,parameters)
            print("Case configured: "+filename, file=sys.stdout, flush=True)

            output_dir = './outputdata'

            if not os.path.exists(output_dir+'/wallPressure'+filename):
                run_case(filename_console, filename, parameters)

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
            filename, filename_console = configure_case(config,parameters)
            print("Case configured: "+filename, file=sys.stdout, flush=True)
            
            output_dir = './outputdata'

            if not os.path.exists(output_dir+'/yPlus'+filename):
                run_case(filename_console, filename, parameters)

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

def replace_res_tol_u(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'RES_TOL_U' with 'replacement_value'
    modified_data = file_data.replace('RES_TOL_U', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'RES_TOL_U' with '{replacement_value}' in '{
          input_file}'. Output saved to '{output_file}'.", file=sys.stdout, flush=True)

def replace_res_tol_p(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'RES_TOL_P' with 'replacement_value'
    modified_data = file_data.replace('RES_TOL_P', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'RES_TOL_P' with '{replacement_value}' in '{
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

def replace_xforyslice(input_file, output_file, replacement_value):
    # Read input file
    with open(input_file, 'r') as f:
        file_data = f.read()

    # Replace all occurrences of 'XFORYSLICE' with 'replacement_value'
    modified_data = file_data.replace('XFORYSLICE', replacement_value)

    # Write modified content to output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    print(f"Replaced 'XFORYSLICE' with '{replacement_value}' in '{input_file}'. Output saved to '{output_file}'.", file=sys.stdout, flush=True)

# Define UM-BRIDGE Models and serve
reattachment_model = ReattachmentModel()
cf_model = CfModel()
cp_model = CpModel()
pWall_model = PwallModel()
yplus_model = yPlusModel()
nasa2dwmh = Nasa2DWMHModel()

umbridge.serve_models([reattachment_model,cf_model,cp_model,yplus_model,pWall_model,nasa2dwmh], 4242)
