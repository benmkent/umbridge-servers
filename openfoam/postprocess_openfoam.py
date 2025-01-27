import fluidfoam as fluidfoam
from scipy import interpolate
from numpy import linspace
import numpy as np
import os
import re
import argparse

def extract_reattachment_point(filename, final_time=None):
    """
    Extracts the reattachment point from wall shear stress data at the bottom wall.

    Parameters:
        filename (str): The case directory containing the mesh and field data.
        final_time (float, optional): The final simulation time (not used).

    Returns:
        tuple: Contains the reattachment point (x), the X-coordinates of the bottom wall, 
               and the wall shear stress (Tx).
    """
    # Read the bottom wall mesh
    X, Y, Z = fluidfoam.readmesh(filename, boundary="bottomWall")
    try:
        X_hump, Y_hump, Z_hump = fluidfoam.readmesh(filename, boundary="hump")
    except Exception as e:
        print("Hump not found, setting it to 0")
        X_hump = []
        Y_hump = []
        Z_hump = []
    
    # Identify the latest available time directory
    largest_subdir = get_largest_number_subdirectory(filename)

    # Read the wall shear stress field data
    Tx, Ty, Tz = fluidfoam.readfield(
        filename,
        name="wallShearStress",
        time_name=largest_subdir,
        boundary="bottomWall"
    )

    try:
        Tx_hump, Ty_hump, Tz_hump = fluidfoam.readfield(
            filename,
            name="wallShearStress",
            time_name=largest_subdir,
            boundary="hump"
        )
    except Exception as e:
        print("Hump not found, setting it to 0")
        Tx_hump = []
        Ty_hump = []
        Tz_hump = []

    if len(X_hump) > 0:
        X = X + X_hump
        Tx = Tx + X_hump

        X, sort_pattern = zip(*sorted((x, i) for i, x in enumerate(X)))
        Tx = [Tx[i] for i in sort_pattern]

    # Extract the reattachment point based on the data series
    try:
        x = extract_reattachment_point_from_dataseries(X, Tx)
    except Exception as e:
        raise ValueError("Error while extracting the reattachment point.") from e

    return x, X, Tx


def extract_reattachment_point_from_dataseries(X, Tx):
    """
    Extracts the reattachment point from the wall shear stress.

    Parameters:
        X : X-coordinates.
        Tx : wall shear stress values.

    Returns:
        float: Reattachment point .
    """
    # Create a cubic spline interpolation for the data
    spline = interpolate.CubicSpline(X, Tx, extrapolate=False)

    # Find the roots
    roots = spline.roots()

    # Return the last root or issue a warning if none found
    if len(roots) > 0:
        reattachment_point = float(roots[-1])
    else:
        print("Unable to detect reattachment point. Setting it to 0.66")
        reattachment_point = 0.66

    return reattachment_point

def extract_cf(filename, final_time, rhoinf, uinf):
    """
    Extracts the skin-friction coefficient (Cf) from the wall shear stress data.

    Parameters:
        filename (str): Case directory
        final_time (float):Final simulation time (not used).
        rhoinf (float): Density
        uinf (float): Inflow velocity

    Returns:
        tuple Skin-friction coefficient (cf), the X-coordinates of the bottom wall,
               and the wall shear stress (Tx).
    """
    # Read mesh
    X, Y, Z = fluidfoam.readmesh(filename, boundary="bottomWall")

    # Latest available time directory
    largest_subdir = get_largest_number_subdirectory(filename)

    # Read the wall shear stress field data
    Tx, Ty, Tz = fluidfoam.readfield(
        filename,
        name="wallShearStress",
        time_name=largest_subdir,
        boundary="bottomWall"
    )

    # Extract the skin-friction coefficient
    try:
        cf = extract_cf_from_dataseries(X, Tx, rhoinf, uinf)
    except Exception as e:
        raise ValueError("Failed to extract the skin-friction coefficient.") from e

    return cf, X, Tx

import numpy as np

def extract_cf_from_dataseries(X, Tx, rhoinf, uinf):
    """
    Computes the skin-friction coefficient (Cf) from wall shear stress data.

    Parameters:
        X: X-coordinates.
        Tx: wall shear stress values.
        rhoinf : freestream density.
        uinf : freestream velocity.

    Returns:
       Skin-friction coefficient.
    """
    Tx = np.array(Tx, dtype=float)

    # Compute the skin-friction coefficient
    cf_at_x = Tx / (0.5 * rhoinf * uinf**2)

    return cf_at_x

def extract_cp(filename, final_time, rhoinf, uinf):
    """
    Extracts the pressure coefficient (Cp) from wall pressure data.

    Parameters:
        filename: Case directory
        final_time: Final simulation time (not used).
        rhoinf: Freestream density.
        uinf: Freestream velocity.

    Returns:
        Pressure coefficient (cp), X-coordinates of the bottom wall,
               and the wall pressure (P).
    """
    # Read the bottom wall mesh
    X, Y, Z = fluidfoam.readmesh(filename, boundary="bottomWall")

    # Identify the latest available time directory
    largest_subdir = get_largest_number_subdirectory(filename)

    # Read the wall pressure field data
    P = fluidfoam.readfield(
        filename,
        name="pWall",
        time_name=largest_subdir,
        boundary="bottomWall"
    )

    # Compute the pressure coefficient
    try:
        cp = extract_cp_from_dataseries(X, P, rhoinf, uinf)
    except Exception as e:
        raise ValueError("Failed to extract the pressure coefficient.") from e

    return cp, X, P


def extract_cp_from_dataseries(X, P, rhoinf, uinf):
    """
    Computes the pressure coefficient (Cp) from pressure data.

    Parameters:
        X : X-coordinates.
        P : Pressure
        rhoinf :Freestream density.
        uinf : Freestream velocity.

    Returns:
        Pressure coefficient (Cp).
    """
    # Convert P to a NumPy array for consistency
    P = np.array(P, dtype=float)

    # Compute the pressure coefficient
    cp_at_x = P / (0.5 * rhoinf * uinf**2)

    return cp_at_x

def extract_yplus(filename, final_time=None):
    """
    Extracts the y+ values.

    Parameters:
        filename : Case directory.
        final_time : Final simulation time (not used).

    Returns:
        X-coordinates of the bottom wall and the y+ values.
    """

    # Identify the latest available time directory
    largest_subdir = get_largest_number_subdirectory(filename)

    # Read the bottom wall mesh
    X1, Y, Z = fluidfoam.readmesh(filename, boundary="bottomWall")

    # Read the y+ field data
    yPlus1 = fluidfoam.readfield(
        filename,
        name="yPlus",
        time_name=largest_subdir,
        boundary="bottomWall"
    )

    try:
        # Read the hump mesh
        X2, Y, Z = fluidfoam.readmesh(filename, boundary="hump")
        # Read the y+ field data
        yPlus2 = fluidfoam.readfield(
            filename,
            name="yPlus",
            time_name=largest_subdir,
            boundary="hump"
        )
    except:
        # If the hump is not defined
        X2 = []
        yPlus2 = []

    dataset1 = list(zip(X1, yPlus1))
    dataset2 = list(zip(X2, yPlus2))
    
    # Interleave and combine the datasets
    combined_dataset = dataset1 + dataset2
    
    # Sort the combined dataset (first by x, then by y)
    sorted_dataset = sorted(combined_dataset)
    X, yPlus = zip(*sorted_dataset)

    return X, yPlus

def extract_pWall(filename, final_time=None):
    """
    Extracts wall pressure (p) values from the bottom wall.

    Parameters:
        filename : Case directory
        final_time : Final simulation time (not used).

    Returns:
       X-coordinates of the bottom wall and the wall pressure (p).
    """
    # Read the bottom wall mesh
    X, Y, Z = fluidfoam.readmesh(filename, boundary="bottomWall")

    # Identify the latest available time directory
    largest_subdir = get_largest_number_subdirectory(filename)

    p = fluidfoam.readfield(
        filename,
        name="p",
        time_name=largest_subdir,
        boundary="bottomWall"
    )

    return X, p


def extract_integrals(filename):
    """
    Extracts the integrals of pressure before and after the jet

    Parameters:
        filename : Case directory

    Returns:
        Before integral and after integral values.
    """
    beforeIntegral = None
    afterIntegral = None
    
    # Open the file and process line by line
    with open(filename+'/postProcessing/integrationBeforeJet/0/surfaceFieldValue.dat', 'r') as file:
        for line in file:
            # Split the line into columns and check if it matches the data format
            parts = line.split()
            if len(parts) == 2:  # Assuming data lines have exactly 2 columns
                try:
                    # Attempt to convert the second column to a float
                    beforeIntegral = float(parts[1])
                except ValueError:
                    pass  # Skip lines that don't match the expected format
    # Open the file and process line by line
    with open(filename+'/postProcessing/integrationAfterJet/0/surfaceFieldValue.dat', 'r') as file:
        for line in file:
            # Split the line into columns and check if it matches the data format
            parts = line.split()
            if len(parts) == 2:  # Assuming data lines have exactly 2 columns
                try:
                    # Attempt to convert the second column to a float
                    afterIntegral = float(parts[1])
                except ValueError:
                    pass  # Skip lines that don't match the expected format                
    return [beforeIntegral,afterIntegral]

def extract_forces(filename):
    """
    Extracts total, pressure, and viscous forces.

    Parameters:
        filename : Case directory

    Returns:
        list: Total forces, pressure forces, and viscous forces. Each is a list of three components.
    """
    
    total_forces = None
    pressure_forces = None
    viscous_forces = None
    
    # Construct the file path
    file_path = f"{filename}/postProcessing/forces1/0/force.dat"
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                components = line.split()
                
                # Ensure the line has the expected number of components
                if len(components) == 10:
                    try:
                        # Extract forces from the components
                        time = float(components[0])
                        total_forces = [float(components[i]) for i in range(1, 4)]
                        pressure_forces = [float(components[i]) for i in range(4, 7)]
                        viscous_forces = [float(components[i]) for i in range(7, 10)]
                    except ValueError:
                        print(f"Warning: Malformed data line (could not convert to float): {line}")
                        continue  # Skip malformed lines

    except FileNotFoundError:
        print(f"Warning: File {file_path} not found.")
        return [None, None, None]

    return [total_forces, pressure_forces, viscous_forces]


def extract_linesamples(filename, final_time=None):
    """
    Line sample data of position (y), pressure (p),
    and velocity components (ux, uy, uz).

    Parameters:
        filename : Case directory
        final_time (optional): Not Used

    Returns:
        tuple: A tuple containing five lists: y, p, ux, uy, uz.
               - y: y coordinate
               - p: pressure
               - ux: x-component of velocity
               - uy: y-component of velocity
               - uz: z-component of velocity
    """
    
    # Find the largest subdirectory (latest time step)
    largest_subdir = get_largest_number_subdirectory(filename)
    
    # Initialize an empty list to store the data
    data = []

    # Construct the file path
    file_path = f'{filename}/postProcessing/lineSample/{largest_subdir}/ySlice_p_U.xy'
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split()

                # Ensure the line has 5 values (y, p, ux, uy, uz)
                if len(columns) == 5:
                    try:
                        # Convert columns to floats and append to data
                        row = [float(value) for value in columns]
                        data.append(row)
                    except ValueError:
                        print(f"Warning: Skipping malformed line: {line.strip()}")
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        return [], [], [], [], []

    # Unpack the data into separate lists for each variable
    if data:
        y, p, ux, uy, uz = zip(*data)
        return y, p, ux, uy, uz
    else:
        print("Warning: No valid data found in the file.")
        return [], [], [], [], []

def extract_residuals(filename):
    # Define the file path
    file_path = f'{filename}/postProcessing/solverInfo/0/solverInfo.dat'

    # Initialize lists to store time and residuals
    time = []
    initial_residuals = {
        "Ux": [],
        "Uy": [],
        "nuTilda": [],
        "p": []
    }
    final_residuals = {
        "Ux": [],
        "Uy": [],
        "nuTilda": [],
        "p": []
    }
    
    print("Open file "+file_path)
    # Read and process the file
    with open(file_path, "r") as file:
        lines = file.readlines()[2:]  # Skip the header lines
        for line in lines:
            if not line.strip():  # Skip empty lines
                continue
            print(line)

            columns = line.split()  # Split line into columns based on whitespace
            
            # Extract time
            time.append(float(columns[0]))
            
            # Extract initial and final residuals for each field
            initial_residuals["Ux"].append(float(columns[2]))
            final_residuals["Ux"].append(float(columns[3]))
            
            initial_residuals["Uy"].append(float(columns[5]))
            final_residuals["Uy"].append(float(columns[6]))

            if len(columns) < 19:
                initial_residuals["nuTilda"].append(0.0)
                final_residuals["nuTilda"].append(0.0)
                
                initial_residuals["p"].append(float(columns[10]))
                final_residuals["p"].append(float(columns[11]))
            else:
                initial_residuals["nuTilda"].append(float(columns[10]))
                final_residuals["nuTilda"].append(float(columns[11]))
                
                initial_residuals["p"].append(float(columns[15]))
                final_residuals["p"].append(float(columns[16]))

    n = len(time)
    time.extend([0.0] * (5000 - n))
    initial_residuals["Ux"].extend([0.0] * (5000 - n))
    initial_residuals["Uy"].extend([0.0] * (5000 - n))
    initial_residuals["p"].extend([0.0] * (5000 - n))
    initial_residuals["nuTilda"].extend([0.0] * (5000 - n))
    final_residuals["Ux"].extend([0.0] * (5000 - n))
    final_residuals["Uy"].extend([0.0] * (5000 - n))
    final_residuals["p"].extend([0.0] * (5000 - n))
    final_residuals["nuTilda"].extend([0.0] * (5000 - n))


    return (time, initial_residuals, final_residuals)

import os
import re

def get_largest_number_subdirectory(path):
    """
    Finds the largest numerical value subdirectory.

    Parameters:
        path : Case directory 

    Returns:
        str: The largest number, or None if no such subdirectory exists.
    """
    
    # Check if the directory exists
    if not os.path.exists(path):
        print(f"Error: The directory {path} does not exist.")
        return None
    
    # List all subdirectories in the specified path
    subdirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

    # If no subdirectories are found, return None
    if not subdirs:
        print(f"Warning: No subdirectories found in {path}.")
        return None
    
    # Regex pattern to match numbers (including decimal points)
    number_pattern = re.compile(r'(\d+(\.\d+)?)')

    largest_number = None
    largest_subdir = None
    
    # Iterate over subdirectories to find the largest numerical value
    for subdir in subdirs:
        match = number_pattern.fullmatch(subdir)
        if match:
            try:
                num = float(match.group(0))
                # Update largest number and corresponding subdirectory
                if largest_number is None or num > largest_number:
                    largest_number = num
                    largest_subdir = subdir
            except ValueError:
                continue  # Skip non-numeric subdirectory names

    if largest_subdir is None:
        print("Warning: No valid numeric subdirectories found.")
        return None

    return largest_subdir