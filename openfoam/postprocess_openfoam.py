import fluidfoam as fluidfoam
from scipy import interpolate
import os
import re
import argparse

def extract_reattachment_point(filename, final_time):
    X, Y, Z = fluidfoam.readmesh(filename,boundary="bottomWall")

    # Find latestTime (may not be final_time if converges early)
    largest_subdir = get_largest_number_subdirectory(filename)

    Tx,Ty,Tz = fluidfoam.readfield(filename,name="wallShearStress",time_name=largest_subdir,boundary="bottomWall")
    # Px,Py,Pz = fluidfoam.readfield(filename,name="p",time_name="5000",boundary="bottomWall")

    spline = interpolate.CubicSpline(X,Tx, extrapolate=False)

    r = spline.roots()

    # Return last root
    return (float(r[-1]) , X, Tx)

def get_largest_number_subdirectory(path):
    # List all subdirectories in the specified path
    subdirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

    # Use regex to extract numbers from subdirectory names
    number_pattern = re.compile(r'(\d+(\.\d+)?)')
    subdir_numbers = []
    
    for subdir in subdirs:
        match = number_pattern.fullmatch(subdir)
        if match:
            subdir_numbers.append(float(match.group(0)))

    if not subdir_numbers:
        return None

    # Find the largest number
    largest_number = max(subdir_numbers)
    largest_subdir = str(largest_number)

    # Ensure the format matches the original subdirectory name
    if largest_subdir not in subdirs:
        largest_subdir = f"{largest_number:.6g}"

    return largest_subdir

# def main():
#     # Set up the argument parser
#     parser = argparse.ArgumentParser(description="A script to run different functions based on command-line arguments")

#     # Add arguments
#     parser.add_argument("--case", type=str, help="This is the case directory")

#     # Parse the arguments
#     args = parser.parse_args()

#     if args.case == None:
#         print("Please specify a case")
#     else:
#         extract_reattachment_point(args.case, final_time)

# if __name__ == "__main__":
#     main()