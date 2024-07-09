import fluidfoam as fluidfoam
from scipy import interpolate

def extract_reattachment_point(filename, final_time):
    X, Y, Z = fluidfoam.readmesh(filename,boundary="bottomWall")
    Tx,Ty,Tz = fluidfoam.readfield(filename,name="wallShearStress",time_name=str(final_time),boundary="bottomWall")
    # Px,Py,Pz = fluidfoam.readfield(filename,name="p",time_name="5000",boundary="bottomWall")

    spline = interpolate.CubicSpline(X,Tx, extrapolate=False)

    r = spline.roots()

    # Return last root
    return (float(r[-1]) , X, Tx)