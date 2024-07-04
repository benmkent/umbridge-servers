import fluidfoam as fluidfoam
from scipy import interpolate

def extract_reattachment_point(filename):
    X, Y, Z = fluidfoam.readmesh(filename,boundary="bottomWall")
    Tx,Ty,Tz = fluidfoam.readfield(filename,name="wallShearStress",time_name="5000",boundary="bottomWall")
    # Px,Py,Pz = fluidfoam.readfield(filename,name="p",time_name="5000",boundary="bottomWall")

    spline = interpolate.CubicSpline(X,Tx)

    r = spline.roots()

    # Return last root
    return r[-1]