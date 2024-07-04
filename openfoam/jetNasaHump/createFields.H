// Reading dictionary - - - - - - - - - - - - - - - - - - - - //
IOdictionary controlDict
(
    IOobject
    (
        "controlDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

word jetName = controlDict.getOrDefault<word>("jet_patch_name", "jet");
scalar inletMag = readScalar(controlDict.lookup("inlet_mag")); 
scalar jetMag = readScalar(controlDict.lookup("jet_mag")); 

// Reading the BC initialisations - - - - - - - - - - - - - - //
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField nuTilda
(
    IOobject
    (
        "nuTilda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

scalarField Cx=mesh.C().component(0);
scalarField Cy=mesh.C().component(1);
