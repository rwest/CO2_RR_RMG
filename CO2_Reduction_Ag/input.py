# Data sources
database(
    thermoLibraries=['CO2RR_Adsorbates_Ag111', 'surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo', 'electrocatThermo'],
    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006_adjusted', False)],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['electrochem',
                        # 'surface',
                        'Surface_Abstraction',
                        'Surface_Abstraction_vdW',
                        'Surface_Abstraction_Single_vdW',
                        'Surface_Abstraction_Beta_double_vdW',
                        'Surface_Adsorption_Dissociative',
                        'Surface_Adsorption_Dissociative_Double',
                        'Surface_Adsorption_vdW',
                        'Surface_Dissociation',
                        'Surface_Dissociation_Double_vdW',
                        'Surface_Dissociation_vdW',
                        'Surface_EleyRideal_Addition_Multiple_Bond',
                        'Surface_Migration',
                        ],
    kineticsEstimator = 'rate rules',

)

catalystProperties(
    metal = 'Ag111'
)

# List of species
species(
    label='CO2',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 O u0 p2 c0 {2,D}
"""),
)


species(
   label='proton',
   reactive=True,
   structure=adjacencyList(
       """
1 H u0 p0 c+1
"""),
)

species(
    label='vacantX',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

species(
   label='H',
   reactive=True,
   structure=adjacencyList(
       """
1 H u1 p0 c0
"""),
)

# species(
#    label='H2O',
#    reactive=False,
#    structure=adjacencyList(
#        """
# 1 H u0 p0 c0 {3,S}
# 2 H u0 p0 c0 {3,S}
# 3 O u0 p2 c0 {1,S} {2,S}
# """),
# )

species(
    label='CO2X',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
4 X u0 p0 c0
"""),
)

species(
    label='CHO2X',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 X u0 p0 c0 {1,S}
"""),
)

species(
    label='CO2HX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,S} {4,S}
2 C u0 p0 c0 {1,S} {3,D} {5,S}
3 O u0 p2 c0 {2,D}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {2,S}

"""),
)

species(
    label='OCX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 X u0 p0 c0 {2,D}
"""),
)

species(
    label='OX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
"""),
)

species(
    label='CH2O2X',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
6 X u0 p0 c0
"""),
)

species(
    label='CHOX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 X u0 p0 c0 {2,S}
"""),
)

species(
    label='CH2OX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 X u0 p0 c0
"""),
)

species(
    label='CO',
    reactive=True,
    structure=SMILES('[C-]#[O+]')
)


forbidden(
        label='CO2-bidentate',
        structure=adjacencyList(
                """
                1 O u0 p2 c0 {2,D}
                2 C u0 p0 c0 {1,D} {3,S} {4,S}
                3 X u0 p0 c0 {2,S}
                4 O u0 p2 c0 {2,S} {5,S}
                5 X u0 p0 c0 {4,S}
                """
        )
)

liquidSurfaceReactor(
    temperature=(300,'K'),
    liqPotential=(0,'V'),
    surfPotential=(-1.914,'V'),
    initialConcentrations={
        # "H2O": (55.6e3, 'mol/m^3'),
        "CO2": (1.0e1,'mol/m^2'),
        "proton": (1.0e-4,'mol/m^3'),
    },
	initialSurfaceCoverages={
        # "HX": 0.5,
        # # "CXO2": 0.0,
        "CHO2X": 0.1,
        "CO2HX": 0.1,
        "vacantX": 0.1,
        "CO2X": 0.4,
        'OX': 0.1,
        'OCX': 0.1,
        'CH2O2X': 0.05,
        'CHOX': 0.04,
        'CH2OX': 0.01
    },
    surfaceVolumeRatio=(36, 'm^-1'),
    terminationTime=(1.0e3,'sec'),
    # terminationConversion={'CO2': 0.90},
    constantSpecies=["proton","CO2"],
 )

liquidSurfaceReactor(
    temperature=(300,'K'),
    liqPotential=(0,'V'),
    surfPotential=(-1.414,'V'),
    initialConcentrations={
        # "H2O": (55.6e3, 'mol/m^3'),
        "CO2": (1.0e1,'mol/m^2'),
        "proton": (1.0e-4,'mol/m^3'),
    },
	initialSurfaceCoverages={
        # "HX": 0.5,
        # # "CXO2": 0.0,
        "CHO2X": 0.1,
        "CO2HX": 0.1,
        "vacantX": 0.1,
        "CO2X": 0.4,
        'OX': 0.1,
        'OCX': 0.1,
        'CH2O2X': 0.05,
        'CHOX': 0.04,
        'CH2OX': 0.01
    },
    surfaceVolumeRatio=(36, 'm^-1'),
    terminationTime=(1.0e3,'sec'),
    # terminationConversion={'CO2': 0.90},
    constantSpecies=["proton","CO2"],
 )

liquidSurfaceReactor(
    temperature=(300,'K'),
    liqPotential=(0,'V'),
    surfPotential=(-0.914,'V'),
    initialConcentrations={
        # "H2O": (55.6e3, 'mol/m^3'),
        "CO2": (1.0e1,'mol/m^2'),
        "proton": (1.0e-4,'mol/m^3'),
    },
	initialSurfaceCoverages={
        # "HX": 0.5,
        # # "CXO2": 0.0,
        "CHO2X": 0.1,
        "CO2HX": 0.1,
        "vacantX": 0.1,
        "CO2X": 0.4,
        'OX': 0.1,
        'OCX': 0.1,
        'CH2O2X': 0.05,
        'CHOX': 0.04,
        'CH2OX': 0.01
    },
    surfaceVolumeRatio=(36, 'm^-1'),
    terminationTime=(1.0e3,'sec'),
    # terminationConversion={'CO2': 0.90},
    constantSpecies=["proton","CO2"],
 )

liquidSurfaceReactor(
    temperature=(300,'K'),
    liqPotential=(0,'V'),
    surfPotential=(-0.614,'V'),
    initialConcentrations={
        # "H2O": (55.6e3, 'mol/m^3'),
        "CO2": (1.0e1,'mol/m^2'),
        "proton": (1.0e-4,'mol/m^3'),
    },
	initialSurfaceCoverages={
        # "HX": 0.5,
        # # "CXO2": 0.0,
        "CHO2X": 0.1,
        "CO2HX": 0.1,
        "vacantX": 0.1,
        "CO2X": 0.4,
        'OX': 0.1,
        'OCX': 0.1,
        'CH2O2X': 0.05,
        'CHOX': 0.04,
        'CH2OX': 0.01
    },
    surfaceVolumeRatio=(36, 'm^-1'),
    terminationTime=(1.0e3,'sec'),
    # terminationConversion={'CO2': 0.90},
    constantSpecies=["proton","CO2"],
 )

solvation(
	solvent='dummy solvent'
)

simulator(
    atol=1e-22,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=1E-16,
    toleranceMoveToCore=1E-1,
    toleranceRadMoveToCore=1E-6,
    toleranceInterruptSimulation=1E3,
    filterReactions=False,
    maximumEdgeSpecies=5000,
    # toleranceBranchReactionToCore=1E-6,
    # branchingIndex=0.5,
    # branchingRatioMax=1.0,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=True,
    saveEdgeSpecies=True,
    saveSimulationProfiles=False,
)

generatedSpeciesConstraints(
    allowed=['input species','reaction libraries'],
    maximumSurfaceSites=2,
    maximumCarbonAtoms=2,
    maximumOxygenAtoms=2,
    maximumRadicalElectrons=1,
)

