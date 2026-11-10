"""Physical constants and atomic weights. Units: mol, J, K, Pa unless a name says otherwise."""

R = 8.31446261815324  # gas constant [J/mol/K]
P0 = 100000.0  # standard pressure [Pa]
T_REF = 298.15  # reference temperature of formation data [K]
P_REF = 1e5  # reference pressure of SUPCRT mineral volumes [Pa]
CAL_TO_J = 4.184
MW_WATER_KG = 0.01801528  # molar mass of water used in molalities [kg/mol]

# Standard atomic weights [g/mol]. Cases may override individual values
# (Case.atomic_weights); this changes feed and balance arithmetic only.
ATOMIC_WEIGHTS = {
    'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007, 'S': 32.06, 'Na': 22.989769,
    'He': 4.002602, 'Li': 6.94, 'Be': 9.0121831, 'B': 10.81, 'F': 18.998403163,
    'Ne': 20.1797, 'Mg': 24.305, 'Al': 26.9815385, 'Si': 28.085, 'P': 30.973761998,
    'Cl': 35.45, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908, 'Ti': 47.867,
    'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938044, 'Fe': 55.845, 'Co': 58.933194,
    'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.630, 'As': 74.921595,
    'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90584,
    'Zr': 91.224, 'Nb': 92.90637, 'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.90550,
    'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.710,
    'Sb': 121.760, 'Te': 127.60, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.90545196,
    'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 'Pr': 140.90766, 'Nd': 144.242,
    'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.500,
    'Ho': 164.93033, 'Er': 167.259, 'Tm': 168.93422, 'Yb': 173.045, 'Lu': 174.9668,
    'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84, 'Re': 186.207, 'Os': 190.23,
    'Ir': 192.217, 'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.592, 'Tl': 204.38,
    'Pb': 207.2, 'Bi': 208.98040, 'Th': 232.0377, 'U': 238.02891,
}

PHASES = ('gas', 'aqueous', 'organic_liquid', 'solid')
