"""Shared constants used by the sequence utilities and assembly planning code."""

import math

# Thermodynamic parameters for nearest-neighbor Tm estimation.
NN_THERMO_PARAMS = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

NN_INIT_AT = (2.3, 4.1)
NN_INIT_GC = (0.1, -2.8)

GAS_CONSTANT_CAL = 1.987

DEFAULT_OLIGO_CONC_M = 250e-9
DEFAULT_MONOVALENT_CONC_M = 50e-3
