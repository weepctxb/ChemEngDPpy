import numpy as np
import dsg

# CEPCI Index
# To access, e.g. CEPCI[2019]
CEPCI = {
    2019: 607.5,
    2018: 603.1,
    2001: 394.3
}

# USD to SGD forex rate (annualised average)
# To access, e.g. USSG[2019]
USSG = {
    2019: 1.3493,
    2018: 1.3912,
    2001: 1.7912
}

# Consumer price index
# To access, e.g. CPI['SG'][2019]
CPI = {
    'SG': {  # SG benchmark 2010 = 100
        2019: 115.0,
        2018: 113.8,
        2001: 86.05
    },
    'US': {  # US benchmark 1983 = 100
        2019: 257.0,
        2018: 251.2,
        2001: 176.7
    }
}


# Equipment cost correlation parameters
# A tuple of (K1, K2, K3, Amin, Amax, n)
# where (K1, K2, K3) = cost correlation params, (Amin, Amax) = min/max capacity (range of validity), n = cost exponent
# To access, e.g. eqptcostlib['pump']['centrifugal']
eqptcostlib = {
    'compressor': {
        'centrifugal': (2.2897, 1.3604, -0.1027, 450., 3000., 0.67),
        'axial': (2.2897, 1.3604, -0.1027, 450., 3000., 0.67),
        'reciprocating': (2.2897, 1.3604, -0.1027, 450., 3000., 0.84),
        'rotary': (5.0355, -1.8002, 0.8253, 18., 950., 0.6)
    },
    'pump': {
        'reciprocating': (3.8696, 0.3161, 0.1220, .1, 200., 0.6),
        'positivedisp': (3.4771, 0.1350, 0.1438, 1., 100., 0.6),
        'centrifugal': (3.3892, 0.0536, 0.1538, 1., 300., 0.67),
    },
    'heatexc': {
        'fixedtube': (4.3247, -0.3030, 0.1634, 10., 1000., 0.62),
        'Utube': (4.1884, -0.2503, 0.1974, 10., 1000., 0.53),
        'kettle': (4.4646, -0.5277, 0.3955, 10., 1000., 0.59),
        'doublepipe': (3.3444, 0.2745, -0.0472, 1., 10., 0.59),
        'multipipe': (2.7652, 0.7282,0.0783, 10., 100., 0.6)
    },
    'vessel': {
        'horizontal': (3.5565, 0.3776, 0.0905, 0.1, 628., 0.5),
        'vertical': (3.4974, 0.4485, 0.1074, 0.3, 520., 0.6)
    },
    'trays': {
        'sieve': (2.9949, 0.4465, 0.3961, 0.7, 12.3, 0.86),
        'valve': (3.3322, 0.4838, 0.3434, 0.7, 10.5, 1.0),
        'demister': (3.2353, 0.4838, 0.3434, 0.7, 10.5, 1.0)
    },
    'mixer': {
        'impeller': (3.8511, 0.7009, -0.0003, 5., 150., 0.6),
        'propeller': (4.3207, 0.0359, 0.1346, 5., 500., 0.5),
        'turbine': (3.4092, 0.4896, 0.0030, 5., 150., 0.3)
    }
}


# Equipment pressure factor correlation parameters
# A tuple of (C1, C2, C3, Pmin, Pmax)
# where (C1, C2, C3) = cost correlation params, (Pmin, Pmax) = min/max pressure (range of validity)
# To access, e.g. pressurefaclib['pump']['centrifugal']
pressurefaclib = {
    'compressor': {
        'centrifugal': (0., 0., 0., -np.inf, np.inf),
        'axial': (0., 0., 0., -np.inf, np.inf),
        'reciprocating': (0., 0., 0., -np.inf, np.inf),
        'rotary': (0., 0., 0., -np.inf, np.inf)
    },
    'pump': {
        'reciprocating': (-0.245382, 0.259016, -0.01363, 10., 100.),
        'positivedisp': (-0.245382, 0.259016, -0.01363, 10., 100.),
        'centrifugal': (-0.3935, 0.3957, -0.00226, 10., 100.),
    },
    'heatexc': {
        'fixedtube': (0.03881, -0.11272, 0.08183, 5., 140.),
        'Utube': (0.03881, -0.11272, 0.08183, 5., 140.),
        'kettle': (0.03881, -0.11272, 0.08183, 5., 140.),
        'doublepipe': (0.6072, -0.9120, 0.3327, 40., 100.),
        'multipipe': (0.6072, -0.9120, 0.3327, 40., 100.)
    },  # use the pressure factor equation for vessels instead
    'trays': {
        'sieve': (0., 0., 0., -np.inf, np.inf),
        'valve': (0., 0., 0., -np.inf, np.inf),
        'demister': (0., 0., 0., -np.inf, np.inf)
    },
    'mixer': {
        'impeller': (0., 0., 0., -np.inf, np.inf),
        'propeller': (0., 0., 0., -np.inf, np.inf),
        'turbine': (0., 0., 0., -np.inf, np.inf)
    }
}


# Equipment material factors
# To access, e.g. baremodlib['pump']['centrifugal']['SS']
matfaclib = {
    'compressor': {
        'centrifugal': {
            'CS': 2.8,  # CS = carbon steel
            'SS': 5.8 / 2.8,  # SS = stainless steel
            'Ni': 11.5 / 2.8  # Ni = nickel alloy
        },
        'axial': {
            'CS': 3.8,
            'SS': 8.0 / 3.8,
            'Ni': 15.9 / 3.8
        },
        'reciprocating': {
            'CS': 3.4,
            'SS': 7.0 / 3.4,
            'Ni': 13.9 / 3.4
        },
        'rotary': {
            'CS': 2.4,
            'SS': 5.0 / 2.4,
            'Ni': 9.9 / 2.4
        }
    },
    'pump': {
        'reciprocating': {
            'Fe': 1.0,  # Fe = cast iron
            'CS': 1.5,
            'SS': 2.4,
            'Ni': 4.0,
            'Ti': 6.5  # Ti = titanium alloy
        },
        'positivedisp': {
            'Fe': 1.0,
            'CS': 1.4,
            'SS': 2.7,
            'Ni': 4.7,
            'Ti': 10.7
        },
        'centrifugal': {
            'Fe': 1.0,
            'CS': 1.6,
            'SS': 2.3,
            'Ni': 4.4
        }
    },
    'heatexc': {
        HXtype: {
            'CS/CS': 1.0,
            'CS/SS': 1.8,
            'SS/SS': 2.9,
            'CS/Ni': 2.8,
            'Ni/Ni': 3.8,
            'CS/Ti': 4.6,
            'Ti/Ti': 11.4
        }
        for HXtype in ['fixedtube', 'Utube', 'kettle', 'doublepipe', 'multipipe']
    },
    'vessel': {
        vestype: {
            'CS': 1.0,
            'SS': 3.1,
            'Ni': 7.1,
            'Ti': 9.4
        }
        for vestype in ['horizontal', 'vertical']
    },
    'trays': {
        'sieve': {
            'CS': 1.0,
            'SS': 1.8,
            'Ni': 5.6
        },
        'valve': {
            'CS': 1.0,
            'SS': 1.8,
            'Ni': 5.6
        },
        'demister': {
            'SS': 1.0,
            'FC': 1.8,  # FC = fluorocarbon
            'Ni': 5.6
        }
    }
}

# Equipment bare module correlation parameters
# A tuple of (B1, B2)
# where (B1, B2) = bare module correlation params
# To access, e.g. baremodlib['pump']['centrifugal']
baremodlib = {
    'compressor': {  # FBM = (B2 for CS) * FM = (B2 for CS) * (B2 for material / B2 for CS)
        'centrifugal': (0., 2.8),
        'axial': (0., 3.8),
        'reciprocating': (0., 3.4),
        'rotary': (0., 2.4)
    },
    'pump': {  # FBM = B1 + B2 * FM * FP
        'reciprocating': (1.89, 1.35),
        'positivedisp': (1.89, 1.35),
        'centrifugal': (1.89, 1.35),
    },
    'heatexc': {  # FBM = B1 + B2 * FM * FP
        'fixedtube': (1.63, 1.66),
        'Utube': (1.63, 1.66),
        'kettle': (1.63, 1.66),
        'doublepipe': (1.74, 1.55),
        'multipipe': (1.74, 1.55)
    },
    'vessel': {  # FBM = B1 + B2 * FM * FP
        'horizontal': (1.49, 1.52),
        'vertical': (2.25, 1.82)
    },
    'trays': {  # FBM = FM for trays. Assuming tray quantity factor Fq = 1.
        'sieve': (0., 1.),
        'valve': (0., 1.),
        'demister': (0., 1.)
    },
    'mixer': {  # FBM = 1.38 (constant)
        'impeller': (1.38, 0.),
        'propeller': (1.38, 0.),
        'turbine': (1.38, 0.)
    }
}


def eqptpurcost(A, Ktuple):

    """
    Calculate equipment purchased cost (Cp^o) cost at ambient pressure and using carbon steel as MOC
    :param A: equipment capacity (various units)
    :param Ktuple: tuple of cost correlation factors (K1, K2, K3)
    :return: Cpo: equipment purchased cost ($)
    """

    Cpo = pow(10., Ktuple[0] + Ktuple[1] * np.log10(A) + Ktuple[2] * (np.log10(A)) ** 2)

    return Cpo


def pressurefacves(D, ts, P):

    """
    Calculate pressure factor (F_P) for vessels
    :param D: vessel diameter (m)
    :param ts: vessel thickness (in)
    :param P: pressure (barg)
    :return: FP: amplification factor for pressure
    """

    if P < -0.5:
        FP = 1.25
    elif P > -0.5 and ts < dsg.tmin:
        FP = 1
    else:
        FP = max(((P+1)*D / (2*(850-0.6*(P+1))) + 0.00315) / 0.0063, 1)

    return FP


def pressurefacanc(P, Ctuple):

    """
    Calculate pressure factor (F_P) for ancillary equipment (e.g. pumps and exchangers)
    at specified elevated pressure and MOC
    :param P: pressure (barg)
    :param Ctuple: tuple of pressure correlation factors (C1, C2, C3)
    :return: FP: amplification factor for pressure
    """

    FP = eqptpurcost(P, Ctuple)  # borrowing quadratic-exponential relation

    return FP


def baremodfac(Btuple, FM, FP):

    """
    Calculate bare module factor (F_BM) at specified elevated pressure and MOC
    :param Btuple: tuple of bare module correlation factors (B1, B2)
    :param FM: amplification factor for material of construction (MOC)
    :param FP: amplification factor for pressure
    :return: FBM: bare module factor (dimensionless)
    """

    FBM = Btuple[0] + Btuple[1] * FM * FP

    return FBM


def baremodcost(Cpo, FBM):

    """
    Calculate bare module cost (C_BM) at specified elevated pressure and MOC
    :param Cpo: equipment purchased cost ($)
    :param FBM: bare module factor (dimensionless)
    :return: CBM: bare module cost ($)
    """

    CBM = FBM * Cpo

    return CBM


def totmodcost(CBM):

    """
    Calculate total module cost (C_TM) at specified elevated pressure and MOC
    :param CBM: bare module cost ($)
    :return: CTM = total module cost ($)
    """

    CTM = 1.18 * CBM

    return CTM


def grasscost(CTM, Cpo):

    """
    Calculate grassroots cost (CGR)
    :param CTM: total module cost ($)
    :param Cpo: purchased equipment cost at ambient pressure and carbon steel MOC ($)
    :return: CGR: grassroots cost ($)
    """

    CGR = CTM + 0.5 * Cpo

    return CGR


def annualcapex(FCI, pbp=3):

    """
    Estimate total annualised capital cost based on assumed payback period
    :param FCI: fixed capital investment (= CTM or total module cost for brownfield projects, or =CGR or grassroots cost
    for greenfield projects) ($)
    :param pbp: payback period estimate (yr, default = 3). If a value of pbp is assumed, note that this should only be
    used for optimisation purposes! Alternatively, calculate pbp based on projected revenue estimates.
    :return: ACC: annualised capital cost estimate ($/yr)
    """

    ACC = FCI / pbp

    return ACC