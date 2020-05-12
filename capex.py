import numpy as np
import dsg
import warnings
import time
from typing import Tuple, Any, List

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
        2016: 112.6,
        2001: 86.05
    },
    'US': {  # US benchmark 1983 = 100
        2019: 257.0,
        2018: 251.2,
        2016: 240.0,
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
        'utube': (4.1884, -0.2503, 0.1974, 10., 1000., 0.53),
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
# where (C1, C2, C3) = pressure factor correlation params, (Pmin, Pmax) = min/max pressure (range of validity)
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
        'utube': (0.03881, -0.11272, 0.08183, 5., 140.),
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
# To access, e.g. matfaclib['pump']['centrifugal']['SS']
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
            'SS/CS': 1.8,  # duplicate
            'SS/SS': 2.9,
            'CS/Ni': 2.8,
            'Ni/CS': 2.8,  # duplicate
            'Ni/Ni': 3.8,
            'CS/Ti': 4.6,
            'Ti/CS': 4.6,  # duplicate
            'Ti/Ti': 11.4
        }
        for HXtype in ['fixedtube', 'utube', 'kettle', 'doublepipe', 'multipipe']
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
        'utube': (1.63, 1.66),
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


def eqptpurcost(A: float=None, Ktuple: Tuple[float]=None, eqpt: Any=None) -> (float, Any):

    """
    Calculate equipment purchased cost (Cp^o) cost at ambient pressure and using carbon steel as MOC
    Two methods of calculation:
    Method 1 - Specify A and Ktuple manually:
    :param A: equipment capacity (various units)
    :param Ktuple: tuple of cost correlation factors (K1, K2, K3, Amin [optional], Amax [optional], n [optional])
    Optional inputs within Ktuple refer to minimum/maximum capacity and exponential factor respectively
    :return: Cpo: equipment purchased cost ($)
    Method 2 - Specify the equipment object directly
    :param eqpt: equipment object as generated by the dsg.size(...) or dsg.design(...) functions
    :return: Cpo: equipment purchased cost ($)
    :return: eqpt: the same equipment object with Cpo updated
    """

    # Method 1 - Specify A and Ktuple manually
    if A is not None and Ktuple is not None:
        pass

    # Method 2 - Specify the equipment object directly
    elif eqpt is not None:
        # For vessel convert in^3 to m^3
        try:
            A = eqpt.comppower if eqpt.category is 'compressor' \
                else eqpt.pumppower if eqpt.category is 'pump' \
                else eqpt.area if eqpt.category is 'heatexc' \
                else (eqpt.Vi*1.639e-5) if eqpt.category is 'vessel' \
                else eqpt.area if eqpt.category is 'trays' \
                else eqpt.mixerpower if eqpt.category is 'mixer' \
                else None
            Ktuple = eqptcostlib[eqpt.category][eqpt.etype]
        except KeyError:
            raise KeyError('Equipment category and/or type (eqpt.category and/or eqpt.etype) not supported!')

    else:
        raise ValueError('Specify either (A + Ktuple) or eqpt!')

    if len(Ktuple) == 6:
        if A < Ktuple[3]:
            warnings.warn('Extrapolating {}={} below minimum capacity of {}! Switching to exponential rule!'.format(eqpt.id, A, Ktuple[3]))
            Cpo = pow(10., Ktuple[0] + Ktuple[1] * np.log10(Ktuple[3]) + Ktuple[2] * (np.log10(Ktuple[3])) ** 2)
            Cpo *= pow(A / Ktuple[3], Ktuple[5])
        elif A > Ktuple[4]:
            warnings.warn('Extrapolating {}={} above maximum capacity of {}! Switching to exponential rule!'.format(eqpt.id, A, Ktuple[4]))
            Cpo = pow(10., Ktuple[0] + Ktuple[1] * np.log10(Ktuple[4]) + Ktuple[2] * (np.log10(Ktuple[4])) ** 2)
            Cpo *= pow(A / Ktuple[4], Ktuple[5])
        else:
            Cpo = pow(10., Ktuple[0] + Ktuple[1] * np.log10(A) + Ktuple[2] * (np.log10(A)) ** 2)
    else:
        Cpo = pow(10., Ktuple[0] + Ktuple[1] * np.log10(A) + Ktuple[2] * (np.log10(A)) ** 2)

    if eqpt is None:
        return Cpo
    else:
        eqpt.Cpo = Cpo
        return Cpo, eqpt


def pressurefacves(D: float=None, ts: float=None, P: float=None, eqpt: Any=None) -> (float, Any):

    """
    Calculate pressure factor (F_P) for vessels
    Two methods of calculation:
    Method 1 - Specify D, ts and P manually:
    :param D: vessel internal diameter (m)
    :param ts: vessel thickness (in)
    :param P: pressure (barg)
    :return: FP: amplification factor for pressure
    Method 2 - Specify the equipment object directly:
    :param eqpt: equipment object as generated by the dsg.design(...) functions
    :return: FP: amplification factor for pressure
    :return: eqpt: the same equipment object with FP updated
    """

    # Method 1 - Specify A and Ktuple manually
    if D is not None and ts is not None and P is not None:
        pass

    # Method 2 - Specify the equipment object directly
    elif eqpt is not None:

        if eqpt.category is not 'vessel':
            raise ValueError('Use pressurefacanc(P, Ctuple, eqpt) for ancillary equipment instead!')

        D = eqpt.Di / 39.37  # convert inch to m
        ts = eqpt.ts  # in inch
        P = eqpt.Pd * 0.06895  # convert psig to barg

    if P < -0.5:
        FP = 1.25
    elif P > -0.5 and ts < dsg.tmin:
        FP = 1
    else:
        FP = max(((P+1)*D / (2*(850-0.6*(P+1))) + 0.00315) / 0.0063, 1)

    if eqpt is None:
        return FP
    else:
        eqpt.FP = FP
        return FP, eqpt


def pressurefacanc(P: float=None, Ctuple: Tuple[float]=None, eqpt: Any=None) -> (float, Any):

    """
    Calculate pressure factor (F_P) for ancillary equipment (e.g. pumps and exchangers)
    at specified elevated pressure and MOC
    Two methods of calculation:
    Method 1 - Specify P and Ctuple manually:
    :param P: pressure (barg)
    :param Ctuple: tuple of pressure correlation factors (C1, C2, C3, Pmin [optional], Pmax [optional])
    Optional inputs within Ctuple refer to minimum/maximum pressure respectively
    :return: FP: amplification factor for pressure
    Method 2 - Specify the equipment object directly:
    :param eqpt: equipment object as generated by the dsg.design(...) functions
    :return: FP: amplification factor for pressure
    :return: eqpt: the same equipment object with FP updated
    """

    # Method 1 - Specify P and Ctuple manually
    if P is not None and Ctuple is not None:
        pass

    # Method 2 - Specify the equipment object directly
    elif eqpt is not None:

        if eqpt.category is 'vessel':
            raise ValueError('Use pressurefacves(D, ts, P, eqpt) for vessels instead!')

        try:
            # Compressor: P2 in bar, but assume pressure doesn't affect bare module factor anyways
            # Pump: P2 in kPa
            # Heat Exchanger: P in bar
            # Vessel: Pd in psig, but vessel pressure calculations is not done here anyways
            # Trays: Assume 1 atm, dummy pressure because pressure doesn't affect bare module factor
            # Mixer: Assume 1 atm, dummy pressure because pressure doesn't affect bare module factor
            P = (eqpt.P2 - dsg.Patmb) if eqpt.category is 'compressor' \
                else (eqpt.P2 - dsg.Patmb*100.)/100. if eqpt.category is 'pump' \
                else (eqpt.P - dsg.Patmb) if eqpt.category is 'heatexc' \
                else (eqpt.Pd * 0.06895) if eqpt.category is 'vessel' \
                else dsg.Patmb if eqpt.category is 'trays' \
                else dsg.Patmb if eqpt.category is 'mixer' \
                else None
            Ctuple = pressurefaclib[eqpt.category][eqpt.etype]
        except KeyError:
            raise KeyError('Equipment category and/or type (eqpt.category and/or eqpt.etype) not supported!')

    else:
        raise ValueError('Specify either (P + Ctuple) or eqpt!')

    if len(Ctuple) == 5:
        if P < Ctuple[3]:
            warnings.warn('Equipment {} pressure={} below minimum pressure of {} for pressure factor correlation! Using FP = 1 instead!'.format(eqpt.id, P, Ctuple[3]))
            FP = 1.  # borrowing quadratic-exponential relation
        elif P > Ctuple[4]:
            warnings.warn('Equipment {} pressure={} above maximum pressure of {} for pressure factor correlation! Using max P instead!'.format(eqpt.id, P, Ctuple[3]))
            FP = max(1., eqptpurcost(A=Ctuple[4], Ktuple=tuple(Ctuple[0:3])))  # borrowing quadratic-exponential relation
        else:
            FP = max(1., eqptpurcost(A=P, Ktuple=Ctuple[0:3]))  # borrowing quadratic-exponential relation
    else:
        FP = max(1., eqptpurcost(A=P, Ktuple=Ctuple[0:3]))  # borrowing quadratic-exponential relation

    if eqpt is None:
        return FP
    else:
        eqpt.FP = FP
        return FP, eqpt


def baremodfac(Btuple: Tuple[float]=None, FM: float=None, FP: float=None, eqpt: Any=None) -> (float, Any):

    """
    Calculate bare module factor (F_BM) at specified elevated pressure and MOC
    Two methods of calculation:
    Method 1 - Specify Btuple, FM and FP manually:
    :param Btuple: tuple of bare module correlation factors (B1, B2)
    :param FM: amplification factor for material of construction (MOC)
    :param FP: amplification factor for pressure
    :return: FBM: bare module factor (dimensionless)
    Method 2 - Specify the equipment object directly:
    :param eqpt: equipment object as generated by the dsg.design(...) or dsg.size(...) functions
    :return: FBM: bare module factor (dimensionless)
    :return: eqpt: the same equipment object with FBM updated
    """

    # Method 1 - Specify Btuple, FM and FP manually
    if Btuple is not None and FM is not None and FP is not None:
        pass

    # Method 2 - Specify the equipment object directly
    elif eqpt is not None:
        try:
            Btuple = baremodlib[eqpt.category][eqpt.etype]
            FM = matfaclib[eqpt.category][eqpt.etype][eqpt.mat]
            if eqpt.category is 'vessel':
                FP, eqpt = pressurefacves(eqpt=eqpt)
            else:
                FP, eqpt = pressurefacanc(eqpt=eqpt)
        except KeyError:
            raise KeyError('Equipment category and/or type and/or material ' +
                           '(eqpt.category and/or eqpt.etype and/or eqpt.mat) not supported!')

    else:
        raise ValueError('Specify either (Btuple + FM + FP) or eqpt!')

    FBM = max(1., Btuple[0] + Btuple[1] * FM * FP)

    if eqpt is None:
        return FBM
    else:
        eqpt.FBM = FBM
        eqpt.FM = FM
        return FBM, eqpt


def baremodcost(Cpo=None, FBM=None, eqpt=None):

    """
    Calculate bare module cost (CBM) at specified elevated pressure and MOC
    Two methods of calculation:
    Method 1 - Specify Cpo and FBM manually:
    :param Cpo: equipment purchased cost ($)
    :param FBM: bare module factor (dimensionless)
    :return: CBM: bare module cost ($)
    Method 2 - Specify the equipment object directly:
    :param eqpt: equipment object as generated by the dsg.design(...) or dsg.size(...) functions
    :return: CBM: bare module cost ($)
    :return: eqpt: the same equipment object with CBM updated
    """

    # Method 1 - Specify Cpo and FBM manually:
    if Cpo is not None and FBM is not None:
        pass

    # Method 2 - Specify the equipment object directly:
    elif eqpt is not None:
        Cpo, eqpt = eqptpurcost(eqpt=eqpt)
        FBM, eqpt = baremodfac(eqpt=eqpt)

    else:
        raise ValueError('Specify either (Cpo + FBM) or eqpt!')

    CBM = FBM * Cpo

    if eqpt is None:
        return CBM
    else:
        eqpt.CBM = CBM
        return CBM, eqpt


def totmodcost(CBM: float=None, eqpt: float=None) -> (float, Any):

    """
    Calculate total module cost (CTM) at specified elevated pressure and MOC
    Two methods of calculation:
    Method 1 - Specify CBM manually:
    :param CBM: bare module cost ($)
    :return: CTM = total module cost ($)
    Method 2 - Specify the equipment object directly:
    :param eqpt: equipment object as generated by the dsg.design(...) or dsg.size(...) functions
    :return: CTM = total module cost ($)
    :return: eqpt: the same equipment object with CTM updated
    """

    # Method 1 - Specify CBM manually:
    if CBM is not None:
        pass

    # Method 2 - Specify the equipment object directly:
    elif eqpt is not None:
        CBM, eqpt = baremodcost(eqpt=eqpt)

    else:
        raise ValueError('Specify either CBM or eqpt!')

    CTM = 1.18 * CBM

    if eqpt is None:
        return CTM
    else:
        eqpt.CTM = CTM
        return CTM, eqpt


def grasscost(CTM: float=None, Cpo: float=None, eqpt: Any=None) -> (float, Any):

    """
    Calculate grassroots cost (CGR)
    Two methods of calculation:
    Method 1 - Specify CTM and Cpo manually:
    :param CTM: total module cost ($)
    :param Cpo: purchased equipment cost at ambient pressure and carbon steel MOC ($)
    :return: CGR: grassroots cost ($)
    Method 2 - Specify the equipment object directly:
    :param eqpt: equipment object as generated by the dsg.design(...) or dsg.size(...) functions
    :return: CGR: grassroots cost ($)
    :return: eqpt: the same equipment object with CTM updated
    """

    # Method 1 - Specify CTM and Cpo manually:
    if CTM is not None and Cpo is not None:
        pass

    # Method 2 - Specify the equipment object directly:
    elif eqpt is not None:
        CTM, eqpt = totmodcost(eqpt=eqpt)  # this will update eqpt.Cpo too if all goes well
        Cpo = eqpt.Cpo

    else:
        raise ValueError('Specify either CBM or eqpt!')

    CGR = CTM + 0.5 * Cpo

    if eqpt is None:
        return CGR
    else:
        eqpt.CGR = CGR
        return CGR, eqpt


def annualcapex(FCI: float=None, pbp: float=None, eqpt: Any=None, planttype: str='brown') -> (float, Any):

    """
    Estimate total annualised capital cost based on assumed payback period
    Two methods of calculation:
    Method 1 - Specify FCI and pbp manually:
    :param FCI: fixed capital investment (= CTM or total module cost for brownfield projects, or =CGR or grassroots cost
    for greenfield projects) ($)
    :param pbp: payback period estimate (yr, default = 3). If a value of pbp is assumed, note that this should only be
    used for optimisation purposes! Alternatively, calculate pbp based on projected revenue estimates.
    :return: ACC: annualised capital cost estimate ($/yr)
    Method 2 - Specify eqpt directly, then pbp and planttype:
    :param eqpt: equipment object as generated by the dsg.design(...) or dsg.size(...) functions
    :param pbp: payback period estimate, as described in Method 1
    :param planttype: 'brown' for brownfield project (using CTM) or 'green' for greenfield project (using CGR)
    :return: ACC: annualised capital cost estimate ($/yr)
    """

    if pbp is None:
        pbp = 3
        warnings.warn('Payback period (pbp) not specified - 3 years is assumed')

    if FCI is not None:
        pass

    elif eqpt is not None:
        if 'brown' in str.lower(planttype):
            _, eqpt = grasscost(eqpt=eqpt)
            FCI = eqpt.CTM
        elif 'green' in str.lower(planttype):
            FCI, eqpt = grasscost(eqpt=eqpt)
        else:
            raise ValueError('planttype should be \'brown\' for brownfield project (default) ' +
                             'or \'green\' for greenfield project!')

    ACC = FCI / pbp

    if eqpt is None:
        return ACC
    else:
        eqpt.ACC = ACC
        return ACC, eqpt


def econreport(eqptlist: List[Any], planttype: str='green', reporttype: str='numpy',
               pbp: float=3., year: int=2019, currency: str='SGD', verbose: bool=False) -> (float, Any):

    """
    Generates an economic capex report of the plant
    :param eqptlist: List of equipment objects as generated by the dsg.design(...) or dsg.size(...) functions
    :param planttype: Type of project ('green' for greenfield [default] or 'brown' for brownfield)
    :param reporttype: Data structure of report ('list' for 2D list, 'dict' for dictionary, or numpy for numpy array [default])
    :param pbp: payback period estimate (yr, default = 3). If a value of pbp is assumed, note that this should only be
    used for optimisation purposes! Alternatively, calculate pbp based on projected revenue estimates.
    :param year: Year for CEPCI updating (integer, either 2001, 2018 or 2019 [default])
    :param currency: Currency (string, either 'USD' or 'SGD' [default])
    :param verbose: True to print economic capex report, False to print nothing [default]
    :return: report: The capex report formatted according to the reporttype input ('list', 'dict' or 'numpy')
    """
    report_list = [['Equipment'], ['Purchased eqpt. cost (Cpo)'], ['Bare mod. cost (CBM)'], ['Total mod. cost (CTM)'],
                   ['Grassroots cost (CGR)'], ['Annualised capital cost (ACC)']]
    report_dict = {eqpt.id: {'Cpo': None, 'CBM': None, 'CTM': None, 'CGR': None, 'ACC': None} for eqpt in eqptlist}
    report_dict['Total'] = {'Cpo': None, 'CBM': None, 'CTM': None, 'CGR': None, 'ACC': None}

    yearcurrfac = CEPCI[year] / CEPCI[2001] * (USSG[year] if currency is 'SGD' else 1.)

    for eqpt in eqptlist:
        ACC, eqpt = annualcapex(pbp=pbp, eqpt=eqpt, planttype=planttype)

        eqpt.Cpo = round(eqpt.Cpo * yearcurrfac, 2)
        eqpt.CBM = round(eqpt.CBM * yearcurrfac, 2)
        eqpt.CTM = round(eqpt.CTM * yearcurrfac, 2)
        eqpt.CGR = round(eqpt.CGR * yearcurrfac, 2)
        eqpt.ACC = round(eqpt.ACC * yearcurrfac, 2)

        i = list.index(eqptlist, eqpt)
        report_list[0].append(eqpt.id)
        report_list[1].append(eqpt.Cpo)
        report_list[2].append(eqpt.CBM)
        report_list[3].append(eqpt.CTM)
        report_list[4].append(eqpt.CGR)
        report_list[5].append(eqpt.ACC)

        report_dict[eqpt.id]['Cpo'] = eqpt.Cpo
        report_dict[eqpt.id]['CBM'] = eqpt.CBM
        report_dict[eqpt.id]['CTM'] = eqpt.CTM
        report_dict[eqpt.id]['CGR'] = eqpt.CGR
        report_dict[eqpt.id]['ACC'] = eqpt.ACC

        if verbose:
            time.sleep(0.1)
            print(eqpt.spec())
            print(eqpt.econ())

    report_list[0].append('Total')
    for i in range(1, 5+1):
        report_list[i].append(round(sum(report_list[i][1:len(eqptlist)+1]), 2))

    report_numpy = np.array(report_list)

    report_dict['Total']['Cpo'] = round(sum([eqpt.Cpo for eqpt in eqptlist]), 2)
    report_dict['Total']['CBM'] = round(sum([eqpt.CBM for eqpt in eqptlist]), 2)
    report_dict['Total']['CTM'] = round(sum([eqpt.CTM for eqpt in eqptlist]), 2)
    report_dict['Total']['CGR'] = round(sum([eqpt.CGR for eqpt in eqptlist]), 2)
    report_dict['Total']['ACC'] = round(sum([eqpt.ACC for eqpt in eqptlist]), 2)

    FCI = round(report_dict['Total']['CGR'] if planttype is 'green' else report_dict['Total']['CTM'], 2)

    if verbose:
        time.sleep(0.1)
        print('----------------------------')
        print('TOTAL PLANT COST (' + str.upper(planttype) + 'FIELD): $' + str(FCI))
        print('----------------------------')
        print('CAPEX REPORT:')
        print(report_list if reporttype is 'list' else report_dict if reporttype is 'dict' else report_numpy)
        print('----------------------------')

    return FCI, report_list if reporttype is 'list' else report_dict if reporttype is 'dict' else report_numpy
