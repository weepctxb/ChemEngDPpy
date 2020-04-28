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