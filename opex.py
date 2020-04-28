import numpy as np

runtime = 8000  # Operational runtime (h/yr)
shiftdur = 8  # Duration per workshift (h/shift)
shiftperwk = 5  # Number of workshifts per week
yearww = 49  # Number of work weeks per year
SF = runtime / (365*24)  # Stream factor


# Utility Costs (per GJ basis)
# From "Analysis, Synthesis & Design of Chemical Processes, 5th Ed. by Turton et al."
util = {
    'LPS': 4.54,  # Utility cost for LPS (5 barg, 160 degC) ($/GJ)
    'MPS': 4.77,  # Utility cost for MPS (10 barg, 184 degC) ($/GJ)
    'HPS': 5.66,  # Utility cost for HPS (41 barg, 254 degC) ($/GJ)
    'CW': 0.378,  # Utility cost for cooling water (30-45 degC) ($/GJ)
    'ChW': 4.77,  # Utility cost for chilled water (5 degC) ($/GJ)
    'LTR': 8.49,  # Utility cost for low temperature refrigerant (-20 degC) ($/GJ)
    'VLTR': 14.12,  # Utility cost for very low temperature refrigerant (-50 degC) ($/GJ)
    'elec': 18.72  # Utility cost for electricity (110-440 V) ($/GJ)
}


def operatorspershift(P=0, Nnp=0):

    """
    Calculate number of operators per shift
    :param P: number of processing steps involving particulate solids (P=0 for fluid-processing plants)
    :param Nnp: number of non-particulate/fluid handling equipment/steps (include compressors, towers, reactors,
    heaters and exchangers; exclude pumps, vessels and tanks)
    :return: NOL: number of operators required per shift
    """

    NOL = round(np.sqrt(6.29 + 31.7 * P ** 2 + 0.23 * Nnp))

    return NOL


def labourcost(P=0, Nnp=0, wage=0):

    """
    Calculate annualised labour cost
    :param P: number of processing steps involving particulate solids (P=0 for fluid-processing plants)
    :param Nnp: number of non-particulate/fluid handling equipment/steps (include compressors, towers, reactors,
    heaters and exchangers; exclude pumps, vessels and tanks)
    :param wage: annualised per-operator wage ($/yr)
    :return: COL: annualised labour cost ($/yr)
    :return: NOL: total number of operators required
    """

    NOL = operatorspershift(P, Nnp)
    shiftperyr = runtime / shiftdur
    shiftperopperyr = yearww * shiftperwk
    Nop = round(shiftperyr / shiftperopperyr * NOL)
    COL = Nop * wage

    return COL, NOL


def costofutil(HPS=0, MPS=0, LPS=0, CW=0, ChW=0, LTR=0, VLTR=0, elec=0):

    """
    Calculates the annualised cost of utilities
    :param HPS: annual consumption of high-pressure steam (GJ/yr)
    :param MPS: annual consumption of medium-pressure steam (GJ/yr)
    :param LPS: annual consumption of low-pressure steam (GJ/yr)
    :param CW: annual consumption of cooling water (GJ/yr)
    :param ChW: annual consumption of chilled water (GJ/yr)
    :param LTR: annual consumption of low-temperature refrigerant (GJ/yr)
    :param VLTR: annual consumption of very low-temperature refrigerant (GJ/yr)
    :param elec: annual consumption of electricity (GJ/yr)
    :return: CUT: annualised cost of utilities ($/yr)
    """

    a = np.array([util['HPS'], util['MPS'], util['LPS'],
         util['CW'], util['ChW'],
         util['LTR'], util['VLTR'], util['elec']])

    b = np.array([HPS, MPS, LPS, CW, ChW, LTR, VLTR, elec])

    CUT = a @ b.T

    return CUT


def costofraw(rawmaterials, unitprices):

    """
    Calculates the annualised cost of raw materials
    :param rawmaterials: annual flow of raw materials, in a tuple (flow/yr)
    :param unitprices: per-flow unit cost price of raw materials, in a tuple in the same order as rawmaterials ($/flow)
    :return: CRM: annualised cost of raw materials ($/yr)
    """

    if type(rawmaterials) is tuple:
        a = np.array(rawmaterials)
    if type(unitprices) is tuple:
        b = np.array(unitprices)

    if a.size != b.size:
        raise ValueError('Number of unit prices do not match number of raw materials!')

    CRM = a @ b.T

    return CRM


def costofmanfc(FCI, COL, CRM, CWT, CUT):

    """
    Estimate annualised cost of manufacture.
    :param FCI: fixed capital investment (= CTM or total module cost for brownfield projects, or =CGR or grassroots cost
    for greenfield projects) ($)
    :param COL: annualised labour cost ($/yr)
    :param CRM: annualised cost of raw materials ($/yr)
    :param CWT: annualised cost of waste treatment ($/yr)
    :param CUT: annualised cost of utilities ($/yr)
    :return:
    """

    COMd = 0.18 * FCI + 2.73 * COL + 1.23 * (CRM + CWT + CUT)
    d = 0.1 * FCI
    COM = COMd + d
    DMC = CRM + CWT + CUT + 1.33 * COL + 0.069 * FCI + 0.03 * COM
    FMC = 0.708 * COL + 0.068 * FCI
    GE = 0.177 * COL + 0.009 * FCI + 0.16 * COM

    return COMd, d, COM, DMC, FMC, GE