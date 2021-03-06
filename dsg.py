import numpy as np
import warnings
import capex

Patm = 14.696  # Standard atmospheric pressure (psi)
Patmb = 1.01325  # Standard atmospheric pressure (bar)
Troom = 77.  # Ambient temperature (degF)
tmin = 1/4  # Universal minimum allowable vessel thickness (in)
tc = 0.125  # Corrosion allowance (in) for both corrosive and non-corrosive conditions (default is 1/8)
rhosteel = 0.2836  # Density of SA-285C/SA-387B/carbon/low-alloy steels (lb/in^3)

g = 9.80665  # standard Earth gravitational acceleration (m/s^2)
R = 8.31446261815324  # universal ideal gas constant (J/(K.mol))
Ta = 10.  # minimum heat exchanger temperature approach (K)


class MechDesign(object):
    def __init__(self, Po=Patm, To=Troom, Di=None, L=None, rho=rhosteel,
                 Pd=None, Td=None, MOC=None, Smax=None, E=0.85, tp=tmin, tc=None, ts=None, tsfinal=None,
                 tv=None, tw=None,
                 Do=None, W=None, V=None, Vi=None,
                 EM=None, tE=None, tEC=None,
                 category='vessel', etype=None, mat=None, id=None,
                 Cpo=None, FP=None, FM=None, FBM=None, CBM=None, CTM=None, CGR=None, ACC=None):
        self.Po = Po
        self.To = To
        self.Di = Di
        self.L = L
        self.rho = rho
        self.Pd = Pd
        self.Td = Td
        self.MOC = MOC
        self.Smax = Smax
        self.E = E
        self.tp = tp
        self.tc = tc
        self.ts = ts
        self.tsfinal = tsfinal
        self.tv = tv
        self.tw = tw
        self.Do = Do
        self.W = W
        self.V = V
        self.Vi = Vi
        self.EM = EM
        self.tE = tE
        self.tEC = tEC
        self.category = category
        self.etype = etype
        self.mat = mat
        self.id = id
        self.Cpo = Cpo
        self.FP = FP
        self.FM = FM
        self.FBM = FBM
        self.CBM = CBM
        self.CTM = CTM
        self.CGR = CGR
        self.ACC = ACC

    def __repr__(self):
        return '%s: MechDesign(Pd=%spsig, Td=%sdegF, tsfinal=%sin, L=%sin, Do=%sin, W=%slb, V=%sin^3, etype=%s, mat=%s)' \
               % (self.id, round(self.Pd, 2), self.Td, self.tsfinal, self.L, self.Do, int(self.W), int(self.V), self.etype, self.mat)

    def spec(self):
        return '---Design specs for {}:---\n{}\n----------------------------'.format(self.id, vars(self))

    def econ(self):
        return '---Econ report for {}:---\nCpo=${}\nCBM=${}\nCTM=${}\nCGR=${}\nACC=${}\n----------------------------'.format(self.id, round(self.Cpo, 2), round(self.CBM, 2), round(self.CTM, 2), round(self.CGR, 2), round(self.ACC, 2))


class Mixer(object):
    def __init__(self, mixerpower=None,
                 category='mixer', etype=None, mat=None, id=None,
                 Cpo=None, FP=None, FM=None, FBM=None, CBM=None, CTM=None, CGR=None, ACC=None):
        self.mixerpower = mixerpower
        self.category = category
        self.etype = etype
        self.mat = mat
        self.id = id
        self.Cpo = Cpo
        self.FP = FP
        self.FM = FM
        self.FBM = FBM
        self.CBM = CBM
        self.CTM = CTM
        self.CGR = CGR
        self.ACC = ACC

    def __repr__(self):
        return '%s: Mixer(mixerpower=%skW, etype=%s, mat=%s)' \
               % (self.id, self.mixerpower, self.etype, self.mat)

    def spec(self):
        return '---Design specs for {}:---\n{}\n----------------------------'.format(self.id, vars(self))

    def econ(self):
        return '---Econ report for {}:---\nCpo=${}\nCBM=${}\nCTM=${}\nCGR=${}\nACC=${}\n----------------------------'.format(self.id, round(self.Cpo, 2), round(self.CBM, 2), round(self.CTM, 2), round(self.CGR, 2), round(self.ACC, 2))


class Trays(object):
    def __init__(self, numtrays=None, area=None,
                 category='trays', etype=None, mat=None, id=None,
                 Cpo=None, FP=None, FM=None, FBM=None, CBM=None, CTM=None, CGR=None, ACC=None):
        self.numtrays = numtrays
        self.area = area
        self.category = category
        self.etype = etype
        self.mat = mat
        self.id = id
        self.Cpo = Cpo
        self.FP = FP
        self.FM = FM
        self.FBM = FBM
        self.CBM = CBM
        self.CTM = CTM
        self.CGR = CGR
        self.ACC = ACC

    def __repr__(self):
        return '%s: Trays(numtrays=%s, area=%s, etype=%s, mat=%s)' \
               % (self.id, self.numtrays, self.area, self.etype, self.mat)

    def spec(self):
        return '---Design specs for {}:---\n{}\n----------------------------'.format(self.id, vars(self))

    def econ(self):
        return '---Econ report for {}:---\nCpo=${}\nCBM=${}\nCTM=${}\nCGR=${}\nACC=${}\n----------------------------'.format(self.id, round(self.Cpo, 2), round(self.CBM, 2), round(self.CTM, 2), round(self.CGR, 2), round(self.ACC, 2))


class Reactor(MechDesign, Mixer):
    def __init__(self):
        super().__init__(self)

    def __repr__(self):
        return super().__repr__(self)


class Distillation(MechDesign, Trays):
    def __init__(self):
        super().__init__(self)

    def __repr__(self):
        return super().__repr__(self)


class Compressor(object):
    def __init__(self, m=None, P1=Patm, P2=None, T1=Troom, T2=None, cp=None, cv=None, Z=None,
                 compeff=None, comppower=None,
                 category='compressor', etype=None, mat=None, id=None,
                 Cpo=None, FP=None, FM=None, FBM=None, CBM=None, CTM=None, CGR=None, ACC=None):
        self.m = m
        self.P1 = P1
        self.P2 = P2
        self.T1 = T1
        self.T2 = T2
        self.cp = cp
        self.cv = cv
        self.Z = Z
        self.compeff = compeff
        self.comppower = comppower
        self.category = category
        self.etype = etype
        self.mat = mat
        self.id = id
        self.Cpo = Cpo
        self.FP = FP
        self.FM = FM
        self.FBM = FBM
        self.CBM = CBM
        self.CTM = CTM
        self.CGR = CGR
        self.ACC = ACC

    def __repr__(self):
        return '%s: Compressor(P1=%sbar, P2=%sbar, compeff=%s, comppower=%skW, etype=%s, mat=%s)' \
               % (self.id, self.P1, self.P2, round(self.compeff, 3), round(self.comppower, 3), self.etype, self.mat)

    def spec(self):
        return '---Design specs for {}:---\n{}\n----------------------------'.format(self.id, vars(self))

    def econ(self):
        return '---Econ report for {}:---\nCpo=${}\nCBM=${}\nCTM=${}\nCGR=${}\nACC=${}\n----------------------------'.format(self.id, round(self.Cpo, 2), round(self.CBM, 2), round(self.CTM, 2), round(self.CGR, 2), round(self.ACC, 2))


class Pump(object):
    def __init__(self, Q=None, P1=None, P2=None, dP=None, rho=1000,
                 pumpeff=0.75, pumppower=None,
                 category='pump', etype=None, mat=None, id=None,
                 Cpo=None, FP=None, FM=None, FBM=None, CBM=None, CTM=None, CGR=None, ACC=None):
        self.Q = Q
        self.P1 = P1
        self.P2 = P2
        self.dP = dP
        self.rho = rho
        self.pumpeff = pumpeff
        self.pumppower = pumppower
        self.category = category
        self.etype = etype
        self.mat = mat
        self.id = id
        self.Cpo = Cpo
        self.FP = FP
        self.FM = FM
        self.FBM = FBM
        self.CBM = CBM
        self.CTM = CTM
        self.CGR = CGR
        self.ACC = ACC

    def __repr__(self):
        return '%s: Pump(P1=%skPa, P2=%skPa, pumpeff=%s, pumppower=%skW, etype=%s, mat=%s)' \
               % (self.id, self.P1, self.P2, round(self.pumpeff, 3), round(self.pumppower, 3), self.etype, self.mat)

    def spec(self):
        return '---Design specs for {}:---\n{}\n----------------------------'.format(self.id, vars(self))

    def econ(self):
        return '---Econ report for {}:---\nCpo=${}\nCBM=${}\nCTM=${}\nCGR=${}\nACC=${}\n----------------------------'.format(self.id, round(self.Cpo, 2), round(self.CBM, 2), round(self.CTM, 2), round(self.CGR, 2), round(self.ACC, 2))


class HeatExc(object):
    def __init__(self, mh=None, mc=None, cph=None, cpc=None, Thin=None, Thout=None, Tcin=None, Tcout=None,
                 U=None, F=0.9, Ns=1, area=None, P=None,
                 category='heatexc', etype=None, mat=None, id=None,
                 Cpo=None, FP=None, FM=None, FBM=None, CBM=None, CTM=None, CGR=None, ACC=None):
        self.mh = mh
        self.mc = mc
        self.cph = cph
        self.cpc = cpc
        self.Thin = Thin
        self.Thout = Thout
        self.Tcin = Tcin
        self.Tcout = Tcout
        self.U = U
        self.F = F
        self.Ns = Ns
        self.area = area
        self.P = P
        self.category = category
        self.etype = etype
        self.mat = mat
        self.id = id
        self.Cpo = Cpo
        self.FP = FP
        self.FM = FM
        self.FBM = FBM
        self.CBM = CBM
        self.CTM = CTM
        self.CGR = CGR
        self.ACC = ACC

    def __repr__(self):
        return '%s: HeatExc(Thin=%sdegC, Thout=%sdegC, Tcin=%sdegC, Tcout=%sdegC, F=%s, Ns=%s, area=%sm^2, etype=%s, mat=%s)' \
               % (self.id, self.Thin, self.Thout, self.Tcin, self.Tcout, round(self.F, 3), self.Ns, round(self.area, 2), self.etype, self.mat)

    def spec(self):
        return '---Design specs for {}:---\n{}\n----------------------------'.format(self.id, vars(self))

    def econ(self):
        return '---Econ report for {}:---\nCpo=${}\nCBM=${}\nCTM=${}\nCGR=${}\nACC=${}\n----------------------------'.format(self.id, round(self.Cpo, 2), round(self.CBM, 2), round(self.CTM, 2), round(self.CGR, 2), round(self.ACC, 2))


def stepwise_leq(a, b, x):

    if len(b) != len(a) + 1:
        raise ValueError('len(b) should be len(a) + 1 !')

    a = (-np.inf,) + a + (np.inf,)

    for i in range(0, len(a)):
        if a[i] <= x < a[i + 1]:
            y = b[i]
            break

    if y == 'error':
        raise ValueError('Input out of supported range!')

    return y


def stepwise_req(a, b, x):

    if len(b) != len(a) + 1:
        raise ValueError('len(b) should be len(a) + 1 !')

    a = (-np.inf,) + a + (np.inf,)

    for i in range(0, len(a)):
        if a[i] < x <= a[i + 1]:
            y = b[i]
            break

    if y == 'error':
        raise ValueError('Input out of supported range!')

    return y


def designP(Po: float) -> float:

    """
    Calculate design pressure for pressure and vacuum vessels
    :param Po: most deviated operating pressure (psig)
    :return: Pd design pressure (psig)
    """

    expo = np.exp(0.60608 + 0.91615 * np.log(Po) + 0.0015655 * pow((np.log(Po)), 2.))
    a = (0., 5., 10., 10.e3)
    b = (Po, 10., max(10., expo), expo, 'error')
    Pd = stepwise_leq(a, b, Po)

    return Pd


def designT(To: float, heuristic: str='towler') -> float:

    """
    Calculate design temperature for pressure and vacuum vessels
    :param To: most deviated operating temperature (degF)
    :param heuristic: either 'Towler' or 'Turton' (optional)
    :return: Td design temperature (degF)
    """

    if 'towler' in str.lower(heuristic):
        if To < Troom:
            Td = To - 25
        else:
            Td = To + 50
    elif 'turton' in str.lower(heuristic):
        if -22 <= To <= 644:
            Td = To + 45
        else:
            raise ValueError('To temperature input out of supported range using Turton heuristic!')
    else:
        raise ValueError('Heuristic not supported! Please check heuristic input!')
    return Td


def maxstress(Td: float, MOC: str='387B') -> (float, float):

    """
    Calculate maximum allowable stress for pressure vessel material
    :param Td: design temperature (degF)
    :param MOC: user-specified material of construction (optional input)
    :return: Smax maximum allowable stress for pressure vessel MOC (psi)
    :return: MOC prescribed material of construction which is in stainless steel family (string, optional).
    If MOC is not user-specified, the returned MOC will be a default value (SA-285C or SA-387B).
    """

    if '317L' in str.upper(MOC):

        MOC = '317L'
        a = (-20., 68., 200., 400., 600., 800., 1000., 1200., 1400., 1600.)
        b = ('error', 25286., 22957., 20957., 19400., 17633., 16733., 15767., 12857., 8300., 'error')

    elif '316Ti' in str.upper(MOC):

        MOC = '316Ti'
        a = (-22., 302., 392., 482., 572., 617., 662., 707., 752., 797., 842., 887.,
             932., 977., 1022., 1067., 1112.)
        b = ('error', 20015., 19435., 18130., 16969., 16824., 16534., 16244., 16099., 15954., 15809., 15664., 15519.,
             15374., 15229., 14475., 11647., 'error')

    elif '316L' in str.upper(MOC):

        MOC = '316L'
        a = (-22., 302., 392., 482., 572., 617., 662., 707., 752., 797., 842., 887.)
        b = ('error', 16679., 15809., 14939., 14214., 13880., 13648., 13460., 13184., 12908., 12734., 12560., 'error')

    elif '304' in str.upper(MOC):

        MOC = '304'
        a = (-22., 149., 212., 257., 302., 392., 482., 572., 617., 662., 707., 752., 797., 842., 887.,
             932., 977., 1022., 1067., 1112.)
        b = ('error', 20015., 19870., 19435., 18855., 18275., 17695., 16824., 16534., 16099., 15809., 15519.,
             15229., 14939., 14649., 14402., 14214., 13532., 11545., 9485., 'error')

    else:  # use default MOCs

        a = (-20., 650., 750., 800., 850., 900.)
        b = ('error', 13750., 15000., 14750., 14200., 13100., 'error')
        c = ('error', '285C', '387B', '387B', '387B', '387B', 'error')
        MOC = stepwise_leq(a, c, Td)

    Smax = stepwise_leq(a, b, Td)

    return Smax, MOC


def elasmod(Td: float, MOC :str='carbon') -> (float, float):

    """
    Calculate modulus of elasticity for vacuum vessel material
    :param Td: design pressure (degF)
    :param MOC: material of construction (only either 'carbon' or 'low-alloy', string)
    :return: EM modulus of elasticity for vacuum vessel MOC (psi)
    :return: MOC returns the input MOC for consistency with the equivalent computation for pressure vessels (optional string)
    """

    if 'carbon' in str.lower(MOC):
        MOC = 'carbon'
        a = (-20., 200., 400., 650.)
        b = (30.2e6, 29.5e6, 28.3e6, 26.0e6, 'error')
        EM = stepwise_req(a, b, Td)

    elif 'low' in str.lower(MOC) and 'alloy' in str.lower(MOC):
        MOC = 'low-alloy'
        a = (-20., 200., 400., 650., 700., 800., 900.)
        b = (30.2e6, 29.5e6, 28.6e6, 27.0e6, 26.6e6, 25.7e6, 24.5e6, 'error')
        EM = stepwise_req(a, b, Td)

    else:
        raise ValueError('Specified MOC not found! Please check MOC input!')

    return EM, MOC


def wallthk(Pd: float, Di: float, Smax: float) -> (float, float):

    """
    Calculate cylindrical shell wall thickness for pressure vessels, including minimum thickness check for structural rigidity
    :param Pd: design pressure (psig)
    :param Di: internal diameter (in)
    :param Smax: maximum allowable stress (psi)
    :return: tp: cylindrical shell wall thickness for pressure
    :return: E: fractional weld efficiency used (string, optional)
    """

    E = 0.85  # first assume 10% X-ray spot check
    tp = Pd * Di / (2*Smax*E - 1.2*Pd)

    if tp > 1.25:  # tp not large enough, 100% X-ray check needed
        E = 1
        tp = Pd * Di / (2*Smax*E - 1.2*Pd)

    # Check if minimum wall thickness to provide rigidity satisfied

    a = (48., 72., 96., 120., 144.)
    b = (max(1/4, tp), max(5/16, tp), max(3/8, tp), max(7/16, tp), max(1/2, tp), 'error')
    tp = stepwise_leq(a, b, tp)

    return tp, E


def wallthkvac(Pd: float, Do: float, Di: float, L: float, EM: float) -> (float, float, float):

    """
    Calculate cylindrical shell wall thickness for vacuum vessels
    :param Pd: design pressure (psig)
    :param Do: external diameter (in)
    :param Di: internal diameter (in)
    :param L: vessel length (in)
    :param EM: modulus of elasticity (psi)
    :return: tp: cylindrical shell wall thickness for vacuum vessels (in)
    :return: tE: necessary thickness for vacuum vessels (in, optional)
    :return: tEC: correction factor for vacuum vessels (in, optional)
    """

    tE = pow(1.3 * Do * (Pd * L / (EM * Do)), 0.4)
    if tE / Do > 0.05:
        warnings.warn('tE is > 0.05*Do, which is' +
                      ' outside the validity range for tE computation!' +
                      ' Nevertheless carrying on with calculation - beware!')

    tEC = L * (0.18*Di - 2.2)*1e-5 - 0.19
    tp = tE + tEC

    return tp, tE, tEC


def shellthkhorz(tp: float) -> float:

    """
    Calculate shell thickness for horizontal vessels
    :param tp: wall thickness (in)
    :return: ts: shell thickness with corrosion allowance for horizontal vessels (in)
    """

    ts = tp + tc

    return ts


def windalw(Do: float, L: float, Smax: float) -> float:

    """
    Calculate wind/earthquake allowance for vertical vessels
    Caution: Using WINDALW requires an assumed value of Do
    which is dependent on tw. If Do is unknown, use
    SHELLTHKVERT directly instead which internally calls WINDALW.
    :param Do: external diameter (in)
    :param L: internal tangent-to-tangent height (in)
    :param Smax: maximum allowable stress (psi)
    :return: tw: wind/earthquake allowance for vertical vessels (in)
    """

    tw = 0.22 * (Do + 18.) * (L ** 2.) / (Smax * Do ** 2.)

    return tw


def shellthkvert(tp: float, Di: float, L: float, Smax: float) -> (float, float, float):

    """
    Calculate shell thickness for vertical vessels
    :param tp: wall thickness (in)
    :param Di: internal diameter (in)
    :param L: internal tangent-to-tangent height (in)
    :param Smax: maximum allowable stress of MOC (psi)
    :return: ts: shell thickness with wind allowance after adding corrosion allowance for vertical vessels (in)
    :return: tv: shell thickness with wind allowance before adding corrosion allowance for vertical vessels (in)
    :return: tw: wind allowance (in, optional)
    """

    ts0 = 2. * tp  # dummy initialisation
    Do = Di + 2. * ts0
    tw = windalw(Do, L, Smax)
    tv = (tp + (tp+tw)) / 2.
    ts1 = tv + tc  # add corrosion allowance

    reltol = 1e-9
    i = 0
    while abs(ts1 - ts0) / ts0 > reltol and i < 1e3:
        ts0 = ts1
        i += 1
        Do = Di + 2. * ts0
        tw = windalw(Do, L, Smax)
        tv = (tp + (tp+tw)) / 2.
        ts1 = tv + tc  # add corrosion allowance

    if i == 1e3:
        warnings.warn('Vertical vessel thickness failed to converge!' +
                      ' Nevertheless carrying on with calculation - beware!')

    ts = ts1

    return ts, tv, tw


def ceilplatethk(ts: float) -> float:

    """
    Round up metal plate thickness to nearest increment
    :param ts: shell wall thickness before before rounding to nearest increment (in)
    :return: tsfinal: final shell wall thickness after rounding to nearest increment (in)
    """

    if ts > 3.:
        warnings.warn('Calculated ts not in supported range.' +
                      ' Assuming metal plate thickness in increments' +
                      ' of 1/4 inches above 3 inches!')

    a = (3/16, 1/2, 2., 3.)
    b = ('error', 1/16, 1/8, 1/4, 1/4)
    acc = stepwise_req(a, b, ts)

    tsfinal = np.ceil(ts / acc) * acc

    return tsfinal


def vesselweight(Di: float, tsfinal: float, L: float, rho: float=rhosteel) -> float:

    """
    Calculate final weight of vessel of the vessel with the shell and two 2:1 elliptical heads
    :param Di: internal diameter (in)
    :param tsfinal: shell thickness with corrosion allowance, rounded to nearest thickness increment for metal plates (in)
    :param L: internal tangent-to-tangent length/height (in)
    :param rho: density of material of construction (MOC) (lb/in^3)
    :return: W: weight of vessel (lb)
    """

    W = np.pi * (Di+tsfinal) * (L+0.8*Di) * tsfinal * rho

    return W


def vesselvol(Do: float, L: float) -> float:

    """
    Calculate final external volume of vessel with the shell and two 2:1 elliptical heads
    :param Do: external diameter (in)
    :param L: internal tangent-to-tangent length/height (in)
    :return: volume of vessel (in^3)
    """

    Vcyl = np.pi * pow(Do, 2.) / 4. * L
    H = Do / 4.
    Vheads = 4. / 3. * np.pi * H * pow((Do / 2.), 2.)
    V = Vcyl + Vheads

    return V


def designhorzpres(Di: float, L: float, Po: float=Patm, To: float=Troom, rho: float=rhosteel, MOC: str='387B',
                   mat: str='SS', id: str='UnnamedVessel') -> MechDesign():

    """
    The main function to be called for designing horizontal pressure vessels
    Example implementation:
    md = designhorzpres(Di=78, L=480, Po=470, To=850)
    :param Di: internal diameter (in)
    :param L: tangent-to-tangent horizontal length (in)
    :param Po: most deviated operating pressure from ambient pressure (psig)
    :param To: most deviated operating temperature from ambient temperature (degF)
    :param rho: density of material of construction (lb/in^3, optional)
    :param MOC: material of construction (string e.g. '387B' [default] or '317L' etc., optional)
    :param mat: category of material of construction (string e.g. 'CS' or 'SS' [default] etc., optional)
    :param id: id/name of equipment (string, e.g. V100, optional)
    :return: md: MechDesign object (optional) consisting of:
    Pd = design pressure (psig)
    Td = design pressure (degF)
    MOC = material of construction to use (string)
    Smax = maximum allowable stress of MOC used (psi)
    E = weld efficiency to use (dimensionless)
    tp = wall thickness (in)
    tc = corrosion allowance used (= 1/8 in)
    ts = tp with tc (in)
    tsfinal = ts rounded up to next increment in metal plate thickness (in)
    Do = external diameter (in)
    W = total vessel weight (lb)
    V = total vessel external volume (in^3)
    Vi = total vessel internal volume (in^3)
    """

    md = MechDesign()
    md.Po = Po
    md.To = To
    md.Di = Di
    md.L = L
    md.rho = rho

    md.Pd = designP(Po)
    md.Td = designT(To)
    md.Smax, md.MOC = maxstress(md.Td, MOC)
    md.tp, md.E = wallthk(md.Pd, Di, md.Smax)
    md.tc = tc
    md.ts = shellthkhorz(md.tp)
    md.tsfinal = ceilplatethk(md.ts)
    md.Do = Di + 2. * md.tsfinal
    md.W = vesselweight(Di, md.tsfinal, L, rho)
    md.V = vesselvol(md.Do, L)
    md.Vi = np.pi * (Di ** 2.) / 4. * L

    md.id = id
    md.category = 'vessel'
    md.etype = 'horizontal'
    if mat is not None:
        md.mat = mat if mat in capex.matfaclib['vessel']['horizontal'].keys() else None
    else:
        md.mat = 'SS' if MOC in ['317L', '316L', '304'] else 'CS' if MOC in ['285C', '387B', 'low-alloy', 'carbon'] else 'Ti' if MOC in ['316Ti'] else None
    if md.mat is None:
        md.mat = 'SS'
        warnings.warn('Type of MOC (mat variable) cannot be identified and is assumed to be stainless steel! ' +
                      'You can specify a mat input (mat=' + capex.matfaclib['vessel']['horizontal'].keys())

    return md


def designvertpres(Di: float, L: float, Po: float=Patm, To: float=Troom, rho: float=rhosteel,
                   MOC: str='387B', mat: str='SS', id: str='UnnamedVessel') -> MechDesign():

    """
    The main function to be called for designing vertical pressure vessels
    Example implementation:
    md = designhorzpres(Di=120, L=2544, Po=95.5, To=150)
    :param Di: internal diameter (in)
    :param L: tangent-to-tangent horizontal height (in)
    :param Po: most deviated operating pressure from ambient pressure (psig)
    :param To: most deviated operating temperature from ambient temperature (degF)
    :param rho: density of material of construction (lb/in^3, optional)
    :param MOC: material of construction (string e.g. '387B' [default] or '317L' etc., optional)
    :param mat: category of material of construction (string e.g. 'CS' or 'SS' [default] etc., optional)
    :param id: id/name of equipment (string, e.g. V100, optional)
    :return: md: MechDesign object (optional) consisting of:
    Pd = design pressure (psig)
    Td = design pressure (degF)
    MOC = material of construction to use (string)
    Smax = maximum allowable stress of MOC used (psi)
    E = weld efficiency to use (dimensionless)
    tp = wall thickness (in)
    tc = corrosion allowance used (= 1/8 in)
    tw = wind/earthquake allowance for vertical vessels (in)
    tv = tp with tw without tc (in)
    ts = tp with tc (in)
    tsfinal = ts rounded up to next increment in metal plate thickness (in)
    Do = external diameter (in)
    W = total vessel weight (lb)
    V = total vessel external volume (in^3)
    Vi = total vessel internal volume (in^3)
    """

    md = MechDesign()
    md.Po = Po
    md.To = To
    md.Di = Di
    md.L = L
    md.rho = rho

    md.Pd = designP(Po)
    md.Td = designT(To)
    md.Smax, md.MOC = maxstress(md.Td, MOC)
    md.tp, md.E = wallthk(md.Pd, Di, md.Smax)
    md.tc = tc
    md.ts, md.tv, md.tw = shellthkvert(md.tp, Di, L, md.Smax)
    md.tsfinal = ceilplatethk(md.ts)
    md.Do = Di + 2. * md.tsfinal
    md.W = vesselweight(Di, md.tsfinal, L, rho)
    md.V = vesselvol(md.Do, L)
    md.Vi = np.pi * (Di ** 2.) / 4. * L

    md.id = id
    md.category = 'vessel'
    md.etype = 'vertical'
    if mat is not None:
        md.mat = mat if mat in capex.matfaclib['vessel']['horizontal'].keys() else None
    else:
        md.mat = 'SS' if MOC in ['317L', '316L', '304'] else 'CS' if MOC in ['285C', '387B', 'low-alloy', 'carbon'] else 'Ti' if MOC in ['316Ti'] else None
    if md.mat is None:
        md.mat = 'SS'
        warnings.warn('Type of MOC (mat variable) cannot be identified and is assumed to be stainless steel! ' +
                      'You can specify a mat input (mat=' + capex.matfaclib['vessel']['vertical'].keys())

    return md


def designvac(Di: float, L: float, Po: float=Patm, To: float=Troom, rho: float=rhosteel,
              MOC: str='carbon', etype: str=None, mat: str='CS', id: str='UnnamedVessel') -> MechDesign():

    """
    The main function to be called for designing vacuum vessels
    Example implementation:
    md = dsg.designvac(Di=168., L=1080., Po=7.977, To=257.)
    :param Di: internal diameter (in)
    :param L: tangent-to-tangent horizontal length/height (in)
    :param Po: most deviated operating pressure from ambient pressure (psig)
    :param To: most deviated operating temperature from ambient temperature (degF)
    :param rho: density of material of construction (lb/in^3, optional)
    :param MOC: material of construction (string e.g. 'carbon' as default, optional)
    :param etype: type of vessel ('horizontal' or 'vertical')
    :param mat: category of material of construction (string e.g. 'CS' [default] or 'SS' etc., optional)
    :param id: id/name of equipment (string, e.g. V100, optional)
    :return: md: MechDesign object (optional) consisting of:
    Pd = design pressure (psig)
    Td = design pressure (degF)
    MOC = material of construction to use (string)
    EM = modulus of elasticity of MOC used (psi)
    tE = vacuum wall thickness (in)
    tEC = vacuum wall correction factor (in)
    tp = tE with tEC (in)
    tc = corrosion allowance used (= 1/8 in)
    ts = tp with tc (in)
    tsfinal = ts rounded up to next increment in metal plate thickness (in)
    Do = external diameter (in)
    W = total vessel weight (lb)
    V = total vessel external volume (in^3)
    Vi = total vessel internal volume (in^3)
    """

    md = MechDesign()
    md.Po = Po
    md.To = To
    md.Di = Di
    md.L = L
    md.rho = rho

    md.Pd = Patm - Po
    md.Td = designT(To)
    if -20 < md.Td <= 650:
        [md.EM, md.MOC] = elasmod(md.Td, 'carbon')
    elif 650 < md.Td < 900:
        [md.EM, md.MOC] = elasmod(md.Td, 'low-alloy')
    else:
        raise ValueError('Td out of supported range for both carbon and low-alloy steel!')

    ts0 = 1  # dummy initialisation
    md.Do = Di + 2. * ts0
    md.tp, md.tE, md.tEC = wallthkvac(md.Pd, md.Do, Di, L, md.EM)
    md.tc = tc
    ts1 = shellthkhorz(md.tp)  # horz/vert orientation does not matter for vacuum

    reltol = 1.e-9
    i = 0
    while abs(ts1 - ts0) / ts0 > reltol and i < 1e3:
        ts0 = ts1
        i += 1
        md.Do = Di + 2. * ts0
        md.tp, md.tE, md.tEC = wallthkvac(md.Pd, md.Do, Di, L, md.EM)
        md.tc = tc
        ts1 = shellthkhorz(md.tp)

    if i == 1e3:
        warnings.warn('Vacuum vessel thickness failed to converge! ' +
                      'Nevertheless carrying on with calculation - beware!')

    md.ts = ts1
    md.tsfinal = ceilplatethk(md.ts)
    md.Do = Di + 2. * md.tsfinal
    md.W = vesselweight(Di, md.tsfinal, L, rho)
    md.V = vesselvol(md.Do, L)
    md.Vi = np.pi * (Di ** 2.) / 4. * L

    md.id = id
    md.category = 'vessel'
    if etype is None:
        etype = 'vertical'
        warnings.warn('Assuming vacuum vessel is vertical! ' +
                      'You can specify a etype input (etype=' + str(capex.eqptcostlib['vessel'].keys()))
    md.etype = str.lower(etype) if (str.lower(etype) in capex.eqptcostlib['vessel'].keys()) else None

    if mat is not None:
        md.mat = mat if mat in capex.matfaclib['vessel']['horizontal'].keys() else None
    else:
        md.mat = 'SS' if MOC in ['317L', '316L', '304'] else 'CS' if MOC in ['285C', '387B', 'low-alloy', 'carbon'] else 'Ti' if MOC in ['316Ti'] else None
    if md.mat is None:
        md.mat = 'CS'
        warnings.warn('Type of MOC (mat variable) cannot be identified and is assumed to be carbon steel! ' +
                      'You can specify a mat input (mat=' + str(capex.matfaclib['vessel']['vertical'].keys()))

    return md


def sizecompressor(m: float, P1: float, P2: float, T1: float, cp: float, cv: float, Z: float=1.,
                   etype: str=None, mat: str=None, id: str='UnnamedCompressor') -> (float, float, float, Compressor()):

    """
    Conducts compressor sizing by determining required compressor power
    based on its flow rate and inlet/outlet pressures
    Example implementation:
    comppower, compeff, T2, compressor = dsg.sizecompressor(m=1e5, P1=100, P2=300, T1=323.15, cp=1.02, cv=0.72, Z=0.99)
    :param m: mass flow rate through compressor (kg/h)
    :param P1: gas inlet pressure (bar)
    :param P2: gas inlet pressure (bar)
    :param T1: gas inlet temperature (K)
    :param cp: constant-pressure heat capacity of gas
    :param cv: constant-volume heat capacity of gas
    :param Z: gas compressibility factor (optional, default = 1)
    :param etype: type of equipment (string, e.g. 'centrifugal' [default] or 'axial' etc., optional)
    :param mat: category of material of construction (string e.g. 'CS' or 'SS' [default] etc., optional)
    :param id: id/name of equipment (string, e.g. K100, optional)
    :return: comppower: required compressor power (kW)
    :return: compeff: compressor efficiency (optional, dimensionless)
    :return: gas outlet temperature (optional, K)
    :return: compressor: Compressor object (optional)
    """

    if P2/P1 > 4.:
        warnings.warn('Compression ratio > 4 is too large -' +
                      ' check that outlet temperature is not too high!' +
                      ' Nevertheless continuing calculation...')
    elif P2/P1 < 1.:
        raise ValueError('Outlet pressure smaller than inlet pressure!')

    m = m / 3600.
    k = cp / cv
    a = (k - 1) / k
    power = (m * Z * R * T1) * (pow((P2 / P1), a) - 1.) / a  # useful power
    power /= 1000.  # convert Pa to kPa

    compeff = np.interp(P2 / P1, [1., 1.5, 2., 3., 6., 10.],
                        [0.65-np.spacing(1), 0.65, 0.75, 0.8, 0.85, 0.85+np.spacing(1)])

    comppower = power / compeff

    T2 = T1 * pow(P2 / P1, a)

    if T2 > 273.15 + 200:
        warnings.warn('Gas outlet temperature too high! ' +
                      'Consider reducing compression ratio P2/P1! ' +
                      'Nevertheless continuing calculation...')

    compressor = Compressor()
    compressor.m = m
    compressor.P1 = P1
    compressor.P2 = P2
    compressor.T1 = T1
    compressor.T2 = T2
    compressor.cp = cp
    compressor.cv = cv
    compressor.Z = Z
    compressor.compeff = compeff
    compressor.comppower = comppower

    compressor.id = id
    compressor.category = 'compressor'
    if etype is None:
        etype = 'centrifugal'
        warnings.warn('Assuming compressor is centrifugal! ' +
                      'You can specify a etype input (etype=' + str(capex.eqptcostlib['compressor'].keys()))
    compressor.etype = str.lower(etype) if (str.lower(etype) in capex.eqptcostlib['compressor'].keys()) else None

    if mat is None:
        mat = 'SS'
        warnings.warn('Assuming compressor material is stainless steel! ' +
                      'You can specify a mat input (mat=' + str(capex.matfaclib['compressor']['centrifugal'].keys()))
    compressor.mat = mat if (mat in capex.matfaclib['compressor']['centrifugal'].keys()) else None

    return comppower, compeff, T2, compressor


def sizepump(Q: float, dP: float=None, P1: float=None, P2: float=None, rho: float=1000., pumpeff: float=None,
             etype: str=None, mat: str=None, id: str='UnnamedPump') -> (float, float, Pump()):

    """
    Conducts pump sizing by determining required pump power
    based on its flow rate and pressure differential (discharge - suction pressure)
    Example implementation:
    pumppower, pumpeff, pump = dsg.sizepump(Q=35, dP=500)
    :param Q: volumetric flow rate through pump (m^3/h)
    :param P1: suction/inlet pressure (kPa)
    :param P2: discharge/outlet pressure (kPa)
    :param rho: stream density (kg/m^3) (optional, default = 1000)
    :param pumpeff: pump efficiency (optional, default = 0.75)
    :param etype: type of equipment (string, e.g. 'centrifugal' [default] or 'reciprocating' etc., optional)
    :param mat: category of material of construction (string e.g. 'CS' or 'SS' [default] etc., optional)
    :param id: id/name of equipment (string, e.g. P100, optional)
    :return: pumppower = required pump power (kW)
    :return: pumpeff  = pump efficiency (optional output - if it is not specified in input, it will be calculated)
    :return: pump: Pump object (optional)
    """

    if dP is None:
        if P1 is None or P2 is None:
            raise ValueError('dP, P1 or P2 not specified!')
        elif P2 > P1:
            dP = P2 - P1
        else:
            raise ValueError('Outlet pressure lower than inlet pressure!')

    power = (Q / 3600.) * dP  # useful power in kW

    H = dP / (rho * g)  # required head in m
    H_ft = H * 3.281  # required head in ft
    Q_gpm = Q * 4.403  # flowrate in gal/min (gpm)

    if pumpeff is None:
        if 50 <= H_ft <= 300 and 100 <= Q_gpm <= 1000:
            a = np.array([80., -0.2855, 3.78e-4, -2.38e-7, 5.39e-4, -6.39e-7, 4.e-10])
            b = np.array([1, H_ft, H_ft*Q_gpm, H_ft*pow(Q_gpm, 2),
                           pow(H_ft, 2), pow(H_ft, 2)*Q_gpm, pow(H_ft, 2)*pow(Q_gpm, 2)])
            pumpeff = (a @ b.T) / 100.
        elif 0 <= power <= 300:
            # Maximum useful power for centrifugal pumps = 300 kW
            pumpeff = np.interp(power, [0., 2., 5., 10., 30., 55., 300.],
                                [0.55-np.spacing(1), 0.55, 0.6, 0.65, 0.7, 0.75, 0.75+np.spacing(1)])
        else:
            pumpeff = 0.75

    pumppower = power / pumpeff

    pump = Pump()
    pump.Q = Q
    pump.P1 = P1
    pump.P2 = P2
    pump.dP = dP
    pump.rho = rho
    pump.pumpeff = pumpeff
    pump.pumppower = pumppower

    pump.id = id
    pump.category = 'pump'
    if etype is None:
        etype = 'centrifugal'
        warnings.warn('Assuming pump is centrifugal! ' +
                      'You can specify a type einput (etype=' + str(capex.eqptcostlib['pump'].keys()))
    pump.etype = str.lower(etype) if (str.lower(etype) in capex.eqptcostlib['pump'].keys()) else None

    if mat is None:
        mat = 'SS'
        warnings.warn('Assuming pump material is stainless steel! ' +
                      'You can specify a mat input (mat=' + str(capex.matfaclib['pump']['centrifugal'].keys()))
    pump.mat = mat if (mat in capex.matfaclib['pump']['centrifugal'].keys()) else None

    return pumppower, pumpeff, pump


def sizeHE_heater(mc: float, cpc: float, Tcin: float, Tcout: float, Thin: float, Thout: float,
                  U: float, F: float=None, Ns: int=1, etype: str=None, mat: str=None,
                  P: float=None, id: str='UnnamedHX') -> (float, float, HeatExc()):

    """
    Conducts shell-and-tube heat exchanger sizing (counterflow arrangement), where cold process stream is heated,
    by determining required heat exchange area
    Example implementation:
    area, F, HX = dsg.sizeHE_heater(mc=31715, cpc=3246, Tcin=89, Tcout=101, Thin=160, Thout=156, U=850)
    :param mc: cold stream mass flow rate (kg/h)
    :param cpc: heat capacity of cold stream % J/(kg.K)
    :param Tcin: cold stream inlet temperature (degC)
    :param Tcout: cold stream outlet temperature (degC)
    :param Thin: hot stream inlet temperature (degC)
    :param Thout: hot stream outlet temperature (degC)
    :param U: heat transfer coefficient (W/(m^2.degC))
    :param F: user-specified correction factor (if not specified, F will be calculated)
    :param Ns: number of shell passes (default = 1)
    :param etype: type of equipment (string, e.g. 'utube' [default] or 'doublepipe' etc., optional)
    :param mat: category of material of construction (string e.g. 'CS/CS' or 'SS/CS' [default] etc., optional)
    :param id: id/name of equipment (string, e.g. HX100, optional)
    :param P: operating pressure, for cost calculation purposes only (bar, optional)
    :return: area: required heat exchange area (m^2)
    :return: F: correction factor (optional output - if F is not specified in input, F will be calculated)
    :return: HX: HeatExc object (optional)
    """

    if Thout > Thin:
        raise ValueError('Hot stream outlet cannot be hotter than inlet!')
    elif Tcout < Tcin:
        raise ValueError('Cold stream outlet cannot be colder than inlet!')
    elif Tcout > Thout:
        warnings.warn('Potential temperature cross - Cold stream outlet is hotter than hot stream inlet! Nevertheless continuing with calculations...')
    elif Thout - Tcin < Ta:
        raise ValueError('Minimum temperature not fulfilled for hot outlet / cold inlet side!')
    elif Thin - Tcout < Ta:
        raise ValueError('Minimum temperature not fulfilled for hot inlet / cold outlet side!')

    mc /= 3600.  # convert kg/h to kg/s

    Q = mc * cpc * (Tcout - Tcin)  # calculate heat transfer rate in W

    if (Thout - Tcin) == (Thin - Tcout):
        LMTD = Thout - Tcin
    else:
        LMTD = ((Thin-Tcout) - (Thout-Tcin)) / np.log((Thin-Tcout) / (Thout-Tcin))

    if F is None:
        r = (Thin - Thout) / (Tcout - Tcin)
        p = (Tcout - Tcin) / (Thin - Tcin)
        if r == 1:
            w = (Ns - Ns * p) / (Ns - Ns * p + p)
            F = (np.sqrt(2) * (1 - w) / w) / \
                (np.log(w / (1 - w) + 1 / np.sqrt(2)) / np.log(w / (1 - w) - 1 / np.sqrt(2)))
        else:
            w = pow((1 - p * r) / (1 - p), 1 / Ns)
            s = np.sqrt(r ** 2 + 1) / (r - 1)
            F = s * np.log(w) / np.log((1 + w - s + s * w) / (1 + w + s - s * w))

    area = Q / (U * F * LMTD)

    HX = HeatExc()
    HX.mh = None
    HX.mc = mc
    HX.cph = None
    HX.cpc = cpc
    HX.Thin = Thin
    HX.Thout = Thout
    HX.Tcin = Tcin
    HX.Tcout = Tcout
    HX.U = U
    HX.F = F
    HX.Ns = Ns
    HX.area = area
    HX.P = P

    HX.id = id
    HX.category = 'heatexc'
    if etype is None:
        etype = 'utube'
        warnings.warn('Assuming heat exchanger is U-tube! ' +
                      'You can specify a etype input (etype=' + str(capex.eqptcostlib['heatexc'].keys()))
    HX.etype = str.lower(etype) if (str.lower(etype) in capex.eqptcostlib['heatexc'].keys()) else None

    if mat is None:
        mat = 'CS/SS'
        warnings.warn('Assuming hext exchanger material is carbon steel/stainless steel (or vice versa)! ' +
                      'You can specify a mat input (mat=' + str(capex.matfaclib['heatexc']['utube'].keys()))
    HX.mat = mat if (mat in capex.matfaclib['heatexc']['utube'].keys()) else None

    return area, F, HX


def sizeHE_cooler(mh: float, cph: float, Thin: float, Thout: float, Tcin: float, Tcout: float,
                  U: float, F: float=None, Ns: int=1, etype: str=None, mat: str=None,
                  P: float=None, id: str='UnnamedHX') -> (float, float, HeatExc()):

    """
    Conducts shell-and-tube heat exchanger sizing (counterflow arrangement), where hot process stream is cooled,
    by determining required heat exchange area
    Example implementation:
    area, F, HX = dsg.sizeHE_cooler(mh=31715, cph=3246, Thin=89, Thout=60, Tcin=5, Tcout=10, U=850)
    :param mh: hot stream mass flow rate (kg/h)
    :param cph: heat capacity of hot stream % J/(kg.K)
    :param Thin: hot stream inlet temperature (degC)
    :param Thout: hot stream outlet temperature (degC)
    :param Tcin: cold stream inlet temperature (degC)
    :param Tcout: cold stream outlet temperature (degC)
    :param U: heat transfer coefficient (W/(m^2.degC))
    :param F: user-specified correction factor (if not specified, F will be calculated)
    :param Ns: number of shell passes (default = 1)
    :param etype: type of equipment (string, e.g. 'utube' [default] or 'doublepipe' etc., optional)
    :param mat: category of material of construction (string e.g. 'CS/CS' or 'SS/CS' [default] etc., optional)
    :param id: id/name of equipment (string, e.g. HX100, optional)
    :param P: operating pressure, for cost calculation purposes only (bar, optional)
    :return: area: required heat exchange area (m^2)
    :return: F: correction factor (optional output - if F is not specified in input, F will be calculated)
    :return: HX: HeatExc object (optional)
    """

    if Thout > Thin:
        raise ValueError('Hot stream outlet cannot be hotter than inlet!')
    elif Tcout < Tcin:
        raise ValueError('Cold stream outlet cannot be colder than inlet!')
    elif Tcout > Thout:
        warnings.warn('Potential temperature cross - Cold stream outlet is hotter than hot stream inlet! Nevertheless continuing with calculations...')
    elif Thout - Tcin < Ta:
        raise ValueError('Minimum temperature not fulfilled for hot outlet / cold inlet side!')
    elif Thin - Tcout < Ta:
        raise ValueError('Minimum temperature not fulfilled for hot inlet / cold outlet side!')

    mh /= 3600.  # convert kg/h to kg/s

    Q = mh * cph * (Thin - Thout)  # calculate heat transfer rate

    if (Thout - Tcin) == (Thin - Tcout):
        LMTD = Thout - Tcin
    else:
        LMTD = ((Thin - Tcout) - (Thout - Tcin)) / np.log((Thin - Tcout) / (Thout - Tcin))

    if F is None:
        r = (Thin - Thout) / (Tcout - Tcin)
        p = (Tcout - Tcin) / (Thin - Tcin)
        if r == 1:
            w = (Ns - Ns * p) / (Ns - Ns * p + p)
            F = (np.sqrt(2) * (1 - w) / w) / \
                (np.log(w / (1 - w) + 1 / np.sqrt(2)) / np.log(w / (1 - w) - 1 / np.sqrt(2)))
        else:
            w = pow((1 - p * r) / (1 - p), 1 / Ns)
            s = np.sqrt(r ** 2 + 1) / (r - 1)
            F = s * np.log(w) / np.log((1 + w - s + s * w) / (1 + w + s - s * w))

    area = Q / (U * F * LMTD)

    HX = HeatExc()
    HX.mh = mh
    HX.mc = None
    HX.cph = cph
    HX.cpc = None
    HX.Thin = Thin
    HX.Thout = Thout
    HX.Tcin = Tcin
    HX.Tcout = Tcout
    HX.U = U
    HX.F = F
    HX.Ns = Ns
    HX.area = area
    HX.P = P

    HX.id = id
    HX.category = 'heatexc'
    if etype is None:
        etype = 'utube'
        warnings.warn('Assuming heat exchanger is U-tube! ' +
                      'You can specify a etype input (etype=' + str(capex.eqptcostlib['heatexc'].keys()))
    HX.etype = str.lower(etype) if (str.lower(etype) in capex.eqptcostlib['heatexc'].keys()) else None

    if mat is None:
        HX.mat = 'CS/SS'
        warnings.warn('Assuming hext exchanger material is carbon steel/stainless steel (or vice versa)! ' +
                      'You can specify a mat input (mat=' + str(capex.matfaclib['heatexc']['utube'].keys()))
    HX.mat = mat if (mat in capex.matfaclib['heatexc']['utube'].keys()) else None

    return area, F, HX
