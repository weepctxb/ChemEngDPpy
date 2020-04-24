import numpy as np
import warnings

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
                 Do=None, W=None, V=None,
                 EM=None, tE=None, tEC=None):
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
        self.EM = EM
        self.tE = tE
        self.tEC = tEC


class Compressor(object):
    def __init__(self, m=None, P1=Patm, P2=None, T1=Troom, T2=None, cp=None, cv=None, Z=None,
                 compeff=None, comppower=None):
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


class Pump(object):
    def __init__(self, Q=None, P1=None, P2=None, dP=None, rho=1000,
                 pumpeff=0.75, pumppower=None):
        self.Q = Q
        self.P1 = P1
        self.P2 = P2
        self.dP = dP
        self.rho = rho
        self.pumpeff = pumpeff
        self.pumppower = pumppower


class HeatExc(object):
    def __init__(self, mh=None, mc=None, cph=None, cpc=None, Thin=None, Thout=None, Tcin=None, Tcout=None,
                 U=None, F=0.9, Ns=1, A=None):
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
        self.A = A


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


def designP(Po):

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


def designT(To, heuristic='towler'):

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


def maxstress(Td, MOC='387B'):

    """
    Calculate maximum allowable stress for pressure vessel material
    :param Td: design temperature (degF)
    :param MOC: user-specified material of construction (optional input)
    :return: Smax maximum allowable stress for pressure vessel MOC (psi)
    :return: MOC prescribed material of construction which is in stainless steel family (string, optional).
    If MOC is not user-specified, the returned MOC will be a default value (SA-285C or SA-387B).
    """

    if '317L' in str.upper(MOC):

        MOC = 'SA-317L'
        a = (-20., 68., 200., 400., 600., 800., 1000., 1200., 1400., 1600.)
        b = ('error', 25286., 22957., 20957., 19400., 17633., 16733., 15767., 12857., 8300., 'error')
        Smax = stepwise_leq(a, b, Td)

    else:

        a = (-20., 650., 750., 800., 850., 900.)
        b = ('error', 13750., 15000., 14750., 14200., 13100., 'error')
        c = ('error', 'SA-285C', 'SA-387B', 'SA-387B', 'SA-387B', 'SA-387B', 'error')
        Smax = stepwise_leq(a, b, Td)
        MOC = stepwise_leq(a, c, Td)

    return Smax, MOC


def elasmod(Td, MOC='carbon'):

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


def wallthk(Pd, Di, Smax):

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


def wallthkvac(Pd, Do, Di, L, EM):

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


def shellthkhorz(tp):

    """
    Calculate shell thickness for horizontal vessels
    :param tp: wall thickness (in)
    :return: ts: shell thickness with corrosion allowance for horizontal vessels (in)
    """

    ts = tp + tc

    return ts


def windalw(Do, L, Smax):

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

    tw = 0.22 * (Do + 18) * L^2 / (Smax * Do ** 2)

    return tw


def shellthkvert(tp, Di, L, Smax):

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

    ts0 = 2 * tp  # dummy initialisation
    Do = Di + 2 * ts0
    tw = windalw(Do, L, Smax)
    tv = (tp + (tp+tw)) / 2
    ts1 = tv + tc  # add corrosion allowance

    reltol = 1e-9
    i = 0
    while abs(ts1 - ts0) / ts0 > reltol and i < 1e3:
        ts0 = ts1
        i += 1
        Do = Di + 2 * ts0
        tw = windalw(Do, L, Smax)
        tv = (tp + (tp+tw)) / 2
        ts1 = tv + tc  # add corrosion allowance

    if i == 1e3:
        warnings.warn('Vertical vessel thickness failed to converge!' +
                      ' Nevertheless carrying on with calculation - beware!')

    ts = ts1

    return ts, tv, tw


def ceilplatethk(ts):

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


def vesselweight(Di, tsfinal, L, rho=rhosteel):

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


def vesselvol(Do, L):

    """
    Calculate final volume of vessel with the shell and two 2:1 elliptical heads
    :param Do: external diameter (in)
    :param L: internal tangent-to-tangent length/height (in)
    :return: volume of vessel (in^3)
    """

    Vcyl = np.pi * pow(Do, 2.) / 4. * L
    H = Do / 4.
    Vheads = 4. / 3. * np.pi * H * pow((Do / 2.), 2.)
    V = Vcyl + Vheads

    return V


def designhorzpres(Di, L, Po=Patm, To=Troom, rho=rhosteel, MOC='SA-387B'):

    """
    The main function to be called for designing horizontal pressure vessels
    Example implementation:
    mechdesign1 = designhorzpres(Di=78, L=480, Po=470, To=850)
    :param Di: internal diameter (in)
    :param L: tangent-to-tangent horizontal length (in)
    :param Po: most deviated operating pressure from ambient pressure (psig)
    :param To: most deviated operating temperature from ambient temperature (degF)
    :param rho: density of material of construction (lb/in^3, optional)
    :param MOC: material of construction (string, optional)
    :return: md: MechDesign class consisting of:
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
    V = total vessel volume (in^3)
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
    md.Do = Di + 2 * md.tsfinal
    md.W = vesselweight(Di, md.tsfinal, L, rho)
    md.V = vesselvol(md.Do, L)

    return md


def designvertpres(Di, L, Po=Patm, To=Troom, rho=rhosteel, MOC='SA-387B'):

    """
    The main function to be called for designing vertical pressure vessels
    Example implementation:
    mechdesign2 = designhorzpres(Di=120, L=2544, Po=95.5, To=150)
    :param Di: internal diameter (in)
    :param L: tangent-to-tangent horizontal height (in)
    :param Po: most deviated operating pressure from ambient pressure (psig)
    :param To: most deviated operating temperature from ambient temperature (degF)
    :param rho: density of material of construction (lb/in^3, optional)
    :param MOC: material of construction (string, optional)
    :return: md: MechDesign class consisting of:
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
    V = total vessel volume (in^3)
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
    md.Do = Di + 2 * md.tsfinal
    md.W = vesselweight(Di, md.tsfinal, L, rho)
    md.V = vesselvol(md.Do, L)

    return md


def designvac(Di, L, Po=Patm, To=Troom, rho=rhosteel, MOC='carbon'):

    """
    The main function to be called for designing vacuum vessels
    :param Di: internal diameter (in)
    :param L: tangent-to-tangent horizontal length/height (in)
    :param Po: most deviated operating pressure from ambient pressure (psig)
    :param To: most deviated operating temperature from ambient temperature (degF)
    :param rho: density of material of construction (lb/in^3, optional)
    :param MOC: material of construction (string, optional)
    :return: md: MechDesign class consisting of:
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
    V = total vessel volume (in^3)
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
    md.Do = Di + 2 * ts0
    md.tp, md.tE, md.tEC = wallthkvac(md.Pd, md.Do, Di, L, md.EM)
    md.tc = tc
    ts1 = shellthkhorz(md.tp)  # horz/vert orientation does not matter for vacuum

    reltol = 1.e-9
    i = 0
    while abs(ts1 - ts0) / ts0 > reltol and i < 1e3:
        ts0 = ts1
        i += 1
        md.Do = Di + 2 * ts0
        md.tp, md.tE, md.tEC = wallthkvac(md.Pd, md.Do, Di, L, md.EM)
        md.tc = tc
        ts1 = shellthkhorz(md.tp)

    if i == 1e3:
        warnings.warn('Vacuum vessel thickness failed to converge!' +
                      ' Nevertheless carrying on with calculation - beware!')

    md.ts = ts1
    md.tsfinal = ceilplatethk(md.ts)
    md.Do = Di + 2 * md.tsfinal
    md.W = vesselweight(Di, md.tsfinal, L, rho)
    md.V = vesselvol(md.Do, L)

    return


def sizecompressor(m, P1, P2, T1, cp, cv, Z=1.):

    """
    Conducts compressor sizing by determining required compressor power
    based on its flow rate and inlet/outlet pressures
    Example implementation:
    [comppower, compeff, T2] = dsg.sizecompressor(1e5,1,4,323.15,1.5,1.4,0.96)
    :param m: mass flow rate through compressor (kg/h)
    :param P1: gas inlet pressure (any pressure units)
    :param P2: gas inlet pressure (same pressure unit as P1)
    :param T1: gas inlet temperature (K)
    :param cp: constant-pressure heat capacity of gas
    :param cv: constant-volume heat capacity of gas
    :param Z: gas compressibility factor (optional, default = 1)
    :return: comppower: required compressor power (kW)
    :return: compeff: compressor efficiency (optional, dimensionless)
    :return: gas outlet temperature (optional, K)
    """

    if P2/P1 > 4.:
        warnings.warn('Compression ratio > 4 is too large -' +
                      ' check that outlet temperature is not too high!' +
                      ' Nevertheless continuing calculation...')
    elif P2/P1 < 1.:
        raise ValueError('Outlet pressure smaller than inlet pressure!')

    m = m / 3600
    k = cp / cv
    a = (k - 1) / k
    power = (m * Z * R * T1) * (pow((P2 / P1), a) - 1.) / a  # useful power
    power /= 1000  # convert Pa to kPa

    compeff = np.interp(P2 / P1, [1., 1.5, 2., 3., 6., 10.] ,
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

    return comppower, compeff, T2, compressor


def sizepump(Q, dP=None, P1=None, P2=None, rho=1000, pumpeff=None):

    """
    Conducts pump sizing by determining required pump power
    based on its flow rate and pressure differential (discharge - suction pressure)
    :param Q: volumetric flow rate through pump (m^3/h)
    :param dP: pressure differential (kPa)
    :param P1: suction/inlet pressure (kPa)
    :param P2: discharge/outlet pressure (kPa)
    :param rho: stream density (kg/m^3) (optional, default = 1000)
    :param pumpeff: pump efficiency (optional, default = 0.75)
    :return: pumppower = required pump power (kW)
    :return: pumpeff  = pump efficiency (optional output - if it is not specified in input, it will be calculated)
    """

    if dP is None:
        if P1 is None or P2 is None:
            raise ValueError('dP, P1 or P2 not specified!')
        elif P2 > P1:
            dP = P2 - P1
        else:
            raise ValueError('Outlet pressure lower than inlet pressure!')

    power = (Q / 3600) * dP  # useful power in kW

    H = dP / (rho * g)  # required head in m
    H_ft = H * 3.281  # required head in ft
    Q_gpm = Q * 4.403  # flowrate in gal/min (gpm)

    if pumpeff is None:
        if 50 <= H_ft <= 300 and 100 <= Q_gpm <= 1000:
            a = np.matrix([80., -0.2855, 3.78e-4, -2.38e-7, 5.39e-4, -6.39e-7, 4.e-10])
            b = np.matrix([1, H_ft, H_ft*Q_gpm, H_ft*pow(Q_gpm, 2),
                           pow(H_ft, 2), pow(H_ft, 2)*Q_gpm, pow(H_ft, 2)*pow(Q_gpm, 2)])
            pumpeff = (a @ b.T) / 100
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

    return pumppower, pumpeff, pump


def sizeHE_heater(mc, cpc, Tcin, Tcout, Thin, Thout, U, F=None, Ns=1):

    """
    Conducts shell-and-tube heat exchanger sizing (counterflow arrangement), where cold process stream is heated,
    by determining required heat exchange area
    Example implementation:
    A, F = dsg.sizeHE_heater(31715,3246,89,101,160,156,850)
    :param mc: cold stream mass flow rate (kg/h)
    :param cpc: heat capacity of cold stream % J/(kg.K)
    :param Tcin: cold stream inlet temperature (degC)
    :param Tcout: cold stream outlet temperature (degC)
    :param cph: hot stream inlet temperature (degC)
    :param Thin: hot stream inlet temperature (degC)
    :param Thout: hot stream outlet temperature (degC)
    :param U: heat transfer coefficient (W/(m^2.degC))
    :param F: user-specified correction factor (if not specified, F will be calculated)
    :param Ns: number of shell passes (default = 1)
    :return: A: required heat exchange area (m^2)
    :return: F: correction factor (optional output - if F is not specified in input, F will be calculated)
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

    mc /= 3600  # convert kg/h to kg/s

    Q = mc * cpc * (Tcout - Tcin)  # calculate heat transfer rate

    if (Thout - Tcin) == (Thin - Tcout):
        LMTD = Thout - Tcin
    else:
        LMTD = ((Thin-Tcout) - (Thout-Tcin)) / np.log((Thin-Tcout) / (Thout-Tcin))

    if F is None:
        R = (Thin - Thout) / (Tcout - Tcin)
        P = (Tcout - Tcin) / (Thin - Tcin)
        if R == 1:
            W = (Ns - Ns * P) / (Ns - Ns * P + P)
            F = (np.sqrt(2) * (1 - W) / W) / \
                (np.log(W / (1 - W) + 1 / np.sqrt(2)) / np.log(W / (1 - W) - 1 / np.sqrt(2)))
        else:
            W = pow((1 - P * R) / (1 - P), 1 / Ns)
            S = np.sqrt(R ** 2 + 1) / (R - 1)
            F = S * np.log(W) / np.log((1 + W - S + S * W) / (1 + W + S - S * W))

    A = Q / (U * F * LMTD)

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
    HX.A = A

    return A, F, HX


def sizeHE_cooler(mh, cph, Thin, Thout, Tcin, Tcout, U, F=None, Ns=1):

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

    mh /= 3600  # convert kg/h to kg/s

    Q = mh * cph * (Thin - Thout)  # calculate heat transfer rate

    if (Thout - Tcin) == (Thin - Tcout):
        LMTD = Thout - Tcin
    else:
        LMTD = ((Thin - Tcout) - (Thout - Tcin)) / np.log((Thin - Tcout) / (Thout - Tcin))

    if F is None:
        R = (Thin - Thout) / (Tcout - Tcin)
        P = (Tcout - Tcin) / (Thin - Tcin)
        if R == 1:
            W = (Ns - Ns * P) / (Ns - Ns * P + P)
            F = (np.sqrt(2) * (1 - W) / W) / \
                (np.log(W / (1 - W) + 1 / np.sqrt(2)) / np.log(W / (1 - W) - 1 / np.sqrt(2)))
        else:
            W = pow((1 - P * R) / (1 - P), 1 / Ns)
            S = np.sqrt(R ** 2 + 1) / (R - 1)
            F = S * np.log(W) / np.log((1 + W - S + S * W) / (1 + W + S - S * W))

    A = Q / (U * F * LMTD)

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
    HX.A = A

    return A, F, HX
