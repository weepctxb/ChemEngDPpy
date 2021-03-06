import dsg
import capex
import opex
import time

"""Note: It is best to always use keyword (named) arguments for this ChemEngDPpy library to avoid ambiguity"""

"""Step 1: Design and size all equipment"""

# Example 1: Horizontal pressure vessel V100
# Di = 78 in, L = 480 in, Po = 470 psig, To = 800 degF, MOC = '317L',
V100 = dsg.designhorzpres(Di=78., L=480., Po=470., To=800., MOC='317L', id='V100')
print(V100)

# Example 2: Vertical pressure vessel V200
# Di = 120 in, L = 2544 in, Po = 95.5 psig, To = 150 degF, auto-select a suitable MOC
V200 = dsg.designhorzpres(Di=120., L=2544., Po=95.5, To=150, id='V200')
print(V200)

# Example 3: Vacuum vessel V300
# Di = 168 in, L = 1080 in, Po = 7.977 psig, To = 257 degF, specify vertical for costing purposes
V300 = dsg.designvac(Di=168., L=1080., Po=7.977, To=257., etype='vertical', id='V300')
print(V300)

# Example 4: Compressor K400 sizing for required power and outlet temperature
# m = 1e5 kg/h, P1 = 2 bar, P2 = 6 bar, T1 = 323.15 K,
# cp = 1.02, cv = 0.72, Z = 0.99
# Specify rotary compressor using carbon steel for costing purposes
K400_comppower, K400_compeff, K400_T2, K400 = \
    dsg.sizecompressor(m=1e5, P1=2, P2=6, T1=323.15, cp=1.02, cv=0.72, Z=0.99, etype='rotary', mat='CS', id='K400')
print(K400)

# Example 5: Compressor K400 sizing for required power and outlet temperature
# m = 4e5 kg/h, P1 = 1.5 bar, P2 = 6.5 bar, T1 = 290 K,
# cp = 1.02, cv = 0.72, Z = 0.99
# Specify rotary compressor using carbon steel for costing purposes
K500_comppower, K500_compeff, K500_T2, K500 = \
    dsg.sizecompressor(m=4e5, P1=1.5, P2=6., T1=290, cp=1.02, cv=0.72, Z=0.99, etype='centrifugal', mat='SS', id='K500')
print(K500)

# Example 6: Pump sizing for required power
# Q = 40 m^3/h, P1 = 200 kPa, P2 = 700 kPa, rho = 1200 kg/m^3
# Specify positive displacement pump using stainless steel for costing purposes
P600_pumppower, P600_pumpeff, P600 = \
    dsg.sizepump(Q=40, P1=200, P2=700, rho=1200, etype='positivedisp', mat='SS', id='P600')
print(P600)

# Example 7: Pump sizing for required power
# Q = 100 m^3/h, P1 = 180 kPa, P2 = 800 kPa, rho = 990 kg/m^3
# Specify positive displacement pump using stainless steel for costing purposes
P700_pumppower, P700_pumpeff, P700 = \
    dsg.sizepump(Q=100, P1=180, P2=800, rho=990, etype='centrifugal', mat='Ni', id='P700')
print(P700)

# Example 8: Heat exchanger sizing for required heat exchange area (for
# heating process stream)
# mc = 31715 kg/h, cpc = 3246 J/(kg.K), Tcin = 89 degC, Tcout = 101 degC,
# Thin = 160 degC, Thout = 156 degC, U = 850 W/(m^2.degC), Ns = 2
# Specify pressure = 2 bar, double pipe HX using CS (shell) and Ni (tube) for costing purposes
HX800_area, HX800_F, HX800 = dsg.sizeHE_heater(mc=31715, cpc=3246, Tcin=89, Tcout=101, Thin=160, Thout=156, U=850, P=2, \
                                Ns=2, etype='doublepipe', mat='CS/Ti', id='HX800')
print(HX800)

# Example 9: Heat exchanger sizing for required heat exchange area (for
# cooling process stream)
# mc = 31715 kg/h, cph = 3246 J/(kg.K), Thin = 89 degC, Thout = 60 degC,
# Tcin = 5 degC, Tcout = 10 degC, U = 850 W/(m^2.degC), Ns = 1 (default)
# Specify pressure = 4 bar, kettle reboiler HX using CS (shell) and SS (tube) for costing purposes
HX900_area, HX900_F, HX900 = dsg.sizeHE_cooler(mh=31715, cph=3246, Thin=89, Thout=60, Tcin=5, Tcout=10, U=850, P=4, \
                                etype='kettle', mat='CS/SS', id='HX900')
print(HX900)

"""Step 2: Calculate all equipment capital"""

time.sleep(0.1)

# Example 8: Calculating capital of equipment, assuming greenfield project
eqptlist = [V100, V200, V300, K400, K500, P600, P700, HX800, HX900]
FCI, capexreport = capex.econreport(eqptlist, planttype='green', pbp=3, year=2019, currency='SGD', \
                                    reporttype='numpy', verbose=True)

"""Step 3: Calculate all manufacturing costs"""

time.sleep(0.1)

# Example 9: Calculating cost of manufacturing

COL, Nop = opex.labourcost(wage=42750., eqptlist=eqptlist)

print('Number of workers required : ' + str(int(Nop)))

CRM = opex.costofraw(rawmaterialtuple=(19779., 12606., 325., 240.), unitpricetuple=(89.94, 0.25, 477.74, 1512.39))

CUT = opex.costofutil(utiltuple=(0., 671069., 30815., 889354., 105723., 0., 0., 30643.), year=2019, currency='SGD')

CWT = 0.  # Waste treatment cost has to be calculated manually

COMd, d, COM, DMC, FMC, GE, report_dict = opex.costofmanfc(FCI=FCI, COL=COL, CRM=CRM, CWT=CWT, CUT=CUT, verbose=True)