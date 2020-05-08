import dsg

# Example 1: Horizontal pressure vessel
# Po = 470 psig, To = 800 degF, Di = 78 in, L = 480 in
V100 = dsg.designhorzpres(Di=78., L=480., Po=470., To=800.)
print(V100)

# Example 2: Vertical pressure vessel
# Po = 95.5 psig, To = 150 degF, Di = 120 in, L = 2544 in
V200 = dsg.designhorzpres(Di=120., L=2544., Po=95.5, To=150)
print(V200)

# Example 3: Vacuum vessel
# Po = 7.977 psig, To = 257 degF, Di = 168 in, L = 1080 in
V300 = dsg.designvac(Di=168., L=1080., Po=7.977, To=257.)
print(V300)

# Example 4: Compressor sizing for required power and outlet temperature
# m = 1e5 kg/h, P1 = 100 kPa, P2 = 300 kPa, T1 = 323.15 K,
# cp = 1.02, cv = 0.72, Z = 0.99
_, _, _, K400 = dsg.sizecompressor(m=1e5, P1=100, P2=300, T1=323.15, cp=1.02, cv=0.72, Z=0.99)
print(K400)

# Example 5: Pump sizing for required power
# Q = 35 m^3/h, dP = 500 kPa, rho = 1000 kg/m^3 (default)
_, _, P500 = dsg.sizepump(Q=35, dP=500)
print(P500)

# Example 6: Heat exchanger sizing for required heat exchange area (for
# heating process stream)
# mc = 31715 kg/h, cpc = 3246 J/(kg.K), Tcin = 89 degC, Tcout = 101 degC,
# Thin = 160 degC, Thout = 156 degC, U = 850 W/(m^2.degC), Ns = 1 (default)
_, _, HX600 = dsg.sizeHE_heater(mc=31715, cpc=3246, Tcin=89, Tcout=101, Thin=160, Thout=156, U=850)
print(HX600)

# Example 7: Heat exchanger sizing for required heat exchange area (for
# cooling process stream)
# mc = 31715 kg/h, cph = 3246 J/(kg.K), Thin = 89 degC, Thout = 60 degC,
# Tcin = 5 degC, Tcout = 10 degC, U = 850 W/(m^2.degC), Ns = 1 (default)
_, _, HX700 = dsg.sizeHE_cooler(mh=31715, cph=3246, Thin=89, Thout=60, Tcin=5, Tcout=10, U=850)
print(HX700)