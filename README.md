# ChemEngDPpy
Python code for mechanical design, sizing &amp; capex/opex calculations

The intention is to use these quickly as "black-box" functions for the purposes for TAC optimization in the
detailed design of the main equipment (e.g. PFR wall thickness requirement, sizing of pumps/compressors etc.)

- ***QUICK START***:

  - Download the required `.py` files directly to your project path; or

  - Use Git or checkout with SVN using the git URL.
  
  - Then, see `exampleruns.py` for usage examples (it's very intuitive)
  
To use, download the required `.py` files into the same folder as your existing code, and call the following
functions depending on your situation.

Quick tips:

- Always use the docstrings to check the required units to avoid unit conversion errors!

- Use named arguments as far as possible to avoid ambiguities!

------------------------------------------------

## Project Structure

- dsg: For Mechanical Design & Ancillary Equipment Sizing

- capex: For CAPEX Calculation and Reporting

- opex: For OPEX Calculation and Reporting

------------------------------------------------

## Implementation Methods

Quick tip: Always use the docstrings to check the required units to avoid unit conversion errors, and use
 named arguments as far as possible to avoid ambiguities!

Some definitions first:

- ***Design inputs***: Properties that are critical and required to conduct equipment design/sizing

- ***Design variables***: Properties that are determined after conducting equipment design/sizing

- ***Additional tags***: Properties that are not critical to equipment design/sizing,
but are important for the purposes of cost estimation.
These can be supplied in the equipment object creation for ease of cost estimation later

For example:

- *Pressure* is a *design input* for vessels (to determine final shell thickness),
a *design variable* for pumps/compressor (to determine power requirement),
and is an *additional tag* for heat exchangers (to determine pressure factor)

- *Volume* is a *design variable* for vessels, but is not relevant for other equipment categories

- Operating temperature is a *design input* for vessels (to determine final shell thickness based on MOC properties),
a *design input* for heat exchangers (to determine heat exchange area)
a *design variable* for compressors (specifically the outlet temperature),
and is not relevant to pumps (pump sizing in this library assumes negligible temperature change)

There are two major methods for using this library:

1. ***Method 1 (scalar framework)***: Using scalar input - scalar output to retrieve only the key design variables; or

2. ***Method 2 (OOP framework)***: Using the library's object-oriented programming (OOP) framework.
Generally, common chemical plant equipment (vessels, pumps, compressors, heat exchangers etc.)
are also created as optional outputs and alternative inputs for the various functions.

For example in equipment design, use:

```python
# Method 1 (scalar framework)
comppower, compeff, T2, _ = dsg.sizecompressor(m=1e5, P1=100, P2=300, T1=323.15, cp=1.02, cv=0.72, Z=0.99)
```

to just retrieve the power rating and temperature of the sized compressor as well as its estimated adiabatic efficiency
(as in Method 1), or use:

```python
# Method 2 (OOP framework)
_, _, _, K100 = dsg.sizecompressor(m=1e5, P1=100, P2=300, T1=323.15, cp=1.02, cv=0.72, Z=0.99,
etype='rotary', mat='CS', id='K400')
```

to create a `Compressor()` object containing all relevant design inputs (i.e. the inputs), design variables
(i.e. `comppower`, `compeff` and `T2`) and additional tags.

Additional tags can be supplied in the equipment object creation for ease of cost estimation later, such as:

- `category` for **equipment category** (e.g. pumps, compressors etc. - automatically created when the respective
`dsg.size(...)` or `dsg.design(...)` functions are called - see below or docstrings)

- `etype` for **equipment type** (e.g. centrifugal pumps, rotary pumps etc. - see `capex` documentation below for
list of supported chemical plant equipments)

- `mat` for **material type** (e.g. carbon steel, stainless steel, cast iron etc. - see `capex` documentation
below for list of supported material types)

- `id` to specify an **equipment name** for semantic purposes

- `P` to specify a designed operating pressure for equipment where pressure is neither a required design input
nor a design variable. Most notably, this include heat exchangers whereby the pressure specification is only to
determine the pressure factor `FP` for capital cost estimation of heat exchangers.

Key exceptions to this include:

- Mechanical design for pressure/vacuum vessels is *only* conducted using the OOP framework (Method 2), due the 
large number of critical design variables in the mechanical design.

- For reactors and distillation columns, design the pressure/vacuum vessel and internals (e.g. distillation trays, 
mixers/impellers, packings) *separately* (i.e. as separate objects), as combined object support is not available yet (#TODO).

- CAPEX reporting can *only* be conducted using the OOP framework (Method 2).

As another example in capital cost estimation for equipment, use

1. *Method 1 (scalar framework)*:

```python
# Method 1 (scalar framework)
# manually supply cost coefficients, min/max capacity for validity range (optional),
# exponential factor for extrapolation (optional)...
capex.eqptpurcost(A=20, Ktuple=(3.5565, 0.3776, 0.0905, 0.1, 628., 0.5))
# ... or just retrieve from capex.eqptcostlib
capex.eqptpurcost(A=20, Ktuple=capex.eqptcostlib['vessel']['horizontal'])
```

to size a horizontal vessel of volume `A` = 20 m^3, and then calling all the subsequent functions in order (see `capex`
documentation). However, using

2. *Method 2 (OOP framework)*:

```python
# Method 2 (OOP framework)
capex.eqptpurcost(eqpt=V100)
```

, where `V100` is the output (a `MechDesign()` object in this context) retrieved from `dsg.designvertpres`,
is so much easier.

In fact, using the OOP framework (Method 2), we can just create a list of all the relevant
equipment objects and perform the entire capex estimation directly:

```python
eqptlist = [V100, V200, V300, K400, K500, P600, P700, HX800, HX900]
FCI, capexreport = capex.econreport(eqptlist, planttype='green', pbp=3, year=2019, currency='SGD', \
                                    reporttype='numpy', verbose=True)
```

See `exampleruns.py` for more sample implementations.

Further possible uses of the OOP framework (Method 2) could be to:

- Perform optimization on equipment type or material selection,
using the tags `etype`, `mat` etc. as categorical decision variables for capex minimization

- Determine the best type of heat exchanger and utilities to use by minizing both capex and opex. This is similar
to above, with utilities selection being an additional categorical decision variable which influences both
heat exchange area and cost of utilities.

- Find the optimal configuration for multi-stage compression with minimal electricity consumption and stream cooling,
using the number of compressors, compression ratios, heat exchange area for cooling etc. as decision variables.

However, this implementation is left to the user to do. (#TODO add examples)

------------------------------------------------

## Mechanical Design & Ancillary Equipment Sizing - dsg

Constants:

- `dsg.Patm` = 14.696 - standard atmospheric pressure (psi)

- `dsg.Patmb` = 1.01325 - standard atmospheric pressure (bar)

- `dsg.Troom` = 77 - ambient temperature (degF)

- `dsg.tmin` = 1/4 - universal minimum allowable vessel thickness (in)

- `dsg.tc` = 0.125 - corrosion allowance (in) for both corrosive and non-corrosive conditions (default is 1/8)

- `dsg.rhosteel` = 0.2836 - density of SA-285C/SA-387B/carbon/low-alloy steels (lb/in^3)

- `dsg.g` = 9.80665 - standard Earth gravitational acceleration (m/s^2)

- `dsg.R` = 8.31446261815324 - universal ideal gas constant (J/(K.mol))

- `dsg.Ta` = 10. - minimum heat exchanger temperature approach (K)

Functions:

- `dsg.designhorzpres` - perform entire mechanical design for horizontal pressure vessels

- `dsg.designvertpres` - perform entire mechanical design for vertical pressure vessels

- `dsg.designvac` - perform entire mechanical design for vacuum vessels

- `dsg.sizecompressor` - conducts compressor sizing

- `dsg.sizepump` - conducts pump sizing

- `dsg.sizeHE_heater` - conducts heat exchanger sizing for heating stream

- `dsg.sizeHE_cooler` - conducts heat exchanger sizing for cooling stream

To call intermediate functions (e.g. calculate shell thickness, calculate max. allowable stress, calculate wind
allowance etc.), refer to documentation within code.

------------------------------------------------

## CAPEX Calculation - capex

Constants:

- `capex.CEPCI[20yy]` - retrieves annual CEPCI index for 20yy (where yy = 01, 18 or 19)

- `capex.USSG[20yy]` - retrieves USD:SGD forex rate for 20yy year-average (where yy = 01, 18 or 19)

- `capex.CPI['zz'][20yy]` - retrieves country zz's consumer price index (CPI) for 20yy year-average
(where yy = 01, 16, 18 or 19 and zz = 'SG' or 'US')

  - Reference year for CAPCOST = 2001 (as of Turton et al. 5th Ed.)

  - Reference year for utilities cost = 2016 (as of Turton et al. 5th Ed.)
  
- `capex.eqptcostlib['eqptcategory']['eqpttype']` - retrieves a tuple of equipment cost correlation parameters
`(K1, K2, K3, Amin, Amax, n)` where `(K1, K2, K3)` = cost correlation params, `(Amin, Amax)` = min/max capacity
(range of validity), and `n` = cost exponent used in the exponential costing rule, from which the
purchased equipment cost `Cpo` can then be calculated.

- `capex.pressurefaclib['eqptcategory']['eqpttype']` - retrieves a tuple of equipment pressure factor correlation parameters
`(C1, C2, C3, Pmin, Pmax)` where `(C1, C2, C3)` = pressure factor correlation params, and `(Pmin, Pmax)` = min/max pressure
(range of validity), from which the pressure factor `FP` can then be calculated.

- `capex.matfaclib['eqptcategory']['eqpttype']['mat']` - retrieves equipment material factor `FM` (a float)

- `capex.baremodlib['eqptcategory']['eqpttype']` - retrieves a tuple of equipment bare module correlation parameters
`(B1, B2)`, such that the bare module factor `FBM = B1 + B2 * FP * FM`.

Supported equipment categories, types and materials (`mat`):
 
 - `compressor`
   - `centrifugal`, `axial`, `reciprocating`, `rotary`
       - `CS` (carbon steel)
       - `SS` (stainless steel)
       - `Ni` (nickel alloy)
   
 - `pump`
    - `reciprocating`, `positivedisp`
        - `Fe` (cast iron)
        - `CS`
        - `SS`
        - `Ni`
        - `Ti` (titanium alloy)
    - `centrifugal`
        - `Fe`
        - `CS`
        - `SS`
        - `Ni`
       
 - `heatexc`
    - `fixedtube`, `utube`, `kettle`, `doublepipe`, `multipipe`
       - `CS/CS`
       - `CS/SS` and `SS/CS` (shell/tube order does not matter)
       - `SS/SS`
       - `CS/Ni` and `Ni/CS`
       - `Ni/Ni`
       - `CS/Ti` and `Ti/CS`
       - `Ti/Ti`
       
 - `vessel`
    - `horizontal`, `vertical`
        - `CS`
        - `SS`
        - `Ni`
        - `Ti`
        
 - `trays`
    - `sieve`, `valve`
        - `CS`
        - `SS`
        - `Ni`
    - `demister`
        - `SS`
        - `FC` (fluorocarbon)
        - `Ni`

Functions:

- `capex.eqptpurcost` - calculates purchased equipment cost (`Cpo`)

- `capex.pressurefacves` - calculates pressure factor for vessels (`FP`)

- `capex.pressurefacanc` - calculates pressure factor for ancillary equipment (`FP`)

- `capex.baremodfac` - calculates bare module factor (`FBM`)

- `capex.baremodcost` - calculates bare module cost (`CBM`)

- `capex.totmodcost` - calculates total module cost (`CTM`)

- `capex.grasscost` - calculates grassroots cost (`CGR`)

- `capex.annualcapex` - calculates total annualised capital cost estimated (`ACC`),
based on an assumed payback period (`pbp`). If a value of pbp is assumed, note that this should only be
used for ACC estimation for optimization purposes! Alternatively, calculate pbp based on projected revenue estimates.

------------------------------------------------

## OPEX Calculation - opex

Constants:

- `opex.SF` - retrieves stream factor for plant operation

- `opex.runtime` - retrieves operational runtime per annum

- `opex.shiftdur` - retrieves duration of one workshift

- `opex.shiftperweek` - retrieves number of shifts per year

- `opex.yearww` - retrieves number of work weeks per year after leave entitlements

- `opex.util['xxx']` - utility cost for xxx (USD/GJ basis), where xxx =

  - `"HPS"` for high-pressure steam (41 barg, 254 degC)

  - `"MPS"` for medium-pressure steam (10 barg, 184 degC)

  - `"LPS"` for low-pressure steam (5 barg, 160 degC)

  - `"CW"` for cooling water (30-45 degC)

  - `"ChW"` for chilled water (5 degC)

  - `"LTR"` for low temperature refrigerant (-20 degC)

  - `"VLTR"` for very low temperature refrigerant (-50 degC)

  - `"elec"` for electricity (110-440 V)

Functions:

- `opex.labourcost` - calculates annualised labour cost (`COL`)

- `opex.costofraw` - calculates annualised cost of raw materials (`CRM`)

- `opex.costofutil` - calculates annualised cost of utilities (`CUT`)

- `opex.costofwaste` - (placeholder function for user's customised waste treatment calculation)

- `opex.costofmanfc` - calculates all components of annualised total cost of manufacture (`COM`)

To call intermediate functions (e.g. calculate number of operators per shift etc.), refer to documentation within code.