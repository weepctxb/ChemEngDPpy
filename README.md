# ChemEngDPpy
Python code for mechanical design, sizing &amp; capex/opex calculations

The intention is to use these quickly as "black-box" functions for the purposes for TAC optimization in the second half of the detail design of the main equipment (e.g. PFR length, selection of pumps/compressors etc.)

- ***QUICK START***:

  - Download the required .py files directly to your project path; or

  - Use Git or checkout with SVN using the git URL.
  
  - Then, see exampleruns.py for usage examples (it's very intuitive)
  
To use, download the required .py files into the same folder as your existing code, and call the following functions depending on your situation.

Quick tip: Always use the function documentation to check the required units to avoid unit conversion errors!

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
   - ``CS`` (carbon steel)
   - ``SS`` (stainless steel)
   - ``Ni`` (nickel alloy)
   
 - `pump`
    - `reciprocating`, `positivedisp`
        - `Fe` (cast iron)
        - `CS`
        - `SS`
        - `Ni`
        - ti (titanium alloy)
    - `centrifugal`
        - `Fe`
        - `CS`
        - `SS`
        - `Ni`
       
 - `heatexc`
    - `fixedtube`, `utube`, `kettle`, `doublepipe`, `multipipe`
       - `CS`/`CS`
       - `CS`/`SS` and `SS`/`CS` (shell/tube order does not matter)
       - `SS`/`SS`
       - `CS`/`Ni` and `Ni`/`CS`
       - `Ni`/`Ni`
       - `CS`/`Ti` and `Ti`/`CS`
       - `Ti`/`Ti`
       
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

- `capex.annualcapex` - calculates total annualised capital cost estimated (`ACC`), based on an assumed payback period (`pbp`)

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

- `opex.costofmanfc` - calculates all components of annualised total cost of manufacture (`COM`)

To call intermediate functions (e.g. calculate number of operators per shift etc.), refer to documentation within code.