
.. incomp_solns

Incompressible Solutions
========================


The `CoolProp <http://www.coolprop.org/dev/index.html>`_ project not only supports 
`Pure Fluids <http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids>`_
wrapped by the EngCoolProp **EC_Fluid** object
:ref:`link_ec_fluid_functions`, 
but also  `Incompressible Solutions <http://www.coolprop.org/fluid_properties/Incompressibles.html#massmix>`_
(e.g. Brines and Solutions) wrapped by the EngCoolProp **EC_Incomp_Soln** object
:ref:`link_ec_incomp_soln_functions`.

For incompressible solutions, EngCoolProp uses units of primarily inch, lbm, lbf, sec, BTU (some use of ft and hour).::

    The following are the default units for each property.

    T    = Temperature = degR
    P    = Pressure = psia
    D    = Density = lbm/cu ft
    rho  = Density = lbm/cu inch
    E    = Internal Energy = BTU/lbm
    H    = Enthalpy = BTU/lbm
    S    = Entropy = BTU/lbm degR
    Cp   = Heat Capacity (const. P) = BTU/lbm degR
    V    = Viscosity = 1.0E5 * lb/ft-sec
    C    = Thermal Conductivity = BTU/ft-hr-R
    

Default State Point
-------------------

Create a listing of properties at the default state point.
(i.e. T=(Tmax+Tmin)/2, P=Pmax/10)::
    
    from engcoolprop.ec_incomp_soln import EC_Incomp_Soln

    # Create incompressible solution object. (without specifying state point)
    ec_soln = EC_Incomp_Soln(symbol="MEG-30%")

    # Print state point
    ec_soln.printProps() # Print state point

Because the freezing point of MEG-30% is above the published minimum temperature 
for generic MEG 0% to 60%, the minimum permissible temperature for MEG-30%
was raised to just above the freezing point.

Output::

    NOTICE: Tmin=311.7 degR has been increased to T_freeze + 1 = 466.4 degR
    State Point for fluid INCOMP::MEG-30% (MEG-30%)
    T =      569  degR,                              Range( 466.434 -   671.67) degR
    P =     1000  psia                               Range(       0 -    10000) psia
    D =  64.1291  lbm/cuft                           Range( 61.7671 -  65.4448) lbm/cuft
    E =  36.5211  BTU/lbm                            Range(-55.0162 -  131.546) BTU/lbm
    H =  39.4067  BTU/lbm                            Range(-53.6702 -  152.698) BTU/lbm
    S =0.0667365  BTU/lbm degR                       Range(-0.110884 - 0.220262) BTU/lbm degR
    Cp= 0.903676  BTU/lbm degR                       Range(0.863298 - 0.936644) BTU/lbm degR
    V =  80.8991  viscosity [1.0E5 * lbm/ft-sec]     Range( 33.4415 -  524.521)
    C = 0.280781  thermal conductivity [BTU/ft-hr-R] Range( 0.24988 - 0.307045)
        T_freeze = 465.434 degR
        rho      = 0.0371118  lbm/cuin               Range(0.035745 - 0.037873) lbm/cuin
        mass%    =        30 base mass percent       Range(0% - 60%)


Can also print short forms of properties as::

    ec_soln.printTPD()

    ec_soln.printTransport()


Output::    

    INCOMP::MEG-30% T= 569.0 P=1000.0 D=64.1291 E= 36.52 H= 39.41 S=0.067
    INCOMP::MEG-30% Cp=0.903676 Visc=80.8991 ThCond=0.280781


State Point
-----------

Create a listing of properties at a given T and P. Note that Pmax is specified.

Pmax is the highest pressure considered in any iterative calcs. 
The default value for Pmax is 10,000 psia(see Range of P above).
It is usually best to keep Pmax above the max pressure being analyzed.::

    from engcoolprop.ec_incomp_soln import EC_Incomp_Soln

    # Create incompressible soln object at T=500 degR, P=500 psia and max pressure = 5000 psia
    ec_soln = EC_Incomp_Soln(symbol="MEG-30%", T=500, P=500, Pmax=5000) # T=degR, P=psia

    # OR... After ec_soln has been crated
    # ec_soln.setTP( 500, 500)

    # Print state point
    ec_soln.printProps()

Output::

    NOTICE: Tmin=311.7 degR has been increased to T_freeze + 1 = 466.4 degR
    State Point for fluid INCOMP::MEG-30% (MEG-30%)
    T =       500  degR,                              Range( 466.434 -   671.67) degR
    P =       500  psia                               Range(       0 -     5000) psia
    D =   65.1491  lbm/cuft                           Range( 61.7671 -  65.4448) lbm/cuft
    E =  -24.5781  BTU/lbm                            Range(-54.3432 -  131.546) BTU/lbm
    H =  -23.1579  BTU/lbm                            Range(-53.6702 -  142.122) BTU/lbm
    S =-0.0477664  BTU/lbm degR                       Range(-0.109441 - 0.220262) BTU/lbm degR
    Cp=  0.877102  BTU/lbm degR                       Range(0.863298 - 0.936644) BTU/lbm degR
    V =   242.449  viscosity [1.0E5 * lbm/ft-sec]     Range( 33.4415 -  524.521)
    C =    0.2604  thermal conductivity [BTU/ft-hr-R] Range( 0.24988 - 0.307045)
        T_freeze = 465.434 degR
        rho      =   0.037702  lbm/cuin               Range(0.035745 - 0.037873) lbm/cuin
        mass%    =         30 base mass percent       Range(0% - 60%)


