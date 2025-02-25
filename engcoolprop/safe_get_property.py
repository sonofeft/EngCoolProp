import traceback
from CoolProp.CoolProp import PropsSI
from engcoolprop.find_exception_threshold import find_exception_limit
# from engcoolprop.parameter_units import si_unitsD
from engcoolprop.ec_fluid import  Peng_fromSI, PSI_fromEng

prop_descD = {} # key:eng prop desc, value:coolprop prop desc
prop_descD ['T'] = 'T'
prop_descD ['P'] = 'P'
prop_descD ['D'] = 'D'
prop_descD ['E'] = 'U'
prop_descD ['H'] = 'H'
prop_descD ['S'] = 'S'
prop_descD ['Cp'] = 'C'
prop_descD ['Visc'] = 'V'
prop_descD ['Cond'] = 'L'

# NOTE: 'P' is always an independent variable
LEGAL_PROP_CHR = set( ['T', 'D', 'U', 'H', 'S', 'C', 'V', 'L', 'P'] )


def get_si_prop( Psi_val, targ_prop='H', ind_name='T', ind_si_val=1000.0, symbol='Water'):
    """
    Psi_val (float): Pressure in SI units (Pa)
    target_prop (string): property to be calculated (e.g. 'T', 'D', 'U', 'H', 'S', 'C', 'V', 'L')
    ind_name (string): independent variable name (e.g. 'T', 'D', 'H', 'S')
    ind_si_val (float): independent variable value in SI units (e.g. K, kg/m³, J/kg, J/kg.K)
    symbol (string): symbol of incompressible fluid

    NOTE: this will often throw an Exception if the pressure is below saturation pressure
    """
    # print( 'Psi_val, targ_prop, ind_name, ind_si_val =', Psi_val, targ_prop, ind_name, ind_si_val )
    
    SI_prop = PropsSI(targ_prop, ind_name,ind_si_val,'P',Psi_val,'INCOMP::%s'%symbol)
    return SI_prop


def safe_get_INCOMP_prop( prop_desc, Psi_val=100000, ind_name='T', ind_si_val=1000.0, 
                          symbol='Water', show_warnings=True, Pmax=10000 ):
    """
    Return desired property from PropsSI for INCOMP fluid with Exception protection
    from an invalid Pressure. (e.g. below saturation pressure)

    Args:
        prop_desc (string): property to calculate: can be T, D, E, H, S, Cp, Visc, Cond    
        Psi_val (float): SI pressure value in "Pa"
        ind_name (string): independent variable other than P. (i.e. D, H, S, T)
        ind_si_val (float): independent variable value 
            D = "kg/m^3"
            H = "J/kg"
            S = "J/kg/K"
            T = "K"
        show_warnings (bool): will show warnings if == True
        Pmax (float): psia max pressure searched 
            
    ...CoolProp Incompressible Properties...
    Temperature (T): Kelvin (K)
    Pressure (P): Pascal (Pa)
    Density (D): Kilograms per cubic meter (kg/m³)
    Specific Heat Capacity (C): Joules per kilogram per Kelvin (J/kg.K)
    Internal Energy (U): Joules per kilogram (J/kg)
    Enthalpy (H): Joules per kilogram (J/kg)
    Entropy (S): Joules per kilogram per Kelvin (J/kg.K)
    Viscosity (V): Pascal-seconds (Pa.s)
    Thermal Conductivity (L): Watts per meter per Kelvin (W/m.K)
    """

    # if Psi_val == 0.0:
    #     print( 'Entered safe_get... at Psi_val=',Psi_val, "Peng_fromSI(Psi_val)=", Peng_fromSI(Psi_val))

    # get legal CoolProp description from a more generous prop_descD dictionary
    if prop_desc in prop_descD:
        prop_chr = prop_descD[ prop_desc ] # e.g. prop_descD ['Cond'] = 'L'
    else:
        prop_chr = prop_desc

    # demand that prop_chr is in legal CoolProp values
    if prop_chr not in  LEGAL_PROP_CHR: # i.e. in set( ['T', 'D', 'U', 'H', 'S', 'C', 'V', 'L', 'P'] )
        raise ValueError( 'prop_chr "%s" must be in '%prop_chr + repr(LEGAL_PROP_CHR) )
    
    
    try:
        prop_val =  get_si_prop( Psi_val, targ_prop=prop_desc, ind_name=ind_name, ind_si_val=ind_si_val, symbol=symbol )
        good_Psi_val = Psi_val
    except:
        # If there is an exception, assume it is a pressure issue.
        # Find the pressure that just works (i.e. right on the Exception threshold)
        # print( 'Exception at Psi_val=',Psi_val, "Peng_fromSI(Psi_val)=", Peng_fromSI(Psi_val))
        # tb_str = traceback.format_exc()
        # print( tb_str )
        
        f = lambda P: get_si_prop( P, targ_prop=prop_desc, ind_name=ind_name, ind_si_val=ind_si_val, symbol=symbol )

        good_Psi_val = find_exception_limit(f, tolerance=1e-4, lower_bound=1.0e-9, upper_bound=PSI_fromEng(Pmax),
                                             show_warnings=show_warnings )

        # start interpreting the output of find_exception_limit
        psia_in = Peng_fromSI(Psi_val)
        if good_Psi_val is None:
            if show_warnings:
                print( 'REALLY bad error: good_Psi_val is None')
                # print( traceback.format_exc() )
            psia_2 = 10000.0 # psia
        else:
            psia_2  = Peng_fromSI(good_Psi_val)

        if show_warnings and Peng_fromSI(good_Psi_val - Psi_val) > 0.1: # only show warning if error is significant
            s = ind_name + 'P'
            print( s,'WARNING: P override (%s'%Psi_val + '-->%s Pa)'%good_Psi_val, 
                '(%.1f-->%.1f psia)'%(psia_in, psia_2), 'calculating:',prop_desc )

        prop_val =  get_si_prop( good_Psi_val, targ_prop=prop_desc, ind_name=ind_name, ind_si_val=ind_si_val, symbol=symbol )

    return prop_val, good_Psi_val

if __name__ == "__main__":
    from engcoolprop.ec_fluid import (toEng_callD,  DSI_fromEng ,  
                                      PropsSI ,   SSI_fromEng , TSI_fromEng ,  
                                      UHSI_fromEng )
    from engcoolprop.parameter_units import si_unitsD

    Psi_val = 0 # PSI_fromEng(14.696) # run through all possible calcs at 1 atm
    # print( 'Psi_val =', Psi_val)

    ind_valD = {} # key:T, D, H or S, value: SI value
    ind_valD['T'] = TSI_fromEng( 800.0) # 527.67 )
    ind_valD['D'] =  DSI_fromEng( 55 ) # 62.4273 )
    ind_valD['H'] = UHSI_fromEng( 0.0 )
    ind_valD['S'] = SSI_fromEng( 0.0 )

    for prop_desc in ['T', 'D', 'U', 'H', 'S', 'C', 'V', 'L']:
        print()
        for ind_name in ['T', 'D', 'H', 'S']: # TP, HP, SP, DP
            ind_si_val = ind_valD[ ind_name ]

            prop_val, good_Psi_val = safe_get_INCOMP_prop( prop_desc, Psi_val=Psi_val, ind_name=ind_name, 
                                                          ind_si_val=ind_si_val, symbol='Water', show_warnings=False )

            print( '%s:'%prop_desc, 'at P=%g Pa'%good_Psi_val, 'and %s=%g'%(ind_name, ind_si_val), '(%s)'%si_unitsD[ind_name], 
                   prop_desc + '=' + '%8g'%toEng_callD[prop_desc](prop_val), '(%s)'%si_unitsD[prop_desc] )
