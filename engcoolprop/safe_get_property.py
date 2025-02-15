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


def safe_get_INCOMP_prop( prop_desc, Psi_val=100000, ind_name='T', ind_si_val=1000.0, symbol='Water' ):
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
    # get legal CoolProp description from a more generous prop_descD dictionary
    if prop_desc in prop_descD:
        prop_chr = prop_descD[ prop_desc ]
    else:
        prop_chr = prop_desc

    # demand tht prop_chr is in legal CoolProp values
    if prop_chr not in  LEGAL_PROP_CHR:
        raise ValueError( 'prop_chr "%s" must be in '%prop_chr + repr(LEGAL_PROP_CHR) )
    
    
    try:
        prop_val =  get_si_prop( Psi_val, targ_prop=prop_desc, ind_name=ind_name, ind_si_val=ind_si_val, symbol=symbol )
        good_Psi_val = Psi_val
    except:
        # If there is an exception, assume it is a pressure issue.
        # Find the pressure that just works (i.e. right on the Exception threshold)
        
        f = lambda P: get_si_prop( P, targ_prop=prop_desc, ind_name=ind_name, ind_si_val=ind_si_val, symbol=symbol )

        good_Psi_val = find_exception_limit(f, tolerance=1e-4, lower_bound=1.0e-9, upper_bound=PSI_fromEng(5000) )
        print( 'WARNING: P changed from %8g Pa to'%Psi_val, '%8g Pa'%good_Psi_val, '(%g psia)'%Peng_fromSI(good_Psi_val) )

        prop_val =  get_si_prop( good_Psi_val, targ_prop=prop_desc, ind_name=ind_name, ind_si_val=ind_si_val, symbol=symbol )

    return prop_val, good_Psi_val

if __name__ == "__main__":
    from engcoolprop.ec_fluid import (toEng_callD,  DSI_fromEng ,  
                                      PropsSI ,   SSI_fromEng , TSI_fromEng ,  
                                      UHSI_fromEng )

    ind_valD = {} # key:T, D, H or S, value: SI value
    ind_valD['T'] = TSI_fromEng( 527.67 )
    ind_valD['D'] =  DSI_fromEng( 55 ) # 62.4273 )
    ind_valD['H'] = UHSI_fromEng( 0.0 )
    ind_valD['S'] = SSI_fromEng( 0.0 )

    for prop_desc in ['T', 'D', 'U', 'H', 'S', 'C', 'V', 'L']:
        print()
        for ind_name in ['T', 'D', 'H', 'S']: # TP, HP, SP, DP
            ind_si_val = ind_valD[ ind_name ]

            Psi_val = PSI_fromEng(14.696)
            # print( 'Psi_val =', Psi_val)

            prop_val, good_Psi_val = safe_get_INCOMP_prop( prop_desc, Psi_val=Psi_val, ind_name=ind_name, ind_si_val=ind_si_val, symbol='Water' )

            print( prop_desc, 'Psi_val=%g'%Psi_val, '%s=%g'%(ind_name, ind_si_val), 'prop_val =', '%8g'%toEng_callD[prop_desc](prop_val), prop_desc )
