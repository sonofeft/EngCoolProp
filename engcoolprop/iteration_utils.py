import traceback
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from engcoolprop.ec_fluid import ( Peng_fromSI, PropsSI , TSI_fromEng, PSI_fromEng, UHeng_fromSI,
                                   Teng_fromSI, Deng_fromSI, UHSI_fromEng, DSI_fromEng  )
# from engcoolprop.find_exception_threshold import find_exception_limit

from engcoolprop.ec_fluid import toSI_callD

from engcoolprop.root_solver import brentq

# ==================================================================================
def calc_Psat_psia( TdegR, fluid ):
    Psat_psia = Peng_fromSI( PropsSI('P','T', TSI_fromEng( TdegR ) ,'Q',0, fluid) )
    return Psat_psia

# ==================================================================================
def find_T_from_Dterp(ec_inc, Dtarg, tol=1e-12, max_iter=1000):
    """
    Using the Interpolator ec_inc.Dterp (D=f(T)), calculate T for a target D
    Dtarg = Lbm/cuft
    Return T = degR
    """


    D_low = ec_inc.Dmin
    D_high = ec_inc.Dmax

    if Dtarg < D_low:
        print( 'WARNING: Dtarg=%g too low in find_T_from_D. Returning Tmax=%g R'%( Dtarg, ec_inc.Tmax ))
        return ec_inc.Tmax, 0 # all good (returning TdegR)
    if Dtarg > D_high:
        print( 'WARNING: Dtarg=%g too high in find_T_from_D. Returning Tmin=%g R'%( Dtarg, ec_inc.Tmin ))
        return ec_inc.Tmin, 0 # all good (returning TdegR)
    # ............................................................

    def get_d_param_at_T( T ):
            D = ec_inc.Dterp( T )
            return D - Dtarg # finds 0 point

    sol = brentq(get_d_param_at_T, ec_inc.Tmin, ec_inc.Tmax, xtol=tol )

    TdegR = sol.root

    if sol.converged:
        return TdegR, 0 # all good
    else:
        return TdegR, 1 # error


# ==================================================================================
def find_T_from_D(ec_inc, D, tol=1e-7, max_iter=1000):
    # Define the search range for temperature in Kelvin
    # T_low_si = ec_inc.Tmin_si  # Lower bound of temperature in K
    # T_high_si = ec_inc.Tmax_si  # Upper bound of temperature in K

    Dtarg_si = DSI_fromEng( D  ) # Do calcs in SI

    D_low_si = ec_inc.Dmin_si
    D_high_si = ec_inc.Dmax_si

    if Dtarg_si < D_low_si:
        print( 'WARNING: Dtarg_si=%g too low in find_T_from_D. Returning Tmax=%g R'%( Dtarg_si, ec_inc.Tmax ))
        return ec_inc.Tmax, 0 # all good (returning TdegR)
    if Dtarg_si > D_high_si:
        print( 'WARNING: Dtarg_si=%g too high in find_T_from_D. Returning Tmin=%g R'%( Dtarg_si, ec_inc.Tmin ))
        return ec_inc.Tmin, 0 # all good (returning TdegR)
    # ............................................................

    def get_d_param_at_T( Tsi ):
        if 1:#try:
            Dsi = CP.PropsSI('D', 'T', Tsi, 'P', ec_inc.Pmax_si, ec_inc.fluid)
            return Dsi - Dtarg_si # finds 0 point
        # except:
        #     return float('inf')

    sol = brentq(get_d_param_at_T, ec_inc.Tmin_si, ec_inc.Tmax_si, xtol=tol )

    TdegR = Teng_fromSI(sol.root)

    if sol.converged:
        return TdegR, 0 # all good
    else:
        return TdegR, 1 # error




# ==================================================================================
def find_T_at_P( ec_inc, P, dep_name='H', dep_val=0, tol=1.0E-6):
    """
    Iterate on T to find the value of the dependent variable for given P
    """

    TargSI = toSI_callD[dep_name]( dep_val ) # get target from Eng units to SI

    def get_opt_targ_of_T( T ): # T in degR
        Psat = ec_inc.get_Psat( T )
        Psi = PSI_fromEng( max(P, Psat + 0.0000001) ) # force fluid to Liquid if below Psat line
        try:
            SI_prop = PropsSI(dep_name, 'T',TSI_fromEng(T),'P',Psi,'INCOMP::%s'%ec_inc.symbol)
            return SI_prop - TargSI # finds 0 point
        except:
            return float('inf')

    sol = brentq(get_opt_targ_of_T, ec_inc.Tmin, ec_inc.Tmax, xtol=tol )

    # print( "sol.function_calls=%s, sol.iterations=%s"%(sol.function_calls, sol.iterations),
    #        'for', "dep_name=%s, dep_val=%s"%(dep_name, dep_val) )

    TdegR = sol.root

    if sol.converged:
        return TdegR, 0 # all good
    else:
        return TdegR, 1 # error


# ==================================================================================
def calc_T_freeze( ec_inc  ):
    """ 
    Use brentq method to solve for T_freeze
    (code taken from scipy)
    see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html
    or: https://en.wikipedia.org/wiki/Brent%27s_method
    for background info.
    """

    # find Tnbp where Psat = 1 atm (14.6959 psia)

    Tmin_si =   PropsSI('Tmin','T',0,'P',0,'INCOMP::%s'%ec_inc.symbol) * 1.0000000000000002
    Tmax_si =   PropsSI('Tmax','T',0,'P',0,'INCOMP::%s'%ec_inc.symbol) * 0.9999999999999999
    Tmin = Teng_fromSI( Tmin_si )
    Tmax = Teng_fromSI( Tmax_si )

    T = int( Tmin + 1)
    Psi = PSI_fromEng( 5000 )

    target_str = 'below the freezing point of'

    while T < Tmax:
        try:
            PropsSI('D', 'T',TSI_fromEng(T),'P',Psi, ec_inc.fluid)
        except:
            tb_str = traceback.format_exc()
            if target_str in tb_str:
                sL = tb_str.split( target_str )
                sL = sL[1].split()
                s = sL[0][:-1]
                try:
                    T_freeze = float( s.strip() )
                    return T_freeze
                except:
                    pass
        T += 1

    return 0



# ==================================================================================
def calc_Tnbp( ec_inc, method=None  ):
    """ 
    Use brentq method to solve for Tnbp
    (code taken from scipy)
    see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html
    or: https://en.wikipedia.org/wiki/Brent%27s_method
    for background info.
    """
    # only calc Tnbp if Psat is supported for ec_inc
    if not  ec_inc.get_Psat( ec_inc.Tmax ) > 0:
        return 0.0, 1 # error

    # find Tnbp where Psat = 1 atm (14.6959 psia)

    f = lambda T: calc_Psat_psia( T, ec_inc.fluid) - 14.6959 # find t(T) == 0

    tol=1.0e-6
    sol = brentq(f, ec_inc.Tmin, ec_inc.Tmax, xtol=tol )

    Tnbp = sol.root

    # psia = calc_Psat_psia( Tnbp, ec_inc.fluid)
    # print( "Tnbp=%.1f, psia=%.1f, sol.method=%s"%(Tnbp, psia, sol.method) )

    if sol.converged:
        return Tnbp, 0 # all good
    else:
        return Tnbp, 1 # error

# ==================================================================================
if __name__ == "__main__":

    from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid

    ec_inc = EC_Incomp_Fluid(symbol="DowQ", P=100.0, show_warnings=1 ) # T=degR, P=psia

    Tnbp, error_code = calc_Tnbp( ec_inc, method=None )

    if not error_code:
        print( 'For fluid = %s, Tnbp = %g degR'%( ec_inc.fluid, Tnbp ) )
    else:
        print( 'For fluid = %s, Tnbp calculation FAILED.'%ec_inc.fluid )
    print( '='*22, "Tnbp Calcs", '='*22 )

    EPSILON = 1.0e-5

    N = 10
    t_range = ec_inc.Tmax - ec_inc.Tmin
    tL = [ec_inc.Tmin + (i/N)*t_range for i in range(N+1)] 

    p_range = ec_inc.Pmax - ec_inc.Pmin
    pL = [ec_inc.Pmin + (i/N)*p_range for i in range(N+1)] 

    print( '='*22, "T from H Calcs", '='*22 )
    for T in tL:
        for P in pL:
            
            ec_inc.setTP( T, P)

            Htarg = ec_inc.H 
            TdegR, err_flag =find_T_at_P( ec_inc, P, dep_name='H', dep_val=Htarg)

            error = T - TdegR
            if abs(error) > EPSILON:
                print( 'T=%g'%TdegR, 'for Htarg=%g'%Htarg, ' and P=%g'%P, "T=%g"%T , 'error=', error, 'err_flag=', err_flag)


    print( '='*22, "T from D Calcs", '='*22 )
    for T in tL:
        for P in pL:
            
            ec_inc.setTP( T, P)

            TdegR, err_flag = find_T_from_Dterp(ec_inc, ec_inc.D, max_iter=1000)

            error = T - TdegR
            if abs(error) > EPSILON:
                print( 'T=%g'%TdegR, 'for Dtarg=%g lbm/cuft'%ec_inc.D, ' and P=%g'%P, "T=%g"%T , 'error=', error, 'err_flag=', err_flag)

    # Dtarg = ec_inc.D
    # P = 70
    # TdegR, err_flag =find_T_at_P( ec_inc, P, dep_name='D', dep_val=Dtarg)
    # print( 'T=%g'%TdegR, 'for Dtarg=%g'%Dtarg, ' and P=%g'%P)

    # ec_inc.setTP( TdegR, P)
    # ec_inc.printProps()




