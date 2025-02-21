from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from engcoolprop.ec_fluid import ( Peng_fromSI, PropsSI , TSI_fromEng, PSI_fromEng, UHeng_fromSI,
                                   Teng_fromSI, Deng_fromSI, UHSI_fromEng, DSI_fromEng  )
# from engcoolprop.find_exception_threshold import find_exception_limit
from engcoolprop.safe_get_property import safe_get_INCOMP_prop as safe_get_prop
from engcoolprop.ec_fluid import toSI_callD
from scipy import optimize

# ==================================================================================
def calc_Psat_psia( TdegR, fluid ):
    Psat_psia = Peng_fromSI( PropsSI('P','T', TSI_fromEng( TdegR ) ,'Q',0, fluid) )
    return Psat_psia

# ==================================================================================
def find_T_from_D(ec_inc, D, tol=1e-7, max_iter=1000):
    # Define the search range for temperature in Kelvin
    T_low_si = ec_inc.Tmin_si  # Lower bound of temperature in K
    T_high_si = ec_inc.Tmax_si  # Upper bound of temperature in K

    Dtarg_si = DSI_fromEng( D  ) # Do calcs in SI

    D_low_si = ec_inc.Dmin_si
    D_high_si = ec_inc.Dmax_si

    if Dtarg_si < D_low_si:
        print( 'WARNING: Dtarg_si=%g too low in find_T_from_D. Returning Tmax=%g R'%( Dtarg_si, ec_inc.Tmax ))
        return ec_inc.Tmax, 0 # all good (returning TdegR)
    if Dtarg_si > D_high_si:
        print( 'WARNING: Dtarg_si=%g too high in find_T_from_D. Returning Tmin=%g R'%( Dtarg_si, ec_inc.Tmin ))
        return ec_inc.Tmin, 0 # all good (returning TdegR)


    for i in range(max_iter):
        T_mid_si = (T_low_si + T_high_si) / 2
        try:
            density = CP.PropsSI('D', 'T', T_mid_si, 'P', ec_inc.Pmax_si, ec_inc.fluid)
        except ValueError as e:
            if "below the freezing point" in str(e):
                print(T_mid_si, "below the freezing point")
                T_low_si = T_mid_si
                continue
            else:
                raise e
        
        if abs(density - Dtarg_si) < tol:
            # print( 'Number of iterations in find_T_from_D =', i)
            return Teng_fromSI(T_mid_si), 0 # all good (returning TdegR)
        elif density < Dtarg_si:
            T_high_si = T_mid_si
        else:
            T_low_si = T_mid_si

    # print( 'Number of iterations in find_T_from_D =', i)
    return Teng_fromSI(T_mid_si), 1 # Failed (returning TdegR)


# ==================================================================================
def find_T_at_P( ec_inc, P, dep_name='H', dep_val=0, tol=1.0E-7):
    """
    Iterate on T to find the value of the dependent variable for given P
    """

    P = max( P, ec_inc.get_Psat(T) )

    Psi = PSI_fromEng( P )
    TargSI = toSI_callD[dep_name]( dep_val )
    def get_opt_targ_of_T( T ):
        try:
            SI_prop = PropsSI(dep_name, 'T',TSI_fromEng(T),'P',Psi,'INCOMP::%s'%ec_inc.symbol)
            return SI_prop - TargSI # finds 0 point
        except:
            return float('inf')
    
    try:
        opt_targ = get_opt_targ_of_T(ec_inc.Tmin  )
        # print( "Tmin opt_targ=%s at Tmin"%opt_targ,
        #        "dep_name=%s, dep_val=%s"%(dep_name, dep_val) )
        if abs(opt_targ) <= tol:
            # print( "opt_targ=%g at Tmin"%opt_targ )
            return ec_inc.Tmin, 0 # all good
    except:
        pass

    try:
        opt_targ = get_opt_targ_of_T(ec_inc.Tmax  )
        # print( "Tmax opt_targ=%s at Tmin"%opt_targ,
        #        "dep_name=%s, dep_val=%s"%(dep_name, dep_val) )
        if abs(opt_targ) <= tol:
            # print( "opt_targ=%g at Tmax"%opt_targ )
            return ec_inc.Tmax, 0 # all good
    except:
        pass

    sol = optimize.root_scalar(get_opt_targ_of_T, bracket=[ec_inc.Tmin, ec_inc.Tmax], xtol=tol )

    # print( "sol.function_calls=%s, sol.iterations=%s"%(sol.function_calls, sol.iterations),
    #        'for', "dep_name=%s, dep_val=%s"%(dep_name, dep_val) )

    TdegR = sol.root

    if sol.converged:
        return TdegR, 0 # all good
    else:
        return TdegR, 1 # error


# ==================================================================================
def calc_Tnbp( ec_inc, method=None  ):
    """ 
    see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html
    for method options.
    If method == None, uses "brentq"
    """

    # find Tnbp where Psat = 1 atm (14.6959 psia)

    f = lambda T: calc_Psat_psia( T, ec_inc.fluid) - 14.6959 # find t(T) == 0

    tol=1.0e-6
    sol = optimize.root_scalar(f, bracket=[ec_inc.Tmin, ec_inc.Tmax], xtol=tol )

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

            TdegR, err_flag = find_T_from_D(ec_inc, ec_inc.D, max_iter=1000)

            error = T - TdegR
            if abs(error) > EPSILON:
                print( 'T=%g'%TdegR, 'for Dtarg=%g lbm/cuft'%ec_inc.D, ' and P=%g'%P, "T=%g"%T , 'error=', error, 'err_flag=', err_flag)

    # Dtarg = ec_inc.D
    # P = 70
    # TdegR, err_flag =find_T_at_P( ec_inc, P, dep_name='D', dep_val=Dtarg)
    # print( 'T=%g'%TdegR, 'for Dtarg=%g'%Dtarg, ' and P=%g'%P)

    # ec_inc.setTP( TdegR, P)
    # ec_inc.printProps()




