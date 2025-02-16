from CoolProp.CoolProp import PropsSI
from engcoolprop.ec_fluid import ( Peng_fromSI, PropsSI , TSI_fromEng, PSI_fromEng, UHeng_fromSI,
                                   Teng_fromSI, Deng_fromSI,  )
from engcoolprop.find_exception_threshold import find_exception_limit

from scipy import optimize

def calc_Psat_psia( TdegR, fluid ):
    Psat_psia = Peng_fromSI( PropsSI('P','T', TSI_fromEng( TdegR ) ,'Q',0, fluid) )
    return Psat_psia

def calc_Tnbp( ec_inc, method=None  ):
    """ 
    see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html
    for method options.
    If method == None, uses "brentq"
    """

    # find Tnbp where Psat = 1 atm (14.6959 psia)

    f = lambda T: calc_Psat_psia( T, ec_inc.fluid) - 14.6959 # find t(T) == 0

    tol=1.0e-4
    sol = optimize.root_scalar(f, bracket=[ec_inc.Tmin, ec_inc.Tmax], xtol=tol )

    Tnbp = sol.root

    # psia = calc_Psat_psia( Tnbp, ec_inc.fluid)
    # print( "Tnbp=%.1f, psia=%.1f, sol.method=%s"%(Tnbp, psia, sol.method) )

    if sol.converged:
        return Tnbp, 0 # all good
    else:
        return Tnbp, 1 # error

if __name__ == "__main__":

    from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid

    ec_inc = EC_Incomp_Fluid(symbol="DowQ", P=100.0 ) # T=degR, P=psia

    Tnbp, error_code = calc_Tnbp( ec_inc, method=None )

    if not error_code:
        print( 'For fluid = %s, Tnbp = %g degR'%( ec_inc.fluid, Tnbp ) )
    else:
        print( 'For fluid = %s, Tnbp calculation FAILED.'%ec_inc.fluid )



