from CoolProp.CoolProp import PropsSI
from engcoolprop.conv_funcs import ( Peng_fromSI, TSI_fromEng, PSI_fromEng, UHeng_fromSI,
                                   Teng_fromSI, Deng_fromSI,  )
from engcoolprop.find_exception_threshold import find_exception_limit

def run_experiments( symbol ):
    """Start iterating pressure values to find the reference H=0 point"""
    Tmin_si =   PropsSI('Tmin','T',0,'P',0,'INCOMP::%s'%symbol)  
    Tmax_si =   PropsSI('Tmax','T',0,'P',0,'INCOMP::%s'%symbol) 

    Tref_si = (Tmin_si + Tmax_si) / 2.0

    print( 'Tref = %g degR'%Teng_fromSI( Tref_si ), '  degK = %g'%Tref_si )

    # density is invariant with pressure for INCOMP
    Psi = PSI_fromEng( 1000.0 )
    Dmax_si =  PropsSI('D','T',Tmin_si,'P',Psi,'INCOMP::%s'%symbol) # Dmax at Tmin
    Dmin_si =  PropsSI('D','T',Tmax_si,'P',Psi,'INCOMP::%s'%symbol) # Dmin at Tmax

    Dmax = Deng_fromSI( Dmax_si )
    Dmin = Deng_fromSI( Dmin_si )
    print( 'Dmin = %g, Dmax = %g lbm/cuft'%(Dmin, Dmax) )

    # check reverse T calc
    Tmin_calc = Teng_fromSI( PropsSI('T','D',Dmax_si,'P',Psi, 'INCOMP::%s'%symbol ) )
    Tmax_calc = Teng_fromSI( PropsSI('T','D',Dmin_si,'P',Psi, 'INCOMP::%s'%symbol ) )

    print( 'Tmin = %g degR'%Teng_fromSI( Tmin_si ), '   Tmin_calc = %g'%Tmin_calc  )
    print( 'Tmax = %g degR'%Teng_fromSI( Tmax_si ), '   Tmin_calc = %g'%Tmax_calc  )

    Hmax = UHeng_fromSI( PropsSI('H','D',Dmin_si,'P',Psi, 'INCOMP::%s'%symbol ) )
    Hmin = UHeng_fromSI( PropsSI('H','D',Dmax_si,'P',0, 'INCOMP::%s'%symbol ) )
    print( '%s: '%symbol, 'Hmin = %g, Hmax = %g BTU/lbm'%(Hmin, Hmax) )

    Smax = UHeng_fromSI( PropsSI('S','D',Dmin_si,'P',Psi, 'INCOMP::%s'%symbol ) )
    Smin = UHeng_fromSI( PropsSI('S','D',Dmax_si,'P',Psi, 'INCOMP::%s'%symbol ) )
    print( '%s: '%symbol,  'Smin = %g, Smax = %g BTU/lbm degR'%(Smin, Smax) )
    
    print( '%s: '%symbol,  'Tmin_si = %g,  Tmax_si= %g'%( Tmin_si, Tmax_si) )
    t_range = Tmax_si - Tmin_si
    N = 9
    tL = [Tmin_si + (i/N)*t_range for i in range(N+1)]
    Psi2 = PSI_fromEng( 100.0 )
    for Tsi in tL:
        if 1:#try:
            Dsi =  PropsSI('D','T',Tsi,'P',Psi,'INCOMP::%s'%symbol)
            print(  '%s: '%symbol, 'D =', Deng_fromSI(Dsi) , 'lbm/cuft', '   at T = %g degK'%Tsi, ' = %g degR'%Teng_fromSI( Tsi )  )

            Hsi =  PropsSI('H','T',Tsi,'P',Psi,'INCOMP::%s'%symbol)
            print( '%s: '%symbol,  '    H =', UHeng_fromSI(Hsi) , 'BTU/lbm', '   at T = %g degK'%Tsi, ' = %g degR'%Teng_fromSI( Tsi )  )

            # Hsi2 =  PropsSI('H','T',Tsi,'P',Psi2,'INCOMP::%s'%symbol)
            # print( '    H =', UHeng_fromSI(Hsi2) , 'BTU/lbm', '   at T = %g degK'%Tsi, ' = %g degR'%Teng_fromSI( Tsi )  )

        # except:
        #     print( 'FAILED: P = %g'%P )

    def get_H( P, Tsi):
        Psi = PSI_fromEng( P )
        Hsi = PropsSI('H','T',Tsi,'P',Psi,'INCOMP::%s'%symbol)

    N = 9
    tL = [Tmin_si + (i/N)*t_range for i in range(N+1)]
    for Tsi in tL:
        f = lambda P: get_H( P, Tsi )
        good_limit_Psi = find_exception_limit(f, tolerance=1e-4, lower_bound=1.0e-9, upper_bound=10000)
        if good_limit_Psi is None:
            print( "==> FAILED:",  '%s: '%symbol, '   at T = %g degK'%Tsi, '(%g degR)'%Teng_fromSI( Tsi )  )
        else:
            print( '%s: '%symbol,  'good_limit_Psi = %8g Pa'%good_limit_Psi, '   at T = %g degK'%Tsi, '(%g degR)'%Teng_fromSI( Tsi )  )

if __name__ == "__main__":

    run_experiments( 'DowQ')



