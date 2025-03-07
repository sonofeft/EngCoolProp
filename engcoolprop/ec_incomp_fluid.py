#!/usr/bin/env python
# -*- coding: ascii -*-


r"""
CoolProp has a set of incompressible fluids
EC_Incomp_Fluid is a wrapper of those fluids using Engineering Units.

EngCoolProp uses units of primarily inch, lbm, lbf, sec, BTU (some use of ft and hour).::

    #:c     The following are the default units for each property of Incompressible Fluids.
    #:c 
    #:c     T = Temperature = degR
    #:c     P = Pressure = psia
    #:c     D = Density = lbm/cu ft
    #:c     rho = Density = lbm/cu inch
    #:c     E = Internal Energy = BTU/lbm
    #:c     H = Enthalpy = BTU/lbm
    #:c     S = Entropy = BTU/lbm degR
    #:c     Cp = Heat Capacity (const. P) = BTU/lbm degR
    #:c     V = Viscosity = 1.0E5 * lbm/ft-sec
    #:c     C = Thermal Conductivity = BTU/ft-hr-R

"""
import traceback
from engcoolprop.ec_fluid import (CPeng_fromSI ,  CondEng_fromSI ,   DSI_fromEng ,  
                                  Deng_fromSI ,   PSI_fromEng, Peng_fromSI,
                                  PropsSI ,   SSI_fromEng ,  Seng_fromSI ,  TSI_fromEng ,  
                                  Teng_fromSI ,  UHSI_fromEng ,  UHeng_fromSI ,   Veng_fromSI )
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from engcoolprop.safe_get_property import safe_get_INCOMP_prop
from engcoolprop.find_exception_threshold import find_exception_limit

from engcoolprop.InterpProp_scipy import get_density_interpolator
from engcoolprop.iteration_utils import find_T_from_Dterp, find_T_at_P
from engcoolprop.iteration_utils import calc_Tnbp
from engcoolprop.utils import Same_g_len, print_avoid_valerr, incomp_pure_fluidL
from engcoolprop.banner import banner


# create simple look-up that is order-independent for input pairs

call_tuplesD = {}  # key:tuple indep vars (e.g. ("T","P")), value: dependent property (e.g. "S", "H", etc.)
call_tuplesD[('D', 'P')] = {'L', 'T', 'S', 'U', 'H', 'V', 'C'}
call_tuplesD[('H', 'P')] = {'L', 'T', 'S', 'U', 'V', 'C', 'D'}
call_tuplesD[('P', 'D')] = {'L', 'T', 'S', 'U', 'H', 'V', 'C'}
call_tuplesD[('P', 'H')] = {'L', 'T', 'S', 'U', 'V', 'C', 'D'}
call_tuplesD[('P', 'S')] = {'L', 'T', 'U', 'H', 'V', 'C', 'D'}
call_tuplesD[('P', 'T')] = {'L', 'S', 'U', 'H', 'V', 'C', 'D'}
call_tuplesD[('S', 'P')] = {'L', 'T', 'U', 'H', 'V', 'C', 'D'}
call_tuplesD[('T', 'P')] = {'L', 'S', 'U', 'H', 'V', 'C', 'D'}          


# how to make a list of all incompressible fluids in coolprop
# incomp_pure_fluidL = CP.get_global_param_string('incompressible_list_pure').split(',')

class EC_Incomp_Fluid(object):
    
    def __init__(self,symbol="DowQ", T=None ,P=None, Pmax=10000.0,
                 show_warnings=2, auto_fix_value_errors=True):
        """        Init generic Incompressible Fluid

        Args:
            symbol (str): Name and percent mass. (e.g. "Water", "DowQ")
            T (None or float): Temperature degR (None sets T to (Tmin+Tmax)/2)
            P (None or float): Pressure psia. (None set P to Pmax/10)
            Pmax (float): Max expected pressure psia. Defaults to 10000.0.
            show_warnings (int): Sets warning level(0=None, 1=Only serious, 2=all). Defaults to 2.
            auto_fix_value_errors (bool): Action when ValueError occurs. 
                                        True=correct problem with a warning, 
                                        False throws Exception Defaults to False.
        """

        if symbol not in incomp_pure_fluidL:
            raise ValueError( '"%s" is NOT in coolprop incompressible list\n%s'%(symbol, repr(incomp_pure_fluidL) ) )
        
        self.symbol = symbol
        self.Pmax = Pmax # highest pressure considered in any iterative calcs (can still input P > Pmax)
        self.Pmax_si = PSI_fromEng( Pmax )
        if P is None:
            P = int( Pmax / 10 )

        self.show_warnings = show_warnings # some calcs will issue warning statements if show_warnings >0 severe>1
        self.auto_fix_value_errors = auto_fix_value_errors

        if show_warnings>1 and auto_fix_value_errors:
            banner( 'NOTICE: any input violations on limits of T, D, H or S\nwill be automatically corrected (set to min or max).'+\
                   '\nTo change this behavior set "auto_fix_value_errors" to False'+\
                    '\nTo suppress this banner set "show_warnings" to 0 or 1')

        self.fluid = 'INCOMP::%s'%symbol

        # get temperature limits (Note attempt to avoid CoolProp round-off problems)
        self.Tmin_si =   PropsSI('Tmin','T',0,'P',0,'INCOMP::%s'%symbol) * 1.0000000000000002
        self.Tmax_si =   PropsSI('Tmax','T',0,'P',0,'INCOMP::%s'%symbol) * 0.9999999999999999

        # ==== some fluids have min or max point in density curve ( screws up state = f(P, D) ) ====
        if self.symbol == 'PMR':
            self.Tmax_si = 500.0 # degK # has minimum in D=f(T) curve.
            if show_warnings > 0:
                print( 'WARNING: PMR Tmax reduced to 500 degK (900 degR) due to error in CoolProp equation for Density' )
        elif self.symbol == 'NBS':
            self.Tmin_si = 275.0 # degK # has maximum in D=f(T) curve.
            if show_warnings > 0:
                print( 'WARNING: NBS Tmin increased to 275 degK (495 degR) due to error in CoolProp equation for Density' )
        elif self.symbol == 'FoodWater':
            self.Tmin_si =  280.0 # degK # has maximum in D=f(T) curve.
            if show_warnings > 0:
                print( 'WARNING: FoodWater Tmin increased to 280 degK (504 degR) due to error in CoolProp equation for Density' )
        elif self.symbol == 'LiqNa':
            self.Tmax_si = 2455.0 # degK # has maximum in H=f(T,P) curve.
            if show_warnings > 0:
                print( 'WARNING: LiqNa Tmax reduced to 2455 degK (4419 degR) due to error in CoolProp equation for Enthalpy' )



        self.Tmin =  Teng_fromSI( self.Tmin_si  )
        self.Tmax =  Teng_fromSI( self.Tmax_si )

        self.Tmid = (self.Tmin + self.Tmax) / 2.0

        self.Dterp = get_density_interpolator( self )

        self.check_visc_support()

        # Note that Psat_min and Psat_max are calculated in calc_min_max_props
        self.Pmin = self.get_Psat( self.Tmin ) #+ 0.0001

        if T is None:
            T = int((self.Tmin + self.Tmax) / 2.0)
            # else: T = input T
        
        self.P = P
        self.Pinput = P # in setTP, if P is below Psat, P will be increased to Psat

        try:
            Tnbp, error_code = calc_Tnbp( self )
            if not error_code:
                self.Tnbp = Tnbp 
            else:
                self.Tnbp = None 
        except:
            self.Tnbp = None 

        # if input T is in range, use it... otherwise
        if T<self.Tmin:
            if self.auto_fix_value_errors:
                self.T = self.Tmin + 0.001
                print( 'NOTICE: input T too low. Changed from: %g to: %g degR'%(T, self.T) )
                print( "        INCOMP: Tmin=%g, Tmax=%g"%(self.Tmin, self.Tmax) )
            else:
                print_avoid_valerr()
                raise ValueError("Input T too low. T=%g,  Tmin=%g"%(T, self.Tmin))

        elif T>self.Tmax:
            if self.auto_fix_value_errors:                            
                self.T = self.Tmax - 0.001
                print( 'NOTICE: input T too high. Changed from: %g to: %g degR'%(T, self.T) )
                print( "        INCOMP: Tmin=%g, Tmax=%g"%(self.Tmin, self.Tmax) )
            else:
                print_avoid_valerr()
                raise ValueError("Input T too high. T=%g,  Tmax=%g"%(T, self.Tmax))
        else:
            self.T = T

        self.calc_min_max_props()

        # set properties to input T and P
        # print( 'First call to setTP')
        self.setTP(self.T, self.P)
        # print( 'Back from: First call to setTP')
            
            
    def check_visc_support(self):
        try:
            _ =  PropsSI('V','T',self.Tmin_si,'P',self.Pmax_si, self.fluid ) # Vmax at Tmin
            _ =  PropsSI('V','T',self.Tmax_si,'P',self.Pmax_si, self.fluid ) # Vmin at Tmax
            self.visc_is_supported = True
        except:
            self.visc_is_supported = False
            if self.show_warnings:
                print( 'check_visc_support: Visc NOT supported for fluid:', self.symbol)

        # print( 'Found visc_is_supported =', self.visc_is_supported)



    def set_warnings(self, show_warnings=2):
        """Set the show_warnings flag"""
        if show_warnings is True:
            show_warnings = 1
        elif show_warnings is False:
            show_warnings = 0

        self.save_show_warnings = self.show_warnings
        self.show_warnings = show_warnings

    def pause_warnings(self):
        """Pause all warnings until self.show_warnings is changed"""
        self.save_show_warnings = self.show_warnings
        self.show_warnings = 0

    def restore_warnings(self):
        """Restore the self.show_warnings flag to the last pause_warnings condition"""
        try:
            self.show_warnings = self.save_show_warnings
        except:
            pass
    

    def setProps(self, **inpD):
        '''Generic call using any P with supported inputs T,D,S,H::
        
        #: for example
        #: ec.setProps(T=100, P=200)
        #: ec.setProps(D=0.1, P=100)
        #: ec.setProps(P=100, H=20)
        #:
        #: Made so both:
        #: setProps(T=530, P=100)  AND  
        #: setProps(P=100, T=530) will work

        '''
        
        key = tuple( inpD.keys() )
        if key not in call_tuplesD:
            print('ERROR... setProps called with unknown inputs', inpD)
            print('Allowable inputs are: T,P,D,H,S')
            if len(inpD) != 2:
                print('ONLY 2 inputs are allowed.  (e.g. T=530.0,P=100.0)')
            raise ValueError('setProps called with illegal inputs: %s'%str(inpD))
        
        if 'P' not in inpD:
            raise ValueError('"P" must be one of the setProps inputs.')
        
        if len(inpD) != 2:
            raise ValueError('ONLY 2 inputs are allowed.  (e.g. T=530.0,P=100.0), not' + repr(inpD))

        if "T" in inpD:
            self.setTP( ** inpD )
        elif "D" in inpD:
            self.setPD( **inpD )
        elif "H" in inpD:
            self.setPH( **inpD )
        elif "S" in inpD:
            self.setPS( **inpD )
        else:
            raise ValueError( 'Invalid input to setProps ' + repr(inpD) )
        
    def get_D_at_T(self, T): # degR
        Psi = PSI_fromEng( 5000 )
        
        Tsi = TSI_fromEng( T )
        
        D = Deng_fromSI( PropsSI('D', 'T',Tsi,'P',Psi, self.fluid) )
        return D # lbm/cuft


    def get_Psat(self, T): # degR
        """
        Given T in degR (Engineering units) return Psat in psia (if supported)
        NOTE: if not supported, return 0 for Psat
        """
        try:
            Psat = Peng_fromSI( PropsSI('P','T', TSI_fromEng( T ) ,'Q',0, self.fluid) )
        except:
            Psat = 0.0
        return Psat # psia


    def get_Tsat(self, P ): # psia
        """
        Calculate Tsat at given P (psia).
        If P is greater than max Psat, simply return Tmax.
        """
        raise Exception( 'INCOMP Psat data not complete enough for Tsat Calculation')


    def setTP(self,T=530.0,P=1000.0):
        '''Calc props from T and P'''


        if T < self.Tmin:
            if self.auto_fix_value_errors:
                if self.show_warnings and self.Tmin - T > 0.01:
                    print( 'T too low in setTP. Changed T=%g to Tmin=%g'%( T, self.Tmin ))
                T = self.Tmin
            else:
                print_avoid_valerr()
                raise ValueError("Input T too low in setTP. T=%g,  Tmin=%g"%(T, self.Tmin))
        if T > self.Tmax:
            if self.auto_fix_value_errors:
                # if self.show_warnings and T - self.Tmax  > 0.01:
                if self.show_warnings and T > self.Tmax + 0.1 :
                    print( 'T too high in setTP. Changed T=%g to Tmax=%g'%( T, self.Tmax ))
                T = self.Tmax
            else:
                print_avoid_valerr()
                raise ValueError("Input T too high in setTP. T=%g,  Tmin=%g"%(T, self.Tmax))

        self.Pinput = P # save input P in case Psat changes it.

        # === Don't allow Vapor phase... increase P to phase line ===
        # P = max( P,  self.get_Psat( T ))

        self.Psat = self.get_Psat( self.T )

        # if P > self.Pinput:
        if self.Psat > P:
            # print( 'Psat - P =', self.Psat - P)
            if self.auto_fix_value_errors:
                P = self.Psat
                if self.show_warnings:
                    print( 'P too low in setTP. Changed P=%g to Psat=%g'%( self.Pinput, self.Psat ))
            else:
                print_avoid_valerr()
                raise ValueError( 'P too low in setTP. P=%g should be >= Psat=%g'%( self.Pinput, self.Psat ))


        
        Tsi = TSI_fromEng( T )
        Psi = PSI_fromEng( P )
        def get_prop( prop_desc='H' ):
            try:
                prop, good_Psi = safe_get_INCOMP_prop( prop_desc, Psi_val=Psi, ind_name='T', ind_si_val=Tsi, 
                                                       symbol=self.symbol, show_warnings=self.show_warnings>1,
                                                       Pmax=self.Pmax )
                return prop
            except:
                return float('inf')

        self.T = T
        self.P = P

        # self.D = Deng_fromSI( get_prop('D') )
        self.D = self.get_D_at_T( self.T )

        self.rho = self.D / 1728.0

        self.E = UHeng_fromSI( get_prop('U') )
        self.H = UHeng_fromSI( get_prop('H') )
        self.S = Seng_fromSI( get_prop('S') )
        self.Cp = CPeng_fromSI( get_prop('C') )

        # Some INCOMP fluids have no viscosity data
        if self.visc_is_supported:
            self.Visc = Veng_fromSI( get_prop('V') ) # multiply only for display * 1.0E5
        else:
            self.Visc = float('inf')

        self.Cond = CondEng_fromSI( get_prop('L') )


        # self.Tsat = self.get_Tsat( self.P )

        # ============== Debug prints ================
        # if self.H > self.Hmax:
        #     print( '...... WARNING: in setTP got H=%g > Hmax=%g ..................'%(self.H, self.Hmax))

        # print( '.......................EXITING setTP..........................')
        
    def constP_newH(self,H):
        '''Calc properties at new H with same P'''

        self.setPH( self.P, H)
        
    def constH_newP(self,P=1000.0):
        '''Calc properties at new P with same H'''

        self.setPH( P, self.H)
        
    def setPH(self,P,H):
        '''Calc properties at P and H'''

        if self.symbol == 'Air':
            raise Exception( 'PH data for Air will not work properly in setPH.' )
  
        if H < self.Hmin:
            if self.auto_fix_value_errors:
                if self.show_warnings:
                    print( 'H too low in setPH. Changed H=%g to H=%g'%( H, self.Hmin ))
                H = self.Hmin
            else:
                print_avoid_valerr()
                raise ValueError("Input H too low in setPH. H=%g,  Hmin=%g"%(H, self.Hmin))
        if H > self.Hmax:
            if self.auto_fix_value_errors:
                if self.show_warnings:
                    print( 'H too high in setPH. Changed H=%g to Hmax=%g'%( H, self.Hmax ))
                H = self.Hmax
            else:
                print_avoid_valerr()
                raise ValueError("Input H too high in setPH. H=%g,  Hmax=%g"%(H, self.Hmax))


        self.Pinput = P # save input P in case Psat changes it.
        # P = self.adjust_P_for_Psat( P, ind_param="H", ind_val=H )
        # if P > self.Pinput:
        #     if self.show_warnings:
        #         print( 'P too low in setPH. Changed P=%g to Psat=%g'%( self.Pinput, P ))
        
        self.P = P 
        self.H = H
        
        Psi = PSI_fromEng( P )
        Hsi = UHSI_fromEng( H )

        # This approach worked for PS, but not PH
        # self.T, err_flag = find_T_at_P(self, P, dep_name='H', dep_val=H, tol=1.0E-12)
        # self.setTP( self.T, self.P)
        # return

        def get_prop( prop_desc='T' ):
            try:
                prop, good_Psi = safe_get_INCOMP_prop( prop_desc, Psi_val=Psi, ind_name='H', ind_si_val=Hsi, 
                                                       symbol=self.symbol, show_warnings=self.show_warnings>1,
                                                       Pmax=self.Pmax )
                return prop
            except:
                return float('inf')

        self.T = Teng_fromSI( get_prop('T') )

        
        self.setTP( self.T, self.P)
        
    def setPS(self,P,S):
        '''Calc properties at P and H'''

        if self.symbol == 'Air':
            raise Exception( 'PS data for Air will not work properly in setPS.' )

        if S < self.Smin:
            if self.auto_fix_value_errors:
                if self.show_warnings:
                    print( 'S too low in setPS. Changed S=%g to S=%g'%( S, self.Smin ))
                S = self.Smin
            else:
                print_avoid_valerr()
                raise ValueError("Input D too low in setPS. S=%g,  Smin=%g"%(S, self.Smin))
            
        if S > self.Smax:
            if self.auto_fix_value_errors:
                if self.show_warnings:
                    print( 'S too high in setPS. Changed S=%g to Smax=%g'%( S, self.Smax ))
                S = self.Smax
            else:
                print_avoid_valerr()
                raise ValueError("Input D too low in setPS. S=%g,  Smin=%g"%(S, self.Smin))

        self.Pinput = P # save input P in case Psat changes it.
        # P = self.adjust_P_for_Psat( P, ind_param="S", ind_val=S )
        # if P > self.Pinput:
        #     if self.show_warnings:
        #         print( 'P too low in setPS. Changed P=%g to Psat=%g'%( self.Pinput, P ))

        self.P = P 
        self.S = S

        self.T, err_flag = find_T_at_P(self, P, dep_name='S', dep_val=S, tol=1.0E-12)
        self.setTP( self.T, self.P)

    def constS_newP(self,P=1000.0):
        '''Calc properties at new P with same S'''

        self.setPS( P, self.S )
        
    def setPD(self,P=1000.0,D=0.01):
        '''
        Calc props from P and D
        NOTE: The pressure has NO EFFECT on calculated temperature for incompressible density.
        '''

        if D < self.Dmin:
            if self.auto_fix_value_errors:
                if self.show_warnings:
                    print( 'D too low in setDP. Changed D=%g to D=%g'%( D, self.Dmin ))
                D = self.Dmin
            else:
                print_avoid_valerr()
                raise ValueError("Input D too low in setPD. D=%g,  Dmin=%g"%(D, self.Dmin))

        if D > self.Dmax:
            if self.auto_fix_value_errors:
                if self.show_warnings:
                    print( 'D too high in setDP. Changed D=%g to Dmax=%g'%( D, self.Dmax ))
                D = self.Dmax
            else:
                print_avoid_valerr()
                raise ValueError("Input D too high in setPD. D=%g,  Dmax=%g"%(D, self.Dmax))

        self.Pinput = P # save input P in case Psat changes it.
        # P = self.adjust_P_for_Psat( P, ind_param="D", ind_val=D )

        if P > self.Pinput:
            if self.show_warnings:
                print( 'P too low in setPD. Changed P=%g to Psat=%g'%( self.Pinput, P ))

        self.P = P 
        self.D = D

        self.T, err_flag = find_T_from_Dterp(self, D)
        if err_flag and self.show_warnings:
            print( 'WARNING: in setPD, find_T_from_Dterp has an error.')

        self.setTP( self.T, self.P )
        
        
    def getStrTransport(self):
        '''create a string from the Transport properties'''
        Visc = self.Visc * 1.0E5
        return  "%s Cp=%6g Visc=%6g ThCond=%6g" %\
        (self.fluid,self.Cp, Visc, self.Cond)

    def getStrTPD(self):
        '''create a string from the TPDEHS properties'''
        return  "%s T=%6.1f P=%6.1f D=%.4f E=%6.2f H=%6.2f S=%.3f" %\
        (self.fluid,self.T, self.P, self.D, self.E, self.H, self.S)

    # def getStrTPDphase(self):
    #     '''create a string from the TPDEHS properties'''
    #     return  "%s T=%6.1f P=%6.1f D=%.4f E=%6.2f H=%6.2f S=%.3f" %\
    #     (self.fluid,self.T, self.P, self.D, self.E, self.H, self.S)

    def printTPD(self):
        '''print a string from the TPDEHS properties'''
        print(self.getStrTPD())

    def printTransport(self):
        '''print a string of Transport properties'''
        print(self.getStrTransport())

    def printProps(self):
        '''print a multiline property summary with units'''

        # get formatted floats that are same length to help the look of table
        SGL = Same_g_len(self, ['Cond', "Cp", 'D', 'E', 'H', 'P', 'S', 'T', 'Visc', 'rho', 'Tnbp', 'Psat'] ) # , 'Tsat'


        print("State Point for fluid",self.fluid,"("+ self.symbol +")")
        print("T =%s"%SGL.T," degR,",    '                             Range(%8g - %8g) degR'%(self.Tmin, self.Tmax))

        print("P =%s"%SGL.P," psia",    '                              Range(%8g - %8g) psia'%(self.Pmin, self.Pmax))
        if self.Pinput != self.P:
            print("    Pinput =%8g"%self.Pinput, '(Adjusted due to Psat)')

        print("D =%s"%SGL.D," lbm/cuft",    '                          Range(%8g - %8g) lbm/cuft'%(self.Dmin, self.Dmax))

        print("E =%s"%SGL.E," BTU/lbm",    '                           Range(%8g - %8g) BTU/lbm'%(self.Emin, self.Emax))
        print("H =%s"%SGL.H," BTU/lbm",    '                           Range(%8g - %8g) BTU/lbm'%(self.Hmin, self.Hmax))
        print("S =%s"%SGL.S," BTU/lbm degR",    '                      Range(%.6f - %.6f) BTU/lbm degR'%(self.Smin, self.Smax))
        print("Cp=%s"%SGL.Cp," BTU/lbm degR",   '                      Range(%8g - %8g) BTU/lbm degR'%(self.Cpmin, self.Cpmax))

        if self.Visc < float('inf'):
            print("V =%s"%SGL.Visc," viscosity [1.0E5 * lbm/ft-sec]", 
                  '    Range(%8g - %8g)'%(self.Viscmin * 1.0E5, self.Viscmax * 1.0E5) )
        else:
            print("V =UNDEFINED","viscosity [1.0E5 * lbm/ft-sec]" )

        print("C =%s"%SGL.Cond," thermal conductivity [BTU/ft-hr-R]", 'Range(%8g - %8g)'%(self.Condmin, self.Condmax) )

        if self.Tnbp is not None:
            print("    Tnbp =%s"%SGL.Tnbp," degR,")
        print("    rho  =%s"%SGL.rho," lbm/cuin",   '                   Range(%.6f - %.6f) lbm/cuin'%(self.rho_min, self.rho_max))
        
        if self.psat_is_supported:
            print("    Psat =%s"%SGL.Psat," psia", '                       Range(%g - %g) psia'%(self.Psat_min, self.Psat_max))
            # print("    Tsat =%s"%SGL.Tsat," degR", '                       Range(%g - %g) degR'%(self.Tsat_min, self.Tsat_max))

    def calc_min_max_props(self, do_print=False):
        # print( '.......................Entered calc_min_max_props ..........................')
        T_save = self.T
        P_save = self.P

        self.setTP(self.T, self.P)

        self.Psat_min = self.get_Psat( self.Tmin )
        self.Psat_max = self.get_Psat( self.Tmax )
        if self.Psat_max > 0:
            self.psat_is_supported = True
        else:
            self.psat_is_supported = False

        # self.Tsat_min = self.get_Tsat( self.Pmin )
        # self.Tsat_max = self.get_Tsat( self.Pmax )


        # make a list of the four corners of the T,P space
        tpL = [(self.Tmin, 0), (self.Tmin, self.Pmax), (self.Tmax, 0), (self.Tmax, self.Pmax)]

        resultD = {} # key:min/max name (e.g. Dmin, Cpmax), value:float value

        # set starting values for min/max of these properties
        for prop_name in [ 'D', 'E', 'H', 'S', 'Cp', 'Visc', 'Cond' ]:
            resultD[prop_name+'min'] = float('inf')
            resultD[prop_name+'max'] = float('-inf')

        self.pause_warnings()
        for T,P in tpL:
            self.setTP( T, P )

            for prop_name in [ 'D', 'E', 'H', 'S', 'Cp', 'Visc', 'Cond' ]:
                if getattr(self, prop_name) < resultD[prop_name+'min']:
                    resultD[prop_name+'min'] = getattr(self, prop_name)

                if getattr(self, prop_name) > resultD[prop_name+'max']:
                    resultD[prop_name+'max'] = getattr(self, prop_name)
        self.restore_warnings()

        # set min/max values found in resultD
        for prop_name, value in resultD.items():
            # print( 'setting ', prop_name, value)
            setattr( self, prop_name, value )

        self.rho_min = self.Dmin / 1728.0
        self.rho_max = self.Dmax / 1728.0
        
        self.Dmin_si = DSI_fromEng( self.Dmin )
        self.Dmax_si =  DSI_fromEng( self.Dmax )

        if do_print:
            minmaxL = ['Cond', "Cp", 'D', 'E', 'H', 'P', 'S', 'T', 'Visc']

            for name in minmaxL:
                max_name = name + 'max'
                print( '%10s'%max_name, '%9g'%getattr(self, max_name) )
                min_name = name + 'min'
                print( '%10s'%min_name, '%9g'%getattr(self, min_name) )
        # print( '.......................Finished calc_min_max_props ..........................')

        # Restore state to T and P from calling method
        self.setTP( T_save, P_save)        



def dev_tests():
    """
    __incompressibles_pure__ = ['Acetone', 'Air', 'AS10', 'AS20', 'AS30', 'AS40', 'AS55', 'DEB', 
    'DowJ', 'DowJ2', 'DowQ', 'DowQ2', 'DSF', 'Ethanol', 'ExampleDigitalPure', 'ExamplePure',
    'FoodAsh', 'FoodCarbohydrate', 'FoodFat', 'FoodFiber', 'FoodIce', 'FoodProtein', 'FoodWater', 
    'HC10', 'HC20', 'HC30', 'HC40', 'HC50', 'HCB', 'HCM', 'Hexane', 'HFE', 'HFE2', 'HY20', 'HY30',
    'HY40', 'HY45', 'HY50', 'LiqNa', 'NaK', 'NBS', 'PBB', 'PCL', 'PCR', 'PGLT', 'PHE', 'PHR', 
    'PLR', 'PMR', 'PMS1', 'PMS2', 'PNF', 'PNF2', 'S800', 'SAB', 'T66', 'T72', 'TCO', 'TD12', 
    'TVP1', 'TVP1869', 'TX22', 'TY10', 'TY15', 'TY20', 'TY24', 'Water', 'XLT', 'XLT2', 'ZS10', 
    'ZS25', 'ZS40', 'ZS45', 'ZS55']
    """

    symbol = 'Water'
    
    print( '='*22, "%s at Pressure = 0"%symbol, '='*22 )
    C = EC_Incomp_Fluid( symbol=symbol, T=None, P=0, auto_fix_value_errors=True )
    
    # C.setTP(T= (C.Tmin+C.Tmax)/2.0, P=100) # this temperature throws an exception ???
    
    C.printProps()

    if C.Tnbp is not None:
        print( '='*22, "%s at Tnbp"%symbol, 'Tnbp=%g degR'%C.Tnbp, '='*22 )
        C.setTP( C.Tnbp, 15.0)
        C.printProps()
    
        print( '='*22, "%s short prints at Tnbp"%symbol, '='*22 )
        print(C.getStrTransport())
        C.printTPD()



    print()
    
    print( '='*22, "%s over full Temp Range"%symbol, '='*22 )
    t_range = C.Tmax - C.Tmin
    N = 9
    tL = [C.Tmin + (i/N)*t_range for i in range(N+1)]
    print( ['%.1f'%T for T in tL])
    for T in tL:
        try:
            C.setTP(T, 5000.0)
            C.printTPD()
        except:
            print( 'T=%.1f FAILED'%T )
            Tsi = TSI_fromEng( T )
            Psi = PSI_fromEng( 50000 )

            Dsi = PropsSI('D','T',Tsi,'P',Psi,'INCOMP::%s'%symbol)
            print( '           Dsi =', Dsi)
        
    print()
    print( '='*22, "%s test constP_newH"%symbol, '='*22 )
    C.setTP(C.Tmid, 500.0)
    C.printTPD()
    C.constP_newH( C.H * 0.9 )
    C.printTPD()


    print()
    print( '='*22, "%s test constH_newP"%symbol, '='*22 )
    C.setTP(C.Tmid, 500.0)
    C.printTPD()
    C.constH_newP( 530.0 )
    C.printTPD()
    
    print()
    print( '='*22, "%s test setPH"%symbol, '='*22 )
    C.setTP(C.Tmid, 500.0)
    C.printTPD()
    C.setPH( 500., C.H * 0.95 )
    C.printTPD()
    
    print()
    print( '='*22, "%s test constS_newP"%symbol, '='*22 )
    C.setTP(C.Tmid, 500.0)
    C.printTPD()
    C.constS_newP(P=1000.0)
    C.printTPD()
    
    print()
    print( '='*22, "%s test setPD"%symbol, '='*22 )
    C.setTP(C.Tmid, 500.0)
    C.printTPD()
    C.setPD(P=800.0,D=C.D)
    C.printTPD()

    print()
    print( '='*22, "%s test calc_min_max_props"%symbol, '='*22 )
    C.calc_min_max_props( do_print=True )    
    print()
    print( '='*22, "print full %s properties"%symbol, '='*22 )
    C.printProps()
    print()

    banner( 'This should fail due to "auto_fix_value_errors = False"' )
    print( '='*22, "print full %s properties"%symbol, '='*22 )
    try:
        C_bad = EC_Incomp_Fluid( symbol='Water', T=None, P=0 )
    except:
        tb_str = traceback.format_exc()
        print( tb_str.split('raise ValueError')[-1])
        # print( tb_str )
    


if __name__ == '__main__':
    dev_tests()
