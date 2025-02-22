#!/usr/bin/env python
# -*- coding: ascii -*-


r"""
CoolProps has a set of incompressible fluids
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
import os
from engcoolprop.ec_fluid import (CPeng_fromSI ,  CondEng_fromSI ,   DSI_fromEng ,  
                                  Deng_fromSI ,   PSI_fromEng, Peng_fromSI,
                                  PropsSI ,   SSI_fromEng ,  Seng_fromSI ,  TSI_fromEng ,  
                                  Teng_fromSI ,  UHSI_fromEng ,  UHeng_fromSI ,   Veng_fromSI )
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from engcoolprop.safe_get_property import safe_get_INCOMP_prop as safe_get_prop
from engcoolprop.find_exception_threshold import find_exception_limit

from engcoolprop.iteration_utils import find_T_from_D
from engcoolprop.iteration_utils import calc_Tnbp
from engcoolprop.utils import Same_g_len


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


# make a list of all incompressible fluids in coolprop
incomp_pure_fluidL = CP.get_global_param_string('incompressible_list_pure').split(',')
# print( 'incomp_pure_fluidL =', incomp_pure_fluidL)

class EC_Incomp_Fluid(object):
    
    fluidNameL = incomp_pure_fluidL

    def __init__(self,symbol="DowQ", T=None ,P=None, Pmax=10000.0,
                 show_warnings=2, child=1):
        '''Init generic Fluid'''

        if symbol not in self.fluidNameL:
            raise ValueError( '"%s" is NOT in coolprop incompressible list\n%s'%(symbol, repr(self.fluidNameL) ) )
        
        self.symbol = symbol
        self.Pmax = Pmax # highest pressure considered in any iterative calcs (can still input P > Pmax)
        self.Pmax_si = PSI_fromEng( Pmax )
        if P is None:
            P = int( Pmax / 10 )

        self.show_warnings = show_warnings # some calcs will issue warning statements if show_warnings >0 severe>1

        self.fluid = 'INCOMP::%s'%symbol
        self.child = child

        # get temperature limits
        self.Tmin_si =   PropsSI('Tmin','T',0,'P',0,'INCOMP::%s'%symbol)  
        self.Tmax_si =   PropsSI('Tmax','T',0,'P',0,'INCOMP::%s'%symbol) 

        self.Tmin =  Teng_fromSI( self.Tmin_si  ) + 0.0000001 # adjust slightly for round off
        self.Tmax =  Teng_fromSI( self.Tmax_si ) - 0.0000001 # adjust slightly for round off

        self.Tmid = (self.Tmin + self.Tmax) / 2.0

        self.Pmin = self.get_Psat( self.Tmin ) #+ 0.0001

        if T is None:
            T = int((self.Tmin + self.Tmax) / 2.0)
        
        self.P = P
        self.Pinput = P 

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
            # Each fluid has a different T range, so just set T to a mid point
            self.T = self.Tmin + 0.001
            print( 'NOTICE: input T too low. Changed from: %g to: %g degR'%(T, self.T) )
            print( "        INCOMP: Tmin=%g, Tmax=%g"%(self.Tmin, self.Tmax) )
        elif T>self.Tmax:
            # Each fluid has a different T range, so just set T to a mid point
            self.T = self.Tmax - 0.001
            print( 'NOTICE: input T too high. Changed from: %g to: %g degR'%(T, self.T) )
            print( "        INCOMP: Tmin=%g, Tmax=%g"%(self.Tmin, self.Tmax) )
        else:
            self.T = T

        self.calc_min_max_props()

        # set properties to input T and P
        self.setTP(self.T, self.P)
            
        if child==1: 
            self.dup = EC_Incomp_Fluid(symbol=self.symbol, T=self.T, P=self.P, child=0, Pmax=self.Pmax, 
                                       show_warnings=self.show_warnings)

    def pause_warnings(self):
        self.save_show_warnings = self.show_warnings
        self.show_warnings = 0

    def restore_warnings(self):
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
        if P >= self.Psat_max:
            return self.Tmax # 
        else:
            Psi = PSI_fromEng( P )
            def get_Dsi(T):
                Dsi = PropsSI('D','P', Psi ,'T', TSI_fromEng(T) , self.fluid)
                # print( 'Dsi =', Dsi, 'at T=', T)

            # # need to iterate on T to find Tsat            
            f = lambda T: get_Dsi( T )

            Tsat = find_exception_limit(f, tolerance=1e-4, 
                                        lower_bound=self.Tmin, upper_bound=self.Tmax,
                                        no_exception='upper_bound',  # return upper bound if no exception
                                        show_warnings=self.show_warnings )
            return Tsat


    def adjust_P_for_Psat(self, P, ind_param="T", ind_val=530.0): # psia, eng units for ind_param
        """
        If Psat is greater than P, then to maintain liquid state, use Psat
        
        Given independent variable name and value, calc T and then Psat
        NOTE: P is in Engineering units (psia)
        NOTE: ind_val is in Engineering units (degR, lbm/cuft, etc.)
        """

        if ind_param == 'T':
            Psat = self.get_Psat( ind_val ) # ind_val == T
        else:
            # to stay liquid, P will need to be above Psat
            Tsat = self.get_Tsat( P )
            Psat = self.get_Psat( Tsat )

        return max( P, Psat )# Warning shown below if Psat > P

    def setTP(self,T=530.0,P=1000.0):
        '''Calc props from T and P'''

        if T < self.Tmin:
            if self.show_warnings and self.Tmin - T > 0.01:
                print( 'T too low in setTP. Changed T=%g to T=%g'%( T, self.Tmin ))
            T = self.Tmin
        if T > self.Tmax:
            # if self.show_warnings and T - self.Tmax  > 0.01:
            if self.show_warnings and T > self.Tmax + 0.1 :
                print( 'T too high in setTP. Changed T=%g to T=%g'%( T, self.Tmax ))
            T = self.Tmax

        self.Pinput = P # save input P in case Psat changes it.
        P = self.adjust_P_for_Psat( P, ind_param="T", ind_val=T ) # ind_param="T", ind_val=530.0

        if P > self.Pinput:
            if self.show_warnings:
                print( 'P too low in setTP. Changed P=%g to Psat=%g'%( self.Pinput, P ))

        
        Tsi = TSI_fromEng( T )
        Psi = PSI_fromEng( P )
        def get_prop( prop_desc='H' ):
            try:
                prop, good_Psi = safe_get_prop( prop_desc, Psi_val=Psi, ind_name='T', ind_si_val=Tsi, 
                                                symbol=self.symbol, show_warnings=self.show_warnings>1 )
                return prop
            except:
                return float('inf')

        self.T = T
        self.P = P

        self.D = Deng_fromSI( get_prop('D') )

        self.rho = self.D / 1728.0

        self.E = UHeng_fromSI( get_prop('U') )
        self.H = UHeng_fromSI( get_prop('H') )
        self.S = Seng_fromSI( get_prop('S') )
        self.Cp = CPeng_fromSI( get_prop('C') )

        # Some INCOMP fluids have no viscosity data
        if self.Viscmin < float('inf'):
            self.Visc = Veng_fromSI( get_prop('V') ) * 1.0E5
        else:
            self.Visc = self.Viscmin

        self.Cond = CondEng_fromSI( get_prop('L') )


        self.Psat = self.get_Psat( self.T )
        self.Tsat = self.get_Tsat( self.P )

        # ============== Debug prints ================
        if self.H > self.Hmax:
            print( '...... WARNING: in setTP got H=%g > Hmax=%g ..................'%(self.H, self.Hmax))

        # print( '.......................EXITING setTP..........................')
        
    def constP_newH(self,H):
        '''Calc properties at new H with same P'''

        self.setPH( self.P, H)
        
    def constH_newP(self,P=1000.0):
        '''Calc properties at new P with same H'''

        self.setPH( P, self.H)
        
    def setPH(self,P,H):
        '''Calc properties at P and H'''

        if H < self.Hmin:
            if self.show_warnings:
                print( 'H too low in setPH. Changed H=%g to H=%g'%( H, self.Hmin ))
            H = self.Hmin
        if H > self.Hmax:
            if self.show_warnings:
                print( 'H too high in setPH. Changed H=%g to H=%g'%( H, self.Hmax ))
            H = self.Hmax


        self.Pinput = P # save input P in case Psat changes it.
        P = self.adjust_P_for_Psat( P, ind_param="H", ind_val=H )

        if P > self.Pinput:
            if self.show_warnings:
                print( 'P too low in setTP. Changed P=%g to Psat=%g'%( self.Pinput, P ))
        
        self.P = P 
        self.H = H
        
        Psi = PSI_fromEng( P )
        Hsi = UHSI_fromEng( H )

        def get_prop( prop_desc='T' ):
            try:
                prop, good_Psi = safe_get_prop( prop_desc, Psi_val=Psi, ind_name='H', ind_si_val=Hsi, 
                                                symbol=self.symbol, show_warnings=self.show_warnings>1 )
                return prop
            except:
                return float('inf')

        self.T = Teng_fromSI( get_prop('T') )
        # self.setTP( self.T, self.P)
        # return

        self.D = Deng_fromSI( get_prop('D') )
        self.rho = self.D / 1728.0

        self.E = UHeng_fromSI( get_prop('U') )
        self.S = Seng_fromSI( get_prop('S') )
        self.Cp = CPeng_fromSI( get_prop('C') )

        
        # Some INCOMP fluids have no viscosity data
        if self.Viscmin < float('inf'):
            self.Visc = Veng_fromSI( get_prop('V') ) * 1.0E5
        else:
            self.Visc = self.Viscmin

        self.Cond = CondEng_fromSI( get_prop('L') )

        self.Psat = self.get_Psat( self.T )
        self.Tsat = self.get_Tsat( self.P )
        
    def setPS(self,P,S):
        '''Calc properties at P and H'''

        if S < self.Smin:
            if self.show_warnings:
                print( 'S too low in setPS. Changed S=%g to S=%g'%( S, self.Smin ))
            S = self.Smin
        if S > self.Smax:
            if self.show_warnings:
                print( 'S too high in setPS. Changed S=%g to S=%g'%( S, self.Smax ))
            S = self.Smax

        self.Pinput = P # save input P in case Psat changes it.
        P = self.adjust_P_for_Psat( P, ind_param="S", ind_val=S )

        if P > self.Pinput:
            if self.show_warnings:
                print( 'P too low in setTP. Changed P=%g to Psat=%g'%( self.Pinput, P ))

        self.P = P 
        self.S = S
        
        Psi = PSI_fromEng( P )
        Ssi = SSI_fromEng( S )

        def get_prop( prop_desc='H' ):
            try:
                prop, good_Psi = safe_get_prop( prop_desc, Psi_val=Psi, ind_name='S', ind_si_val=Ssi, 
                                                symbol=self.symbol, show_warnings=self.show_warnings>1 )
                return prop
            except:
                return float('inf')


        self.T = Teng_fromSI( get_prop('T') )

        self.D = Deng_fromSI( get_prop('D') )
        self.rho = self.D / 1728.0

        self.H = UHeng_fromSI( get_prop('H') )
        self.E = UHeng_fromSI( get_prop('U') )
        
        self.Cp = CPeng_fromSI( get_prop('C') )

        # Some INCOMP fluids have no viscosity data
        if self.Viscmin < float('inf'):
            self.Visc = Veng_fromSI( get_prop('V') ) * 1.0E5
        else:
            self.Visc = self.Viscmin

        self.Cond = CondEng_fromSI( get_prop('L') )

        self.Psat = self.get_Psat( self.T )
        self.Tsat = self.get_Tsat( self.P )

    def constS_newP(self,P=1000.0):
        '''Calc properties at new P with same S'''

        self.setPS( P, self.S )
        
    def setPD(self,P=1000.0,D=0.01):
        '''
        Calc props from P and D
        NOTE: The pressure has NO EFFECT on incompressible density calc.
        '''

        if D < self.Dmin:
            if self.show_warnings:
                print( 'D too low in setDP. Changed D=%g to D=%g'%( D, self.Dmin ))
            D = self.Dmin
        if D > self.Dmax:
            if self.show_warnings:
                print( 'D too high in setDP. Changed D=%g to D=%g'%( D, self.Dmax ))
            D = self.Dmax


        self.Pinput = P # save input P in case Psat changes it.
        P = self.adjust_P_for_Psat( P, ind_param="D", ind_val=D )

        if P > self.Pinput:
            if self.show_warnings:
                print( 'P too low in setTP. Changed P=%g to Psat=%g'%( self.Pinput, P ))

        self.P = P 
        self.D = D

        self.T, err_flag = find_T_from_D(self, D)
        if err_flag and self.show_warnings:
            print( 'WARNING: in setPD, find_T_from_D has an error.')

        self.setTP( self.T, self.P )
        
        

    def restoreFromDup(self):
        '''restore properties from duplicate n_fluid'''
        self.T = self.dup.T
        self.P = self.dup.P
        self.D = self.dup.D
        self.rho = self.D / 1728.0
        self.E = self.dup.E
        self.H = self.dup.H
        self.S = self.dup.S
        self.Cp = self.dup.Cp
        self.Visc = self.dup.Visc
        self.Cond = self.dup.Cond
        self.Psat = self.dup.Psat

    def saveToDup(self):
        '''save properties to duplicate n_fluid'''
        self.dup.T = self.T
        self.dup.P = self.P
        self.dup.D = self.D
        self.dup.E = self.E
        self.dup.H = self.H
        self.dup.S = self.S
        self.dup.Cp = self.Cp
        self.dup.Visc = self.Visc
        self.dup.Cond = self.Cond
        self.dup.Psat = self.Psat

    def initFromObj(self, obj):
        '''initialize properties from another n_fluid object'''
        if  self.symbol == obj.symbol:
            
            self.T = obj.T
            self.P = obj.P
            self.D = obj.D
            self.rho = self.D / 1728.0
            self.E = obj.E
            self.H = obj.H
            self.S = obj.S
            self.Cp = obj.Cp
            self.Visc = obj.Visc
            self.Cond = obj.Cond
            self.Psat = obj.Psat
        else:
            raise Exception('Wrong fluid for initializing')

    def getStrTransport(self):
        '''create a string from the Transport properties'''
        return  "%s Cp=%6g Visc=%6g ThCond=%6g" %\
        (self.symbol,self.Cp, self.Visc, self.Cond)

    def getStrTPD(self):
        '''create a string from the TPDEHS properties'''
        return  "%s T=%6.1f P=%6.1f D=%.4f E=%6.2f H=%6.2f S=%.3f" %\
        (self.symbol,self.T, self.P, self.D, self.E, self.H, self.S)

    def getStrTPDphase(self):
        '''create a string from the TPDEHS properties'''
        return  "%s T=%6.1f P=%6.1f D=%.4f E=%6.2f H=%6.2f S=%.3f" %\
        (self.symbol,self.T, self.P, self.D, self.E, self.H, self.S)

    def printTPD(self):
        '''print a string from the TPDEHS properties'''
        print(self.getStrTPD())


    def printProps(self):
        '''print a multiline property summary with units'''

        # get formatted floats that are same length to help the look of table
        SGL = Same_g_len(self, ['Cond', "Cp", 'D', 'E', 'H', 'P', 'S', 'T', 'Visc', 'rho', 'Tnbp', 'Psat', 'Tsat'] )


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
            print("V =%s"%SGL.Visc," viscosity [1.0E5 * lbm/ft-sec]", '    Range(%8g - %8g)'%(self.Viscmin, self.Viscmax) )
        else:
            print("V =UNDEFINED","viscosity [1.0E5 * lbm/ft-sec]" )

        print("C =%s"%SGL.Cond," thermal conductivity [BTU/ft-hr-R]", 'Range(%8g - %8g)'%(self.Condmin, self.Condmax) )

        if self.Tnbp is not None:
            print("    Tnbp =%s"%SGL.Tnbp," degR,")
        print("    rho  =%s"%SGL.rho," lbm/cuin",   '                   Range(%.6f - %.6f) lbm/cuin'%(self.rho_min, self.rho_max))
        
        if self.Psat is not None and self.Psat_max > 0:
            print("    Psat =%s"%SGL.Psat," psia", '                       Range(%g - %g) psia'%(self.Psat_min, self.Psat_max))
            print("    Tsat =%s"%SGL.Tsat," degR", '                       Range(%g - %g) degR'%(self.Tsat_min, self.Tsat_max))

    def calc_min_max_props(self, do_print=False):
        # print( '.......................Entered calc_min_max_props ..........................')

        self.Psat_min = self.get_Psat( self.Tmin )
        self.Psat_max = self.get_Psat( self.Tmax )

        self.Tsat_min = self.get_Tsat( self.Pmin )
        self.Tsat_max = self.get_Tsat( self.Pmax )

        def get_prop( T, P, prop_desc='H' ):
            Tsi = TSI_fromEng( T )
            Psi = PSI_fromEng( P )
            try:
                prop, good_Psi = safe_get_prop( prop_desc, Psi_val=Psi, ind_name='T', ind_si_val=Tsi, 
                                                symbol=self.symbol, show_warnings=self.show_warnings>1 )
                return prop
            except:
                return float('inf')

        Dmin_si = get_prop(self.Tmax, self.Psat_max,'D')
        Dmax_si =get_prop(self.Tmin, self.Pmax,'D')

        self.Dmin = Deng_fromSI( Dmin_si )
        self.Dmax = Deng_fromSI( Dmax_si )

        self.rho_min = self.Dmin / 1728.0
        self.rho_max = self.Dmax / 1728.0

        self.Emin = UHeng_fromSI( get_prop(self.Tmin, self.Pmax,'U') )
        self.Emax = UHeng_fromSI(  get_prop(self.Tmax, self.Psat_max, 'U') )

        self.Hmin = UHeng_fromSI( get_prop(self.Tmin, self.Psat_min, 'H') )
        self.Hmax = UHeng_fromSI( get_prop(self.Tmax, self.Pmax, 'H') )
        
        self.Smin = Seng_fromSI( get_prop(self.Tmin, self.Pmax, 'S') )
        self.Smax = Seng_fromSI( get_prop(self.Tmax, self.Psat_max, 'S') )
        
        self.Cpmin = CPeng_fromSI( get_prop(self.Tmin, self.Pmax, 'C') )
        self.Cpmax = CPeng_fromSI( get_prop(self.Tmax, self.Pmax, 'C') )

        self.pause_warnings()
        self.Viscmin = Veng_fromSI( get_prop(self.Tmax, self.Pmax, 'V') ) * 1.0E5
        self.Viscmax = Veng_fromSI( get_prop(self.Tmin, self.Pmax, 'V') ) * 1.0E5
        self.restore_warnings()

        if self.Viscmin is None:
            self.Viscmin = float('inf')
        if self.Viscmax is None:
            self.Viscmax = float('inf')
        
        self.Condmin = CondEng_fromSI( get_prop(self.Tmin, self.Pmax, 'L') )
        self.Condmax = CondEng_fromSI( get_prop(self.Tmax, self.Pmax, 'L') )
        
        if do_print:
            minmaxL = ['Cond', "Cp", 'D', 'E', 'H', 'P', 'S', 'T', 'Visc']

            for name in minmaxL:
                max_name = name + 'max'
                print( '%10s'%max_name, '%9g'%getattr(self, max_name) )
                min_name = name + 'min'
                print( '%10s'%min_name, '%9g'%getattr(self, min_name) )
        # print( '.......................Finished calc_min_max_props ..........................')


if __name__ == '__main__':
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
    C = EC_Incomp_Fluid( symbol=symbol, T=None, P=0 )
    
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
