#!/usr/bin/env python
# -*- coding: ascii -*-

r"""
CoolProps has a set of incompressible fluids
EC_Incomp_Fluid is a wrapper of those fluids using Engineering Units.

EngCoolProp uses units of primarily inch, lbm, lbf, sec, BTU (some use of ft and hour).

    The following are the default units for each property of Incompressible Fluids.

    T = Temperature = degR
    P = Pressure = psia
    D = Density = lbm/cu ft
    rho = Density = lbm/cu inch
    E = Internal Energy = BTU/lbm
    H = Enthalpy = BTU/lbm
    S = Entropy = BTU/lbm degR
    Cp = Heat Capacity (const. P) = BTU/lbm degR
    V = Viscosity = 1.0E5 * lbm/ft-sec
    C = Thermal Conductivity = BTU/ft-hr-R

"""
import os
from engcoolprop.ec_fluid import (CPeng_fromSI ,  CondEng_fromSI ,   DSI_fromEng ,  
                                  Deng_fromSI ,   PSI_fromEng, Peng_fromSI,
                                  PropsSI ,   SSI_fromEng ,  Seng_fromSI ,  TSI_fromEng ,  
                                  Teng_fromSI ,  UHSI_fromEng ,  UHeng_fromSI ,   Veng_fromSI )
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from engcoolprop.iteration_utils import calc_Tnbp
from engcoolprop.safe_get_property import safe_get_INCOMP_prop as SC
from engcoolprop.safe_get_property import get_si_prop


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

    def __init__(self,symbol="DowQ", T=500 ,P=1000.0, Pmax=10000.0,
                 show_warnings=True, child=1):
        '''Init generic Fluid'''

        if symbol not in self.fluidNameL:
            raise ValueError( '"%s" is NOT in coolprop incompressible list\n%s'%(symbol, repr(self.fluidNameL) ) )
        
        self.symbol = symbol
        self.Pmax = Pmax # highest pressure considered in any iterative calcs (can still input P > Pmax)
        self.Pmin = 0
        self.show_warnings = show_warnings # some calcs will issue warning statements if show_warnings == True

        self.fluid = 'INCOMP::%s'%symbol
        self.child = child

        # get temperature limits
        self.Tmin_si =   PropsSI('Tmin','T',0,'P',0,'INCOMP::%s'%symbol)  
        self.Tmax_si =   PropsSI('Tmax','T',0,'P',0,'INCOMP::%s'%symbol) 

        self.Tmin =  Teng_fromSI( self.Tmin_si  )
        self.Tmax =  Teng_fromSI( self.Tmax_si )


        # H = 0.0 is the arbitrary ref for the fluid (http://www.coolprop.org/fluid_properties/Incompressibles.html#)
        P_h_ref = 14.696 # 1 atm
        T_h_ref = 527.67 # degR == 20 C

        # desired ref temperature may not be in fluids temperature range
        if T_h_ref <   self.Tmin:
            T_h_ref =  self.Tmin + 0.001 # give a slight move into range
        elif T_h_ref > self.Tmax:
            T_h_ref =  self.Tmax - 0.001 # give a slight move into range

        

        Psi_1_atm = PSI_fromEng( 14.696 ) # SI for 1 atm (Pa)
        try:
            # try calculating T and P for H=0 (i.e. Tref and Pref)
            prop, good_Psi = SC( T, Psi_val=Psi_1_atm, ind_name='H', ind_si_val=0, 
                                             symbol=self.symbol, show_warnings=self.show_warnings )
            self.Tref = Teng_fromSI( prop )
            self.Pref = Peng_fromSI(good_Psi)
        except:
            # set ref P and T to best estimate of desired ref point
            self.Pref = P_h_ref
            self.Tref = T_h_ref


        # try to get Psat at Tref, Pref
        try:
            self.Psat_ref = Peng_fromSI( PropsSI('P','T', TSI_fromEng( self.Tref ) ,'Q',0, self.fluid) )
            if self.Pref < self.Psat_ref:
                # can't use 1 atm as Pref
                self.Pref = 1.001 * self.Psat_ref # give slight nudge into liquid P range
                print( 'Changed Pref from 1 atm to: %g psia'%self.Pref )
        except:
            self.Psat_ref = None # Psat may not be supported by this INCOMP fluid
        # print( 'self.Psat_ref =', self.Psat_ref)

        # print( "EC) self.Tref=%g, self.Pref=%g"%(self.Tref, self.Pref) )
        self.setTP(self.Tref,self.Pref)
        self.Href = self.H # CoolProp uses this Tref, Pref for setting Href and Sref to 0
        self.Sref = self.S
        self.P = P
        self.Pinput = P 

        # If Psat is supported, try to get Tnbp
        if self.Psat_ref is not None:
            try:
                Tnbp, error_code = calc_Tnbp( self )
                if not error_code:
                    self.Tnbp = Tnbp 
                else:
                    self.Tnbp = None 
            except:
                self.Tnbp = None 
        else:
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

        # set properties to input T and P
        self.setTP(self.T, self.P)

        self.calc_min_max_props()
            
        if child==1: 
            self.dup = EC_Incomp_Fluid(symbol=self.symbol, T=self.T, P=self.P, child=0, Pmax=self.Pmax, 
                                       show_warnings=self.show_warnings)

    def setProps(self, **inpD):
        '''Generic call using any P with supported inputs T,D,S,H
        
        #: for example
        #: ec.setProps(T=100, P=200)
        #: ec.setProps(D=0.1, P=100)
        #: ec.setProps(P=100, H=20)

        # Make so both:
        # setProps(T=530, P=100)  AND  
        # setProps(P=100, T=530) will work

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
            self.set_PS( **inpD )
        else:
            raise ValueError( 'Invalid input to setProps ' + repr(inpD) )


    def dH_FromHref(self):
        """Calculate Enthalpy difference from self.Href"""
        return self.H - self.Href

    def get_Psat(self, T):
        """
        Given T in degR (Engineering units) return Psat in psia (if supported)
        NOTE: if not supported, return 0 for Psat
        """
        try:
            Psat = Peng_fromSI( PropsSI('P','T', TSI_fromEng( T ) ,'Q',0, self.fluid) )
        except:
            Psat = 0
        return Psat # psia


    def adjust_P_for_Psat(self, P, ind_name='T', ind_val=None):
        """
        If Psat is greater than P, then to maintain liquid state, use Psat
        
        Given independent variable name and value, calc T and then Psat
        NOTE: P is in Engineering units (psia)
        NOTE: ind_val is in Engineering units (degR, lbm/cuft, etc.)
        """

        if ind_name == 'T':
            Psat = self.get_Psat( ind_val ) # psia
        elif ind_name == 'H':
            T = Teng_fromSI( get_si_prop( PSI_fromEng( P ), targ_prop='T', 
                                          ind_name='H', ind_si_val=UHSI_fromEng( ind_val ), symbol=self.symbol) )
            Psat = self.get_Psat( T )
        elif ind_name == 'D':
            T = Teng_fromSI( get_si_prop( PSI_fromEng( P ), targ_prop='T', 
                                          ind_name='D', ind_si_val=DSI_fromEng( ind_val ), symbol=self.symbol) )
            Psat = self.get_Psat( T )
        elif ind_name == 'S':
            T = Teng_fromSI( get_si_prop( PSI_fromEng( P ), targ_prop='T', 
                                          ind_name='S', ind_si_val=SSI_fromEng( ind_val ), symbol=self.symbol) )
            Psat = self.get_Psat( T )
        else:
            raise ValueError( 'ind_name must be T, H, D or S')
        
        if Psat > P:
            self.Padjusted = Psat + 0.001 # slight increase in Psat to stay in liquid region
        else:
            self.Padjusted = None

        return max( P, Psat )

    def setTP(self,T=530.0,P=1000.0):
        '''Calc props from T and P'''

        self.Pinput = P # save input P in case Psat changes it.
        P = self.adjust_P_for_Psat( P, ind_name='T', ind_val=T)
        
        Tsi = TSI_fromEng( T )
        Psi = PSI_fromEng( P )
        def get_prop( prop_desc='H' ):
            try:
                prop, good_Psi = SC( prop_desc, Psi_val=Psi, ind_name='T', ind_si_val=Tsi, 
                                    symbol=self.symbol, show_warnings=self.show_warnings )
                return prop
            except:
                return float('inf')

        self.T = T
        self.P = P
        self.adjusted = None

        self.D = Deng_fromSI( get_prop('D') )

        self.rho = self.D / 1728.0

        self.E = UHeng_fromSI( get_prop('U') )
        self.H = UHeng_fromSI( get_prop('H') )
        self.S = Seng_fromSI( get_prop('S') )
        self.Cp = CPeng_fromSI( get_prop('C') )
        self.Visc = Veng_fromSI( get_prop('V') ) * 1.0E5
        self.Cond = CondEng_fromSI( get_prop('L') )


        try:
            self.Psat = Peng_fromSI( PropsSI('P','T', TSI_fromEng( self.T ) ,'Q',0, self.fluid) )
        except:
            self.Psat = None        
        
    def constP_newH(self,H):
        '''Calc properties at new H with same P'''

        self.setPH( self.P, H)
        
    def constH_newP(self,P=1000.0):
        '''Calc properties at new P with same H'''

        self.setPH( P, self.H)
        
    def setPH(self,P,H):
        '''Calc properties at P and H'''


        self.Pinput = P # save input P in case Psat changes it.
        P = self.adjust_P_for_Psat( P, ind_name='H', ind_val=H)
        
        self.P = P 
        self.H = H
        
        Psi = PSI_fromEng( P )
        Hsi = UHSI_fromEng( H )

        def get_prop( prop_desc='T' ):
            try:
                prop, good_Psi = SC( prop_desc, Psi_val=Psi, ind_name='H', ind_si_val=Hsi, 
                                    symbol=self.symbol, show_warnings=self.show_warnings )
                return prop
            except:
                return float('inf')

        self.T = Teng_fromSI( get_prop('T') )

        self.D = Deng_fromSI( get_prop('D') )
        self.rho = self.D / 1728.0

        self.E = UHeng_fromSI( get_prop('U') )
        self.S = Seng_fromSI( get_prop('S') )
        self.Cp = CPeng_fromSI( get_prop('C') )
        self.Visc = Veng_fromSI( get_prop('V') ) * 1.0E5
        self.Cond = CondEng_fromSI( get_prop('L') )

        try:
            self.Psat = Peng_fromSI( PropsSI('P','T', TSI_fromEng( self.T ) ,'Q',0, self.fluid) )
        except:
            self.Psat = None
        
    def setPS(self,P,S):
        '''Calc properties at P and H'''

        self.Pinput = P # save input P in case Psat changes it.
        P = self.adjust_P_for_Psat( P, ind_name='S', ind_val=S)

        self.P = P 
        self.S = S
        
        Psi = PSI_fromEng( P )
        Ssi = SSI_fromEng( S )

        def get_prop( prop_desc='H' ):
            try:
                prop, good_Psi = SC( prop_desc, Psi_val=Psi, ind_name='S', ind_si_val=Ssi, 
                                    symbol=self.symbol, show_warnings=self.show_warnings )
                return prop
            except:
                return float('inf')


        self.T = Teng_fromSI( get_prop('T') )

        self.D = Deng_fromSI( get_prop('D') )
        self.rho = self.D / 1728.0

        self.H = UHeng_fromSI( get_prop('H') )
        self.E = UHeng_fromSI( get_prop('U') )
        
        self.Cp = CPeng_fromSI( get_prop('C') )
        self.Visc = Veng_fromSI( get_prop('V') ) * 1.0E5
        self.Cond = CondEng_fromSI( get_prop('L') )

        try:
            self.Psat = Peng_fromSI( PropsSI('P','T', TSI_fromEng( self.T ) ,'Q',0, self.fluid) )
        except:
            self.Psat = None

    def constS_newP(self,P=1000.0):
        '''Calc properties at new P with same S'''

        self.setPS( P, self.S )
        
    def setPD(self,P=1000.0,D=0.01):
        '''
        Calc props from P and D
        NOTE: The pressure has NO EFFECT on incompressible density calc.
        '''

        self.Pinput = P # save input P in case Psat changes it.
        P = self.adjust_P_for_Psat( P, ind_name='D', ind_val=D)

        self.P = P 
        self.D = D
        
        Psi = PSI_fromEng( P )
        Dsi = DSI_fromEng( D )
        self.rho = self.D / 1728.0

        def get_prop( prop_desc='H' ):
            try:
                prop, good_Psi = SC( prop_desc, Psi_val=Psi, ind_name='D', ind_si_val=Dsi, 
                                    symbol=self.symbol, show_warnings=self.show_warnings )
                return prop
            except:
                return float('inf')

        self.T = Teng_fromSI( get_prop('T') )

        self.H = UHeng_fromSI( get_prop('H') )
        self.E = UHeng_fromSI( get_prop('U') )
        self.S = Seng_fromSI( get_prop('S') )

        self.Cp = CPeng_fromSI( get_prop('C') )
        self.Visc = Veng_fromSI( get_prop('V') ) * 1.0E5
        self.Cond = CondEng_fromSI( get_prop('L') )


        try:
            self.Psat = Peng_fromSI( PropsSI('P','T', TSI_fromEng( self.T ) ,'Q',0, self.fluid) )
        except:
            self.Psat = None
        

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
        print("State Point for fluid",self.fluid,"("+ self.symbol +")")
        print("T =%8g"%self.T," degR,",    '                             Range(%8g - %8g) degR'%(self.Tmin, self.Tmax))

        print("P =%8g"%self.P," psia",    '                              Range(%8g - %8g) psia'%(self.Pmin, self.Pmax))
        if self.Padjusted is not None:
            print("    Pinput =%8g"%self.Pinput, '(Adjusted due to Psat)')

        print("D =%8g"%self.D," lbm/cuft",    '                          Range(%8g - %8g) lbm/cuft'%(self.Dmin, self.Dmax))
        print("   rho =%8g"%self.rho," lbm/cuin",   '                    Range(%.6f - %.6f) lbm/cuin'%(self.rho_min, self.rho_max))

        print("E =%8g"%self.E," BTU/lbm",    '                           Range(%8g - %8g) BTU/lbm'%(self.Emin, self.Emax))
        print("H =%8g"%self.H," BTU/lbm",    '                           Range(%8g - %8g) BTU/lbm'%(self.Hmin, self.Hmax))
        print("S =%8g"%self.S," BTU/lbm degR",    '                      Range(%.6f - %.6f) BTU/lbm degR'%(self.Smin, self.Smax))
        print("Cp=%8g"%self.Cp," BTU/lbm degR",   '                      Range(%8g - %8g) BTU/lbm degR'%(self.Cpmin, self.Cpmax))
        print("V =%8g"%self.Visc," viscosity [1.0E5 * lbm/ft-sec]", '    Range(%8g - %8g)'%(self.Viscmin, self.Viscmax) )
        print("C =%8g"%self.Cond," thermal conductivity [BTU/ft-hr-R]", 'Range(%8g - %8g)'%(self.Condmin, self.Condmax) )

        if self.Tnbp is not None:
            print("    Tnbp =%8g"%self.Tnbp," degR,")
        
        if self.Psat is not None:
            print("    Psat =%8g"%self.Psat," psia")

    def dH_FromHZero(self):
        return self.H - self.Href

    def calc_min_max_props(self, do_print=False):

        try:
            Psat_Tmax = Peng_fromSI( PropsSI('P','T', TSI_fromEng( self.Tmax ) ,'Q',0, self.fluid) )
        except:
            Psat_Tmax = 0

        # try:
        #     Psat_Tmin = Peng_fromSI( PropsSI('P','T', TSI_fromEng( self.Tmin ) ,'Q',0, self.fluid) )
        # except:
        #     Psat_Tmin = 0


        def get_prop( T, P, prop_desc='H' ):
            Tsi = TSI_fromEng( T )
            Psi = PSI_fromEng( P )
            try:
                prop, good_Psi = SC( prop_desc, Psi_val=Psi, ind_name='T', ind_si_val=Tsi, 
                                    symbol=self.symbol, show_warnings=self.show_warnings )
                return prop
            except:
                return float('inf')


        self.Dmin = Deng_fromSI( get_prop(self.Tmax, self.Pmin,'D') )
        self.Dmax = Deng_fromSI( get_prop(self.Tmin, self.Pmax,'D') )
        self.rho_min = self.Dmin / 1728.0
        self.rho_max = self.Dmax / 1728.0

        self.Emin = UHeng_fromSI( get_prop(self.Tmin, self.Pmax,'U') )
        self.Emax = UHeng_fromSI(  get_prop(self.Tmax, Psat_Tmax, 'U') )

        self.Hmin = UHeng_fromSI( get_prop(self.Tmin, self.Pmin, 'H') )
        self.Hmax = UHeng_fromSI( get_prop(self.Tmax, self.Pmax, 'H') )
        
        self.Smin = Seng_fromSI( get_prop(self.Tmin, self.Pmax, 'S') )
        self.Smax = Seng_fromSI( get_prop(self.Tmax, Psat_Tmax, 'S') )
        
        self.Cpmin = CPeng_fromSI( get_prop(self.Tmin, self.Pmax, 'C') )
        self.Cpmax = CPeng_fromSI( get_prop(self.Tmax, self.Pmax, 'C') )

        self.Viscmin = Veng_fromSI( get_prop(self.Tmax, self.Pmax, 'V') ) * 1.0E5
        self.Viscmax = Veng_fromSI( get_prop(self.Tmin, self.Pmax, 'V') ) * 1.0E5
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


if __name__ == '__main__':
    
    C = EC_Incomp_Fluid( symbol='Water', T=851, P=0 )
    
    # C.setTP(T= (C.Tmin+C.Tmax)/2.0, P=100) # this temperature throws an exception ???
    
    C.printProps()
    print('='*55)
    
    C = EC_Incomp_Fluid( symbol='Water' )
    
    
    print(C.getStrTransport())
    print(C.getStrTPDphase())
    C.printTPD()
    print('.'*55)
    
    if 1:
        print('='*55)
        t_range = C.Tmax - C.Tmin
        N = 9
        tL = [C.Tmin + (i/N)*t_range for i in range(N)]
        print( ['%.1f'%T for T in tL])
        for T in tL:
            try:
                C.setTP(T, 5000.0)
                C.printTPD()
            except:
                print( 'T=%.1f FAILED'%T )
                Tsi = TSI_fromEng( T )
                Psi = PSI_fromEng( 50000 )

                Dsi = PropsSI('D','T',Tsi,'P',Psi,'INCOMP::Water')
                print( '           Dsi =', Dsi)
            
    print('/'*55)

    C.setTP(500.0, 500.0)
    
    C.constP_newH( C.H * 0.9 )
    C.printTPD()
    C.constH_newP( 500.0 )
    C.printTPD()
    
    C.setPH( 300., C.H * 0.95 )
    C.printTPD()
    
    C.constS_newP(P=1000.0)
    C.printTPD()
    
    C.setPD(P=800.0,D=C.D*0.9)
    C.printTPD()

    print('='*55)
    C.calc_min_max_props( do_print=True )    
    C.printProps()
