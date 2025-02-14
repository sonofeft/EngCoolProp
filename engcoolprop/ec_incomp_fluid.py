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

  Working methods of AbstractState
Prandtl d() = 6.51198382628316
...........................................
Q d() = -inf
...........................................
T d() = 500.0
...........................................
Tmax d() = 633.15
...........................................
Tmin d() = 238.15
...........................................
backend_name d() = IncompressibleBackend
...........................................
conductivity d() = 0.09255528717032524
...........................................
cpmass d() = 2288.1643758645673
...........................................
cvmass d() = 2288.1643758645673
...........................................
has_melting_line d() = False
...........................................
hmass d() = 408385.32146991
...........................................
name d() = DowQ
...........................................
p d() = 101325.0
...........................................
rhomass d() = 809.0654931667821
...........................................
rhomolar d() = -inf
...........................................
smass d() = 1039.0927361604681
...........................................
umass d() = 408260.0843897777
...........................................
viscosity d() = 0.00026340700845078847
...........................................    
"""
import os
from engcoolprop.ec_fluid import (  AbstractState ,  Aeng_fromSI ,  CP ,  CPSI_fromEng ,  
                                  CPeng_fromSI ,  CondEng_fromSI ,  CondSI_fromEng ,  DSI_fromEng ,  
                                  Deng_fromSI ,  EchoInput ,  PSI_fromEng ,  Peng_fromSI , 
                                  PropsSI ,   SSI_fromEng ,  Seng_fromSI ,  TSI_fromEng ,  
                                  Teng_fromSI ,  UHSI_fromEng ,  UHeng_fromSI ,  VSI_fromEng ,  Veng_fromSI ,  
                                  here ,   toSI_callD  )
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP



# ========== build up generic call logic for "setProps" =============
# (i.e. discover which input pairs are supported in T,P,D,Q,S,H,E)

# Make so both:
# setProps(T=530, P=100)  AND  
# setProps(P=100, T=530) will work

# map property names between EngCoolProp and CoolProp
cool_varD = {} # index=eng, value=coolprop  (e.g.  cool_varD['H'] = 'Hmass')
cool_varD['H'] = 'Hmass'
cool_varD['D'] = 'Dmass'
cool_varD['S'] = 'Smass'
cool_varD['E'] = 'Umass'
cool_varD['T'] = 'T'
cool_varD['P'] = 'P'
cool_varD['Q'] = 'Q'

# map Eng to SI conversion functions for each of the fluid properties
toSI_callD = {} # index=eng, value=conversion func (e.g. SSI_fromEng)
toSI_callD['H'] = UHSI_fromEng
toSI_callD['D'] = DSI_fromEng
toSI_callD['S'] = SSI_fromEng
toSI_callD['E'] = UHSI_fromEng
toSI_callD['T'] = TSI_fromEng
toSI_callD['P'] = PSI_fromEng
toSI_callD['Q'] = EchoInput

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

call_InputD = {}  # key:tuple indep vars (e.g. ("T","P")), value: AbstractState input (PT_INPUTS, DmassP_INPUTS, etc.)
call_InputD[('D', 'P')] = CP.DmassP_INPUTS
call_InputD[('H', 'P')] = CP.HmassP_INPUTS
call_InputD[('P', 'D')] = CP.DmassP_INPUTS
call_InputD[('P', 'H')] = CP.HmassP_INPUTS
call_InputD[('P', 'S')] = CP.PSmass_INPUTS
call_InputD[('P', 'T')] = CP.PT_INPUTS
call_InputD[('S', 'P')] = CP.PSmass_INPUTS
call_InputD[('T', 'P')] = CP.PT_INPUTS


# make a list of all incompressible fluids in coolprop
incomp_pure_fluidL = CP.get_global_param_string('incompressible_list_pure').split(',')
# print( 'incomp_pure_fluidL =', incomp_pure_fluidL)

class EC_Incomp_Fluid(object):
    
    fluidNameL = incomp_pure_fluidL

    def __init__(self,symbol="DowQ",T=530.0,P=1000.0, child=1):
        '''Init generic Fluid'''

        if symbol not in self.fluidNameL:
            raise ValueError( '"%s" is NOT in coolprop incompressible list\n%s'%(symbol, repr(self.fluidNameL) ) )
        
        self.symbol = symbol #.upper()
        
        # AS = CP.AbstractState("INCOMP", symbol)
        # self.AS = AS

        self.fluid = 'INCOMP::%s'%symbol
        self.T = T
        self.P = P
        self.child = child

        # get temperature limits
        self.Tmin =  Teng_fromSI( PropsSI('Tmin','T',0,'P',0,'INCOMP::%s'%symbol)  )
        self.Tmax =  Teng_fromSI( PropsSI('Tmax','T',0,'P',0,'INCOMP::%s'%symbol) )
        
        # set Href by using SATP
        self.Tref = 536.67 # 536.67R = (SATP, Standard Ambient T P) = 77F, 25C 
        self.Pref = 14.7
        self.setTP(self.Tref,self.Pref)
        self.Href = self.H

        # set properties to input T and P
        self.setTP(T,P)
            
        if child==1: 
            self.dup = EC_Incomp_Fluid(symbol=self.symbol, T=self.T, P=self.P, child=0)

    def setProps(self, **inpD):
        '''Generic call using any two of supported inputs T,P,D,Q,S,H,E::
        
        #: for example
        #: ec.setProps(T=100, D=0.1)
        #: ec.setProps(D=0.1, T=100)
        #: ec.setProps(P=100, H=20)
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

    def setTP(self,T=530.0,P=1000.0):
        '''Calc props from T and P'''
        
        # CPeng_fromSI ,  CondEng_fromSI ,  CondSI_fromEng ,  DSI_fromEng ,  
        # Deng_fromSI ,  Peng_fromSI ,  Seng_fromSI , 
        # Teng_fromSI ,  UHeng_fromSI ,   Veng_fromSI

        Tsi = TSI_fromEng( T )
        Psi = PSI_fromEng( P )

        self.T = T
        self.P = P

        self.D = Deng_fromSI( PropsSI('D','T',Tsi,'P',Psi, self.fluid ) )
        self.rho = self.D / 1728.0

        self.E = UHeng_fromSI( PropsSI('U','T',Tsi,'P',Psi, self.fluid ) )
        self.H = UHeng_fromSI( PropsSI('H','T',Tsi,'P',Psi, self.fluid ) )
        self.S = Seng_fromSI( PropsSI('S','T',Tsi,'P',Psi, self.fluid ) )
        self.Cp = CPeng_fromSI( PropsSI('C','T',Tsi,'P',Psi, self.fluid ) )
        self.Visc = Veng_fromSI( PropsSI('V','T',Tsi,'P',Psi, self.fluid ) )
        self.Cond = CondSI_fromEng( PropsSI('L','T',Tsi,'P',Psi, self.fluid ) )


    def newDE(self,D=0.1,E=50.0):
        '''Calc properties at new D and E'''
        Dsi = DSI_fromEng( D )
        Esi = UHSI_fromEng( E )

        #print('Failed in newDE for D=%g, E=%g'%(D, E))
        '''iterate on temperature until internal energy is correct
        tRbegin--beginning search temperature [R]'''
        
        Tbegin = self.T
        
        tolr = 1.0E-8
        tR = Tbegin
        #print '---'
        for i in range(48): # limit number of iterations
            self.setTD(T=tR, D=D)
            dedt = self.Cv
            eError = E-self.E
            
            #print 'dedt=',dedt,'  eError=',eError,'  tR=',tR
            
            if self.Cv==0.0:
                print( '==> ERROR in wrap_dll.  CV=0 for tR=',tR )
            
            tR = tR + eError / dedt
            
            if abs(eError) <= tolr:
                break
            
        
        
        
    def constP_newH(self,H):
        '''Calc properties at new H with same P'''

        self.setPH( self.P, H)
        
    def constH_newP(self,P=1000.0):
        '''Calc properties at new P with same H'''

        self.setPH( P, self.H)
        
    def setPH(self,P,H):
        '''Calc properties at P and H'''

        self.P = P 
        self.H = H
        
        Psi = PSI_fromEng( P )
        Hsi = UHSI_fromEng( H )

        self.T = Teng_fromSI( PropsSI('T','H',Hsi,'P',Psi, self.fluid ) )

        self.D = Deng_fromSI( PropsSI('D','H',Hsi,'P',Psi, self.fluid ) )
        self.rho = self.D / 1728.0

        self.E = UHeng_fromSI( PropsSI('U','H',Hsi,'P',Psi, self.fluid ) )
        self.S = Seng_fromSI( PropsSI('S','H',Hsi,'P',Psi, self.fluid ) )
        self.Cp = CPeng_fromSI( PropsSI('C','H',Hsi,'P',Psi, self.fluid ) )
        self.Visc = Veng_fromSI( PropsSI('V','H',Hsi,'P',Psi, self.fluid ) )
        self.Cond = CondSI_fromEng( PropsSI('L','H',Hsi,'P',Psi, self.fluid ) )
        
    def setPS(self,P,S):
        '''Calc properties at P and H'''

        self.P = P 
        self.S = S
        
        Psi = PSI_fromEng( P )
        Ssi = UHSI_fromEng( S )

        self.T = Teng_fromSI( PropsSI('T','S',Ssi,'P',Psi, self.fluid ) )

        self.D = Deng_fromSI( PropsSI('D','S',Ssi,'P',Psi, self.fluid ) )
        self.rho = self.D / 1728.0

        self.H = UHeng_fromSI( PropsSI('H','S',Ssi,'P',Psi, self.fluid ) )
        self.E = UHeng_fromSI( PropsSI('U','S',Ssi,'P',Psi, self.fluid ) )
        self.Cp = CPeng_fromSI( PropsSI('C','S',Ssi,'P',Psi, self.fluid ) )
        self.Visc = Veng_fromSI( PropsSI('V','S',Ssi,'P',Psi, self.fluid ) )
        self.Cond = CondSI_fromEng( PropsSI('L','S',Ssi,'P',Psi, self.fluid ) )

    def constS_newP(self,P=1000.0):
        '''Calc properties at new P with same S'''

        self.setPS( P, self.S )

    def setTD(self,T=530.0,D=0.01):
        '''Calc P from T and D'''
        
        Tsi = TSI_fromEng( T )
        Dsi = DSI_fromEng( D )

        raise Exception( 'setTD not ready' )
        
    def setPD(self,P=1000.0,D=0.01):
        '''Calc props from P and D'''

        self.P = P 
        self.D = D
        
        Psi = PSI_fromEng( P )
        Dsi = DSI_fromEng( D )
        self.rho = self.D / 1728.0

        self.T = Teng_fromSI( PropsSI('T','D',Dsi,'P',Psi, self.fluid ) )
        self.H = UHeng_fromSI( PropsSI('H','D',Dsi,'P',Psi, self.fluid ) )
        self.E = UHeng_fromSI( PropsSI('U','D',Dsi,'P',Psi, self.fluid ) )
        self.S = Seng_fromSI( PropsSI('S','D',Dsi,'P',Psi, self.fluid ) )
        self.Cp = CPeng_fromSI( PropsSI('C','D',Dsi,'P',Psi, self.fluid ) )
        self.Visc = Veng_fromSI( PropsSI('V','D',Dsi,'P',Psi, self.fluid ) )
        self.Cond = CondSI_fromEng( PropsSI('L','D',Dsi,'P',Psi, self.fluid ) )

        

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

    def initFromObj(self, obj):
        '''initialize properties from another n_fluid object'''
        if  self.symbol.upper() == obj.symbol.upper():
            #self.setTD(T=obj.T,D=obj.D)
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

    def html_desc(self):
        '''html output a multiline property summary with units'''
        return '<table><th colspan="3" align="left">&nbsp;&nbsp;&nbsp;for fluid "'+ str(self.fluid)+ " (" + str(self.symbol) + ')"</th>' \
            + "<tr><td>E = </td><td>" + "%5g"%self.E + "</td><td> BTU/lbm"+ "</td></tr>" \
            + "<tr><td>H = </td><td>" + "%5g"%self.H + "</td><td> BTU/lbm"+ "</td></tr>" \
            + "<tr><td>S = </td><td>" + "%5g"%self.S + "</td><td> BTU/lbm degR"+ "</td></tr>" \
            + "<tr><td>Cp= </td><td>" + "%5g"%self.Cp + "</td><td>BTU/lbm degR"+ "</td></tr></table>"


    def printProps(self):
        '''print a multiline property summary with units'''
        print("State Point for fluid",self.fluid,"("+self.symbol+")", sep=' ')
        print("T =%8g"%self.T," degR,", sep=' ')
            
        print("P =%8g"%self.P," psia", sep=' ')
        print("D =%8g"%self.D," lbm/cu ft", sep=' ')
        print("E =%8g"%self.E," BTU/lbm", sep=' ')
        print("H =%8g"%self.H," BTU/lbm", sep=' ')
        print("S =%8g"%self.S," BTU/lbm degR", sep=' ')
        print("Cp=%8g"%self.Cp," BTU/lbm degR", sep=' ')
        print("V =%8g"%self.Visc," viscosity [1.0E5 * lbm/ft-sec]", sep=' ')
        print("C =%8g"%self.Cond," thermal conductivity [BTU/ft-hr-R]", sep=' ')

    def dH_FromHZero(self):
        return self.H - self.Href

if __name__ == '__main__':
    
    C = EC_Incomp_Fluid( symbol='DowQ' )
    
    C.setProps(T=C.Tref, P=C.Pref)
    #C.setProps(P=5., Q=0.5)
    #C.setProps(H=120., E=0.5) # illegal inputs
    
    C.printProps()
    print('='*55)
    
    C = EC_Incomp_Fluid( symbol='Water' )
    
    
    print(C.getStrTransport())
    print(C.getStrTPDphase())
    C.printTPD()
    print('.'*55)
    print( C.html_desc() )
    
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
            
    print('$'*55)
    # C.newDE(D=1,E=50.0)
    # C.printTPD()

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
    
    C.setTD(T=500.0,D=5.0)
    C.printTPD()
    
