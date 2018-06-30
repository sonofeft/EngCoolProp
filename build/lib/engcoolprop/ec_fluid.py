#!/usr/bin/env python
# -*- coding: ascii -*-
from __future__ import print_function

r"""
??? Tfreeze ???
EngCoolProp is a python Engineering Units wrapper around the CoolProp project.

<Paragraph description see docstrings at http://www.python.org/dev/peps/pep-0257/>
EngCoolProp uses units of primarily inch, lbm, lbf, sec, BTU (some use of ft and hour).

    The following are the default units for each property.

    T = Temperature = degR
    Tc= Critical Temperature = degR 
    Tnbp = Normal Boiling Point = degR
    P = Pressure = psia
    Pc = Critical Pressure = psia
    D = Density = lbm/cu ft
    rho = Density = lbm/cu inch
    Dc = Critical Density = lbm/cu ft
    E = Internal Energy = BTU/lbm
    H = Enthalpy = BTU/lbm
    S = Entropy = BTU/lbm degR
    Cv = Heat Capacity (const. V) = BTU/lbm degR
    Cp = Heat Capacity (const. P) = BTU/lbm degR
    g = Ratio of Specific Heats (Cp/Cv) = (-)
    A = Sonic Velocity = ft/sec
    V = Viscosity = 1.0E5 * lb/ft-sec
    C = Thermal Conductivity = BTU/ft-hr-R
    MW = Molecular Weight = lbm/lbmmole
    Q = Quality (mass fraction gas) = (-)
    Z = Compressibility() = (-)



EngCoolProp
Copyright (C) 2018  Applied Python

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

-----------------------

"""
import os
here = os.path.abspath(os.path.dirname(__file__))


# for multi-file projects see LICENSE file for authorship info
# for single file projects, insert following information
__author__ = 'Charlie Taylor'
__copyright__ = 'Copyright (c) 2018 Charlie Taylor'
__license__ = 'GPL-3'
exec( open(os.path.join( here,'_version.py' )).read() )  # creates local __version__ variable
__email__ = "cet@appliedpython.com"
__status__ = "4 - Beta" # "3 - Alpha", "4 - Beta", "5 - Production/Stable"

#
import string, os, sys
from math import exp
from CoolProp.CoolProp import PropsSI
from CoolProp import AbstractState
import CoolProp.CoolProp as CP

# Q values to indicat all liquid or gas
Q_LIQUID = 0  # all liquid
Q_GAS = 1     # all gas

def Deng_fromSI( D ):
    return D  / 16.01843417

def DSI_fromEng( D ):
    return D * 16.01843417

def Teng_fromSI( T ):
    return T * 1.8

def TSI_fromEng( T ):
    return T / 1.8

def Peng_fromSI( P ):
    return P / 6894.76

def PSI_fromEng( P ):
    return P * 6894.76

def UHeng_fromSI( UH ):
    return UH * 0.000429923

def UHSI_fromEng( UH ):
    return UH / 0.000429923
    
def Seng_fromSI( S ):
    return S * 0.0002388461111111111

def SSI_fromEng( S ):
    return S / 0.0002388461111111111

def CPeng_fromSI( S ):
    return S * 0.0002388461111111111

def CPSI_fromEng( S ):
    return S / 0.0002388461111111111

def Aeng_fromSI( A ):
    return A * 3.28084
    
def ASI_fromEng( A ):
    return A / 3.28084

def Veng_fromSI( V ):
    return V * 0.6719689751
    
def VSI_fromEng( V ):
    return V / 0.6719689751

def CondEng_fromSI( Cond ):
    return Cond * 0.578176

def CondSI_fromEng( Cond ):
    return Cond / 0.578176

def EchoInput( x ):
    return x

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
call_tuplesD = {} # index=(ec1, ec2), value=(XX_INPUTS, c1, c2)

for s in dir(CP):
    if s.endswith('_INPUTS'):
        if s.find('molar')<0:
            name = s[:-7]
            for ec1, c1 in list(cool_varD.items()):
                if name.startswith(c1):
                    for ec2, c2 in list(cool_varD.items()):
                        if name.endswith(c2):
                            # put in both orderings of input pairs
                            call_tuplesD[(ec1, ec2)] = (getattr(CP, s), ec1, ec2)
                            call_tuplesD[(ec2, ec1)] = (getattr(CP, s), ec1, ec2)
            

class EC_Fluid(object):

    def __init__(self,symbol="N2",T=530.0,P=1000.0, child=1):
        '''Init generic Fluid'''
        
        self.symbol = symbol.upper()
        
        # http://www.coolprop.org/_static/doxygen/html/class_cool_prop_1_1_abstract_state.html
        AS = AbstractState("HEOS", self.symbol)
        self.AS = AS
        
        self.name = AS.name()
        self.T = T
        self.P = P
        self.child = child
        
        Tsi = TSI_fromEng( T )
        Psi = PSI_fromEng( P )
        AS.update(CP.PT_INPUTS, Psi, Tsi)

        self.WtMol = AS.molar_mass() * 1000.0
        self.Tc = Teng_fromSI( AS.T_critical() )
        self.Pc = Peng_fromSI( AS.p_critical() )
        self.Dc = Deng_fromSI( AS.rhomass_critical() )
        self.Ttriple = Teng_fromSI( AS.Ttriple() )


        self.Tmin =  Teng_fromSI( AS.Tmin() )
        self.Tmax =  Teng_fromSI( AS.Tmax() )
        self.Pmax = Peng_fromSI( AS.pmax() )
        
        dcIdeal = self.Pc*self.WtMol/self.Tc/10.729
        self.Zc = dcIdeal / self.Dc
        
        TnbpSI = PropsSI("T", "P", 101325, "Q", Q_GAS, self.symbol)
        self.Tnbp = Teng_fromSI( TnbpSI )
        
        if self.Tnbp < 536.67:
            self.Tref = self.Tnbp
        else:
            self.Tref = 536.67
        self.Pref = 14.7
        
        #print( 'About to call setTP #1 with Tref, Pref=',self.Tref,self.Pref )
        self.setTP(self.Tref,self.Pref)
        #print( 'Back from call setTP')
        self.Href = self.H
        #print( 'About to call setTP #2')
        self.setTP(T,P)
        #print( 'Back from call setTP')
            
        if child==1: 
            self.dup = EC_Fluid(symbol=self.symbol,T=self.T,P=self.P, child=0)
            
        self.calcdFreezePt = 0

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
            print('Allowable inputs are: T,P,D,Q,S,H,E')
            if len(inpD) != 2:
                print('ONLY 2 inputs are allowed.  (e.g. T=530.0,P=100.0)')
            raise ValueError('setProps called with illegal inputs: %s'%str(inpD))
        
        
        INPUTS, ec1, ec2 = call_tuplesD[key]
        P1si = toSI_callD[ec1]( inpD[ec1] )
        P2si = toSI_callD[ec2]( inpD[ec2] )
        #print('Calling',key,'with',call_tuplesD[key],'P1si=%g,  P2si=%g'%(P1si, P2si))
        
        self.AS.update(INPUTS, P1si, P2si)
        
        self.setPropsFromAS()

    def dH_FromHref(self):
        """Calculate Enthalpy difference from self.Href"""
        return self.H - self.Href

    def setPropsFromAS(self):
        """Assume that self.AS (AbstractState) has been set elsewhere."""
        
        AS = self.AS

        self.T = Teng_fromSI( AS.T() )
        self.P = Peng_fromSI( AS.p() )

        self.H = UHeng_fromSI( AS.hmass() )
        self.D = Deng_fromSI( AS.rhomass() )
        self.rho = self.D / 1728.0
        #print( "rho=",rho,"in setTP")
        self.E = UHeng_fromSI( AS.umass() )
        self.S = Seng_fromSI( AS.smass() )
        self.Cv = CPeng_fromSI( AS.cvmass() )
        self.Cp = CPeng_fromSI( AS.cpmass() )
        self.sonicV = Aeng_fromSI( AS.speed_sound() )
        self.Visc = Veng_fromSI( AS.viscosity() ) * 1.0E5
        self.Cond = CondEng_fromSI( AS.conductivity() )
        self.Q = AS.Q()  # maybe use AS.PIP() ???
        
        dIdeal = self.P*self.WtMol/self.T/10.729
        self.Z = dIdeal / self.D

    def setTP(self,T=530.0,P=1000.0):
        '''Calc props from T and P'''
        
        Tsi = TSI_fromEng( T )
        Psi = PSI_fromEng( P )
        self.AS.update(CP.PT_INPUTS, Psi, Tsi)
        
        self.setPropsFromAS()
    
    def newDE(self,D=0.1,E=50.0):
        '''Calc properties at new D and E'''
        Dsi = DSI_fromEng( D )
        Esi = UHSI_fromEng( E )

        self.AS.update(CP.DmassUmass_INPUTS, Dsi, Esi)
        
        self.setPropsFromAS()
        
    def constP_newH(self,H):
        '''Calc properties at new H with same P'''
        
        Psi = PSI_fromEng( self.P )
        Hsi = UHSI_fromEng( H )
        self.AS.update(CP.HmassP_INPUTS, Hsi, Psi)
        
        self.setPropsFromAS()
        
    def constH_newP(self,P=1000.0):
        '''Calc properties at new P with same H'''
        
        Psi = PSI_fromEng( P )
        Hsi = UHSI_fromEng( self.H )
        
        self.AS.update(CP.HmassP_INPUTS, Hsi, Psi)
        self.setPropsFromAS()
        
    def setPH(self,P,H):
        '''Calc properties at P and H'''
        
        Psi = PSI_fromEng( P )
        Hsi = UHSI_fromEng( H )
        
        self.AS.update(CP.HmassP_INPUTS, Hsi, Psi)
        self.setPropsFromAS()

    def constS_newP(self,P=1000.0):
        '''Calc properties at new P with same S'''
        
        Psi = PSI_fromEng( P )
        Ssi = SSI_fromEng( self.S )
        self.AS.update(CP.PSmass_INPUTS, Psi, Ssi)
        self.setPropsFromAS()

    def setTD(self,T=530.0,D=0.01):
        '''Calc P from T and D'''
        
        Tsi = TSI_fromEng( T )
        Dsi = DSI_fromEng( D )
        self.AS.update(CP.DmassT_INPUTS, Dsi, Tsi)
        self.setPropsFromAS()
        
    def setPD(self,P=1000.0,D=0.01):
        '''Calc props from P and D'''
        
        Psi = PSI_fromEng( P )
        Dsi = DSI_fromEng( D )
        self.AS.update(CP.DmassP_INPUTS, Dsi, Psi)
        self.setPropsFromAS()

    def getSatP(self, T=None):
        """Get Psat for Liquid at input T or self.T"""
        if not(T is None) and (T < self.Tc):
            Tsi = TSI_fromEng( T )
            
            self.dup.AS.update(CP.QT_INPUTS, Q_LIQUID, Tsi)
            Psat = Peng_fromSI( self.dup.AS.p() )
        else:
            print('ERROR... getSatP called with T>Tc.   T=%s,   Tc=%g'%(T, self.Tc))
            Psat = self.P
            
        return Psat

    def getSatPandDens(self, T=None):
        """Get Psat and rho for both Liquid & Gas at input T or self.T"""
        if not(T is None) and (T < self.Tc):
            Tsi = TSI_fromEng( T )
            
            
            self.dup.AS.update(CP.QT_INPUTS, Q_GAS, Tsi)
            Dgas = self.D
            
            self.dup.AS.update(CP.QT_INPUTS, Q_LIQUID, Tsi)
            Psat = Peng_fromSI( self.dup.AS.p() )
            Dliq = Deng_fromSI( self.dup.AS.rhomass() )
            
            self.getSatP(T)
        else:
            print('ERROR... getSatP called with T>Tc.   T=%s,   Tc=%g'%(T, self.Tc))
            Psat = self.P
            Dgas = Dliq = self.D
            
            
        return Psat, Dliq, Dgas
        
    def getSatT(self, P=None):
        """Get Tsat for Liquid at input P or self.P"""
        if not(P is None) and (P < self.Pc):
            Psi = PSI_fromEng( P )
            self.dup.AS.update(CP.PQ_INPUTS, Psi, 1)
            Tsat = Teng_fromSI( self.dup.AS.T() )
        else:
            print('ERROR... getSatT called with P>Pc.   P=%s,   Pc=%g'%(P, self.Pc))
            Tsat = self.T
        
        return Tsat
            
    #for s in dir(CP):
    #    if s.endswith('_INPUTS'):
    #        print s,

    def getSatTandDens(self,P=None):
        """Get Tsat and Density for both Liquid & Gas at input P or self.P"""
        if not(P is None) and (P < self.Pc):
            Psi = PSI_fromEng( P )
            
            self.dup.AS.update(CP.PQ_INPUTS, Psi, Q_GAS)
            Dgas = Deng_fromSI( self.dup.AS.rhomass() )

            self.dup.AS.update(CP.PQ_INPUTS, Psi, Q_LIQUID)
            Tsat = Teng_fromSI( self.dup.AS.T() )
            Dliq = Deng_fromSI( self.dup.AS.rhomass() )
        else:
            print('ERROR... getSatT called with P>Pc.   P=%s,   Pc=%g'%(P, self.Pc))
            Tsat = self.T
            Dgas = Dliq = self.D
        
        return Tsat, Dliq, Dgas

    def dHvap(self, P=None):
        '''calculate heat of vaporization at input P or self.P (BTU/lbm)'''
        
        if P is None:
            P = self.P
        
        if (P < self.Pc):
            
            Psi = PSI_fromEng( P )
            self.dup.AS.update(CP.PQ_INPUTS, Psi, Q_LIQUID)
            Hliq = UHeng_fromSI( self.dup.AS.hmass() )
            
            self.dup.AS.update(CP.PQ_INPUTS, Psi, Q_GAS)
            Hgas = UHeng_fromSI( self.dup.AS.hmass() )

            return Hgas - Hliq
        else:
            print('ERROR in dHvap P=%g, Pc=%g'%(P, self.Pc))
            return 0.0
    
    def surfTen(self):
        '''calculate surface tension, lbf / inch'''
        
        st = self.AS.surface_tension() # units = N/m
        try:
            st_eng = st * 0.00571014709755764 # convert to lbf / inch
        except:
            st_eng = st
        #print('st = %s N/m = %s'%(st, st_eng))
        return st_eng # return lbf / inch
        

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
        self.Cv = self.dup.Cv
        self.sonicV = self.dup.sonicV
        self.Visc = self.dup.Visc
        self.Cond = self.dup.Cond
        self.Q = self.dup.Q

    def saveToDup(self):
        '''save properties to duplicate n_fluid'''
        self.dup.T = self.T
        self.dup.P = self.P
        self.dup.D = self.D
        self.dup.E = self.E
        self.dup.H = self.H
        self.dup.S = self.S
        self.dup.Cp = self.Cp
        self.dup.Cv = self.Cv
        self.dup.sonicV = self.sonicV
        self.dup.Visc = self.Visc
        self.dup.Cond = self.Cond
        self.dup.Q = self.Q

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
            self.Cv = obj.Cv
            self.sonicV = obj.sonicV
            self.Visc = obj.Visc
            self.Cond = obj.Cond
            self.Q = obj.Q
        else:
            raise Exception('Wrong fluid for initializing')


    def gamma(self):
        '''calculate ratio of specific heats from Cp and Cv'''
        try:
            if self.Cp <= 0.0 or self.Cv <=0.0:
                g = 1.0 / 0.0
            g = self.Cp / self.Cv
            return g
        except:
            return 1.000001

    def getStrTransport(self):
        '''create a string from the Transport properties'''
        return  "%s Cp=%6g Cv=%6g gamma=%.4f Visc=%6g ThCond=%6g" %\
        (self.symbol,self.Cp, self.Cv, self.gamma(), self.Visc, self.Cond)

    def getStrTPD(self):
        '''create a string from the TPDEHS properties'''
        return  "%s T=%6.1f P=%6.1f D=%.4f E=%6.2f H=%6.2f S=%.3f Q=%.2f" %\
        (self.symbol,self.T, self.P, self.D, self.E, self.H, self.S, self.Q)

    def getStrTPDphase(self):
        '''create a string from the TPDEHS properties'''
        return  "%s T=%6.1f P=%6.1f D=%.4f E=%6.2f H=%6.2f S=%.3f %s" %\
        (self.symbol,self.T, self.P, self.D, self.E, self.H, self.S, self.Qdescription())

    def printTPD(self):
        '''print a string from the TPDEHS properties'''
        print(self.getStrTPD())
        
    def Qdescription(self):
        '''CoolProp Q values areDifferent from RefProp::
        
        #:c == RefProp ==
        #:c        q--vapor quality on a MOLAR basis [moles vapor/total moles]
        #:c           q < 0 indicates subcooled (compressed) liquid
        #:c           q = 0 indicates saturated liquid
        #:c           q = 1 indicates saturated vapor
        #:c           q > 1 indicates superheated vapor
        #:c           q = -998 subcooled liquid, but quality not defined (p > Pc)
        #:c           q = 999 indicates supercritical state (t > Tc) and (p > Pc)
        #:
        #:c == CoolProp ==
        #:c        q--vapor quality on a MASS basis [mass vapor/total mass]
        #:c           q  indicates subcooled (compressed) liquid
        #:c           q  indicates saturated liquid
        #:c           q  indicates saturated vapor
        #:c           q < 0 indicates superheated vapor OR (p > Pc)
        #:c           q  subcooled liquid, but quality not defined (p > Pc)
        #:c           q  indicates supercritical state (t > Tc) and (p > Pc)
        #:
        #:Phase String            Phase Region
        #:liquid                  p < pcrit & T < Tcrit ; above saturation
        #:gas                     p < pcrit & T < Tcrit ; below saturation
        #:twophase                p < pcrit & T < Tcrit ; mixed liquid/gas
        #:
        #:supercritical_liquid    p > pcrit & T < Tcrit
        #:supercritical_gas       p < pcrit & T > Tcrit
        #:supercritical           p > pcrit & T > Tcrit        
        
        '''
        if self.P >= self.Pc: # supercritical pressure
            if self.T < self.Tc:
                return 'Supercritical Liquid'
            else:
                return 'Supercritical Fluid'
        else:
            # P < Pc
            if self.T >= self.Tc:
                return 'Supercritical Gas'
            else:
                # both T<Tc  AND  P<Pc
                if self.Q <= 0.0:
                    return "Liquid"
                elif self.Q > 900.0:
                    return "Supercrit"
                elif self.Q >= 1.0:
                    return "Gas"
                else:
                    pcent = self.Q * 100.0
                    return "%.2f "%pcent + " %gas"

    def html_desc(self):
        '''html output a multiline property summary with units'''
        return '<table><th colspan="3" align="left">&nbsp;&nbsp;&nbsp;for fluid "'+ str(self.name)+ " (" + str(self.symbol) + ')"</th>' \
            + "<tr><td>T = </td><td>" + "%5g"%self.T + "</td><td> degR &nbsp;&nbsp;(Tc=" + "%5g"%self.Tc + ")"+ "</td></tr>" \
            + "<tr><td>P = </td><td>" + "%5g"%self.P + "</td><td> psia &nbsp;&nbsp;&nbsp;(Pc=" + "%5g"%self.Pc + ")"+ "</td></tr>" \
            + "<tr><td>D = </td><td>" + "%5g"%self.D + "</td><td> lbm/cu ft &nbsp;&nbsp;(Dc=" + "%5g"%self.Dc + ")"+ "</td></tr>" \
            + "<tr><td>E = </td><td>" + "%5g"%self.E + "</td><td> BTU/lbm"+ "</td></tr>" \
            + "<tr><td>H = </td><td>" + "%5g"%self.H + "</td><td> BTU/lbm"+ "</td></tr>" \
            + "<tr><td>S = </td><td>" + "%5g"%self.S + "</td><td> BTU/lbm degR"+ "</td></tr>" \
            + "<tr><td>Cv= </td><td>" + "%5g"%self.Cv + "</td><td> BTU/lbm degR"+ "</td></tr>" \
            + "<tr><td>Cp= </td><td>" + "%5g"%self.Cp + "</td><td>BTU/lbm degR"+ "</td></tr>" \
            + "<tr><td>A = </td><td>" + "%5g"%self.sonicV + "</td><td> ft/sec"+ "</td></tr>" \
            + "<tr><td>MW= </td><td>" + "%5g"%self.WtMol + "</td><td> lbm/lbmmole"+ "</td></tr>" \
            + "<tr><td>Q = </td><td>" + "%s"%self.Qdescription() + "</td><td> phase"+ "</td></tr>" \
            + "<tr><td>Z = </td><td>" + "%5g"%self.Z + "</td><td>(-)"+ "</td></tr></table>"

    def printCriticalProps(self):
        '''print a multiline property summary with units'''
        print("Critical Properties for fluid",self.name,"("+self.symbol+")")
        print("Tc=%8g"%self.Tc,"degR")
        print("Tnbp=%8g"%self.Tnbp,"degR")
        print("Pc =%8g"%self.Pc,"psia")
        print("Dc =%8g"%self.Dc,"lbm/cu ft")
        print("MW=%8g"%self.WtMol," lbm/lbmmole")
        print("Zc =%8g"%self.Zc," (-)")
        
        try:
            
            Tsi = TSI_fromEng( self.Tc )
            Dsi = DSI_fromEng( self.Dc )
            self.dup.AS.update(CP.DmassT_INPUTS, Dsi, Tsi)
            
            H = UHeng_fromSI( self.dup.AS.hmass() )
            print("Hc =%8g"%H," BTU/lbm", sep=' ')
            E = UHeng_fromSI( self.dup.AS.umass() )
            print("Ec =%8g"%E," BTU/lbm", sep=' ')
            S = Seng_fromSI( self.dup.AS.smass() )
            print("Sc =%8g"%S," BTU/lbm degR", sep=' ')
            #Cv = CPeng_fromSI( self.dup.AS.cvmass() )
            #print("Cvc=%8g"%Cv," BTU/lbm degR", sep=' ')
            #Cp = CPeng_fromSI( self.dup.AS.cpmass() )
            #print("Cpc=%8g"%Cp," BTU/lbm degR", sep=' ')
            sonicV = Aeng_fromSI( self.dup.AS.speed_sound() )
            print("Ac =%8g"%sonicV," ft/sec", sep=' ')
            Visc = Veng_fromSI( self.dup.AS.viscosity() ) * 1.0E5
            print("Vc =%8g"%Visc," viscosity [1.0E5 * lb/ft-sec]", sep=' ')
            Cond = CondEng_fromSI( self.dup.AS.conductivity() )
            print("Cc =%8g"%Cond," thermal conductivity [BTU/ft-hr-R]", sep=' ')
        except:
            pass

    def printProps(self):
        '''print a multiline property summary with units'''
        print("State Point for fluid",self.name,"("+self.symbol+")", sep=' ')
        print("T =%8g"%self.T," degR (Tc=%8g"%self.Tc,", Tnbp=%8g"%self.Tnbp, "Ttriple=%8g"%self.Ttriple,")", sep=' ')
        print("P =%8g"%self.P," psia (Pc=%8g"%self.Pc,")", sep=' ')
        print("D =%8g"%self.D," lbm/cu ft (Dc=%8g"%self.Dc,")", sep=' ')
        print("E =%8g"%self.E," BTU/lbm", sep=' ')
        print("H =%8g"%self.H," BTU/lbm", sep=' ')
        print("S =%8g"%self.S," BTU/lbm degR", sep=' ')
        print("Cv=%8g"%self.Cv," BTU/lbm degR", sep=' ')
        print("Cp=%8g"%self.Cp," BTU/lbm degR", sep=' ')
        print("g =%8g"%self.gamma()," Cp/Cv (-)", sep=' ')
        print("A =%8g"%self.sonicV," ft/sec", sep=' ')
        print("V =%8g"%self.Visc," viscosity [1.0E5 * lb/ft-sec]", sep=' ')
        print("C =%8g"%self.Cond," thermal conductivity [BTU/ft-hr-R]", sep=' ')
        print("MW=%8g"%self.WtMol," lbm/lbmmole", sep=' ')
        print("Q =%8g"%self.Q," Vapor Quality (mass fraction gas)", sep=' ')
        print("Z =%8g"%self.Z," (-)", sep=' ')


if __name__ == '__main__':
    
    C = EC_Fluid()
    
    C.setProps(T=120., Q=0.5)
    #C.setProps(P=5., Q=0.5)
    #C.setProps(H=120., E=0.5) # illegal inputs
    
    C.printProps()
    print('-'*55)
    
    
    
    print(C.getStrTransport())
    print(C.getStrTPDphase())
    C.printTPD()
    print('.'*55)
    C.printCriticalProps()
    print()
    #print( C.html_desc() )
    
    if 1:
        print('='*55)
        for T in [115., 150., 200., 250., 300., 350., 400., 500., 550.]:
            C.setTP(T, 100.0)
            C.printTPD()
            
    print('$'*55)
    C.newDE(D=1,E=50.0)
    C.printTPD()
    
    C.constP_newH( 80.0 )
    C.printTPD()
    C.constH_newP( 500.0 )
    C.printTPD()
    
    C.setPH( 300., 75.0 )
    C.printTPD()
    
    C.constS_newP(P=1000.0)
    C.printTPD()
    
    C.setTD(T=500.0,D=5.0)
    C.printTPD()
    
    C.setPD(P=800.0,D=4.0)
    C.printTPD()
    
    print('------- Psat Test ----------')
    Psat = C.getSatP( T=200.0 )
    print( 'Psat = %g psia'%Psat )
    C.setProps(Q=Q_LIQUID, P=Psat)
    C.printTPD()
    C.setProps(Q=Q_GAS, P=Psat)
    C.printTPD()
    print()
    print('Psat=%g, rhoLiq=%g, rhoGas==%g'%C.getSatPandDens(T=190.0))
    C.printTPD()
    
    Tsat = C.getSatT( P=150.0 )
    print('Tsat = %g degR'%Tsat)
    C.setProps(Q=Q_LIQUID, T=Tsat)
    C.printTPD()
    C.setProps(Q=Q_GAS, T=Tsat)
    C.printTPD()
    
    print()
    print('Tsat=%g, rhoLiq=%g, rhoGas==%g'%C.getSatTandDens( P=200.0 ))
    C.printTPD()
    C.getSatT( P=196.1 )
    C.printTPD()
    
    print('_'*55)
    Tsat = C.getSatT( P=14.7 )
    print('Tnpb =',C.Tnbp,'   Tsat @1atm =',Tsat)
    #C.setTP(C.Tnbp, 14.7)
    for P in [5.0, 14.7, 50.0]:
        dH = C.dHvap( P=P )
        print('dH=%g at P=%g'%(dH,P) )
    
    
    