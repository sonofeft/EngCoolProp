
"""
# =========== Use these import statements for conversion functions ============

from engcoolprop.conv_funcs import (Aeng_fromSI, CondEng_fromSI, CPeng_fromSI, 
                                    Deng_fromSI, Peng_fromSI, Seng_fromSI, 
                                    Teng_fromSI, UHeng_fromSI, Veng_fromSI)
from engcoolprop.conv_funcs import (ASI_fromEng, CondSI_fromEng, CPSI_fromEng, 
                                    DSI_fromEng, PSI_fromEng, SSI_fromEng, 
                                    TSI_fromEng, UHSI_fromEng, VSI_fromEng)
from engcoolprop.conv_funcs import  EchoInput


"""
def Aeng_fromSI( A ):
    """Return ft/sec from m/sec"""
    return A * 3.28084
    
def ASI_fromEng( A ):
    """Return m/sec from ft/sec"""
    return A / 3.28084

def CondEng_fromSI( Cond ):
    """Return BTU/ft-hr-R from W/m/K"""
    return Cond * 0.578176

def CondSI_fromEng( Cond ):
    """Return W/m/K from BTU/ft-hr-R"""
    return Cond / 0.578176

def CPeng_fromSI( S ):
    """Return BTU/lbm degR from J/kg/K"""
    return S * 0.0002388461111111111

def CPSI_fromEng( S ):
    """Return J/kg/K from BTU/lbm degR"""
    return S / 0.0002388461111111111

def Deng_fromSI( D ): # 1 lbm/ft^3 = 16.01843417 kg/m^3
    """Return lbm/cu ft from kg/m^3"""
    return D  / 16.01843417

def DSI_fromEng( D ):
    """Return kg/m^3 from lbm/cu ft"""
    return D * 16.01843417

def EchoInput( x ):
    """Return input unchanged (can be float, string, int, etc.)"""
    return x

def Peng_fromSI( P ):
    """Return psia from Pa"""
    return P / 6894.76

def PSI_fromEng( P ):
    """Return Pa from psia"""
    return P * 6894.76
    
def Seng_fromSI( S ):
    """Return BTU/lbm degR from J/kg/K"""
    return S * 0.0002388461111111111

def SSI_fromEng( S ):
    """Return J/kg/K from BTU/lbm degR"""
    return S / 0.0002388461111111111

def Teng_fromSI( T ):
    """Return degR from degK"""
    return T * 1.8

def TSI_fromEng( T ):
    """Return degK from degR"""
    return T / 1.8

def UHeng_fromSI( UH ): # 1 BTU/lbm = 2326 J/kg
    """Return BTU/lbm from J/kg"""
    return UH / 2326.0

def UHSI_fromEng( UH ):
    """Return J/kg from BTU/lbm"""
    return UH * 2326.0

def Veng_fromSI( V ):
    """Return lbm/ft/sec from Pa-s"""
    # AS.viscosity() returns Pa*s
    return V * 0.6719689751 # convert from Pa*s to lbm/ft/sec
    
def VSI_fromEng( V ):
    """Return Pa-s from lbm/ft/sec"""
    return V / 0.6719689751
