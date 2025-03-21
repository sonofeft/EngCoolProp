
def Aeng_fromSI( A ):
    return A * 3.28084
    
def ASI_fromEng( A ):
    return A / 3.28084

def CondEng_fromSI( Cond ):
    return Cond * 0.578176

def CondSI_fromEng( Cond ):
    return Cond / 0.578176

def CPeng_fromSI( S ):
    return S * 0.0002388461111111111

def CPSI_fromEng( S ):
    return S / 0.0002388461111111111

def Deng_fromSI( D ): # 1 lbm/ft^3 = 16.01843417 kg/m^3
    return D  / 16.01843417

def DSI_fromEng( D ):
    return D * 16.01843417

def EchoInput( x ):
    return x

def Peng_fromSI( P ):
    return P / 6894.76

def PSI_fromEng( P ):
    return P * 6894.76
    
def Seng_fromSI( S ):
    return S * 0.0002388461111111111

def SSI_fromEng( S ):
    return S / 0.0002388461111111111

def Teng_fromSI( T ):
    return T * 1.8

def TSI_fromEng( T ):
    return T / 1.8

def UHeng_fromSI( UH ): # 1 BTU/lbm = 2326 J/kg
    return UH / 2326.0

def UHSI_fromEng( UH ):
    return UH * 2326.0

def Veng_fromSI( V ):
    # AS.viscosity() returns Pa*s
    return V * 0.6719689751 # convert from Pa*s to lbm/ft/sec
    
def VSI_fromEng( V ):
    return V / 0.6719689751
