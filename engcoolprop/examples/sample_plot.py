import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots.Common import BasePlot, PropertyDict, BaseDimension
from engcoolprop.ec_fluid import EC_Fluid

symbol = 'N2'

class ENGunits(PropertyDict):
    def __init__(self):
        self._D = BaseDimension(add_SI=0.0, mul_SI=1./16.01843417, off_SI=0.0, label='Density',       symbol=u'd', unit=u'lbm/ft3')
        self._H = BaseDimension(add_SI=0.0, mul_SI=0.000429923, off_SI=0.0, label='Specific Enthalpy',symbol=u'h', unit=u'BTU/lbm')
        self._P = BaseDimension(add_SI=0.0, mul_SI=1.0/6894.76, off_SI=0.0, label='Pressure',         symbol=u'p', unit=u'psia')
        self._S = BaseDimension(add_SI=0.0, mul_SI=0.0002388461,off_SI=0.0, label='Specific Entropy', symbol=u's', unit=u'BTU/lbm/R')
        self._T = BaseDimension(add_SI=0.0, mul_SI=1.8,         off_SI=0.0, label='Temperature',      symbol=u'T', unit=u'R')
        self._U = BaseDimension(add_SI=0.0, mul_SI=0.000429923, off_SI=0.0, label='Specific Internal Energy', symbol=u'u', unit=u'BTU/lbm')
        self._Q = BaseDimension(add_SI=0.0, mul_SI=1.0,         off_SI=0.0, label='Vapour Quality',   symbol=u'x', unit=u'')
        
BasePlot.UNIT_SYSTEMS['ENG'] = ENGunits() # ['EUR','KSI','SI']

ec = EC_Fluid(symbol)
ec.printCriticalProps()

# Plot Types: 'PS', 'PT', 'HS', 'TS', 'PD', 'TD', 'PH'
# tp_limits:  'NONE', 'DEF','ACHP','ORC'
plot = PropertyPlot(symbol, 'TD', unit_system='ENG', tp_limits='DEF')
plot.calc_isolines(CP.iHmass, num=11)
plot.calc_isolines(CP.iQ, num=11)
plot.title('%s (%s)'%(ec.name, ec.symbol))
plot.grid()
plot.show()

