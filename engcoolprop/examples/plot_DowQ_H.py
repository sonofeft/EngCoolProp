import matplotlib.pyplot as plt
import numpy as np
from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid

# Create incompressible solution object. (without specifying state point)
ec_inc = EC_Incomp_Fluid(symbol="DowQ", auto_fix_value_errors=True, show_warnings=0 ) # T=degR, P=psia

tArr = np.linspace(500, 800, 50)

# Use a few different pressures
for P in [10000, 5000, 2000, 15]:

    densL = []
    for T in tArr:
        ec_inc.setTP( T, P)
        densL.append( ec_inc.H )

    plt.plot( tArr, densL, label="P = %g"%P)

plt.grid( True )
plt.title( 'DowQ Enthalpy')
plt.xlabel( 'Temperature (degR)')
plt.ylabel( 'Enthalpy (BTU/lbm)')
plt.legend( loc='best' )

plt.savefig( 'DowQ_enthalpy.png', dpi=200)
plt.show()
