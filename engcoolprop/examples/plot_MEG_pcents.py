import matplotlib.pyplot as plt
import numpy as np
from engcoolprop.ec_incomp_soln import EC_Incomp_Soln

P = 1000
for pcent in [60, 40, 20, 0]:

    # Create incompressible solution object. (without specifying state point)
    ec_soln = EC_Incomp_Soln(symbol="MEG-%i%%"%pcent, auto_fix_value_errors=True, show_warnings=0)

    tArr = np.linspace(ec_soln.Tmin, ec_soln.Tmax, 50)

    densL = []
    for T in tArr:
        ec_soln.setTP( T, P)
        densL.append( ec_soln.D )

    plt.plot( tArr, densL, label="MEG-%i%%"%pcent)

plt.grid( True )
plt.title( 'MEG Solution Densities')
plt.xlabel( 'Temperature (degR)')
plt.ylabel( 'Density (lbm/cuft)')
plt.legend( loc='best' )

plt.savefig( 'MEG_pcent_D.png', dpi=200)
plt.show()
