import matplotlib.pyplot as plt
import numpy as np
from engcoolprop.ec_humid_air import EC_Humid_Air, DEG_F_R_OFFSET

HA = EC_Humid_Air()

relhumArr = np.linspace(.05, 1, 20)
drybulbL = [100, 90, 80, 70] # degF

for drybulb in drybulbL:

    dewptL = []
    wetbulbL = []
    for RelHum in relhumArr:
        HA.setProps( TdegF=drybulb, RelHum=RelHum)
        
        dewptL.append( HA.Tdp - DEG_F_R_OFFSET )
        wetbulbL.append( HA.Twb - DEG_F_R_OFFSET )

    plt.plot( dewptL, wetbulbL, label='Dry Bulb=%g F'%drybulb)
plt.grid( True )
plt.title( 'Humid Air Temperatures')
plt.xlabel( 'Dew Point (degF)')
plt.ylabel( 'Wet Bulb Temperature (degF)')
plt.legend( loc='best' )

plt.savefig( 'dew_point.png', dpi=200)
plt.show()
