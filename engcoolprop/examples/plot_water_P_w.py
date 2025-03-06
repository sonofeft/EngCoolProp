import matplotlib.pyplot as plt
from engcoolprop.ec_humid_air import EC_Humid_Air

HA = EC_Humid_Air()

for RelHum in [1.0, 0.7, 0.5, 0.3, 0.1]:

    tL = [500 + i for i in range(61)]
    hL = []
    for T in tL:
        HA.setProps( Tdb=T, RelHum=RelHum)
        hL.append( HA.P_w )

    plt.plot( tL, hL, label='RelHum=%g'%RelHum)
plt.grid( True )
plt.title( 'Humid Air Water Partial Pressure')
plt.xlabel( 'Dry Bulb Temperature (degR)')
plt.ylabel( 'Water Partial Pressure (psia)')
plt.legend( loc='best' )

plt.savefig( 'water_partial_pressure.png', dpi=200)
plt.show()
