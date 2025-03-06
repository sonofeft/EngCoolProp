import matplotlib.pyplot as plt
from engcoolprop.ec_humid_air import EC_Humid_Air

HA = EC_Humid_Air()

for RelHum in [1.0, 0.7, 0.5, 0.3, 0.1]:

    tL = [500 + i for i in range(61)]
    hL = []
    for T in tL:
        HA.setProps( Tdb=T, RelHum=RelHum)
        hL.append( HA.WetBulb )

    plt.plot( tL, hL, label='RelHum=%g'%RelHum)
plt.grid( True )
plt.title( 'Humid Air Wet Bulb Temperature')
plt.xlabel( 'Dry Bulb Temperature (degR)')
plt.ylabel( 'Wet Bulb Temperature (degR)')
plt.legend( loc='best' )

plt.savefig( 'wetbulb.png', dpi=200)
plt.show()
