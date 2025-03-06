import matplotlib.pyplot as plt
from engcoolprop.ec_humid_air import EC_Humid_Air

HA = EC_Humid_Air()

for RelHum in [1.0, 0.7, 0.5, 0.3, 0.1]:

    tL = [500 + i for i in range(61)]
    hL = []
    for T in tL:
        HA.setProps( Tdb=T, RelHum=RelHum)
        hL.append( HA.Hha )

    plt.plot( tL, hL, label='RelHum=%g'%RelHum)
plt.grid( True )
plt.title( 'Humid Air Enthalpy (BTU/lbm)')
plt.xlabel( 'Dry Bulb Temperature (degR)')
plt.ylabel( 'Humid Air Enthalpy (BTU/lbm)')
plt.legend( loc='best' )

plt.savefig( 'humid_air_H.png', dpi=200)
plt.show()
