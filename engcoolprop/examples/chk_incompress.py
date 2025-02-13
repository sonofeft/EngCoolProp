from CoolProp.CoolProp import PropsSI

#Specific heat capacity of Downtherm Q at 500 K and 1 atm
print( ' '*18, PropsSI('C','T',500,'P',101325,'INCOMP::DowQ') )
print( 'expected: Out[2]: 2288.1643758645673 ' )

#Density of Downtherm Q at 500 K and 1 atm.
print( ' '*18, PropsSI('D','T',500,'P',101325,'INCOMP::DowQ') )
print( 'expected: Out[3]: 809.0654931667821' )

#Round trip in thermodynamic properties
T_init = 500.0

P_init = 101325

D_init = PropsSI('D','T',T_init,'P',P_init,'INCOMP::DowQ')

S_init = PropsSI('S','D',D_init,'P',P_init,'INCOMP::DowQ')

H_init = PropsSI('H','S',S_init,'P',P_init,'INCOMP::DowQ')

T_init = PropsSI('T','H',H_init,'P',P_init,'INCOMP::DowQ')

print( ' '*18, 'T_init =', T_init )
print( 'expected: Out[10]: 500.0000000000001' )

#Saturation pressure of Downtherm Q at 500 K
print( ' '*18, PropsSI('P','T',500,'Q',0,'INCOMP::DowQ') )
print( 'expected: Out[11]: 38091.37403658103' )

#Minimum temperature for Downtherm Q
print( ' '*18, PropsSI('Tmin','T',0,'P',0,'INCOMP::DowQ') )
print( 'expected: Out[12]: 238.15' )

#Maximum temperature for Downtherm Q
print( ' '*18, PropsSI('Tmax','T',0,'P',0,'INCOMP::DowQ') )
print( 'expected: Out[13]: 633.15' )
