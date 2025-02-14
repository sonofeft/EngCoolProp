from engcoolprop.ec_fluid import EC_Fluid

# Print state point
ec = EC_Fluid(symbol="N2", T=530.0,P=100.0 ) # T=degR, P=psia
ec.setProps(T=500., D=0.1)

ec.printProps() # Print state point at given T,P

import engcoolprop.ec_fluid
pL = dir( engcoolprop.ec_fluid )

print( 'from engcoolprop.ec_fluid import (', end=' ')         
for p in pL:
    if not p.startswith( '__'):
     print( p, ', ', end=' ')
print( ')' )
print()

print( engcoolprop.ec_fluid.call_tuplesD )
