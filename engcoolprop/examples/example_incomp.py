from engcoolprop.ec_fluid import EC_Fluid

print( "NOTE: This is EXPECTED TO FAIL!!!")
# Print state point
ec = EC_Fluid(symbol="DowQ", T=530.0,P=100.0 ) # T=degR, P=psia
ec.setProps(T=500., D=0.1)

ec.printProps() # Print state point at given T,P

