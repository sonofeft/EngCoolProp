from engcoolprop.ec_fluid import EC_Fluid

# Print state point
ec = EC_Fluid(symbol="N2", T=530.0,P=100.0 ) # T=degR, P=psia
ec.setProps(T=500., D=0.1)

ec.printProps() # Print state point at given T,P

