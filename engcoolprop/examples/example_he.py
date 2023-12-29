from engcoolprop.ec_fluid import EC_Fluid

# Print state point ( using 298.15K and 1 bar )
ec = EC_Fluid(symbol="Helium", T=298.15*1.8,P=14.5038 ) # T=degR, P=psia

ec.printProps() # Print state point at given T,P

