from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid

ec_inc = EC_Incomp_Fluid(symbol="DowQ" ) # T=degR, P=psia

# Print state point
ec_inc.printProps() # Print state point at given T,P
