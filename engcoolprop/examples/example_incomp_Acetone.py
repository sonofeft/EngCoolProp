from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid

# Create incompressible object. (without specifying state point)
ec_inc = EC_Incomp_Fluid(symbol="Acetone" )

# Print state point
ec_inc.printProps()
