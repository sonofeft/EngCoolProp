
from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid

# Create incompressible object. (without specifying state point)
ec_inc = EC_Incomp_Fluid(symbol="Water")

# Print state point
ec_inc.printProps() # Print state point at given T,P


from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid

# Create incompressible object at T=500 degR, P=500 psia and max pressure = 5000 psia
ec_inc = EC_Incomp_Fluid(symbol="Water", T=500, P=500, Pmax=5000)

# Print state point
ec_inc.printProps()

