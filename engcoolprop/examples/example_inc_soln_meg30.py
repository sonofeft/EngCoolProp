from engcoolprop.ec_incomp_soln import EC_Incomp_Soln

# Create incompressible solution object. (without specifying state point)
ec_soln = EC_Incomp_Soln(symbol="MEG-30%", auto_fix_value_errors=True, show_warnings=2)

# Print state point
ec_soln.printProps() # Print state point

print( '\n\n\n\n' )

ec_soln.printTPD()

ec_soln.printTransport()


print( '='*66 )
# printSIUnits
from engcoolprop.ec_incomp_soln import EC_Incomp_Soln

# Create incompressible solution object. (without specifying state point)
ec_soln = EC_Incomp_Soln(symbol="MEG-30%", auto_fix_value_errors=True, show_warnings=0)

# Print state point
ec_soln.printSIProps() # Print state point with SI units
