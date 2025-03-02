from engcoolprop.ec_incomp_soln import EC_Incomp_Soln

# Create incompressible soln object at T=500 degR, P=500 psia and max pressure = 5000 psia
ec_soln = EC_Incomp_Soln(symbol="MEG-30%", T=500, P=500, Pmax=5000) # T=degR, P=psia

# OR... After ec_soln has been crated
# ec_soln.setTP( 500, 500)

# Print state point
ec_soln.printProps()
