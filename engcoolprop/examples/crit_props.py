from engcoolprop.ec_fluid import EC_Fluid

ec = EC_Fluid(symbol="O2", T=530.0,P=100.0 ) # T=degR, P=psia

ec.printCriticalProps() # Print critical properties
