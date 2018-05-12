from engcoolprop.ec_fluid import EC_Fluid

ec = EC_Fluid(symbol="AR")
print('Properties for: %s (%s)'%(ec.name, ec.symbol) )
for T in range(460, 570, 10):
    ec.setProps(T=T, P=14.7)
    ec.printTPD()

