from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid
from engcoolprop.ec_fluid import ( Peng_fromSI, TSI_fromEng )
from CoolProp.CoolProp import PropsSI

# Print state point
ec_inc = EC_Incomp_Fluid(symbol="DowQ", T=630.0,P=100.0 ) # T=degR, P=psia

# print( "NOTE: This is EXPECTED TO FAIL!!!")
# ec_inc.setProps(T=500., D=0.1)

ec_inc.printProps() # Print state point at given T,P

t_range = (ec_inc.Tmax - ec_inc.Tmin) / 1.1
N = 9
tL = [ec_inc.Tmin + (i/N)*t_range for i in range(N)]

for T in tL:
    Psat = Peng_fromSI( PropsSI('P','T', TSI_fromEng( T ) ,'Q',0, ec_inc.fluid) )
    print( 'T = %g'%T, 'Psat = %f'%Psat)

Psat_max = Peng_fromSI( PropsSI('P','T', TSI_fromEng( ec_inc.Tmax ) ,'Q',0, ec_inc.fluid) )
print( 'Tmax = %g'%ec_inc.Tmax, 'Psat_max = %g'%Psat_max)
