import CoolProp.CoolProp as CP
import traceback

N2 = CP.AbstractState("HEOS", 'N2')

DowQ = CP.AbstractState("INCOMP", 'DowQ')
DowQ.update(CP.PT_INPUTS, 101325, 500)

n2L = dir( N2 )
dowqL = dir( DowQ )

# print( 'n2L =', n2L )
print()
# print( 'dowqL =', dowqL )

for d_str in dowqL:
    d = getattr( DowQ, d_str )
    if not d_str.startswith('__'):
        try:
            print( d_str, 'd() =', d() )
            print( '...........................................')
        except:
            tb_str = traceback.format_exc()
            if tb_str.find('not implement') < 0 :
                print( d_str, type(d))
                print( 'FAILED', d_str)
                print( tb_str )
                print( '...........................................')

# print( '... .molar_mass() ...')
# print( N2.molar_mass() )
# print( DowQ.molar_mass() )

