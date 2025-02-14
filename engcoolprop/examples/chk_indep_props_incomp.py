import CoolProp.CoolProp as CP

# Example fluid: Water
fluid = 'INCOMP::Water'

# Specifying Temperature and Pressure
T = 298.15  # Temperature in K
P = 101325  # Pressure in Pa

density = CP.PropsSI('D', 'T', T, 'P', P, fluid)
specific_heat_capacity = CP.PropsSI('C', 'T', T, 'P', P, fluid)
internal_energy = CP.PropsSI('U', 'T', T, 'P', P, fluid)
enthalpy = CP.PropsSI('H', 'T', T, 'P', P, fluid)
entropy = CP.PropsSI('S', 'T', T, 'P', P, fluid)
viscosity = CP.PropsSI('V', 'T', T, 'P', P, fluid)
thermal_conductivity = CP.PropsSI('L', 'T', T, 'P', P, fluid)

propertyD = {} # key:property letter (e.g. D, P, T), value:numeric SI value
propertyD['D'] = density
propertyD['C'] = specific_heat_capacity
propertyD['U'] = internal_energy
propertyD['H'] = enthalpy
propertyD['S'] = entropy
propertyD['V'] = viscosity
propertyD['L'] = thermal_conductivity
propertyD['T'] = T
propertyD['P'] = P


all_propL = sorted( propertyD.keys() )

call_tuplesD = {} # key:tuple indep vars (e.g. ("T","P")), value: set of dependent properties (e.g. "S", "H", etc.)

num_workedD = {} # key:independent property pair, value:number of successful calcs
# Iterate through all pairs of properties
for prop1 in all_propL:
    for prop2 in all_propL:
        if prop1 != prop2:
            for prop in all_propL:
                if prop not in [prop1, prop2]:
                    try:
                        val = CP.PropsSI( prop, prop1, propertyD[prop1], prop2, propertyD[prop2], fluid)
                        print( 'success: %s%s for property:%s = %g'%(prop1, prop2, prop, val))
                        
                        t = (prop1, prop2)
                        if t in call_tuplesD:
                            call_tuplesD[t].add( prop )
                        else:
                            call_tuplesD[t] = set( [prop] )

                        if prop2 > prop1:
                            pair = prop1 + prop2
                        else:
                            pair = prop2 + prop1

                        num_workedD[pair] = num_workedD.get( pair, 0) + 1
                    except:
                        print( 'FAILED: %s%s for property:%s'%(prop1, prop2, prop))

print( num_workedD )
for k,v in num_workedD.items():
    print( k, v)

"""
Temperature (T)
Pressure (P)
Density (D)
Specific Heat Capacity (C)
Internal Energy (U)
Enthalpy (H)
Entropy (S)
Viscosity (V)
Thermal Conductivity (L)
"""
print( 'call_tuplesD = {}  # key:tuple indep vars (e.g. ("T","P")), value: dependent property (e.g. "S", "H", etc.)')
for t,v in call_tuplesD.items():
    print( 'call_tuplesD[%s] = %s'%(t,v))
