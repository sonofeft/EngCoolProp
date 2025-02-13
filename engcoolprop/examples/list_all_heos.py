import CoolProp.CoolProp as CP

def list_heos_fluids():
    heos_fluids = []
    for fluid in CP.FluidsList():
        try:
            CP.AbstractState("HEOS", fluid)
            heos_fluids.append(fluid)
        except ValueError:
            continue
    return heos_fluids

# List all fluids in HEOS backend
heos_fluids = list_heos_fluids()
print("Fluids available in HEOS backend:")
print( heos_fluids )
