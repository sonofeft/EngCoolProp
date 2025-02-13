import CoolProp.CoolProp as CP
from CoolProp import __fluids__

def list_heos_fluids():
    heos_fluids = []
    for fluid in CP.FluidsList():
        try:
            CP.AbstractState("HEOS", fluid)
            heos_fluids.append(fluid)
        except ValueError:
            continue
    return heos_fluids

if __name__ == "__main__":

    # List all fluids in HEOS backend
    heos_fluids = list_heos_fluids()
    print("Fluids available in HEOS backend:")
    print( heos_fluids )

    print( '.'*66 )
    print( '__fluids__ =', __fluids__)

