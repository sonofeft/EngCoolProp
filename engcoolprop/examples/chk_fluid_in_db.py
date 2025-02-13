import CoolProp.CoolProp as CP

def is_heos_fluid(fluid_name):
    try:
        CP.AbstractState("HEOS", fluid_name)
        return True
    except ValueError:
        return False

# Example usage
for fluid_name in ['n2', 'N2', 'XXX', 'DOWQ']:
    if is_heos_fluid(fluid_name):
        print(f"{fluid_name} is available in HEOS fluids.")
    else:
        print(f"{fluid_name} is not available in HEOS fluids.")

print( '='*66)


def is_incomp_fluid(fluid_name):
    try:
        CP.AbstractState("INCOMP", fluid_name)
        return True
    except ValueError:
        return False

# Example usage
for fluid_name in [ 'DOWQ', 'DowQ']:
    if is_incomp_fluid(fluid_name):
        print(f"{fluid_name} is available in INCOMP fluids.")
    else:
        print(f"{fluid_name} is not available in INCOMP fluids.")
