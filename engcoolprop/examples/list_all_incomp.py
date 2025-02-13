import CoolProp.CoolProp as CP
from CoolProp import __incompressibles_pure__, __incompressibles_solution__

def list_incompressible_fluids():
    # Retrieve the list of all available fluids
    incompressible_fluids = CP.get_global_param_string('incompressible_list_pure').split(',')
    return incompressible_fluids

def list_incompressible_solutions():
    # Retrieve the list of all available fluids
    incompressible_solutions = CP.get_global_param_string('incompressible_list_solution').split(',')
    return incompressible_solutions

if __name__ == "__main__":
    fluids = list_incompressible_fluids()
    print( fluids )

    solutions = list_incompressible_solutions()
    print( solutions )

    print( '.'*66 )
    print( '__incompressibles_pure__ =', __incompressibles_pure__)

    print( '.'*66 )
    print( '__incompressibles_solution__ =', __incompressibles_solution__)
