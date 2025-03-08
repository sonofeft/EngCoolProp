
from engcoolprop.ec_humid_air import EC_Humid_Air
ha = EC_Humid_Air()
ha.print_input_params()

print()
ha.print_output_params()

print()
ha.printProps()


print()
ha.printProps( eng_units=False )
