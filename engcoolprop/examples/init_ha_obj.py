from engcoolprop.ec_humid_air import EC_Humid_Air
from engcoolprop.banner import banner
        
def dict_to_string(D):
    return ', '.join([f"{key}={value:g}" for key, value in D.items()])

def make_ha_obj( init_str ):
    HA = eval( init_str )
    print( '%36s'%init_str, '#', dict_to_string(HA.eng_inputD), end='\n\n' )

banner(  '------ Results of different calls to EC_Humid_Air ------', leftMargin=10)

margin = ' ' * 16
print( margin, 'Can initialize with no input parameters')
make_ha_obj( 'EC_Humid_Air()' )

print( margin, 'Can omit pressure (P) and 1 atm will be assumed')
make_ha_obj( 'EC_Humid_Air(T=600, R=1)' )

print( margin, 'Can append "degF" to any temperature and degR will result')
make_ha_obj( 'EC_Humid_Air(TdegF=70, RelHum=.6)' )

print( margin, 'Can use all synonyms for inputs (e.g. T_db, Tdb, T)')
make_ha_obj( 'EC_Humid_Air(P=14.7, RH=.3, Tdb=530)' )


