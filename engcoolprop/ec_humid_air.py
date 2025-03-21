""" 
Humid Air Properties
http://www.coolprop.org/fluid_properties/HumidAir.html#table-of-inputs-outputs-to-hapropssi

"""
import traceback
import CoolProp.CoolProp as CP
from engcoolprop.humid_air_params import (preferred_list, param_synonymD, param_eng_unitsD, 
                                          param_descD, param_si_unitsD, preferred_nameD, 
                                          input_set)

from engcoolprop.conv_funcs import (CPeng_fromSI ,  CondEng_fromSI ,   DSI_fromEng ,  CPSI_fromEng,
                                  Deng_fromSI ,   PSI_fromEng, Peng_fromSI, EchoInput,
                                  SSI_fromEng ,  Seng_fromSI ,  TSI_fromEng ,  
                                  Teng_fromSI ,  UHSI_fromEng ,  UHeng_fromSI ,   Veng_fromSI,
                                  CondSI_fromEng, VSI_fromEng  )
from engcoolprop.banner import banner

def Voleng_fromSI( Vol ): # 1 lbm/ft^3 = 16.01843417 kg/m^3
    return Vol * 16.01843417

def VolSI_fromEng( Vol ): # 1 lbm/ft^3 = 16.01843417 kg/m^3
    return Vol / 16.01843417

PSIA_PER_ATM = 14.6959 # used to convert between psia and atm

DEG_F_R_OFFSET = 459.67


print_orderL = ["Tdb","WetBulb","DewPoint","P","P_w","",
"Vda","Vha","","cp","cp_ha","CV","CVha","","Hda","Hha","Sda","Sha",
"",     "Cond","Visc","","RelHum","HumRat","Y","Z"]


# map Eng to SI conversion functions for each of the fluid properties
toSI_callD = {} # index=eng, value=conversion func (e.g. SSI_fromEng)
toSI_callD['Conductivity'] = CondSI_fromEng
toSI_callD['cp'] = CPSI_fromEng
toSI_callD['cp_ha'] = CPSI_fromEng
toSI_callD['CV'] = CPSI_fromEng
toSI_callD['CVha'] = CPSI_fromEng
toSI_callD['DewPoint'] = TSI_fromEng
toSI_callD['Hda'] = UHSI_fromEng
toSI_callD['Hha'] = UHSI_fromEng
toSI_callD['HumRat'] = EchoInput
toSI_callD['P'] = PSI_fromEng
toSI_callD['P_w'] = PSI_fromEng
toSI_callD['RelHum'] = EchoInput
toSI_callD['Sda'] = SSI_fromEng
toSI_callD['Sha'] = SSI_fromEng
toSI_callD['Tdb'] = TSI_fromEng
toSI_callD['Vda'] = VolSI_fromEng
toSI_callD['Vha'] = VolSI_fromEng
toSI_callD['Visc'] = VSI_fromEng
toSI_callD['WetBulb'] = TSI_fromEng
toSI_callD['Y'] = EchoInput
toSI_callD['Z'] = EchoInput

toSI_callD['Cond'] = CondSI_fromEng

# need to add to toSI_callD with all synonyms of preferred names
for pref_name in preferred_list:
    syn_set = param_synonymD[pref_name]
    for syn in syn_set:
        if syn not in toSI_callD:
            toSI_callD[syn] = toSI_callD[pref_name]

# print( 'toSI_callD.keys =', sorted(toSI_callD.keys(), key=str.lower))
# ==========================================================
# map SI to Eng conversion functions for each of the fluid properties
toEng_callD = {} # index=si char, value=conversion func (e.g. Seng_fromSI)

toEng_callD['Conductivity'] = CondEng_fromSI
toEng_callD['cp'] = CPeng_fromSI
toEng_callD['cp_ha'] = CPeng_fromSI
toEng_callD['CV'] = CPeng_fromSI
toEng_callD['CVha'] = CPeng_fromSI
toEng_callD['DewPoint'] = Teng_fromSI
toEng_callD['Hda'] = UHeng_fromSI
toEng_callD['Hha'] = UHeng_fromSI
toEng_callD['HumRat'] = EchoInput
toEng_callD['P'] = Peng_fromSI
toEng_callD['P_w'] = Peng_fromSI
toEng_callD['RelHum'] = EchoInput
toEng_callD['Sda'] = Seng_fromSI
toEng_callD['Sha'] = Seng_fromSI
toEng_callD['Tdb'] = Teng_fromSI
toEng_callD['Vda'] = Voleng_fromSI
toEng_callD['Vha'] = Voleng_fromSI
toEng_callD['Visc'] = Veng_fromSI
toEng_callD['WetBulb'] = Teng_fromSI
toEng_callD['Y'] = EchoInput
toEng_callD['Z'] = EchoInput

toEng_callD['Cond'] = CondEng_fromSI

# may need to show viscosity in lbf-s/ft^2 ???
def convert_viscosity_to_lbf_s_per_ft2(lbm_per_s_per_ft):
    # Conversion factor
    conversion_factor = 0.03108095
    # Convert to lbf-s/ft^2
    lbf_s_per_ft2 = lbm_per_s_per_ft * conversion_factor
    return lbf_s_per_ft2

    # Example usage
    # viscosity_lbm_per_s_per_ft = 1.23324e-5
    # viscosity_lbf_s_per_ft2 = convert_viscosity_to_lbf_s_per_ft2(viscosity_lbm_per_s_per_ft)
    # print(f"Dynamic Viscosity: {viscosity_lbf_s_per_ft2} lbf-s/ft²")


class EC_Humid_Air(object):
    
    def __init__(self,  **kwargs):

        self.si_propD = {} # key:preferred name, value: value in SI units
        self.eng_propD = {} # key:preferred name, value: value in Engineering units
        self.setProps( **kwargs )

    def setProps(self, **kwargs):
        """
        Calculates various properties of humid air using CoolProp's HAPropsSI function.

        Args:
        - kwargs (dict): Dictionary with three keys corresponding to the input properties and their values.

        Note the values in kwargs are in Engineering units.

        Results in object properties corresponding to state properties in engineering units.
        """
        # Allow user to call setProps with none or partial parameters
        if len(kwargs) == 0:
            kwargs = {'T':536.4 ,'P':PSIA_PER_ATM, 'RelHum':0.5}
        
        # Allow any temperature to be input in degF
        self.eng_inputD = {} # key:CoolProp input name: value=value in Engineering units
        for k,v in kwargs.items():
            if k.endswith('degF'):
                self.eng_inputD[ k[:-4] ] = v + DEG_F_R_OFFSET
            else:
                self.eng_inputD[k] = v

        if len(self.eng_inputD) == 2 and 'P' not in kwargs:
            self.eng_inputD['P'] = PSIA_PER_ATM

        # Ensure the kwargs dictionary has exactly three items
        if len(self.eng_inputD) != 3:
            raise ValueError("Exactly three input properties must be provided")
        
        if 'P' not in self.eng_inputD:
            raise ValueError('"P" must be specified')

        self.si_inputD = {}
        for k,v in self.eng_inputD.items():
            self.si_inputD[k] = toSI_callD[k](v)

        # print( 'eng_inputD =', self.eng_inputD)
        # print( 'si_inputD =', self.si_inputD)

        # Extract keys and values from the self.si_inputD dictionary
        input_keys = list(self.si_inputD.keys())
        input_vals = list(self.si_inputD.values())

        # Iterate through the list of properties and calculate each one
        self.success = True
        for prop in preferred_list:
            try:
                # Place SI results into self.si_propD
                self.si_propD[prop] = CP.HAPropsSI(prop, input_keys[0], input_vals[0], input_keys[1], input_vals[1], input_keys[2], input_vals[2])
            except:
                self.success = False

                print( 'Failed to calc:', prop)
                self.si_propD[prop] = float('inf')

                tb_str = traceback.format_exc()
                if 'ValueError' in tb_str:
                    print( 'ValueError::', tb_str.split('ValueError')[-1])
                else:
                    print( tb_str )

                break # all props will be set to float('inf')
                
        # if not success, set everything to infinity
        if not self.success:
            for prop in preferred_list:
                self.si_propD[prop] = float('inf')


        # include shorter, non-CoolProp name for Conductivity
        self.si_propD['Cond'] = self.si_propD['Conductivity']
                

        # convert properties from SI to Engineering
        for k,v in self.si_propD.items():
            self.eng_propD[k] = toEng_callD[k](v)

        # include shorter, non-CoolProp name for Conductivity
        self.eng_propD['Cond'] = self.eng_propD['Conductivity']

        # assign Engineering properties to self
        for name, value in self.eng_propD.items():
            setattr(self, name, value)
            for syn in param_synonymD[name]:
                if syn != name:
                    setattr(self, syn, value)


        # if 'Y' in self.eng_propD and self.eng_propD['Y'] < 1.0:
        #     # Constants
        #     M_dry_air = 28.96  # Molecular weight of dry air in g/mol
        #     M_water_vapor = 18.0153  # Molecular weight of water vapor in g/mol
        #     Y = self.eng_propD['Y']

        #     MolWt = (1 - Y) * M_dry_air + Y * M_water_vapor
        #     self.eng_propD['MolWt'] = MolWt
        #     self.si_propD['MolWt'] = MolWt

    def get_eng_fmt_size(self):
        """Make a g format sized for largest number of property in propL"""
        # sL = ['%g'%self.eng_propD[name] for name in propL]
        max_len = max([len('%g'%self.eng_propD[name]) for name in preferred_list])
        return '%' + '%ig'%max_len, max_len

    def get_si_fmt_size(self):
        """Make a g format sized for largest number of property in propL"""
        # sL = ['%g'%self.eng_propD[name] for name in propL]
        max_len = max([len('%g'%self.si_propD[name]) for name in preferred_list])
        return '%' + '%ig'%max_len, max_len

    def printProps(self, eng_units=True):
        '''
        Print a multiline property summary with Engineering units.
        
        If eng_units=False, SI units will be printed.
        '''

        
        def dict_to_string(D):
            return ', '.join([f"{key}={value:g}" for key, value in D.items()])

        if eng_units:
            s = dict_to_string( self.eng_inputD )
            fmt, max_len = self.get_eng_fmt_size( )
            propD = self.eng_propD.copy()
            unitsD = param_eng_unitsD
            propD['Visc'] *= 1.0E5
        else:
            s = dict_to_string( self.si_inputD )
            fmt, max_len = self.get_si_fmt_size( )
            propD = self.si_propD.copy()
            unitsD = param_si_unitsD
        
        banner(  '-'*4 + ' State Point for Humid Air (%s) '%s +'-'*4, leftMargin=3 )

        for name in print_orderL:
            if name:
                if unitsD[name] == 'degR':
                    alt_units = '%.1f'%(propD[name] - DEG_F_R_OFFSET,)
                    alt_units = '(' + alt_units.strip() + ' degF)'
                elif unitsD[name] == 'psia':
                    alt_units = fmt%(propD[name]/PSIA_PER_ATM,)
                    alt_units = '(' + alt_units.strip() + ' atm)'
                else:
                    alt_units = ''
                print("%10s ="%name, fmt%propD[name], unitsD[name],
                      " :: %s"%param_descD[name], alt_units)
                
                # if eng_units and name=='Visc':
                #     print( '           =', 
                #           fmt%convert_viscosity_to_lbf_s_per_ft2(self.eng_propD['Visc']),
                #            'lbf-s/ft^2', " :: %s"%param_descD[name] )
            else:
                print()


    def print_input_params(self):
        """Print all the legal inputs for Humid Air"""
        banner(  '-'*22 + ' Humid Air Input Parameters ' +'-'*22, leftMargin=3 )
        for name in print_orderL:
            if name:
                if name == 'Cond':
                    desc = preferred_nameD['Conductivity']
                else:
                    desc = preferred_nameD[name]

                # for name, desc in preferred_nameD.items():
                if name in input_set:
                    syn_set = param_synonymD[name]-{name}
                    if syn_set:
                        print( '%10s'%name, '%23s'%param_eng_unitsD[name], '%-34s'%desc, '::AKA', syn_set)
                    else:
                        print( '%10s'%name, '%23s'%param_eng_unitsD[name], '%-34s'%desc)

    def print_output_params(self):
        """Print all the legal outputs for Humid Air"""
        banner(  '-'*22 + ' Humid Air Output Parameters ' +'-'*22, leftMargin=3 )
        for name in print_orderL:
            if name:
                if name == 'Cond':
                    desc = preferred_nameD['Conductivity']
                else:
                    desc = preferred_nameD[name]
            
                syn_set = param_synonymD[name]-{name}
                # if name=='Conductivity': 
                #     name='Cond'
                if syn_set:
                    print( '%10s'%name, '%23s'%param_eng_unitsD[name], '%-34s'%desc, '::AKA', syn_set)
                else:
                    print( '%10s'%name, '%23s'%param_eng_unitsD[name], '%-34s'%desc)

def dev_tests():
    ha = EC_Humid_Air( T=536.4 , RelHum=0.5 )

    # kL = sorted( ha.eng_propD.keys(), key=str.lower)
    # for k in kL:
    #     print( '%10s = %12g'%(k,ha.eng_propD[k]), '= %12G'%ha.si_propD[k])

    print()

    ha.printProps()
    # print('.'*66)
    # ha.printProps(eng_units=False)
    print( 'ha.Cond =', ha.Cond)
    print( 'ha.Visc =', ha.Visc)

    ha.print_input_params()
    ha.print_output_params()

    print( '-'*66 )
    HA = ha = EC_Humid_Air( TdegF=70, RelHum=0.5 )
    HA.printProps()

if __name__ == "__main__":

    dev_tests()
