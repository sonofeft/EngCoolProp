""" 
Humid Air Properties
http://www.coolprop.org/fluid_properties/HumidAir.html#table-of-inputs-outputs-to-hapropssi

"""
import CoolProp.CoolProp as CP
from engcoolprop.humid_air_params import (preferred_list, param_synonymD, param_eng_unitsD, 
                                          param_descD, param_si_unitsD)

from engcoolprop.ec_fluid import (CPeng_fromSI ,  CondEng_fromSI ,   DSI_fromEng ,  CPSI_fromEng,
                                  Deng_fromSI ,   PSI_fromEng, Peng_fromSI, EchoInput,
                                  SSI_fromEng ,  Seng_fromSI ,  TSI_fromEng ,  
                                  Teng_fromSI ,  UHSI_fromEng ,  UHeng_fromSI ,   Veng_fromSI,
                                  CondSI_fromEng, VSI_fromEng  )


def Voleng_fromSI( Vol ): # 1 lbm/ft^3 = 16.01843417 kg/m^3
    return Vol * 16.01843417

def VolSI_fromEng( Vol ): # 1 lbm/ft^3 = 16.01843417 kg/m^3
    return Vol / 16.01843417


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

# need to add to toEng_callD with all synonyms for preferred names
for pref_name in preferred_list:
    syn_set = param_synonymD[pref_name]
    for syn in syn_set:
        if syn not in toSI_callD:
            toSI_callD[syn] = toSI_callD[pref_name]


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
                         The values in kwargs are in Engineering units

        Returns:
        - dict: A dictionary with property names as keys and their calculated values as values.
        """
        # Allow user to call setProps with none or partial parameters
        if len(kwargs) == 0:
            kwargs = {'T':536.4 ,'P':14.6959, 'RelHum':0.5}
        elif len(kwargs) == 2 and 'P' not in kwargs:
            kwargs['P'] = 14.6959

        # Ensure the kwargs dictionary has exactly three items
        if len(kwargs) != 3:
            raise ValueError("Exactly three input properties must be provided")
        
        if 'P' not in kwargs:
            raise ValueError('"P" must be specified')
        
        self.eng_inputD = kwargs.copy()

        self.si_inputD = {}
        for k,v in kwargs.items():
            self.si_inputD[k] = toSI_callD[k](v)


        # Extract keys and values from the self.si_inputD dictionary
        input_keys = list(self.si_inputD.keys())
        input_vals = list(self.si_inputD.values())

        # Iterate through the list of properties and calculate each one
        for prop in preferred_list:
            try:
                # Place SI results into self.si_propD
                self.si_propD[prop] = CP.HAPropsSI(prop, input_keys[0], input_vals[0], input_keys[1], input_vals[1], input_keys[2], input_vals[2])
            except:
                print( 'Failed to calc:', prop)
                self.si_propD[prop] = float('inf')
                

        # convert properties from SI to Engineering
        for k,v in self.si_propD.items():
            self.eng_propD[k] = toEng_callD[k](v)

        # include shorter, non-CoolProp name for Conductivity
        self.si_propD['Cond'] = self.si_propD['Conductivity']
        self.eng_propD['Cond'] = self.eng_propD['Conductivity']

        self.eng_propD['Visc'] *= 1.0E5

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
        return '%' + '%ig'%max_len

    def get_si_fmt_size(self):
        """Make a g format sized for largest number of property in propL"""
        # sL = ['%g'%self.eng_propD[name] for name in propL]
        max_len = max([len('%g'%self.si_propD[name]) for name in preferred_list])
        return '%' + '%ig'%max_len

    def printProps(self, eng_units=True):
        '''print a multiline property summary with units'''

        
        def dict_to_string(D):
            return ', '.join([f"{key}={value:g}" for key, value in D.items()])


        print_orderL = ["Tdb","DewPoint","WetBulb","P","P_w","",
        "Vda","Vha","","cp","cp_ha","CV","CVha","","Hda","Hha","Sda","Sha",
        "",     "Visc","Cond","","RelHum","HumRat","Y","Z"]

        if eng_units:
            s = dict_to_string( self.eng_inputD )
            fmt = self.get_eng_fmt_size( )
            propD = self.eng_propD
            unitsD = param_eng_unitsD
        else:
            s = dict_to_string( self.si_inputD )
            fmt = self.get_si_fmt_size( )
            propD = self.si_propD
            unitsD = param_si_unitsD
        print("============ State Point for Humid Air (%s) ============"%s)

        for name in print_orderL:
            if name:
                print("%10s ="%name, fmt%propD[name], unitsD[name],
                      " :: %s"%param_descD[name])
            else:
                print()


if __name__ == "__main__":

    ha = EC_Humid_Air( T=536.4 , RelHum=0.5 )

    kL = sorted( ha.eng_propD.keys(), key=str.lower)
    for k in kL:
        print( '%10s = %12g'%(k,ha.eng_propD[k]), '= %12G'%ha.si_propD[k])

    print()

    ha.printProps()
    print('.'*66)
    ha.printProps(eng_units=False)
