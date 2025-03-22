
from engcoolprop.conv_funcs import (ASI_fromEng, CondSI_fromEng, CPSI_fromEng, 
                                    DSI_fromEng, PSI_fromEng, SSI_fromEng, 
                                    TSI_fromEng, UHSI_fromEng, VSI_fromEng)
from engcoolprop.conv_funcs import  EchoInput
from engcoolprop.parameter_units import si_unitsD

def rhoSI_fromRhoEng( rho ):
    """Return g/cm^3 from lbm/cuin"""
    return rho * 27.679905
    
# map SI from Engineering conversion functions for each of the fluid properties
toSI_callD = {} # index=eng prop name, value=conversion func (e.g. SSI_fromEng)

class SI_obj:
    def __init__(self, eng_obj, propL):
        """
        Create an object of SI units from one of Engineering units
        eng_obj = the object having properties in Engineering units
        propL = names of properties to be converted to SI and made properties of SI_obj
        """
        self.propL = propL
        self.max_g_len = 8 # find max %g string length

        self.si_valueD = {} # key:prop name, value:SI value

        for prop in propL:
            # get Engineering value from eng_obj
            val = getattr(eng_obj, prop)

            
            # Convert Engineering value into SI value
            if callable(val):# might be a method name, i.e. callable, so check for that
                si_val = toSI_callD[prop](val())
            else:                
                si_val = toSI_callD[prop](val)
            self.si_valueD[prop] = si_val

            # find longest %g string representing SI value
            try:
                si_val_str = '%g'%si_val
                self.max_g_len = max(self.max_g_len, len(si_val_str) )
            except:
                pass

        # create formats for g and s 
        self.g_fmt = '%' + '%ig'%self.max_g_len


    def __getattr__(self, name):
        """
        Override __getattr__ to return right justified string of the maximum length.

        Parameters:
        name (str): The name of the attribute to access.

        Returns:
        str: The right justified string of the maximum length.
        """
        if name in self.si_valueD:
            try:
                s = self.g_fmt%self.si_valueD[name]
            except:
                s = self.si_valueD[name].rjust(self.max_g_len)
            return s
        else:
            return '?'*self.max_g_len


    def show_si_values(self):
        for prop in self.propL:
            print( '%20s'%prop, getattr(self, prop), si_unitsD.get(prop, ''))

toSI_callD['Cond']     = CondSI_fromEng
toSI_callD['Condmax']  = CondSI_fromEng
toSI_callD['Condmin']  = CondSI_fromEng
toSI_callD['Cp']       = CPSI_fromEng
toSI_callD['Cpmax']    = CPSI_fromEng
toSI_callD['Cpmin']    = CPSI_fromEng
toSI_callD['Cv']       = CPSI_fromEng
toSI_callD['D']        = DSI_fromEng
toSI_callD['Dc']       = DSI_fromEng
toSI_callD['Dmax']     = DSI_fromEng
toSI_callD['Dmin']     = DSI_fromEng
toSI_callD['E']        = UHSI_fromEng
toSI_callD['Emax']     = UHSI_fromEng
toSI_callD['Emin']     = UHSI_fromEng
toSI_callD['fluid']    = EchoInput
toSI_callD['gamma']    = EchoInput # NOTE: gamma is a function
toSI_callD['good_nbp'] = EchoInput
toSI_callD['H']        = UHSI_fromEng
toSI_callD['Hmax']     = UHSI_fromEng
toSI_callD['Hmin']     = UHSI_fromEng
toSI_callD['name']     = EchoInput
toSI_callD['P']        = PSI_fromEng
toSI_callD['Pc']       = PSI_fromEng
toSI_callD['percentage']    = EchoInput
toSI_callD['pcent_max_str'] = EchoInput
toSI_callD['pcent_min_str'] = EchoInput
toSI_callD['Pinput']   = PSI_fromEng
toSI_callD['Pmax']     = PSI_fromEng
toSI_callD['Pmin']     = PSI_fromEng
toSI_callD['psat_is_supported'] = EchoInput
toSI_callD['Psat']     = PSI_fromEng
toSI_callD['Psat_max'] = PSI_fromEng
toSI_callD['Psat_min'] = PSI_fromEng
toSI_callD['Q']        = EchoInput
toSI_callD['rho']      = rhoSI_fromRhoEng # NOTE: this is g/cm^3 from kg/m^3
toSI_callD['rho_max']  = rhoSI_fromRhoEng # NOTE: this is g/cm^3 from kg/m^3
toSI_callD['rho_min']  = rhoSI_fromRhoEng # NOTE: this is g/cm^3 from kg/m^3
toSI_callD['S']        = SSI_fromEng
toSI_callD['Smax']     = SSI_fromEng
toSI_callD['Smin']     = SSI_fromEng
toSI_callD['sonicV']   = ASI_fromEng
toSI_callD['symbol']   = EchoInput
toSI_callD['T']        = TSI_fromEng
toSI_callD['T_freeze'] = TSI_fromEng
toSI_callD['Tc']       = TSI_fromEng
toSI_callD['Tmax']     = TSI_fromEng
toSI_callD['Tmin']     = TSI_fromEng
toSI_callD['Tnbp']     = TSI_fromEng
toSI_callD['Ttriple']  = TSI_fromEng
toSI_callD['Visc']     = VSI_fromEng
toSI_callD['Viscmax']  = VSI_fromEng
toSI_callD['Viscmin']  = VSI_fromEng
toSI_callD['WtMol']    = EchoInput
toSI_callD['Z']        = EchoInput


if __name__ == "__main__":

    from engcoolprop.ec_fluid import EC_Fluid
    eng_obj = EC_Fluid( symbol='H2O' )
    eng_obj.printProps()

    ec_fluidL = ['Cond', 'Cp', 'Cv', 'D', 'Dc', 'E', 'gamma', 'good_nbp', 'H', 'name', 
         'P', 'Pc', 'Q', 'S', 'sonicV', 'symbol', 'T', 'Tc', 'Tnbp', 'Ttriple', 
         'Visc', 'WtMol', 'Z']

    si_obj = SI_obj( eng_obj, ec_fluidL)
    si_obj.show_si_values()
