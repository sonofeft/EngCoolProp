
input_set = set() # set of input params

output_set = set() # set of output params

param_si_unitsD = {} # key:name, value:SI units
param_eng_unitsD= {} # key:name, value:Engineering units

param_descD = {} # key:name, value:description
param_synonymD = {} # key:name, value: set of all synonyms

preferred_nameD = {} # key:preferred name, value:description


# ============= input_set ===================
input_set = {'DewPoint', 'Sha', 'T_dp', 'H', 'P', 'W', 'psi_w', 'R', 'Tdb', 'Entropy', 'RelHum', 
             'WetBulb', 'Enthalpy', 'T_db', 'Hda', 'D', 'T', 'V', 'B', 'Y', 'Vha', 'Omega', 'Sda', 
             'S', 'Hha', 'HumRat', 'Twb', 'P_w', 'RH', 'Vda', 'T_wb', 'Tdp'}

input_list = sorted( input_set, key=str.lower)

# ============= output_set ===================
output_set = {'DewPoint', 'Visc', 'mu', 'Sha', 'Z', 'T_dp', 'M', 'cv_ha', 'H', 'P', 'W', 'psi_w', 
              'CV', 'R', 'Tdb', 'Entropy', 'RelHum', 'CVha', 'WetBulb', 'Enthalpy', 'k', 'T_db', 
              'Hda', 'D', 'cp', 'T', 'V', 'B', 'Y', 'Vha', 'Omega', 'Sda', 'S', 'Hha', 'HumRat', 
              'Twb', 'Conductivity', 'P_w', 'RH', 'Cha', 'C', 'Vda', 'T_wb', 'cp_ha', 'Tdp', 'K'}

output_list = sorted( output_set, key=str.lower)

# ============= param_si_unitsD ===================
param_si_unitsD["B"] = "K"
param_si_unitsD["Twb"] = "K"
param_si_unitsD["T_wb"] = "K"
param_si_unitsD["WetBulb"] = "K"
param_si_unitsD["C"] = "J/kg dry air/K"
param_si_unitsD["cp"] = "J/kg dry air/K"
param_si_unitsD["Cha"] = "J/kg humid air/K"
param_si_unitsD["cp_ha"] = "J/kg humid air/K"
param_si_unitsD["CV"] = "J/kg dry air/K"
param_si_unitsD["CVha"] = "J/kg humid air/K"
param_si_unitsD["cv_ha"] = "J/kg humid air/K"
param_si_unitsD["D"] = "K"
param_si_unitsD["Tdp"] = "K"
param_si_unitsD["DewPoint"] = "K"
param_si_unitsD["T_dp"] = "K"
param_si_unitsD["H"] = "J/kg dry air"
param_si_unitsD["Hda"] = "J/kg dry air"
param_si_unitsD["Enthalpy"] = "J/kg dry air"
param_si_unitsD["Hha"] = "J/kg humid air"
param_si_unitsD["K"] = "W/m/K"
param_si_unitsD["k"] = "W/m/K"
param_si_unitsD["Conductivity"] = "W/m/K"
param_si_unitsD["M"] = "Pa-s"
param_si_unitsD["Visc"] = "Pa-s"
param_si_unitsD["mu"] = "Pa-s"
param_si_unitsD["psi_w"] = "mol water/mol humid air"
param_si_unitsD["Y"] = "mol water/mol humid air"
param_si_unitsD["P"] = "Pa"
param_si_unitsD["P_w"] = "Pa"
param_si_unitsD["R"] = ""
param_si_unitsD["RH"] = ""
param_si_unitsD["RelHum"] = ""
param_si_unitsD["S"] = "J/kg dry air/K"
param_si_unitsD["Sda"] = "J/kg dry air/K"
param_si_unitsD["Entropy"] = "J/kg dry air/K"
param_si_unitsD["Sha"] = "J/kg humid air/K"
param_si_unitsD["T"] = "K"
param_si_unitsD["Tdb"] = "K"
param_si_unitsD["T_db"] = "K"
param_si_unitsD["V"] = "m^3 /kg dry air"
param_si_unitsD["Vda"] = "m^3 /kg dry air"
param_si_unitsD["Vha"] = "m^3 /kg humid air"
param_si_unitsD["W"] = "kg water/kg dry air"
param_si_unitsD["Omega"] = "kg water/kg dry air"
param_si_unitsD["HumRat"] = "kg water/kg dry air"
param_si_unitsD["Z"] = ""

# Non-CoolProp names
param_si_unitsD["MolWt"] = "g/gmole"
param_si_unitsD["Cond"] = "W/m/K"

# ============= param_eng_unitsD ===================
param_eng_unitsD["B"] = "degR"
param_eng_unitsD["Twb"] = "degR"
param_eng_unitsD["T_wb"] = "degR"
param_eng_unitsD["WetBulb"] = "degR"
param_eng_unitsD["C"] = "BTU/lbm dry air/degR"  
param_eng_unitsD["cp"] = "BTU/lbm dry air/degR"
param_eng_unitsD["Cha"] = "BTU/lbm humid air/degR"
param_eng_unitsD["cp_ha"] = "BTU/lbm humid air/degR"
param_eng_unitsD["CV"] = "BTU/lbm dry air/degR"
param_eng_unitsD["CVha"] = "BTU/lbm humid air/degR"
param_eng_unitsD["cv_ha"] = "BTU/lbm humid air/degR"
param_eng_unitsD["D"] = "degR"
param_eng_unitsD["Tdp"] = "degR"
param_eng_unitsD["DewPoint"] = "degR"
param_eng_unitsD["T_dp"] = "degR"
param_eng_unitsD["H"] = "BTU/lbm dry air"  
param_eng_unitsD["Hda"] = "BTU/lbm dry air"
param_eng_unitsD["Enthalpy"] = "BTU/lbm dry air"
param_eng_unitsD["Hha"] = "BTU/lbm humid air"
param_eng_unitsD["K"] = "BTU/ft-hr-R"
param_eng_unitsD["k"] = "BTU/ft-hr-R"
param_eng_unitsD["Conductivity"] = "BTU/ft-hr-R"
param_eng_unitsD["M"] = "[1.0E5 * lbm/ft-sec]"
param_eng_unitsD["Visc"] = "[1.0E5 * lbm/ft-sec]"
param_eng_unitsD["mu"] = "[1.0E5 * lbm/ft-sec]"
param_eng_unitsD["psi_w"] = "mol water/mol humid air"
param_eng_unitsD["Y"] = "mol water/mol humid air"
param_eng_unitsD["P"] = "psia"
param_eng_unitsD["P_w"] = "psia"
param_eng_unitsD["R"] = ""
param_eng_unitsD["RH"] = ""
param_eng_unitsD["RelHum"] = ""
param_eng_unitsD["S"] = "BTU/lbm dry air/degR"
param_eng_unitsD["Sda"] = "BTU/lbm dry air/degR"
param_eng_unitsD["Entropy"] = "BTU/lbm dry air/degR"
param_eng_unitsD["Sha"] = "BTU/lbm humid air/degR"
param_eng_unitsD["T"] = "degR"
param_eng_unitsD["Tdb"] = "degR"
param_eng_unitsD["T_db"] = "degR"
param_eng_unitsD["V"] = "ft^3/lbm dry air"
param_eng_unitsD["Vda"] = "ft^3/lbm dry air"
param_eng_unitsD["Vha"] = "ft^3/lbm humid air"
param_eng_unitsD["W"] = "lbm water/lbm dry air"
param_eng_unitsD["Omega"] = "lbm water/lbm dry air"
param_eng_unitsD["HumRat"] = "lbm water/lbm dry air"
param_eng_unitsD["Z"] = ""

# Non-CoolProp names
param_eng_unitsD["MolWt"] = "lbm/lbmmole"
param_eng_unitsD["Cond"] = "BTU/ft-hr-R"


# ============= param_descD ===================
param_descD["B"] = "Wet-Bulb Temperature"
param_descD["Twb"] = "Wet-Bulb Temperature"
param_descD["T_wb"] = "Wet-Bulb Temperature"
param_descD["WetBulb"] = "Wet-Bulb Temperature"
param_descD["C"] = "Mixture Cp per unit dry air"
param_descD["cp"] = "Mixture Cp per unit dry air"
param_descD["Cha"] = "Mixture Cp per unit humid air"
param_descD["cp_ha"] = "Mixture Cp per unit humid air"
param_descD["CV"] = "Mixture Cv per unit dry air"
param_descD["CVha"] = "Mixture Cv per unit humid air"
param_descD["cv_ha"] = "Mixture Cv per unit humid air"
param_descD["D"] = "Dew-Point Temperature"
param_descD["Tdp"] = "Dew-Point Temperature"
param_descD["DewPoint"] = "Dew-Point Temperature"
param_descD["T_dp"] = "Dew-Point Temperature"
param_descD["H"] = "Mixture enthalpy per dry air"
param_descD["Hda"] = "Mixture enthalpy per dry air"
param_descD["Enthalpy"] = "Mixture enthalpy per dry air"
param_descD["Hha"] = "Mixture enthalpy per humid air"
param_descD["K"] = "Mixture thermal conductivity"
param_descD["k"] = "Mixture thermal conductivity"
param_descD["Conductivity"] = "Mixture thermal conductivity"
param_descD["M"] = "Mixture viscosity"
param_descD["Visc"] = "Mixture viscosity"
param_descD["mu"] = "Mixture viscosity"
param_descD["psi_w"] = "Water mole fraction"
param_descD["Y"] = "Water mole fraction"
param_descD["P"] = "Pressure"
param_descD["P_w"] = "Partial pressure of water vapor"
param_descD["R"] = "Relative humidity in range [0, 1]"
param_descD["RH"] = "Relative humidity in range [0, 1]"
param_descD["RelHum"] = "Relative humidity in range [0, 1]"
param_descD["S"] = "Mixture entropy per unit dry air"
param_descD["Sda"] = "Mixture entropy per unit dry air"
param_descD["Entropy"] = "Mixture entropy per unit dry air"
param_descD["Sha"] = "Mixture entropy per unit humid air"
param_descD["T"] = "Dry-Bulb Temperature"
param_descD["Tdb"] = "Dry-Bulb Temperature"
param_descD["T_db"] = "Dry-Bulb Temperature"
param_descD["V"] = "Mixture volume per unit dry air"
param_descD["Vda"] = "Mixture volume per unit dry air"
param_descD["Vha"] = "Mixture volume per unit humid air"
param_descD["W"] = "Humidity Ratio"
param_descD["Omega"] = "Humidity Ratio"
param_descD["HumRat"] = "Humidity Ratio"
param_descD["Z"] = "Compressibility factor (Z=pv/(RT))"

# Non-CoolProp names
param_descD["MolWt"] = "Mixture Molecular Weight"
param_descD["Cond"] = "Mixture thermal conductivity"

# ============= param_synonymD ===================
param_synonymD["B"] = "{'Twb', 'B', 'T_wb', 'WetBulb'}"
param_synonymD["Twb"] = "{'Twb', 'B', 'T_wb', 'WetBulb'}"
param_synonymD["T_wb"] = "{'Twb', 'B', 'T_wb', 'WetBulb'}"
param_synonymD["WetBulb"] = "{'Twb', 'B', 'T_wb', 'WetBulb'}"
param_synonymD["C"] = "{'C', 'cp'}"
param_synonymD["cp"] = "{'C', 'cp'}"
param_synonymD["Cha"] = "{'cp_ha', 'Cha'}"
param_synonymD["cp_ha"] = "{'cp_ha', 'Cha'}"
param_synonymD["CV"] = "{'CV'}"
param_synonymD["CVha"] = "{'CVha', 'cv_ha'}"
param_synonymD["cv_ha"] = "{'CVha', 'cv_ha'}"
param_synonymD["D"] = "{'Tdp', 'T_dp', 'D', 'DewPoint'}"
param_synonymD["Tdp"] = "{'Tdp', 'T_dp', 'D', 'DewPoint'}"
param_synonymD["DewPoint"] = "{'Tdp', 'T_dp', 'D', 'DewPoint'}"
param_synonymD["T_dp"] = "{'Tdp', 'T_dp', 'D', 'DewPoint'}"
param_synonymD["H"] = "{'Hda', 'Enthalpy', 'H'}"
param_synonymD["Hda"] = "{'Hda', 'Enthalpy', 'H'}"
param_synonymD["Enthalpy"] = "{'Hda', 'Enthalpy', 'H'}"
param_synonymD["Hha"] = "{'Hha'}"
param_synonymD["K"] = "{'Conductivity', 'k', 'K'}"
param_synonymD["k"] = "{'Conductivity', 'k', 'K'}"
param_synonymD["Conductivity"] = "{'Conductivity', 'k', 'K'}"
param_synonymD["M"] = "{'mu', 'M', 'Visc'}"
param_synonymD["Visc"] = "{'mu', 'M', 'Visc'}"
param_synonymD["mu"] = "{'mu', 'M', 'Visc'}"
param_synonymD["psi_w"] = "{'Y', 'psi_w'}"
param_synonymD["Y"] = "{'Y', 'psi_w'}"
param_synonymD["P"] = "{'P'}"
param_synonymD["P_w"] = "{'P_w'}"
param_synonymD["R"] = "{'RH', 'R', 'RelHum'}"
param_synonymD["RH"] = "{'RH', 'R', 'RelHum'}"
param_synonymD["RelHum"] = "{'RH', 'R', 'RelHum'}"
param_synonymD["S"] = "{'Sda', 'Entropy', 'S'}"
param_synonymD["Sda"] = "{'Sda', 'Entropy', 'S'}"
param_synonymD["Entropy"] = "{'Sda', 'Entropy', 'S'}"
param_synonymD["Sha"] = "{'Sha'}"
param_synonymD["T"] = "{'T_db', 'Tdb', 'T'}"
param_synonymD["Tdb"] = "{'T_db', 'Tdb', 'T'}"
param_synonymD["T_db"] = "{'T_db', 'Tdb', 'T'}"
param_synonymD["V"] = "{'Vda', 'V'}"
param_synonymD["Vda"] = "{'Vda', 'V'}"
param_synonymD["Vha"] = "{'Vha'}"
param_synonymD["W"] = "{'W', 'Omega', 'HumRat'}"
param_synonymD["Omega"] = "{'W', 'Omega', 'HumRat'}"
param_synonymD["HumRat"] = "{'W', 'Omega', 'HumRat'}"
param_synonymD["Z"] = "{'Z'}"

param_synonymD["Cond"] = "{'Conductivity', 'k', 'K'}"

# ============= preferred_nameD ===================
preferred_nameD["WetBulb"] = "Wet-Bulb Temperature"
preferred_nameD["cp"]    = "Mixture Cp per unit dry air"
preferred_nameD["cp_ha"] = "Mixture Cp per unit humid air"
preferred_nameD["CV"]    = "Mixture Cv per unit dry air"
preferred_nameD["CVha"]  = "Mixture Cv per unit humid air"
preferred_nameD["DewPoint"] = "Dew-Point Temperature"
preferred_nameD["Hda"] = "Mixture enthalpy per dry air"
preferred_nameD["Hha"] = "Mixture enthalpy per humid air"
preferred_nameD["Conductivity"] = "Mixture thermal conductivity"
preferred_nameD["Visc"] = "Mixture viscosity"
preferred_nameD["Y"] = "Water mole fraction"
preferred_nameD["P"] = "Pressure"
preferred_nameD["P_w"] = "Partial pressure of water vapor"
preferred_nameD["RelHum"] = "Relative humidity in range [0, 1]"
preferred_nameD["Sda"] = "Mixture entropy per unit dry air"
preferred_nameD["Sha"] = "Mixture entropy per unit humid air"
preferred_nameD["Tdb"] = "Dry-Bulb Temperature"
preferred_nameD["Vda"] = "Mixture volume per unit dry air"
preferred_nameD["Vha"] = "Mixture volume per unit humid air"
preferred_nameD["HumRat"] = "Humidity Ratio"
preferred_nameD["Z"] = "Compressibility factor (Z=pv/(RT))"

preferred_list = sorted( [k for k in preferred_nameD.keys()], key=str.lower)
preferred_set  = set( preferred_list )

if __name__ == "__main__":

    # make sure that all the descriptions are covered by preferred names
    desc_set = set()
    for v in param_descD.values():
        desc_set.add( v )

    pref_desc_set = set()
    for v in preferred_nameD.values():
        pref_desc_set.add( v )

    print( 'len(desc_set)      =', len(desc_set))
    print( 'len(pref_desc_set) =', len(pref_desc_set))

