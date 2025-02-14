"""
In the CoolProp module, properties available for incompressible fluids include:

Temperature (T)
Pressure (P)
Density (D)
Specific Heat Capacity (C)
Internal Energy (U)
Enthalpy (H)
Entropy (S)
Viscosity (V)
Thermal Conductivity (L)
Minimum Temperature (Tmin)
Maximum Temperature (Tmax)

These properties can be accessed using the PropsSI function. 
CoolProp supports a variety of incompressible fluids, including pure fluids and binary mixtures.

Working independent pairs: TP, HP, SP, DP
from chk_indep_props_incomp.py:
    ONLY WORKS FOR: DP, HP, PS, PT  <-- NEED Pressure <--

NOT working: DT, HD, HS, UP, UD

"""
import CoolProp.CoolProp as CP

# Example fluid: Water
fluid = 'INCOMP::Water'

# Specifying Temperature and Pressure
T = 298.15  # Temperature in K
P = 101325  # Pressure in Pa

density = CP.PropsSI('D', 'T', T, 'P', P, fluid)
specific_heat_capacity = CP.PropsSI('C', 'T', T, 'P', P, fluid)
internal_energy = CP.PropsSI('U', 'T', T, 'P', P, fluid)
enthalpy = CP.PropsSI('H', 'T', T, 'P', P, fluid)
entropy = CP.PropsSI('S', 'T', T, 'P', P, fluid)
viscosity = CP.PropsSI('V', 'T', T, 'P', P, fluid)
thermal_conductivity = CP.PropsSI('L', 'T', T, 'P', P, fluid)

print(f"Density: {density} kg/m³")
print(f"Specific Heat Capacity: {specific_heat_capacity} J/kg.K")
print(f"Internal Energy: {internal_energy} J/kg")
print(f"Enthalpy: {enthalpy} J/kg")
print(f"Entropy: {entropy} J/kg.K")
print(f"Viscosity: {viscosity} Pa.s")
print(f"Thermal Conductivity: {thermal_conductivity} W/m.K")

# Specifying Density and Temperature
D = density  # Density in kg/m³

# DT NOT supported
# pressure_from_DT = CP.PropsSI('P', 'D', D, 'T', T, fluid)
# print(f"Pressure from Density and Temperature: {pressure_from_DT} Pa")

# Specifying Enthalpy and Pressure
H = enthalpy  # Enthalpy in J/kg

temperature_from_HP = CP.PropsSI('T', 'H', H, 'P', P, fluid)
print(f"Temperature from Enthalpy and Pressure: {temperature_from_HP} K")

# Specifying Entropy and Pressure
S = entropy  # Entropy in J/kg.K

temperature_from_SP = CP.PropsSI('T', 'S', S, 'P', P, fluid)
print(f"Temperature from Entropy and Pressure: {temperature_from_SP} K")

# Minimum and Maximum Temperature
Tmin = CP.PropsSI('Tmin', fluid)
Tmax = CP.PropsSI('Tmax', fluid)

print(f"Minimum Temperature: {Tmin} K")
print(f"Maximum Temperature: {Tmax} K")

# =============================================================


# Specifying Density and Pressure
D = 997  # Density in kg/m³
P = 101325  # Pressure in Pa

temperature_from_DP = CP.PropsSI('T', 'D', D, 'P', P, fluid)
print(f"Temperature from Density and Pressure: {temperature_from_DP} K")

# Specifying Enthalpy and Density
H = 419000  # Enthalpy in J/kg

# HD NOT supported
# temperature_from_HD = CP.PropsSI('T', 'H', H, 'D', D, fluid)
# print(f"Temperature from Enthalpy and Density: {temperature_from_HD} K")

# Specifying Enthalpy and Entropy
S = 1500  # Entropy in J/kg.K

# HS NOT supported
# temperature_from_HS = CP.PropsSI('T', 'H', H, 'S', S, fluid)
# print(f"Temperature from Enthalpy and Entropy: {temperature_from_HS} K")

# Specifying Internal Energy and Pressure
U = 210000  # Internal Energy in J/kg

# UP NOT supported
# temperature_from_UP = CP.PropsSI('T', 'U', U, 'P', P, fluid)
# print(f"Temperature from Internal Energy and Pressure: {temperature_from_UP} K")

# UD NOT supported
# Specifying Internal Energy and Density
# temperature_from_UD = CP.PropsSI('T', 'U', U, 'D', D, fluid)
# print(f"Temperature from Internal Energy and Density: {temperature_from_UD} K")

