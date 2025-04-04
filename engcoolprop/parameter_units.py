

eng_unitsD = {} # key:engineering param name, value: Engineering units
eng_unitsD["A"] = "ft/sec"
eng_unitsD["Ac"] = "ft/sec"
eng_unitsD["C"] = "BTU/ft-hr-R"
eng_unitsD["Cc"] = "BTU/ft-hr-R"
eng_unitsD["Cond"] = "BTU/ft-hr-R"
eng_unitsD["Cp"] = "BTU/lbm degR"
eng_unitsD["Cv"] = "BTU/lbm degR"
eng_unitsD["D"] = "lbm/cu ft"
eng_unitsD["Dc"] = "lbm/cu ft"
eng_unitsD["Dmax"] = "lbm/cu ft"
eng_unitsD["Dmin"] = "lbm/cu ft"
eng_unitsD["E"] = "BTU/lbm"
eng_unitsD["Ec"] = "BTU/lbm"
eng_unitsD["g"] = "Cp/Cv (-)"
eng_unitsD["H"] = "BTU/lbm"
eng_unitsD["Hc"] = "BTU/lbm"
eng_unitsD["MW"] = "lbm/lbmmole"
eng_unitsD["P"] = "psia"
eng_unitsD["Pc"] = "psia"
eng_unitsD["Pref"] = "psia"
eng_unitsD["Psat"] = "psia"
eng_unitsD["Q"] = ""
eng_unitsD["rho"] = "lbm/cu in"
eng_unitsD["S"] = "BTU/lbm degR"
eng_unitsD["Sc"] = "BTU/lbm degR"
eng_unitsD["T"] = "degR"
eng_unitsD["Tc"] = "degR"
eng_unitsD["Tmax"] = "degR"
eng_unitsD["Tmin"] = "degR"
eng_unitsD["Tnbp"] = "degR"
eng_unitsD["Tref"] = "degR"
eng_unitsD["Tsat"] = "degR"
eng_unitsD["U"] = "BTU/lbm" # U is SI for E
eng_unitsD["V"] = "[1.0E5 * lbm/ft-sec]"
eng_unitsD["Vc"] = "[1.0E5 * lbm/ft-sec]"
eng_unitsD["Visc"] = "[1.0E5 * lbm/ft-sec]"
eng_unitsD["Z"] = ""
eng_unitsD["Zc"] = ""
# eng_unitsD["xxx"] = "xxx"
# eng_unitsD["xxx"] = "xxx"
# eng_unitsD["xxx"] = "xxx"
# eng_unitsD["xxx"] = "xxx"
# eng_unitsD["xxx"] = "xxx"
# eng_unitsD["xxx"] = "xxx"


# coolprop units from: http://www.coolprop.org/coolprop/HighLevelAPI.html#fluid-information
#  http://www.coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function
si_unitsD = {} # key:coolprop param name, value: SI units
si_unitsD["A"] = "m/s"
si_unitsD["sonicV"] = "m/s"
si_unitsD["ACENTRIC"] = ""
si_unitsD["acentric"] = ""
si_unitsD["ALPHA0"] = ""
si_unitsD["alpha0"] = ""
si_unitsD["ALPHAR"] = ""
si_unitsD["alphar"] = ""
si_unitsD["BVIRIAL"] = ""
si_unitsD["Bvirial"] = ""
si_unitsD["C"] = "J/kg/K"
si_unitsD["CONDUCTIVITY"] = "W/m/K"
si_unitsD["conductivity"] = "W/m/K"
si_unitsD["Cond"] = "W/m/K"
si_unitsD["CP0MASS"] = "J/kg/K"
si_unitsD["Cp0mass"] = "J/kg/K"
si_unitsD["CP0MOLAR"] = "J/mol/K"
si_unitsD["Cp0molar"] = "J/mol/K"
si_unitsD["CPMASS"] = "J/kg/K"
si_unitsD["Cpmass"] = "J/kg/K"
si_unitsD["Cp"] = "J/kg/K"
si_unitsD["CPMOLAR"] = "J/mol/K"
si_unitsD["Cpmolar"] = "J/mol/K"
si_unitsD["CVIRIAL"] = ""
si_unitsD["Cvirial"] = ""
si_unitsD["CVMASS"] = "J/kg/K"
si_unitsD["Cvmass"] = "J/kg/K"
si_unitsD["Cv"] = "J/kg/K"
si_unitsD["CVMOLAR"] = "J/mol/K"
si_unitsD["Cvmolar"] = "J/mol/K"
si_unitsD["D"] = "kg/m^3"
si_unitsD["Dc"] = "kg/m^3"
si_unitsD["D2ALPHA0_DDELTA2_CONSTTAU"] = ""
si_unitsD["d2alpha0_ddelta2_consttau"] = ""
si_unitsD["D3ALPHA0_DDELTA3_CONSTTAU"] = ""
si_unitsD["d3alpha0_ddelta3_consttau"] = ""
si_unitsD["DALPHA0_DDELTA_CONSTTAU"] = ""
si_unitsD["dalpha0_ddelta_consttau"] = ""
si_unitsD["DALPHA0_DTAU_CONSTDELTA"] = ""
si_unitsD["dalpha0_dtau_constdelta"] = ""
si_unitsD["DALPHAR_DDELTA_CONSTTAU"] = ""
si_unitsD["dalphar_ddelta_consttau"] = ""
si_unitsD["DALPHAR_DTAU_CONSTDELTA"] = ""
si_unitsD["dalphar_dtau_constdelta"] = ""
si_unitsD["DBVIRIAL_DT"] = ""
si_unitsD["dBvirial_dT"] = ""
si_unitsD["DCVIRIAL_DT"] = ""
si_unitsD["dCvirial_dT"] = ""
si_unitsD["DELTA"] = ""
si_unitsD["Delta"] = ""
si_unitsD["DIPOLE_MOMENT"] = "C m"
si_unitsD["dipole_moment"] = "C m"
si_unitsD["DMASS"] = "kg/m^3"
si_unitsD["Dmass"] = "kg/m^3"
si_unitsD["DMOLAR"] = "mol/m^3"
si_unitsD["Dmolar"] = "mol/m^3"
si_unitsD["FH"] = ""
si_unitsD["FRACTION_MAX"] = ""
si_unitsD["fraction_max"] = ""
si_unitsD["FRACTION_MIN"] = ""
si_unitsD["fraction_min"] = ""
si_unitsD["FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS"] = ""
si_unitsD["fundamental_derivative_of_gas_dynamics"] = ""
si_unitsD["gamma"] = ""
si_unitsD["G"] = "J/kg"
si_unitsD["GAS_CONSTANT"] = "J/mol/K"
si_unitsD["gas_constant"] = "J/mol/K"
si_unitsD["GMASS"] = "J/kg"
si_unitsD["Gmass"] = "J/kg"
si_unitsD["GMOLAR"] = "J/mol"
si_unitsD["Gmolar"] = "J/mol"
si_unitsD["GMOLAR_RESIDUAL"] = "J/mol/K"
si_unitsD["Gmolar_residual"] = "J/mol/K"
si_unitsD["GWP100"] = ""
si_unitsD["GWP20"] = ""
si_unitsD["GWP500"] = ""
si_unitsD["H"] = "J/kg"
si_unitsD["HELMHOLTZMASS"] = "J/kg"
si_unitsD["Helmholtzmass"] = "J/kg"
si_unitsD["HELMHOLTZMOLAR"] = "J/mol"
si_unitsD["Helmholtzmolar"] = "J/mol"
si_unitsD["HH"] = ""
si_unitsD["HMASS"] = "J/kg"
si_unitsD["Hmass"] = "J/kg"
si_unitsD["HMOLAR"] = "J/mol"
si_unitsD["Hmolar"] = "J/mol"
si_unitsD["HMOLAR_RESIDUAL"] = "J/mol/K"
si_unitsD["Hmolar_residual"] = "J/mol/K"
si_unitsD["I"] = "N/m"
si_unitsD["ISENTROPIC_EXPANSION_COEFFICIENT"] = ""
si_unitsD["isentropic_expansion_coefficient"] = ""
si_unitsD["ISOBARIC_EXPANSION_COEFFICIENT"] = "1/K"
si_unitsD["isobaric_expansion_coefficient"] = "1/K"
si_unitsD["ISOTHERMAL_COMPRESSIBILITY"] = "1/Pa"
si_unitsD["isothermal_compressibility"] = "1/Pa"
si_unitsD["L"] = "W/m/K"
si_unitsD["M"] = "kg/mol"
si_unitsD["WtMol"] = "kg/mol"
si_unitsD["MOLAR_MASS"] = "kg/mol"
si_unitsD["molar_mass"] = "kg/mol"
si_unitsD["MOLARMASS"] = "kg/mol"
si_unitsD["molarmass"] = "kg/mol"
si_unitsD["MOLEMASS"] = "kg/mol"
si_unitsD["molemass"] = "kg/mol"
si_unitsD["O"] = "J/kg/K"
si_unitsD["ODP"] = ""
si_unitsD["P"] = "Pa"
si_unitsD["Pc"] = "Pa"
si_unitsD["P_CRITICAL"] = "Pa"
si_unitsD["p_critical"] = "Pa"
si_unitsD["P_MAX"] = "Pa"
si_unitsD["P_max"] = "Pa"
si_unitsD["P_MIN"] = "Pa"
si_unitsD["P_min"] = "Pa"
si_unitsD["P_REDUCING"] = "Pa"
si_unitsD["p_reducing"] = "Pa"
si_unitsD["P_TRIPLE"] = "Pa"
si_unitsD["p_triple"] = "Pa"
si_unitsD["PCRIT"] = "Pa"
si_unitsD["Pcrit"] = "Pa"
si_unitsD["pcrit"] = "Pa"
si_unitsD["PH"] = ""
si_unitsD["PHASE"] = ""
si_unitsD["Phase"] = ""
si_unitsD["PIP"] = ""
si_unitsD["PMAX"] = "Pa"
si_unitsD["pmax"] = "Pa"
si_unitsD["PMIN"] = "Pa"
si_unitsD["pmin"] = "Pa"
si_unitsD["PRANDTL"] = ""
si_unitsD["Prandtl"] = ""
si_unitsD["PTRIPLE"] = "Pa"
si_unitsD["ptriple"] = "Pa"
si_unitsD["Q"] = "mol/mol"
si_unitsD["RHOCRIT"] = "kg/m^3"
si_unitsD["rhocrit"] = "kg/m^3"
si_unitsD["RHOMASS_CRITICAL"] = "kg/m^3"
si_unitsD["rhomass_critical"] = "kg/m^3"
si_unitsD["RHOMASS_REDUCING"] = "kg/m^3"
si_unitsD["rhomass_reducing"] = "kg/m^3"
si_unitsD["RHOMOLAR_CRITICAL"] = "mol/m^3"
si_unitsD["rhomolar_critical"] = "mol/m^3"
si_unitsD["RHOMOLAR_REDUCING"] = "mol/m^3"
si_unitsD["rhomolar_reducing"] = "mol/m^3"
si_unitsD["S"] = "J/kg/K"
si_unitsD["SMASS"] = "J/kg/K"
si_unitsD["Smass"] = "J/kg/K"
si_unitsD["SMOLAR"] = "J/mol/K"
si_unitsD["Smolar"] = "J/mol/K"
si_unitsD["SMOLAR_RESIDUAL"] = "J/mol/K"
si_unitsD["Smolar_residual"] = "J/mol/K"
si_unitsD["SPEED_OF_SOUND"] = "m/s"
si_unitsD["speed_of_sound"] = "m/s"
si_unitsD["SURFACE_TENSION"] = "N/m"
si_unitsD["surface_tension"] = "N/m"
si_unitsD["T"] = "K"
si_unitsD["Tc"] = "K"
si_unitsD["Tnbp"] = "K"
si_unitsD["Ttriple"] = "K"
si_unitsD["T_CRITICAL"] = "K"
si_unitsD["T_critical"] = "K"
si_unitsD["T_FREEZE"] = "K"
si_unitsD["T_freeze"] = "K"
si_unitsD["T_MAX"] = "K"
si_unitsD["T_max"] = "K"
si_unitsD["T_MIN"] = "K"
si_unitsD["T_min"] = "K"
si_unitsD["T_REDUCING"] = "K"
si_unitsD["T_reducing"] = "K"
si_unitsD["T_TRIPLE"] = "K"
si_unitsD["T_triple"] = "K"
si_unitsD["TAU"] = ""
si_unitsD["Tau"] = ""
si_unitsD["TCRIT"] = "K"
si_unitsD["Tcrit"] = "K"
si_unitsD["TMAX"] = "K"
si_unitsD["Tmax"] = "K"
si_unitsD["TMIN"] = "K"
si_unitsD["Tmin"] = "K"
si_unitsD["TTRIPLE"] = "K"
si_unitsD["Ttriple"] = "K"
si_unitsD["U"] = "J/kg"
si_unitsD["E"] = "J/kg"
si_unitsD["UMASS"] = "J/kg"
si_unitsD["Umass"] = "J/kg"
si_unitsD["UMOLAR"] = "J/mol"
si_unitsD["Umolar"] = "J/mol"
si_unitsD["V"] = "Pa s"
si_unitsD["Visc"] = "Pa s"
si_unitsD["VISCOSITY"] = "Pa s"
si_unitsD["viscosity"] = "Pa s"
si_unitsD["Z"] = ""
