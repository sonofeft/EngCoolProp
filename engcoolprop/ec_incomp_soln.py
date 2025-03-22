#!/usr/bin/env python
# -*- coding: ascii -*-


r"""
CoolProp has a set of incompressible solutions (e.g. Brines and Solutions)
EC_Incomp_Soln is a wrapper of those solutions using Engineering Units.

To get a list of these solutions use::

    incompressible_list_solution = CP.get_global_param_string('incompressible_list_solution')

EngCoolProp uses units of primarily inch, lbm, lbf, sec, BTU (some use of ft and hour).::

    #:c     The following are the default units for each property of Incompressible Fluids.
    #:c 
    #:c     T = Temperature = degR
    #:c     P = Pressure = psia
    #:c     D = Density = lbm/cu ft
    #:c     rho = Density = lbm/cu inch
    #:c     E = Internal Energy = BTU/lbm
    #:c     H = Enthalpy = BTU/lbm
    #:c     S = Entropy = BTU/lbm degR
    #:c     Cp = Heat Capacity (const. P) = BTU/lbm degR
    #:c     V = Viscosity = 1.0E5 * lbm/ft-sec
    #:c     C = Thermal Conductivity = BTU/ft-hr-R

"""
import traceback
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP

from engcoolprop.iteration_utils import calc_T_freeze
from engcoolprop.InterpProp_scipy import get_density_interpolator
from engcoolprop.conv_funcs import (CPeng_fromSI ,  CondEng_fromSI ,   DSI_fromEng ,  
                                  Deng_fromSI ,   PSI_fromEng, Peng_fromSI,
                                  SSI_fromEng ,  Seng_fromSI ,  TSI_fromEng ,  
                                  Teng_fromSI ,  UHSI_fromEng ,  UHeng_fromSI ,   Veng_fromSI )

from engcoolprop.safe_get_property import safe_get_INCOMP_prop, is_frac_max_check_inclusive
from engcoolprop.utils import Same_g_len, parse_coolprop_mixture, incomp_pure_solnL, print_avoid_valerr
from engcoolprop.iteration_utils import find_T_from_Dterp
from engcoolprop.banner import banner
from engcoolprop.si_object import SI_obj

# how to make a list of all incompressible solutions in coolprop
# incomp_pure_solnL = CP.get_global_param_string('incompressible_list_solution').split(',')


class EC_Incomp_Soln(object):
    
    def __init__(self,symbol="MEG-20%", T=None ,P=None, Pmax=10000.0,
                 show_warnings=2, auto_fix_value_errors=True):
        """        Init generic Incompressible Solution

        Args:
            symbol (str): Name and percent mass. ("LiBr[0.23]" or "LiBr-23%") Defaults to "MEG-20%".
            T (None or float): Temperature degR (None sets T to (Tmin+Tmax)/2)
            P (None or float): Pressure psia. (None set P to Pmax/10)
            Pmax (float): Max expected pressure psia. Defaults to 10000.0.
            show_warnings (int): Sets warning level(0=None, 1=Only serious, 2=all). Defaults to 2.
            auto_fix_value_errors (bool): Action when ValueError occurs. 
                                        True=correct problem with a warning, 
                                        False throws Exception Defaults to False.
        """
        self.auto_fix_value_errors = auto_fix_value_errors

        self.basename, self.percentage = parse_coolprop_mixture(symbol)
        if self.basename is None:
            raise ValueError( 'Input symbol "" not recognized.'%symbol )
        
        # If only basename was given, create a proper symbol with mass fraction
        if symbol == self.basename:
            symbol = symbol + '-' + '%g'%( self.percentage ) + '%'
        
        # print( "basename=%s, percentage=%s"%(self.basename, self.percentage) )

        frac_min = PropsSI('fraction_min','INCOMP::%s'%self.basename)
        frac_max = PropsSI('fraction_max','INCOMP::%s'%self.basename)
        # print( "CoolProp DB frac_min=%s, frac_max=%s"%(frac_min, frac_max) )

        frac_max_hi = frac_max  * 1.0000000000002 # needed to handle CoolProp non-inclusive range check


        self.pcent_min_str = '%g'%( frac_min*100 ) + '%'
        self.pcent_max_str = '%g'%( frac_max*100 ) + '%'

        self.fraction = self.percentage / 100.0
        # print( "self.fraction=%g"%self.fraction )

        # Make sure that solution base name is in CoolProp DB
        if self.basename not in incomp_pure_solnL:
            raise ValueError( '"%s" is NOT in coolprop incompressible list\n%s'%(symbol, repr(incomp_pure_solnL) ) )

        if show_warnings>1 and auto_fix_value_errors:
            banner( 'NOTICE: any input violations on limits of T, D or mass fraction\nwill be automatically corrected (set to min or max).'+\
                   '\nTo change this behavior set "auto_fix_value_errors" to False'+\
                    '\nTo suppress this banner set "show_warnings" to 0 or 1')


        # ==========  Account for errors in CoolProp DB ======================
        # Feb 2025 non-inclusive fluids ['ExampleSolution', 'IceEA', 'IceNA', 'IcePG', 'MAM2', 'MKA2', 
        #                                'MMG2', 'MPG2', 'VMG', 'ZLC', 'ZMC' ]
        frac_max_is_inclusive = is_frac_max_check_inclusive( self.basename )

        frac_max_adjusted = frac_max # might be changed if not frac_max_is_inclusive

        if not frac_max_is_inclusive:
            
            frac_max_adjusted = frac_max - 0.00001

            if show_warnings and self.fraction >= frac_max:
                print( 'WARNING: CoolProp has a round-off problem with frac_max=%g for %s'%(frac_max, self.basename) )
                print( "   CoolProp DB frac_min=%s, frac_max=%s"%(frac_min, frac_max) )
                print( '   Changing frac_max from %g to %g'%(frac_max, frac_max_adjusted))


        # make sure that fraction is in range
        if self.fraction < frac_min:
            if self.auto_fix_value_errors:
                new_pcent = '%g'%( frac_min*100 ) + '%'
                if show_warnings:
                    print( 'WARNING: input percentage(%s) is too low. Reset to %s'%( self.percentage, new_pcent ))
                self.percentage = new_pcent
                self.fraction = frac_min
                symbol = self.basename + '-' + new_pcent
            else:
                print_avoid_valerr()
                raise ValueError('Input percentage(%s) is too low. min percent = %g'%( self.percentage, frac_min*100 ))

        elif self.fraction >= frac_max_adjusted and self.fraction <= frac_max_hi: 
            # Always fix the CoolProp problem of non-inclusive range check
            # Need to reduce self.fraction below frac_max with non-inclusive range check.
            if not frac_max_is_inclusive:
                new_pcent = '%.3f'%( frac_max_adjusted*100 ) + '%'

                if show_warnings:
                    print( 'WARNING: input percentage(%s) is adjusted for CoolProp issue. Reset to %s'%( self.percentage, new_pcent ))
                self.percentage = new_pcent
                self.fraction = frac_max
                symbol = self.basename + '-' + new_pcent

        elif self.fraction > frac_max:
            if self.auto_fix_value_errors:
                # Need to reduce self.fraction below frac_max with non-inclusive range check.
                if not frac_max_is_inclusive:
                    new_pcent = '%.3f'%( frac_max_adjusted*100 ) + '%'
                else:
                    new_pcent = '%g'%( frac_max*100 ) + '%'

                if show_warnings:
                    print( 'WARNING: input percentage(%s) is too high. Reset to %s'%( self.percentage, new_pcent ))
                self.percentage = new_pcent
                self.fraction = frac_max
                symbol = self.basename + '-' + new_pcent
            else:
                print_avoid_valerr()
                raise ValueError( 'Input percentage(%s) is too high. max percent = %g'%( self.percentage, frac_max*100 ) )


        
        self.symbol = symbol


        self.Pmax = Pmax # highest pressure considered in any iterative calcs (can still input P > Pmax)
        self.Pmax_si = PSI_fromEng( Pmax )
        self.Pmin = 0

        if P is None:
            P = int( Pmax / 10 )
        self.P = P

        self.show_warnings = show_warnings # some calcs will issue warning statements if show_warnings >0 severe>1

        self.fluid = 'INCOMP::%s'%symbol

        # get temperature limits (Note attempt to avoid CoolProp round-off problems)
        self.Tmin_si =   PropsSI('Tmin','T',0,'P',0,'INCOMP::%s'%symbol) * 1.0000000000000002
        self.Tmax_si =   PropsSI('Tmax','T',0,'P',0,'INCOMP::%s'%symbol) * 0.9999999999999999

        # print( "self.Tmin_si=%g degK, self.Tmax_si=%g degK"%(self.Tmin_si, self.Tmax_si) )

        try:
            self.T_freeze_si = PropsSI('T_freeze','INCOMP::%s'%self.symbol)
            # print( "T_freeze = %.1f degR"%Teng_fromSI(self.T_freeze_si) )
            
            if self.T_freeze_si > self.Tmin_si:
                if self.show_warnings:
                    print( 'NOTICE: Tmin=%.1f degR has been increased to T_freeze + 1 = %.1f degR'%\
                          ( Teng_fromSI(self.Tmin_si), Teng_fromSI(self.T_freeze_si + 5.0/9.0)) )
                self.Tmin_si = self.T_freeze_si + 5.0/9.0
        except:
            # If no freezing point found, set it to 0
            self.T_freeze_si = 0

        if self.T_freeze_si <= 0:
            if self.show_warnings:
                print('NOTICE: T_freeze NOT in CoolProp DB. Trying to find it by traceback Exception.')
            self.T_freeze_si = calc_T_freeze(self)
            if self.T_freeze_si > self.Tmin_si:
                if self.show_warnings:
                    print( 'NOTICE: Tmin=%.1f degR has been increased to T_freeze + 1 = %.1f degR'%\
                          ( Teng_fromSI(self.Tmin_si), Teng_fromSI(self.T_freeze_si + 5.0/9.0)) )
                self.Tmin_si = self.T_freeze_si  + 5.0/9.0

        self.T_freeze = Teng_fromSI(self.T_freeze_si)
        # print( "self.T_freeze=%.1f"%self.T_freeze )



        self.Tmin =  Teng_fromSI( self.Tmin_si  )
        self.Tmax =  Teng_fromSI( self.Tmax_si )
        # print( "Tmin=%.1f, Tmax=%.1f"%(self.Tmin, self.Tmax) )

        self.Tmid = (self.Tmin + self.Tmax) / 2.0

        self.Dterp = get_density_interpolator( self )

        self.check_visc_support()

        if T is None:
            T = int((self.Tmin + self.Tmax) / 2.0)
            # else: T = input T

        # if input T is in range, use it... otherwise
        if T<self.Tmin:
            if self.auto_fix_value_errors:                
                self.T = self.Tmin + 0.001
                if show_warnings:
                    print( 'NOTICE: input T too low. Changed from: %g to: %g degR'%(T, self.T) )
                    print( "        INCOMP: Tmin=%g, Tmax=%g"%(self.Tmin, self.Tmax) )
            else:
                print_avoid_valerr()
                raise ValueError("Input T too low. T=%g,  Tmin=%g"%(T, self.Tmin))
                
        elif T>self.Tmax:
            if self.auto_fix_value_errors:                
                self.T = self.Tmax - 0.001
                if show_warnings:
                    print( 'NOTICE: input T too high. Changed from: %g to: %g degR'%(T, self.T) )
                    print( "        INCOMP: Tmin=%g, Tmax=%g"%(self.Tmin, self.Tmax) )
            else:
                print_avoid_valerr()
                raise ValueError("Input T too high. T=%g,  Tmax=%g"%(T, self.Tmax))
                
        else:
            self.T = T

        self.calc_min_max_props()

        # set properties to input T and P
        # print( 'First call to setTP')
        self.setTP(self.T, self.P)
        

        # Note that Psat is NOT supported
            
    def check_visc_support(self):
        try:
            _ =  PropsSI('V','T',self.Tmin_si,'P',self.Pmax_si, self.fluid ) # Vmax at Tmin
            _ =  PropsSI('V','T',self.Tmax_si,'P',self.Pmax_si, self.fluid ) # Vmin at Tmax
            self.visc_is_supported = True
        except:
            self.visc_is_supported = False
            if self.show_warnings:
                print( 'check_visc_support: Visc NOT supported for fluid:', self.symbol)

        # print( 'Found visc_is_supported =', self.visc_is_supported)


    def set_warnings(self, show_warnings=2):
        """Set the show_warnings flag"""
        if show_warnings is True:
            show_warnings = 1
        elif show_warnings is False:
            show_warnings = 0

        self.save_show_warnings = self.show_warnings
        self.show_warnings = show_warnings

    def pause_warnings(self):
        """Pause all warnings until self.show_warnings is changed"""
        self.save_show_warnings = self.show_warnings
        self.show_warnings = 0

    def restore_warnings(self):
        """Restore the self.show_warnings flag to the last pause_warnings condition"""
        try:
            self.show_warnings = self.save_show_warnings
        except:
            pass
    

    def getStrTransport(self):
        '''create a string from the Transport properties'''
        Visc = self.Visc  * 1.0E5
        return  "%s Cp=%6g Visc=%6g ThCond=%6g" %\
        (self.fluid,self.Cp, Visc, self.Cond)

    def getStrTPD(self):
        '''create a string from the TPDEHS properties'''
        return  "%s T=%6.1f P=%6.1f D=%.4f E=%6.2f H=%6.2f S=%.3f" %\
        (self.fluid,self.T, self.P, self.D, self.E, self.H, self.S)

    # def getStrTPDphase(self):
    #     '''create a string from the TPDEHS properties'''
    #     return  "%s T=%6.1f P=%6.1f D=%.4f E=%6.2f H=%6.2f S=%.3f" %\
    #     (self.fluid,self.T, self.P, self.D, self.E, self.H, self.S)

    def printTPD(self):
        '''print a string from the TPDEHS properties'''
        print(self.getStrTPD())

    def printTransport(self):
        '''print a string of Transport properties'''
        print(self.getStrTransport())
        
    def get_D_at_T(self, T): # degR
        Psi = PSI_fromEng( 5000 )
        
        Tsi = TSI_fromEng( T )
        
        D = Deng_fromSI( PropsSI('D', 'T',Tsi,'P',Psi, self.fluid) )
        return D # lbm/cuft


    def setTP(self,T=530.0,P=1000.0):
        '''Calc props from T and P'''


        if T < self.Tmin:
            if self.auto_fix_value_errors:
                if self.show_warnings and self.Tmin - T > 0.01:
                    print( 'T too low in setTP. Changed T=%g to Tmin=%g'%( T, self.Tmin ))
                T = self.Tmin
            else:
                print_avoid_valerr()
                raise ValueError("Input T too low in setTP. T=%g,  Tmin=%g"%(T, self.Tmin))
        if T > self.Tmax:
            if self.auto_fix_value_errors:
                # if self.show_warnings and T - self.Tmax  > 0.01:
                if self.show_warnings and T > self.Tmax + 0.1 :
                    print( 'T too high in setTP. Changed T=%g to Tmax=%g'%( T, self.Tmax ))
                T = self.Tmax
            else:
                print_avoid_valerr()
                raise ValueError("Input T too high in setTP. T=%g,  Tmin=%g"%(T, self.Tmax))

        self.Pinput = P # save input P in case Psat changes it.

        if P > self.Pmax: # Not a bad infraction so just warn
            if self.show_warnings:
                print( 'P=%g is greater than Pmax=%g'%(P, self.Pmax))


        self.T = T
        self.P = P
        
        Tsi = TSI_fromEng( T )
        Psi = PSI_fromEng( P )
        def get_prop( prop_desc='H' ):
            try:
                prop, good_Psi = safe_get_INCOMP_prop( prop_desc, Psi_val=Psi, ind_name='T', ind_si_val=Tsi, 
                                                       symbol=self.symbol, show_warnings=self.show_warnings>1,
                                                       Pmax=self.Pmax )
                return prop
            except:
                return float('inf')

        # self.D = Deng_fromSI( get_prop('D') )
        self.D = self.get_D_at_T( self.T )

        self.rho = self.D / 1728.0

        self.E = UHeng_fromSI( get_prop('U') )
        self.H = UHeng_fromSI( get_prop('H') )
        self.S = Seng_fromSI( get_prop('S') )
        self.Cp = CPeng_fromSI( get_prop('C') )

        # Some INCOMP fluids have no viscosity data
        if self.visc_is_supported:
            self.Visc = Veng_fromSI( get_prop('V') )# multiply only for display * 1.0E5
        else:
            self.Visc = float('inf')

        self.Cond = CondEng_fromSI( get_prop('L') )
        
    def setPD(self,P=1000.0,D=0.01):
        '''
        Calc props from P and D
        NOTE: The pressure has NO EFFECT on calculated temperature for incompressible density.
        '''

        if D < self.Dmin:
            if self.auto_fix_value_errors:
                if self.show_warnings:
                    print( 'D too low in setDP. Changed D=%g to Dmin=%g'%( D, self.Dmin ))
                D = self.Dmin
            else:
                print_avoid_valerr()
                raise ValueError("Input D too low in setPD. D=%g,  Dmin=%g"%(D, self.Dmin))

        if D > self.Dmax:
            if self.auto_fix_value_errors:
                if self.show_warnings:
                    print( 'D too high in setDP. Changed D=%g to Dmax=%g'%( D, self.Dmax ))
                D = self.Dmax
            else:
                print_avoid_valerr()
                raise ValueError("Input D too high in setPD. D=%g,  Dmax=%g"%(D, self.Dmax))


        self.Pinput = P # save input P in case Psat changes it.
        # P = self.adjust_P_for_Psat( P, ind_param="D", ind_val=D )

        if P > self.Pinput:
            if self.show_warnings:
                print( 'P too low in setPD. Changed P=%g to Psat=%g'%( self.Pinput, P ))

        self.P = P 
        self.D = D

        self.T, err_flag = find_T_from_Dterp(self, D)
        if err_flag and self.show_warnings:
            print( 'WARNING: in setPD, find_T_from_Dterp has an error.')

        self.setTP( self.T, self.P )
        


    def printProps(self):
        '''print a multiline property summary with units'''

        # get formatted floats that are same length to help the look of table
        SGL = Same_g_len(self, ['Cond', "Cp", 'D', 'E', 'H', 'P', 'S', 'T', 'Visc', 'rho', 'percentage', 'T_freeze'] ) # , 'Tsat'


        print("State Point for fluid",self.fluid,"("+ self.symbol +")")
        print("T =%s"%SGL.T," degR,",    '                             Range(%8g - %8g) degR'%(self.Tmin, self.Tmax))

        print("P =%s"%SGL.P," psia",    '                              Range(%8g - %8g) psia'%(self.Pmin, self.Pmax))
        if self.Pinput != self.P:
            print("    Pinput =%8g"%self.Pinput, '(Adjusted due to Psat)')

        print("D =%s"%SGL.D," lbm/cuft",    '                          Range(%8g - %8g) lbm/cuft'%(self.Dmin, self.Dmax))

        print("E =%s"%SGL.E," BTU/lbm",    '                           Range(%8g - %8g) BTU/lbm'%(self.Emin, self.Emax))
        print("H =%s"%SGL.H," BTU/lbm",    '                           Range(%8g - %8g) BTU/lbm'%(self.Hmin, self.Hmax))
        print("S =%s"%SGL.S," BTU/lbm degR",    '                      Range(%.6f - %.6f) BTU/lbm degR'%(self.Smin, self.Smax))
        print("Cp=%s"%SGL.Cp," BTU/lbm degR",   '                      Range(%8g - %8g) BTU/lbm degR'%(self.Cpmin, self.Cpmax))

        if self.Visc < float('inf'):
            # print( "Visc=%s, Viscmin=%s, Viscmax=%s"%(SGL.Visc, self.Viscmin, self.Viscmax) )
            print("V =%s"%SGL.Visc," viscosity [1.0E5 * lbm/ft-sec]", 
                  '    Range(%8g - %8g)'%(self.Viscmin * 1.0E5, self.Viscmax * 1.0E5) )
        else:
            print("V =UNDEFINED","viscosity [1.0E5 * lbm/ft-sec]" )

        print("C =%s"%SGL.Cond," thermal conductivity [BTU/ft-hr-R]", 'Range(%8g - %8g)'%(self.Condmin, self.Condmax) )

        if self.T_freeze > 0:
            print( "    T_freeze = %g degR"%self.T_freeze)

        print("    rho      = %s"%SGL.rho," lbm/cuin",   '              Range(%.6f - %.6f) lbm/cuin'%(self.rho_min, self.rho_max))
        print("    mass%%    = %s"%SGL.percentage,"base mass percent", '      Range(%s - %s)'%(self.pcent_min_str, self.pcent_max_str))
        

    def printSIProps(self):
        '''print a multiline property summary with SI units'''

        ec_incomp_solnL = ['Cond','Condmax', 'Condmin', 'Cp', 'Cpmax', 'Cpmin', 'D', 'Dmax', 'Dmin', 'E', 'Emax', 'Emin', 
                           'fluid', 'H', 'Hmax', 'Hmin', 'P', 'percentage', 'pcent_max_str', 'pcent_min_str', 'Pinput', 
                           'Pmax', 'Pmin', 'rho', 'rho_max', 'rho_min', 'S', 'Smax', 'Smin', 'symbol', 'T', 'T_freeze', 
                           'Tmax', 'Tmin', 'Visc', 'Viscmax', 'Viscmin']
        si_obj = SI_obj( self, ec_incomp_solnL)


        print("State Point for fluid",self.fluid,"("+ self.symbol +")")
        print("T =%s"%si_obj.T," degK,",    '                       Range(%s - %s) degK'%(si_obj.Tmin.strip(), si_obj.Tmax.strip()))

        print("P =%s"%si_obj.P," Pa",    '                          Range(%s - %s) Pa'%(si_obj.Pmin.strip(), si_obj.Pmax.strip()))
        if self.Pinput != self.P:
            print("    Pinput =%s"%si_obj.Pinput, '(Adjusted due to Psat)')

        print("D =%s"%si_obj.D," kg/m^3",    '                      Range(%s - %s) lbm/cuft'%(si_obj.Dmin.strip(), si_obj.Dmax.strip()))

        print("E =%s"%si_obj.E," J/kg",    '                        Range(%s - %s) J/kg'%(si_obj.Emin.strip(), si_obj.Emax.strip()))
        print("H =%s"%si_obj.H," J/kg",    '                        Range(%s - %s) J/kg'%(si_obj.Hmin.strip(), si_obj.Hmax.strip()))
        print("S =%s"%si_obj.S," J/kg/K",    '                      Range(%s - %s) J/kg/K'%(si_obj.Smin.strip(), si_obj.Smax.strip()))
        print("Cp=%s"%si_obj.Cp," J/kg/K",   '                      Range(%s - %s) J/kg/K'%(si_obj.Cpmin.strip(), si_obj.Cpmax.strip()))

        if self.Visc < float('inf'):
            # print( "Visc=%s, Viscmin=%s, Viscmax=%s"%(si_obj.Visc, si_obj.Viscmin.strip(), si_obj.Viscmax.strip()) )
            print("V =%s"%si_obj.Visc," viscosity Pa s", 
                  '              Range(%s - %s) Pa s'%(si_obj.Viscmin.strip(), si_obj.Viscmax.strip()) )
        else:
            print("V =UNDEFINED","viscosity Pa s" )

        print("C =%s"%si_obj.Cond," thermal conductivity W/m/K", '  Range(%s - %s) W/m/K'%(si_obj.Condmin.strip(), si_obj.Condmax.strip()) )

        if self.T_freeze > 0:
            print( "    T_freeze = %s degK"%si_obj.T_freeze)

        print("    rho      = %s"%si_obj.rho," g/cm^3",   '          Range(%s - %s) g/cm^3'%(si_obj.rho_min.strip(), si_obj.rho_max.strip()))
        print("    mass%%    = %s"%si_obj.percentage,"base mass percent", 'Range(%s - %s)'%(si_obj.pcent_min_str.strip(), si_obj.pcent_max_str.strip()))
        
        
        
    def calc_min_max_props(self, do_print=False):
        # print( '.......................Entered calc_min_max_props ..........................')
        T_save = self.T
        P_save = self.P

        self.psat_is_supported = False

        # self.Tsat_min = self.get_Tsat( self.Pmin )
        # self.Tsat_max = self.get_Tsat( self.Pmax )


        # make a list of the four corners of the T,P space
        tpL = [(self.Tmin, 0), (self.Tmin, self.Pmax), (self.Tmax, 0), (self.Tmax, self.Pmax)]

        resultD = {} # key:min/max name (e.g. Dmin, Cpmax), value:float value

        # set starting values for min/max of these properties
        for prop_name in [ 'D', 'E', 'H', 'S', 'Cp', 'Visc', 'Cond' ]:
            resultD[prop_name+'min'] = float('inf')
            resultD[prop_name+'max'] = float('-inf')

        self.pause_warnings()
        for T,P in tpL:
            self.setTP( T, P )

            for prop_name in [ 'D', 'E', 'H', 'S', 'Cp', 'Visc', 'Cond' ]:
                if getattr(self, prop_name) < resultD[prop_name+'min']:
                    resultD[prop_name+'min'] = getattr(self, prop_name)

                if getattr(self, prop_name) > resultD[prop_name+'max']:
                    resultD[prop_name+'max'] = getattr(self, prop_name)
        self.restore_warnings()

        # set min/max values found in resultD
        for prop_name, value in resultD.items():
            # print( 'setting ', prop_name, value)
            setattr( self, prop_name, value )

        self.rho_min = self.Dmin / 1728.0
        self.rho_max = self.Dmax / 1728.0
        
        self.Dmin_si = DSI_fromEng( self.Dmin )
        self.Dmax_si =  DSI_fromEng( self.Dmax )

        if do_print:
            minmaxL = ['Cond', "Cp", 'D', 'E', 'H', 'P', 'S', 'T', 'Visc']

            for name in minmaxL:
                max_name = name + 'max'
                print( '%10s'%max_name, '%9g'%getattr(self, max_name) )
                min_name = name + 'min'
                print( '%10s'%min_name, '%9g'%getattr(self, min_name) )
        # print( '.......................Finished calc_min_max_props ..........................')

        # Restore state to T and P from calling method
        self.setTP( T_save, P_save)        



def dev_tests():
    """
    === incompressible_list_solution ===
    AEG,AKF,AL,AN,APG,ExampleDigital,ExampleMelinder,ExampleSecCool,ExampleSolution,
    FRE,GKN,IceEA,IceNA,IcePG,LiBr,MAM,MAM2,MCA,MCA2,MEA,MEA2,MEG,MEG2,MGL,MGL2,
    MITSW,MKA,MKA2,MKC,MKC2,MKF,MLI,MMA,MMA2,MMG,MMG2,MNA,MNA2,MPG,MPG2,PK2,PKL,
    VCA,VKC,VMA,VMG,VNA,ZAC,ZFC,ZLC,ZM,ZMC
    """

    symbol = 'MEG-99%'    
    print( '='*22, "%s with %% too high"%symbol, '='*22 )
    C = EC_Incomp_Soln( symbol=symbol, T=None, P=None, 
                        auto_fix_value_errors=True, show_warnings=1 )
    print()

    symbol = 'IceNA-1%'
    print( '='*22, "%s with %% too low"%symbol, '='*22 )
    C = EC_Incomp_Soln( symbol=symbol, T=None, P=None, 
                        auto_fix_value_errors=True, show_warnings=1 )
    print()

    
    symbol = 'MEG-20%'
    print( '='*22, "%s at T=None P=None"%symbol, '='*22 )
    C = EC_Incomp_Soln( symbol=symbol, T=None, P=None )

    C.printTPD()
    C.printTransport()
    print()

    C.printProps()
    
    # C.setTP(T= (C.Tmin+C.Tmax)/2.0, P=100) # this temperature throws an exception ???

    
    symbol = 'MPG2-57%'
    print( '='*22, "%s Has frac_max round-off problem"%symbol, '='*22 )
    C = EC_Incomp_Soln( symbol=symbol, T=None, P=None )

    C.printTPD()
    C.printTransport()
    print()

    C.printProps()

    print( '='*22, "%s Check setPD"%symbol, '='*22 )
    C.setPD( C.P, C.D )
    C.printProps()


    banner( 'This should fail due to "auto_fix_value_errors = False"\nand max mass fraction = "MEG-60%"' )
    print( '='*22, "print full %s properties"%symbol, '='*22 )
    try:
        C_bad = EC_Incomp_Soln( symbol='MEG-99%', T=None, P=0 )
    except:
        tb_str = traceback.format_exc()
        print( tb_str.split('raise ValueError')[-1])
        # print( tb_str )
    

    banner( 'This should fail due to "auto_fix_value_errors = False"\nand T=100 below Tmin=428.68' )
    print( '='*22, "print full %s properties"%symbol, '='*22 )
    try:
        C_bad = EC_Incomp_Soln( symbol='MEG-50%', T=100, P=0 )
        C_bad.printTPD()
    except:
        tb_str = traceback.format_exc()
        print( tb_str.split('raise ValueError')[-1])
        # print( tb_str )


    print( '='*22, "Check omission of mass fraction", '='*22 )
    C = EC_Incomp_Soln( symbol='MEG', T=None, P=None )
    C.printProps()

    print( '='*22, "Check printSIProps", '='*22 )
    C.printSIProps()
    

if __name__ == '__main__':
    dev_tests()
