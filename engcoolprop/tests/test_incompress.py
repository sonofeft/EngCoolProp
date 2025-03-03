
import unittest
# import unittest2 as unittest # for versions of python < 2.7

"""
        Method                            Checks that
self.assertEqual(a, b)                      a == b   
self.assertNotEqual(a, b)                   a != b   
self.assertTrue(x)                          bool(x) is True  
self.assertFalse(x)                         bool(x) is False     
self.assertIs(a, b)                         a is b
self.assertIsNot(a, b)                      a is not b
self.assertIsNone(x)                        x is None 
self.assertIsNotNone(x)                     x is not None 
self.assertIn(a, b)                         a in b
self.assertNotIn(a, b)                      a not in b
self.assertIsInstance(a, b)                 isinstance(a, b)  
self.assertNotIsInstance(a, b)              not isinstance(a, b)  
self.assertAlmostEqual(a, b, places=5)      a within 5 decimal places of b
self.assertNotAlmostEqual(a, b, delta=0.1)  a is not within 0.1 of b
self.assertGreater(a, b)                    a is > b
self.assertGreaterEqual(a, b)               a is >= b
self.assertLess(a, b)                       a is < b
self.assertLessEqual(a, b)                  a is <= b

for expected exceptions, use:

with self.assertRaises(Exception):
    blah...blah...blah

with self.assertRaises(KeyError):
    blah...blah...blah

Test if __name__ == "__main__":


    def test__main__(self):
        # Change bottom of source file to call "dev_tests"
        
         def dev_tests():
            pass

         if __name__ == "__main__":
            dev_tests()
            
        # then test by calling <name>.dev_tests()


See:
      https://docs.python.org/2/library/unittest.html
         or
      https://docs.python.org/dev/library/unittest.html
for more assert options
"""

import sys, os

from engcoolprop.ec_incomp_fluid import EC_Incomp_Fluid, dev_tests
from engcoolprop.utils import incomp_pure_fluidL

here = os.path.abspath(os.path.dirname(__file__)) # Needed for py.test
up_one = os.path.split( here )[0]  # Needed to find models development version
if here not in sys.path[:2]:
    sys.path.insert(0, here)
if up_one not in sys.path[:2]:
    sys.path.insert(0, up_one)

class MyTest(unittest.TestCase):


    def test_should_always_pass_cleanly(self):
        """Should always pass cleanly."""
        pass

    def test_myclass_existence(self):
        """Check that myclass exists"""

        ec_inc = EC_Incomp_Fluid( symbol='DowJ', auto_fix_value_errors=True )
        
        # See if the self.EC_Incomp_Fluid object exists
        self.assertIsInstance(ec_inc, EC_Incomp_Fluid, msg=None)

    def test_make_all_fluids(self):
        """test make all fluids"""

        for symbol in incomp_pure_fluidL:
            ec_inc = EC_Incomp_Fluid( symbol=symbol, auto_fix_value_errors=True )
            self.assertIsInstance(ec_inc, EC_Incomp_Fluid, msg=None)

    def test_setPropsDP_all_fluids(self):
        """test setPropsDP for all fluids"""

        for symbol in incomp_pure_fluidL:
            ec_inc = EC_Incomp_Fluid( symbol=symbol, auto_fix_value_errors=True )
            T = ec_inc.T
            ec_inc.setProps( D=ec_inc.D, P=ec_inc.P)

            self.assertAlmostEqual(T, ec_inc.T, delta=0.1)

    def test_setPropsHP_all_fluids(self):
        """test setPropsHP for all fluids"""

        for symbol in incomp_pure_fluidL:
            if symbol == 'Air':
                continue
            ec_inc = EC_Incomp_Fluid( symbol=symbol, auto_fix_value_errors=True )
            T = ec_inc.T
            ec_inc.setProps( H=ec_inc.H, P=ec_inc.P)

            self.assertAlmostEqual(T, ec_inc.T, delta=0.1)
            

    def test_setPropsSP_all_fluids(self):
        """test setPropsSP for all fluids"""

        for symbol in incomp_pure_fluidL:
            if symbol == 'Air':
                continue
            ec_inc = EC_Incomp_Fluid( symbol=symbol, auto_fix_value_errors=True )
            T = ec_inc.T
            ec_inc.setProps( S=ec_inc.S, P=ec_inc.P)

            self.assertAlmostEqual(T, ec_inc.T, delta=0.1)
            

    def test_setPropsTP_Tmin_Tmax_all_fluids(self):
        """test setPropsDP for all fluids"""

        for symbol in incomp_pure_fluidL:
            ec_inc = EC_Incomp_Fluid( symbol=symbol, auto_fix_value_errors=True )
            
            # check that safeguards for Tmin/Tmax are working
            ec_inc.setProps( T=ec_inc.Tmin-1, P=ec_inc.P)
            ec_inc.setProps( T=ec_inc.Tmax+1, P=ec_inc.P)

            

    def test__main__(self):
        old_sys_argv = list(sys.argv)
        sys.argv = list(sys.argv)
        sys.argv.append('suppress_show')
        
        try:
            dev_tests()
        except:
            raise Exception('ERROR... failed in __main__ routine')
        finally:
            sys.argv = old_sys_argv
        

if __name__ == '__main__':
    # Can test just this file from command prompt
    #  or it can be part of test discovery from nose, unittest, pytest, etc.
    unittest.main()

