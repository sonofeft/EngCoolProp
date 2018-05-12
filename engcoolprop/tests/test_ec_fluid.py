import unittest
# import unittest2 as unittest # for versions of python < 2.7

"""
        Method                  Checks that
self.assertEqual(a, b)           a == b   
self.assertNotEqual(a, b)        a != b   
self.assertTrue(x)               bool(x) is True  
self.assertFalse(x)              bool(x) is False     
self.assertIs(a, b)              a is b
self.assertIsNot(a, b)           a is not b
self.assertIsNone(x)             x is None 
self.assertIsNotNone(x)          x is not None 
self.assertIn(a, b)              a in b
self.assertNotIn(a, b)           a not in b
self.assertIsInstance(a, b)      isinstance(a, b)  
self.assertNotIsInstance(a, b)   not isinstance(a, b)  

See:
      https://docs.python.org/2/library/unittest.html
         or
      https://docs.python.org/dev/library/unittest.html
for more assert options
"""

import sys, os

here = os.path.abspath(os.path.dirname(__file__)) # Needed for py.test
up_one = os.path.split( here )[0]  # Needed to find engcoolprop development version
if here not in sys.path[:2]:
    sys.path.insert(0, here)
if up_one not in sys.path[:2]:
    sys.path.insert(0, up_one)

from engcoolprop.ec_fluid import EC_Fluid

class MyTest(unittest.TestCase):


    def test_should_always_pass_cleanly(self):
        """Should always pass cleanly."""
        self.assertAlmostEqual(0.0, 0.0000005, places=5)
        pass

    def test_myclass_existence(self):
        """Check that myclass exists"""
        ec = EC_Fluid( symbol="N2" )

        # See if the self.myclass object exists
        self.assertTrue( ec )
        del( ec )

    def test_setProps_N2(self):
        '''Check call to setProps with fluid N2'''
        ec = EC_Fluid( symbol="N2", T=200., P=300. )
        ec.setProps( T=530.0, P=1000.0)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_setProps_O2(self):
        '''Check call to setProps with fluid O2'''
        ec = EC_Fluid( symbol="O2" , T=200., P=300.)
        ec.setProps( T=530.0, P=1000.0)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_dH_FromHref_N2(self):
        '''Check call to dH_FromHref with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        dH = ec.dH_FromHref()
        
        self.assertAlmostEqual(dH, 177.5413109961724 , places=5)

        del( ec )


    def test_dH_FromHref_O2(self):
        '''Check call to dH_FromHref with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        dH = ec.dH_FromHref()
        
        self.assertAlmostEqual(dH, 165.11167057591305, places=5)

        del( ec )


    def test_setPropsFromAS_N2(self):
        '''Check call to setPropsFromAS with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.setPropsFromAS()
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_setPropsFromAS_O2(self):
        '''Check call to setPropsFromAS with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.setPropsFromAS()
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_setTP_N2(self):
        '''Check call to setTP with fluid N2'''
        ec = EC_Fluid( symbol="N2", T=200., P=300. )
        ec.setTP(T=530.0,P=1000.0)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_setTP_O2(self):
        '''Check call to setTP with fluid O2'''
        ec = EC_Fluid( symbol="O2", T=200., P=300. )
        ec.setTP(T=530.0,P=1000.0)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_newDE_N2(self):
        '''Check call to newDE with fluid N2'''
        ec = EC_Fluid( symbol="N2" , T=200., P=300.)
        
        ec.newDE(D=4.94659381374,E=87.6731930702)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_newDE_O2(self):
        '''Check call to newDE with fluid O2'''
        ec = EC_Fluid( symbol="O2" , T=200., P=300.)
        ec.newDE(D=5.86052712206,E=76.1981078385)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_constP_newH_N2(self):
        '''Check call to constP_newH with fluid N2'''
        ec = EC_Fluid( symbol="N2" , T=200., P=1000.)
        ec.constP_newH(125.082832491)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_constP_newH_O2(self):
        '''Check call to constP_newH with fluid O2'''
        ec = EC_Fluid( symbol="O2" , T=200., P=1000.)
        ec.constP_newH(107.773815464 )
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_constH_newP_N2(self):
        '''Check call to constH_newP with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.setProps(P=500., H=125.082832491)
        ec.constH_newP(P=1000.0)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_constH_newP_O2(self):
        '''Check call to constH_newP with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.setProps(P=500., H=107.773815464)
        ec.constH_newP(P=1000.0)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_setPH_N2(self):
        '''Check call to setPH with fluid N2'''
        ec = EC_Fluid( symbol="N2", T=200., P=1000. )
        ec.setPH(1000.0,125.082832491)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_setPH_O2(self):
        '''Check call to setPH with fluid O2'''
        ec = EC_Fluid( symbol="O2" , T=200., P=1000.)
        ec.setPH(1000.0,107.773815464)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_constS_newP_N2(self):
        '''Check call to constS_newP with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.setProps(T=500., S=1.31927325217)
        ec.constS_newP(P=1000.0)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_constS_newP_O2(self):
        '''Check call to constS_newP with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.setProps(T=500., S=1.25461825929)
        ec.constS_newP(P=1000.0)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_setTD_N2(self):
        '''Check call to setTD with fluid N2'''
        ec = EC_Fluid( symbol="N2" , T=200., P=300.)
        ec.setTD(T=530.0,D=4.94659381374)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_setTD_O2(self):
        '''Check call to setTD with fluid O2'''
        ec = EC_Fluid( symbol="O2" , T=200., P=300. )
        ec.setTD(T=530.0,D=5.86052712206)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_setPD_N2(self):
        '''Check call to setPD with fluid N2'''
        ec = EC_Fluid( symbol="N2"  , T=200., P=300.)
        ec.setPD(P=1000.0,D=4.94659381374)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_setPD_O2(self):
        '''Check call to setPD with fluid O2'''
        ec = EC_Fluid( symbol="O2"  , T=200., P=300.)
        ec.setPD(P=1000.0,D=5.86052712206)
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_getSatP_N2(self):
        '''Check call to getSatP with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        Psat = ec.getSatP( T=None)
        Psat = ec.getSatP( T=200.0)
        
        self.assertAlmostEqual(Psat, 226.61149712125732, places=5)

        del( ec )


    def test_getSatP_O2(self):
        '''Check call to getSatP with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        Psat = ec.getSatP( T=None)
        Psat = ec.getSatP( T=200.0)
        
        self.assertAlmostEqual(Psat, 85.02309450018438, places=5)


        del( ec )


    def test_getSatPandDens_N2(self):
        '''Check call to getSatPandDens with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        Psat, Dliq, Dgas = ec.getSatPandDens( T=None)
        Psat, Dliq, Dgas = ec.getSatPandDens( T=200.0)
        
        self.assertAlmostEqual(Psat, 226.61149712125732, places=5)
        self.assertAlmostEqual(Dliq, 38.2470012845856, places=5)
        self.assertAlmostEqual(Dgas, 4.946593813736714, places=5)

        del( ec )


    def test_getSatPandDens_O2(self):
        '''Check call to getSatPandDens with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        Psat, Dliq, Dgas = ec.getSatPandDens( T=None)
        Psat, Dliq, Dgas = ec.getSatPandDens( T=200.0)
        
        self.assertAlmostEqual(Psat, 85.02309450018438, places=5)
        self.assertAlmostEqual(Dliq, 64.23715727948591, places=5)
        self.assertAlmostEqual(Dgas, 5.860527122060697, places=5)

        del( ec )


    def test_getSatT_N2(self):
        '''Check call to getSatT with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        Tsat = ec.getSatT( P=None)
        Tsat = ec.getSatT( P=80.0)
        
        self.assertAlmostEqual(Tsat, 171.4764377703347, places=5)

        del( ec )


    def test_getSatT_O2(self):
        '''Check call to getSatT with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.getSatT( P=None)
        Tsat = ec.getSatT( P=None)
        Tsat = ec.getSatT( P=80.0)
        
        self.assertAlmostEqual(Tsat, 198.39085041700642, places=5)

        del( ec )


    def test_getSatTandDens_N2(self):
        '''Check call to getSatTandDens with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        Tsat, Dliq, Dgas = ec.getSatTandDens(P=None)
        Tsat, Dliq, Dgas = ec.getSatTandDens(P=80.0)
        
        self.assertAlmostEqual(Tsat, 171.4764377703347, places=5)
        self.assertAlmostEqual(Dliq, 44.747923887617794, places=5)
        self.assertAlmostEqual(Dgas, 1.4181440790822393, places=5)

        del( ec )


    def test_getSatTandDens_O2(self):
        '''Check call to getSatTandDens with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.getSatTandDens(P=None)
        Tsat, Dliq, Dgas = ec.getSatTandDens(P=80.0)
        
        self.assertAlmostEqual(Tsat, 198.39085041700642, places=5)
        self.assertAlmostEqual(Dliq, 64.56368984711504, places=5)
        self.assertAlmostEqual(Dgas, 1.3475849732116867, places=5)

        del( ec )


    def test_dHvap_N2(self):
        '''Check call to dHvap with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        dH = ec.dHvap( P=None)
        dH = ec.dHvap( P=15.0)
        
        self.assertAlmostEqual(dH, 85.5336483246067, places=5)

        del( ec )


    def test_dHvap_O2(self):
        '''Check call to dHvap with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        dH = ec.dHvap( P=None)
        dH = ec.dHvap( P=15.0)
        
        self.assertAlmostEqual(dH, 91.51510849550952, places=5)

        del( ec )


    def test_surfTen_N2(self):
        '''Check call to surfTen with fluid N2'''
        ec = EC_Fluid( symbol="N2", P=80, T=200 )
        st = ec.surfTen()        
        self.assertAlmostEqual(st, 1.1726812809126956e-05, places=5)

        del( ec )


    def test_surfTen_O2(self):
        '''Check call to surfTen with fluid O2'''
        ec = EC_Fluid( symbol="O2", P=80, T=200 )
        st = ec.surfTen()
        self.assertAlmostEqual(st, 4.638551601014212e-05, places=5)

        del( ec )


    def test_restoreFromDup_N2(self):
        '''Check call to restoreFromDup with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.restoreFromDup()
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_restoreFromDup_O2(self):
        '''Check call to restoreFromDup with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.restoreFromDup()
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=3)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_saveToDup_N2(self):
        '''Check call to saveToDup with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.saveToDup()
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_saveToDup_O2(self):
        '''Check call to saveToDup with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.saveToDup()
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_initFromObj_N2(self):
        '''Check call to initFromObj with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.initFromObj( ec )
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 4.94659381374, places=5)
        self.assertAlmostEqual(ec.S, 1.31927325217, places=5)

        del( ec )


    def test_initFromObj_O2(self):
        '''Check call to initFromObj with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.initFromObj( ec )
        
        self.assertAlmostEqual(ec.T, 530.0, places=5)
        self.assertAlmostEqual(ec.P, 1000.0, places=5)
        self.assertAlmostEqual(ec.D, 5.86052712206, places=5)
        self.assertAlmostEqual(ec.S, 1.25461825929, places=5)

        del( ec )


    def test_gamma_N2(self):
        '''Check call to gamma with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        gam = ec.gamma()
        
        self.assertAlmostEqual(gam, 1.5214910882792643, places=5)

        del( ec )


    def test_gamma_O2(self):
        '''Check call to gamma with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        gam = ec.gamma()
        
        self.assertAlmostEqual(gam, 1.539276126838285, places=5)

        del( ec )


    def test_getStrTransport_N2(self):
        '''Check call to getStrTransport with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.getStrTransport()
        
        del( ec )


    def test_getStrTransport_O2(self):
        '''Check call to getStrTransport with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.getStrTransport()
        
        del( ec )


    def test_getStrTPD_N2(self):
        '''Check call to getStrTPD with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.getStrTPD()
        
        del( ec )


    def test_getStrTPD_O2(self):
        '''Check call to getStrTPD with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.getStrTPD()
        
        del( ec )


    def test_getStrTPDphase_N2(self):
        '''Check call to getStrTPDphase with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.getStrTPDphase()

        del( ec )


    def test_getStrTPDphase_O2(self):
        '''Check call to getStrTPDphase with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.getStrTPDphase()

        del( ec )


    def test_printTPD_N2(self):
        '''Check call to printTPD with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.printTPD()

        del( ec )


    def test_printTPD_O2(self):
        '''Check call to printTPD with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.printTPD()

        del( ec )


    def test_Qdescription_N2(self):
        '''Check call to Qdescription with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        s = ec.Qdescription()
        
        self.assertEqual(s, 'Supercritical Fluid')

        del( ec )


    def test_Qdescription_O2(self):
        '''Check call to Qdescription with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        s = ec.Qdescription()
        
        self.assertEqual(s, 'Supercritical Fluid')

        del( ec )


    def test_html_desc_N2(self):
        '''Check call to html_desc with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.html_desc()

        del( ec )


    def test_html_desc_O2(self):
        '''Check call to html_desc with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.html_desc()

        del( ec )


    def test_printCriticalProps_N2(self):
        '''Check call to printCriticalProps with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.printCriticalProps()

        del( ec )


    def test_printCriticalProps_O2(self):
        '''Check call to printCriticalProps with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.printCriticalProps()

        del( ec )


    def test_printProps_N2(self):
        '''Check call to printProps with fluid N2'''
        ec = EC_Fluid( symbol="N2" )
        ec.printProps()

        del( ec )


    def test_printProps_O2(self):
        '''Check call to printProps with fluid O2'''
        ec = EC_Fluid( symbol="O2" )
        ec.printProps()

        del( ec )


if __name__ == '__main__':
    # Can test just this file from command prompt
    #  or it can be part of test discovery from nose, unittest, pytest, etc.
    unittest.main()

