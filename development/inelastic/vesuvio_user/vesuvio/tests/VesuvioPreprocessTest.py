"""
Unit test for Vesuvio preprocessing

Assumes that mantid can be imported and the data paths
are configured to find the Vesuvio data
"""
import unittest

from mantid.api import AlgorithmManager
from mantid.simpleapi import LoadVesuvio
import vesuvio

class VesuvioPreprocessTest(unittest.TestCase):

    _test_ws = None

    def setUp(self):
        if self._test_ws is not None:
            return
        # Cache a TOF workspace
        self.__class__._test_ws = \
            vesuvio.workflow.load_and_crop_data(runs="15039-15045",
                                                ip_file="IP0004_10.par",
                                                spectra="135-136")

    # -------------- Success cases ------------------

    def test_smooth_uses_requested_number_of_points(self):
        alg = self._create_algorithm(InputWorkspace=self._test_ws,
                                     Smoothing="Neighbour", SmoothingOptions="NPoints=3")
        alg.execute()
        output_ws = alg.getProperty("OutputWorkspace").value

        self.assertEqual(2, output_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, output_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, output_ws.readX(0)[-1])

        self.assertAlmostEqual(0.0071192100574770656, output_ws.readY(0)[0])
        self.assertAlmostEqual(0.011063685963343062, output_ws.readY(0)[-1])
        self.assertAlmostEqual(0.002081985833021438, output_ws.readY(1)[0])
        self.assertAlmostEqual(-0.01209794157313937, output_ws.readY(1)[-1])

    def xtest_single_run_produces_correct_output_workspace_index1_kfixed_no_background(self):
        profiles = "function=GramCharlier,width=[2, 5, 7],hermite_coeffs=[1, 0, 0],k_free=0,sears_flag=1;"\
                   "function=Gaussian,width=10;function=Gaussian,width=13;function=Gaussian,width=30;"

        alg = self._create_algorithm(InputWorkspace=self._test_ws, WorkspaceIndex=1,
                                     Masses=[1.0079, 16, 27, 133],
                                     MassProfiles=profiles,
                                     IntensityConstraints="[0,1,0,-4]")
        alg.execute()
        output_ws = alg.getProperty("OutputWorkspace").value

        self.assertEqual(7, output_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, output_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, output_ws.readX(0)[-1])

        self.assertAlmostEqual(-0.005852648610523481, output_ws.readY(0)[0])
        self.assertAlmostEqual(-0.013112461599666836, output_ws.readY(0)[-1])
        self.assertAlmostEqual(1.5165126628163728e-05, output_ws.readY(1)[0])
        self.assertAlmostEqual(7.6619346342727019e-05, output_ws.readY(1)[-1])


    def xtest_single_run_produces_correct_output_workspace_index0_kfixed_including_background(self):
        profiles = "function=GramCharlier,width=[2, 5, 7],hermite_coeffs=[1, 0, 0],k_free=0,sears_flag=1;"\
                   "function=Gaussian,width=10;function=Gaussian,width=13;function=Gaussian,width=30;"
        background = "function=Polynomial,order=3"

        alg = self._create_algorithm(InputWorkspace=self._test_ws, WorkspaceIndex=0,
                                     Masses=[1.0079, 16.0, 27.0, 133.0],
                                     MassProfiles=profiles,
                                     Background=background,
                                     IntensityConstraints="[0,1,0,-4]")

        alg.execute()
        output_ws = alg.getProperty("OutputWorkspace").value

        self.assertEqual(8, output_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, output_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, output_ws.readX(0)[-1])

        self.assertAlmostEqual(0.000928695463881635, output_ws.readY(0)[0])
        self.assertAlmostEqual(0.00722948549525415, output_ws.readY(0)[-1])
        self.assertAlmostEqual(-0.00756178413274695, output_ws.readY(1)[0])
        self.assertAlmostEqual(0.00355843687365601, output_ws.readY(1)[-1])

    # -------------- Failure cases ------------------

    def test_invalid_smooth_opt_raises_error_on_validate(self):
        alg = self._create_algorithm(InputWorkspace=self._test_ws,
                                     Smoothing="Neighbour", SmoothingOptions="npts=3")
        self.assertRaises(RuntimeError, alg.execute)


    # -------------- Helpers --------------------

    def _create_algorithm(self, **kwargs):
        alg = AlgorithmManager.createUnmanaged("VesuvioPreprocess")
        alg.initialize()
        alg.setChild(True)
        alg.setProperty("OutputWorkspace", "__unused")
        for key, value in kwargs.iteritems():
            alg.setProperty(key, value)
        return alg

if __name__ == "__main__":
    unittest.main()

