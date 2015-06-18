"""
Unit test for Vesuvio corrections steps

Assumes that mantid can be imported and the data paths
are configured to find the Vesuvio data
"""
import unittest

from mantid.api import AlgorithmManager
from mantid.simpleapi import LoadVesuvio
import vesuvio

#TODO: There is a lot of missing tests here that are to be added once the
#      correction scripts are accepted by the scientists to be doing the correct
#      thing

class VesuvioCorrections(unittest.TestCase):

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

    def xtest_smooth_uses_requested_number_of_points(self):
        alg = self._create_algorithm(InputWorkspace=self._test_ws,
                                     Smoothing="Neighbour", SmoothingOptions="NPoints=3",
                                     BadDataError=-1)
        alg.execute()
        output_ws = alg.getProperty("OutputWorkspace").value

        self.assertEqual(2, output_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, output_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, output_ws.readX(0)[-1])

        self.assertAlmostEqual(0.0071192100574770656, output_ws.readY(0)[0])
        self.assertAlmostEqual(0.011063685963343062, output_ws.readY(0)[-1])
        self.assertAlmostEqual(0.002081985833021438, output_ws.readY(1)[0])
        self.assertAlmostEqual(-0.01209794157313937, output_ws.readY(1)[-1])

    def xtest_mask_only_masks_over_threshold(self):
        err_start = self._test_ws.readE(1)[-1]
        self._test_ws.dataE(1)[-1] = 1.5e6

        alg = self._create_algorithm(InputWorkspace=self._test_ws,
                                     Smoothing="None", BadDataError=1.0e6)
        alg.execute()
        self._test_ws.dataE(1)[-1] = err_start
        output_ws = alg.getProperty("OutputWorkspace").value

        self.assertEqual(2, output_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, output_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, output_ws.readX(0)[-1])

        self.assertAlmostEqual(0.00092869546388163471, output_ws.readY(0)[0])
        self.assertAlmostEqual(0.007229485495254151, output_ws.readY(0)[-1])
        self.assertAlmostEqual(-0.005852648610523481, output_ws.readY(1)[0])
        # Masked
        self.assertAlmostEqual(0.0, output_ws.readY(1)[-1])

    # -------------- Failure cases ------------------

    def test_gamma_correction_without_params_raises_error(self):
        alg = self._create_algorithm(InputWorkspace=self._test_ws,
                                     GammaBackground=True)
        self.assertRaises(RuntimeError, alg.execute)

    def test_ms_correction_without_params_raises_error(self):
        alg = self._create_algorithm(InputWorkspace=self._test_ws,
                                     MultipleScattering=True)
        self.assertRaises(RuntimeError, alg.execute)

    # -------------- Helpers --------------------

    def _create_algorithm(self, **kwargs):
        alg = AlgorithmManager.createUnmanaged("VesuvioCorrections")
        alg.initialize()
        alg.setChild(True)
        alg.setProperty("OutputWorkspace", "__unused")
        for key, value in kwargs.iteritems():
            alg.setProperty(key, value)
        return alg

if __name__ == "__main__":
    unittest.main()

