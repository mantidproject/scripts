"""
Unit test for Vesuvio reduction

Assumes that mantid can be imported and the data paths
are configured to find the Vesuvio data
"""
import unittest

from mantid.api import AlgorithmManager

class VesuvioReductionTest(unittest.TestCase):

    # -------------- Success cases ------------------

    def test_single_run_produces_correct_output_workspace_index0(self):
        profiles = "function=GramCharlier,width=[2, 5, 7],hermite_coeffs=[1, 0, 1],k_free=1,sears_flag=1;"\
                   "function=Gaussian,width=10;function=Gaussian,width=13;function=Gaussian,width=30;"

        alg = self._create_algorithm(Runs="15039", IPFilename="IP0004_10.par",
                                     Masses=[1.0079, 33, 27, 133],
                                     MassProfiles=profiles)
        alg.execute()
        output_ws = alg.getProperty("FittedWorkspace").value

        self.assertEqual(7, output_ws.getNumberHistograms())
        self.assertAlmostEqual(50.0, output_ws.readX(0)[0])
        self.assertAlmostEqual(562.0, output_ws.readX(0)[-1])

        self.assertAlmostEqual(0.012819595690617192, output_ws.readY(0)[0])
        self.assertAlmostEqual(0.020656898709105809, output_ws.readY(0)[-1])


    # -------------- Failure cases ------------------

    def test_empty_runs_raises_error(self):
        alg = self._create_algorithm()

        self.assertRaises(ValueError, alg.setProperty,
                          "Runs", "")

    def test_empty_masses_raises_error(self):
        alg = self._create_algorithm()

        self.assertRaises(ValueError, alg.setProperty,
                          "Masses", [])

    def test_empty_functions_raises_error(self):
        alg = self._create_algorithm()

        self.assertRaises(ValueError, alg.setProperty, "MassProfiles", "")

    def test_number_functions_in_list_not_matching_length_masses_throws_error(self):
        alg = self._create_algorithm(Runs="15039-15045", IPFilename="IP0004_10.par",
                                     Masses=[1.0079, 33],
                                     MassProfiles="function=Gaussian,width=5;")

        self.assertRaises(RuntimeError, alg.execute)

    # -------------- Helpers --------------------

    def _create_algorithm(self, **kwargs):
        alg = AlgorithmManager.createUnmanaged("VesuvioReduction")
        alg.initialize()
        alg.setChild(True)
        alg.setProperty("FittedWorkspace", "__unused")
        alg.setProperty("FittedParameters", "__unused")
        for key, value in kwargs.iteritems():
            alg.setProperty(key, value)
        return alg

if __name__ == "__main__":
    unittest.main()

