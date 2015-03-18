"""
Unit test for Vesuvio reduction

Assumes that mantid can be imported and the data paths
are configured to find the Vesuvio data
"""
import unittest

from mantid.api import AlgorithmManager

class VesuvioReductionTest(unittest.TestCase):

    # -------------- Success cases ------------------

    def test_single_run_produces_correct_output(self):
        alg = self._create_algorithm(Runs="15039",
                                     Masses=[1.0079, 33, 27, 133],
                                     FixedWidths=[0, 10, 13, 30],
                                     WidthConstraints=[2, 5, 7],
                                     MassProfiles=["GramCharlier",
                                                   "Gaussian", "Gaussian", "Gaussian"])
        alg.execute()
        output_ws = alg.getProperty("OutputWorkspace").value

        self.assertEqual(64, output_ws.getNumberHistograms())

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

        self.assertRaises(ValueError, alg.setProperty, "MassProfiles", [])

    def test_functions_list_not_matching_length_masses_throws_error(self):
        alg = self._create_algorithm(Runs="15039-15045",
                                     Masses=[1.0079, 33],
                                     FixedWidths=[5,30],
                                     MassProfiles=["GramCharlier"])

        self.assertRaises(RuntimeError, alg.execute)

    def test_fixedwidth_list_not_matching_length_masses_throws_error(self):
        alg = self._create_algorithm(Runs="15039-15045",
                                     Masses=[1.0079, 33],
                                     FixedWidths=[5],
                                     MassProfiles=["GramCharlier", "Gaussian"])

        self.assertRaises(RuntimeError, alg.execute)

    def test_widthconstraints_not_3_times_length_of_number_of_non_fixed_widths_throws_error(self):
        alg = self._create_algorithm(Runs="15039-15045",
                                     Masses=[1.0079, 33],
                                     FixedWidths=[0, 5],
                                     WidthConstraints=[2, 5],
                                     MassProfiles=["GramCharlier", "Gaussian"])

        self.assertRaises(RuntimeError, alg.execute)

    # -------------- Helpers --------------------

    def _create_algorithm(self, **kwargs):
        alg = AlgorithmManager.createUnmanaged("VesuvioReduction")
        alg.initialize()
        alg.setChild(True)
        alg.setProperty("Outputworkspace", "__unused")
        for key, value in kwargs.iteritems():
            alg.setProperty(key, value)
        return alg

if __name__ == "__main__":
    unittest.main()

