"""
Unit test for Vesuvio reduction

Assumes that mantid can be imported and the data paths
are configured to find the Vesuvio data
"""
import unittest

from mantid.api import AlgorithmManager

class VesuvioReductionTest(unittest.TestCase):

    # -------------- Success cases ------------------

    # -------------- Failure cases ------------------

    def test_runs_greater_two_length_raise_error(self):
        alg = self._create_initialized_algorithm()

        self.assertRaises(ValueError, alg.setProperty,
                          "Runs", [14188,14189,14190])

    def test_empty_masses_raises_error(self):
        alg = self._create_initialized_algorithm()

        self.assertRaises(ValueError, alg.setProperty,
                          "Masses", [])

    def test_empty_functions_raises_error(self):
        alg = self._create_initialized_algorithm()

        self.assertRaises(ValueError, alg.setProperty, "MassProfiles", [])

    def test_functions_list_not_matching_length_masses_throws_error(self):
        alg = self._create_initialized_algorithm()
        alg.setProperty("Runs", [14188, 14189])
        alg.setProperty("Masses", [1.0079, 33])
        alg.setProperty("Widths", [5, 30])
        alg.setProperty("MassProfiles", ["GramCharlier"])

        self.assertRaises(RuntimeError, alg.execute)

    # -------------- Helpers --------------------

    def _create_initialized_algorithm(self):
        alg = AlgorithmManager.createUnmanaged("VesuvioReduction")
        alg.initialize()
        return alg

if __name__ == "__main__":
    unittest.main()

