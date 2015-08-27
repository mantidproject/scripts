#pylint: disable=too-many-public-methods

import unittest

from mantid.api import ITableWorkspace
from mantid.simpleapi import *


class VesuvioLt0ResolutionTest(unittest.TestCase):

    def tearDown(self):
        mtd.clear()


    def _validate_result_shape(self, forward, back, mean):
        """
        Validates the shape of the result tables.
        """

        self.assertTrue(isinstance(forward, ITableWorkspace))
        self.assertEqual(forward.columnCount(), 5)
        self.assertEqual(forward.rowCount(), 64)

        self.assertTrue(isinstance(back, ITableWorkspace))
        self.assertEqual(back.columnCount(), 5)
        self.assertEqual(back.rowCount(), 132)

        self.assertTrue(isinstance(mean, ITableWorkspace))
        self.assertEqual(mean.rowCount(), 4)


    def test_run_with_par_file(self):
        """
        Tests running using distances from an instrument parameter file.
        """
        forward, back, mean = VesuvioLt0Resolution(InstrumentParFile='IP0005.par')
        self._validate_result_shape(forward, back, mean)

        # Validate mean values
        mean_values = mean.column('Mean')
        TOLERANCE = 7
        self.assertAlmostEqual(mean_values[0], 0.0298101, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[1], 0.2581127, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[2], 0.0206698, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[3], 0.2972636, places=TOLERANCE)


if __name__ == '__main__':
    unittest.main()
