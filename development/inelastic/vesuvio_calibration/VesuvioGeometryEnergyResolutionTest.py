#pylint: disable=too-many-public-methods

import unittest

from mantid.api import (ITableWorkspace, WorkspaceGroup)
from mantid.simpleapi import *


class VesuvioGeometryEnergyResolutionTest(unittest.TestCase):

    def tearDown(self):
        mtd.clear()


    def _execute_resolution_algorithm(self, **argv):
        default_args = {
            'InstrumentParFile': 'IP0005.par',
            'MonteCarloEvents': 1000
        }
        default_args.update(argv)
        output, mean = VesuvioGeometryEnergyResolution(**default_args)
        return (mean, output)


    def _validate_result_shape(self, mean, resolution):
        """
        Validates the shape of the result tables.
        """

        self.assertTrue(isinstance(mean, ITableWorkspace))
        self.assertEqual(mean.columnCount(), 6)
        self.assertEqual(mean.rowCount(), 24)

        self.assertTrue(isinstance(resolution, ITableWorkspace))
        self.assertEqual(resolution.columnCount(), 17)
        self.assertEqual(resolution.rowCount(), 196)


    def test_resolution(self):
        """
        Check values calculated by resolution algorithm match those expected.
        """
        mean, resolution = self._execute_resolution_algorithm()
        self._validate_result_shape(mean, resolution)

        # Validate mean values
        mean_values = mean.column('Mean')
        TOLERANCE = 7
        self.assertAlmostEqual(mean_values[0],    0.5023026, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[1],    0.9461753, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[2],    0.7906631, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[3],    0.0298165, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[4],    0.0206698, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[5],    0.2581127, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[6],    0.2972636, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[7],    0.0300434, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[8],    0.3971947, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[9],    7.1871166, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[10],   0.4046330, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[11],   7.6269999, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[12],  55.4417038, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[13],  55.4496880, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[14], 140.3843994, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[15],  53.2056618, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[16],  53.2166023, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[17],  31.4454365, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[18],  72.5857315, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[19],  72.5914993, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[20],  45.2145004, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[21], 143.4886322, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[22], 143.4915924, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[23],  97.8484039, places=TOLERANCE)


if __name__ == '__main__':
    unittest.main()
