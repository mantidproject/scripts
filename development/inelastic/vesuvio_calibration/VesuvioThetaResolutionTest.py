#pylint: disable=too-many-public-methods

import unittest

from mantid.api import ITableWorkspace
from mantid.simpleapi import *


class VesuvioThetaResolutionTest(unittest.TestCase):

    def setUp(self):
        self._default_args = {}
        self._default_args['Samples'] = '17083-17084'
        self._default_args['Background'] = '17086'
        self._default_args['DSpacings'] = '1.489,1.750,2.475,2.858'
        self._default_args['DeltaL0'] = 2.2e-3
        self._default_args['DeltaL1'] = 2.3e-3
        self._default_args['DeltaT0'] = 3.7e-1
        self._default_args['InstrumentParFile'] = 'IP0005.par'


    def tearDown(self):
        mtd.clear()


    def _execute_resolution_algorithm(self, **argv):
        self._default_args.update(argv)
        forward, back, mean = VesuvioThetaResolution(**self._default_args)
        return (forward, back, mean)


    def _validate_result_shape(self, forward, back, mean):
        """
        Validates the shape of the result tables.
        """

        self.assertTrue(isinstance(forward, ITableWorkspace))
        self.assertEqual(forward.columnCount(), 3)
        self.assertEqual(forward.rowCount(), 64)

        self.assertTrue(isinstance(back, ITableWorkspace))
        self.assertEqual(back.columnCount(), 3)
        self.assertEqual(back.rowCount(), 132)

        self.assertTrue(isinstance(mean, ITableWorkspace))
        self.assertEqual(mean.rowCount(), 4)


    def test_run_with_par_file(self):
        """
        Tests running using distances from an instrument parameter file.
        """
        forward, back, mean = self._execute_resolution_algorithm()

        self._validate_result_shape(forward, back, mean)

        # Validate mean values
        mean_values = mean.column('Mean')
        TOLERANCE = 7
        self.assertAlmostEqual(mean_values[0], 0.9915794, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[1], 4.4192872, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[2], 0.9501636, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[3], 5.0413828, places=TOLERANCE)


if __name__ == '__main__':
    unittest.main()
