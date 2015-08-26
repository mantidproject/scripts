#pylint: disable=too-many-public-methods

import unittest

from mantid.api import (ITableWorkspace, WorkspaceGroup)
from mantid.simpleapi import *


class VesuvioEnergyResolutionTest(unittest.TestCase):

    def setUp(self):
        self._default_args = {}
        self._default_args['Samples'] = '17083-17084'
        self._default_args['DeltaL0'] = 2.2e-3
        self._default_args['DeltaL1'] = 2.3e-3
        self._default_args['DeltaT0'] = 3.7e-1
        self._default_args['InstrumentParFile'] = 'IP0005.par'


    def tearDown(self):
        mtd.clear()


    def _execute_resolution_algorithm(self, **argv):
        self._default_args.update(argv)
        output, mean = VesuvioEnergyResolution(**self._default_args)
        return (mean, output)


    def _validate_result_shape(self, mean, output):
        """
        Validates the shape of the result tables.
        """

        self.assertTrue(isinstance(mean, ITableWorkspace))
        self.assertEqual(mean.columnCount(), 6)
        self.assertEqual(mean.rowCount(), 8)

        self.assertTrue(isinstance(output, WorkspaceGroup))
        self.assertEqual(len(output), 3)

        for idx in range(len(output)):
            wks = output.getItem(idx)
            self.assertTrue(isinstance(wks, ITableWorkspace))
            self.assertEqual(wks.columnCount(), 3)
            self.assertEqual(wks.rowCount(), 196)


    def test_run_with_par_file(self):
        """
        Tests running using distances from an instrument parameter file.
        """
        mean, output = self._execute_resolution_algorithm(InstrumentParFile='IP0005.par')

        self._validate_result_shape(mean, output)

        # Validate mean values
        mean_values = mean.column('Mean')
        TOLERANCE = 7
        self.assertAlmostEqual(mean_values[0],  55.3837471, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[1], 138.2033997, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[2],  53.2623253, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[3],  30.9145622, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[4],  74.1999588, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[5],  61.5803642, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[6], 143.5410919, places=TOLERANCE)
        self.assertAlmostEqual(mean_values[7], 108.9285660, places=TOLERANCE)


if __name__ == '__main__':
    unittest.main()
