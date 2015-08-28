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


    def test_run_with_par_file(self):
        """
        Tests running using distances from an instrument parameter file.
        """
        mean, resolution = self._execute_resolution_algorithm()
        self._validate_result_shape(mean, resolution)

        # Validate mean values
        #TODO


if __name__ == '__main__':
    unittest.main()
