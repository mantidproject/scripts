import unittest

from vesuvio.fitting import FittingOptions
from vesuvio.profiles import GramCharlierMassProfile, GaussianMassProfile

class FittingOptionsTest(unittest.TestCase):

    def test_function_str_with_no_given_params_looks_as_expected(self):
        fit_opts = self._create_test_fitting_opts()

        expected = \
            "composite=ComptonScatteringCountRate,NumDeriv=1,IntensityConstraints=\"Matrix(1|2)1.000000|-4.000000\";"\
            "name=GramCharlierComptonProfile,Mass=1.0079,HermiteCoeffs=1 0 0,Width=5;"\
            "name=GaussianComptonProfile,Mass=16,Width=10"
        self.assertEqual(expected, fit_opts.create_function_str())

    def test_function_str_with_given_params_looks_as_expected(self):
        fit_opts = self._create_test_fitting_opts()

        param_vals = {"f0.Width": 7.5, "f0.FSECoeff": 0.1, "f0.C_0": 0.25,
                      "f0.C_2": 0.5, "f0.C_4": 0.75}
        param_vals.update({"f1.Width": 11.0, "f1.Intensity": 4.5})

        expected = \
            "composite=CompositeFunction,NumDeriv=1;"\
            "name=GramCharlierComptonProfile,Mass=1.0079,HermiteCoeffs=1 0 0,Width=7.5,FSECoeff=0.100000,C_0=0.250000;"\
            "name=GaussianComptonProfile,Mass=16,Width=11.0,Intensity=4.5"
        self.assertEqual(expected, fit_opts.create_function_str(param_vals))

    def _create_test_fitting_opts(self):
        gramc = GramCharlierMassProfile([2, 5, 7], 1.0079, [1, 0, 0], 1, 1)
        gauss = GaussianMassProfile(10, 16)
        constraints = ([1, -4],)

        return FittingOptions([gramc, gauss], constraints)


if __name__ == '__main__':
    unittest.main()
