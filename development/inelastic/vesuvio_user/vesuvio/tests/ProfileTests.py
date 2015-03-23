import unittest

from vesuvio.profiles import (create_profile_from_str, GramCharlierMassProfile,
                          GaussianMassProfile)

# --------------------------------------------------------------------------------
# Gaussian
# --------------------------------------------------------------------------------

class GaussianMassProfileTest(unittest.TestCase):

    # ---------------- Success cases ---------------------------

    def test_string_with_fixed_width_produces_valid_object(self):
        function_str = "function=Gaussian,width=10"
        mass = 16.0

        profile = create_profile_from_str(function_str, mass)
        self.assertTrue(isinstance(profile, GaussianMassProfile))
        self.assertAlmostEqual(mass, profile.mass)
        self.assertAlmostEqual(10.0, profile.width)

    def test_string_with_constrained_width_produces_valid_object(self):
        function_str = "function=Gaussian,width=[2, 5, 7]"
        mass = 16.0

        profile = create_profile_from_str(function_str, mass)
        self.assertTrue(isinstance(profile, GaussianMassProfile))
        self.assertAlmostEqual(mass, profile.mass)
        self.assertEqual([2, 5, 7], profile.width)

    def test_function_string_has_expected_form_with_no_defaults(self):
        test_profiles = GaussianMassProfile(10, 16)

        expected = "name=GaussianComptonProfile,Mass=16,Width=10;"
        self.assertEqual(expected, test_profiles.create_fitting_str())

    def test_function_string_has_expected_form_with_defaults_given(self):
        test_profiles = GaussianMassProfile(10, 16)
        param_prefix = "f1."
        param_vals = {"f1.Width": 11.0, "f1.Intensity": 4.5}

        expected = "name=GaussianComptonProfile,Mass=16,Width=11.0,Intensity=4.5;"
        self.assertEqual(expected, test_profiles.create_fitting_str(param_vals, param_prefix))

    # ---------------- Failure cases ---------------------------

    def test_string_not_starting_with_function_equals_name_gives_error(self):
        function_str = "function=Gaussia,width=[2, 5, 7]"
        mass = 16.0

        self.assertRaises(TypeError, GaussianMassProfile.from_str,
                          function_str, mass)

    def test_string_not_starting_with_function_gives_error(self):
        function_str = "Gaussian,width=[2, 5, 7]"
        mass = 16.0

        self.assertRaises(TypeError, GaussianMassProfile.from_str,
                          function_str, mass)

    def test_string_with_wrong_function_gives_error(self):
        function_str = "function=GramCharlier,width=[2, 5, 7]"
        mass = 16.0

        self.assertRaises(TypeError, GaussianMassProfile.from_str,
                          function_str, mass)

# --------------------------------------------------------------------------------
# GramCharlier
# --------------------------------------------------------------------------------

class GramCharlierMassProfileTest(unittest.TestCase):

    def test_string_with_fixed_width_produces_valid_object(self):
        function_str = "function=GramCharlier,width=[2, 5,7],k_free=1,hermite_coeffs=[1,0,1],sears_flag=0,"
        mass = 16.0

        profile = create_profile_from_str(function_str, mass)
        self.assertTrue(isinstance(profile, GramCharlierMassProfile))
        self.assertAlmostEqual(mass, profile.mass)
        self.assertEqual([2, 5, 7], profile.width)
        self.assertEqual([1, 0, 1], profile.hermite_co)
        self.assertEqual(0, profile.sears_flag)
        self.assertEqual(1, profile.k_free)

    def test_function_string_has_expected_form_with_no_defaults(self):
        test_profile = GramCharlierMassProfile(10, 16,[1,0,1],1,1)

        expected = "name=GramCharlierComptonProfile,Mass=16,HermiteCoeffs=1 0 1,Width=10;"
        self.assertEqual(expected, test_profile.create_fitting_str())

    def test_function_string_has_expected_form_with_given_values(self):
        test_profile = GramCharlierMassProfile(10, 16,[1,0,1],1,1)
        param_prefix = "f1."
        param_vals = {"f1.Width": 11.0, "f1.FSECoeff": 0.1, "f1.C_0": 0.25,
                      "f1.C_2": 0.5, "f1.C_4": 0.75}

        expected = "name=GramCharlierComptonProfile,Mass=16,HermiteCoeffs=1 0 1,"\
                   "Width=11.0,FSECoeff=0.100000,C_0=0.250000,C_4=0.750000;"
        self.assertEqual(expected, test_profile.create_fitting_str(param_vals, param_prefix))

    # ---------------- Failure cases ---------------------------

    def test_string_not_starting_with_function_equals_name_gives_error(self):
        function_str = "function=GramCharlie,width=[2, 5, 7]"
        mass = 16.0

        self.assertRaises(TypeError, GramCharlierMassProfile.from_str,
                          function_str, mass)

    def test_string_not_starting_with_function_gives_error(self):
        function_str = "GramCharlie,width=[2, 5, 7]"
        mass = 16.0

        self.assertRaises(TypeError, GramCharlierMassProfile.from_str,
                          function_str, mass)

    def test_string_with_wrong_function_gives_error(self):
        function_str = "function=Gaussian,width=[2, 5, 7]"
        mass = 16.0

        self.assertRaises(TypeError, GramCharlierMassProfile.from_str,
                          function_str, mass)


if __name__ == '__main__':
    unittest.main()
