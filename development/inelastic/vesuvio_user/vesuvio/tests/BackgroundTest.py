import unittest

from vesuvio.backgrounds import PolynomialBackground

# --------------------------------------------------------------------------------
# Polynomial
# --------------------------------------------------------------------------------

class PolynomialBackgroundTest(unittest.TestCase):

    def test_create_function_str_for_nth_order_with_no_values(self):
        background = PolynomialBackground(order=2)

        expected = "name=Polynomial,n=2"
        self.assertEqual(expected, background.create_fit_function_str())

    def test_create_function_str_for_nth_order_given_fixed_values(self):
        background = PolynomialBackground(order=2)
        param_values = {"f1.A0": 2.0, "f1.A1": 3.0, "f1.A2": 4.0}

        expected = "name=Polynomial,n=2,A0=2.0,A1=3.0,A2=4.0"
        self.assertEqual(expected, background.create_fit_function_str(param_values, param_prefix="f1."))

if __name__ == '__main__':
    unittest.main()
