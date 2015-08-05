#!/usr/bin/env python

"""
Usage: ./run_test.py [mantid bin directory] [test module name]
"""

import sys
import unittest

sys.path.append(sys.argv[1])

test_loader = unittest.defaultTestLoader
suite = unittest.TestSuite()
suite.addTests(test_loader.loadTestsFromName(sys.argv[2]))
unittest.TextTestRunner(verbosity=2).run(suite)
