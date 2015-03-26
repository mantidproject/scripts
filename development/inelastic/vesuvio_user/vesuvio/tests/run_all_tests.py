"""Run all unit tests
"""
import os
import unittest

# unittest discovery was only introduced in Python 2.7 so we'll do it ourselves
module_dir = os.path.dirname(__file__)
test_files = os.listdir(module_dir)
test_files.remove("__init__.py")
test_files.remove("run_all_tests.py")
test_modules = []
for filename in test_files:
    if filename.endswith(".py"):
        test_modules.append(filename.rstrip(".py"))

test_loader = unittest.defaultTestLoader
suite = unittest.TestSuite()
for name in test_modules:
    suite.addTests(test_loader.loadTestsFromName("vesuvio.tests.{0}".format(name)))

unittest.TextTestRunner(verbosity=2).run(suite)
