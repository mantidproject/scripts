"""Supports the Vesuvio instrument at ISIS

The module requires Mantid and importing it registers additional algorithms into
the framework
"""
from __future__ import absolute_import

import vesuvio.workflow as workflow

import mantid
# Register algorithms
import vesuvio.algorithms