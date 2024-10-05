# __init__.py

# Import everything from the main RNA_pysoforms module
from . import *

# Create an alias for RNApysoforms
import sys
sys.modules['RNApysoforms'] = sys.modules['RNA_pysoforms']