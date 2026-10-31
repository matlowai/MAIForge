"""Compatibility import for the Norway-spruce calibration fixture."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from development.annual_woody import WoodyAnnualBackend, digest
Growth=WoodyAnnualBackend
