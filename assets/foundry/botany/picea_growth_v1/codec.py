"""Compatibility import for the shared constrained trait-codec sampler."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from development.trait_codec import sample_genome
