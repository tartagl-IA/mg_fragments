"""Settings module for the ChEMBL database and fragments generation."""

import os

# --- database settings ---
DATABASE_NAME = "chembl.db"

# --- fragments generation ---
FRAGMENTS_MIN_ATOMS = 5
FRAGMENTS_MAX_ATOMS = 100
FRAGMENTS_FLEXIBILITY = "rigid"  # 'rigid' or 'flexible'
FRAGMENTS_MAX_ROTABLE_BONDS = 1
FRAGMENTS_OUTPUT_DIR = os.path.join("outputs", "fragments")
FRAGMENTS_TOP_RES = 20
REACTIVE_PATTERN_FILE = os.path.join("chem", "reactive_patterns.json")
