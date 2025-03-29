"""Settings module for the ChEMBL database and fragments generation."""

# --- database settings ---
DATABASE_NAME = "chembl.db"

# --- fragments generation ---
FRAGMENTS_MIN_ATOMS = 5
FRAGMENTS_MAX_ATOMS = 100
FRAGMENTS_FLEXIBILITY = "rigid"  # 'rigid' or 'flexible'
FRAGMENTS_MAX_ROTABLE_BONDS = 1
