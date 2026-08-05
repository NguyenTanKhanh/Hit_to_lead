#!/usr/bin/env bash
set -euo pipefail

# Get the script directory (works no matter where you run it from)
SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
RESULTS_DIR="$SCRIPT_DIR/results/lead_generation_output"

# Fixed input/output
IN="$RESULTS_DIR/generated_smiles_pass_pH7_4.smi"          # input: output from pH 7.4 step
OUT="$RESULTS_DIR/generated_smiles_pass_pH6_3.smi"         # output after pH 6.3 filtering
OUT_PROPS="$RESULTS_DIR/generated_smiles_pass_pH6_3_with_props.smi"  # report ALL input + props

# Create directory if it doesn't exist
mkdir -p "$RESULTS_DIR"

# ---- Filtering thresholds ----
MW_LO=200; MW_HI=650
LOGP_LO=0; LOGP_HI=5
TPSA_LO=0;  TPSA_HI=200
HBD_LO=0;    HBD_HI=5
HBA_LO=0;    HBA_HI=10
ROTB_LO=0;   ROTB_HI=10
LOGS_LO=-10;  LOGS_HI=0  # logS range (ESOL model)
# ----------------------------

# Check input file exists
if [[ ! -f "$IN" ]]; then
  echo "ERROR: Input not found: $IN" >&2
  exit 1
fi

# --- Main filtering step (unchanged structure) ---
python utils/compound_properties_filter.py \
  --mode filter \
  --in "$IN" \
  --out "$OUT" \
  --ph 6.3 \
  --mw   "$MW_LO" "$MW_HI" \
  --logp "$LOGP_LO" "$LOGP_HI" \
  --tpsa "$TPSA_LO" "$TPSA_HI" \
  --hbd  "$HBD_LO" "$HBD_HI" \
  --hba  "$HBA_LO" "$HBA_HI" \
  --rotb "$ROTB_LO" "$ROTB_HI" \
  --logs "$LOGS_LO" "$LOGS_HI"

echo "[pH 6.4] Passed: $(wc -l < "$OUT" | tr -d ' ') molecules -> $OUT"

# === ADD-ON ONLY: report ALL INPUT molecules with properties computed by the SAME code as filtering ===
# Output includes ALL molecules from $IN (not only those that passed),
# and uses the same protonation + rdkit_features + logS implementation as the filter.
python - <<'PY' "$SCRIPT_DIR" "$IN" "$OUT_PROPS"
import sys
from pathlib import Path

script_dir, inp, outp = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, str(Path(script_dir).resolve()))

from utils.compound_properties_filter import rdkit_features, protonate_smiles

def parse_line(line: str):
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    low = line.lower()
    if low.startswith("smiles") or low.startswith("mol_id") or low.startswith("id"):
        return None
    parts = line.split()
    smi = parts[0].strip()
    if smi.lower() == "smiles":
        return None
    mol_id = parts[1].strip() if len(parts) > 1 else "mol"
    return smi, mol_id

PH = 6.3  # must match the filter pH above

with open(inp, "r", encoding="utf-8") as f_in, open(outp, "w", encoding="utf-8") as f_out:
    f_out.write("#SMILES\tmol_ID\tMW\tLogP\tTPSA\tHBD\tHBA\tRotBonds\tLogS\n")

    for line in f_in:
        parsed = parse_line(line)
        if parsed is None:
            continue
        smi, mol_id = parsed

        try:
            prot = protonate_smiles(smi, ph=PH)
        except Exception:
            continue

        feats = rdkit_features(prot)
        if not feats:
            continue

        f_out.write(
            f"{prot}\t{mol_id}\t"
            f"{feats['MW']:.3f}\t{feats['LogP']:.3f}\t{feats['TPSA']:.3f}\t"
            f"{int(feats['HBD'])}\t{int(feats['HBA'])}\t{int(feats['RotBonds'])}\t"
            f"{feats['LogS']:.3f}\n"
        )
PY

echo "[pH 6.4] Props file (ALL input) -> $OUT_PROPS"
