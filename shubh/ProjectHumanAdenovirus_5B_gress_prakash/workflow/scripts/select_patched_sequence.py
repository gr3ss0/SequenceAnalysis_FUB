
import re
import sys

from Bio import SeqIO

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

records = list(SeqIO.parse(snakemake.input.patched, "fasta"))

if not records:
    raise RuntimeError(f"No sequences found in {snakemake.input.patched}")

scaffold = next((r for r in records if re.match(r"^scf\d+", r.id)), None)

if scaffold is not None:
    print(
        f"Found ragtag patch scaffold '{scaffold.id}' ({len(scaffold.seq)} bp) "
        "- using it as the reference."
    )
    selected = scaffold
else:
    selected = max(records, key=lambda r: len(r.seq))
    print(
        "No ragtag patch scaffold (id starting with 'scf') found - the "
        "de-novo assembly was probably already a single contig. Falling "
        f"back to the longest sequence: '{selected.id}' "
        f"({len(selected.seq)} bp)."
    )



with open(snakemake.output[0], "w") as out_fh:
    SeqIO.write(selected, out_fh, "fasta")
