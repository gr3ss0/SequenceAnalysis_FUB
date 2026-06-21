
import sys
from collections import defaultdict

from Bio import SeqIO

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]

summed_bitscore = defaultdict(float)
n_hits = defaultdict(int)

with open(snakemake.input.blast) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        fields = line.split("\t")
        row = dict(zip(BLAST_COLUMNS, fields))
        sseqid = row["sseqid"]
        summed_bitscore[sseqid] += float(row["bitscore"])
        n_hits[sseqid] += 1

if not summed_bitscore:
    raise RuntimeError(
        f"No BLAST hits found in {snakemake.input.blast}. None of the "
        f"candidate references in {snakemake.input.candidates} are similar "
        "enough to this sample's de-novo contigs. Consider adding a more "
        "closely related candidate reference, or relaxing "
        "assembly_params.blast.evalue in the config."
    )

best_id = max(summed_bitscore, key=summed_bitscore.get)

print("Candidate reference summed bitscores (sseqid: bitscore, n_hits):")
for sseqid, score in sorted(summed_bitscore.items(), key=lambda kv: kv[1], reverse=True):
    print(f"  {sseqid}\t{score:.1f}\t{n_hits[sseqid]} hit(s)")
print(f"\nSelected best reference: {best_id}")

best_record = None
for record in SeqIO.parse(snakemake.input.candidates, "fasta"):
    if record.id == best_id:
        best_record = record
        break

if best_record is None:
    raise RuntimeError(
        f"Could not find a sequence with id '{best_id}' in "
        f"{snakemake.input.candidates}. Check that the BLAST sseqid matches "
        "the fasta header (first whitespace-delimited token) of the "
        "candidate reference file."
    )

with open(snakemake.output.fasta, "w") as out_fa:
    SeqIO.write(best_record, out_fa, "fasta")

with open(snakemake.output.info, "w") as out_info:
    out_info.write(f"sample\t{snakemake.wildcards.sample}\n")
    out_info.write(f"selected_reference\t{best_id}\n")
    out_info.write(f"summed_bitscore\t{summed_bitscore[best_id]:.1f}\n")
    out_info.write(f"n_hits\t{n_hits[best_id]}\n")
