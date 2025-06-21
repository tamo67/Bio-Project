# main.py

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter

# Setup
START_CODON = "TGA"
STOP_CODONS = {"TAA", "TTA", "TAG"}
VALID_AAS = set("ARNDCEQGHILKMFPSTWYV")
table_11 = CodonTable.unambiguous_dna_by_id[11]

three_letter = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
    'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
    'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
    'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Stop'
}

def find_orfs_with_translation(seq):
    seq = str(seq).upper()
    orfs = []
    for frame in range(3):
        i = frame
        while i <= len(seq) - 3:
            codon = seq[i:i+3]
            if codon == START_CODON:
                for j in range(i+3, len(seq)-2, 3):
                    stop = seq[j:j+3]
                    if stop in STOP_CODONS:
                        orf_nt = seq[i:j+3]
                        if len(orf_nt) % 3 == 0:
                            try:
                                protein_1 = str(Seq(orf_nt).translate(table=11, to_stop=False))
                                protein_3 = '-'.join(three_letter.get(aa, '???') for aa in protein_1)
                                orfs.append({
                                    "frame": frame + 1,
                                    "start": i+1,
                                    "end": j+3,
                                    "length": len(orf_nt),
                                    "nt_seq": orf_nt,
                                    "aa_seq_1letter": protein_1,
                                    "aa_seq_3letter": protein_3
                                })
                            except:
                                pass
                        break
            i += 3
    return orfs

# Run
fasta_path = "sequence.fasta"
output_folder = "orf_outputs"
os.makedirs(output_folder, exist_ok=True)

record = SeqIO.read(fasta_path, "fasta")
orfs = find_orfs_with_translation(record.seq)

# Save PDF
pdf_path = os.path.join(output_folder, "orf_report.pdf")
c = canvas.Canvas(pdf_path, pagesize=letter)
y = 750
line_height = 14

c.drawString(40, y, f"Sequence ID: {record.id}")
y -= line_height

for idx, orf in enumerate(orfs, 1):
    lines = [
        f"ORF {idx}",
        f"Frame: {orf['frame']}",
        f"Start: {orf['start']}, End: {orf['end']}",
        f"Length: {orf['length']} bases",
        f"Nucleotide sequence: {orf['nt_seq']}",
        f"Protein (1-letter): {orf['aa_seq_1letter']}",
        f"Protein (3-letter): {orf['aa_seq_3letter']}",
        ""
    ]
    for line in lines:
        if y <= 50:
            c.showPage()
            y = 750
        c.drawString(40, y, line)
        y -= line_height
c.save()
print(f"✅ ORF report saved to {pdf_path}")

# Save longest ORF to FASTA
if orfs:
    longest = max(orfs, key=lambda o: len(o['aa_seq_1letter']))
    cleaned_seq = ''.join([aa for aa in longest['aa_seq_1letter'] if aa in VALID_AAS])

    fasta_out = os.path.join(output_folder, "longest_orf_protein.fasta")
    with open(fasta_out, "w") as f:
        f.write(f">Longest_ORF_Protein\n{cleaned_seq}\n")
    print(f"✅ Longest ORF saved to {fasta_out}")
else:
    print("⚠️ No ORFs found.")
