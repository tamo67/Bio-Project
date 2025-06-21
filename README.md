# 🧬 ORF 3D Protein Viewer

A lightweight **Streamlit** app to upload protein FASTA files, detect **open reading frames (ORFs)**, and visualize their predicted **3D structures** using **ESMFold** and **py3Dmol**, specifically for a *Prochlorococcus Marinus* genome.

---

## 📁 Files Overview

- `main.py` – Handles FASTA uploads, ORF detection, and interface logic.  
  Generates `orf_report.pdf` and `longest_orf_protein.fasta` from `sequence.fasta`.

- `viewer.py` – Predicts and displays **tertiary structure** via **ESMFold API** on a Streamlit web app.

- `sequence.fasta` – Full *Prochlorococcus Marinus* genome input.

- `random_protein.fasta` – Test protein sequence in FASTA format.

- `longest_orf_protein.fasta` – Auto-generated FASTA file of the **longest ORF** found in the genome.

- `orf_report.pdf` – Downloadable PDF report with every ORF found, including positions, translations, and sequence details.

---
