🧬 ORF 3D Protein Viewer
A lightweight Streamlit app to upload protein FASTA files, detect open reading frames (ORFs), and visualize their predicted 3D structures using ESMFold and py3Dmol for a Prochlorococcus Marinus Genome

📁 Files Overview
main.py – Handles FASTA uploads, ORF detection, and interface logic and genreates orf_report.pdf / longest_orf_protein.fasta from sequence.fasta 

viewer.py – Predicts and displays tertiary structure via ESMFold API on a streamlit web app

sequence.fasta - full Prochlorococcus Marinus Genome

random_protein.fasta – Test fasta sequence

longest_orf_protein.fasta – Auto-generated FASTA file from longest ORF found in the pdf by main.py 

orf_report.pdf – Downloadable PDF report with every ORF found and other details in the marinus sequence
