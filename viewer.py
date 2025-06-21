import streamlit as st
import io
from Bio import SeqIO
import requests
import py3Dmol

st.set_page_config(page_title="🧬 ORF 3D Protein Viewer Tertiary Structure")
st.title("🧬 ORF 3D Protein Viewer Tertiary Structure")

st.subheader("Step 1: Upload a Protein FASTA File")
uploaded_file = st.file_uploader("Choose a .fasta file containing a protein sequence", type=["fasta", "fa"])

if uploaded_file:
    try:
        record = SeqIO.read(io.TextIOWrapper(uploaded_file, encoding="utf-8"), "fasta")
        seq = str(record.seq).upper()
        valid_aas = set("ARNDCEQGHILKMFPSTWYV")
        cleaned_seq = ''.join([aa for aa in seq if aa in valid_aas])

        st.success(f"✅ Sequence loaded ({len(cleaned_seq)} AAs)")
        st.code(cleaned_seq)

        if len(cleaned_seq) < 20:
            st.error("❌ Sequence is too short for prediction. Minimum 20 amino acids required.")
        else:
            st.info("⏳ Predicting structure using ESMFold...")
            response = requests.post(
                "https://api.esmatlas.com/foldSequence/v1/pdb/",
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                data=cleaned_seq
            )

            if response.ok:
                pdb = response.text
                view = py3Dmol.view(width=700, height=600)
                view.addModel(pdb, "pdb")
                view.setStyle({"cartoon": {"color": "spectrum"}})
                view.zoomTo()
                st.components.v1.html(view._make_html(), height=600, width=700)

                st.subheader("🧬 Protein Sequence")
                st.code(cleaned_seq)
            else:
                st.error("❌ Structure prediction failed. Try a different sequence.")
    except Exception as e:
        st.error(f"❌ Error: {str(e)}")
else:
    st.warning("📁 Please upload a FASTA file to begin.")
