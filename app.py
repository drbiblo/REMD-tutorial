import streamlit as st

st.set_page_config(
    page_title="Alanine Dipeptide Free Energy Landscape Tutorial",
    page_icon="🧬",
    layout="wide",
)

# ================== HEADER / AUTHORS ==================
st.title("Alanine Dipeptide Free Energy Landscape Tutorial")

st.markdown("""
**TALIPOV COMPUTATIONAL GROUP, NEW MEXICO STATE UNIVERSITY**  
*Tutorial by:* **Oluwatosin Saibu, Arslan Talipov, Caroline Ayo-Olaojo, Samuel Brown**
""")

st.markdown("---")

# ================== INTRODUCTION ==================
st.subheader("Why this tutorial?")

st.markdown("""
Understanding **free energy landscapes** is key to explaining:

- Conformational preferences  
- Folding/unfolding  
- Binding and recognition  
- Rare transitions in biomolecules  

Plain MD at a single temperature often:

- Gets trapped in local minima  
- Samples transitions too slowly  
- Needs very long simulations to visit all relevant states  

Here we use **alanine dipeptide** as a minimal, clean system to:

- Learn enhanced sampling concepts  
- Set up **Replica Exchange MD (REMD)** and 2D metadynamics (φ, ψ)  
- Reconstruct **1D and 2D free energy surfaces** along φ and ψ
""")

st.markdown("""
**Why enhanced sampling (REMD / MetaD) instead of only plain MD?**

- 🚀 **Better barrier crossing:** high-T replicas / bias help escape local minima.  
- 🎯 **More reliable FES:** multiple basins are sampled → smoother, converged profiles.  
- 📉 **Less hysteresis:** outcome depends less on the starting structure.  
- 🧪 **Generalizable:** the same ideas extend to real peptides, proteins, and ligands.
""")

st.markdown("---")

# ================== ALANINE DIPEPTIDE OVERVIEW ==================
st.subheader("The model system: Alanine dipeptide")

col1, col2 = st.columns([1.2, 1])

with col1:
    st.markdown("""
Alanine dipeptide (Ace–Ala–Nme) is the **“hello world”** of conformational sampling:

- Small and cheap to simulate  
- Backbone described by two key dihedrals:
  - **φ (phi)**: C(i-1)–N–Cα–C
  - **ψ (psi)**: N–Cα–C–N(i+1)
- These define a 2D landscape (Ramachandran-style) with well-known basins:
  - α, β, etc.
- Perfect for:
  - testing force fields,
  - validating enhanced sampling,
  - comparing solvent vs vacuum behavior.

In this tutorial you will:

1. Build alanine dipeptide
2. Run REMD / 2D MetaD (on your own machine/cluster)
3. Load outputs here to visualize φ/ψ distributions and free energy surfaces
""")

with col2:
    st.markdown("**Alanine dipeptide structure**")
    st.markdown("""
You can download the starting structure used in this tutorial:

👉 **[Download dipeptide PDB](https://raw.githubusercontent.com/drbiblo/FEMD-tutorial/main/dipep.pdb)**

View it locally in PyMOL / VMD, or integrate directly into your workflow.
""")

st.markdown("---")

# ================== REQUIREMENTS ==================
st.subheader("What you need to follow this tutorial")

st.markdown("""
To **run the simulations** and reproduce the analysis, you will need:

1. **GROMACS with MPI support**

   For parallel and REMD runs.  
   Official installation guide (build + MPI):
   - https://manual.gromacs.org/current/install-guide/index.html

2. **PLUMED patched into GROMACS** (for metadynamics & CVs)

   - PLUMED install & patching docs:
     - https://www.plumed.org/doc-v2.9/user-doc/html/_installation.html
   - Example notes for using PLUMED with GROMACS:
     - https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-5.html

3. **Alanine dipeptide structure**

   - Use the provided Ace–Ala–Nme:
     - **dipep.pdb** from this link:  
       👉 https://raw.githubusercontent.com/drbiblo/FEMD-tutorial/main/dipep.pdb

4. **Python environment for analysis**

   - This app (or your local env) uses:
     - `streamlit`, `numpy`, `pandas`, `matplotlib`

5. **Basic command-line/HPC familiarity**

   - Editing `.mdp` files  
   - Running `gmx` commands  
   - Submitting jobs on a cluster

---

Below, we will:

- Show ready-to-use REMD and MetaD setup blocks  
- Let you **upload your own** `COLVAR` / FES / dihedral files  
- Generate clean 1D and 2D FES plots suitable for slides and papers
""")
