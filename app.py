import streamlit as st

st.set_page_config(
    page_title="Free Energy Profiles with Replica Exchange MD",
    page_icon="🧬",
    layout="wide",
)

# ================== HEADER ==================
st.markdown(
    """
    <div style="background: linear-gradient(90deg, #1f77b4, #2ca02c); padding: 14px; border-radius: 10px; color: white;">
        <h1 style="margin-bottom: 0.2rem;">Free Energy Profiles with Replica Exchange Molecular Dynamics</h1>
        <p style="margin: 0.2rem 0 0;">
            Using Alanine Dipeptide as a Minimal Model System
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <div style="margin-top: 10px; font-size: 0.95rem;">
        <b>TALIPOV COMPUTATIONAL GROUP, NEW MEXICO STATE UNIVERSITY</b><br>
        Tutorial by: <b>Oluwatosin Saibu, Arslan Talipov, Caroline Ayo-Olaojo, Samuel Brown</b>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("---")

# ================== INTRO: MOTIVATION ==================
left, right = st.columns([1.6, 1.4])

with left:
    st.subheader("Why this tutorial?")

  st.markdown(
    """
    In a standard single-temperature MD run (e.g. 300 K), alanine dipeptide tends to rattle
    inside one or two basins for a long time. Crossings between metastable states in the
    (φ, ψ) landscape can become **rare events**, so any free energy profile you estimate
    from such a trajectory is at risk of being incomplete or biased.

    **Replica Exchange Molecular Dynamics (REMD)** tackles this by simulating multiple
    replicas of the *same* system at different temperatures:
    \\(T_1 < T_2 < \\dots < T_N\\). At regular intervals, neighboring replicas attempt to
    **exchange configurations**, with a Metropolis acceptance rule chosen so that each
    temperature still samples its correct Boltzmann ensemble.

    A well-designed temperature ladder ensures:
    - sufficient **overlap of potential energy distributions** between neighboring replicas,
    - **reasonable exchange probabilities** (typically on the order of 20–40%),
    - and the ability for configurations to **diffuse in temperature space**:
      a structure can move to a higher \\(T\\), cross a barrier in φ/ψ, then return to lower
      \\(T\\) carrying that new conformation with it.

    Practically, this means the lowest-temperature replica (e.g. 300 K) benefits from
    barrier-crossing events that would be extremely rare in plain MD, leading to
    **smoother, better converged 1D and 2D free energy surfaces** along φ and ψ.

    **Metadynamics** provides a complementary route: instead of changing temperature, we
    introduce a time-dependent bias along chosen collective variables (here φ and ψ).
    By depositing Gaussians in visited regions of CV space, metadynamics gradually fills
    free-energy wells and forces exploration of new basins. In the long-time limit,
    the accumulated bias approximates \\(-F(\\phi, \\psi)\\), giving direct access to the
    free energy landscape.

    In this tutorial, alanine dipeptide serves as a clean playground to:
    - construct a sensible REMD temperature ladder,
    - observe how replica exchanges enhance sampling in φ/ψ space,
    - and compare those results with 2D metadynamics-based free energy reconstructions.
    """
)

with right:
    st.subheader("A compact, visual model system")

    st.markdown(
        """
        Alanine dipeptide is small, fast to simulate, and rich enough to display:
        - multiple metastable states,
        - characteristic Ramachandran basins,
        - and clear sensitivity to force field and environment.

        This makes it ideal both for teaching and for testing new methods.
        """,
        unsafe_allow_html=True,
    )

st.markdown("---")

# ================== ALANINE DIPEPTIDE: STRUCTURE ==================
st.subheader("Meet the alanine dipeptide (Ace–Ala–Nme)")

c1, c2 = st.columns([1.3, 1.7])

with c1:
    # Display static structure image (make sure dipep.png is in your repo)
    st.image(
        "dipep.png",
        caption="Alanine dipeptide: a minimal backbone with φ and ψ dihedrals.",
        use_container_width=True,
    )

with c2:
    st.markdown(
        """
        The alanine dipeptide backbone is governed by two dihedral angles:

        - **φ (phi)**: rotation around the N–Cα bond (C(i−1)–N–Cα–C)
        - **ψ (psi)**: rotation around the Cα–C bond (N–Cα–C–N(i+1))

        Together, (φ, ψ) define a **2D conformational landscape** with familiar regions such as
        the α and β basins. Because the molecule is so small, we can:

        - sample this landscape thoroughly,
        - compare vacuum vs solvent,
        - and test different enhanced sampling approaches.

        In the sections below (you will add them next), we will:
        - set up REMD for alanine dipeptide in vacuum,
        - copy a python code to generate and plot the free energy plots
        - or (optionally) load the resulting data into this app to generate clean, publication-ready free energy plots.
        """,
        unsafe_allow_html=True,
    )

st.markdown("---")

# ================== REQUIREMENTS ==================
st.subheader("What you need to follow and reproduce this tutorial")

st.markdown(
    """
    To actually run the simulations on your own machine or cluster and then use this app
    for analysis, you will need a few components. We keep them explicit so everything is reproducible.
    """,
)

st.markdown(
    """
    ### 1. GROMACS with MPI support

   You will use GROMACS for energy minimization, REMD, metadynamics-compatible runs, and basic analysis.

    Please install a recent **GROMACS** version with MPI enabled, following the official guide:
    """,
)
st.markdown(
    """
    ▶️ **GROMACS installation guide**  
    <https://manual.gromacs.org/current/install-guide/index.html>
    """
)

st.markdown(
    """
    ### 2. PLUMED patched into GROMACS

    To bias φ/ψ and analyze complex collective variables, you will need **PLUMED** compiled with GROMACS.

    Useful resources:
    - PLUMED installation: <https://www.plumed.org/doc-v2.9/user-doc/html/_installation.html>  
    - Example of PLUMED + GROMACS usage: <https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-5.html>
    """
)

st.markdown(
    """
    ### 3. Alanine dipeptide starting structure

    We use a capped Ace–Ala–Nme dipeptide (`dipep.pdb`) as the starting point.

    Click below to download the exact structure used in this tutorial:
    """,
)

st.markdown(
    """
    <a href="https://raw.githubusercontent.com/drbiblo/FEMD-tutorial/main/dipep.pdb"
       download="dipep.pdb"
       style="
         display: inline-block;
         padding: 8px 14px;
         margin: 4px 0 10px 0;
         background-color: #1f77b4;
         color: white;
         text-decoration: none;
         border-radius: 6px;
         font-weight: 500;">
       ⬇ Download dipep.pdb
    </a>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    You can inspect this in PyMOL, VMD, ChimeraX, or any viewer to connect the geometry you see
    with the φ/ψ angles and the landscapes we will compute.
    """
)

st.markdown(
    """
    ### 4. Python analysis environment

    This Streamlit app assumes:
    - `streamlit`
    - `numpy`
    - `pandas`
    - `matplotlib`

    (Already listed in your `requirements.txt`.)
    """
)

st.markdown(
    """
    ### 5. What this app will help you do

    In the next sections (which you’ll add under this header), this app will:

    - Show ready-to-use **REMD** and **metadynamics** setup snippets for alanine dipeptide.
    - Let you **upload**:
      - dihedral time series (φ, ψ),
      - PLUMED `fes_psi_phi` outputs,
      - or COLVAR-style data.
    - Automatically generate:
      - 2D FES(φ, ψ) heatmaps,
      - 1D profiles F(φ) and F(ψ),
      - figures clean enough for lab meetings and slides.

    This page is your **landing + context**. Below it, you can now start adding:
    - a “Simulation setup” section,
    - an “Upload & visualize” section,
    - and side-by-side comparisons (vacuum vs solvent, force field changes, etc.).
    """
)
