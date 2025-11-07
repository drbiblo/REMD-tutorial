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

# ================== INTRO & BACKGROUND ==================
left, right = st.columns([1.7, 1.3])

with left:
    st.subheader("Why this tutorial?")

    st.markdown(
        """
        Molecular dynamics (MD) is a powerful tool for exploring conformational dynamics,
        but **plain single-temperature MD** has a well-known limitation:
        the system rapidly falls into one or two low-energy basins and then spends most of
        the trajectory vibrating there.  
        Transitions that require crossing even moderate free-energy barriers become rare
        on accessible timescales. As a result:

        - free energy profiles \(F(\phi)\), \(F(\psi)\), or \(F(\phi,\psi)\) built from a single
          300 K trajectory are often **noisy, biased, or incomplete**;
        - important conformational states may be severely undersampled or completely missed.
        """,
        unsafe_allow_html=True,
    )

with right:
    st.subheader("What this tutorial delivers")

    st.markdown(
        """
        This tutorial walks through:

        - constructing a minimal alanine dipeptide model,
        - setting up **Replica Exchange MD (REMD)** in vacuum,
        - applying **2D metadynamics** along φ/ψ,
        - and reconstructing **1D & 2D free energy surfaces** suitable for analysis and slides.

        All using a system small enough to understand in detail, but rich enough
        to show real enhanced sampling behavior.
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---- Enhanced sampling background ----
st.subheader("Why Replica Exchange MD and Metadynamics?")

st.markdown(
    """
    In a standard fixed-temperature simulation, the sampling problem is fundamentally kinetic:
    barrier crossings are rare, so your estimate of the free energy surface is built from a
    trajectory that has not seen enough of phase space.

    **Replica Exchange Molecular Dynamics (REMD)** tackles this by running multiple replicas
    of the *same* system at different temperatures
    \\(T_1 < T_2 < \\dots < T_N\\) simultaneously.

    - High-temperature replicas explore conformational space broadly and cross barriers more easily.
    - Low-temperature replicas retain sharp, physically relevant distributions.
    - At fixed intervals, neighboring replicas attempt to **exchange configurations**, using a
      Metropolis criterion chosen so that each temperature maintains its correct Boltzmann ensemble.

    If the temperature ladder is chosen so that neighboring replicas have overlapping potential
    energy distributions, the exchange probabilities (typically on the order of a few tenths)
    are high enough that configurations can perform a **random walk in temperature space**:

    - a configuration heats up → crosses a barrier in \\((\\phi, \\psi)\\),
    - then cools back down → bringing new conformations into the 300 K ensemble.

    The payoff is that your **lowest-temperature replica** inherits transitions that would be
    extremely rare in plain MD, yielding **smoother, better-converged 1D and 2D free energy
    surfaces** along \\(\\phi\\) and \\(\\psi\\).

    **Metadynamics** approaches the same problem from a different direction. Instead of changing
    temperature, we:

    - identify slow **collective variables** (here \\(\\phi\\) and \\(\\psi\\)),
    - periodically deposit repulsive Gaussian hills where the system has visited in CV space,
    - progressively fill free-energy wells, forcing exploration of new regions.

    Under appropriate settings (e.g. well-tempered metadynamics), the accumulated bias
    approximates \\(-F(\\phi, \\psi)\\), so you can reconstruct the free energy surface directly.
    For alanine dipeptide, biasing \\(\\phi\\) and \\(\\psi\\) is natural and highly targeted.

    In this tutorial, alanine dipeptide provides a controlled playground to:

    - build and test a sensible REMD temperature ladder,
    - observe how replica exchanges enhance sampling in \\(\\phi/\\psi\\) space,
    - and compare these results with metadynamics-based reconstructions of the same landscape.
    """
)

st.markdown("---")

# ================== ALANINE DIPEPTIDE: STRUCTURE ==================
st.subheader("The model system: Alanine dipeptide (Ace–Ala–Nme)")

c1, c2 = st.columns([1.1, 1.9])

with c1:
    # Visible figure in the body (ensure dipep.png exists in the repo root)
    st.image(
        "dipep.png",
        caption="Alanine dipeptide: minimal backbone with φ and ψ dihedrals.",
        use_container_width=True,
    )

with c2:
    st.markdown(
        """
        Alanine dipeptide is the canonical **testbed for conformational sampling**.

        Its backbone is governed by two dihedral angles:

        - **φ (phi)**: rotation about the N–Cα bond (C(i−1)–N–Cα–C)
        - **ψ (psi)**: rotation about the Cα–C bond (N–Cα–C–N(i+1))

        Together, \\((\\phi, \\psi)\\) define a compact yet non-trivial conformational landscape
        with distinct minima (α, β, etc.) visible on a **Ramachandran-like 2D free energy surface**.

        Why this system?

        - It is small enough that you can afford enhanced sampling methods.
        - The φ/ψ landscape is well-characterized in the literature, so you can check if your
          force field + sampling protocol behave sensibly.
        - Subtle differences (e.g. vacuum vs solvent, different force fields, different biasing
          schemes) show up clearly in the free energy profiles.

        All of this makes alanine dipeptide an ideal didactic and diagnostic system for learning
        **Replica Exchange MD** and **metadynamics** before scaling to larger biomolecules.
        """,
        unsafe_allow_html=True,
    )

st.markdown("---")

# ================== REQUIREMENTS ==================
st.subheader("Reproducing this tutorial: what you need")

st.markdown(
    """
    To follow the workflow on your own machine or cluster and then use this app for visualization,
    you will need the following components.
    """
)

st.markdown(
    """
    #### 1. GROMACS with MPI support

    You will use **GROMACS** for building the system, minimization, REMD, and trajectory handling.
    Install a recent version with MPI enabled using the official guide:

    <https://manual.gromacs.org/current/install-guide/index.html>
    """
)

st.markdown(
    """
    #### 2. PLUMED patched into GROMACS

    To perform metadynamics and handle collective variables (φ, ψ), compile **PLUMED**
    with your GROMACS build.

    Useful references:

    - PLUMED installation:
      <https://www.plumed.org/doc-v2.9/user-doc/html/_installation.html>  
    - Example notes for PLUMED + GROMACS:
      <https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-5.html>
    """
)

st.markdown(
    """
    #### 3. Alanine dipeptide starting structure

    This tutorial uses a capped Ace–Ala–Nme dipeptide (`dipep.pdb`) as the starting point.

    Click below to download exactly the structure used here:
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
    You can inspect this file in PyMOL, VMD, ChimeraX, or any molecular viewer to connect the
    geometry you see with the φ/ψ angles and the landscapes we will compute.
    """
)

st.markdown(
    """
    #### 4. What this app will do (next sections)

    Below this introduction, this app is intended to:

    - present annotated, copy-pasteable REMD and simulation setup blocks for alanine dipeptide,
    - you can compare your results for each step to the figures/screenshots here
    - perform the tutorial step by step in your own machine
    - let you upload (optionally) your own outputs (e.g. `COLVAR`, `fes_psi_phi.dat`, etc),
    - automatically generate:
      - FES(φ, ψ) heatmaps,
      - 1D profiles F(φ) and F(ψ),

    This page is your conceptual and practical starting point. From here, you can extend the app with:
    - a **Simulation Setup** tab (REMD + MetaD scripts),
    - an **Upload & Plot** tab (interactive analysis),
    - and comparison views (vacuum vs solvent, different force fields, etc.).
    """
)
