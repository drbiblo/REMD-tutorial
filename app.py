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

st.markdown("---")

# ================== STEP-BY-STEP: REMD SETUP ==================
st.subheader("Step-by-step: REMD setup in vacuum (Alanine dipeptide)")

st.markdown(
    """
    Below is the exact sequence we use to prepare a 4-replica REMD in vacuum
    with **AMBER99SB-ILDN** and PLUMED-ready inputs.
    """
)

# ---------- STEP 1 ----------
st.markdown("### 1. Build topology and coordinates from `dipep.pdb` (vacuum, AMBER99SB-ILDN)")

st.markdown(
    """
    Use `pdb2gmx` to generate the topology and starting coordinates.  
    When prompted:
    - choose **AMBER99SB-ILDN** (option 6 in this specific environment),
    - set **water model: none** (since we are in vacuum).
    """
)

st.code(
    "gmx_mpi pdb2gmx -f dipep.pdb -o conf.gro -p topol.top",
    language="bash",
)

st.image(
    "screenshots/screenshot1.png",
    caption="Example pdb2gmx selection: AMBER99SB-ILDN, no water.",
    use_container_width=True,
)

# ---------- STEP 2 ----------
st.markdown("### 2. Define a finite vacuum box")

st.markdown(
    """
    Even in vacuum, GROMACS requires a periodic box. We choose a cubic box with 2 nm padding.
    """
)

st.code(
    "gmx_mpi editconf -f conf.gro -o conf_box.gro -bt cubic -d 2",
    language="bash",
)

st.image(
    "screenshots/screenshot2.png",
    caption="Defining a cubic vacuum box around alanine dipeptide.",
    use_container_width=True,
)

# ---------- STEP 3 ----------
st.markdown("### 3. Create the minimization `.mdp` file")

st.markdown("Steepest descent minimization in vacuum with short-range cutoffs:")

st.code(
    r"""cat > min.mdp << 'EOF'
integrator       = steep
nsteps           = 50000
emtol            = 1000
emstep           = 0.01
cutoff-scheme    = Verlet
coulombtype      = Cut-off
rvdw             = 1.0
rcoulomb         = 1.0
rlist            = 1.0
nstlist          = 20
constraints      = none
EOF""",
    language="bash",
)

# ---------- STEP 4 ----------
st.markdown("### 4. Run energy minimization")

st.code(
    """gmx_mpi grompp -f min.mdp -p topol.top -c conf_box.gro -o min.tpr
gmx_mpi mdrun  -v -deffnm min""",
    language="bash",
)

st.image(
    "screenshots/screenshot4.png",
    caption="Minimization run: checking for convergence and stability.",
    use_container_width=True,
)

# ---------- STEP 5 ----------
st.markdown("### 5. Save minimized structure for REMD")

st.code("cp min.gro mini.gro", language="bash")

# ---------- STEP 6 ----------
st.markdown("### 6. Create an NVT template for REMD replicas")

st.markdown(
    """
    We now define a minimal NVT `.mdp` template for each replica.
    Here we set:
    - `dt = 0.002 ps`
    - `nsteps = 1000000` (2 ns; adjust for your needs)
    - no pressure coupling (vacuum),
    - no constraints (simple small system).
    """
)

st.code(
    r"""cat > nvt_template.mdp << 'EOF'
integrator            = md
dt                    = 0.002
nsteps                = 1000000
; outputs
nstxout-compressed    = 5000
nstenergy             = 1000
; cutoffs
cutoff-scheme         = Verlet
rvdw                  = 1.0
rcoulomb              = 1.0
rlist                 = 1.0
nstlist               = 20
; thermostat
tcoupl                = V-rescale
tc-grps               = System
tau_t                 = 0.1
ref_t                 = XXX   ; to be set per replica
; vacuum
pcoupl                = no
constraints           = none
EOF""",
    language="bash",
)

# ---------- STEP 7 ----------
st.markdown("### 7. Create per-replica folders and copy inputs")

st.code(
    """for i in 0 1 2 3; do
  mkdir -p $i
  cp topol.top $i/
  cp mini.gro  $i/
done""",
    language="bash",
)

# ---------- STEP 8 ----------
st.markdown("### 8. Create replica-specific `.mdp` files with distinct temperatures")

st.markdown(
    """
    Here we set a 4-replica ladder (example): 300, 366, 547, 996 K.  
    Adjust later if optimizing exchange probabilities.
    """
)

st.code(
    r"""cp nvt_template.mdp nvt0.mdp; sed -i 's/^ref_t.*/ref_t                 = 300/' nvt0.mdp
cp nvt_template.mdp nvt1.mdp; sed -i 's/^ref_t.*/ref_t                 = 366/' nvt1.mdp
cp nvt_template.mdp nvt2.mdp; sed -i 's/^ref_t.*/ref_t                 = 547/' nvt2.mdp
cp nvt_template.mdp nvt3.mdp; sed -i 's/^ref_t.*/ref_t                 = 996/' nvt3.mdp""",
    language="bash",
)

# ---------- STEP 9 ----------
st.markdown("### 9. Copy each `.mdp` into its replica folder")

st.code(
    """cp nvt0.mdp 0/nvt.mdp
cp nvt1.mdp 1/nvt.mdp
cp nvt2.mdp 2/nvt.mdp
cp nvt3.mdp 3/nvt.mdp""",
    language="bash",
)

# ---------- STEP 10 ----------
st.markdown("### 10. Build `.tpr` files for each replica")

st.code(
    """gmx_mpi grompp -f 0/nvt.mdp -p 0/topol.top -c 0/mini.gro -o 0/remd.tpr -maxwarn 10
gmx_mpi grompp -f 1/nvt.mdp -p 1/topol.top -c 1/mini.gro -o 1/remd.tpr -maxwarn 10
gmx_mpi grompp -f 2/nvt.mdp -p 2/topol.top -c 2/mini.gro -o 2/remd.tpr -maxwarn 10
gmx_mpi grompp -f 3/nvt.mdp -p 3/topol.top -c 3/mini.gro -o 3/remd.tpr -maxwarn 10""",
    language="bash",
)

st.image(
    "screenshots/screenshot10.png",
    caption="Successful generation of replica TPR files for REMD.",
    use_container_width=True,
)

# ---------- STEP 11 ----------
st.markdown("### 11. Add PLUMED input (φ/ψ) to each replica")

st.markdown(
    """
    Download the `plumed.dat` used in this tutorial:
    """
)

st.markdown(
    """
    <a href="https://raw.githubusercontent.com/drbiblo/FEMD-tutorial/main/plumed.dat"
       download="plumed.dat"
       style="
         display: inline-block;
         padding: 8px 14px;
         margin: 4px 0 10px 0;
         background-color: #2ca02c;
         color: white;
         text-decoration: none;
         border-radius: 6px;
         font-weight: 500;">
       ⬇ Download plumed.dat
    </a>
    """,
    unsafe_allow_html=True,
)

st.markdown("Copy it into each replica directory:")

st.code(
    r"""if [[ -f plumed.dat ]]; then
  for d in 0 1 2 3; do cp plumed.dat $d/; done
  echo "PLUMED file copied to all replicas."
else
  echo "WARNING: plumed.dat not found in this folder."
fi""",
    language="bash",
)

st.markdown(
    """
    At this point, your system is ready for a REMD + PLUMED run (to be detailed in the next step).
    """
)

st.markdown("---")
st.subheader("12. Launching the REMD Simulation")

st.markdown(
    """
With the four replica directories (`0/`, `1/`, `2/`, `3/`) prepared and `plumed.dat`
copied into each, we now run Replica Exchange Molecular Dynamics (REMD) in vacuum.

We use **4 replicas at 300, 366, 547, and 996 K**. For a tiny system like alanine dipeptide
in vacuum, a wide temperature ladder is intentional:

- Higher replicas (547 K, 996 K) can cross otherwise prohibitive φ/ψ barriers.
- Lower replicas (300 K, 366 K) preserve physically relevant sampling.
- Frequent exchanges between neighboring replicas allow each *conformation* to visit
  multiple temperatures, improving exploration of the Ramachandran landscape and
  enhancing convergence of 1D/2D free energy profiles.
"""
)

st.markdown("**Run the REMD job:**")

st.code(
    """export OMP_NUM_THREADS=2
mpirun -np 4 gmx_mpi mdrun -v \\
  -multidir 0 1 2 3 \\
  -deffnm remd \\
  -replex 100 \\
  -plumed plumed.dat \\
  -ntomp $OMP_NUM_THREADS \\
  -pin on""",
    language="bash",
)

st.markdown(
    """
### What this command is doing

**`export OMP_NUM_THREADS=2`**  
Sets 2 OpenMP threads per MPI rank. Together with `-ntomp $OMP_NUM_THREADS`, this gives each replica 2 threads.

**`mpirun -np 4`**  
Starts 4 MPI processes → one for each replica. This matches our 4 directories: `0/`, `1/`, `2/`, `3/`.

**`gmx_mpi mdrun`**  
Runs the MPI-enabled `mdrun`, required for multi-replica (`-multidir`) execution.

**`-multidir 0 1 2 3`**  
Tells GROMACS:
- simulation 0 uses `0/remd.tpr` (300 K),
- simulation 1 uses `1/remd.tpr` (366 K),
- simulation 2 uses `2/remd.tpr` (547 K),
- simulation 3 uses `3/remd.tpr` (996 K).
All run simultaneously as one composite REMD job.

**`-deffnm remd`**  
Sets the base name for outputs in each directory:
`0/remd.log`, `0/remd.xtc`, ..., `1/remd.log`, etc.

**`-replex 100`**  
Attempts exchanges between neighboring replicas every 100 MD steps.  
With `dt = 0.002 ps`, that’s every **0.2 ps**, giving frequent mixing for this small system.

**`-plumed plumed.dat`**  
Activates PLUMED in each replica using the local `plumed.dat` to monitor (and, if desired,
bias) φ/ψ. This is what will later let us reconstruct the free energy surfaces.

**`-ntomp $OMP_NUM_THREADS`**  
Ensures each MPI rank uses the 2 OpenMP threads we exported.

**`-pin on`**  
Pins threads to physical cores, improving stability and reproducibility of timings and sampling.

---
### How to read the console output (matches the screenshot)

When things are set up correctly, your terminal will show lines like:

- `This is simulation 0 out of 4 running as a composite GROMACS multi-simulation job.`
- `This is simulation 1 out of 4 ...`
- `This is simulation 2 out of 4 ...`
- `This is simulation 3 out of 4 ...`

This means:

- GROMACS has recognized **4 coupled simulations** (REMD replicas).
- Each is running with `Using 1 MPI process` and `Using 2 OpenMP threads`
  — exactly what we requested with `-np 4` and `OMP_NUM_THREADS=2`.
- For each replica you will see e.g.:
  `1000000 steps,   2000.0 ps.`
  which confirms: **1,000,000 steps × 0.002 ps = 2000 ps (2 ns)** per replica.

Near the end, lines like:

- `step 999900, remaining wall clock time: ...`

indicate that the replicas are approaching completion. In the next steps of the tutorial,
we will use the `remd.log` files to inspect exchange probabilities and then demultiplex
the trajectories for φ/ψ free energy analysis.
"""
)

st.image(
    "screenshots/screenshot11.png",
    caption="Example REMD output: four simulations (replicas) detected and propagated as a composite multi-simulation job.",
    use_container_width=True,
)
