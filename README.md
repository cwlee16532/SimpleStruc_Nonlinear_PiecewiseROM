# Rapid Postearthquake Fragility Assessment Using a Parametric Reduced-Order Model

MATLAB demonstration of damage-parameterized piecewise-linear seismic analysis using a proper orthogonal decomposition reduced-order model (POD-ROM).

The main script compares a trilinear reference model, a piecewise-linear full-order model (FOM), and a POD-ROM using a 100-degree-of-freedom (DOF) shear-building example.

---

## 📌 Related Publication

This repository accompanies the following paper:

**Chanwoo Lee**, Namsu Jeon, and Hyung-Jo Jung  
*Rapid Postearthquake Fragility Assessment of Operating Nuclear Structures Using a Parametric Reduced-Order Model*  
Earthquake Spectra, **42**, e70128 (2026)  
DOI: [10.1002/esp4.70128](https://doi.org/10.1002/esp4.70128)

> The published Version of Record is available at the journal website via the DOI link above.

> **Code scope:** The shared main script covers the illustrative 100-DOF example in Section 3.3. The nuclear containment FE model, monitoring-based model updating, multiple stripe analysis (MSA), and fragility-curve fitting described elsewhere in the paper are outside this script's scope.

---

## 🎯 Motivation

Postearthquake decisions require seismic assessments that reflect the current condition of an operating structure. However, repeatedly performing high-fidelity nonlinear analyses can make condition-specific fragility reassessment computationally demanding.

The paper addresses this challenge by combining:

- **Parametric reduced-order modeling** to accelerate structural response calculations.
- **Piecewise-linear dynamic analysis** to represent damage-dependent stiffness changes.
- **Condition-informed fragility reassessment** for operational early-damage evaluation within a slight-damage regime.

The shared example demonstrates the response-analysis component by approximating a trilinear structural model with a damage-parameterized piecewise-linear model and its POD-ROM.

---

## ⚙️ Algorithm Overview

The broader framework in the paper follows an offline-online workflow:

1. **Offline preparation** → Construct frequency-domain POD-based parametric ROMs and define damage-dependent stiffness regimes.
2. **Online model updating** → Estimate stiffness-related parameters from monitoring data.
3. **Accelerated seismic analysis** → Perform piecewise-linear response analyses across ground motions and intensity levels.
4. **Fragility reassessment** → Estimate condition-specific fragility curves from limit-state exceedance results.

### Workflow of the Shared Example

1. **Build the reference structure** → Assemble the mass and elastic stiffness matrices and determine Rayleigh damping coefficients.
2. **Analyze the trilinear model** → Apply a recorded ground motion at increasing target PGA levels using the Newmark-beta solver.
3. **Calibrate the piecewise-linear model** → Identify two effective stiffness reduction factors, `eta1` and `eta2`, by matching peak top-floor displacements at 0.35 g and 0.45 g.
4. **Construct the POD basis** → Call the frequency-response-based POD basis-building function.
5. **Compare FOM and ROM responses** → Evaluate intensity-displacement curves, displacement time histories, and stiffness-state progression.

---

## 📌 Case Study: 100-DOF Shear-Building Model

The example uses a one-dimensional shear-building model subjected to horizontal ground acceleration.

| Parameter | Default setting |
| --- | --- |
| Number of stories / DOFs | 100 |
| Lumped mass per story | 10,000 kg |
| Initial story stiffness | 2.5857 × 10¹⁰ N/m |
| Fundamental frequency | Approximately 4 Hz |
| Rayleigh damping | 5% in the first two modes |
| Trilinear displacement boundaries | 0.1 mm and 0.2 mm |
| Trilinear stiffness levels | `k1 = k_story`, `k2 = 0.5*k_story`, `k3 = 0.1*k_story` |
| Parameterized damage zone | Stories 1-20 |
| Top-displacement state thresholds | 6.5 mm and 9 mm |
| Calibration PGA levels | 0.35 g and 0.45 g |
| Stiffness-factor search grid | 500 candidates from 0.02 to 0.5 |
| Target PGA grid in the script | 0-0.50 g in 0.025 g increments |
| POD dimension used in the paper | 10 modes |

The default input is the horizontal component of the 2017 Pohang earthquake recorded at the PHA2 station. An alternative Gyeongju-record filename is also provided in the script.

The paper reports calibrated values of `eta1 = 0.436` and `eta2 = 0.153` for this example. The script recalibrates these parameters from the supplied record and solver implementation.

---

## 🚀 Running the Code

### Main Script and Required Supporting Functions

The supplied file is `main_100dof_trilinear.m`. It calls the following supporting functions, which must be available in the current folder or on the MATLAB path:

| File | Role |
| --- | --- |
| `main_100dof_trilinear.m` | Main analysis, calibration, and comparison workflow |
| `assembleK_linear.m` | Elastic stiffness-matrix assembly |
| `newmark_100dof_trilinear.m` | Trilinear reference-model time-history analysis |
| `newmark_100dof_piecewise_eta.m` | Piecewise-linear FOM time-history analysis |
| `build_POD_basis_piecewise_eta.m` | Frequency-response-based POD basis construction |
| `newmark_100dof_piecewise_eta_ROM.m` | Piecewise-linear ROM time-history analysis |

**The five supporting function definitions and the ground-motion CSV are separate dependencies; they are not embedded in the supplied main script.**

### Prepare the Ground-Motion Record

Provide a numeric CSV with a consistent time increment and the following column order:

| Column | Content | Expected unit |
| --- | --- | --- |
| 1 | Time | s |
| 2 | First horizontal acceleration component, H1 | m/s² |
| 3 onward, if present | Other numeric components; unused by this example | As applicable |

The script uses only the first two columns. It subtracts the first acceleration sample and then normalizes the record to each target PGA. The input must have a nonzero peak acceleration.

Select the input file in the main script:

```matlab
% eqFile = 'GY_USN_3axes.csv';  % Alternative Gyeongju record
eqFile = 'Po_PHA2_3axes.csv';   % Default Pohang record
```

Place the selected CSV in the MATLAB working folder or set `eqFile` to its full path.

### Execute

With the supporting functions and input record available, set the MATLAB current folder to the project folder and run:

```matlab
main_100dof_trilinear
```

The script performs the reference analyses, calibrates the piecewise-linear model, constructs the POD basis, and compares the three models.

---

## 📊 Response Comparison and Outputs

The example compares:

- **Baseline trilinear model** — Reference response with three story-stiffness regimes.
- **Piecewise-linear FOM** — Parameterized stiffness changes in the lower 20 stories.
- **POD-ROM** — Reduced-order simulation of the piecewise-linear model.

Generated figures include:

- Trilinear story force-displacement backbone.
- Original ground-acceleration time history.
- PGA versus maximum absolute top-floor displacement.
- Story-wise stiffness-regime histories for the trilinear model.
- Piecewise-linear system-state histories.
- Top-floor displacement time-history comparisons at 0.35 g and 0.45 g.
- Combined intensity-displacement curves for all three models.

---

## ⚖️ Copyright Notice

This repository contains original MATLAB implementation code developed by the author. The published journal article is not redistributed here. Please access the Version of Record through the DOI link above.
