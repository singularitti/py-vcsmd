# Scientific and architectural guide

This guide describes the implemented classical scope and the relationship between the scientific model, the numerical core, and application adapters. It uses $H$ for the cell matrix, with the direct cell vectors as its columns. Fractional particle coordinates are stored as rows of $S$, so Cartesian positions are

$$
R = S H^T.
$$

The implementation does not infer crystal symmetry and does not transform an evolving cell through a standardization operation.

## Scope and boundaries

The solver covers classical Lennard–Jones particles under hydrostatic external pressure. It supports fixed-cell dynamics and minimization, Parrinello–Rahman cell dynamics and minimization, modified-metric dynamics and minimization, and reference-cell strain dynamics and minimization. The active interaction is an unshifted, inclusive-cutoff Lennard–Jones potential with one global $\sigma$ and $\epsilon$ for all particles. The original argon values are available as a named preset.

Tabulated potentials, electronic-structure or DFT calls, tensor external stress, crystal-symmetry detection, standardization, disabled experiments, and historical output writers are outside the implementation boundary. The compatibility adapter reads old data; it does not import historical names or bugs into the numerical package.

## Units and data models

Pint quantities are required at the configuration boundary. Normalization converts all numerical work to float64 in a coherent atomic system:

| Quantity | Internal unit |
| --- | --- |
| length and cell vectors | bohr |
| energy | Rydberg |
| time | `rydberg_time = ℏ / Ry` |
| ordinary mass and cell inertia | `rydberg_mass` |
| modified-metric or reference-strain cell inertia | `rydberg_mass / bohr**4` |
| hydrostatic pressure | `Ry / bohr**3` |

Here `rydberg_mass = Ry * rydberg_time**2 / bohr**2`, equal to two electron masses. The native unit registry converts its own quantities directly; quantities from another Pint registry pass through SI base units for interchange. Keeping native conversions direct prevents repeated configuration reads from accumulating rounding changes.

The public layer returns Pint quantities with dimensional names. Fractional-coordinate rates have units $1/t$; cell rates have units length/time; peculiar Cartesian particle velocities are $H\dot S$. The full instantaneous Cartesian derivative of a particle position also includes the deformation term $\dot Hs_i$:

$$
\dot r_i = H\dot s_i + \dot Hs_i.
$$

The extended-system kinetic stress and atomic kinetic energy use the peculiar velocity $H\dot s_i$. The public state exposes both peculiar and full Cartesian velocities explicitly.

The main models are `Structure`, `SimulationConfig`, normalized `NumericalModel`, `InitialConditions`, `SimulationState`, `Observables`, and `StepResult`. Configuration models retain Pint quantities and enum choices; file adapters serialize the enum values. Normalized models contain resolved data and no file paths. Frozen dataclasses own read-only arrays. A step returns a new state and does not modify caller-owned arrays.

## Geometry and periodic interactions

Let

$$
V=\det H,
\qquad
g=H^T H,
\qquad
C=V H^{-T}.
$$

`C` is the cofactor matrix and is the source’s reciprocal-cell cross-product matrix. Its columns are the cross products of the other two direct cell vectors. The inverse metric is

$$
g^{-1}=\frac{C^T C}{V^2}=H^{-1}H^{-T}.
$$

For an image shift $n\in\mathbb Z^3$, a fractional displacement is $s_j-s_i+n$, and the Cartesian displacement is $(s_j-s_i+n)H^T$ in the row-vector convention. The central cell counts each distinct pair once. Every nonzero image shift counts ordered particle pairs, including self-image interactions; energy, virial, and cell-virial terms receive half weight for those ordered image terms, while the central particle receives the full force. The cutoff condition is $r^2\le r_c^2$, with no energy shift.

For the active potential,

$$
U(r)=4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-
\left(\frac{\sigma}{r}\right)^6\right],
$$

and the force on particle $j$, with $d=r_j-r_i$, is

$$
f_{ji}=\frac{24\epsilon}{r^2}
\left[2\left(\frac{\sigma^2}{r^2}\right)^6-
\left(\frac{\sigma^2}{r^2}\right)^3\right]d.
$$

The periodic implementation evaluates images in bounded vectorized batches. It rejects zero-distance collisions and validates that the requested image range is sufficient for the cutoff.

## Particle equations

Forces are first evaluated in Cartesian coordinates. In row notation, the fractional force acceleration is

$$
\ddot S_{\rm force}=M^{-1}F H^{-T},
$$

where $M^{-1}$ scales each particle row by its inverse mass. In variable-cell modes, the metric coupling is

$$
\dot g=\dot H^T H+H^T\dot H,
\qquad
\ddot S=\ddot S_{\rm force}-\dot S(g^{-1}\dot g)^T.
$$

Fixed-cell modes use only the force term. The implementation stores fractional rates and reconstructs peculiar Cartesian velocities as $\dot S H^T$.

## Cell equations

The kinetic and virial stress entering the cell force is

$$
\Pi=\frac{1}{V}\left(\sum_i m_i v_i v_i^T+\mathcal V\right),
\qquad
A_0=\frac{1}{W}(\Pi-P I)C,
$$

where $v_i=H\dot s_i$, $\mathcal V$ is the configurational virial matrix, $P$ is hydrostatic external pressure, and $W$ is the cell inertia. Ordinary cell dynamics uses $W$ with mass units. Modified-metric and reference-strain formulations use $W$ with mass/length⁴ units and multiply by an inverse metric.

For modified-metric modes, define

$$
F=C^TC,
\qquad F^{-1}=\frac{H^T H}{V^2},
\qquad E=\dot H^T\dot H,
$$

$$
\dot C=\operatorname{tr}(H^{-1}\dot H)C-C\dot H^T H^{-T},
\qquad
\dot F=\dot C^T C+C^T\dot C.
$$

The cell acceleration is

$$
\ddot H=\left[A_0+G_{\rm geom}-\dot H\dot F\right]F^{-1},
$$

with

$$
G_{\rm geom}=\operatorname{tr}(EF)H^{-T}-H^{-T}EF.
$$

The absence of a leading one-half is deliberate. Differentiating $F=C^TC$ produces two equal contractions against the symmetric $E$, which cancel the one-half in the kinetic-metric derivative.

For reference-cell strain modes, with $H_0$, $V_0$, and $C_0=V_0H_0^{-T}$ fixed at initialization,

$$
\ddot H=A_0\frac{H_0^T H_0}{V_0^2}.
$$

After the formulation-specific acceleration, every variable-cell mode applies the source strain projection:

$$
D=\ddot H H_0^{-1},
\qquad
D\leftarrow\frac12(D+D^T),
\qquad
\ddot H\leftarrow D H_0.
$$

This is a continuous tensor projection. It is not crystal-symmetry detection.

The cofactor derivative used for the modified metric can also be expressed without an explicitly stored four-index tensor. For a basis perturbation $E^{ij}$,

$$
\frac{\partial C}{\partial H_{ij}}
=\operatorname{tr}(H^{-1}E^{ij})C-C(E^{ij})^TH^{-T},
$$

followed by

$$
\frac{\partial F}{\partial H_{ij}}
=\left(\frac{\partial C}{\partial H_{ij}}\right)^TC
+C^T\frac{\partial C}{\partial H_{ij}}.
$$

This is the analytical form of the source’s cofactor derivative table.

## Beeman integration and state history

For any generalized coordinate $x$, with acceleration $a$, the predictor is

$$
x_{n+1}^{p}=x_n+h\dot x_n+
\frac{h^2}{6}(4a_n-a_{n-1}).
$$

After forces and cell equations are evaluated at the predicted state, the corrector is

$$
\dot x_{n+1}=\frac{x_{n+1}^{p}-x_n}{h}+
\frac{h}{6}(2a_{n+1}+a_n)
=\dot x_n+\frac{h}{6}(2a_{n+1}+5a_n-a_{n-1}).
$$

The state stores both current and previous particle and cell accelerations. A fresh initialization sets the previous acceleration equal to the current acceleration, which gives a defined second-order self-start. Native checkpoints preserve these histories for exact continuation.

Fractional positions are wrapped into $[0,1)$ after prediction. Floating-point remainders that round to exactly 1 are mapped to 0 so that checkpoints retain the same half-open interval. Cell vectors are never standardized or rotated by the numerical core.

## Energies, pressure, and temperature

Atomic kinetic energy is

$$
K_a=\frac12\sum_i m_i\lVert H\dot s_i\rVert^2.
$$

The cell kinetic energies are

$$
K_c=\frac W2\lVert\dot H\rVert_F^2
\quad\text{(ordinary cell modes)},
$$

$$
K_c=\frac W2\operatorname{tr}(\dot H F\dot H^T)
=\frac W2\lVert\dot H C^T\rVert_F^2
\quad\text{(modified metric)},
$$

and the same expression with $F_0=C_0^TC_0$ for reference strain. The orientation $\lVert\dot H C^T\rVert_F^2$ is required by the metric in the equations; the historical expression $\lVert C^T\dot H\rVert_F^2$ is not equivalent in general.

For variable-cell modes, the external-pressure contribution is $U_P=PV$, so $U=U_{\rm LJ}+PV$. Fixed-cell modes report the Lennard–Jones energy without an external-volume term. The internal hydrostatic pressure is

$$
P_{\rm int}=\frac{2K_a+\operatorname{tr}(\mathcal V)}{3V}.
$$

Atomic temperature uses $3(N-1)$ degrees of freedom after global center-of-mass removal. The historical extended-system controller uses six cell degrees of freedom and therefore $3(N+1)$ total degrees of freedom in variable-cell modes. These are exposed separately as atomic, instantaneous extended, and running extended temperatures. Zero-temperature and single-atom cases are handled explicitly.

At a controller interval, the running average kinetic energy is compared with the target temperature. If the tolerance is exceeded, particle rates are rescaled and dynamic cell rates are rescaled for the dynamic formulations. Minimization modes quench velocity components whose current and previous accelerations have opposite signs. `StepResult.observables` describes the step before these events; `StepResult.state` contains the state after them.

## Configuration, execution, and checkpoints

`vcsmd.config` defines the domain configuration and `prepare` performs quantity conversion and formulation-specific inertia validation. `vcsmd.io.config` is the only native text-format adapter; it accepts JSON, YAML, and TOML through one schema. `vcsmd.execution` creates documented run folders, streams CSV data, records events and provenance, and saves checkpoints. The computational modules do not import these adapters.

The public checkpoint functions are:

```python
from pathlib import Path

from vcsmd.io import load_checkpoint, save_checkpoint

save_checkpoint(Path("checkpoint.npz"), model, state)
model, state = load_checkpoint(Path("checkpoint.npz"))
```

The archive uses named numeric NumPy arrays and a versioned JSON metadata scalar with `allow_pickle=False`. It contains all model parameters, positions, rates, current and previous accelerations, reference geometry, time and step counters, running accumulators, and remaining controller rescalings. The public adapter wraps the raw numerical state in quantity-bearing fields. Exact continuation is verified within the same numerical software environment; different libraries or hardware may change floating-point rounding.

## Compatibility conversion

`vcsmd.compat` is an isolated, one-way reader for historical run directories and inputs. The [compatibility guide](legacy-conversion.md) documents its historical filenames, codes, units, and imported schemas. Dependencies point from the adapter toward native configuration validation; no computational module imports compatibility code.

Conversion combines available energy, pressure, temperature, and cell-parameter records by step, preserves source values and provenance, and writes CSV import datasets plus a conversion report. These datasets describe available source fields and are not interchangeable with complete native trajectory or cell-matrix streams. Header-only historical trajectory files are recorded as having no trajectory. Lengths and angles alone are not turned into a cell-vector history. A damaged or incomplete restart remains explicitly partial: missing reference geometry, integrator history, and damaged velocity columns prevent it from being advertised as an exact native checkpoint. Unused historical potential files remain outside the native computational configuration.

## Corrected historical behavior and known differences

The Python implementation deliberately corrects several demonstrable defects in the source implementation: fresh Beeman histories no longer start from invented zero acceleration; corrected velocity reconstruction does not add a second copy of the reconstructed peculiar velocity; modified-metric and strain cell kinetic energies use the metric orientation in their equations; restart histories are preserved in native checkpoints; single-atom and zero-temperature paths are explicit; and velocity initialization uses a local NumPy generator with the application default seed 119. The public API distinguishes peculiar velocity from cell-deformation motion.

The original source and accompanying project notes are useful provenance but are secondary to the equations implemented here. Historical code identifiers and filenames are confined to this compatibility and discrepancy discussion so that the native API remains organized around scientific concepts.

The implementation was reconciled against the source and the 2017 notes, rather than treating either as an infallible specification:

| Topic | Evidence in the historical material | Resolution |
| --- | --- | --- |
| Applied pressure | The input notes describe scalar pressure. The tensor equation in `invariant-md.tex`, lines 181–209, writes the pressure subtraction without an explicit identity matrix; `celq.f`, lines 1088–1106, subtracts the scalar only on the diagonal. | Interpret external pressure as $P I$. The internal virial remains a tensor; external tensor stress is outside this interface. |
| Cell kinetic energy | `invariant-md.tex`, lines 257–284, uses $\operatorname{tr}(\dot H F\dot H^T)$. `init.tex`, lines 44–52, and the source at lines 1161–1192 use the opposite matrix orientation. | Follow the kinetic metric in the Lagrangian. The two expressions agree in special geometries, which can conceal this defect. |
| Metric derivative | The source at lines 2694–2738 includes both cofactor-derivative terms and a factor of one-half. `sigsp.tex`, line 36, labels the derivative contraction as an inverse metric. | Combine the two equal contractions to obtain the closed geometric term above. The actual inverse metric is $H^T H/V^2$. |
| Strain projection | `init.tex`, lines 26–42, derives the projection; the source at lines 1111–1145 applies it after every variable-cell branch. | Preserve the projection across all variable-cell modes. No crystal-symmetry operations are introduced. |
| Fresh and restarted histories | The source at lines 1047–1054 zeros previous particle acceleration even after reading a restart. Its restart writer at lines 317–326 repeats cell-state columns. | Define fresh Beeman histories and store complete native continuation state. Historical restart data remains partial. |
| Particle velocities | The source at lines 1302–1361 removes a per-species arithmetic mean and rescales. At lines 2259–2293 it adds another reconstructed Cartesian velocity during correction. | Use a global mass-weighted center-of-mass constraint and reconstruct the corrected peculiar velocity once. The latter defect concerns the immediate stored velocity; the next source force pass reconstructs it again. |
| Paper scope | The attached 1993 paper combines strain dynamics with converged electronic energies, forces, and stresses. Later notes also discuss electronic-structure packages. | Reproduce the active classical Lennard–Jones solver. Electronic degrees of freedom and the paper's material-specific calculations are not part of this rewrite. |

These source line numbers refer to the supplied historical files. The numerical checks and preserved implementation snapshots in the [validation report](validation.md) identify the executed Python implementation.

## References and validation boundaries

The variable-cell equations are grounded in the formal literature:

1. M. Parrinello and A. Rahman, “Crystal Structure and Pair Potentials: A Molecular-Dynamics Study,” *Physical Review Letters* 45, 1196 (1980), DOI [10.1103/PhysRevLett.45.1196](https://doi.org/10.1103/PhysRevLett.45.1196).
2. R. M. Wentzcovitch, “Invariant molecular-dynamics approach to structural phase transformations,” *Physical Review B* 44, 2358 (1991), DOI [10.1103/PhysRevB.44.2358](https://doi.org/10.1103/PhysRevB.44.2358).
3. D. Beeman, “Some multistep methods for use in molecular dynamics calculations,” *Journal of Computational Physics* 20, 130–139 (1976), DOI [10.1016/0021-9991(76)90059-0](https://doi.org/10.1016/0021-9991(76)90059-0).
4. J. E. Jones, “On the determination of molecular fields. II. From the equation of state of a gas,” *Proceedings of the Royal Society A* 106, 463–477 (1924), DOI [10.1098/rspa.1924.0082](https://doi.org/10.1098/rspa.1924.0082).
5. R. M. Wentzcovitch, J. L. Martins, and G. D. Price, “Ab Initio Molecular Dynamics with Variable Cell Shape: Application to MgSiO₃,” *Physical Review Letters* 70, 3947–3950 (1993), DOI [10.1103/PhysRevLett.70.3947](https://doi.org/10.1103/PhysRevLett.70.3947). This is the attached paper; its electronic-structure component is broader than the active classical source.

The project’s `MSAE-E6237` notes and the historical source are retained as secondary implementation documentation. Full-run validation claims belong in the executed validation report. In particular, fixed-cell, modified-metric dynamics, and modified-metric minimization examples do not by themselves establish full-trajectory validation for the other formulations, and direct Fortran parity requires an executed reference run.
