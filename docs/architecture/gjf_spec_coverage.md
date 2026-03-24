<!--
 * @Author: TMJ
 * @Date: 2026-03-23 23:07:40
 * @LastEditors: TMJ
 * @LastEditTime: 2026-03-24 15:01:53
 * @Description: 请填写简介
-->

# GJF Spec Coverage (Current Status)

This document compares current implementation against the Gaussian docs baseline:

- https://gaussian.com/input/
- https://gaussian.com/molspec/
- https://gaussian.com/zmat/
- https://gaussian.com/gic/
- https://gaussian.com/opt/
- https://gaussian.com/freq/
- https://gaussian.com/population/
- https://gaussian.com/geom/
- https://gaussian.com/capabilities/
- https://gaussian.com/basissets/

Implementation scope inspected:

- `src/molop/io/logic/qminput_frame_models/GJFFileFrame.py`
- `src/molop/io/logic/qminput_frame_parsers/GJFFileFrameParser.py`
- `src/molop/io/logic/gaussian_route_models.py`

## Coverage Summary

- Core input framing: **partial (improved)**
- Molecule specification: **partial (substantial coverage)**
- Z-matrix variants: **partial (substantial coverage)**
- GIC/ModRedundant modeling: **partial (improved)**
- Gaussian route semantics: **partial (substantial coverage)**

Estimated overall coverage for currently-targeted rules: **about 91-96%**.

## Detailed Matrix

| Spec Area    | Requirement                                               | Status      | Notes                                                                                                                                            |
| ------------ | --------------------------------------------------------- | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| input        | Link0 + Route + Title + Molecule + Additional split       | Partial     | Implemented in parser, including comment stripping and selected route-semantic exceptions; still heuristic blank-line flow                       |
| input        | Route as multi-line `#` section                           | Implemented | Supported by section slicing to terminating blank line                                                                                           |
| input        | Link1 strict semantics and pre-separator blank validation | Implemented | File parser enforces required blank line before `--Link1--`                                                                                      |
| input        | `!` comments stripping in parser                          | Implemented | Frame parser strips comment suffix/content before section detection                                                                              |
| input        | `@filename` include handling                              | Implemented | File parser expands relative includes, detects missing files and include cycles                                                                  |
| input        | Title constraints (required, <=5 lines)                   | Implemented | Title is required in normal flow, max 5 lines enforced                                                                                           |
| molspec      | Charge/multiplicity parsing                               | Implemented | Supports 2-field and multi-fragment forms                                                                                                        |
| molspec      | Fragment assignment by `Fragment=n`                       | Implemented | Multi-fragment inputs now require explicit contiguous `Fragment=n` assignments matching declared fragment charge/spin pairs                      |
| molspec      | Cartesian atom records                                    | Implemented | Basic parsing/rendering supported, with free-format delimiters tolerated                                                                         |
| molspec      | Free-format delimiters (space/tab/comma/slash)            | Implemented | Charge/multiplicity, atom lines, and z-matrix variable assignments tolerate multiple separators                                                  |
| molspec      | Z-matrix atom records                                     | Partial     | Internal-coordinate parsing/rendering exists, including working trailing `0/1` semantics; broader validation remains incomplete                  |
| molspec      | `Geom=Checkpoint/AllCheck` exceptions                     | Partial     | Parser supports missing title/molecule variants for selected route cases; broader route-driven omission semantics still incomplete               |
| zmat         | Alternate trailing format selector (`... 0/1`)            | Implemented | Accepted in internal-coordinate parsing                                                                                                          |
| zmat         | Preserve trailing `0/1` marker in rendering               | Implemented | Internal-coordinate renderer keeps marker token and emits Gaussian-safe numeric fields                                                           |
| zmat         | Render internal coordinates without unit suffixes         | Implemented | Internal-coordinate renderer outputs numeric tokens only (no `angstrom`/`degree`)                                                                |
| zmat         | Alternate `1` row geometry semantics                      | Implemented | Trailing `1` rows are interpreted as two-angle constructions during internal-to-cartesian conversion                                             |
| zmat         | Variables/Constants labeled blocks                        | Implemented | Variable assignments are parsed and applied to atom lines                                                                                        |
| zmat         | Unlabeled positional var/const inference                  | Implemented | Unlabeled assignments after atom block are accepted and applied                                                                                  |
| zmat         | Mixed Cartesian + internal in one block                   | Partial     | Some atom-level forms parse, full section semantics incomplete                                                                                   |
| zmat         | Line-number and label references                          | Partial     | Numeric references parse; symbolic atom-label references and validation remain incomplete                                                        |
| zmat         | Angle-domain and graph validity checks                    | Partial     | Positive bond-length, angle-domain, and basic backward-reference/distinct-reference validation now exist, but no complete graph validator exists |
| gic          | Detect/store GIC section as typed model                   | Implemented | `GJFGICSection` stores typed parsed lines and preserves raw section text                                                                         |
| gic          | Parse expression grammar and AST                          | Partial     | GIC lines now parse labels, options, function-like expressions, and nested top-level args, but not a full evaluable AST                          |
| gic          | Standalone/global options (`FreezeAll`, etc.)             | Partial     | Standalone options are structurally parsed into keyword/target/action/state fields, but semantic execution/state transitions are not implemented |
| gic          | Enforce non-mixing GIC vs ModRedundant                    | Implemented | Mixed section falls back to `unknown` and emits explicit diagnostic                                                                              |
| modredundant | Detect/store as typed section                             | Partial     | `GJFModRedundantSection` exists; syntax parsing is permissive                                                                                    |
| population   | Route-level Pop options (`full`, `nbo`, `nboRead`, `hirshfeld`, `cm5`, `mk`, `chelpg`, `orbitals`, radii options) | Implemented | `GaussianPopOptions` now provides strict typed extraction for common documented Pop suboptions                                                   |
| population   | Trailing NBO input block after molecule specification     | Implemented | `$NBO ... $END` is modeled as `GJFNBOSection` and rendered round-trip                                                                          |
| route        | Shared Gaussian route semantic model for both GJF and G16 logs | Implemented | `GaussianRouteSemantic` is now the single route semantic source for GJF and G16 log consumers                                                   |
| route        | Model chemistry parsing (`method/basis`, `method2//method1`) | Implemented | Supports single-layer and layered high//low model chemistry, including default opt+sp interpretation for double-slash routes                   |
| route        | Basis-set grammar summary (`family`, diffuse, polarization) | Implemented | Shared route model classifies major basis families and common diffuse/polarization markers                                                      |
| route        | Common job types (`sp`, `opt`, `freq`, `force`, `irc`, `scan`, `td`, `nmr`, `pop`, `stable`, `volume`, `admp`, `bomd`, `oniom`) | Implemented | Common route job families are explicitly recognized and added to capabilities                                                                   |
| route        | Common typed route option families (`Opt`, `Freq`, `TD`, `SCRF`, `Pop`, `Geom`) | Partial | High-frequency options are strict-typed; long-tail options still fall back to `option_maps` / `extra_options`                                 |
| route        | `Geom` keyword options from official geom page            | Implemented | Documented `Geom` options are represented in `GaussianGeomOptions`, including output-related options such as `PrintInputOrient`               |
| route        | `Freq` keyword options from official freq page            | Partial     | Many common options are typed (`anharmonic`, `raman`, `readisotopes`, `hpmodes`, etc.), but not every documented niche option is modeled     |
| route        | `Opt` keyword options from official opt page              | Partial     | Common optimization controls are typed (`ts`, `calcfc`, `restart`, `saddle`, `rcfc`, etc.), but some niche options remain fallback           |
| route        | `Pop` keyword options from official population page       | Partial     | Common population/NBO options are typed, but some specialized population-analysis variants remain fallback                                      |
| route        | Dieze tag parsing                                         | Implemented | `#`, `#N`, `#P`, `#T` are normalized only when present in the first route token; `# p` spacing is tolerated                                   |
| diagnostics  | Structured parse diagnostics by category/line span        | Partial     | Section-level diagnostics exist; fine-grained categories incomplete                                                                              |

## Prioritized Next Work

1. Add explicit molspec mode parser (`cartesian` / `zmatrix` / `mixed`) with stricter section termination and semantic validation.
2. Strengthen internal-coordinate validation (complete reference-graph validation and ambiguous alternate-row orientation diagnostics).
3. Implement full GIC evaluable AST and semantic option/state-transition handling.
4. Extend strict typed route parsing to more long-tail `Opt` / `Freq` / `Pop` / `TD` suboptions while keeping dict fallback for rarer cases.
5. Add rule-based diagnostics with line-level categories for all section parsers and route semantic extraction.
