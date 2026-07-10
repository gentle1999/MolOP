# Gaussian-Like Renderer

<!-- format-support:fakeg -->

| Item | Value |
| ---- | ----- |
| Format ID | `fakeg` |
| Extensions | `.fakeg` |
| Read | No |
| Write | Yes |
| Registry role | File writer |
| Data level | Parsed Gaussian output data |

`fakeg` renders parsed Gaussian output data as Gaussian-like text for human
inspection, compatibility testing, and downstream-tool handoff. It is not a
byte-for-byte Gaussian log reproducer.

`format_transform("fakeg")` follows the general transform default `frame=-1`;
pass `frame="all"` when full-file Gaussian-like rendering is required.

| Capability | Support | Rendered content | Boundaries |
| ---------- | ------- | ---------------- | ---------- |
| <!-- feature-area:File-level Gaussian-like writer -->File-level Gaussian-like writer | Partial | `.fakeg` files from parsed Gaussian output data. | File-level rendering only; frame-level fakeG writing is not provided, the selected frames follow `frame`, and byte-for-byte Gaussian log reproduction is not claimed. |
| <!-- feature-area:Structure and SCF energy rendering -->Structure and SCF energy rendering | Partial | Normalized orientation and SCF-cycle text from coordinate and energy fields. | Rendered sections are semantic Gaussian-like text, not original log text. |
| <!-- feature-area:Vibrational frequency rendering -->Vibrational frequency rendering | Partial | Frequency, reduced-mass, force-constant, IR-intensity, and per-mode displacement sections. | Current coverage is limited to frequency and IR-related fields; Raman/VCD sections are not advertised. |
| <!-- feature-area:Thermochemistry rendering -->Thermochemistry rendering | Partial | Temperature, correction, energy, entropy, heat-capacity, mass, inertia, and rotational/vibrational metadata. | This is a normalized thermochemistry summary, not a complete Gaussian thermochemistry pretty-printer. |
| <!-- feature-area:Reparseable frequency and thermochemistry output -->Reparseable frequency and thermochemistry output | Supported | Rendered frequency and thermochemistry content can be parsed back into frame fields. | The round trip proves supported semantic fields, not byte-for-byte equivalence with Gaussian output. |
