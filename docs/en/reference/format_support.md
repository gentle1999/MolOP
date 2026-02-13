# Supported Formats

MolOP supports a variety of chemical file formats for reading and writing. The following table summarizes the built-in support:

| Format ID | Typical Extensions                     | Read | Write | Notes                          |
| --------- | -------------------------------------- | ---- | ----- | ------------------------------ |
| `g16log`  | `.log`, `.out`, `.g16`, `.gal`, `.irc` | Yes  | No    | Gaussian 16 output files.      |
| `xyz`     | `.xyz`                                 | Yes  | Yes   | Standard XYZ coordinate files. |
| `sdf`     | `.sdf`, `.sd`, `.mol`                  | Yes  | Yes   | MDL Structure-Data File.       |
| `smi`     | `.smi`, `.txt`                         | Yes  | Yes   | SMILES strings.                |
| `gjf`     | `.gjf`, `.gif`, `.com`, `.gau`, `.gjc` | Yes  | Yes   | Gaussian input files.          |
| `orcainp` | `.inp`                                 | Yes  | Yes   | ORCA input files.              |
| `cml`     | `.cml`                                 | No   | Yes   | Chemical Markup Language.      |

!!! note
This format support checklist will be updated over time as new formats are added.
