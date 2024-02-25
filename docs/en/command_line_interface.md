# Command line interface

MolOP also provides a simple command line interface for basic information extraction and format conversion. It is worth noting that, thanks to the chaining function calls provided by the [fire library](https://github.com/google/python-fire), the various commands of the MolOP CLI can be used in concert.

A complex example of using molop to read files and extract information. This script completed 5 tasks:

1. Read specific files and extract information with each last frame.
2. Transform the structures to GJF format.
3. Save the summary datasheet of the files.
4. Transform the structures to cdxml format.
5. Print out the smiles of the structures.


```python
! molop read "../../tests/test_files/g16log/*.log" --only_last_frame - gjf "../../tests/test_files/temp" --chk --template="../../tests/test_files/g16gjf/test.gjf" - summary "../../tests/test_files/temp" - chemdraw "../../tests/test_files/temp" - smiles - end
```

    MolOP parsing with 16 jobs: 100%|███████████████| 41/41 [00:01<00:00, 32.14it/s]
    0 files failed to parse, 41 successfully parsed
    gjf files saved to /home/tmj/proj/MolOP/tests/test_files/temp
    summary csv saved to /home/tmj/proj/MolOP/tests/test_files/temp/summary.csv
    chemdraw files saved to /home/tmj/proj/MolOP/tests/test_files/temp
    CC[C@H](C)[C@H](NC(C)=O)C(=O)NCC(=O)N[C@@H](CC([O])=O)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@H](C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CC([O])=O)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](CCCCNC(C)=O)C(=O)NC)[C@@H](C)CC.O.[O]C=O.[O]C=O.[Rh][Rh]
    C=C[C@H]1[C@@H]2[C@H](C[C@@]1(C)O)OC[C@@H]2C
    CCC[C@H]1CO[C@H]2C[C@@](C)(O)C[C@@H]12
    CC(C)(C)[C@@H]1COC(=[C-]C2[N-][C@H](C(C)(C)C)CO2)[N-]1.C[CH-]C(=O)[N-]c1ccccc1.[Cu@TB7+5].[O-]c1ccccc1
    [Br-].[Br-].[Ni@OH26+3]c1ccccc1.[S-]Cc1ccccc1.c1c[nH]cn1.c1c[nH]cn1
    CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F
    CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F
    C=[N+](C)[N-]C.COC(=O)[C@@]1(OC)C#CC(Br)(Br)CCCC1
    CNC(=O)C#[N+]/C(C(C)=O)=C(\C)[O-].COC(=O)C#CC(F)(F)F
    C/[NH+]=C(\C(=O)OC)c1ccccc1
    [CH3-]
    [C-]#CCC#C
    CC/C(C)=N/[O-]
    [C+]#CC#CC#C
    [C+]1=C[C@H]2O[C@H]2C1
    C1[C@@H]2[C]3CN2[C@@H]13
    [H]/N=c1/[cH-]onn1
    O=C1C=C=NC=[NH+]1
    Oc1nc[c]cn1
    Oc1n[c-]ncn1
    N/C=N/C=C(/[O-])O
    COC1=NCC[N-]1
    O=CN[CH][C@H]1CO1
    N#N.O=C=CC=O
    C[NH2+]C/C(C[O-])=N/[O-]
    Cc1[nH]c([O-])nc1O
    Nc1[nH]c([O-])nc1O
    O=CC1=C[N]C(O)=N1
    C#Cc1nc(O)[c]o1
    NC1(N)OCC(=O)O1
    N#C[C-](C#N)/[NH+]=C/N
    N[CH-]NC1=NC(=[OH+])N=N1
    [O-]c1onc2c1CNC2
    [O-]c1onc2c1C[NH2+]C2
    CN1[C+]=CC([N+](=O)[O-])=C1
    [NH2+]=c1[nH]ccc(F)c1[O-]
    [C-]#C[C@@H]1N2C#[N+]C[C@]12C
    C[C@H]1C2C[C-]=C3C=[N+]1[C@@H]32
    C[C@@H]1[C@H](C)C2=C(CC2)[C@@H]1C
    C/C([O-])=N/[N+]#N
    CN/C([O-])=N/[O+]=N\c1ccccc1


You can use MolOP CLI to do seprate tasks of course. There are the available commands:


```python
! molop
```

    [1mNAME[0m
        molop - CLI for MolOP.
    
    [1mSYNOPSIS[0m
        molop [4mCOMMAND[0m
    
    [1mDESCRIPTION[0m
        CLI for MolOP.
    
    [1mCOMMANDS[0m
        [1m[4mCOMMAND[0m[0m is one of the following:
    
         auto
           Auto process the current directory.
    
         chemdraw
           Save the cdxml file of specified frames of each file.
    
         end
           End the command chain, stop printing help comments.
    
         gjf
           Save the GJF file of any frame of each file.
    
         read
           Read the files given and set the file batch object.
    
         sdf
           Save the SDF file of all frames of each file.
    
         smiles
           Print the SMILES of last frame of each files.
    
         summary
           Save the summary csv file.
    
         xyz
           Save the XYZ file of all frames of each file.


Note that the `read` command is a mandatory command, and you can then chain function calls after it to continue with other commands. Any two commands need to be separated by the `-` symbol and terminated with the `end` command. For example:


```python
! molop read "../../tests/test_files/g16log/*.log" --only_last_frame - smiles - end
```

    MolOP parsing with 16 jobs: 100%|███████████████| 41/41 [00:01<00:00, 31.26it/s]
    0 files failed to parse, 41 successfully parsed
    CC[C@H](C)[C@H](NC(C)=O)C(=O)NCC(=O)N[C@@H](CC([O])=O)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@H](C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CC([O])=O)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](CCCCNC(C)=O)C(=O)NC)[C@@H](C)CC.O.[O]C=O.[O]C=O.[Rh][Rh]
    C=C[C@H]1[C@@H]2[C@H](C[C@@]1(C)O)OC[C@@H]2C
    CCC[C@H]1CO[C@H]2C[C@@](C)(O)C[C@@H]12
    CC(C)(C)[C@@H]1COC(=[C-]C2[N-][C@H](C(C)(C)C)CO2)[N-]1.C[CH-]C(=O)[N-]c1ccccc1.[Cu@TB7+5].[O-]c1ccccc1
    [Br-].[Br-].[Ni@OH26+3]c1ccccc1.[S-]Cc1ccccc1.c1c[nH]cn1.c1c[nH]cn1
    CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F
    CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F
    C=[N+](C)[N-]C.COC(=O)[C@@]1(OC)C#CC(Br)(Br)CCCC1
    CNC(=O)C#[N+]/C(C(C)=O)=C(\C)[O-].COC(=O)C#CC(F)(F)F
    C/[NH+]=C(\C(=O)OC)c1ccccc1
    [CH3-]
    [C-]#CCC#C
    CC/C(C)=N/[O-]
    [C+]#CC#CC#C
    [C+]1=C[C@H]2O[C@H]2C1
    C1[C@@H]2[C]3CN2[C@@H]13
    [H]/N=c1/[cH-]onn1
    O=C1C=C=NC=[NH+]1
    Oc1nc[c]cn1
    Oc1n[c-]ncn1
    N/C=N/C=C(/[O-])O
    COC1=NCC[N-]1
    O=CN[CH][C@H]1CO1
    N#N.O=C=CC=O
    C[NH2+]C/C(C[O-])=N/[O-]
    Cc1[nH]c([O-])nc1O
    Nc1[nH]c([O-])nc1O
    O=CC1=C[N]C(O)=N1
    C#Cc1nc(O)[c]o1
    NC1(N)OCC(=O)O1
    N#C[C-](C#N)/[NH+]=C/N
    N[CH-]NC1=NC(=[OH+])N=N1
    [O-]c1onc2c1CNC2
    [O-]c1onc2c1C[NH2+]C2
    CN1[C+]=CC([N+](=O)[O-])=C1
    [NH2+]=c1[nH]ccc(F)c1[O-]
    [C-]#C[C@@H]1N2C#[N+]C[C@]12C
    C[C@H]1C2C[C-]=C3C=[N+]1[C@@H]32
    C[C@@H]1[C@H](C)C2=C(CC2)[C@@H]1C
    C/C([O-])=N/[N+]#N
    CN/C([O-])=N/[O+]=N\c1ccccc1


All commands are controlled with additional parameters, you can use the `molop command --help` command to see the help.


```python
! molop read --help
```

    INFO: Showing help with the command 'molop read -- --help'.
    
    [1mNAME[0m
        molop read - Read the files given and set the file batch object.
    
    [1mSYNOPSIS[0m
        molop read [4mFILE_PATH[0m <flags>
    
    [1mDESCRIPTION[0m
        Read the files given and set the file batch object.
    
    [1mPOSITIONAL ARGUMENTS[0m
        [1m[4mFILE_PATH[0m[0m
            Type: str
            use regax to match files.
    
    [1mFLAGS[0m
        -c, --charge=[4mCHARGE[0m
            Type: Optional[]
            Default: None
            forced charge of the molecule, if not given, will use the charge written in the file or 0.
        -m, --multiplicity=[4mMULTIPLICITY[0m
            Type: Optional[]
            Default: None
            forced multiplicity of the molecule, if not given, will use the charge written in the file or 1.
        --only_extract_structure=[4mONLY_EXTRACT_STRUCTURE[0m
            Default: False
            if True, only extract the structure, else extract the whole file.
        --only_last_frame=[4mONLY_LAST_FRAME[0m
            Default: False
            if True, only extract the last frame, else extract all frames.
    
    [1mNOTES[0m
        You can also use flags syntax for POSITIONAL ARGUMENTS



```python
! molop gjf --help
```

    INFO: Showing help with the command 'molop gjf -- --help'.
    
    [1mNAME[0m
        molop gjf - Save the GJF file of any frame of each file.
    
    [1mSYNOPSIS[0m
        molop gjf <flags>
    
    [1mDESCRIPTION[0m
        Save the GJF file of any frame of each file.
    
    [1mFLAGS[0m
        --file_dir=[4mFILE_DIR[0m
            Type: Optional[str]
            Default: None
        --charge=[4mCHARGE[0m
            Type: Optional[int]
            Default: None
        -m, --multiplicity=[4mMULTIPLICITY[0m
            Type: Optional[int]
            Default: None
        -p, --prefix=[4mPREFIX[0m
            Type: str
            Default: '# g16 gjf \n'
        -s, --suffix=[4mSUFFIX[0m
            Type: str
            Default: '\n\n'
        -t, --template=[4mTEMPLATE[0m
            Type: Optional[str]
            Default: None
        --chk=[4mCHK[0m
            Type: bool
            Default: True
        -o, --oldchk=[4mOLDCHK[0m
            Type: bool
            Default: False
        --frameID=[4mFRAMEID[0m
            Type: int
            Default: -1


More useful features will be added in the future.
