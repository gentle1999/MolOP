# 命令行接口

MolOP 还提供了一个简单的命令行接口，用于基本信息提取和格式转换。值得注意的是，得益于[fire library](https://github.com/google/python-fire)提供的链式函数调用，MolOP CLI 的各种命令可以协同使用。

使用 molop 读取文件并提取信息的复杂示例。该脚本完成了 5 项任务：

1. 读取特定文件并提取最后一帧的信息。
2. 将结构转换为 GJF 格式。
3. 保存文件的摘要数据表。
4. 将结构转换为 cdxml 格式。
5. 打印出结构的smiles字符串。


```python
! molop read "../../tests/test_files/g16log/*.log" --only_last_frame - gjf "../../tests/test_files/temp" --template="../../tests/test_files/g16gjf/test.gjf" - summary "../../tests/test_files/temp" - chemdraw "../../tests/test_files/temp" - smiles - end
```

    MolOP parsing with 32 jobs: 100%|███████████████| 39/39 [00:01<00:00, 20.70it/s]
    0 files failed to parse, 39 successfully parsed
    gjf files saved to /home/cathayana/py_project/HONG/MolOP/tests/test_files/temp
    summary csv saved to /home/cathayana/py_project/HONG/MolOP/tests/test_files/temp/summary.csv
    chemdraw files saved to /home/cathayana/py_project/HONG/MolOP/tests/test_files/temp
    C=C[C@H]1[C@@H]2[C@H](C[C@@]1(C)O)OC[C@@H]2C
    CCC[C@H]1CO[C@H]2C[C@@](C)(O)C[C@@H]12
    CC(C)(C)[C@@H]1COC(=[C-]C2[N-][C@H](C(C)(C)C)CO2)[N-]1.C[CH-]C(=O)[N-]c1ccccc1.[Cu@TB7+5].[O-]c1ccccc1
    [Br-].[Br-].[Ni@OH26+3]c1ccccc1.[S-]Cc1ccccc1.c1c[nH]cn1.c1c[nH]cn1
    CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F
    CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F
    C=[N+](C)[N-]C.COC(=O)[C@@]1(OC)C#CC(Br)(Br)CCCC1
    CNC(=O)C#[N+]/C(C(C)=O)=C(\C)[O-].COC(=O)C#CC(F)(F)F
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


当然，您也可以使用 MolOP CLI 分别执行不同的任务。您可以使用以下命令：


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


请注意，`read` 命令是一条必须执行的命令，您可以在它之后通过函数调用链继续执行其他命令。任何两条命令都需要用 `-` 符号分隔，并以 `end` 命令结束。例如：


```python
! molop read "../../tests/test_files/g16log/*.log" --only_last_frame - smiles - end
```

    MolOP parsing with 32 jobs: 100%|███████████████| 39/39 [00:01<00:00, 23.60it/s]
    0 files failed to parse, 39 successfully parsed
    C=C[C@H]1[C@@H]2[C@H](C[C@@]1(C)O)OC[C@@H]2C
    CCC[C@H]1CO[C@H]2C[C@@](C)(O)C[C@@H]12
    CC(C)(C)[C@@H]1COC(=[C-]C2[N-][C@H](C(C)(C)C)CO2)[N-]1.C[CH-]C(=O)[N-]c1ccccc1.[Cu@TB7+5].[O-]c1ccccc1
    [Br-].[Br-].[Ni@OH26+3]c1ccccc1.[S-]Cc1ccccc1.c1c[nH]cn1.c1c[nH]cn1
    CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F
    CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F
    C=[N+](C)[N-]C.COC(=O)[C@@]1(OC)C#CC(Br)(Br)CCCC1
    CNC(=O)C#[N+]/C(C(C)=O)=C(\C)[O-].COC(=O)C#CC(F)(F)F
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


所有命令都由附加参数控制，您可以使用 `molop command --help` 命令查看帮助。


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


More useful features will be added in the future.
