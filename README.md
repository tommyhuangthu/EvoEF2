# INTRODUCTION to EvoEF2 and EvoEF
EvoEF (**Evo**Design physcical **E**nergy **F**unction) is the physical energy function component of EvoDesign, an approach developed for de novo protein sequence design given a protein backbone scaffold. EvoDesign uses a composite score function that combines the position-specific scoring matrix (PSSM, derived from structure-based evolutionary analysis) and EvoEF for sequence/structure decoy evaluation in the Monte Carlo-based sequence design simulations. In the earlier versions of EvoDesign, FoldX was used to compute the physical energy, which is, however, rather slow. To improve the scoring accuracy and speed, we developed EvoEF to replace FoldX. For more information, please refer to reference 2.

The motivation of developing EvoEF is for protein design. However, we find that EvoEF score performs poorly for producing native-like sequences on native protein backbones, suggesting EvoEF may not be appropriate for sequence design. This is partly because EvoEF was optimized by maximizing its ability to predict the experimental thermodynamic change data (ΔΔ<i>G</i>) rather than the native protein sequences. We reason that a score function specifically optimized for ΔΔ<i>G</i> prediction may not do well for sequence design. Thus, we extended EvoEF into EvoEF2 by introducing extra energy terms and reoptimizing energy weights. Computational experiments show that EvoEF2 significantly outperforms EvoEF in the native sequence recovery benchmark. For more information, please refer to reference 1.

Though we report EvoEF2 as an accurate energy function for protein sequence design in the paper, EvoEF2 is actually far more than an energy function for protein sequence design. It has multiple functions, e.g., protein energy computing, repairing incomplete protein side chains, building mutant model, optimizing the position of rotatable hydrogens, and most importantly, protein design.

# USAGE
EvoEF2 can do the following calculations given a protein structure in the PDB format. <b>Note that EvoEF2 works with amino acids/proteins only; it cannot deal with nucleotides, DNA, RNA, water, and/or other molecules</b>.

## 1. ComputeStability

**ComputeStability** computes the stability (or total energy) of a protein or protein complex.

$path/EvoEF2 --command=ComputeStability  --pdb=protein.pdb

## 2. ComputeBinding

**ComputeBinding** computes the binding interaction energy of a protein-protein complex.

$path/EvoEF2 --command=ComputeBinding --pdb=complex.pdb

$path/EvoEF2 --command=ComputeBinding --split=A,BC --pdb=complex.pdb

User should specify how to split chains into two parts for "ComputeBinding". Otherwise, EvoEF will compute the interaction energy between any pair of chains.

## 3. RepairStructure

**RepairStructure** repairs incomplete side chains of a protein. The side chains will be optimized to reduce steric clashes at the best. The hydroxyl hydrogens of Ser, Thr, and Tyr are optimized. Side-chain groups of His, Asn, and Gln may be flipped for optimizing hydrogen-bonding networks.

$path/EvoEF2 --command=RepairStructure --pdb=your.pdb --num_of_runs=3

Running this command will create a new structure model named "your_Repair.pdb" in the current directory. The option "--number_of_runs" specify the number of repeated times of repairing/optimizing the structure sequentially (default: 3).

## 4. BuildMutant

**BuildMutant** builds mutant model.

$path/EvoEF2 --command=BuildMutant --pdb=your.pdb --mutant_file=individual_list.txt  --num_of_runs=10 

Here, the "individual_list.txt" file shows the mutants that you want to build. It has the following format:

CA171A,DB180E;

Each mutant is written in one line ending with “;”, and multiple mutations in a mutant are divided by “,”. Note that there is no gap or space character between single mutations. For each single mutation, the first letter is the reference amino acid, the second letter is the chain identifier of the amino acid followed by its position in the chain, and the last letter is the amino acid after mutation. Running the command successfully should generate a new structure file named “your_Model_1.pdb”. The option "--number_of_runs" specify the number of repeated times of optimizing the rotamers of the mutated and surrounding residues sequentially (default: 10). 

## 5. OptimizeHydrogen

**OptimizeHydrogen** optimizes hydroxyl hydrogens for Ser, Thr, and Tyr.

$path/EvoEF2 --command=OptimizeHydrogen --pdb=your.pdb --num_of_runs=3

Running this command successfully will create a file named "your_OptH.pdb" in the current directory. The option "--number_of_runs" specify the number of repeated times of optimizing the locations of Ser/Thr/Tyr hydroxyl hydrogens (default: 3).

## 6. ProteinDesign

**ProteinDesign** can design de novo protein sequences on a fixed protein backbone.

$path/EvoEF2 --command=ProteinDesign --monomer --pdb=monomer.pdb

This command will design the whole sequence for a monomer protein ("monomer.pdb").

$path/EvoEF2 --command=ProteinDesign --ppint --design_chains=A --pdb=complex.pdb
  
This command will design the whole sequence for Chain A of a protein complex ("complex.pdb").

$path/EvoEF2 --command=ProteinDesign --ppint --design_chains=AB --pdb=ABCD.pdb

This command will design the whole sequence for Chains A and B of a four-chain protein complex ("ABCD.pdb").

Both backbone-depdent and backbone-indepdent rotamer libraries are supported by EvoEF2, users can specify the library that they want to use:

The available backbone-depdent rotamer libraries are:

<i>dun2010bb.lib

dun2010bb1per.lib (exclude the rotamers with probability < 1% from dun2010bb.lib)

dun2010bb3per.lib (exclude the rotamers with probability < 3% from dun2010bb.lib)</i>

The available backbone-indepdent rotamer libraries are:

<i>honig984.lib (984 rotamers for all 20 amino acids)

honig3222.lib (3222 rotamers for all 20 amino acids)

honig7421.lib (7421 rotamers for all 20 amino acids)

honig11810.lib (11810 rotamers for all 20 amino acids)</i>

By default, the dun2010bb3per.lib is used for a tradeoff of accuracy and running speed, users can choose a different backbone-dependent rotamer library, e.g.:

$path/EvoEF2 --command=ProteinDesign --monomer --rotlib=dun2010bb  --pdb=monomer.pdb

$path/EvoEF2 --command=ProteinDesign --monomer --bbdep=disable  --rotlib=honig984  --pdb=monomer.pdb

Note that you do not have to add the suffix ".lib" when you specify the rotamer library. However, we strongly suggest users take the default parameters if you are not familiar with rotamer library.


# INSTALLATION

First, download and unzip the EvoEF2 package.

## Windows

Users can try running the precompiled EvoEF2.exe first, which is built on my Windows 10. If it fails, you can then rebuild the EvoEF2.exe for your system. To do this, you need to install a g++ compilier. You can install the g++ compiler using MinGW. After installation, you can test if the g++ compiler has been successfully installed with "g++ -v" in the Windows CMD control console. If successful, it will show messages similar to:

Using built-in specs.<br>
COLLECT_GCC=g++<br>
COLLECT_LTO_WRAPPER=c:/mingw/bin/../libexec/gcc/mingw32/8.2.0/lto-wrapper.exe<br>
Target: mingw32<br>
Configured with: ../src/gcc-8.2.0/configure --build=x86_64-pc-linux-gnu <br>
--host=mingw32 --target=mingw32 --prefix=/mingw --disable-win32-registry <br>
--with-arch=i586 --with-tune=generic --enable-languages=c,c++,objc,obj-c++,fortran,ada <br>
--with-pkgversion='MinGW.org GCC-8<br>
Thread model: win32<br>
gcc version 8.2.0 (MinGW.org GCC-8.2.0-3)

Otherwise, it will say that g++ cannot be found. You should check your installation. If you try many times but still cannot install the program successfully, please contact me by the email below.

## UNIX/Linux/Mac
If you are working on a Linux-like system, things can be much easier, as the g++ compiler usually has already been installed. Go to the EvoEF2 main directory and run   "<b>./build.sh</b>" or "<b>g++ -O3 --fast-math -o EvoEF2 src/*.cpp</b>". If it fails, try without "--fast-math". If you are working on Mac, you may use "-fast-math".

# COPYRIGHT & CONTACT
Copyright (c) Xiaoqiang Huang. EvoEF2 is free to academic users. For suggestions, please contact xiaoqiah@umich.edu or xiaoqiah@outlook.com.

# REFERENCE
1. Huang X, Pearce R, Zhang Y. EvoEF2: accurate and fast energy function for computational protein design. Bioinformatics (2020), 36:1135-1142

2. Pearce R, Huang X, Setiawan D, Zhang Y. EvoDesign: Designing Protein–Protein Binding Interactions Using Evolutionary Interface Profiles in Conjunction with an Optimized Physical Energy Function. Journal of Molecular Biology (2019) 431: 2467-2476
