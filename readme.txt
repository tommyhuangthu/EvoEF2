                           README for EvoEF2
------------------------------------------------------------------------------


What's EvoEF?
--------------
EvoEF is the abbreviation of EvoDesign physical Energy Function, where EvoDesign is 
an evolution-based approach to de novo protein sequence design. EvoDesign combines 
an evolution- and physics-based scoring function for protein design simulation. In 
the original EvoDesign, FoldX is used as the physics-based score. We found that FoldX 
is very slow and is not well optimized for modeling protein-protein interactions. 
To improve the computational accuracy and speed of the physical score for EvoDesign, 
we developed EvoEF to replace FoldX. Our benchmark test showed that EvoEF is more 
accurate than FoldX and about three to five times faster than FoldX when they were 
benchmarked in the same way to predict thermodynamic change upon amino acid mutations 
(see reference 2 for more details).


What's EvoEF2?
----------------
EvoEF2 is an updated version of EvoEF. The motivation of developing EvoEF is for 
protein design. However, we found that EvoEF alone without evolutionary profile 
showed a poor performance to produce native-like sequence, suggesting the inability 
of EvoEF for de novo protein sequence design. We found the reason was that we 
optimized EvoEF by correlating with the experimental thermodynamic change data 
(i.e. ddG), and a potential specifically optimized for ddG prediction cannot do well 
for de novo protein design. We therefore extended EvoEF into EvoEF2 by introducing 
a few new energy terms and reoptimizing all the energy weights. EvoEF2 yielded a much 
better performance on sequence design than EvoEF. Although here we regard EvoEF2 as 
an energy function, however, it can work as a framework for different tasks of 
protein modeling, such as computing energy, repairing incomplete protein side chains, 
building mutant model, optimizing the position of rotatable hydrogens, and most 
importantly protein design.


What EvoEF2 can do?
---------------------
The following useful functions are supported in EvoEF:

o  ComputeStability
   -- compute the stability of a given protein molecule in PDB format.

o  ComputeBinding
   -- compute the binding affinity of protein-protein dimeric complexes.

o  RepairStructure
   -- repair incomplete side chains of user-provided model and minimize 
      energy of give model to reduce possible steric clashes.
  
o  BuildMutant
   -- build mutant models, which are useful for ddG assessment in combination
      with the wild-type. Before running the BuildMutant module, we suggest 
      users run the RepairStructure module to minimize the structure using 
      the EvoEF2 force field
  
o  OptimizeHydrogen
   -- optimize the positions of hydrogens in the hydroxyl groups (e.g. Ser, 
      Thr, and Tyr).

o  ProteinDesign
   -- de novo protein sequence design using a fixed backbone.
  
o  and some other commands (please see the source code for details)


Installation
------------
First of all, download the EvoEF2 package.

o If you are using a Windows system, you can directly run the 
  executable EvoEF2.exe program.

o If you are working in a Linux system, go to the EvoEF2 main directory and run:
  ./build.sh
  or run:
  g++ -O3 --fast-math -o EvoEF2 src/*.cpp
  to build the EvoEF2 program.
  If it does not work, try without the '--fast-math' option.

o If you want to build a new EvoEF2.exe in Windows, you need to install a g++ 
  compilier first. You can make the installation with MinGW. After you 
  complete the installation of g++ compiler and set the environmental 
  variable. You can test if the g++ compiler has been successfully 
  installed with 'g++ -v' in the cmd.exe control console. If successful, 
  it will show messages similar to:

  "Using built-in specs.
  COLLECT_GCC=g++
  COLLECT_LTO_WRAPPER=c:/mingw/bin/../libexec/gcc/mingw32/8.2.0/lto-wrapper.exe
  Target: mingw32
  Configured with: ../src/gcc-8.2.0/configure --build=x86_64-pc-linux-gnu \
  --host=mingw32 --target=mingw32 --prefix=/mingw --disable-win32-registry \
  --with-arch=i586 --with-tune=generic --enable-languages=c,c++,objc,obj-c++,fortran,ada \
  --with-pkgversion='MinGW.org GCC-8
  Thread model: win32
  gcc version 8.2.0 (MinGW.org GCC-8.2.0-3)"
   
  Otherwise, it will say that g++ cannot be found. You should check your installation.
  If you try many times but still cannot install the program successfully, please 
  contact report problems to the emails listed below.


Usage
-----
o To compute protein stability, you can run:

  ./EvoEF2 --command=ComputeStability  --pdb=protein.pdb


o To compute protein-protein binding energy of a dimer complex, you can run:

  ./EvoEF2 --command=ComputeBinding --pdb=complex.pdb

  Sometimes, you may have a multi-chain complex structure, which can still be 
  handled by EvoEF, but users need to specify how to split the chains for binding calculation.
  For example, say you have a four-chain (A,B,C, and D) complex, you can split the chains by:
    
  ./EvoEF2 --command=ComputeBinding --pdb=multi_chain_complex.pdb --split_chains=AB,CD
    
  which calculates the binding between partner 'AB' and 'CD'.
    
  or:
    
  ./EvoEF2 --command=ComputeBinding --pdb=multi_chain_complex.pdb --split_chains=AC,BD
    
  which calculates the binding between partner 'AC' and 'BD'
    
  In one word, users can split chains in their own ways.

o To repair the structure model and do energy minimization:

  ./EvoEF2 --command=RepairStructure --pdb=xxxx.pdb

  A new structure model name "xxxx_Repair.pdb" will be built in the directory 
  where you run the command.

o To build mutation model, you can run:

  ./EvoEF2 --command=BuildMutant --pdb=xxxx.pdb --mutant_file=individual_list.txt

  where the "individual_list.txt" file shows the mutants that you want to build. 
  It has the following format:
  
  CA171A;
  CA171A,DB180E;

  Each mutation is written in one line ending with “;”, and multiple mutants are 
  divided by “,”. Note that there’s no gap/space between single mutations. For 
  each single mutation, the first alphabet is the wild-type amino acid, the second 
  is the identifier of the chain that the amino acid is attached to, the number is 
  the position of the amino acid in the chain, and the last alphabet is the amino 
  acid after mutation. Running the command successfully should generate mutant models 
  namded as “xxxx_Model_0001.pdb”, "xxxx_Model_0002.pdb", etc. In the mutant model,  
  optimized polar hydrogen coordinates are also shown.

o To de novo design new sequence, you can run:
  
  ./EvoEF2 --command=ProteinDesign --monomer --pdb=monomer.pdb

  to design a monomer.
  
  or run:

  ./EvoEF2 --command=ProteinDesign --ppint --design_chains=A --pdb=complex.pdb
  
  to design the chain A of a complex.

  For a multi-chain complex (e.g., ABCD.pdb), if you want to design the first 
  two chains, you can run:

  ./EvoEF2 --command=ProteinDesign --ppint --design_chains=AB --pdb=ABCD.pdb

  Note that both backbone-depdent and backbone-indepdent rotamer libraries
  are supported by EvoEF2, you can also specify the library that you want to use:

  Backbone-depdent rotamer library:

  >> dun2010bb.lib
  >> dun2010bb1per.lib (exclude the rotamers with probability < 0.01 from dun2010bb.lib)
  >> dun2010bb3per.lib (exclude the rotamers with probability < 0.03 from dun2010bb.lib)

  ------------------------------------
  Backbone-indepdent rotamer library:
  ------------------------------------

  >> honig984.lib (984 rotamers for all 20 amino acids)
  >> honig3222.lib (3222 rotamers for all 20 amino acids)
  >> honig7421.lib (7421 rotamers for all 20 amino acids)
  >> honig11810.lib (11810 rotamers for all 20 amino acids)

  By default, the dun2010bb3per.lib is used for fast speed, if you want to switch
  to another backbone-dependent rotamer library, you can specify:

  ./EvoEF2 --command=ProteinDesign --monomer --rotlib=dun2010bb  --pdb=monomer.pdb

  Note that you do not have to add the suffix '.lib' when you specify the rotamer library

  or swith to a backbone-independent library by using:

  ./EvoEF2 --command=ProteinDesign --monomer --bbdep=disable  --rotlib=honig984  --pdb=monomer.pdb

  we strongly suggest you use the default settings to acheive both good performance and speed.


Cost and Availability
---------------------
The standalone EvoEF2 program and its source code is freely available to users. 
But unauthorized copying of the source code files via any medium is strictly prohibited.


Disclaimer and Copyright
------------------------
Copyright (c) Xiaoqiang Huang


References
----------
If EvoEF2 is important to your work, please cite: 
  
Xiaoqiang Huang, Robin Pearce, Yang Zhang. EvoEF2: accurate and fast energy function 
for computational protein design. Bioinformatics (2020), 36:1135-1142. 
https://doi.org/10.1093/bioinformatics/btz740

Bugs, comments and suggestions
------------------------------
For suggestions or bug reporting, please contact
xiaoqiah@outlook.com or xiaoqiah@umich.edu
