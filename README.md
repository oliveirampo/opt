# Force-field optimization

* This project is part of an integrated scheme for the automated refinement of force-field parameters 
against experimental condensed-phase data,
considering entire classes of organic molecules
constructed using a fragment library via combinatorial isomer enumeration.

* The main steps of the scheme, referred to as CombiFF, are:

  1. definition of a molecule family;
  2. combinatorial enumeration of all isomers;
  3. [query for experimental data](https://github.com/oliveirampo/combiff);
  4. automatic construction of the molecular topologies by fragment assembly;
  5. iterative refinement of the force-field parameters considering the entire family.
  
![](/images/TOC.png)

[J. Chem. Theory Comput., 16, 7525-​7555 (2020).](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00683)

* This project is only concerned about topic (iv).

* There are two key aspects in this scheme:

  1. The force-field design relies on an
     [electronegativity-equalization scheme](https://github.com/oliveirampo/opt/blob/master/scr/ChargeDistribution.py)
for the atomic partial charges [(J. Chem. Phys. 131, 044127 (2009))](https://aip.scitation.org/doi/10.1063/1.3187034).
  2. The optimization procedure uses statistical-mechanical expressions.
  
## Files Description

The script files are found in the scr/ directory.

Required directories:
  - bin/: Directory with SAMOS executable.
  - scr/: Directory with source code.
  - wrkDir/: working directory at the same level of the scr/ directory.
    Run the pipeline from this directory.
    
  - wrkDir/00_inp/: Directory of input files.
    - mol.dat: File with list of molecules.
    - Conf.dat: Configuration file.
    - listAtom.dat: File with list of atoms for each molecule.
    - listAng.dat: File with list of angles for each molecule.
    - listBond.dat: File with list of bonds for each molecule.
    - matrix.dat: File that specifies usage of C12(II) Lennard-Jones parameters.
    - model.sam: Template of SAMOS input file.
    - prm_0.dat: File with initial guess of parameters.
    - symmetry_sig.dat: File that specifies which parameters are optimized simultaneously in terms of sigma.
    - symmetry_eps.dat: File that specifies which parameters are optimized simultaneously in terms of epsilon.
  
  - wrkDir/cfg/: Directory with configurations.
    - 00_{molecule_cod}_liq_ini.cfg: Default configuration in the liquid phase.
    - 00_{molecule_cod}_gas_ini.cfg: Default configuration in the gas phase.
  - wrkDir/top/: Directory with topologies.
  - wrkDir/trc/: Directory with trajectories.
  - wrkDir/prm/: Directory with default IFP file:
    - FILE.ifp: Default interaction parameter file (IFP).

## Installation

This pipeline uses Python 3.7
and default packages for data visualization
and statistical analysis.

Additional packages:
  - scipy

## Usage

These are the options:

  * GEN: generate param.mod files to run simulations.
  * ANA: perform analysis of the simulation results.
  * OPT: optimize force-field parameters.
  * SUB: submit jobs to run on [euler](https://scicomp.ethz.ch/wiki/Main_Page).
  * PLOT: plot main simulation results.
  
Execute the pipeline with the following command:

> python $PATH/run.py $OPTION $IT
  
where
PATH = path to source code
OPTION = (GEN, ANA, OPT, SUB, PLOT)
and IT is the iteration number.
