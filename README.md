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

  1. The force-field design relies on an [electronegativity-equalization scheme](https://github.com/oliveirampo/opt/blob/master/scr/EEM.py)
for the atomic partial charges [(J. Chem. Phys. 131, 044127 (2009))](https://aip.scitation.org/doi/10.1063/1.3187034).
  2. The optimization procedure uses statistical-mechanical expressions.
  
## Files Description

TODO

## Installation

TODO

## Usage

These are the options:

  * GEN: generate param.mod files to run simulations.
  * ANA: perform analysis of the simulation results.
  * OPT: optimize force-field parameters.
  * SUB: submit jobs to run on [euler](https://scicomp.ethz.ch/wiki/Main_Page).
  * PLOT: plot main simulation results.
  
Execute the pipeline with the following command:

python <path to source code>/run.py $OPTION $IT
  
where OPTION = (GEN, ANA, OPT, SUB, PLOT)
and IT is the iteration number.

