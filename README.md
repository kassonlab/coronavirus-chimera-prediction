# coronavirus-chimera-prediction
Code for generating and predicting coronavirus chimeras.

This repository is described in Simpson and Kasson, "Structural prediction of chimeric immunogens to elicit targeted antibodies against betacoronaviruses".

Here, we provide a set of Python functions and scripts to do the following:
  Create arbitrary chimeras between protein sequences
  Predict structure and pLDDT of these chimeras using AlphaFold
  Score the resulting chimeras for stability using the AlphaFold outputs
  Utilize the above functions to create S1/S2 chimeras between SARS-CoV-2 and a large set of betacoronaviruses
  Computationally validate chimera predictions using molecular dynamics simulations.

Functionality of specific python files is summarized below:

S1ChimeraProduction.py creates chimeras between a reference sequence (SARS-CoV-2) and a user-supplied list of accession numbers
example:
python3 S1ChimeraProduction.py listofaccessionandproteinnames.tsv S1
will create a set of chimeric fasta files from listofaccessionandproteinnames.tsv and name them as "S1" spliced chimeras.

MultimerAlphaFold.sh is a shell script to run AlphaFold

TODO: AnalysisSettings.py analyzes AlphaFold outputs to score chimeric sequences

md_sims/setup.py will set up molecular dynamics simulations using Gromacs.  This requires Gromacs as well as the CHARMM36 forcefield package.
As currently formulated, setup.py will take all PDB files in the current directory beginning with "3mer" and prepare production Gromacs run input files (parsing, energy minimization, and equilibration).

The following files contain utility routines
AccessiontoAlignment.py contains functions necesary for: creating fasta files for the sequences attached to the accession numbers you've collected, creating a concattnated fasta file for multiple sequence alignment, and also using that msa to find homologous sequences useful for splicing.

ChimeraGenerator.py contains functions that splice and then recombine sequence segments of your choice 

Analysis.py contains routines for analysis of AlphaFold outputs

