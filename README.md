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

ProductionScript.py creates chimeras between a reference sequence and a user-supplied list of sequences.

MultimerAlphaFold.sh is a shell script to run AlphaFold

AnalysisSettings.py analyzes AlphaFold outputs to score chimeric sequences

md_sims/setup.py will set up molecular dynamics simulations using Gromacs.  This requires Gromacs as well as the CHARMM36 forcefield package.
As currently formulated, setup.py will take all PDB files in the current directory beginning with "3mer" and prepare production Gromacs run input files (parsing, energy minimization, and equilibration).

The following files contain utility routines
AccessiontoAlignment.py contains functions necesary for: creating fasta files for the sequences attached to the accession numbers you've collected, creating a concattnated fasta file for multiple sequence alignment, and also using that msa to find homologous sequences useful for splicing.

ChimeraGenerator.py contains functions that splice and then recombine sequence segments of your choice 

Analysis.py contains routines for analysis of AlphaFold outputs

More detailed information:
ProductionScript.py creates chimeras between a reference sequence (ex:SARS-CoV-2) and a user-supplied list of accession numbers and names preferred names for the corresponding sequences.
example: python3 ProductionScript.py List_of_Coronaviruses.tsv chimera_arguments.json 6vsb_S1.fasta
will create a set of chimeric fasta files along with the fastas for their wild-type parent proteins named as designated in chimera_arguments.json monomer_fasta/multimer_fasta/chimera_fastas argument where the * in /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3mer*.fasta will be replaced with preferred names outlined in second column of List_of_Coronaviruses.tsv. The basename of the monomer fasta files will  be how the aligned sequences are identified in the multiple sequence alignment. If you input your own msa, make sure these are identical. Alphafold will name output folders after the basename of the fasta file that was input.

In both argument json files, chimera_arguments.json and analysis_arguments.json, some variables are enclosed in [] and accompanied by empty "", filling those "" with # will turn the creation of the associated files or commandline function off.
Ex. "emboss_command":["","/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle"] -> "emboss_command":["#","/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle"] Including this # in the empty quotes of the "emboss_command" variable from analysis_arguments.json will prevent emboss from running
Ex. "emboss_names": ["","/gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/*S1vsSARS2.emboss"] -> "emboss_names": ["","/gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/*S1vsSARS2.emboss"] this will prevent Emboss from running as well as negate emboss data or files from being included in analysis or data table output completely
Ex. "muscle_command_for_msa": ["","module load gcc / 9.2.0 & & module load muscle / 3.8.31 & & muscle"] -> "muscle_command_for_msa": ["#","module load gcc / 9.2.0 & & module load muscle / 3.8.31 & & muscle"] will turn off the creation of the muscle msa file, allowing for inclusion of a personally generated msa

