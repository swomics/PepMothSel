# PepMothSel
This repository holds code for the Biston betularia project. Specifically, the project component that aims to estimate the parameters of selection and mutation rate operating at the gene cortex during a period of industrial pollution.

We have two approaches to consider:

# Pascal's trajectory simulation
1) (Bespoke code) Simulate erosion of haplotypes through recombination based on a given frequency trajectory.

Problems: We don't have the same quality of frequency trajectory data for Europe.

See Linglong's overleaf document: https://www.overleaf.com/project/64bf9a4929d5c439e3af5feb

# A brute force simulation approach
1) (msprime/tsinfer) Perform a backwards simulation, to reach a starting point pre-industrialisation.
2) (Slim) Construct a modelto simulate an onset of selection of a specific strength and with a specific mutation rate.
3) (SBI) Simulate over a parameter space using some priors for selection strength and mutation rate and compare to empirical summary statistics.

Problems: Model mis-specification. Priors - we have can get reasonable priors for selection based on U.K. studies, however TE induced mutation rates may be more challenging.

![Screenshot 2023-09-19 at 11 07 02](https://github.com/swomics/PepMothSel/assets/44398000/3538656c-8e86-4fce-b32a-1b508cf6262e)

Figure 1: An agreed population model consisting of polluted and neighbouring unpolluted demes, roughly corresponding to the areas of industrial development from West to East.

<img width="706" alt="Screenshot 2023-09-15 at 10 38 01" src="https://github.com/swomics/PepMothSel/assets/44398000/5ac7a0a7-f7cd-4d56-a738-3d64f6003881">

Figure 2: The same population model, as encoded and depicted in Slim:

As an initial experiment we constructed a simple two deme model (polluted/unpolluted). Melanism is also simplified to typica vs melanic with selection coefficients for each depending on the population/timepoint.

The characteristics of this model, are that most functional melanic mutations that arise, are lost very quickly. When selection begins the number of alleles present in the two demes is stochastic. If a single allele is present we see a classic hard sweep, followed by a subsequent sweep in the opposing direction.

<img width="682" alt="Screenshot 2023-09-19 at 11 20 27" src="https://github.com/swomics/PepMothSel/assets/44398000/57da9976-bde2-44c5-ab23-503bb612e1d6">

Figure 3: The frequency trajectories of melanic alleles in the two deme simulation.

When multiple functional alleles are present we see a soft sweep favouring the multiple melanic alleles, followed by a hard sweep in the opposing direction favouring the ancestral non-melanic background.

<img width="682" alt="Screenshot 2023-09-19 at 11 20 55" src="https://github.com/swomics/PepMothSel/assets/44398000/e85ae9c4-6e49-456c-822f-5e796647a38a">

Figure 4: An independent simulation run, showing alternate frequency trajectories of melanic alleles in the two deme simulation.

