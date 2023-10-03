# PepMothSel
This repository holds code for the Biston betularia project. Specifically, the project component that aims to estimate the parameters of selection and mutation rate operating at the gene cortex during a period of industrial pollution.

We have two approaches to consider:

# Pascal's/Linglong's trajectory simulation
1) (Bespoke code) Simulate erosion of haplotypes through recombination based on a given frequency trajectory.

Problems: We don't have the same quality of frequency trajectory data for Europe.

See Linglong's overleaf document: https://www.overleaf.com/project/64bf9a4929d5c439e3af5feb

# A brute force simulation approach
1) (msprime/tsinfer) Perform a backwards simulation, to reach a starting point pre-industrialisation.
2) (Slim) Construct a model to simulate an onset of selection of a specific strength and with a specific mutation rate.
3) (SBI) Simulate over a parameter space using some priors for selection strength and mutation rate and compare to empirical summary statistics.

Problems: Model mis-specification. Priors - we have can get reasonable priors for selection based on U.K. studies, however TE induced mutation rates may be more challenging.

As an initial experiment we constructed a simple two deme simulation (polluted/unpolluted). Melanism is also simplified to typica vs melanic with selection coefficients for each depending on the population/timepoint. The slim file is available in the repository.

The characteristics of this simulation, are that most functional melanic mutations that arise, are lost very quickly. When selection begins the number of alleles present in the two demes is stochastic. If a single allele is present we see a classic hard sweep, followed by a subsequent sweep in the opposing direction.

<img width="682" alt="Screenshot 2023-09-19 at 11 20 27" src="https://github.com/swomics/PepMothSel/assets/44398000/dcea8b24-ed84-47bb-bc81-6a1cc6c93611">

Figure 1: The frequency trajectories of melanic alleles in the two deme simulation.

When multiple functional alleles are present we see a soft sweep favouring the multiple melanic alleles, followed by a hard sweep in the opposing direction favouring the ancestral non-melanic background.

<img width="682" alt="Screenshot 2023-09-19 at 11 20 55" src="https://github.com/swomics/PepMothSel/assets/44398000/7df9269d-25ba-431b-ab98-81a604ab9e6c">

Figure 2: An independent simulation run, showing alternate frequency trajectories of melanic alleles in the two deme simulation.




# Specifying a reaonsable demographic model in Slim

We need our model to be able to capture differences in industrialisation onset and enaction of clean air legislation, and also assess whether melanic alleles could spread alongside this anthropogenic activity or are more likely to have arisen/swept from indpendent events.


<img width="1026" alt="Industrial_rev_diffusion" src="https://github.com/swomics/PepMothSel/assets/44398000/bb582425-5b33-4875-b556-e5809206a690">
Figure 3: Our rough understanding of the onset of industrialisation.

![Screenshot 2023-09-19 at 11 31 32](https://github.com/swomics/PepMothSel/assets/44398000/68a9069b-d821-4857-bd47-280b04f304b4)

Figure 4: An agreed population model consisting of polluted and neighbouring unpolluted demes, roughly corresponding to the areas of industrial development from West to East.

<img width="706" alt="Screenshot 2023-09-15 at 10 38 01" src="https://github.com/swomics/PepMothSel/assets/44398000/0ad113b6-dce1-475c-ae5e-a0c80fb57cc5">

Figure 5: The same population model, as encoded and depicted in Slim:




##

export XDG_RUNTIME_DIR=/scratch/

conda activate msprime2

python msprime_initialise.py

python sbi_inference.py






