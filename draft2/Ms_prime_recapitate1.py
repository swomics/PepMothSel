import tskit
import pyslim
import numpy as np

ts = tskit.load("P_struct_out.trees")

# allele frequencies
p = ts.sample_count_stat(
                [ts.samples()], lambda x: x/20000, 1, windows='sites',
                span_normalise=False, polarised=True, strict=False)
print(f"There are {ts.num_sites} segregating sites, of which {np.sum(p > 0.25)}")
print(f"are at frequency above 25%, and {np.sum(p > 0.05)} are above 5%.")


rts = pyslim.recapitate(ts, ancestral_Ne=5000000, recombination_rate=4.4397e-8, random_seed=6)

# check type m0 is not used:
mut_types = set([md['mutation_type']
                for mut in ts.mutations()
                for md in mut.metadata['mutation_list']])
print(f"Keeping {rts.num_mutations} existing mutations of type(s) {mut_types}.")
assert 0 not in mut_types

# add type m0 mutations
next_id = pyslim.next_slim_mutation_id(rts)
rts = msprime.sim_mutations(
            rts, rate=1e-8, random_seed=7, keep=True,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id)
)

p = rts.sample_count_stat(
                [rts.samples()], lambda x: x/20000, 1, windows='sites',
                span_normalise=False, polarised=True, strict=False)
print(f"After mutation, there are {rts.num_sites} segregating sites, of which {np.sum(p > 0.25)}")
print(f"are at frequency above 25%, and {np.sum(p > 0.05)} are above 5%.")

