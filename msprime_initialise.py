import json
import cyvcf2
import tsinfer
import os
import multiprocessing
import cyvcf2
import tsinfer
import numpy as np
import pyslim
import msprime

def add_populations(vcf, samples):
    """
    Add tsinfer Population objects and returns a list of IDs corresponding to the VCF samples.
    """
    # In this VCF, the first letter of the sample name refers to the population

    samples_first_letter = [sample_name.split("-")[0] for sample_name in vcf.samples]
    return samples_first_letter




def chromosome_length(vcf):
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]


def add_diploid_sites(vcf, samples, chromosome=None, start=None, stop=None):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
    pos = 0
    if chromosome is None or start is None or stop is None:
        vcf_iter = vcf
        start = pos = 0  # VCFs are 1-based, so 0 is not in the first file
    else:
        # This is a partial file: use the standard HTSlib string spec
        vcf_iter = vcf(f"{chromosome}:{start}-{stop}")
        pos = start - 1
    for variant in vcf_iter:  # Loop over variants, each assumed at a unique site
        if variant.POS < start:
            # This could be a long variant from a previous section of the VCF that
            # extends into this file. It should have been capture previously, so we skip
            continue
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS

        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF] + variant.ALT
        ancestral = variant.INFO.get('AA', variant.REF)
        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        allele_index = {old_index: ordered_alleles.index(allele)
            for old_index, allele in enumerate(alleles)}
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [allele_index[old_index]
            for row in variant.genotypes for old_index in row[0:2]]
        samples.add_site(pos, genotypes=genotypes, alleles=ordered_alleles)

def read_vcf_part(params):
    vcf = cyvcf2.VCF(params[0])
    vcf.set_index(params[1])
    assert len(vcf.seqnames) == 1  # Check only one chromosome present
    assert len(vcf.seqlens) == 1
    sequence_length = vcf.seqlens[0]
    part = params[2]
    num_parts = params[3]
    assert sequence_length > num_parts
    chunk_size = sequence_length // num_parts + 1
    start = part*chunk_size
    stop = (part+1)*chunk_size - 1  # ranges are inclusive
    stop = min(stop, sequence_length)


    path = f"partial_file{start}-{stop}.samples"
    with tsinfer.SampleData(
        path=path,
        sequence_length=sequence_length) as samples:
            add_diploid_sites(vcf, samples, vcf.seqnames[0], start, stop)
    samples.close()  # No need to leave the file open
    return path

if __name__ == '__main__':
    # URL for the VCF
#    url = "Pangenie.1.2.merged_GQ200_mac3.FR989875_cortex.vcf.gz"
#    idx = "Pangenie.1.2.merged_GQ200_mac3.FR989875_cortex.vcf.gz.tbi"
#    num_threads = 12
#    with multiprocessing.Pool(num_threads) as pool:
#        sample_filenames = [r for r in pool.map(
#            read_vcf_part, ((url, idx, i, num_threads) for i in range(num_threads)))]
#    all_samples = [tsinfer.load(name) for name in sample_filenames]
#    samples = all_samples[0].copy()
#    samples.append_sites(*all_samples[1:])
#    samples.finalise()
#    print(
#        "Created full samples file with",
#        samples.num_sites, "sites and",
#        samples.num_samples, "samples")


	vcf_location = "Pangenie.1.2.merged_GQ200_mac3.FR989875_cortex.vcf.gz"
# NB: could also read from an online version by setting vcf_location to
# "https://github.com/tskit-dev/tsinfer/raw/main/docs/_static/P_dom_chr24_phased.vcf.gz"

	vcf = cyvcf2.VCF(vcf_location)
	with tsinfer.SampleData(
    	path="P_dom_chr24_phased.samples", sequence_length=chromosome_length(vcf)
	) as samples:
		add_diploid_sites(vcf, samples)

	print(
    	"Sample file created for {} samples ".format(samples.num_samples)
    	+ "({} individuals) ".format(samples.num_individuals)
    	+ "with {} variable sites.".format(samples.num_sites),
    	flush=True,
	)

# Do the inference
	ts = tsinfer.infer(samples)
	print(
    	"Inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        	ts.num_trees, ts.sequence_length / 1e6, ts.num_edges
    	)
	)


# Check the metadata
	for sample_node_id in ts.samples():
    		individual_id = ts.node(sample_node_id).individual
    		population_id = ts.node(sample_node_id).population
    		#print("Node", sample_node_id,individual_id,population_id)

import matplotlib.pyplot as plt

kb = [0]  # Starting genomic position
mrca_time = []


demog = msprime.Demography()
demog.add_population(initial_size=10000)
ts = msprime.sim_ancestry(
            samples=200,
            end_time=1000,
            demography=demog,
            recombination_rate=2.9e-9,
            sequence_length=14726131,
            random_seed=5)
ts = pyslim.annotate(ts, model_type="WF", tick=1)
assert ts.num_individuals == 200
assert ts.num_samples == 400

ts.dump("initialize_WF.trees")
