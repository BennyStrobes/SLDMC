import numpy as np
import os
import sys
import pdb
import gzip
import pickle
import argparse
from pandas_plink import read_plink
import time


def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')





def create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file, variant_id_to_genotype_sdev, standardize=True):
	f = gzip.open(est_borzoi_effect_size_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		var_id = data[1]
		chrom_num = data[2]
		snp_pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		effect = float(data[6])

		if var_id not in variant_id_to_genotype_sdev:
			continue

		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()

		geno_sdev = variant_id_to_genotype_sdev[var_id]
		if standardize:
			effect = effect*geno_sdev

		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, effect)
	f.close()
	return mapping


def create_mapping_from_gene_id_to_variant_gene_annotations(sim_variant_gene_annotation_file):
	# New format: columns 6+ are one integer column per annotation, holding the category index
	# the variant-gene pair falls in (-1 if the pair falls in no category of that annotation).
	f = gzip.open(sim_variant_gene_annotation_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			anno_names = np.asarray(data[6:])
			continue
		gene_id = data[0]
		var_id = data[1]
		chrom_num = data[2]
		snp_pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		anno = np.asarray(data[6:]).astype(int)
		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, anno)
	f.close()
	return mapping, anno_names


def extract_annotation_categories(annotation_category_file, anno_names):
	# Companion file to the annotation file. Columns (tab-separated, with header):
	# anno_name  source  category_index  category_name
	# Returns, per annotation, the ordered list of its category names.
	anno_name_to_category_names = {}
	anno_name_to_source = {}
	f = open(annotation_category_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		anno_name = data[0]
		source = data[1]
		category_index = int(data[2])
		category_name = data[3]
		if anno_name not in anno_name_to_category_names:
			anno_name_to_category_names[anno_name] = []
			anno_name_to_source[anno_name] = source
		# Categories are written in index order, so appending keeps them aligned to the index
		if category_index != len(anno_name_to_category_names[anno_name]):
			print('assumption eroror: categories not in index order for ' + anno_name)
			pdb.set_trace()
		anno_name_to_category_names[anno_name].append(category_name)
	f.close()

	# Every annotation column in the annotation file needs an entry in the category file
	for anno_name in anno_names:
		if anno_name not in anno_name_to_category_names:
			print('assumption eroror: ' + anno_name + ' missing from category file')
			pdb.set_trace()

	return anno_name_to_category_names, anno_name_to_source


def create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file):	
	variant_id_to_geno_sdev = {}
	f = gzip.open(est_eqtl_effect_size_file,'rt')
	mapping = {}
	head_count = 0
	bad_genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		var_id = data[1]
		chrom_num = data[2]
		pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		effect = float(data[6])
		se = float(data[7])
		eqtl_sample_size = float(data[8])
		maf = float(data[9])
		if len(data) < 12:
			genotype_sdev = np.sqrt(2.0*maf*(1.0-maf))
		else:
			genotype_sdev = float(data[11])

		std_effect = effect*genotype_sdev
		std_se = se*genotype_sdev
		variant_id_to_geno_sdev[var_id] = genotype_sdev
		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('repeat snp error')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, pos, a0, a1, std_effect, std_se)
	f.close()
	return mapping, variant_id_to_geno_sdev

def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
	mapping = {}

	n_snps = len(ordered_snps)
	for snp_iter in range(n_snps):
		snp_name = ordered_snps[snp_iter]

		if snp_name in mapping:
			print('asssumption erororo')
			pdb.set_trace()
		mapping[snp_name] = snp_iter

	return mapping

def create_mapping_from_variant_id_to_snp_info(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
	if len(snp_array) != len(a0_arr):
		print('assumption eorroro')
		pdb.set_trace()
	if len(snp_array) != len(a1_arr):
		print('assumption eorroro')
		pdb.set_trace()

	dicti = {}

	for ii, snp_id in enumerate(snp_array):
		if snp_id in dicti:
			print('assumpationoenroer')
			pdb.set_trace()
		dicti[snp_id] = (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii])
	return dicti

def extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_est_borzoi_effects, var_to_est_eqtl_effects):
	unique_vars = np.unique(np.hstack(([*var_to_est_borzoi_effects],[*var_to_est_eqtl_effects])))
	final_vars = []
	for var in unique_vars:
		if var not in rsid_to_genotype_index:
			continue
		geno_alleles = (rsid_to_snp_info[var][0], rsid_to_snp_info[var][1])

		passing = True
		if var in var_to_est_borzoi_effects:
			borzoi_alleles = var_to_est_borzoi_effects[var][4:6]
			if set(geno_alleles) != set(borzoi_alleles):
				passing = False
		if var in var_to_est_eqtl_effects:
			eqtl_alleles = var_to_est_eqtl_effects[var][4:6]
			if set(geno_alleles) != set(eqtl_alleles):
				passing = False

		if passing == False:
			continue
		final_vars.append(var)
	return np.asarray(final_vars)


def load_in_snp_gene_eqtl_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []
	effect_ses = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.nan)
			effect_ses.append(np.nan)
			alleles.append(('nan', 'nan'))
		else:
			var_info = var_to_est_eqtl_effects[variant_id]
			effects.append(var_info[6])
			effect_ses.append(var_info[7])
			alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles), np.asarray(effect_ses)

def load_in_snp_gene_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.nan)
			alleles.append(('nan', 'nan'))
		else:
			var_info = var_to_est_eqtl_effects[variant_id]
			effects.append(var_info[6])
			alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles)


def load_in_snp_gene_anno_data(ordered_cis_variants, var_to_variant_gene_anno, n_anno):
	# Returns a (n_variant X n_annotation) matrix of category indices. A variant absent from the
	# annotation file gets all -1 (no category in any annotation); such variants are eqtl-only and
	# get subset out on the borzoi end anyways.
	annos = []
	alleles = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_variant_gene_anno:
			annos.append(np.full(n_anno, -1, dtype=int))
			alleles.append(('nan', 'nan'))
		else:
			var_info = var_to_variant_gene_anno[variant_id]
			annos.append(var_info[6])
			alleles.append((var_info[4], var_info[5]))
	return np.vstack(annos), np.asarray(alleles)


def create_one_hot_from_category_indices(category_indices, n_categories):
	# Turn a vector of category indices (-1 = no category) into a one-hot matrix.
	one_hot = np.zeros((len(category_indices), n_categories))
	for var_iter, category_index in enumerate(category_indices):
		if category_index >= 0:
			one_hot[var_iter, category_index] = 1.0
	return one_hot

def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	chrom_num = var_id_to_est_borzoi_effects[var_id][2]
	return chrom_num

def generate_gene_ld_means(gene_id_to_est_borzoi_effects,gene_id_to_est_borzoi_effects_unstandardized, gene_id_to_est_eqtl_effects, gene_id_to_variant_gene_anno, genotype_plink_filestem, anno_names, anno_name_to_n_categories, genotype_sample_indices, weighted):
	# Initialize output object
	gene_to_ld_means = {}
	
	# Loop through chromsomes
	for chrom_num in range(1,23):
		print(chrom_num)

		##################################
		# Load in per-chrom-genotype data
		##################################
		# string of chromosome name
		chrom_string = 'chr' + str(chrom_num)
		# Load in chromosome plink data
		(bim, fam, G) = read_plink(genotype_plink_filestem + str(chrom_num))
		# Create mapping from variant id to index
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
		# Create mapping from rsid to a0, a1
		rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))


		##################################
		# Loop through genes on this chromosome
		# (Analysis done seperately for each gene)
		##################################
		for gene_id in [*gene_id_to_est_borzoi_effects]:

			# Limit to genes on this chromosome
			gene_chrom_num = extract_gene_chrom_num(gene_id_to_est_borzoi_effects[gene_id])
			if str(gene_chrom_num) != str(chrom_num):
				continue

			# Gene needs both borzoi effects AND eQTLs
			if gene_id not in gene_id_to_est_eqtl_effects:
				continue
			if gene_id not in gene_id_to_variant_gene_anno:
				continue

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id], gene_id_to_est_eqtl_effects[gene_id])
			# Sip genes with fewer than 10 variants
			if len(ordered_cis_variants) < 10:
				continue

			# Load in data for gene
			# eQTL
			eqtl_effects, eqtl_variant_alleles, eqtl_effect_ses = load_in_snp_gene_eqtl_data(ordered_cis_variants, gene_id_to_est_eqtl_effects[gene_id])
			# Borzoi
			borzoi_effects, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])
			borzoi_effects_unstandardized, borzoi_variant_alleles_unstandardized = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects_unstandardized[gene_id])

			# Anno
			variant_anno, borzoi_anno_variant_alleles = load_in_snp_gene_anno_data(ordered_cis_variants, gene_id_to_variant_gene_anno[gene_id], len(anno_names))

			# Load in LD
			cis_genotype_indices = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				snp_info = rsid_to_snp_info[cis_variant]
				geno_alleles = snp_info[:2]
				
				# Also flip signs of eqtls to match LD
				if np.isnan(eqtl_effects[var_index]) == False:
					if eqtl_variant_alleles[var_index,:][0] == geno_alleles[1]:
						eqtl_effects[var_index] = -1.0*eqtl_effects[var_index]
				if np.isnan(borzoi_effects[var_index]) == False:
					if borzoi_variant_alleles[var_index,:][0] == geno_alleles[1]:
						borzoi_effects[var_index] = -1.0*borzoi_effects[var_index]
						borzoi_effects_unstandardized[var_index] = -1.0*borzoi_effects_unstandardized[var_index]
				if borzoi_variant_alleles[var_index,:][0] != borzoi_anno_variant_alleles[var_index,:][0]:
					print('annotation alllele assumption erororo')
					pdb.set_trace()
			
			# Extract genotype
			cis_genotype_indices = np.asarray(cis_genotype_indices)
			# Extract genotype matrix
			geno_mat = (G[cis_genotype_indices,:].compute())[:, genotype_sample_indices]
			row_means = np.nanmean(geno_mat, axis=1)
			nan_rows, nan_cols = np.where(np.isnan(geno_mat))
			geno_mat[nan_rows, nan_cols] = row_means[nan_rows]
			LD = np.corrcoef(geno_mat)

			# Subset LD by missingness
			# A. on eQTL end
			observed_eqtl_indices = np.isnan(eqtl_effects) == False
			eqtl_effects = eqtl_effects[observed_eqtl_indices]
			eqtl_effect_ses = eqtl_effect_ses[observed_eqtl_indices]
			LD = LD[observed_eqtl_indices, :]
			# B. on borzoi end
			observed_borzoi_indices = np.isnan(borzoi_effects) == False
			borzoi_effects = borzoi_effects[observed_borzoi_indices]
			borzoi_effects_unstandardized = borzoi_effects_unstandardized[observed_borzoi_indices]
			variant_anno = variant_anno[observed_borzoi_indices, :]
			LD = LD[:, observed_borzoi_indices]

			# Add to dictionary
			if gene_id in gene_to_ld_means:
				print('assumption erororo')
				pdb.set_trace()

			NNN = geno_mat.shape[1]
			sq_LD = LD**2
			adj_sq_ld = sq_LD - (1.0 - sq_LD)/(NNN - 2.0)

			# These pieces are annotation-independent, so compute them once per gene
			eqtl_var_y = np.square(eqtl_effects) - np.square(eqtl_effect_ses)
			if weighted:
				# LDSC redundancy weights: w_i = 1 / (intercept LD score), clamped at 1.
				intercept_LD_scores = np.sum(adj_sq_ld, axis=1)
				intercept_LD_scores[intercept_LD_scores < 1.0] = 1
				weights = 1.0/intercept_LD_scores

			# Run the model seperately for each annotation. Each annotation is reduced to per-gene
			# sufficient statistics (small matrices summable across genes for the bootstrap).
			gene_to_ld_means[gene_id] = {}
			for anno_iter, anno_name in enumerate(anno_names):
				n_categories = anno_name_to_n_categories[anno_name]
				# One-hot the variants over this annotation's categories (-1 index -> all-zero row)
				one_hot = create_one_hot_from_category_indices(variant_anno[:, anno_iter], n_categories)

				# Averager sufficient stats: A^T A and category column sums (over borzoi-observed variants)
				anno_xt_x = one_hot.T @ one_hot
				anno_col_sums = np.sum(one_hot, axis=0)

				# Mean (calibration) regression: X = ld_means, y = eQTL effect sizes.
				ld_means = LD @ (one_hot * borzoi_effects[:, None])
				calibration_valid = np.isfinite(eqtl_effects) & np.all(np.isfinite(ld_means), axis=1)
				calibration_x = ld_means[calibration_valid, :]
				calibration_y = eqtl_effects[calibration_valid]

				# Borzoi variance regression: X = one_hot, y = borzoi effect^2 (X^T X is anno_xt_x).
				borzoi_var_xt_y = one_hot.T @ np.square(borzoi_effects)
				borzoi_var_xt_y_unstandardized = one_hot.T @ np.square(borzoi_effects_unstandardized)

				# Variance (per-SNP eQTL h2) regression: X = ld_scores, y = eQTL effect^2 - se^2.
				ld_scores = adj_sq_ld @ one_hot
				eqtl_var_valid = np.isfinite(eqtl_var_y) & np.all(np.isfinite(ld_scores), axis=1)
				eqtl_var_x = ld_scores[eqtl_var_valid, :]
				eqtl_var_y_valid = eqtl_var_y[eqtl_var_valid]
				if weighted:
					weights_valid = weights[eqtl_var_valid]
					eqtl_var_xt_x = (eqtl_var_x * weights_valid[:, None]).T @ eqtl_var_x
					eqtl_var_xt_y = (eqtl_var_x * weights_valid[:, None]).T @ eqtl_var_y_valid
				else:
					eqtl_var_xt_x = eqtl_var_x.T @ eqtl_var_x
					eqtl_var_xt_y = eqtl_var_x.T @ eqtl_var_y_valid

				gene_to_ld_means[gene_id][anno_name] = {}
				gene_to_ld_means[gene_id][anno_name]['anno_xt_x'] = anno_xt_x
				gene_to_ld_means[gene_id][anno_name]['anno_col_sums'] = anno_col_sums
				gene_to_ld_means[gene_id][anno_name]['calibration_xt_x'] = calibration_x.T @ calibration_x
				gene_to_ld_means[gene_id][anno_name]['calibration_xt_y'] = calibration_x.T @ calibration_y
				gene_to_ld_means[gene_id][anno_name]['borzoi_var_xt_y'] = borzoi_var_xt_y
				gene_to_ld_means[gene_id][anno_name]['borzoi_var_xt_y_unstandardized'] = borzoi_var_xt_y_unstandardized
				gene_to_ld_means[gene_id][anno_name]['eqtl_var_xt_x'] = eqtl_var_xt_x
				gene_to_ld_means[gene_id][anno_name]['eqtl_var_xt_y'] = eqtl_var_xt_y


	return gene_to_ld_means


def run_ld_corr_single_annotation(anno_name, ordered_gene_names, gene_id_to_ld_means):
	# Accumulate this annotation's per-gene sufficient statistics across genes
	anno_xt_x = None
	anno_col_sums = None
	calibration_xt_x = None
	calibration_xt_y = None
	borzoi_var_xt_y = None
	borzoi_var_xt_y_unstandardized = None
	eqtl_var_xt_x = None
	eqtl_var_xt_y = None
	for gene_name in ordered_gene_names:
		gene_data = gene_id_to_ld_means[gene_name][anno_name]
		if anno_xt_x is None:
			anno_xt_x = np.zeros_like(gene_data['anno_xt_x'])
			anno_col_sums = np.zeros_like(gene_data['anno_col_sums'])
			calibration_xt_x = np.zeros_like(gene_data['calibration_xt_x'])
			calibration_xt_y = np.zeros_like(gene_data['calibration_xt_y'])
			borzoi_var_xt_y = np.zeros_like(gene_data['borzoi_var_xt_y'])
			borzoi_var_xt_y_unstandardized = np.zeros_like(gene_data['borzoi_var_xt_y_unstandardized'])
			eqtl_var_xt_x = np.zeros_like(gene_data['eqtl_var_xt_x'])
			eqtl_var_xt_y = np.zeros_like(gene_data['eqtl_var_xt_y'])
		anno_xt_x = anno_xt_x + gene_data['anno_xt_x']
		anno_col_sums = anno_col_sums + gene_data['anno_col_sums']
		calibration_xt_x = calibration_xt_x + gene_data['calibration_xt_x']
		calibration_xt_y = calibration_xt_y + gene_data['calibration_xt_y']
		borzoi_var_xt_y = borzoi_var_xt_y + gene_data['borzoi_var_xt_y']
		borzoi_var_xt_y_unstandardized = borzoi_var_xt_y_unstandardized + gene_data['borzoi_var_xt_y_unstandardized']
		eqtl_var_xt_x = eqtl_var_xt_x + gene_data['eqtl_var_xt_x']
		eqtl_var_xt_y = eqtl_var_xt_y + gene_data['eqtl_var_xt_y']

	# Solve each regression, then average the predicted value within each category.
	# Averaging a per-category prediction is (A^T A @ coefs) / (A^T A column sums).
	calibration_coefs, _, _, _ = np.linalg.lstsq(calibration_xt_x, calibration_xt_y, rcond=None)
	borzoi_var_coefs, _, _, _ = np.linalg.lstsq(anno_xt_x, borzoi_var_xt_y, rcond=None)
	borzoi_var_coefs_unstandardized, _, _, _ = np.linalg.lstsq(anno_xt_x, borzoi_var_xt_y_unstandardized, rcond=None)
	taus, _, _, _ = np.linalg.lstsq(eqtl_var_xt_x, eqtl_var_xt_y, rcond=None)

	avg_calibration_slopes = np.dot(anno_xt_x, calibration_coefs)/anno_col_sums
	avg_borzoi_variances = np.dot(anno_xt_x, borzoi_var_coefs)/anno_col_sums
	avg_unstandardized_borzoi_variances = np.dot(anno_xt_x, borzoi_var_coefs_unstandardized)/anno_col_sums
	avg_per_snp_h2 = np.dot(anno_xt_x, taus)/anno_col_sums

	return avg_calibration_slopes, avg_per_snp_h2, avg_borzoi_variances, avg_unstandardized_borzoi_variances


def run_ld_corr(anno_names, ordered_gene_names, gene_id_to_ld_means):
	# Run the model seperately for each annotation
	anno_name_to_results = {}
	for anno_name in anno_names:
		anno_name_to_results[anno_name] = run_ld_corr_single_annotation(anno_name, ordered_gene_names, gene_id_to_ld_means)
	return anno_name_to_results








##########################
# Command line args
##########################
parser = argparse.ArgumentParser(description='Estimate eQTL-Borzoi LD correlations.')
parser.add_argument('--est-borzoi-effect-size-file', dest='est_borzoi_effect_size_file', required=True, help='Estimated Borzoi effect sizes file.')
parser.add_argument('--est-eqtl-effect-size-file', dest='est_eqtl_effect_size_file', required=True, help='Estimated eQTL effect sizes file.')
parser.add_argument('--sim-variant-gene-annotation-file', dest='sim_variant_gene_annotation_file', required=True, help='Variant-gene annotation file.')
parser.add_argument('--genotype-plink-filestem', dest='genotype_plink_filestem', required=True, help='Genotype plink filestem (per-chromosome number appended).')
parser.add_argument('--genotype-sample-mapping-file', dest='genotype_sample_mapping_file', required=True, help='Genotype sample indices for in-sample LD.')
parser.add_argument('--ld-corr-output-stem', dest='ld_corr_output_stem', required=True, help='Output filestem.')
parser.add_argument('--bootstrapped-gene-set-filestem', dest='bootstrapped_gene_sets_filestem', default=None, required=False, help='Filestem corresponding to a list of files with bootrapped gene names')
parser.add_argument('--weighted', dest='weighted', type=str2bool, default=True, help='Whether to use weighted regression (True/False).')
args = parser.parse_args()

est_borzoi_effect_size_file = args.est_borzoi_effect_size_file
est_eqtl_effect_size_file = args.est_eqtl_effect_size_file
sim_variant_gene_annotation_file = args.sim_variant_gene_annotation_file
genotype_plink_filestem = args.genotype_plink_filestem
genotype_sample_mapping_file = args.genotype_sample_mapping_file
ld_corr_output_stem = args.ld_corr_output_stem
weighted = args.weighted
bootstrapped_gene_sets_filestem = args.bootstrapped_gene_sets_filestem

# Companion file describing the (annotation, category) pairs in the annotation file
annotation_category_file = sim_variant_gene_annotation_file.split('.txt.gz')[0] + '_categories.txt'


##############################
# Load in data
##############################
# Create mapping from gene id to vector of est (standardized) eqtl effect sizes
gene_id_to_est_eqtl_effects, variant_id_to_genotype_sdev = create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file)


# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file, variant_id_to_genotype_sdev,standardize=True)
gene_id_to_est_borzoi_effects_unstandardized = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file, variant_id_to_genotype_sdev,standardize=False)


# Create mapping from gene id to vector of variant-gene annotations
gene_id_to_variant_gene_anno, anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(sim_variant_gene_annotation_file)

# Extract the categories making up each annotation
anno_name_to_category_names, anno_name_to_source = extract_annotation_categories(annotation_category_file, anno_names)
anno_name_to_n_categories = {}
for anno_name in anno_names:
	anno_name_to_n_categories[anno_name] = len(anno_name_to_category_names[anno_name])

# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

# Generate per gene ld-means
gene_id_to_ld_means = generate_gene_ld_means(gene_id_to_est_borzoi_effects, gene_id_to_est_borzoi_effects_unstandardized, gene_id_to_est_eqtl_effects, gene_id_to_variant_gene_anno, genotype_plink_filestem, anno_names, anno_name_to_n_categories, genotype_sample_indices, weighted)
del gene_id_to_est_eqtl_effects
del gene_id_to_variant_gene_anno
del gene_id_to_est_borzoi_effects
del gene_id_to_est_borzoi_effects_unstandardized


##############################
# Run regressions
##############################

# Run standard regression (all annotations at once)
ordered_gene_names = np.asarray([*gene_id_to_ld_means])
observed_results = run_ld_corr(anno_names, ordered_gene_names, gene_id_to_ld_means)

# Bootstrap estimates! The same resampled gene set is shared across annotations each iteration.
n_bs = 100
bs_results = []
for bs_iter in range(n_bs):
	print('bootstrap itera: ' + str(bs_iter))
	if bootstrapped_gene_sets_filestem is None:
		bs_genes = np.random.choice(ordered_gene_names, size=len(ordered_gene_names), replace=True)
	else:
		bootstrapped_gene_set_file = bootstrapped_gene_sets_filestem + str(bs_iter) + '.txt'
		bs_genes = np.loadtxt(bootstrapped_gene_set_file, dtype=str, skiprows=1)
		# Global gene sets are cross-tissue; keep only genes this tissue has data for
		bs_genes = bs_genes[np.isin(bs_genes, ordered_gene_names)]

	bs_results.append(run_ld_corr(anno_names, bs_genes, gene_id_to_ld_means))


##############################
# Organize results per annotation
##############################
# For each annotation, stack observed (length n_categories) and bootstrap (n_bs X n_categories) estimates
anno_name_to_observed = {}
anno_name_to_bootstrap = {}
for anno_name in anno_names:
	avg_calibration_slopes, avg_per_snp_eqtl_h2, avg_borzoi_variances, avg_unstandardized_borzoi_variances = observed_results[anno_name]

	bs_calibration_slopes = np.asarray([bs_results[bs_iter][anno_name][0] for bs_iter in range(n_bs)])
	bs_per_snp_eqtl_h2 = np.asarray([bs_results[bs_iter][anno_name][1] for bs_iter in range(n_bs)])
	bs_borzoi_variances = np.asarray([bs_results[bs_iter][anno_name][2] for bs_iter in range(n_bs)])
	bs_unstandardized_borzoi_variances = np.asarray([bs_results[bs_iter][anno_name][3] for bs_iter in range(n_bs)])

	# Compute correlation information
	n_categories = anno_name_to_n_categories[anno_name]
	observed_corr = np.full(n_categories, np.nan)
	for category_iter in range(n_categories):
		sqrt_term = avg_borzoi_variances[category_iter]/avg_per_snp_eqtl_h2[category_iter]
		if np.isfinite(sqrt_term) and sqrt_term >= 0.0:
			observed_corr[category_iter] = avg_calibration_slopes[category_iter]*np.sqrt(sqrt_term)
	bs_corr = np.full(bs_calibration_slopes.shape, np.nan)
	for bs_iter in range(bs_calibration_slopes.shape[0]):
		for category_iter in range(bs_calibration_slopes.shape[1]):
			sqrt_term = bs_borzoi_variances[bs_iter, category_iter]/bs_per_snp_eqtl_h2[bs_iter, category_iter]
			if np.isfinite(sqrt_term) and sqrt_term >= 0.0:
				bs_corr[bs_iter, category_iter] = bs_calibration_slopes[bs_iter, category_iter]*np.sqrt(sqrt_term)

	anno_name_to_observed[anno_name] = np.vstack((avg_calibration_slopes, avg_per_snp_eqtl_h2, avg_borzoi_variances, avg_unstandardized_borzoi_variances, observed_corr)).T
	anno_name_to_bootstrap[anno_name] = [bs_calibration_slopes, bs_per_snp_eqtl_h2, bs_borzoi_variances, bs_unstandardized_borzoi_variances, bs_corr]


#############################
# Print results to (shared) output file
##############################
output_names = ['calibration_slope', 'per_snp_eqtl_h2', 'borzoi_variance', 'unstandardized_borzoi_variance', 'correlation']

results_output_file = ld_corr_output_stem + '_bootstrap_summary.txt'
with open(results_output_file, 'w') as t:
	t.write('annotation_name\tcategory_name\tsample_label\tsample_iter\tcalibration_slope\tper_snp_eqtl_h2\tborzoi_variance\tunstandardized_borzoi_variance\tcorrelation\n')
	for anno_name in anno_names:
		category_names = anno_name_to_category_names[anno_name]
		observed_mat = anno_name_to_observed[anno_name]
		bs_calibration_slopes, bs_per_snp_eqtl_h2, bs_borzoi_variances, bs_unstandardized_borzoi_variances, bs_corr = anno_name_to_bootstrap[anno_name]
		for category_iter, category_name in enumerate(category_names):
			t.write(str(anno_name) + '\t' + str(category_name) + '\tobserved\t-1\t' + str(observed_mat[category_iter, 0]) + '\t' + str(observed_mat[category_iter, 1]) + '\t' + str(observed_mat[category_iter, 2]) + '\t' + str(observed_mat[category_iter, 3]) + '\t' + str(observed_mat[category_iter, 4]) + '\n')
		for bs_iter in range(n_bs):
			for category_iter, category_name in enumerate(category_names):
				t.write(str(anno_name) + '\t' + str(category_name) + '\tbootstrap\t' + str(bs_iter) + '\t' + str(bs_calibration_slopes[bs_iter, category_iter]) + '\t' + str(bs_per_snp_eqtl_h2[bs_iter, category_iter]) + '\t' + str(bs_borzoi_variances[bs_iter, category_iter]) + '\t' + str(bs_unstandardized_borzoi_variances[bs_iter, category_iter]) + '\t' + str(bs_corr[bs_iter, category_iter]) + '\n')

summary_stats_output_file = ld_corr_output_stem + '_bootstrap_stats.txt'
with open(summary_stats_output_file, 'w') as t:
	t.write('annotation_name\tcategory_name\toutput_name\tmean\tbootstrapped_mean\tbootstrap_se\tgaussian_z_score\tempirical_ci_lower\tempirical_ci_upper\n')
	for anno_name in anno_names:
		category_names = anno_name_to_category_names[anno_name]
		observed_mat = anno_name_to_observed[anno_name]
		bootstrap_mats = anno_name_to_bootstrap[anno_name]
		for category_iter, category_name in enumerate(category_names):
			for output_iter, output_name in enumerate(output_names):
				bootstrap_distribution = bootstrap_mats[output_iter][:, category_iter]
				bootstrap_mean = np.mean(bootstrap_distribution)
				bootstrap_se = np.std(bootstrap_distribution, ddof=1)
				if bootstrap_se == 0.0:
					gaussian_z_score = np.nan
				else:
					gaussian_z_score = observed_mat[category_iter, output_iter]/bootstrap_se
				empirical_ci = np.quantile(bootstrap_distribution, [.025, .975])
				t.write(str(anno_name) + '\t' + str(category_name) + '\t' + output_name + '\t' + str(observed_mat[category_iter, output_iter]) + '\t' + str(bootstrap_mean) + '\t' + str(bootstrap_se) + '\t' + str(gaussian_z_score) + '\t' + str(empirical_ci[0]) + '\t' + str(empirical_ci[1]) + '\n')

print(results_output_file)
print(summary_stats_output_file)
