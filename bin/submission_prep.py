#!/usr/bin/env python3
import os
import shutil
import logging
import pandas as pd
from submission_helper import (
	GetParams,
	SubmissionConfigParser,
	Sample,
	BiosampleSubmission,
	SRASubmission,
	GenbankSubmission,
	get_compound_extension,
	setup_logging
)

# Module-level variable to store log file path for exception handling
_log_file_path = None

def prepare_sra_fastqs(samples, outdir, copy=False):
	for sample in samples:
		if sample.fastq1 and sample.fastq2:
			ext1 = get_compound_extension(sample.fastq1)
			ext2 = get_compound_extension(sample.fastq2)
			dest_fq1 = os.path.join(outdir, f"{sample.sample_id}_R1{ext1}")
			dest_fq2 = os.path.join(outdir, f"{sample.sample_id}_R2{ext2}")
			if not os.path.exists(dest_fq1):
				if copy:
					shutil.copy(sample.fastq1, dest_fq1)
				else:
					os.symlink(sample.fastq1, dest_fq1)
			if not os.path.exists(dest_fq2):
				if copy:
					shutil.copy(sample.fastq2, dest_fq2)
				else:
					os.symlink(sample.fastq2, dest_fq2)

def main_prepare():
	global _log_file_path
	# parse exactly the same CLI args you already have
	params = GetParams().parameters

	os.makedirs(params['outdir'], exist_ok=True)

	# Ensure log file is created early, even if script fails
	# Use relative path (from outdir) as Nextflow expects files relative to work directory
	log_file_path = os.path.join(params['outdir'], 'prep_submission.log')
	_log_file_path = log_file_path  # Store for exception handler
	
	# Create and write to log file directly to ensure it exists
	with open(log_file_path, 'w') as f:
		f.write('[INFO] Log file initialized\n')
		f.flush()
		os.fsync(f.fileno())  # Force write to disk
	
	setup_logging(log_file=log_file_path, level=logging.DEBUG)
	logging.info("Started logging for preparation.")
	# Force flush to ensure log file is written
	for handler in logging.getLogger().handlers:
		if isinstance(handler, logging.FileHandler):
			handler.flush()
	
	# load config & metadata
	config = SubmissionConfigParser(params).load_config()
	batch_id = os.path.splitext(os.path.basename(params['metadata_file']))[0]
	metadata_df = pd.read_csv(params['metadata_file'], sep='\t')
	identifier = params['identifier']
	submission_dir = 'Test' if params['test'] else 'Production'
	output_root = params['outdir']

	# build sample objects
	samples = []
	for s in params['sample']:
		d = dict(item.split('=') for item in s.split(','))
		samples.append(Sample(
			sample_id   = d['sample_id'],
			batch_id	= batch_id,
			fastq1	  = d.get('fq1'),
			fastq2	  = d.get('fq2'),
			nanopore	= d.get('nanopore'),
			species	 = params['species'],
			databases   = [db for db in params if params[db] and db in ['biosample','sra','genbank']],
			fasta_file  = d.get('fasta'),
			annotation_file = d.get('gff')
		))

	# 1) Prepare BioSample XML + submit.ready
	if params['biosample']:
		submission_dir = os.path.join(output_root, 'biosample')
		os.makedirs(submission_dir, exist_ok=True)
		bs = BiosampleSubmission(
			parameters=params,
			submission_config=config,
			metadata_df=metadata_df,
			outdir=submission_dir,
			submission_mode=params['submission_mode'],
			submission_dir=submission_dir,
			type='biosample',
			sample=None,
			accession_id=None,
			identifier=identifier,
			wastewater=params.get('wastewater', False)
		)
		bs.init_xml_root()
		for s in samples:
			md = metadata_df[metadata_df['sample_name'] == s.sample_id]
			bs.add_sample(s, md)
		bs.finalize_xml()
		# write submit.ready
		open(os.path.join(submission_dir,'submit.ready'),'w').close()

	# 2) Prepare SRA XML + submit.ready (per-platform if needed)
	if params['sra']:
		illum = [s for s in samples if s.fastq1 and s.fastq2]
		nano  = [s for s in samples if s.nanopore]
		platforms = (('illumina', illum), ('nanopore', nano)) if illum and nano else [(None, illum or nano)]
		for platform, samp_list in platforms:
			# submission_dir needs to be unique if submitting both illumina and nanopore
			submission_dir = os.path.join(output_root, 'sra', platform) if platform else os.path.join(output_root, 'sra')
			os.makedirs(submission_dir, exist_ok=True)
			sra = SRASubmission(
				parameters=params,
				submission_config=config,
				metadata_df=metadata_df,
				outdir=submission_dir,
				submission_mode=params['submission_mode'],
				submission_dir=submission_dir,
				type='sra',
				samples=samp_list,
				sample=None,
				accession_id=None,
				identifier=identifier,
				wastewater=params.get('wastewater', False)
			)
			sra.init_xml_root()
			for s in samp_list:
				md = metadata_df[metadata_df['sample_name'] == s.sample_id]
				sra.add_sample(s, md, platform)  # existing signature
			sra.finalize_xml()
			# write submit.ready
			open(os.path.join(submission_dir,'submit.ready'),'w').close()
			# copy/Symlink raw files to SRA folder
			prepare_sra_fastqs(samp_list, submission_dir, copy=False)
			
	# 3) Prepare GenBank submission, per-sample
	if params['genbank']:
		for s in samples:
			submission_dir = os.path.join(output_root, 'genbank', s.sample_id)
			os.makedirs(submission_dir, exist_ok=True)
			gb = GenbankSubmission(
				parameters=params,
				submission_config=config,
				metadata_df=metadata_df,
				outdir=submission_dir,
				submission_mode=params['submission_mode'],
				submission_dir=submission_dir,
				type='genbank',
				samples=samples,
				sample=s,               # class expects one sample
				accession_id=None,
				identifier=identifier
			)
			gb.genbank_submission_driver()
	
	# Ensure log file exists and is flushed before script exits
	# Flush all handlers but don't close them (closing might cause issues)
	for handler in logging.getLogger().handlers:
		if isinstance(handler, logging.FileHandler):
			handler.flush()
	
	# CRITICAL: Write directly to file to ensure it exists on disk
	# This is the last thing we do before returning
	# Use the exact relative path that Nextflow expects
	with open(log_file_path, 'a') as f:
		f.write('[INFO] Script completed successfully\n')
		f.flush()
		os.fsync(f.fileno())  # Force write to disk
	
	# Final verification - ensure file exists at the expected path
	if not os.path.exists(log_file_path):
		# Last resort: create the file
		with open(log_file_path, 'w') as f:
			f.write('[INFO] Log file created at end of script\n')
			f.write('[INFO] Script completed successfully\n')
			f.flush()
			os.fsync(f.fileno())
	
	# Verify file is readable and has content
	try:
		with open(log_file_path, 'r') as f:
			content = f.read()
			if not content.strip():
				# File exists but is empty, write something
				with open(log_file_path, 'a') as f2:
					f2.write('[INFO] Log file verified at end\n')
					f2.flush()
					os.fsync(f2.fileno())
	except Exception as e:
		# If we can't read it, recreate it
		with open(log_file_path, 'w') as f:
			f.write(f'[INFO] Log file recreated (read error: {str(e)})\n')
			f.write('[INFO] Script completed successfully\n')
			f.flush()
			os.fsync(f.fileno())

if __name__=="__main__":
	try:
		main_prepare()
		# Final flush of all handlers
		for handler in logging.getLogger().handlers:
			if isinstance(handler, logging.FileHandler):
				handler.flush()
		# Final write to ensure file exists
		if _log_file_path and os.path.exists(_log_file_path):
			with open(_log_file_path, 'a') as f:
				f.write('[INFO] Script exiting normally\n')
				f.flush()
				os.fsync(f.fileno())
	except Exception as e:
		# Log error if logging is already set up
		if logging.getLogger().handlers:
			logging.error(f"Fatal error in submission_prep.py: {str(e)}", exc_info=True)
			# Flush handlers even on error
			for handler in logging.getLogger().handlers:
				if isinstance(handler, logging.FileHandler):
					handler.flush()
		
		# CRITICAL: Ensure log file exists even on error
		# Use module-level variable to avoid scope issues
		if _log_file_path:
			try:
				with open(_log_file_path, 'a') as f:
					f.write(f'[ERROR] Script failed: {str(e)}\n')
					f.flush()
					os.fsync(f.fileno())
			except Exception as write_error:
				# If we can't write to the expected location, try to create it
				try:
					# Try to get outdir from GetParams if available
					try:
						params = GetParams().parameters
						log_file_path = os.path.join(params['outdir'], 'prep_submission.log')
						with open(log_file_path, 'w') as f:
							f.write(f'[ERROR] Script failed: {str(e)}\n')
							f.write(f'[ERROR] Original write error: {str(write_error)}\n')
							f.flush()
							os.fsync(f.fileno())
					except:
						pass
				except:
					pass  # Last resort - if we can't create the file, at least don't crash
		raise
