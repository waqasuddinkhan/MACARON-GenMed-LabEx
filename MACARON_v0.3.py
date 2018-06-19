#!/usr/bin/python

#alpha version 0.6 (15/06/18)

"""
#Script to identify SNPs located within a genetic codon:

Annotation of SNPs (including synonymous and non-synonymous variants) located within the same genetic codon requires attention. This idea conflicts with the annotations observed by traditional annotation software widely used now a days. While looking at the combined effect within the framework of genetic codon, we have new / altered codons that code for new amino acid that can be predicted by using this MACARON python script.

#####------How to run quickly run MACARON-------

# python MACARON -i yourinputfile.vcf

# python MACARON -i /path/to/yourinputfile.vcf -o /path/to/MACARON_output.txt -f (INFO)_FIELD_HEADER --GATK /path/to/GenomeAnalysisTK.jar --HG_REF /path/to/hg.fasta --SNPEFF /path/to/snpeff.jar --SNPEFF_HG hg19 

# python MACARON -i yourinputfile.vcf -o MACARON_output.txt --gatk4 (when using gatk versions >= 4.0)

"""

import sys, os, time
import itertools
import multiprocessing
import re
import subprocess
from argparse import ArgumentParser

## GLOBAL VARIABLES (IMPORTANT: You can set the default values here)
#GATK="/home/wuk/software/GenomeAnalysisTK.jar"
GATK="/home/wuk/software/gatk-4.0.1.2/gatk-package-4.0.1.2-local.jar"
HG_REF="/home/wuk/Working/gnme_refrnces/Homo_sapiens_assembly19.fasta"
SNPEFF="/home/wuk/software/snpEff/snpEff.jar"
SNPEFF_HG="GRCh37.75" ## SnpEff human genome annotation database version

## PRINTINGS, AESTHETICS
str_progress_list = ["\tIndexing VCF file",
		     "\tIdentifying SnpClusters",
		     "\tExtracting SnpClusters",
		     "\tAnnotating SnpClusters",
		     "\tExcluding InDels",
		     "\tGenerating a SnpCluster Table",
		     "\tRe-annotating Codons",
		     "\tRemoving SnpCluster if AA_Change_pcSNV == (AA1 or AA2)",
		     "\tExtracting established SnpClusters - for which two (or three) SNPs are reference heterozygous (or non-reference homozygous)"]
header = ("\n" +
	  "(###############)\n" +
	  "@@@@ MACARON @@@@\n" +
	  "(###############)\n\n" +
	  "Starting....\n")
footer = "\nMACARON Run Completed. Bon Courage with Analysis ...,,,!!!!\n"

def animate(keep_anim, idx):
	c = "|/-\\"
	i=0
	while keep_anim.is_set():
		i += 1
		sys.stdout.write("\r>" + str_progress_list[idx] + ": " + c[i % len(c)] +"\r")
		sys.stdout.flush()
		time.sleep(0.1)
	sys.stdout.write("\r " + str_progress_list[idx] + ": Done!\n")
def print_step(keep_anim, idx):
	if not(ECO):
		keep_anim.set()
		anim_thread = multiprocessing.Process(target=animate, args=(keep_anim, idx))
		anim_thread.start()
		return anim_thread
	else:
		sys.stdout.write(">" + str_progress_list[idx] + ": in progress...\r" )
		sys.stdout.flush()
		return None
def end_print_step(keep_anim, anim_thread, idx):
	if not(ECO):
		keep_anim.clear(); anim_thread.join()
	else:
		sys.stdout.write(" " + str_progress_list[idx] + ": Done!	 \n")
		sys.stdout.flush()


## CLASS DEFINITIONS 
		
class fileHandler:
	def __init__(self):
		self.data = []; 
	def open_file(self, readfl):
		self.rfile = open(readfl, 'r').readlines()
		return self.rfile
	def write_file(self, writefl):
		self.wfile = open(writefl, 'w')
		return self.wfile

class SearchDB(fileHandler):
	def __init__(self):
		self.data = [] 
		from collections import defaultdict
		self.ident_ranges_HMBM = defaultdict(list)
	def Search_CODON(self, vcf_input, workdir, FIELDS):
		"""
		Calling SNPClusters
		
		USAGE INSTRUCTIONS:	Full path to the software directories should be set before compiling.
		"""

		#######################
		Fld_Len = int(len(FIELDS.split(",")))
		FldwithF = " ".join(["-F " + str(x) for x in FIELDS.split(",")])

		## Options compatible with GATK versions >= 4.0: add option --gatk4 when calling MACARON
		if GATK4:
			GATK_v, SNPeff_v = (" ", " -v ") if VERBOSE else (" --QUIET true --verbosity ERROR ", " ")
			cmd0 = "java -Xmx12g -jar " + GATK + " IndexFeatureFile -F " + vcf_input + GATK_v
			cmd1 = "java -Xmx12g -jar " + GATK + " VariantFiltration -R " + HG_REF + " -V "+ vcf_input +" -O " + workdir + "snp_clsters2_ws3.vcf --cluster-size 2 --cluster-window-size 3" + GATK_v
			cmd2 = "java -Xmx12g -jar " + GATK + " SelectVariants -R " + HG_REF + " -V " + workdir + "snp_clsters2_ws3.vcf -O " + workdir + "snp_clsters2_ws3_clstronly.vcf -select 'FILTER == SnpCluster'" + GATK_v
			cmd3 = "java -Xmx12g -jar " + SNPEFF + SNPeff_v + SNPEFF_HG + " -formatEff -lof -classic " + workdir + "snp_clsters2_ws3_clstronly.vcf > " + workdir + "snp_clsters2_ws3_clstronly_annt.vcf"
			cmd4 = "java -Xmx12g -jar " + GATK + " SelectVariants -R " + HG_REF + " -V " + workdir + "snp_clsters2_ws3_clstronly_annt.vcf -O " + workdir + "snp_clsters2_ws3_clstronly_annt_snv.vcf --select-type-to-exclude INDEL" + GATK_v
			cmd5 = "java -Xmx12g -jar " + GATK + " VariantsToTable -R " + HG_REF + " -V " + workdir + "snp_clsters2_ws3_clstronly_annt_snv.vcf -F CHROM -F POS -F ID -F REF -F ALT -F EFF " + FldwithF + " -GF GT --error-if-missing-data --show-filtered -O " + workdir + "snp_clsters2_ws3_clstronly_annt_snv_clstronly.table" + GATK_v
			## GATK4 needs to index the VCF upfront
			thread = print_step(keep_anim, 0)
			subprocess.check_output(cmd0, shell=True)
			end_print_step(keep_anim, thread, 0)

		else: ## Options comptatible with GATK versions < 4
			GATK_v, SNPeff_v = (" ", " -v ") if VERBOSE else (" --logging_level ERROR ", " ")
			cmd1 = "java -Xmx4g -jar " + GATK + " -T VariantFiltration -R " + HG_REF + " -V " + vcf_input + " -o " + workdir + "snp_clsters2_ws3.vcf --clusterSize 2 --clusterWindowSize 3" + GATK_v
			cmd2 = "java -Xmx4g -jar " + GATK + " -T SelectVariants -R " + HG_REF + " -V " + workdir + "snp_clsters2_ws3.vcf -o " + workdir + "snp_clsters2_ws3_clstronly.vcf -select 'FILTER == SnpCluster'" + GATK_v
			cmd3 = "java -Xmx4g -jar " + SNPEFF + SNPeff_v + SNPEFF_HG + " -formatEff -lof -classic " + workdir + "snp_clsters2_ws3_clstronly.vcf > " + workdir + "snp_clsters2_ws3_clstronly_annt.vcf"
			cmd4 = "java -Xmx4g -jar " + GATK + " -T SelectVariants -R " + HG_REF + " -V " + workdir + "snp_clsters2_ws3_clstronly_annt.vcf -o " + workdir + "snp_clsters2_ws3_clstronly_annt_snv.vcf --selectTypeToExclude INDEL" + GATK_v
			cmd5 = "java -Xmx4g -jar " + GATK + " -T VariantsToTable -R " + HG_REF + " -V " + workdir + "snp_clsters2_ws3_clstronly_annt_snv.vcf -F CHROM -F POS -F ID -F REF -F ALT -F EFF " + FldwithF + " -GF GT --showFiltered -o " + workdir + "snp_clsters2_ws3_clstronly_annt_snv_clstronly.table" + GATK_v
		
		thread = print_step(keep_anim, 1)
		subprocess.check_output(cmd1, shell=True)
		end_print_step(keep_anim, thread, 1)
		thread = print_step(keep_anim, 2)
		subprocess.check_output(cmd2, shell=True)
		end_print_step(keep_anim, thread, 2)
		thread = print_step(keep_anim, 3)
		subprocess.check_output(cmd3, shell=True)
		end_print_step(keep_anim, thread, 3)
		thread = print_step(keep_anim, 4)
		subprocess.check_output(cmd4, shell=True)
		end_print_step(keep_anim, thread, 4)
		thread = print_step(keep_anim, 5)
		subprocess.check_output(cmd5, shell=True)
		end_print_step(keep_anim, thread, 5)

		subprocess.check_output("rm snpEff_genes.txt", shell=True)
		subprocess.check_output("rm snpEff_summary.html", shell=True)
                
		####--------------------------------------

		def Change_zygo(ref, alt, zyg):
			"""
			program to convert zygosity code ref/alt to 0/1.
			Input Variables:
			
			ref = "A";alt = "G"
			zyg =  ['A/G', 'G/A', 'G/G', 'A/A', './.', 'G/.', './G', './A', 'A/.']
			chzyg = ['0/1', '1/0', '1/1', '0/0', './.', '1/.', './1', './0', '0/.']
			InputUsage Variables:
			chgz = Change_zygo(ref, alt, zyg)

			"""
			import re
			chg_zyg = {}; i = 1
			for cs in zyg:
				csp = re.split('[|/]', cs)
				if ((ref == csp[0]) and (ref == csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i] = "0/0"
				elif ((ref != csp[0]) and (ref == csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i] = "1/0"
				elif ((ref == csp[0]) and (ref != csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i] = "0/1"
				elif ((ref != csp[0]) and (ref != csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i] = "1/1"
				elif ((csp[0] == ".") and (csp[1] == ".")):
					chg_zyg[i] = cs
				elif ((csp[1] == ".")):
					if (ref == csp[0]):
						chg_zyg[i] = "0/."
					elif (ref != csp[0]):
						chg_zyg[i] = "1/."
				elif ((csp[0] == ".")):
					if (ref == csp[1]):
						chg_zyg[i] = "./0"
					elif (ref != csp[1]):
						chg_zyg[i] = "./1"
				i += 1
			return list(chg_zyg.values())
		####--------------------------------------

		with open(workdir + "snp_clsters2_ws3_clstronly_annt_snv_clstronly.table", 'r') as f1, open(workdir + "temp_file1", 'w') as output:
			first_line = f1.readline().strip()
			zyg_head = '\t'.join(first_line.split()[6+Fld_Len:])
			output.write(first_line + "\t" + zyg_head + "\t" + str("Protein_coding_EFF	AA-Change	REF-codon	ALT-codon") + "\n")
			for line in f1:
				line1 = line.strip()
				line_TAB = line1.split("\t")
				line_EFF = line_TAB[5].split("|")
				if (len(line_EFF) > 1):
					if ((line_EFF[1] == "SILENT") or (line_EFF[1] == "MISSENSE") or (line_EFF[1] == "NONSENSE")):
						True
						linesp = line_EFF[2].split("/")
						ref = line_TAB[3]
						alt = line_TAB[4]	
						zyg = line_TAB[6+Fld_Len:]
						chgz = Change_zygo(ref, alt, zyg)
						chgz_out = '\t'.join(chgz)
						wrt = str(line_EFF[1] + "\t" + line_EFF[3] + "\t" + linesp[0] + "\t" + linesp[1])
						output.write(line1 + "\t" + chgz_out + "\t" + wrt + "\n")
		return None

	##-------------------------------
	import string
	gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T','AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R','CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P','CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R','GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A','GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G','TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L','TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
	
	# a function to translate a single codon
	def translate_codon(self, codon):
		return self.gencode.get(codon.upper(), '#')
		# a function to split a sequence into codons
	def split_into_codons(self, dna, frame):
		codons = []
		for i in range(frame-1, len(dna)-2, 3):
			codon = dna[i:i+3]
			codons.append(codon)
		return codons
	# a function to translate a dna sequence in a single frame
	def translate_dna_single(self, dna, frame=1):
		codons = self.split_into_codons(dna, frame)
		amino_acids = ''
		for codon in codons:
			amino_acids = amino_acids + self.translate_codon(codon)
		return amino_acids
	
	def TWO_VAR(self, workdir):
		"""
		###---Two variants (2VAR) codon changes---
		##----------------------------------
		"""
		lines =	 open(workdir + "temp_file1", "r").read().splitlines()
		writ2 = self.write_file(workdir + "temp_file2")
		#---
		midline_head =	lines[0].strip().split("\t")
		try:
			midline_headcrp = '\t'.join([w.replace(midline_head[5], 'Gene_Name') for w in midline_head])
		except ValueError:
			pass
		writ2.write(str(midline_headcrp + "\t" + "ALT-codon_merge-2VAR" + "\t" + "AA-Change-2VAR") + "\n")
		#---

		i=0; TRcode =""; protcode = ""; midline_crpit = ""
		for i in range(len(lines)):
			#---
			midline_crp =  lines[i].strip().split("\t")
			try:
				if len(midline_crp[5].split("|")) != 1:
					GeneName = midline_crp[5].split("|")[5]
					midline_crpit = '\t'.join([w.replace(midline_crp[5], GeneName) for w in midline_crp])
			except ValueError:
				pass
			#---
			try:
				beforeline = lines[i-1].strip()
				line0 = beforeline.split("\t")[-3]
				beforeline1 = re.findall("\d+", line0)
				
				midline = lines[i].strip()
				line1 = midline.split("\t")[-3]
				midline1 = re.findall("\d+", line1)
				
				nextline = lines[i+1].strip()
				line2 = nextline.split("\t")[-3]
				nextline1 = re.findall("\d+", line2)
				#----Condition to replace empty list ([] to ['000']) produced from the header column "AA-Change". 
				if (line0 == "AA-Change"):
					beforeline1.append('000')
				elif (line1 == "AA-Change"):
					midline1.append('000')
				elif (line2 == "AA-Change"):
					nextline1.append('000')
				#----
				REFbf=[]; line22=""
				if ((beforeline1[0] == midline1[0]) or (midline1[0] == nextline1[0])):
					spREFbf=[]; lscod1=[]; lscod2=[]
					REF= lines[i].strip().split("\t")[-2]
					if (midline1[0] == beforeline1[0]):
						REFbf = lines[i-1].strip().split("\t")[-2]
						
					line11 = lines[i].strip().split("\t")[-1]
					if (midline1[0] == nextline1[0]):
						line22 = lines[i+1].strip().split("\t")[-1]
					if REFbf != []:
						for cod in REFbf:
							spREFbf.append(cod)
					for cod in line11:
						lscod1.append(cod)
					for cod in line22:
						lscod2.append(cod)

					if (((lscod1[0].islower()==True and lscod2[0].islower()==True) or (lscod1[1].islower()==True and lscod2[1].islower()==True) or (lscod1[2].islower()==True and lscod2[2].islower()==True)) and ((lscod1[0] == lscod2[0]) or (lscod1[1] == lscod2[1]) or (lscod1[2] == lscod2[2]))):
						threeltr_code = []
						if ((lscod1[0].isupper()==True and lscod2[0].islower()==True) or (lscod1[0] == lscod2[0])):
							threeltr_code.append(lscod1[0])
						elif((lscod1[0].islower()==True and lscod2[0].isupper()==True) or (lscod1[0] == lscod2[0])):
							threeltr_code.append(lscod2[0])
	
						if((lscod1[1].isupper()==True and lscod2[1].islower()==True) or (lscod1[1] == lscod2[1])):
							threeltr_code.append(lscod1[1])
						elif((lscod1[1].islower()==True and lscod2[1].isupper()==True) or (lscod1[1] == lscod2[1])):
							threeltr_code.append(lscod2[1])

						if((lscod1[2].isupper()==True and lscod2[2].islower()==True) or (lscod1[2] == lscod2[2])):
							threeltr_code.append(lscod1[2])
						elif((lscod1[2].islower()==True and lscod2[2].isupper()==True) or (lscod1[2] == lscod2[2])):
							threeltr_code.append(lscod2[2])
						#-----------------------
						if(len(threeltr_code)==3):
							TRcode = ''.join(threeltr_code)
						else:
							TRcode = 'Multiallelic-t1'
						#-----------------------
					else:
						TRcode = '...'

					if(TRcode != line22):
						protcode = self.translate_dna_single(TRcode)
					else:				
						TRcode = "."

					if (REF == REFbf):
						protcode = "#####";TRcode = "Multiallelic-t2"	
					writ2.write(str(midline_crpit + "\t" + TRcode + "\t" + protcode) + "\n")

			except IndexError:
				if not midline_crpit.split("\t")[0] == "CHROM":
					writ2.write(str(midline_crpit + "\t" + "." + "\t" + protcode) + "\n")
				pass
		return None

	def THREE_VAR(self, workdir):
		"""
		###---Three variants (3VAR) codon changes---
		##----------------------------------
		"""
		lines =	open(workdir + "temp_file2", "r").read().splitlines()
		writ3 = self.write_file(workdir + "temp_file3")
		#---
		midline_head =lines[0].strip()
		writ3.write(str(midline_head + "\t" + "ALT-codon_merge-3VAR" + "\t" + "AA-Change-3VAR") + "\n")
		#---
		i = 0; TRcode = ""; protcode = ""
		for i in range(len(lines)):
			try:
				#---
				midline_crp = lines[i].strip()
				#---
				beforeline = lines[i-1].strip()
				line0 = beforeline.split("\t")[-5]
				beforeline1 = re.findall("\d+", line0)
				midline = lines[i].strip()
				line1 = midline.split("\t")[-5]
				midline1 = re.findall("\d+", line1)
				nextline = lines[i+1].strip()
				line2 = nextline.split("\t")[-5]
				nextline1 = re.findall("\d+", line2)

				line11=[]; line22=[]; writelst= []
				if ((beforeline1[0] == midline1[0]) or (midline1[0] == nextline1[0])):
					lscod1=[]; lscod2=[]
					line11 = lines[i].strip().split("\t")[-2]
					if (midline1[0] == nextline1[0]):
						line22 = lines[i+1].strip().split("\t")[-2]
					for cod in line11:
						lscod1.append(cod)
					for cod in line22:
						lscod2.append(cod)
					#-----------------------
					if(len(lscod1) == 3 and len(lscod2) == 3):
						threeltr_code = [];

						if ((lscod1[0].isupper()==True) and (lscod2[0].islower()==True) or
						    (lscod1[0].isupper()==True and lscod2[0].isupper()==True)):
							threeltr_code.append(lscod1[0])
						elif ((lscod1[0].islower()==True) and (lscod2[0].isupper()==True)):
							threeltr_code.append(lscod2[0])

						if ((lscod1[1].isupper()==True) and (lscod2[1].islower()==True) or
						    (lscod1[1].isupper()==True and lscod2[1].isupper()==True)):
							threeltr_code.append(lscod1[1])
						elif ((lscod1[1].islower()==True) and (lscod2[1].isupper()==True)):
							threeltr_code.append(lscod2[1])
			
						if ((lscod1[2].isupper()==True) and (lscod2[2].islower()==True) or
						    (lscod1[2].isupper()==True and lscod2[2].isupper()==True)):
							threeltr_code.append(lscod1[2])
						elif ((lscod1[2].islower()==True) and (lscod2[2].isupper()==True)):
							threeltr_code.append(lscod2[2])
						#-----------------------
						if(len(threeltr_code) == 3):
							TRcode = ''.join(threeltr_code)
						else:
							TRcode = 'Multiallelic-t1'
						#-----------------------
					else:
						TRcode = '.'

					if(TRcode != line22):
						protcode = self.translate_dna_single(TRcode)
						if len(protcode) == 0:
							protcode = "."; TRcode = "."
					else:
						TRcode = "."

					if protcode == "" or TRcode == "":
						protcode = "."; TRcode = "."
					#-----------------------
					writ3.write(str(midline_crp + "\t" + TRcode + "\t" + protcode) + "\n")

					##----------------------
			except IndexError:
				if not midline_crp.split("\t")[0] == "CHROM":
					if len(protcode) ==0:
						protcode = "."; TRcode = "."
						writ3.write(str(midline_crp + "\t" + TRcode + "\t" + protcode) + "\n")
					else:
						protcode = "."; TRcode = "."
						writ3.write(str(midline_crp + "\t" + TRcode + "\t" + protcode) + "\n")
				pass
		return None

	def PARS_OUT_VAR(self, workdir):
		"""
		###---Pars variants (2VAR _ 3VAR) based on change of protein codons ---
		##----------------------------------
		"""
		#This first records the ones that do match and saves the value in column 3 (less the final matching amino acid residue).
		subprocess.check_output("awk 'index($(NF-6), $(NF-2)) {print}' " + workdir + "temp_file3 | awk '{print $(NF-6)}' | sed s'/.$//' > " + workdir + "matches.list", shell=True)

		#This then searches for the ones that do not match, and then also eliminates any value recorded in matches.list
		subprocess.check_output("awk '!index($(NF-6), $(NF-2)) { print }' " + workdir + "temp_file3 | grep -w -v -f " + workdir + "matches.list | awk -F'\t' '{ print }' > " + workdir + "temp_file4", shell=True)

		return None

	def ZYGO_PAIR(self, macaron_output, workdir, FIELDS):
		"""
		###---Pair of zygosity check: remove pair of codons for which one SNP is reference homozygous ---
		##----------------------------------
		"""

		Fld_Len = int(len(FIELDS.split(",")))
		lines =	open(workdir + "temp_file4", "r").read().splitlines()
		writ4 = self.write_file(macaron_output)
		#---
		midline_head = lines[0].strip()
		writ4.write(str(midline_head) + "\n")
		#---
		nextline_pos_lst = []

		for i in range(len(lines)):
			#---
			midline_crp = lines[i].strip()
			#---
			beforeline = lines[i-1].strip()
			beforeline_pos = re.findall("\d+", beforeline.split("\t")[-7])
			midline = lines[i].strip()
			midline_pos = re.findall("\d+", midline.split("\t")[-7])
			try:
				nextline = lines[i+1].strip()
				nextline_pos =re.findall("\d+", nextline.split("\t")[-7])
				##-------
				if (midline_pos[0] == nextline_pos[0]):
					Start_lns = int(6) + Fld_Len
					End_lns = int(8)
					Start_zyg_lns = int(len(midline.split("\t")) - int(Start_lns))
					midline_zyg = midline.split("\t")[-int(Start_zyg_lns):-int(End_lns)]
					nextline_zyg = nextline.split("\t")[-int(Start_zyg_lns):-int(End_lns)]
					checkList = list(['0/1:1/0', '1/0:0/1', '1/0:1/0', '0/1:0/1', '1/0:1/1', '1/1:1/0', '0/1:1/1', '1/1:0/1', '1/1:1/1'])
					Mergetwozyg = ','.join([str(a) + ":" + b for a,b in zip(midline_zyg, nextline_zyg)])
					Mergetwozygsp = Mergetwozyg.split(",")
					if set(Mergetwozygsp).intersection(checkList) != set([]):
						if ((midline_pos[0] == nextline_pos[0]) and (midline_pos[0] == beforeline_pos[0])):
							writ4.write(str(nextline)+"\n")
						else:
							writ4.write(str(midline+"\n"+nextline)+"\n")
				##-------
			except IndexError:
				#Condition to remove the duplicated line at the end and prints only paired lines.
				nextline_pos_lst.append(nextline_pos[0])
				result = dict((i, nextline_pos_lst.count(i)) for i in nextline_pos_lst)
				if len(result) == 1:
					writ4.write(str(nextline) + "\n")

		return None
		###------


## MACARON MAIN
if __name__ == "__main__":
	## Parsing arguments
	parser = ArgumentParser(description="-Script to identify SnpClusters (SNPs within the same genetic codon)")
	parser.add_argument("-i", "--infile", dest="INPUTFile",default=False, required=True, help="Full path of the input VCF file.")
	parser.add_argument("-o", "--outfile", dest="OUTPUTFile",default="./MACARON_output.txt", required=False, help="Path of the output txt file (Default Output file: MACARON_output.txt)")
	parser.add_argument("-f", "--fields", dest="Fields", default="QUAL", required=False, help=" Single field name or comma-seperated ',' multiple field names can be given. Field name should be given according to the (INFO) field header of the input vcf file. Example: -f Func.refGene,ExonicFunc.refGene,Gene.refGene,1000g2015aug_all,ExAC_ALL,ExAC_EAS,clinvar_20161128,gnomAD_exome_ALL,gnomAD_genome_ALL,EFF,CSQ")
	parser.add_argument("--GATK", dest="GATK_path", default=GATK, required=False, help="Indicate the full path to GATK jar file")
	parser.add_argument("--HG_REF", dest="HG_REF_path", default=HG_REF, required=False, help="Indicate the full path to the reference genome fasta file")
	parser.add_argument("--SNPEFF", dest="SNPEFF_path", default=SNPEFF, required=False, help="Indicate the full path to SnpEff jar file")
	parser.add_argument("--SNPEFF_HG", dest="SNPEFF_HG_version", default=SNPEFF_HG, required=False, help="Indicate SnpEff human genome annotation database version")
	parser.add_argument("--gatk4", dest="GATK4", default=False, required=False, action='store_true', help="Add this option when using GATK versions >= 4.0")
	parser.add_argument("-v", "--verbosity", dest="Verbosity", default=False, required=False, action='store_true', help="Use to print verbosity (Mostly GATK/SNPEFF output)")
	parser.add_argument("-c", "--eco_friendly", dest="ECO", default=False, required=False, action='store_true', help="Save a thread, but you won't be able to stare at the fabulous animation while waiting ...")

	## Assign arguments to global variables
	args = parser.parse_args()
	FIELDS = args.Fields
	GATK4 = args.GATK4
	VERBOSE = args.Verbosity
	GATK = args.GATK_path
	HG_REF = args.HG_REF_path
	SNPEFF = args.SNPEFF_path
	SNPEFF_HG = args.SNPEFF_HG_version
	ECO = args.ECO

	## Inputs / Outputs path & names
	INF = args.INPUTFile
	OUTF = args.OUTPUTFile
	TMPDIR = os.path.dirname(os.path.abspath(OUTF)) + "/macaron_tmp/"
	subprocess.check_output("mkdir -p " + TMPDIR, shell=True)

	########################
	##    MAIN PROCESS    ##
	########################
	print(header)
	
	## Check if global variables point to existing files:
	INF_check = os.path.exists(INF)
	GATK_check = os.path.exists(GATK)
	HG_REF_check = os.path.exists(HG_REF)
	SNPEFF_check = os.path.exists(SNPEFF)
	SNPEFF_HG_non_empty = (SNPEFF_HG != "")
	if not(INF_check and GATK_check and HG_REF_check and SNPEFF_check and SNPEFF_HG_non_empty):
		print(">ERROR : One or several global variable: \n  VCF={} > {}\n  GATK={} > {}\n  HG_REF={} > {}\n  SNPEFF={} > {}\n  SNPEFF_HG={} > {}\n>Please correct and try again!".format(INF, INF_check, GATK, GATK_check, HG_REF, HG_REF_check, SNPEFF, SNPEFF_check, SNPEFF_HG, SNPEFF_HG_non_empty))
		sys.exit(1)
	else: ## If everything checks out, start MACARON
		## Animation event initialization
		keep_anim = multiprocessing.Event()

		## 1)VARIANTS FILTERING, ANNOTATION (GATK,SNPEff)
		clF1 = SearchDB().Search_CODON(INF, TMPDIR, FIELDS)
		## 2)SEARCH MULTI-SNPS CODONS
		thread = print_step(keep_anim, 6)
		clF2 = SearchDB().TWO_VAR(TMPDIR)
		clF3 = SearchDB().THREE_VAR(TMPDIR)
		end_print_step(keep_anim, thread, 6)
		## 3)CHECK IF SNPCLUSTERS IMPACT CODON
		thread = print_step(keep_anim, 7)
		clF4 = SearchDB().PARS_OUT_VAR(TMPDIR)
		end_print_step(keep_anim, thread, 7)
		## 4) EXTRACT SNPCLUSTERS (Keeping SnpCluster if >=1 sample is Ref-Heterozygous or nonRef-Homozygous)
		thread = print_step(keep_anim, 8)
		clF5 = SearchDB().ZYGO_PAIR(OUTF, TMPDIR, FIELDS)
		end_print_step(keep_anim, thread, 8)
		print(footer)
	subprocess.check_output("rm -r " + TMPDIR, shell=True)
