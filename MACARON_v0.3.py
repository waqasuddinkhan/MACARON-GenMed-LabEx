#!/usr/bin/python

#alpha version 0.3

"""
#Script to identify SNPs located within a genetic codon:

Annotation of SNPs (including synonymous and non-synonymous variants) located within the same genetic codon requires attention. This idea conflicts with the annotations observed by traditional annotation software widely used now a days. While looking at the combined effect within the framework of genetic codon, we have new / altered codons that code for new amino acid that can be predicted by using this MACARON python script.

#####------Inputs-------
# python MACARON_v0.3.py INPUT=Inputfile(VCF) OUTPUT=Outputfile(TXT) DIR-PATH=Full_Path_of_Working_Directory FIELDS=Required-Fields

# python MACARON_v0.3.py -i yourinputfile.vcf -o MACARON_output.txt -d ./Example/ -f INFO_HEADER

# python MACARON_v0.3.py -i input_test0.vcf -o MACARON_output.txt -d /Users/varma/Desktop/TEST_MACARON/Example/ -f INFO_HEADER

"""

import sys
import re
import commands
import subprocess
from argparse import ArgumentParser

class fileHandler:
	def __init__(self):
		self.data = []; 
		#print "Calling fileHandler constructor"
	def open_file(self,readfl):
		self.rfile = open(readfl,'r').readlines()
		return self.rfile
	def write_file(self,writefl):
		self.wfile = open(writefl,'w')
		return self.wfile

class SearchDB(fileHandler):

	def __init__(self):
		self.data = []; 
		from collections import defaultdict
		self.ident_ranges_HMBM = defaultdict(list)
	def Search_CODON(self,readfl1,workdir,FIELDS):
		"""
		Calling SNPClusters

		USAGE INSTRUCTIONS:	Full path to the software directories should be set before compiling.
		"""
		GATK="/home/wuk/software/GenomeAnalysisTK.jar"
		HG_REF="/home/wuk/Working/gnme_refrnces/Homo_sapiens_assembly19.fasta"
		SNPEFF="/home/wuk/software/snpEff/snpEff.jar"
		SNPEFF_HG = "GRCh37.75"

		#######################
		Fld_Len = int(len(FIELDS.split(",")))
		FldwithF = " ".join(["-F "+str(x) for x in FIELDS.split(",")])

		cmd1 = "java -Xmx4g -jar "+GATK+" -T VariantFiltration -R "+HG_REF+" -V "+workdir+readfl1+" -o "+workdir+"snp_clsters2_ws3.vcf --clusterSize 2 --clusterWindowSize 3"
		
		cmd2 = "java -Xmx4g -jar "+GATK+" -T SelectVariants -R "+HG_REF+" -V "+workdir+"snp_clsters2_ws3.vcf -o "+workdir+"snp_clsters2_ws3_clstronly.vcf -select 'FILTER == SnpCluster'"
		
		cmd3 = "java -Xmx4g -jar "+SNPEFF+" -v "+SNPEFF_HG+" -formatEff -lof -classic "+workdir+"snp_clsters2_ws3_clstronly.vcf > "+workdir+"snp_clsters2_ws3_clstronly_annt.vcf"
		
		cmd4 = "java -Xmx4g -jar "+GATK+" -T SelectVariants -R "+HG_REF+" -V "+workdir+"snp_clsters2_ws3_clstronly_annt.vcf -o "+workdir+"snp_clsters2_ws3_clstronly_annt_snv.vcf --selectTypeToExclude INDEL"
		
		cmd5 = "java -Xmx4g -jar "+GATK+" -T VariantsToTable -R "+HG_REF+" -V "+workdir+"snp_clsters2_ws3_clstronly_annt_snv.vcf -F CHROM -F POS -F ID -F REF -F ALT -F EFF "+FldwithF+" -GF GT --allowMissingData --showFiltered -o "+workdir+"snp_clsters2_ws3_clstronly_annt_snv_clstronly.table"
		
		subprocess.check_output(cmd1, shell=True)
		subprocess.check_output(cmd2, shell=True)
		subprocess.check_output(cmd3, shell=True)
		subprocess.check_output(cmd4, shell=True)
		subprocess.check_output(cmd5, shell=True)
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
			chg_zyg = {};i=1
			for cs in zyg:
				csp = re.split('[|/]',cs)
				if ((ref == csp[0]) and (ref == csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i]= "0/0"
				elif ((ref != csp[0]) and (ref == csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i]= "1/0"
				elif ((ref == csp[0]) and (ref != csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i]= "0/1"
				elif ((ref != csp[0]) and (ref != csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i]= "1/1"
				elif ((csp[0] == ".") and (csp[1] == ".")):
					chg_zyg[i]= cs
				elif ((csp[1] == ".")):
					if (ref == csp[0]):
						chg_zyg[i]= "0/."
					elif (ref != csp[0]):
						chg_zyg[i]= "1/."
				elif ((csp[0] == ".")):
					if (ref == csp[1]):
						chg_zyg[i]= "./0"
					elif (ref != csp[1]):
						chg_zyg[i]= "./1"
				i+=1
			return list(chg_zyg.values())
		####--------------------------------------

		with open(workdir+"snp_clsters2_ws3_clstronly_annt_snv_clstronly.table",'r') as f1, open(workdir+"temp_file1",'w') as output:
			first_line = f1.readline().strip()
			zyg_head = '\t'.join(first_line.split()[6+Fld_Len:])
			output.write(first_line+"\t"+zyg_head+"\t"+str("Protein_coding_EFF	AA-Change	Ref-codon	Alt-codon")+"\n")
			for line in f1:
				line1 = line.strip()
				line_TAB = line1.split("\t")
				line_EFF = line_TAB[5].split("|")
				if ((line_EFF[1] == "SILENT") or (line_EFF[1] == "MISSENSE") or (line_EFF[1] == "NONSENSE")):
					True
					linesp = line_EFF[2].split("/")
					ref = line_TAB[3]; alt = line_TAB[4]	
					zyg = line_TAB[6+Fld_Len:]
					chgz = Change_zygo(ref, alt, zyg)
					chgz_out = '\t'.join(chgz)
					wrt = str(line_EFF[1]+"\t"+line_EFF[3]+"\t"+linesp[0]+"\t"+linesp[1])
					output.write(line1+"\t"+chgz_out+"\t"+wrt+"\n")
		print "Done with annotation..."
		return None

	##-------------------------------
	import string
	gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T','AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R','CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P','CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R','GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A','GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G','TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L','TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
	
	# a function to translate a single codon
	def translate_codon(self,codon):
		return self.gencode.get(codon.upper(), '#')
		# a function to split a sequence into codons
	def split_into_codons(self,dna,frame):
		codons = []
		for i in range(frame - 1, len(dna)-2, 3):
			codon = dna[i:i+3]
			codons.append(codon)
		return codons
	# a function to translate a dna sequence in a single frame
	def translate_dna_single(self,dna, frame=1):
		codons = self.split_into_codons(dna, frame)
		amino_acids = ''
		for codon in codons:
			amino_acids = amino_acids + self.translate_codon(codon)
		return amino_acids
	
	def TWO_VAR(self,workdir):
		"""
		###---Two variants(2VAR) codon changes---
		##----------------------------------
		"""
		lines =  open(workdir+"temp_file1","r").read().splitlines()
		writ2 = self.write_file(workdir+"temp_file2")
		#---
		midline_head =  lines[0].strip().split("\t")
		try:
			midline_headcrp = '\t'.join([w.replace(midline_head[5], 'Gene_Name') for w in midline_head])
		except ValueError:
			pass
		writ2.write(str(midline_headcrp+"\t"+"Altcodon_merge-2VAR"+"\t"+"AA-change-2VAR")+"\n")
		#---

		i=0; TRcode =""; protcode = ""; midline_crpit = ""
		for i in range(len(lines)):
			True
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
				beforeline1 =  re.findall("\d+", line0)
				
				midline = lines[i].strip()
				line1 = midline.split("\t")[-3]
				midline1 =  re.findall("\d+", line1)
				
				nextline = lines[i+1].strip()
				line2 = nextline.split("\t")[-3]
				nextline1 =  re.findall("\d+", line2)
				#----Condition to replace empty list ([] to ['000']) produced from the header column "AA-Change". 
				if (line0 == "AA-Change"):
					beforeline1.append('000')
				elif (line1 == "AA-Change"):
					midline1.append('000')
				elif (line2 == "AA-Change"):
					nextline1.append('000')
				#----
				REFbf=[];
				if ((beforeline1[0] == midline1[0]) or (midline1[0] == nextline1[0])):
					True 
					spREFbf=[];lscod1=[];lscod2=[]
					REF= lines[i].strip().split("\t")[-2]
					if (midline1[0] == beforeline1[0]):
						REFbf = lines[i-1].strip().split("\t")[-2]
						
					line11 = lines[i].strip().split("\t")[-1]
					if (midline1[0] == nextline1[0]):
						line22 = lines[i+1].strip().split("\t")[-1]
					if REFbf!=[]:
						for cod in REFbf:
							spREFbf.append(cod)
					
					for cod in line11:
						lscod1.append(cod)
					
					try:
						True
						for cod in line22:
							lscod2.append(cod)
					except:
						False
					if (((lscod1[0].islower()==True and lscod2[0].islower()==True) or (lscod1[1].islower()==True and lscod2[1].islower()==True) or (lscod1[2].islower()==True and lscod2[2].islower()==True)) and ((lscod1[0] == lscod2[0]) or (lscod1[1] == lscod2[1]) or (lscod1[2] == lscod2[2]))):
						True
						threeltr_code =[]
						if ((lscod1[0].isupper()==True and lscod2[0].islower()==True) or (lscod1[0] == lscod2[0])):
							True
							threeltr_code.append(lscod1[0])
						elif((lscod1[0].islower()==True and lscod2[0].isupper()==True) or (lscod1[0] == lscod2[0])):
							True
							threeltr_code.append(lscod2[0])
	
						if((lscod1[1].isupper()==True and lscod2[1].islower()==True) or (lscod1[1] == lscod2[1])):
							True
							threeltr_code.append(lscod1[1])
						elif((lscod1[1].islower()==True and lscod2[1].isupper()==True) or (lscod1[1] == lscod2[1])):
							True
							threeltr_code.append(lscod2[1])

						if((lscod1[2].isupper()==True and lscod2[2].islower()==True) or (lscod1[2] == lscod2[2])):
							True
							threeltr_code.append(lscod1[2])
						elif((lscod1[2].islower()==True and lscod2[2].isupper()==True) or (lscod1[2] == lscod2[2])):
							True
							threeltr_code.append(lscod2[2])
						#-----------------------
						if(len(threeltr_code)==3):
							True
							TRcode = ''.join(threeltr_code)
						else:
							True
							TRcode = 'Multiallelic-t1'
						#-----------------------
					else:
						True
						TRcode = '...'

					if(TRcode != line22):
						True
						protcode = self.translate_dna_single(TRcode)
					else:				
						True
						TRcode = ".";

					if (REF == REFbf):
						True
						protcode = "#####";TRcode = "Multiallelic-t2"	
					writ2.write(str(midline_crpit+"\t"+TRcode+"\t"+protcode)+"\n")

			except IndexError:
				if not midline_crpit.split("\t")[0] == "CHROM":
					writ2.write(str(midline_crpit+"\t"+"."+"\t"+protcode)+"\n")
				pass
		print "Done 2VAR of protein coding variants..."
		return None

	def THREE_VAR(self,workdir):
		"""
		###---Three variants(3VAR) codon changes---
		##----------------------------------
		"""
		lines =  open(workdir+"temp_file2","r").read().splitlines()
		writ3 = self.write_file(workdir+"temp_file3")
		#---
		midline_head =  lines[0].strip()
		writ3.write(str(midline_head+"\t"+"Altcodon_merge-3VAR"+"\t"+"AA-change-3VAR")+"\n")
		#---
		i=0; TRcode =""; protcode = ""
		for i in range(len(lines)):
			True
			try:
				#---
				midline_crp =  lines[i].strip()
				#---
				beforeline = lines[i-1].strip()
				line0 = beforeline.split("\t")[-5]
				beforeline1 =  re.findall("\d+", line0)
				midline = lines[i].strip()
				line1 = midline.split("\t")[-5]
				midline1 =  re.findall("\d+", line1)
				nextline = lines[i+1].strip()
				line2 = nextline.split("\t")[-5]
				nextline1 =  re.findall("\d+", line2)

				line11=[];line22=[];writelst= [];
				if ((beforeline1[0] == midline1[0]) or (midline1[0] == nextline1[0])):
					True
					lscod1=[];lscod2=[]
					line11 = lines[i].strip().split("\t")[-2]
					if (midline1[0] == nextline1[0]):
						line22 = lines[i+1].strip().split("\t")[-2]
					for cod in line11:
						lscod1.append(cod)
					for cod in line22:
						lscod2.append(cod)
					#-----------------------
					if(len(lscod1)==3 and len(lscod2)==3):
						True
						threeltr_code =[];

						if ((lscod1[0].isupper()==True) and (lscod2[0].islower()==True) or (lscod1[0].isupper()==True and lscod2[0].isupper()==True) ):
							True
							threeltr_code.append(lscod1[0])
						elif ((lscod1[0].islower()==True) and (lscod2[0].isupper()==True)):
							True
							threeltr_code.append(lscod2[0])

						if ((lscod1[1].isupper()==True) and (lscod2[1].islower()==True) or (lscod1[1].isupper()==True and lscod2[1].isupper()==True)):
							True
							threeltr_code.append(lscod1[1])
						elif ((lscod1[1].islower()==True) and (lscod2[1].isupper()==True)):
							True
							threeltr_code.append(lscod2[1])
			
						if ((lscod1[2].isupper()==True) and (lscod2[2].islower()==True) or (lscod1[2].isupper()==True and lscod2[2].isupper()==True)):
							True
							threeltr_code.append(lscod1[2])
						elif ((lscod1[2].islower()==True) and (lscod2[2].isupper()==True)):
							True
							threeltr_code.append(lscod2[2])
						#-----------------------

						if(len(threeltr_code)==3):
							True
							TRcode = ''.join(threeltr_code)
						else:
							True
							TRcode = 'Multiallelic-t1'
						#-----------------------
					else:
						True
						TRcode = '.'
					
					if(TRcode != line22):
						True
						protcode = self.translate_dna_single(TRcode)
						if len(protcode) ==0:
							protcode = ".";TRcode = "."
					else:				
						True
						TRcode = ".";

					if protcode == "" or TRcode == "":
						protcode = ".";TRcode = "."
					#-----------------------
					writ3.write(str(midline_crp+"\t"+TRcode+"\t"+protcode)+"\n")
					
					##----------------------
			except IndexError:
				if not midline_crp.split("\t")[0] == "CHROM":
					if len(protcode) ==0:
						protcode = ".";TRcode = "."
						writ3.write(str(midline_crp+"\t"+TRcode+"\t"+protcode)+"\n")
					else:
						protcode = ".";TRcode = "."
						writ3.write(str(midline_crp+"\t"+TRcode+"\t"+protcode)+"\n")
				pass
		print "Done 3VAR of protein coding variants..."
		return None
		
	def PARS_OUT_VAR(self,workdir):
		"""
		###---Pars variants(2VAR _ 3VAR) based on change of protein codons ---
		##----------------------------------
		"""
		#This first records the ones that do match and saves the value in column 3 (less the final matching amino acid residue).
		subprocess.check_output("awk 'index($(NF-6), $(NF-2)) {print}' "+workdir+"temp_file3 | awk '{print $(NF-6)}' | sed s'/.$//' > "+workdir+"matches.list", shell=True)

		#This then searches for the ones that do not match, and then also eliminates any value recorded in matches.list
		subprocess.check_output("awk '!index($(NF-6), $(NF-2)) { print }' "+workdir+"temp_file3 | grep -w -v -f "+workdir+"matches.list | awk -F'\t' '{ print }' > "+workdir+"temp_file4", shell=True)
		
		print "Done Pars 2VAR and 3VAR of protein coding variants..."
		return None
		
	def ZYGO_PAIR(self,readfl3,workdir,FIELDS):
		"""
		###---Pair of zygosity check: remove pair of codons for which one is reference homozygous ---
		##----------------------------------
		"""
	
		Fld_Len = int(len(FIELDS.split(",")))
		lines =  open(workdir+"temp_file4","r").read().splitlines()
		writ4 = self.write_file(workdir+readfl3)
		#---
		midline_head =  lines[0].strip()
		writ4.write(str(midline_head)+"\n")
		#---
		nextline_pos_lst=[]

		for i in range(len(lines)):
			#---
			midline_crp =  lines[i].strip()
			#---
			beforeline = lines[i-1].strip()
			beforeline_pos =  re.findall("\d+", beforeline.split("\t")[-7])
			midline = lines[i].strip()
			midline_pos =  re.findall("\d+", midline.split("\t")[-7])
			try:
				nextline = lines[i+1].strip()
				nextline_pos =  re.findall("\d+", nextline.split("\t")[-7])
				##-------
				if (midline_pos[0] == nextline_pos[0]):
					Start_lns = int(6)+Fld_Len
					End_lns = int(8)
					Start_zyg_lns = int(len(midline.split("\t"))-int(Start_lns))
					midline_zyg = midline.split("\t")[-int(Start_zyg_lns):-int(End_lns)]
					nextline_zyg = nextline.split("\t")[-int(Start_zyg_lns):-int(End_lns)]
					checkList = list(['0/1:1/0','1/0:0/1','1/0:1/0','0/1:0/1','1/0:1/1','1/1:1/0','0/1:1/1','1/1:0/1','1/1:1/1'])
					Mergetwozyg = ','.join([str(a)+":"+ b for a,b in zip(midline_zyg, nextline_zyg)])
					Mergetwozygsp =  Mergetwozyg.split(",")
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
					writ4.write(str(nextline)+"\n")

		subprocess.check_output("rm "+workdir+"temp_file*", shell=True)
		subprocess.check_output("rm "+workdir+"matches.list", shell=True)
		subprocess.check_output("rm "+workdir+"snp*", shell=True)

		print "Done checking pairs which are matching at least once with the zygosity list..."
		return None
		###------
		print "Done all pairs of codon variants..."

### Assigning arguments #--------

if __name__ == "__main__":
	parser = ArgumentParser(
	description="-Script to identify SNPClusters (SNPs within the same genetic code).")

	parser.add_argument("-i", "--infile", dest="INPUTFile",default=False, required=True, help="Full name of the input VCF file." )
	parser.add_argument("-o", "--outfile", dest="OUTPUTFile",default="./MACARON_output.txt", required=False, help="Full name of the output txt file (Default Output file: MACARON_output.txt)." )
	#
	parser.add_argument("-f", "--fields", dest="Fields", default="QUAL", required=False, help=" Single field name or comma-seperated ',' multiple field names can be given. Field name should be given according to the INFO header of the input vcf file. Example: -f Func.refGene,ExonicFunc.refGene,Gene.refGene,1000g2015aug_all,ExAC_ALL,ExAC_EAS,clinvar_20161128,gnomAD_exome_ALL,gnomAD_genome_ALL,EFF,CSQ")
	parser.add_argument("-d", "--dir", dest="Directory", default="./", required=False, help="Full path of the working directory; Input VCF file should be located within the same folder as output file will also be generated here. Default: Current directory")

	if len(sys.argv) == 2:
		parser.print_help()
		print "\nUSAGE: python MACARON_v0.3.py -i <inputfile.vcf>"
		print "Example1: python MACARON_v0.3.py -i yourinputfile.vcf"
		print "Example2: python MACARON_v0.3.py -i yourinputfile.vcf -o MACARON_output.txt -d /WORKDIR/FULLPATH/ -f 1000g2015aug_all,ExAC_ALL \n"
		sys.exit()

	args = parser.parse_args()
	INF = args.INPUTFile
	OUTF = args.OUTPUTFile
	DIR = args.Directory
	FIELDS = args.Fields

clF1 = SearchDB().Search_CODON(INF,DIR,FIELDS)
clF2 = SearchDB().TWO_VAR(DIR)
clF3 = SearchDB().THREE_VAR(DIR)
clF4 = SearchDB().PARS_OUT_VAR(DIR)
clF5 = SearchDB().ZYGO_PAIR(OUTF,DIR,FIELDS)
