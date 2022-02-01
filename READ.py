#!/usr/bin/python

Usage =""" 

Current READ version only supports .tped files.

READ.py <InputFile> <normalization> <normalization_value>

A normalization is only required when a user-defined values is used instead of the median, mean or maximum of the test population. All normalization settings are described below:

    median (default) - assuming that most pairs of compared individuals are unrelated, READ uses the median across all pairs for normalization.
    mean - READ uses the mean across all pairs for normalization, this would be more sensitive to outliers in the data (e.g. recent migrants or identical twins)
    max - READ uses the maximum across all pairs for normalization. This should be used to test trios where both parents are supposed to be unrelated but the offspring is a first degree relative to both others.
    value <val> - READ uses a user-defined value for normalization. This can be used if the value from another population should be used for normalization. That would be useful if only two individuals are available for the test population. A value can be obtained from the NonNormalizedP0 column of the file meansP0_AncientDNA_normalized from a previous run of READ.

Optionally, one can add --window_size <value> at the end in order to modify the window size used (default: 1000000).

"""

print("   ===Thank you for using Relationship Estimation from Ancient DNA (READ)===")
print("")
print("")

import re
import random
import sys
import subprocess
import os

path_Rscript='/usr/bin/Rscript'
Arguments = sys.argv[1:]
norm_method="median"
norm_value=''
Filetype='TPED'
window_size=1000000

if len(Arguments)==0:
	sys.exit(Usage)

if len(Arguments)>=2 and Arguments[1] in ["median","mean","max","value"]:
	norm_method=Arguments[1]
if len(Arguments)>=3 and Arguments[1]=="value":
	norm_value=float(Arguments[2])
if "--window_size" in Arguments:
	ws_index=Arguments.index("--window_size")
	if len(Arguments)>(ws_index+1):
		window_size=int(Arguments[ws_index+1])
	else:
		sys.exit("No window size specified!")
	
	
	
InFileName= Arguments[0] + '.tped'
try:
	open(InFileName,'r')
except IOError:
	print("the file %s.tped does not exist" % InFileName)

InFile = open(InFileName, 'r')

OutFileName='Read_intermediate_output'
OutFile=open(OutFileName, 'w')



dictionary_pair_individuals={}
dictionary_pair_individuals_missinginfo={}
list_all_individuals=[]
list_2_all_individuals=[]
previous_window=0
Missing_info=0
IBS2=0
IBS0=0
Total_snps_window=0
full_call_var=0
full_call_var_corrected=0
List_individuals=[]

print("Starting pairwise comparisons in windows")
print("")


if Filetype == "TPED":
	FirstLine = InFile.readline()
	FirstLine = FirstLine.strip("\n")
	ElementsFirstLine = FirstLine.split("\t")
	Genotype = ElementsFirstLine[4:]
	Number_individuals= len(Genotype)/2
	for l in open(Arguments[0] + '.tfam'):
		split=l.split()
		list_all_individuals.append(split[1])
	#for i in range(1,Number_individuals+1):
	#	list_all_individuals.append('Ind%d' % i)
	list_2_all_individuals.extend(list_all_individuals)
	for idx1, j in enumerate(list_all_individuals):
		for idx2, i in enumerate(list_2_all_individuals):
			if (i == j or idx1 >= idx2):
				continue
			else:
				dictionary_pair_individuals["%s%s" % (j,i)]=0
				dictionary_pair_individuals_missinginfo["%s%s" % (j,i)]=0
				List_individuals.append("%s%s" % (j,i))
	InFile.seek(0)
	OutFile.write("PairIndividuals\tChromosome\tWindowIndex\tSNVperWindow\tIBS2\tIBS0\tP1\tP0\tMissing\n")
	previous_window=0
	snp_count=0
	
	for Line in InFile:
		snp_count+=1
		Line = Line.strip("\n")
		ElementList = Line.split()
		Chromosome = int(ElementList[0])
		Position = int(ElementList[3])
		Genotype= ElementList[4:]
		window_index=Position/window_size

		if (snp_count%5000)==0:
			print("Current position: Chromosome %s, bp %s" %(Chromosome,Position))
		
		Genotype="".join(Genotype)
		Alleles_individuals = [Genotype[i:i+2] for i in range(0, len(Genotype), 2)]
		if window_index == previous_window:
			full_call_var+=1
			for idx1, i in enumerate(list_all_individuals):
				for idx2, j in enumerate(list_2_all_individuals): 
						if (i == j or idx1 >= idx2): #I think this comparison is wrong
							continue
						elif (Alleles_individuals[idx2]=="00") or (Alleles_individuals[idx1]=="00"):
							dictionary_pair_individuals_missinginfo["%s%s" % (list_all_individuals[idx1],list_all_individuals[idx2])]+=1
						elif Alleles_individuals[idx1] == Alleles_individuals[idx2]:
							dictionary_pair_individuals["%s%s" % (list_all_individuals[idx1],list_all_individuals[idx2])]+=1

		else: 
			IBS2= [dictionary_pair_individuals[pair] for pair in List_individuals]
			Missing_info = [dictionary_pair_individuals_missinginfo[pair] for pair in List_individuals]
			for valor, missing, pair in zip(IBS2,Missing_info,List_individuals):
				full_call_var_corrected = full_call_var - missing
				IBS0= float(full_call_var_corrected - valor)
				if (full_call_var_corrected!=0):
					P2 = float(valor)/full_call_var_corrected
					P0 = float(IBS0/full_call_var_corrected)
					OutString= "%s\t%s\t%s\t%i\t%i\t%i\t%f\t%f\t%i\t%i" % (pair, Chromosome, window_index, full_call_var, valor, IBS0,P2,P0,missing,full_call_var_corrected)
					OutFile.write(OutString+ '\n')
				
			for key in dictionary_pair_individuals:
				dictionary_pair_individuals[key]=0	
			for key in dictionary_pair_individuals_missinginfo:
				dictionary_pair_individuals_missinginfo[key]=0
			full_call_var=1
			previous_window= window_index	
				
			for idx1, i in enumerate(list_all_individuals):
				for idx2, j in enumerate(list_2_all_individuals):
					if (i == j or idx1 >= idx2): #I think this comparison is wrong
						continue
					elif (Alleles_individuals[idx2]=="00") or (Alleles_individuals[idx1]=="00"):
						dictionary_pair_individuals_missinginfo["%s%s" % (list_all_individuals[idx1],list_all_individuals[idx2])]+=1
					elif Alleles_individuals[idx1] == Alleles_individuals[idx2]:
						dictionary_pair_individuals["%s%s" % (list_all_individuals[idx1],list_all_individuals[idx2])]+=1
else:
	print("Format was not specified or it is not supported.")

OutFile.close()
InFile.close()

print("")

###Order intermediate Output
OutFile_1=open('READ_output_ordered_firstline', 'w') #This will be the ordered Output
OutFile_1.write("PairIndividuals\tChromosome\tWindowIndex\tSNVperWindow_original\tIBS2\tIBS0\tP2\tP0\tMissing\tSNVperWindow\n")
OutFile_1.close()
list_filestocat=[]
list_filestocat.append('READ_output_ordered_firstline')


for a in List_individuals:
	list_filestocat.append("READ_output_ordered_%s.txt" %(a))
	
InFileName="Read_intermediate_output"
InFile = open(InFileName, 'r')
line_number=0
for Line in InFile:
	line_number+=1
	if Line.startswith("PairIndividuals"):
		continue
	Line = Line.strip("\n")
	ElementList = Line.split("\t")
	Pair_individuals= ElementList[0]
	Chromosome = int(ElementList[1])
	WindowIndex = int(float(ElementList[2]))
		
	for a in List_individuals:
		if a == Pair_individuals:
			with open("READ_output_ordered_%s.txt" %(a), "a") as f:
				f.write(Line + "\n")


my_cmd = ['cat'] + list_filestocat
with open('READ_output_ordered', "w") as outfile:
	subprocess.call(my_cmd, stdout=outfile)

InFile.close()

print("Normalization started")

print("")

###Run R script
import subprocess
process = subprocess.Popen("%s READscript.R %s %s" %("Rscript",norm_method,norm_value), stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 
###

out, err = process.communicate()

print(out)
print(err)

print("Estimating degree of relationships")

###Predict relationships
InFileName= "meansP0_AncientDNA_normalized"
InFile = open(InFileName, 'r')
OutFileName="READ_results"
OutFile=open(OutFileName, 'w')

OutFile.write("PairIndividuals\tRelationship\tZ_upper\tZ_lower\n")

for Line in InFile:
	if Line.startswith('Pair'):
		continue
	else:
		Line = Line.strip("\n")
		ElementList = Line.split(" ")
		PairIndividuals=ElementList[0]
		P0=float(ElementList[1])
		err=float(ElementList[2])
		if (P0 >= 0.90625):
			Relationship="Unrelated"
			Zup='NA'
			Zdown=(0.90625-P0)/err
		elif (P0 >= 0.8125):
			Relationship="Second Degree"
			Zup=(0.90625-P0)/err
			Zdown=(0.8125-P0)/err
		elif (P0 >= 0.625):
			Relationship="First Degree"
			Zup=(0.8125-P0)/err
			Zdown=(0.625-P0)/err
		else:
			Relationship="IdenticalTwins/SameIndividual"
			Zup=(0.625-P0)/err
			Zdown='NA'
		
		OutString="%s\t%s\t%s\t%s" % (PairIndividuals,Relationship,Zup,Zdown)
		OutFile.write(OutString+ '\n')

InFile.close()
OutFile.close()

print("")
###

for file in list_filestocat:
	os.remove(file)
#os.remove("Read_intermediate_output")
#os.remove("READ_output_ordered")

print("READ analysis finished. Please check READ_results for results!")
