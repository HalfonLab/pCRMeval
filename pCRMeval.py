#!/usr/bin/env python
from __future__ import absolute_import, division, print_function



####################################################################################------BRIED DESCRIPTION-------###################################################################################################################################
# Halfon Lab																																																													#
#Date: Sept 2018																																																									#
#Purpose: This pipeline is written to evaluate the results of Supervised Enhancer Prediction Algorithm such as SCRMshaw, imogene etc..Specifically it will acess the Training Set Sensitivity, REDfly Recovery and Expression Pattern Precision	#
#Input and Output: To learn the details about the required input files and what this pipeline would generate as its output please see README enclosed in the package																					#
#####################################################################################################################################################################################################################################################
import os
import pybedtools
import statistics
import argparse
import sys
import shutil
from scipy import stats
import scipy.stats
from collections import Counter
import numpy as np
import pandas as pd
import csv
import itertools
import pprint


##################################################-------FUNCTIONS--------#########################################################
# This function will parse file of multiple outputs of scrmshaw to individual unique files for each training set and statisitcal method used. This will return three methods dictionaries with training sets as their keys
def parse_output(outfile,numparse): 
	# creating three separate dictionaries based on methods used  
	d_hexmcd={}
	d_imm={}
	d_pac={}
	x={''}
	global d_hexmcd_val
	global d_imm_val
	global d_pac_val
	d_hexmcd_val=['']
	d_imm_val=['']
	d_pac_val=['']
	
	with open (outfile) as file:
		rows=(line2.split('\t') for line2 in file)

		for row in rows:
		#based on the 14th column(names of different data sets) and 15th column (statistical method used) of scrmshaw_joined_output file giving values to each of the three method`s dictionaries
			if (row[16]=='hexmcd') and (int(row[17]) <= int(numparse)):
				#print(numparse)
				if row[15] not in d_hexmcd:
					myRow = [] # create a new list to use
					myRow.append(row) # add my new row to our new list
					d_hexmcd[row[15]] = myRow  #create a new entry if it isn't in the dictionary already

				else:
					d_hexmcd.get(row[15]).append(row)
					#count_hexmcd=count_hexmcd+1
			elif (row[16]=='imm')and (int(row[17]) <= int(numparse)):
				if row[15] not in d_imm:
					myRow = []
					myRow.append(row)
					d_imm[row[15]] = myRow
				else:
					d_imm.get(row[15]).append(row)
			elif (row[16]=='pac') and (int(row[17]) <= int(numparse) ):
				if row[15] not in d_pac:
					myRow = []
					myRow.append(row)
					d_pac[row[15]] = myRow
				else:
					d_pac.get(row[15]).append(row)
					#count_pac=count_pac+1

	#calculating number of keys(datasets) each dictionary ends up having		
		for key in d_hexmcd.keys():
			d_hexmcd_val.append(key)
		for key in d_imm.keys():
			d_imm_val.append(key)
		for key in d_pac.keys():
			d_pac_val.append(key)
			
	#creating separate files for each method w.r.t datasets, using the above newly created three dictionaries and moving them to tmp (temporary folder).
	
	#These individual unique files(based on methods and their data sets) will be used as input to perform the functions for getting the values to fill up the outfile file
	for key in d_hexmcd.keys():
		noOflines=len(d_hexmcd[key])
		
		with open(os.path.join(subdirectory,'hexmcd_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_hexmcd[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")
	
	for key in d_imm.keys():
		with open(os.path.join(subdirectory,'imm_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_imm[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")

	for key in d_pac.keys():
		with open(os.path.join(subdirectory,'pac_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_pac[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")
				
	return(d_imm_val,d_hexmcd_val,d_pac_val)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#This function will take file of scrmshaw output of unique dataset and method. It will sort and merge them and return number of sorted and merged scrmshaw predictions and also the path of file which will have these, and save it as "scrms.merged.X.bed"

def sort_and_merge_output(dName,mName,path):

	#sorting and merging the intervals of modified crms
	a=pybedtools.BedTool(path)
	b=a.sort().merge().saveas(subdirectory+'/scrms.merged'+dName+'_'+mName+'.bed')
	numOfsortedMergedPredictedCrms=b.count()
	
	for root, dirs, files in os.walk(os.getcwd()):
					for name in files:
						if name=='scrms.merged'+dName+'_'+mName+'.bed':
							p=os.path.abspath(os.path.join(root,name))
	return p,numOfsortedMergedPredictedCrms



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#This function will extract user specified number of predictions(or top good hits calculated from top good predictions) from the given training set's predictions

def parse_output_individual(indoutfile,parse,dataset,method): 
	# creating three separate dictionaries based on methods used  

	count=0
	#print("inside function")
		
	#print(parse)
	
	if nameOfmethod == 'SCRMshaw' or nameOfmethod =='SCRMSHAW' or nameOfmethod =='scrmshaw':
	#
		with open(subdirectory+"/"+method+"_"+dataset+"_fullLength.bed") as file2, open(subdirectory+"/"+method+"_"+dataset+"tmp.bed",'w') as tmp:		
			for line2 in file2:
			
				#print("opened file")
				row=line2.split('\t')
				lastCol=len(row)-1
				if int(row[lastCol]) <= parse:
					#print(row[16])
					tmp.write(line2)
					count+=1
				
	else:
		with open(indoutfile) as file2, open(subdirectory+"/"+method+"_"+dataset+"tmp.bed",'w') as tmp:		
			for line2 in file2:
		
				#print("opened file")
				row=line2.split('\t')
				lastCol=len(row)-1
				if int(row[lastCol]) <= parse:
					#print(row[16])
					tmp.write(line2)
					count+=1
				
		
	outfile=open(os.path.join(subdirectory,method+'_'+dataset+'.bed'),'w')
	#print(outfile)
	with open(subdirectory+"/"+method+"_"+dataset+"tmp.bed") as file6:	
		for line6 in file6:
			outfile.write(line6)
	
# 			
	pathThis=subdirectory+"/"+method+"_"+dataset+".bed"
	#return number of lines extracted
	#print(count)
	return count, pathThis
	
#---------------------------------------------------------------------------------------------------------------------------
#This function taked the tab delimited file and returns py bed version of it to perform the functions of bedtools

def bed_conversion(tab_delimited_path):
	
	bed_tabdelimited=pybedtools.BedTool(tab_delimited_path)
		
	return bed_tabdelimited
	

#---------------------------------------------------------------------------------------------------------------------------

#This function will take the excluded crms file and sort it to another file saved as "sortedExcluded_X_Crms.bed"
def sort_and_merge_tset(tset_Bed,dN,mN):
	a=tset_Bed.sort().merge().saveas(subdirectory+'/sorted_and_merged'+dN+'_'+mN+'_tsetCrms.bed')
	#size_of_sortedMergedTsetCrms=a.count()
	
	for root, dirs, files in os.walk(os.getcwd()):
					for name in files:
						if name=='sorted_and_merged'+dN+'_'+mN+'_tsetCrms.bed':
							p3=os.path.abspath(os.path.join(root,name))
	#return p3,size_of_sortedMergedTsetCrms
	return p3

#----------------------------------------------------------------------------------------------------------------------------

#This function will take the modified crms file of redfly(output of above function) and sort and merge it to another file saved as "sortedMergedModified_X_Crms.bed"
def sort_and_merge_modifiedCrms(ftype,file,dataset,method):
	a=pybedtools.BedTool(file)
	b=a.sort().merge(c=4,o='distinct').saveas(subdirectory+'/sortedMergedModified_'+ftype+'_'+dataset+'_'+method+'_Crms.bed')
	size_of_sortedMergedModifiedCrms=b.count() #Step 5

	for root, dirs, files in os.walk(os.getcwd()):
					for name in files:
						if name=='sortedMergedModified_'+ftype+'_'+dataset+'_'+method+'_Crms.bed':
							p3=os.path.abspath(os.path.join(root,name))
	return p3,size_of_sortedMergedModifiedCrms
#----------------------------------------------------------------------------------------------------------------------------

#This function will sort the the expression mapped crms file excluding the ones used as training sets(called as Modified expressionMapped crms file)
def sort_expressionMappedCrms(ftype,file,dataset,method):
	a=pybedtools.BedTool(file)
	b=a.sort().saveas(subdirectory+'/sortedModified_'+ftype+'_'+dataset+'_'+method+'_Crms.bed')
	size_of_sortedMergedModifiedCrms=b.count() #Step 5

	for root, dirs, files in os.walk(os.getcwd()):
					for name in files:
						if name=='sortedModified_'+ftype+'_'+dataset+'_'+method+'_Crms.bed':
							p3=os.path.abspath(os.path.join(root,name))
	return p3,size_of_sortedMergedModifiedCrms


#----------------------------------------------------------------------------------------------------------------------------
#This function will take two sorted and and merged files(1: Scrmshaw unique output 2: Excluded Crms) and return the number of common crms and save it to files named as "Intersected_excludedcrms & __name of dataset of scrmshaw output_.bed"

def intersect_Scrm_and_Excluded_crms(methName,dtName,sortedMergedScrmsBedV,sortedMergedExcludedCrmsBed):
	count=0
	scrm_training=sortedMergedExcludedCrmsBed.intersect(sortedMergedScrmsBedV,wa=True,u=True,F=0.10).saveas(subdirectory+'/Intersected'+methName+'_'+dtName+'_excludedcrms_and_predictions.bed')
	
	#creating a file which will have scrms ranks with it 
	
	count=scrm_training.count()

	return count


#------------------------------------------------------------------------------------------------

#This function will take FullCrms (redfly file) and exclude those crms which were used as training data for the particular scrmshaw output. It will save this file as "modifiedCrms_X_.bed"
def exclude_training_set_crms(ftype,dName,mName,inp,knCrms):
	#creating file of Modified crms( that does not have training set crms)
	countTotalCrms=0
	countExclCrms=0
	countModifiedCrms=0

	file= open(os.path.join(subdirectory,"modified_crms_"+ftype+'_'+dName+"_"+mName+".txt"),"w")
	file2= open(os.path.join(subdirectory,"excludedcrms"+ftype+'_'+dName+"_"+mName+".txt"),"w")
	#print(len(d[inp]))
	with open(knCrms) as f:
		for line in f:
			#linecount=linecount+1
			countTotalCrms=countTotalCrms+1
			#print(str(line))
			if any(x.rstrip() in line.split() for x in d[inp]):
				#print(line)
				file2.write(str(line))
				countExclCrms=countExclCrms+1
				#next(f)
			else:
				file.write(line)
				countModifiedCrms=countModifiedCrms+1


	file.close()
	file2.close()
	f.close()
	#converting the modified crms to tab delimited (bed file)	
	with open(os.path.join(subdirectory,'modifiedCrms'+ftype+'_'+dName+'_'+mName+'.bed'),'w') as outfile:
		with open(subdirectory+"/modified_crms_"+ftype+'_'+dName+"_"+mName+".txt") as infile:
			for line in infile:
				line=line.strip('chr')
				outfile.write(" ".join(line.split()).replace(' ', '\t'))
				outfile.write("\n")
	#converting excluded to bed file
	with open(os.path.join(subdirectory,'excluded_crms'+ftype+'_'+dName+'_'+mName+'.bed'),'w') as outfile:
		with open(subdirectory+"/excludedcrms"+ftype+'_'+dName+"_"+mName+".txt") as infile:
			for line in infile:
				line=line.strip('chr')
				outfile.write(" ".join(line.split()).replace(' ', '\t'))
				outfile.write("\n")

	#removing earlier modified crms files without tab
	my_d=os.getcwd()+'/'+subdirectory+'/'
	os.remove(my_d+'modified_crms_'+ftype+'_'+dName+'_'+mName+'.txt')
	os.remove(my_d+'excludedcrms'+ftype+'_'+dName+'_'+mName+'.txt')
	
	for root, dirs, files in os.walk(os.getcwd()):
					for name in files:
						if name=='modifiedCrms'+ftype+'_'+dName+'_'+mName+'.bed':
							p2=os.path.abspath(os.path.join(root,name))
	
	for root, dirs, files in os.walk(os.getcwd()):
					for name in files:
						if name=="excluded_crms"+ftype+'_'+dName+"_"+mName+".bed":
							p4=os.path.abspath(os.path.join(root,name))
	
	
	return countTotalCrms,countModifiedCrms,countExclCrms,p2,p4


#----------------------------------------------------------------------------------------------------------------------------
#This function will take FullCrms (redfly file) and exclude those crms which were used as training data(provided in Bed format not list format). It will save this file as "X_modifiedCrms_.bed"
def exclude_tsetBedFormat(ftype,methName,dtName,tsetbed,redflybed):

	#bedtools intersect -a without_chrallredfly_2.5kb.July2017.txt -b tsetCrmsBed.bed -v -f 0.90 -r > modTest.txt
	resultOfintersect=redflybed.intersect(tsetbed,v=True,f=0.90,r=True).saveas(subdirectory+'/'+ftype+'_'+methName+'_'+dtName+'_modifiedCrms.bed')
	
	return(resultOfintersect)
	
#----------------------------------------------------------------------------------------------------------------------------

#This function will take two sorted and and merged files(1: Scrmshaw unique output 2: Modified Crms) and return the number of common crms and saving it to files named as "Intersected_modifiedcrms & __name of dataset of scrmshaw output_.bed"
def intersect_scrms_and_modifiedCrms(ftype,methName,dtName,sortedMergedModifiedCrmsBed,sortedMergedScrmsBed):
	
	a_and_b_count=0
	countOfHits=0
	
	#getting the intersection keeping all data intact and saving the file to look later what crms got intersected  
	sortedMergedModifiedCrmsBed.intersect(sortedMergedScrmsBed,wao=True,F=0.10).saveas(subdirectory+'/Intersected_Scrms_and_modified_'+ftype+'_'+methName+'_'+dtName+'.bed')
	#for counting no of hits discarding lines other than hits : because we need the number of hits  SEE IF ORDER MAKES DIFFerence to NO OF HITs
	a_and_b=sortedMergedModifiedCrmsBed.intersect(sortedMergedScrmsBed,wa=True,F=0.10).saveas(subdirectory+'/IntersectedNotIntact_Scrms_and_modified_'+ftype+'_'+methName+'_'+dtName+'.bed')
	#print("A: Redfly B:SCRMs keeping redfly:Total no of lines "+str(a_and_b.count()))
	##Flipped
	a_and_b2=sortedMergedScrmsBed.intersect(sortedMergedModifiedCrmsBed,wa=True,u=True,f=0.10).saveas(subdirectory+'/FlippedIntersectedNotIntact_Scrms_and_modified_'+ftype+'_'+methName+'_'+dtName+'.bed')	
	#print("A:SCRMs B: Redfly keeping A scrms (Should be the actual one) why is this wrong? may be also not the lines but the names of crms inn each:Total no of lines "+str(a_and_b2.count()))
	
	#based on predictions how many of the predictions overlapped with how many redfly crms
	if ftype=="expressionMapped":
		
		a_and_b_new=sortedMergedScrmsBed.intersect(sortedMergedModifiedCrmsBed,wo=True,f=0.10)
		#print(a_and_b_new)
		if a_and_b_new != '':
			a_and_b_new.merge(c=7,o='distinct').saveas(subdirectory+'/ResultOfIntersectBasedOnScrms_modifiedcrms_and_'+ftype+'_'+methName+'_'+dtName+'.bed')
		else:
			print('No intesect found')
			a_and_b_new.saveas(subdirectory+'/ResultOfIntersectBasedOnScrms_modifiedcrms_and_'+ftype+'_'+methName+'_'+dtName+'.bed')
		#print(a_and_b_new)
	#counting number of hits (crms of scrmshaw output that matched with modified crms)
	a_and_b_count=a_and_b2.count()
	#print(a_and_b_count)
	if ftype=="expressionMapped":
		for root, dirs, files in os.walk(os.getcwd()):
				for name in files:
					if name=='ResultOfIntersectBasedOnScrms_modifiedcrms_and_'+ftype+'_'+methName+'_'+dtName+'.bed':
						p6=os.path.abspath(os.path.join(root,name))
	
	else:
		for root, dirs, files in os.walk(os.getcwd()):
						for name in files:
							if name=='IntersectedNotIntact_Scrms_and_modified_'+ftype+'_'+methName+'_'+dtName+'.bed':
								p6=os.path.abspath(os.path.join(root,name))
							

	return a_and_b_count,p6

#--------------------------------------------------------------------------------------------------------------------------
#This function is counting number of predictions intersected with redfly crms excluding the ones used for training set, from the intersected file created from the function of intersect
def redfly_recovered(IntersectScrmsRedflyResultFile):
	countNoOfHits=0
	redflyRecovered=['']
	with open (IntersectScrmsRedflyResultFile) as file:
		rows=(line.split('\t') for line in file)
		for row in rows:
			if scrmshawOrNot=="True" or scrmshawOrNot=="T" or scrmshawOrNot=="TRUE":
				if row[9]!='0\n':
					countNoOfHits=countNoOfHits+1
			else:
				if row[7]!='0\n':
					countNoOfHits=countNoOfHits+1

	
	return countNoOfHits


#----------------------------------------------------------------------------------------------------------------------------
#This will shuffle through the non coding part of genome (and not drawing out from our predictions file) and calculate all the stats like mean median p, z value etc 
def shuffle(tmp_redfly_merged,scrms_merged,e,ge,a_and_b_count,numOfShuffles):
	no_of_intersects=[]
	i=""
	mean=0
	sd=0
	p=0
	z=0
	median=0
	maxV=0
	minV=0
	for i in range(numOfShuffles):
		shuffled=scrms_merged.shuffle(excl=e, noOverlapping=True, g=ge)
		#and then bedtools intersect -a scrms.merged.bed -b shuffled2.bed -f 0.10 >step8aOutput2.bed 
		intersect_shuffled_and_scrms_merged=tmp_redfly_merged.intersect(shuffled,f=0.10)
		intersect_shuffled_and_scrms_merged.saveas('shuffled.bed')
		# actually need only number of lines 
		no_of_intersects.append(intersect_shuffled_and_scrms_merged.count())
		#print(no_of_intersects[i])

	#get summary stats 
	maxV=max(no_of_intersects)
	minV=min(no_of_intersects)
	mean=statistics.mean(no_of_intersects)
	median=statistics.median(no_of_intersects)
	sd=statistics.stdev(no_of_intersects)

	summary_stats="Mean= "+ str(mean)+"\nMedian= "+str(median)+"\nStandard Deviation= "+str(sd)+"\nMinimum value= "+str(minV)+"\nMaximum Value= "+str(maxV)

	#from a_and_b_count and mean sd of above array
	#print(a_and_b_count)
	try:	
		z=(a_and_b_count-mean)/sd
		#print(z)
		p=stats.norm.sf(abs(z))
		#print(p)
	except ZeroDivisionError:
		np.seterr(divide='ignore', invalid='ignore')
	return mean,median,sd,minV,maxV,summary_stats,z,p
#----------------------------------------------------------------------------------------------------------------------------
#This function will calculate size of sorted and merged full redfly2010 and subset redfly file(expression mapped crms file)
def sort_merge_redfly_find_size(ftype,file,dataset,method):
	
	newFileName='temporaryFileName.bed'
	with open(os.path.join(subdirectory,newFileName),'w') as outfile:
		with open(file) as infile:
			for line in infile:
				line=line.strip('chr')
				outfile.write(" ".join(line.split()).replace(' ', '\t'))
				outfile.write('\n')
			
	for root, dirs, files in os.walk(os.getcwd()):
			for name in files:
				if name==newFileName:
					newFilePath=os.path.abspath(os.path.join(root,name))		
	
	#newFilePath=os.path.abspath(newFileName)
	a=pybedtools.BedTool(newFilePath)
	#a=pybedtools.BedTool(file)
	b=a.sort().merge(c=4,o='distinct').saveas(subdirectory+'/sortedMerged_'+ftype+'_'+dataset+'_'+method+'_Crms.bed')
	size_of_sortedMergedCrms=b.count() #Step 5
	return size_of_sortedMergedCrms

#----------------------------------------------------------------------------------------------------------------------------
#This function will take in the intersected file of SCRMShaw predictions and expression mapped crms, and calculate if the predictions belong to the training's group of expression, if not then to which other group it belongs
def parse_intersected_scrmshaw(outfile,dSet,intersected,meth):
	countHits=0
	numberOfSCrmsMatchedToExpression=0
	numberOfSCrmsMatchedToExpression=0
	numberOfCrmsNotMatchedToOtherExpression=0
	hitsWithExpression=['']
	
	
	file2= open(os.path.join(subdirectory,"Hits_modifiedCrms_Scrms_"+dSet+"_"+meth),"w")
	file3=open(os.path.join(subdirectory,"Hits_NotMatchedWithItsGroup_modifiedCrms_Scrms_"+dSet+"_"+meth),"a")
	fileE=open(outfile+'_ExpressionPatternRecovereyDistribution.bed','a')
	
	with open (intersected) as file:
		rows=(line.split('\t') for line in file)
		for row in rows:
			alreadyDoneSet=['']
			matchedCrmsString=row[3].strip('\n')	
			countHits=countHits+1
			hitsWithExpression=matchedCrmsString.split(',')
			for i in range(len(hitsWithExpression)):
			
				for key in d2:
					if hitsWithExpression[i] in d2[key]:
						if key not in alreadyDoneSet:
							countd[key]=countd[key]+1
							alreadyDoneSet.append(key)
						else:
							continue
							#countd[key]=countd[key]+1
							#alreadyDoneSet.append(key)
							
			for i in range(len(hitsWithExpression)):
				#if pattern recovery is true then find out how many of the predicted crms belong to the correct group of expression
				if patternRecovery=='True' or patternRecovery == "T" or patternRecovery =="TRUE":
					if hitsWithExpression[i] in d2[dSet]:
						numberOfSCrmsMatchedToExpression=numberOfSCrmsMatchedToExpression+1
						
						#print(hitsWithExpression[i]+ " is in "+ matchedCrmsString)
						file2.write(hitsWithExpression[i])
						file2.write('\n')
						break
					else:
						numberOfCrmsNotMatchedToOtherExpression=numberOfCrmsNotMatchedToOtherExpression+1
						file3.write(matchedCrmsString)
						file3.write('\n')

		numberOfSCrmsNotMatchedToOtherExpression=countHits-numberOfSCrmsMatchedToExpression
		#print("Total Number Recovered that matched any expression group "+dSet+"_"+meth+": "+str(countHits)+"\nOut of those number of recovered that matched with its own expression group:"+str(numberOfSCrmsMatchedToExpression)+"\nRemaining number of recovered that matched with other expression groups:"+str(numberOfSCrmsNotMatchedToOtherExpression))
		
		#writing the distribution of crms across all the expression groups to a file		
		if continuous=='TRUE' or continuous=='True' or continuous=='T':
			fileE.write(dSet+"_"+meth+'_'+str(t)+":"+'\t')
		else:
			fileE.write(dSet+"_"+meth+":"+'\t')

		for keys in countd:
			fileE.write(str(countd[keys])+'\t')
			#print("number of occurences in group "+keys+":"+str(countd[keys]))
		
		fileE.write('\n')
		fileE.write("percentageExpressionPatternPrecision"+str(t)+'\t')
		
		for key in countd:

			#fileE.write()
			#print("denom"+str(len(d2[key]))+"nominator:"+str(countd[key]))
			fileE.write(str(countd[key]/countHits)+'\t')
		 	countd[key]=0
		 
		fileE.write('\n')

	file2.close()
	file3.close()

	if patternRecovery=='False' or patternRecovery=="F" or patternRecovery=="FALSE":
		return(numberOfSCrmsMatchedToExpression,len(d2['mapping1.pns']),countHits)
		
	else:
		return(numberOfSCrmsMatchedToExpression,len(d2[dSet]),countHits)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#This function will discard the peaks below the score cutoff determined by score curve
def extract_below_cutoff_vals(indoutfile,parse,dataset,method,numScore1,numAmp1):
	count=0
	with open(subdirectory+"/"+method+"_"+dataset+"tmp.bed") as fileA, open(subdirectory+"/"+method+"_"+dataset+"tmp2.bed",'w') as tmp2:
		for lineA in fileA:
			rowA=lineA.split('\t')

			#if the peak's score value is equal or above the score cutoff only then it will write to file
			if float(rowA[4]) >= numScore1:
				#print(row[16])
				tmp2.write(lineA)
				count+=1
	
	outfile=open(os.path.join(subdirectory,method+'_'+dataset+'.bed'),'w')
	with open(subdirectory+"/"+method+"_"+dataset+"tmp2.bed") as fileB:	
		for lineB in fileB:
			outfile.write(lineB)
		
	pathThis2=subdirectory+"/"+method+"_"+dataset+".bed"
	#print("now")
	#print(count)
	return count, pathThis2
#----------------------------------------------------------------------------------------------------------------------------

#This function will take in the intersected file of predictions and expression mapped crms, and calculate to which of the expression mapped group these predictions belong
def parse_intersected(outfile,intersected,dSet,meth):
	countHits=0
	numberOfSCrmsMatchedToExpression=0
	numberOfCrmsNotMatchedToOtherExpression=0
	hitsWithExpression=['']
	fileE=open(outfile+'_ExpressionPatternRecovereyDistribution.bed','a')
	#print(intersected)
	with open (intersected) as file:
		rows=(line.split('\t') for line in file)
		for row in rows:
			#print(row)
			alreadyDoneSet=['']
			matchedCrmsString=row[3].strip('\n')
			countHits=countHits+1
			hitsWithExpression=matchedCrmsString.split(',')
			
			for i in range(len(hitsWithExpression)):
				for key in d2:
					if hitsWithExpression[i] in d2[key]:
						if key not in alreadyDoneSet:
							countd_random[key]=countd_random[key]+1
							alreadyDoneSet.append(key)
						else:
							continue
							#countd[key]=countd[key]+1
							#alreadyDoneSet.append(key)
			
			
			for i in range(len(hitsWithExpression)):
			
				if hitsWithExpression[i] in d2[dSet]:
					numberOfSCrmsMatchedToExpression=numberOfSCrmsMatchedToExpression+1

					break
				else:
					numberOfCrmsNotMatchedToOtherExpression=numberOfCrmsNotMatchedToOtherExpression+1

		numberOfSCrmsNotMatchedToOtherExpression=countHits-numberOfSCrmsMatchedToExpression
		#print("Total Number Recovered that matched any expression group "+dSet+"_"+meth+": "+str(countHits)+"\nOut of those number of recovered that matched with its own expression group:"+str(numberOfSCrmsMatchedToExpression)+"\nRemaining number of recovered that matched with other expression groups:"+str(numberOfSCrmsNotMatchedToOtherExpression))
		
		#writing the distribution of crms across all the expression groups to a file	
		if continuous=='TRUE' or continuous=='True' or continuous=='T':
			fileE.write(dSet+"_"+meth+'_'+str(t)+":"+'\t')
		else:
			fileE.write(dSet+"_"+meth+":"+'\t')

		for keys in countd_random:
			fileE.write(str(countd_random[keys])+'\t')
			#print('writing')
		fileE.write('\n')
		
		#if the intersected file is not random(to calculate expected Precision) then calcule and write the Precision to file
		if not outfile.endswith("random"):
			fileE.write("percentageExpressionPatternPrecision"+str(t)+'\t')

			for key in countd_random:
					fileE.write(str(countd_random[key]/countHits)+'\t')
					countd_random[key]=0
		 	fileE.write('\n')
		 		
		else:

		 	for keys in countd_random:
		 		countd_random[keys]=0
		 		
	return(numberOfSCrmsMatchedToExpression,countHits)
	

#----------------------------------------------------------------------------------------------------------------------------
#This function will find out the score cutoff and amplitude cutoff based on their cutoffs 
def topGoodPrediction(unsorted_indScrOutputfile):
	row=[]
	valuesScore=[]
	valuesAmp=[]
	
	with open(unsorted_indScrOutputfile) as sfile:
		for line in sfile:
			row=line.split('\t')
			valuesScore.append(float(row[4]))
			valuesAmp.append(float(row[3]))
	#print(values)
	# pull out the list from pandas frame
	valuesScore=list(valuesScore)
	valuesScore=sorted(valuesScore,key=float,reverse=True)
	#valuesAmp=sorted(valuesAmp,key=float,reverse=True)
	
	#for scores cutoff
	
	#get coordinates of all the points
	nPointsScore = len(valuesScore)
	allCoordScore = np.vstack((range(nPointsScore), valuesScore)).T
	#np.array([range(nPoints), values])

	# get the first point
	firstPointScore = allCoordScore[0]
	# get vector between first and last point - this is the line
	lineVecScore = allCoordScore[-1] - allCoordScore[0]
	lineVecNormScore = lineVecScore / np.sqrt(np.sum(lineVecScore**2))

	# find the distance from each point to the line:
	# vector between all points and first point
	vecFromFirstScore = allCoordScore - firstPointScore

	# To calculate the distance to the line, we split vecFromFirst into two 
	# components, one that is parallel to the line and one that is perpendicular 
	# Then, we take the norm of the part that is perpendicular to the line and 
	# get the distance.
	# We find the vector parallel to the line by projecting vecFromFirst onto 
	# the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
	# We project vecFromFirst by taking the scalar product of the vector with 
	# the unit vector that points in the direction of the line (this gives us 
	# the length of the projection of vecFromFirst onto the line). If we 
	# multiply the scalar product by the unit vector, we have vecFromFirstParallel
	scalarProductScore = np.sum(vecFromFirstScore * np.matlib.repmat(lineVecNormScore, nPointsScore, 1), axis=1)
	vecFromFirstParallelScore = np.outer(scalarProductScore, lineVecNormScore)
	vecToLineScore = vecFromFirstScore - vecFromFirstParallelScore
	# distance to line is the norm of vecToLine
	distToLineScore = np.sqrt(np.sum(vecToLineScore ** 2, axis=1))
	# knee/elbow is the point with max distance value
	idxOfBestPointScore = np.argmax(distToLineScore)
	
	#print("value of score at score elbow")
	#print(valuesScore[idxOfBestPointScore])
	scoreAtElbow=valuesScore[idxOfBestPointScore]
	
	#for amplitude cutoff
	
	valuesAmp=list(valuesAmp)
	#print(valuesAmp)
	#get coordinates of all the points
	nPointsAmp = len(valuesAmp)
	allCoordAmp = np.vstack((range(nPointsAmp), valuesAmp)).T
	#np.array([range(nPoints), values])

	# get the first point
	firstPointAmp = allCoordAmp[0]
	# get vector between first and last point - this is the line
	lineVecAmp = allCoordAmp[-1] - allCoordAmp[0]
	lineVecNormAmp = lineVecAmp / np.sqrt(np.sum(lineVecAmp**2))

	# find the distance from each point to the line:
	# vector between all points and first point
	vecFromFirstAmp = allCoordAmp - firstPointAmp
	scalarProductAmp = np.sum(vecFromFirstAmp * np.matlib.repmat(lineVecNormAmp, nPointsAmp, 1), axis=1)
	vecFromFirstParallelAmp = np.outer(scalarProductAmp, lineVecNormAmp)
	vecToLineAmp = vecFromFirstAmp - vecFromFirstParallelAmp
	# distance to line is the norm of vecToLine
	distToLineAmp = np.sqrt(np.sum(vecToLineAmp ** 2, axis=1))

	# knee/elbow is the point with max distance value
	idxOfBestPointAmp = np.argmax(distToLineAmp)
	#print("value at 843")
	#print(valuesAmp[idxOfBestPointAmp])
	ampAtElbow=valuesAmp[idxOfBestPointAmp]
	#print(idxOfBestPointScore,idxOfBestPointAmp)
	
	return idxOfBestPointScore,scoreAtElbow,ampAtElbow
	
	
#----------------------------------------------------------------------------------------------------------------------------
#This function will calculate summary stats for the expected Expression mapped crms recovery using pandas
def calculating_sum_stats_random_file(outfile):
	setNum='mean'+str(t)
	my_file_name=outfile+'_random_ExpressionPatternRecovereyDistribution.bed'
	df=pd.read_csv(my_file_name,sep='\t',skiprows=(0,1),header=None)
	dfSummary=df.describe()
	dfSummary= dfSummary.rename({'mean':setNum},axis='index')
	#dfSummary=dfSummary.dropna()
	dfSummary.to_csv('summaryRadom.txt',sep='\t')
	
	with open('summaryRadom.txt','r') as sumFile, open(outfile+'_ExpressionPatternRecovereyDistribution.bed','a') as outf:
		for line in sumFile:
			#if line.startswith('mean') or line.startswith('std') or line.startswith('min') or line.startswith('max'):
			if line.startswith('mean') or line.startswith('std'):
				outf.write(line)
				

#----------------------------------------------------------------------------------------------------------------------------	
#This function is used to transpose the expression pattern recovery distribution file
def transposeDistributionFile(outfile):
	fileN=outfile+'_ExpressionPatternRecovereyDistribution.bed'
	
	# Read all Tab-delimited rows from stdin.
	all_data = []
	with open(fileN,'r') as file:
		for line in file:
		    all_data.append(line.strip().split('\t'))
	file.close()
	
	# Transpose it.
	transposed = []
	for row_i in range(len(all_data[0])):
	    new_row = []
	    for col_i in range(len(all_data)):
	        new_row.append(all_data[col_i][row_i])
	    transposed.append(new_row)
	
	# Write it back out
	with open(fileN, 'w') as fileT:
		for row in transposed:
			line='\t'.join(row)
			fileT.write(line+'\n')		


#----------------------------------------------------------------------------------------------------------------------------		'
#This function is used to display the summary results of Precision of each of the training set as it goes through the loop on the screen
#Should probably change it to writing to another file,instead of displaying on terminal
def transposeDistributionFile_temp(outfile,x):
	fileN=outfile+'_ExpressionPatternRecovereyDistribution.bed'
	fileO=outfile+'tmp_ExpressionPatternRecovereyDistribution.bed'
	
	# Read all Tab-delimited rows from stdin.
	all_data = []
	with open(fileN,'r') as file:
		for line in file:
		    all_data.append(line.strip().split('\t'))
	file.close()
	
	# Transpose it.
	transposed = []
	for row_i in range(len(all_data[0])):
	    new_row = []
	    for col_i in range(len(all_data)):
	        new_row.append(all_data[col_i][row_i])
	    transposed.append(new_row)
	
	# Write it back out
	with open(fileO, 'w') as fileT:
		for row in transposed:
			line='\t'.join(row)
			fileT.write(line+'\n')	

#one set max FE and max spec result..printing on screen
	spec='percentageExpressionPatternPrecision'+str(t)	
	FE='foldEnriched'+str(t)
	df=pd.read_csv(fileO,delimiter='\t')
	#print(df)
	#print("following are the summary stats for your set"+x)
	#print("Percentage Precision of your set")
	#print(df[df['Tset_Name']==x][spec])
	
	maxSpec=df[spec].max()
	#print("Percentage Precision of max set")
	#print(df[df[spec]==maxSpec])
	
	
	#print("Fold Enrichment of your set")
	#print(df[df['Tset_Name']==x][FE])

	maxFE=df[FE].max()
	#print("Fold Enrichment of max set")
	#print(df[df[FE]==maxFE])
#----------------------------------------------------------------------------------------------------------------------------	
#This function is used to calculate fold enrichment values using the actual expression pattern recovery and expected/random 
def calculateFoldEnrichment(outfile,x,tab):
	setNum='mean'+str(t)
	setName=x+'_'+tab+':'
	if continuous=="True" or continuous=="TRUE" or continuous=="T":	
		setName=x+'_'+tab+'_'+str(t)+':'
	fileN=outfile+'_ExpressionPatternRecovereyDistribution.bed'

	df=pd.read_csv(fileN,sep='\t',header=(0))
	#df=df.set_index([0])
	df=df.set_index('Tset_Name')
	set=setName
	foldEnrichmentColName='foldEnriched'+str(t)
	dfEnrich=df.loc[df.index==set].divide(df.loc[df.index == setNum].values)
	dfEnrich=dfEnrich.rename({set:foldEnrichmentColName},axis='index')
	dfEnrich=dfEnrich.fillna(0)
	#df=df.append(dfEnrich)
	dfEnrich.to_csv('enriched.txt',sep='\t')
	with open('enriched.txt','r') as sumFile, open(outfile+'_ExpressionPatternRecovereyDistribution.bed','a') as outf:
		for line in sumFile:
		#if line.startswith('mean') or line.startswith('std') or line.startswith('min') or line.startswith('max'):
			if not line.startswith('Tset_Name'):
				outf.write(line)


#----------------------------------------------------------------------------------------------------------------------------
#This function is used to calculate maximum fold enrichment across all the data sets
def max_FE(outfile,x):
	
	fileN=outfile+'_ExpressionPatternRecovereyDistribution.bed'
	df=pd.read_csv(fileN,delimiter='\t')
		
	spec='percentageExpressionPatternPrecision'+str(t)	
	#find max 
	#print("following set has maximum percentage Precision")
	maxSpec=df[spec].max()
	#print(df[df[spec]==maxSpec])

	#print("following set has maximum fold enrichment")
	FE='foldEnriched'+str(t)
	maxFE=df[FE].max()
	#print(df[df[FE]==maxFE])

	
	
	
#----------------------------------------------------------------------------------------------------------------------------
#writing the output of evaluation to tab delimited file

def table_file(outfile,myList=[],*args):
	#print('in table file')
	with open(outfile+'.bed','a') as t:
		for items in myList:
			t.write(items+'\t')
		t.write('\n')
	
#################################################################--------------Main function----------------------####################################################################
def main():
	global subdirectory
	global d
	global d2
	global countd,countd_random
	global tsetBedOrNot
	global scrmshawOrNot
	global patternRecovery
	global t,x
	global nameOfmethod
	global continuous
	subdirectory='tmp' #Temporary directory
	t=1 #count of sets done
	global totalNumberOfCrmsKnownToCauseExpression
	totalNumberOfCrmsKnownToCauseExpression=0	
	#command line parsing files from users
	parser=argparse.ArgumentParser() 
	parser.add_argument('-nameOfmet','--nameOfmethod',help='if the CRM predictions were made through SCRMshaw or not',required=True)
	parser.add_argument('-predInBed','--predictionsInBed',help= 'Set of predictions in bed format')
	parser.add_argument('-fullredfly','--fullcrmFile',help='File of of all known crms: downloaded from redfly in bed format',required=True)
	parser.add_argument('-so','--scrmJoinedOutputFile',help='Scrmshaw Output file')
	parser.add_argument('-tset','--tsetBedOrNot',help='If training set information is in the bed format or not')
	parser.add_argument('-bedTset','--bedTsetFormat',help='File of training set CRMs in the bed format')
	parser.add_argument('-listTset','--tsetListFormat',help='File of training set CRMs in the list format')
	parser.add_argument('-pattern','--patternRecovery',help='if your input has nothing to do with how many expression paterrn crms recovered,set this parameter to False. Default value is True',default='True')
	parser.add_argument('-subsetcrmsExpBed','--subsetExpressionMappedCrmsBedFormat',help='Subset of redfly having only crms known to cause expression in Bed Format',required=True)
	parser.add_argument('-finalcrmsExp','--finalExpressionMappedCrms2015',help='Final expression file of 2015 ',required=True)
	parser.add_argument('-cont','--continuous',help='Do you want the output calculated in continuous fashion or not',default=True)
	parser.add_argument('-goodHits','--goodHitsOnly',help='Do you only want to see the evaluation of Top Hits above the inflevtion point or not, if yes specify number of hits for parse: 5000. if no, specify number of hits for parse',default=False)
	parser.add_argument('-p','--numParse',help='if the answer of goodHitsOnly is YES then you should set this to 5000 otherwise specify how much you want to parse your scrmshaw output to certain number of hits',required=True)
	parser.add_argument('-e','--exon',help='Exon file',required=True)	
	parser.add_argument('-drosog','--drosoGenome',help='Genome file',required=True)
	parser.add_argument('-s','--numberOfShuffles',help='You can give number of shuffles here.Higher the number,more time it will take. Default value is 20',default=20)
	parser.add_argument('-o','--outfile',help='Output file',required=True)
	parser.add_argument('-which','--whichSetToLookFor',help='your best guess about the set')
	args = parser.parse_args()
	
	nameOfmethod=args.nameOfmethod
	predictionsInBed=args.predictionsInBed
	scrmJoinedOutputFile=args.scrmJoinedOutputFile
	tsetBedOrNot=args.tsetBedOrNot
	bedTsetFormat=args.bedTsetFormat
	tsetListFormat=args.tsetListFormat
	fullcrmFile=args.fullcrmFile
	patternRecovery=args.patternRecovery
	numberOfShuffles=args.numberOfShuffles
	outfile=args.outfile
	exons=args.exon
	drosoGenome=args.drosoGenome
	numberOfShuffles=int(numberOfShuffles)
	subsetExpressionMappedCrmsBedFormat=args.subsetExpressionMappedCrmsBedFormat
	finalExpressionMappedCrms2015=args.finalExpressionMappedCrms2015
	continuous=args.continuous
	numParse=args.numParse
	goodHitsOnly= args.goodHitsOnly
	whichSetToLookFor=args.whichSetToLookFor
	outfile_random=outfile+'_random'
	
	fullCrms=os.path.abspath(fullcrmFile)
	expressionMappedCrmsBedFormat=os.path.abspath(subsetExpressionMappedCrmsBedFormat)
	
	#creating temporary directory to move stuff there
	my_path=os.getcwd()
	if not os.path.isdir(my_path+'/'+subdirectory):
		os.makedirs(my_path+'/'+subdirectory)
	
	#creating dictionary from the file that  contain expression-mapped CRMs, use to assess training set Precision or expression pattern recovery
	
	path_to_subset_crms=os.path.abspath(finalExpressionMappedCrms2015)
	with open(path_to_subset_crms) as fi:
		rows = (line.split('\t') for line in fi )
		d2 = {row[0].strip(':\:'):row[1:] for row in rows }
		
	#creating dictionary for keeping track of counters of each group to calculate the distribution of expression groups
	with open(path_to_subset_crms) as fi:
		rows = (line.split('\t') for line in fi )
		countd={row5[0].strip(':\:'):0 for row5 in rows}
	#for random data	
	with open(path_to_subset_crms) as fi:
		rows = (line.split('\t') for line in fi )
		countd_random={row5[0].strip(':\:'):0 for row5 in rows}
	
	print("check!")
	#counting all the crms in the expression mapped crms files dictionary
	for keys in d2:
		totalNumberOfCrmsKnownToCauseExpression=totalNumberOfCrmsKnownToCauseExpression+len(d2[keys])
		#print(keys+" values "+str(len(d2[keys])))		

	#Creating a file..............
	fileE=open(outfile+'_ExpressionPatternRecovereyDistribution.bed','w')

	#writing header string to expression distribution files
	stringOfGroupNames='Tset_Name'+'\t'
	for keys in countd:
		stringOfGroupNames=stringOfGroupNames+keys+'\t'	
		
	stringOfSizes='Size'+'\t'
	for key in d2:
		stringOfSizes=stringOfSizes+str(len(d2[key]))+'\t'
	with open(outfile+'_ExpressionPatternRecovereyDistribution.bed','a') as e:
		e.write(stringOfGroupNames+'\n')
		e.write(stringOfSizes+'\n')
	
	headingList=[]
	
	#if this is true it will parse SCRMshaw joined output to each set and each method and loop through all of them to get the individual set's evaluation result
	if nameOfmethod=="SCRMshaw" or nameOfmethod=="scrmshaw" or nameOfmethod=="SCRMSHAW" :
		####creating necessary files for output
		if patternRecovery=='True' or patternRecovery =="T" or patternRecovery == "TRUE" :
			fin2=open(outfile+'.bed','w')
			headingList.extend(["TsetName","Method","TsetSize","TotalREDfly","ModifiedREDfly","TotalREDfly2010","ModifiedREDfly2010","SCRMs","TrainingSetRecovered","PercentageTrainingSetSensitivity","ExpectedValueRandom","PercentageExpectedTrainingSetSensitivity","DifferencesInPercentagesActual&Expected","P_value","REDflyRecovered","PercentageRedflyRecovered","ExpectedOverlapRandom","PercentageExpectedOverlap","DifferencesInPercentagesActual&Expected","P_value","PercentageRecallREDflyRecovery","numberOfKnownCrmsCausesExpressionInTset","numOfSCRMsRecoveredExpMappedCrms","PercentagePatternRecovery","RandomREDflyRecoveredPattern","PercentageOfRandomPatternRecovered","numOfRecoveredSCrmsBelongingToOwnGroup","numOfRecoveredSCrmsBelongingToAnyGroup","percentageExpressionPatternPrecision","ExpectedpercentageExpressionPatternPrecision","P_value","DifferencesInPercentagesActual&ExpectedExpressionRecovered","percentageExpressionPatternRecall"])
			with open(outfile+'.bed','a') as h:
				for things in headingList:
					h.write(things+'\t')
				h.write('\n')
		else:
			fin2=open(outfile+'.bed','w')
			headingList.extend(["TsetName","Method","TsetSize","TotalREDfly","ModifiedREDfly","SCRMs","TrainingSetRecovered","PercentageTrainingSetSensitivity","ExpectedValueRandom","PercentageExpectedTrainingSetSensitivity","DifferencesInPercentagesActual&Expected","P_value","REDflyRecovered","PercentageRedflyRecovered","ExpectedOverlapRandom","PercentageExpectedOverlap","DifferencesInPercentagesActual&Expected","P_value","PercentageRecallREDflyRecovery"])
			with open(outfile+'.bed','a') as h:
				for things in headingList:
					h.write(things+'\t')
				h.write('\n')
				
		
						
		methods=['imm','hexmcd','pac']
		
		#creating dictionary from the file that contain names of all crms against each training set; use to assess training set sensitivity; in the format name of sets being keys and their respective crms as key's values
		path_to_known_crms=os.path.abspath(tsetListFormat) 
		with open(path_to_known_crms) as fin:
			rows = (line.split('\t') for line in fin )
			d = {row[0].strip(':\:'):row[1:] for row in rows }			
		

		#Creating list to iterate through three methods
		methods_val=[None,None,None]
		
		#Parsing the output file into separate files for each training set and each method via creating three dictionaries for each method: keys being the names of data sets associated with that method 
		scrmJoinedOutputFile=os.path.abspath(scrmJoinedOutputFile) 
		methods_val[0],methods_val[1],methods_val[2]=parse_output(scrmJoinedOutputFile,35000)
		
		#This loop is used to iterate through all three method's dictionary:
		for num in range(len(methods)):
			print("Now method:"+methods[num])
			print("sets in this method "+str(methods_val[num]))
			#This loop is used to iterate through all the training sets in the given methods dictionary
			for x in methods_val[num]:
				if x in d.keys():
					#Getting tmp path of methods dictionary file
					for root, dirs, files in os.walk(os.getcwd()):
						for name in files:
							if name==methods[num]+'_'+x+'_fullLength.bed':
								unsorted_indScrmOutputfile=os.path.abspath(os.path.join(root,name))
							
							
					#if user wants to evaluate the result in a continuous fashion, like evaluate at every 500 predcitions..
					if continuous=="True" or continuous=="TRUE" or continuous=="T":		
						numberparse=250
					else:
						numberparse=int(numParse)
				
					#if user wants to evaluate only top good hits based on the composite curve method 
					if goodHitsOnly == "True" or goodHitsOnly== "true" or goodHitsOnly== "TRUE" or goodHitsOnly =="T":
						numParse,numScore,numAmp = topGoodPrediction(unsorted_indScrmOutputfile)
						if continuous=="False" or continuous=="FALSE" or continuous=="F":
							numberparse=numParse
							print("Top pred after first(amplitude cutoff): " + str(numParse))
						
					else:
						numParse=int(numParse)
						print("user specified hits" + str(numParse))
				
					tab1=methods[num]								
					#if user has set it up to find the continuous evaluation it will start fromm 100 and go up to max otherwise it would run only one time
					while numberparse <= numParse:
						print("Numparse "+str(numberparse)+" parse "+str(numParse))
						
						#creating files for keeping track of distribution of expression recovered crms 
						fileE=open(outfile_random+'_ExpressionPatternRecovereyDistribution.bed','w')
						with open(outfile_random+'_ExpressionPatternRecovereyDistribution.bed','a') as e:
							e.write(stringOfSizes+'\n')
							e.write(stringOfGroupNames+'\n')					
					
						#extract number of lines of peaks given by the user
						noOflinesExtracted,unsorted_indScrmOutputfile=parse_output_individual(unsorted_indScrmOutputfile,numberparse,x,tab1)
						#if user wants to evaluate top hits only then extract those and discard the rest of the data
						if goodHitsOnly == "True" or goodHitsOnly== "true" or goodHitsOnly== "TRUE" or goodHitsOnly =="T":
							noOflinesExtracted,unsorted_indScrmOutputfile=extract_below_cutoff_vals(unsorted_indScrmOutputfile,numberparse,x,tab1,numScore,numAmp)
							print('Number of Top Predictions: '+str(noOflinesExtracted))
						
						#sorting and merging the each unique scrmshaw output file (the files that have been created earlier by d_imm dictionary for each different data set > probably unsorted.)
						sortedMergedScrmsFilePath,numOfsortedMergedPredictedCrms=sort_and_merge_output(x,tab1,unsorted_indScrmOutputfile)	
				
						#creating a modified list of redfly files(FULLCRMS and subsetExpressionMappedCrmsBedFormat) will be called as "modifiedCrms2010 and ModifiedCrms" which will not have the crms that were used as training set for scrmshaw
	
						#excluding training set crms from redfly file and expression mapped file..based on names...will have modified and excluded file 
						totalKcrms,countModifiedCrms,excludedCrmsCount,modifiedCrmsFilePath,excludedCrmsFilePath= exclude_training_set_crms("redfly",x,tab1,x,fullCrms)
						totalKSubsetcrms,countSubsetModifiedCrms,excludedSubsetCrmsCount,modifiedSubsetCrmsFilePath,excludedSubsetCrmsFilePath= exclude_training_set_crms("expressionMapped",x,tab1,x,expressionMappedCrmsBedFormat)
							
						#Total known crms in both redfly files and expression mapped file after sorting and merging
						size_of_sortedMergedCrms=sort_merge_redfly_find_size("redfly",fullCrms,x,tab1)
						size_of_sortedMergedSubsetCrms=sort_merge_redfly_find_size("expressionMapped",expressionMappedCrmsBedFormat,x,tab1)

						#converting to bed format
						sortedMergedScrmsBed=bed_conversion(sortedMergedScrmsFilePath)
						excludedCrmsFilePathBed=bed_conversion(excludedCrmsFilePath)
						
																				############THREE MEASURES OF EVALUATION###############
						
						#`-------------------------------------	1: Tset Sensitivity---------------------------
						print("Calculating Training set Sensitivity")
						#sorting excluded crms
						sortedExcludedCrms=sort_and_merge_tset(excludedCrmsFilePathBed,x,tab1)		
						#BED conversion
						sortedExcludedCrmsBed=bed_conversion(sortedExcludedCrms)
						countCommonScrmAndExcluded=intersect_Scrm_and_Excluded_crms(x,tab1,sortedMergedScrmsBed,sortedExcludedCrmsBed)
						#print("TsetRecovered "+str(countCommonScrmAndExcluded))	
						#percentage of sensitivity of result
						#percentageOfSensitivity="{0:.2f}%".format((countCommonScrmAndExcluded/excludedCrmsCount)*100) #or  equal print(len(d[x]))
						percentageOfSensitivity=countCommonScrmAndExcluded/excludedCrmsCount
						#print("Percent True Training set Recovered "+str(percentageOfSensitivity))	
			
						#calculating expected Training set Recovery
						ex=os.path.abspath(exons)
						g=os.path.abspath(drosoGenome)
						#Generating random data from genome and finding Significance of SCRMSHAW training set sensitivity
						mean2V,med2V,sd2V,min2Val,max2Val,stats2,z_value2,p_value2=shuffle(sortedExcludedCrmsBed,sortedMergedScrmsBed,ex,g,countCommonScrmAndExcluded,int(numberOfShuffles))
						#percentage of expected sensitivity of result
						
						percentageOfExpectedSensititvity=(mean2V/excludedCrmsCount)
						#calculating differences in percentages of sensitivity of scrms vs expected 
						#differenceInPercentagesTsetSensitivity= "{0:.2f}".format(((countCommonScrmAndExcluded/excludedCrmsCount)*100)-((mean2V/excludedCrmsCount)*100))
						differenceInPercentagesTsetSensitivity= (countCommonScrmAndExcluded/excludedCrmsCount)-((mean2V/excludedCrmsCount))
						#print("Percent Random Training set Recovered "+str(percentageOfExpectedSensititvity))	
				
						#	--------------------------------2: REDfly Recovery ----------------------------
						print("Calculating REDfly Recovery")
						#sorting and merging the files of modified crm(redfly files excluding tsets crms)
						sortedMergedModifiedCrmsFilePath,sizeOfSortedMergedModifiedCrms=sort_and_merge_modifiedCrms("redfly",modifiedCrmsFilePath,x,tab1)
							
						#converting to bed format
						sortedMergedModifiedCrmsBed=bed_conversion(sortedMergedModifiedCrmsFilePath)
					
						#Finding the "Number of HITS" (common crms between modified crms and the individual parsed sorted+merged scrmshaw output)
						no_of_overlaps2,pathOfIntersectResult=intersect_scrms_and_modifiedCrms("redfly",x,tab1,sortedMergedModifiedCrmsBed,sortedMergedScrmsBed)
						
						#for calculating number of redfly recovered 
						#no_of_overlaps=redfly_recovered(pathOfIntersectResult)
						pathOfIntersectResultBed=bed_conversion(pathOfIntersectResult)
						no_of_overlaps=pathOfIntersectResultBed.count()
						
						
						percentageOfOverlaps=(no_of_overlaps/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded))
						#print("Percent True Redfly Recovered =: "+str(percentageOfOverlaps))
						
						print("Calculating REDfly Recovery Recall")
						#Recall/Sensitivity of REDfly recovery TP/TP+FN .. True positive=SCRMs recovered REDfly CRMs, False Negative= Total REDfly CRMs - True Positives
						percentageRecallREDfly=(no_of_overlaps/(no_of_overlaps+(sizeOfSortedMergedModifiedCrms-no_of_overlaps)))
						
						##calculating expected REDfly Recovery
						#Generating random data from genome and finding Significance of SCRMSHAW no of hits
						meanV,medV,sdV,minVal,maxVal,stats,z_value,p_value=shuffle(sortedMergedModifiedCrmsBed,sortedMergedScrmsBed,ex,g,no_of_overlaps,int(numberOfShuffles))
							
						#print("Exp Redfly Recovered =: "+str(meanV))		
						#Percentage of Expected redfly recovered
						percentageOfExpectedOverlaps=(meanV/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded))
						#print("Percent Random Expectation Redfly Recovered =: "+str(percentageOfExpectedOverlaps))	
						
						#Difference in percentages of scrms and expected redfly recovered
						differenceInPercentagesRedflyRecovered=((no_of_overlaps/numOfsortedMergedPredictedCrms))-((meanV/numOfsortedMergedPredictedCrms))		
						
						
						#	--------------------------------3:Expression Pattern Precision-------------------------------------
						print("Calculating Expression pattern Precision and Recall")
						#basically repeat all the steps that we did with actual redfly like extract tset from the prediction file and then intersect it ..etc
						
						#sorting and merging the files of modified expression mapped crm(expression mapped redfly files excluding tsets crms)
						#sortedMergedModifiedSubsetCrmsFilePath,sizeOfSortedMergedModifiedSubsetCrms=sort_and_merge_modifiedCrms("expressionMapped",modifiedSubsetCrmsFilePath,x,tab1,)
						sortedMergedModifiedSubsetCrmsFilePath,sizeOfSortedMergedModifiedSubsetCrms=sort_expressionMappedCrms("expressionMapped",modifiedSubsetCrmsFilePath,x,tab1,)
						sortedMergedModifiedSubsetCrmsBed=bed_conversion(sortedMergedModifiedSubsetCrmsFilePath)

						
															##-----------expression pattern Recovery
						#for finding out no of hits matched to expected expression we need to intersect file of subset(because that is the file of those crms having some expression)
						no_of_overlapsSubset,pathOfSubsetIntersectResult= intersect_scrms_and_modifiedCrms("expressionMapped",x,tab1,sortedMergedModifiedSubsetCrmsBed,sortedMergedScrmsBed)
						pathOfSubsetIntersectResultBed=bed_conversion(pathOfSubsetIntersectResult)
						no_of_overlapsX=pathOfSubsetIntersectResultBed.count()
						#print("Expression Pattern Recovered =: "+str(no_of_overlapsSubset))							
						#print("recovery defined by : num of scrms recovered any crms from redfly expr mapped file divide by total num of predictions excl ones used for tset")
						percentageOfHitsRecoveredExpressionPattern=no_of_overlapsSubset/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded)
						#print("Percent Expression Pattern Recovered =: "+str(percentageOfHitsRecoveredExpressionPattern))	
						
						#Generating random data from genome and finding Significance of SCRMSHAW no of hits
						meanExp,medExp,sdExp,minVExp,maxVExp,statsExp,zExp,pExp=shuffle(sortedMergedModifiedSubsetCrmsBed,sortedMergedScrmsBed,ex,g,no_of_overlapsX,int(numberOfShuffles))
												
						#Percentage of Expected expression pattern recovered
						#print("Expected Expression Pattern Recovered =: "+str(meanExp))
						percentageOfExpectedHitsRecoveredExpressionPattern=meanExp/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded)
						#print("Percent Expected Expression Pattern Recovered =: "+ str(percentageOfExpectedHitsRecoveredExpressionPattern))	

						
															#---------Overlap with shared annotation (Precision)
						#calculating expected Precision
						# getting number of crms predicted by scrmshaw that matches with expected expression pattern
						hitsMatchedExpression,numberOfKnownCrmsCausesExpressionInTset,numberOfAnyPatternRecovered=parse_intersected_scrmshaw(outfile,x,pathOfSubsetIntersectResult,tab1)
						
						
						#Expression Pattern Precision: TP/TP+FP..... True Positives= SCRMs with correct expression pattern...False Positive=SCRMs with incorrect expression pattern (SCRMs with any pattern - TP)
						percentageExpressionPatternPrecision=hitsMatchedExpression/numberOfAnyPatternRecovered
						#print("Percent Precision =: "+str(percentageExpressionPatternPrecision))	
						
						#Expression Pattern Recall/Sensitivity: TP/TP+FN ..... True Positives= SCRMs with correct expression pattern..False Negative= Total Known CRMs with Expression Pattern in Training set - TP
						percentageExpressionPatternRecall= hitsMatchedExpression/(hitsMatchedExpression+(numberOfKnownCrmsCausesExpressionInTset-hitsMatchedExpression))
						
						#trying to see the expected distribtuion across all the sets..
						#print("Random Expectation Values")
						percentageExpressionPatternPrecisionExp=[]
						for i in range(numberOfShuffles):
							shuffled2=sortedMergedScrmsBed.shuffle(excl=ex, noOverlapping=True, g=g).sort().saveas('shuffled')	
							shuffledpath=os.path.abspath('shuffled')
							shuffledBed=bed_conversion(shuffledpath)		
							no_of_overlapsR,pathOfIntersectResultR=intersect_scrms_and_modifiedCrms("expressionMapped",x,tab1,sortedMergedModifiedSubsetCrmsBed,shuffledBed)		
							
							
							hitsMatchedExpressionExp,numberOfKnownCrmsCausesExpressionInTsetExp,numberOfAnyPatternRecoveredExp=parse_intersected_scrmshaw(outfile_random,x,pathOfIntersectResultR,tab1)
							percentageExpressionPatternPrecisionExp.append(hitsMatchedExpressionExp/numberOfAnyPatternRecoveredExp)
							
						#calculating summary stats for the expected distribution result
					
						#calculating p values for Precision 
						#get summary stats 
						maxRspec=max(percentageExpressionPatternPrecisionExp)
						minRspec=min(percentageExpressionPatternPrecisionExp)
						meanRspec=statistics.mean(percentageExpressionPatternPrecisionExp)
						medianRspec=statistics.median(percentageExpressionPatternPrecisionExp)
						sdRspec=statistics.stdev(percentageExpressionPatternPrecisionExp)
						#print("Median Percent Precision Expected =: "+str(medianRspec))	
						summary_stats="Mean= "+ str(meanRspec)+"\nMedian= "+str(medianRspec)+"\nStandard Deviation= "+str(sdRspec)+"\nMinimum value= "+str(minRspec)+"\nMaximum Value= "+str(maxRspec)

						#from a_and_b_count and mean sd of above array
						#print(a_and_b_count)
						try:	
							zRspec=float((percentageExpressionPatternPrecision-meanRspec)/sdRspec)
							#print(zRspec)
							#pRspec=stats.norm.sf(float(abs(zRspec)))
							pRspec= 1 - scipy.special.ndtr(zRspec)
							#print(pRspec)
						except ZeroDivisionError:
							np.seterr(divide='ignore', invalid='ignore')
							pRspec=0
						#return mean,median,sd,minV,maxV,summary_stats,z,p
						
						differenceInPercentagesExpressionPatternRecovered=percentageExpressionPatternPrecision-meanRspec
						
						### Following lines of code calculates and write fold enrichment along with the Precision of individual sets in the Expression Pattern Distribtuiton File
						
						calculating_sum_stats_random_file(outfile)
						#calculating fold enrichment 
						calculateFoldEnrichment(outfile,x,tab1)
						
						#temporary transpose and displaying on the terminal
						transposeDistributionFile_temp(outfile,x)
						#max_FE(outfile,x)

						#writing to files depending upon _condition applied depending upon if Expression Pattern Recovery step is needed or not
						printList=[]
						
						if patternRecovery=='True' or  patternRecovery == "T" or patternRecovery=="TRUE":
							printList.extend([x,tab1,str(excludedCrmsCount),str(size_of_sortedMergedSubsetCrms),str(sizeOfSortedMergedModifiedSubsetCrms),str(size_of_sortedMergedCrms),str(sizeOfSortedMergedModifiedCrms),str(numOfsortedMergedPredictedCrms),str(countCommonScrmAndExcluded),str(percentageOfSensitivity),str(mean2V),str(percentageOfExpectedSensititvity),str(differenceInPercentagesTsetSensitivity),str(p_value2),str(no_of_overlaps),str(percentageOfOverlaps),str(meanV),str(percentageOfExpectedOverlaps),str(differenceInPercentagesRedflyRecovered),str(p_value),str(percentageRecallREDfly),str(numberOfKnownCrmsCausesExpressionInTset),str(no_of_overlapsSubset),str(percentageOfHitsRecoveredExpressionPattern),str(meanExp),str(percentageOfExpectedHitsRecoveredExpressionPattern),str(hitsMatchedExpression),str(numberOfAnyPatternRecovered),str(percentageExpressionPatternPrecision),str(meanRspec),str(pRspec),str(differenceInPercentagesExpressionPatternRecovered),str(percentageExpressionPatternRecall)]) 
							#table_file(outfile,printList)
							#table_file(outfile,x,tab1,str(excludedCrmsCount),str(size_of_sortedMergedSubsetCrms),str(sizeOfSortedMergedModifiedSubsetCrms),str(size_of_sortedMergedCrms),str(sizeOfSortedMergedModifiedCrms),str(numOfsortedMergedPredictedCrms),str(countCommonScrmAndExcluded),str(percentageOfSensitivity),str(mean2V),str(percentageOfExpectedSensititvity),str(differenceInPercentagesTsetSensitivity),str(p_value2),str(no_of_overlaps),str(percentageOfOverlaps),str(meanV),str(percentageOfExpectedOverlaps),str(differenceInPercentagesRedflyRecovered),str(p_value),str(numberOfKnownCrmsCausesExpressionInTset),str(no_of_overlapsSubset),str(percentageOfHitsRecoveredExpressionPattern),str(meanExp),str(percentageOfExpectedHitsRecoveredExpressionPattern),str(hitsMatchedExpression),str(numberOfAnyPatternRecovered),str(percentageExpressionPatternPrecision),str(meanRspec),str(pRspec),str(differenceInPercentagesExpressionPatternRecovered))						
						else:
							printList.extend([x,tab1,str(excludedCrmsCount),str(size_of_sortedMergedCrms),str(sizeOfSortedMergedModifiedCrms),str(numOfsortedMergedPredictedCrms),str(countCommonScrmAndExcluded),str(percentageOfSensitivity),str(mean2V),str(percentageOfExpectedSensititvity),str(differenceInPercentagesTsetSensitivity),str(p_value2),str(no_of_overlaps),str(percentageOfOverlaps),str(meanV),str(percentageOfExpectedOverlaps),str(differenceInPercentagesRedflyRecovered),str(p_value),str(percentageRecallREDfly)])
							#table_file_excludingExpressionRecovered(outfile,printList)
							#table_file_excludingExpressionRecovered(outfile,x,tab1,str(excludedCrmsCount),str(size_of_sortedMergedCrms),str(sizeOfSortedMergedModifiedCrms),str(numOfsortedMergedPredictedCrms),str(countCommonScrmAndExcluded),str(percentageOfSensitivity),str(mean2V),str(percentageOfExpectedSensititvity),str(differenceInPercentagesTsetSensitivity),str(p_value2),str(no_of_overlaps),str(percentageOfOverlaps),str(meanV),str(percentageOfExpectedOverlaps),str(differenceInPercentagesRedflyRecovered),str(p_value))

						table_file(outfile,printList)											
						print("Done",t,"set of length",numberparse)
						t=t+1
						if continuous=="True" or continuous=="TRUE" or continuous=="T":	
							print("adding 250 to numparse")
						numberparse=numberparse+250
		#final transposing of the file
		transposeDistributionFile(outfile)
		
		print("Done")
		
#-------------------------------------------------------------------------------------------------------------------		
	#if the algorithm used for prediction is not SCRMshaw then it would go through these steps
	else:
		headingList=[]
		#files processing
		#Creating Excel (Single tab delimited) file to work_ condition applied depending upon if Expression Pattern Recovery step is needed or not
		####creating necessary files for output
		if patternRecovery=='True' or patternRecovery =="T" or patternRecovery == "TRUE" :
			fin2=open(outfile+'.bed','w')
			headingList.extend(["TsetName","Method","TsetSize","TotalREDfly","ModifiedREDfly","TotalREDfly2010","ModifiedREDfly2010","SCRMs","TrainingSetRecovered","PercentageTrainingSetSensitivity","ExpectedValueRandom","PercentageExpectedTrainingSetSensitivity","DifferencesInPercentagesActual&Expected","P_value","REDflyRecovered","PercentageRedflyRecovered","ExpectedOverlapRandom","PercentageExpectedOverlap","DifferencesInPercentagesActual&Expected","P_value","PercentageRecallREDflyRecovery","numOfSCRMsRecoveredExpMappedCrms","PercentagePatternRecovery","RandomREDflyRecoveredPattern","PercentageOfRandomPatternRecovered","numOfRecoveredSCrmsBelongingToOwnGroup","numOfRecoveredSCrmsBelongingToAnyGroup","percentageExpressionPatternPrecision","ExpectedpercentageExpressionPatternPrecision","P_value","DifferencesInPercentagesActual&ExpectedExpressionRecovered","percentageExpressionPatternRecall"])
			with open(outfile+'.bed','a') as h:
				for things in headingList:
					h.write(things+'\t')
				h.write('\n')
			
		else:
			fin2=open(outfile+'.bed','w')
			headingList.extend(["TsetName","Method","TsetSize","TotalREDfly","ModifiedREDfly","SCRMs","TrainingSetRecovered","PercentageTrainingSetSensitivity","ExpectedValueRandom","PercentageExpectedTrainingSetSensitivity","DifferencesInPercentagesActual&Expected","P_value","REDflyRecovered","PercentageRedflyRecovered","ExpectedOverlapRandom","PercentageExpectedOverlap","DifferencesInPercentagesActual&Expected","P_value","PercentageRecallREDflyRecovery"])
			with open(outfile+'.bed','a') as h:
				for things in headingList:
					h.write(things+'\t')
				h.write('\n')
			
		
		
		#checking if user wants to evaluate output in continuous fashion meaning top 100 predictions evaluated first then top 200 and so on
		
		if continuous=="True" or continuous=="TRUE" or continuous=="T":		
			numberparse=250
		else:
			numberparse=int(numParse)
		
		
		x='query'
		tab1=nameOfmethod
		print("Numparse"+str(numberparse)+"parse"+str(numParse))
		#if user has set it up to find the continuous evaluation it will start fromm 100 and go up to max otherwise it would run only one time
		while(int(numberparse) <= int(numParse)):
		
			print("Numparse "+str(numberparse)+" parse "+str(numParse))
		
			#creating file for distribution of random shuffle results
			fileE=open(outfile_random+'_ExpressionPatternRecovereyDistribution.bed','w')
			with open(outfile_random+'_ExpressionPatternRecovereyDistribution.bed','a') as e:
				e.write(stringOfSizes+'\n')
				e.write(stringOfGroupNames+'\n')
			#1: PREDICTION post processing
			unsorted_indScrmOutputfile=predictionsInBed
			noOflinesExtracted,unsorted_indScrmOutputfile=parse_output_individual(unsorted_indScrmOutputfile,numberparse,x,tab1)
			
			#sorting and merging the each unique scrmshaw output file (the files that have been created earlier by d_imm dictionary for each different data set > probably unsorted.)
			sortedMergedScrmsFilePath,numOfsortedMergedPredictedCrms=sort_and_merge_output(x,tab1,unsorted_indScrmOutputfile)
			
			#converting to bed format
			sortedMergedScrmsBed=bed_conversion(sortedMergedScrmsFilePath)
			
																			############THREE MEASURES OF EVALUATION###############
			
			#`-------------------------------------	1: Tset Sensitivity---------------------------
			print("Calculating Training Sensitivity")
			#converting REDfly bed file to BED format
			fullredflyBED=bed_conversion(fullcrmFile)
			expressionMappedCrmsBedFormatBED=bed_conversion(expressionMappedCrmsBedFormat)
		
			#if training set provided is not in bed then extract tset coordinates information from REDfly based on the names in the list provided through creating dictionary
			if tsetBedOrNot=="False" or tsetBedOrNot=="F" or tsetBedOrNot=="FALSE":
				#find out method and set name from the user 
				x='mapping1.adult_mesoderm'
				tab1='imm'

				#creating dictionary from the file that contain names of all crms against each training set; use to assess training set sensitivity; in the format name of sets being keys and their respective crms as key's values
				path_to_known_crms=os.path.abspath(tsetListFormat) 
				with open(path_to_known_crms) as fin:
					rows = (line.split('\t') for line in fin )
					d = {row[0].strip(':\:'):row[1:] for row in rows }					
				#excluding TSET crms from redfly using names provided in the list(dictionary d)	
				totalKcrms,countModifiedCrms,excludedCrmsCount,modifiedCrmsFilePath,excludedCrmsFilePath= exclude_training_set_crms("redfly",x,tab1,x,fullcrmFile)
				excludedCrmsFilePathBed=bed_conversion(excludedCrmsFilePath)
			
				###excluding TSET crms from expression mapped crms in bed using names provided in the list(dictionary d)	
				totalKSubsetcrms,countSubsetModifiedCrms,excludedSubsetCrmsCount,modifiedSubsetCrmsFilePath,excludedSubsetCrmsFilePath= exclude_training_set_crms("expressionMapped",x,tab1,x,expressionMappedCrmsBedFormat)
			
			#if the training set is in bed format then using its coordinate information for evaluation purpose, e.g to find out sensitivity and recovery etc
			else: 
				excludedCrmsFilePath=bedTsetFormat
				excludedCrmsFilePathBed=bed_conversion(excludedCrmsFilePath)
				#size of tset
				excludedCrmsCount=excludedCrmsFilePathBed.count()
			
				#expressionMapped
				excludedSubsetCrmsFilePath=bedTsetFormat
				excludedSubsetCrmsFilePathBed=bed_conversion(excludedSubsetCrmsFilePath)
		
				#modified crms = Redfly pybed - tset crms
				modifiedCrmsFilePath=exclude_tsetBedFormat('redfly',x,tab1,excludedCrmsFilePathBed,fullredflyBED)
				#modified subset crms=expression mapped pybed - tset crms
				modifiedSubsetCrmsFilePath=exclude_tsetBedFormat('expressionMapped',x,tab1,excludedSubsetCrmsFilePathBed,expressionMappedCrmsBedFormatBED)

			#sorting and merging training set file , but before that needs BED version of that file
			#sorting excluded crms
			sortedExcludedCrms=sort_and_merge_tset(excludedCrmsFilePathBed,x,tab1)		
			#BED conversion
			sortedExcludedCrmsBed=bed_conversion(sortedExcludedCrms)
			#counting Number of Tsets getting back
			countCommonScrmAndExcluded=intersect_Scrm_and_Excluded_crms(x,tab1,sortedMergedScrmsBed,sortedExcludedCrmsBed)
			#print("TsetRecovered "+str(countCommonScrmAndExcluded))	
			#percentage of sensitivity of result
			
			percentageOfSensitivity=countCommonScrmAndExcluded/excludedCrmsCount
			#print("PercentTsetRecovered "+str(percentageOfSensitivity))	
		
			#calculating expected Training set Sensitivity
			ex=os.path.abspath(exons)
			g=os.path.abspath(drosoGenome)
			#Generating random data from genome and finding Significance of SCRMSHAW training set sensitivity
			mean2V,med2V,sd2V,min2Val,max2Val,stats2,z_value2,p_value2=shuffle(sortedExcludedCrmsBed,sortedMergedScrmsBed,ex,g,countCommonScrmAndExcluded,int(numberOfShuffles))
		
			#percentage of expected sensitivity of result
			#percentageOfExpectedSensititvity="{0:.2f}%".format((mean2V/excludedCrmsCount)*100)
			percentageOfExpectedSensititvity=(mean2V/excludedCrmsCount)
			#calculating differences in percentages of sensitivity of scrms vs expected 
			differenceInPercentagesTsetSensitivity= (countCommonScrmAndExcluded/excludedCrmsCount)-((mean2V/excludedCrmsCount))
			#print("Percent Exp TsetRecovered "+str(percentageOfExpectedSensititvity))	
		
			#	--------------------------------2: REDfly Recovery ----------------------------
			print("Calculating REDfly Recovery")
			#Size of the files
			#Total known crms in both redfly files(2010 and expression mapped file) after sorting and merging
			size_of_sortedMergedCrms=sort_merge_redfly_find_size("redfly",fullCrms,x,tab1)
			size_of_sortedMergedSubsetCrms=sort_merge_redfly_find_size("expressionMapped",expressionMappedCrmsBedFormat,x,tab1)
	
			#sorting and merging the files of modified crm and modified subset crm(redfly files excluding tsets crms)
			sortedMergedModifiedCrmsFilePath,sizeOfSortedMergedModifiedCrms=sort_and_merge_modifiedCrms("redfly",modifiedCrmsFilePath,x,tab1)
			
			#converting to bed format
			sortedMergedModifiedCrmsBed=bed_conversion(sortedMergedModifiedCrmsFilePath)

	
			#Finding the "Number of HITS" (common crms between modified crms and the individual parsed sorted+merged scrmshaw output)
			no_of_overlaps2,pathOfIntersectResult=intersect_scrms_and_modifiedCrms("redfly",x,tab1,sortedMergedModifiedCrmsBed,sortedMergedScrmsBed)
		
			#for calculating number of redfly recovered 
			#no_of_overlaps=redfly_recovered(pathOfIntersectResult)
			pathOfIntersectResultBed=bed_conversion(pathOfIntersectResult)
			no_of_overlaps=pathOfIntersectResultBed.count()
		
			#print("Redfly Recovered =: "+str(no_of_overlaps))		
			#Percentage of redfly recovered....moved after calculating Tset Sensitivity 
			
			percentageOfOverlaps=(no_of_overlaps/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded))
			#print("Perc Redfly Recovered =: "+str(percentageOfOverlaps))	
			print("Calculating REDfly Recovery Recall")
			#Recall/Sensitivity of REDfly recovery TP/TP+FN .. True positive=SCRMs recovered REDfly CRMs, False Negative= Total REDfly CRMs - True Positives
			percentageRecallREDfly=(no_of_overlaps/(no_of_overlaps+(sizeOfSortedMergedModifiedCrms-no_of_overlaps)))
		
			#calculating expected REDfly Recovery
			#Generating random data from genome and finding Significance of SCRMSHAW no of hits
			meanV,medV,sdV,minVal,maxVal,stats,z_value,p_value=shuffle(sortedMergedModifiedCrmsBed,sortedMergedScrmsBed,ex,g,no_of_overlaps,int(numberOfShuffles))
			
			#print("Exp Redfly Recovered =: "+str(meanV))		
			#Percentage of Expected redfly recovered
			
			percentageOfExpectedOverlaps=(meanV/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded))
			#print("Perc Expected Redfly Recovered =: "+str(percentageOfExpectedOverlaps))	
		
			#Difference in percentages of scrms and expected redfly recovered
			differenceInPercentagesRedflyRecovered=((no_of_overlaps/numOfsortedMergedPredictedCrms))-((meanV/numOfsortedMergedPredictedCrms))		
		

		
			#	--------------------------------3:Expression Pattern Precision-------------------------------------
			print("Calculating Expression pattern Precision and Recall")
			#basically repeat all the steps that we did with actual redfly like extract tset from the prediction file and then intersect it ..etc
		
			##expression pattern recovery
			sortedMergedModifiedSubsetCrmsFilePath,sizeOfSortedMergedModifiedSubsetCrms=sort_expressionMappedCrms("expressionMapped",modifiedSubsetCrmsFilePath,whichSetToLookFor,tab1)
			sortedMergedModifiedSubsetCrmsBed=bed_conversion(sortedMergedModifiedSubsetCrmsFilePath)		
			#for finding out no of hits matched to expected expression we need to intersect file of subset(because that is the file of those crms having some expression)
			no_of_overlapsSubset,pathOfSubsetIntersectResult=intersect_scrms_and_modifiedCrms("expressionMapped",whichSetToLookFor,tab1,sortedMergedModifiedSubsetCrmsBed,sortedMergedScrmsBed)
			pathOfSubsetIntersectResultBed=bed_conversion(pathOfSubsetIntersectResult)
			no_of_overlapsX=pathOfSubsetIntersectResultBed.count()
			#print("Expression Pattern Recovered =: "+str(no_of_overlapsSubset))							
			#print("recovery defined by : num of scrms recovered any crms from redfly expr mapped file divide by total num of predictions excl ones used for tset")
			percentageOfHitsRecoveredExpressionPattern=no_of_overlapsSubset/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded)
			#print("Percent Expression Pattern Recovered =: "+str(percentageOfHitsRecoveredExpressionPattern))	
		
		

			#Expected Exprssion Pattern Recovery
			#Generating random data from genome and finding Significance of SCRMSHAW no of hits
			meanExp,medExp,sdExp,minVExp,maxVExp,statsExp,zExp,pExp=shuffle(sortedMergedModifiedSubsetCrmsBed,sortedMergedScrmsBed,ex,g,no_of_overlapsX,int(numberOfShuffles))		
			#Percentage of Expected expression pattern recovered
			#print("Expected Expression Pattern Recovered =: "+str(meanExp))
			percentageOfExpectedHitsRecoveredExpressionPattern=meanExp/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded)
			#calculating in difference among percentages of scrms and expected recovered expression pattern						
			differenceInPercentagesExpressionPatternRecovered=percentageOfHitsRecoveredExpressionPattern-percentageOfExpectedHitsRecoveredExpressionPattern
			#print("Percent Expected Expression Pattern Recovered =: "+ str(meanExp/(numOfsortedMergedPredictedCrms-countCommonScrmAndExcluded)))	
		
			#Overlap with shared annotation (Part III) 
			# getting number of crms predicted by scrmshaw that matches with any of the expression patterngroup
			hitsMatchedExpression,numberOfAnyPatternRecovered=parse_intersected(outfile,pathOfSubsetIntersectResult,whichSetToLookFor,tab1)
			
			#Expression Pattern Precision: TP/TP+FP..... True Positives= SCRMs with correct expression pattern...False Positive=SCRMs with incorrect expression pattern (SCRMs with any pattern - TP)
			percentageExpressionPatternPrecision=hitsMatchedExpression/numberOfAnyPatternRecovered
			#print("Percent Precision =: "+str(percentageExpressionPatternPrecision))	
			
			#Expression Pattern Recall/Sensitivity: TP/TP+FN ..... True Positives= SCRMs with correct expression pattern..False Negative= Total Known CRMs with Expression Pattern in Training set - TP
			percentageExpressionPatternRecall= hitsMatchedExpression/(hitsMatchedExpression+(numberOfKnownCrmsCausesExpressionInTset-hitsMatchedExpression))
			
			#print("testing Random Expectations Spec")
			percentageExpressionPatternPrecisionExp= []
			#trying to see the expected distribtuion across all the sets..
			for i in range(numberOfShuffles):
				shuffled2=sortedMergedScrmsBed.shuffle(excl=ex, noOverlapping=True, g=g).sort().saveas('shuffled')	
				shuffledpath=os.path.abspath('shuffled')
				shuffledBed=bed_conversion(shuffledpath)		
				no_of_overlapsR,pathOfIntersectResultR=intersect_scrms_and_modifiedCrms("expressionMapped",whichSetToLookFor,tab1,sortedMergedModifiedSubsetCrmsBed,shuffledBed)		
				
				hitsMatchedExpressionExp,numberOfAnyPatternRecoveredExp=parse_intersected(outfile_random,pathOfIntersectResultR,whichSetToLookFor,tab1)
				percentageExpressionPatternPrecisionExp.append(hitsMatchedExpressionExp/numberOfAnyPatternRecoveredExp)

			#see if you need these upto 3 #
			#calculating summary stats for the expected distribution result
			#get summary stats 
			maxRspec=max(percentageExpressionPatternPrecisionExp)
			minRspec=min(percentageExpressionPatternPrecisionExp)
			meanRspec=statistics.mean(percentageExpressionPatternPrecisionExp)
			medianRspec=statistics.median(percentageExpressionPatternPrecisionExp)
			sdRspec=statistics.stdev(percentageExpressionPatternPrecisionExp)
			#print("Median Percent Precision Expected =: "+str(medianRspec))	
			summary_stats="Mean= "+ str(meanRspec)+"\nMedian= "+str(medianRspec)+"\nStandard Deviation= "+str(sdRspec)+"\nMinimum value= "+str(minRspec)+"\nMaximum Value= "+str(maxRspec)

			#from a_and_b_count and mean sd of above array
			#print(a_and_b_count)
			try:	
				zRspec=float((percentageExpressionPatternPrecision-meanRspec)/sdRspec)
				#print(zRspec)
				#pRspec=stats.norm.sf(float(abs(zRspec)))
				pRspec= 1 - scipy.special.ndtr(zRspec)
				#print(pRspec)
			except ZeroDivisionError:
				np.seterr(divide='ignore', invalid='ignore')
				pRspec=0
			
			
			differenceInPercentagesExpressionPatternRecovered=percentageExpressionPatternPrecision-meanRspec

			#calculating summary stats for the expected distribution result
			calculating_sum_stats_random_file(outfile)
			#calculating fold enrichment 
			calculateFoldEnrichment(outfile,whichSetToLookFor,tab1)
			#temporary transpose to display the results of every parsed output on the terminal, for instance if user has selected to evaluate continuously, this eill show up the summary results fo every N predictions 
			transposeDistributionFile_temp(outfile,whichSetToLookFor)			
			printList=[]
			
			#writing to files depending upon _condition applied depending upon if Expression Pattern Recovery step is needed or not			
			if patternRecovery=='True' or  patternRecovery == "T" or patternRecovery=="TRUE":
				printList.extend([whichSetToLookFor,tab1,str(excludedCrmsCount),str(size_of_sortedMergedSubsetCrms),str(sizeOfSortedMergedModifiedSubsetCrms),str(size_of_sortedMergedCrms),str(sizeOfSortedMergedModifiedCrms),str(numOfsortedMergedPredictedCrms),str(countCommonScrmAndExcluded),str(percentageOfSensitivity),str(mean2V),str(percentageOfExpectedSensititvity),str(differenceInPercentagesTsetSensitivity),str(p_value2),str(no_of_overlaps),str(percentageOfOverlaps),str(meanV),str(percentageOfExpectedOverlaps),str(differenceInPercentagesRedflyRecovered),str(p_value),str(percentageRecallREDfly),str(no_of_overlapsSubset),str(percentageOfHitsRecoveredExpressionPattern),str(meanExp),str(percentageOfExpectedHitsRecoveredExpressionPattern),str(hitsMatchedExpression),str(numberOfAnyPatternRecovered),str(percentageExpressionPatternPrecision),str(meanRspec),str(pRspec),str(differenceInPercentagesExpressionPatternRecovered),str(percentageExpressionPatternRecall)])
			else:
				printList.extend([whichSetToLookFor,tab1,str(excludedCrmsCount),str(size_of_sortedMergedSubsetCrms),str(sizeOfSortedMergedModifiedSubsetCrms),str(size_of_sortedMergedCrms),str(sizeOfSortedMergedModifiedCrms),str(numOfsortedMergedPredictedCrms),str(countCommonScrmAndExcluded),str(percentageOfSensitivity),str(mean2V),str(percentageOfExpectedSensititvity),str(differenceInPercentagesTsetSensitivity),str(p_value2),str(no_of_overlaps),str(percentageOfOverlaps),str(meanV),str(percentageOfExpectedOverlaps),str(differenceInPercentagesRedflyRecovered),str(p_value),str(percentageRecallREDfly)])
			#	if tsetBedOrNot=="False" or tsetBedOrNot=="F" or tsetBedOrNot=="FALSE":
				#table_file_nonScrmshaw(outfile,x,tab1,str(excludedCrmsCount),str(size_of_sortedMergedSubsetCrms),str(sizeOfSortedMergedModifiedSubsetCrms),str(size_of_sortedMergedCrms),str(sizeOfSortedMergedModifiedCrms),str(numOfsortedMergedPredictedCrms),str(countCommonScrmAndExcluded),str(percentageOfSensitivity),str(mean2V),str(percentageOfExpectedSensititvity),str(differenceInPercentagesTsetSensitivity),str(p_value2),str(no_of_overlaps),str(percentageOfOverlaps),str(meanV),str(percentageOfExpectedOverlaps),str(differenceInPercentagesRedflyRecovered),str(p_value),str(no_of_overlapsSubset),str(percentageOfHitsRecoveredExpressionPattern),str(numberOfAnyPatternRecovered),str(meanExp),str(percentageOfExpectedHitsRecoveredExpressionPattern),str(differenceInPercentagesExpressionPatternRecovered),str(pExp))
			table_file(outfile,printList)
			
			
			print("Done",t,"set of length",numberparse)
			t=t+1
			if continuous=="True" or continuous=="TRUE" or continuous=="T":	
				print("adding 250 to numparse")
			numberparse=numberparse+250
		
		t=t-1
		#final transposing of the file
		transposeDistributionFile(outfile)
		#for calculation of max fold enrichment and max percent Precision
		max_FE(outfile,x)
main()	
#./generic_evaluationPipeline_scrmshawIncl.py -nameOfmet SCRMshaw -so scrmshawOutput_orig_2sets_5000hits.bed -fullredfly without_chrallredfly_2.5kb.July2017.txt -tset False -listTset trainingset_assignments_2010.txt -pattern TRUE -subsetcrmsExpBed without_chrredfly-analysis-set.assignments.2015.REDFly_format.txt -finalcrmsExp redfly_analysis_set.assignments.2015.txt -e exons.bed -drosog genome_chr_lengths_r6_copy.txt -o outfile_test_cont_5k -cont True -goodHits False -p 5000
#./generic_evaluationPipeline_scrmshawIncl.py -nameOfmet imogene -predInBed scrms.merged.bed -fullredfly allredfly_2.5kb.July2017.txt -tset False -listTset trainingset_assignments_2010.txt -pattern TRUE -subsetcrmsExpBed redfly-analysis-set.assignments.2015.REDFly_format.txt -finalcrmsExp redfly_analysis_set.assignments.2015.txt -e exons.bed -drosog genome_chr_lengths_r6_copy.txt -o scrmshawModifiedPipelineOutput_notScrmshawInput_tsetList_adultMesoImm_100shuffle_500preds_writeToFile -cont False -goodHits False -p 500 -s 10 -pattern True
#./evaluationPipeline_withFE.py -nameOfmet pnas -predInBed allPredictions_withoutChr.txt -fullredfly allredfly_2.5kb.July2017.txt -bedTset pnas_Tset_withoutChr.txt -pattern True -subsetcrmsExpBed redfly-analysis-set.assignments.2015.REDFly_format.txt -finalcrmsExp redfly_analysis_set.assignments.2015.txt -e exons.bed -drosog genome_chr_lengths_r6_copy.txt -o testingPnas -cont False -goodHits False -p 7000 -s 10 -tset True -which mapping1.blastoderm
