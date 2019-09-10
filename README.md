

################################################################################																			
   								  			
		`````````Enhancer Predictions Evaluation Pipeline`````````````
				   (pCRM_eval)
				   Halfon Lab			
  				Date: Sept 2018	 										
################################################################################

	1. INPUT
	2.  USAGE
	3. PARAMETERS
	4. OUTPUT


#1. INPUT

Following files are required to run the script. Most of these required files can be downloaded from our github site. All you need to provide is the set of predictions in bed format and the training set used (in bed or list format) in case of Supervised Method.\
i) Set of predictions in bed format. (SCRMshaw-HD or any other method's output)\
ii) genome of the specie\
iii) exons information\
iv) Training set assignment file (List or bed format)\
v) All REDfly CRMs:\
This file includes everything in REDfly annotated as a CRM, including evidence solely from cell culture experiments, that is no longer than 2.5kb in length. Its about 16,000 entities.\
vi) Expression Mapped CRMs: These are the expression-mapped CRMs, use to assess training set specificity\
vii) Expression Mapped CRMs in BED format: This will be the bed formatted version of above file. It will have the coordinates with the names of crms.


#2.  USAGE

Following basic modules are required to run this script.Please make sure these modules have already been properly installed and are recognizable.
pybedtools, statistics, scipy, numpy, pandas, csv and itertools

For evaluating SCRMshaw output, following script could be used (just change the file names according to your data)

>python generic_evaluationPipeline_NO_FE.py  -nameOfmet SCRMshaw  -fullredfly without_chrallredfly_2.5kb.July2017.txt  -subsetcrmsExpBed without_chrredfly-analysis-set.assignments.2015.REDFly_format.txt  -finalcrmsExp redfly_analysis_set.assignments.2015.txt  -tset False  -listTset trainingset_assignments_2010.txt  -pattern TRUE  -e exons.bed  -drosog genome_chr_lengths_r6_copy.txt  -so scrmshawOutput_peaksCalledover5kcrms_allSets_IMM.bed  -o output -cont False  -goodHits True  -p 35000  -s 10  -pattern True 

#3. PARAMETERS

	-nameOfmet	<str>	SCRMshaw/ name of any other method 
	-fullredfly	<str>	This file includes everything in REDfly annotated as a CRM, including evidence solely from cell culture experiments, that is no longer than 2.5kb in length. It'92s about 16,000 entities.  
	-subsetcrmsExpBed <str>	These are the expression-mapped CRMs, use to assess training set specificity, in the form of list
	-finalcrmsExp	<str>	This will be the bed formatted version of expression-mapped CRMs file. It will have the coordinates with the names of crms. 	
	-tset	 <str>	This is the binary variable, which takes the value of True or False depending upon if you have the training set information or not (Supervised vs unsupervised)
	-listTset	<str>	If the training set information is in the form of list, provide its name here.
	-bedTset	<str> If the training set information is in the form of bed, provide the name here.
	-pattern	<str>	This is binary variable which take the value of True or False depending on if user wants to test the specificity of the training set. 
	-e	<str>	Exon information file of Drosophila in bed format.
	-drosog	<str>	 Genome length file of Drosophila.

	-so	<str>	Set of Predictions, user wants to evaluate.
	-o	<str>	Name of output file.
	-cont	<str>	This is binary variable which takes the value of True or False, depending on user'92s choice if they want to evaluate the predictions in a semi continuous fashion (evaluating at every addition of 250 CRMs) 
	-goodHits	<str>	This is binary variable which takes the value of True or False, if user wants to evaluate Top Good Predictions only.
	-p 	<str>	Total number of predictions to start evaluation  from the predictions file. Default value is 35000
	-s	<str>	Number of shuffles you want to for calculating the p-value of your data. Default value is 20

#4. OUTPUT

The output will have following files.

1) output.bed: The output of evaluation in Tab Delimited format.
2) output_ExpressionPatternRecovereyDistributionExcel.bed: File containing the Expression pattern recovery information for the predictions, and the fold enrichment values for specificity.
3) output_random_ExpressionPatternRecovereyDistributionExcel: File containing the Expression pattern recovery information for the random expectations.
4) temp: Directory containing temporary files created during execution.

