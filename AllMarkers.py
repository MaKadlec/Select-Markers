#!/usr/bin/env python
# -*- coding: utf-8 -*-


##############################
##############################
##    Kadlec Malvina        ##
##############################

#-----------------------------------------------------------------------------#

## Select single copy sequences in several transcriptome in fasta format by 
# comparing sequences with Whole Genome Sequence (WGS) data or database of 
# single copy genes in flowering plants (De Smet et al, 2013, 
# doi:10.1073/pnas.1300127110)
# Input : several transcriptomes in a directory. Optional: WGS data.
# Output : selected sequences for each transcriptome in fasta format + homologues
# in WGS or database data
# Can be used with python 2 and python 3

#----------------------------------------------------------------------------#
######## Module #########

import argparse #For arguments
import os 
import glob #To obtain path(s)
from multiprocessing import Pool #Parallel
import shutil #To move files
from Bio.Blast.Applications import NcbiblastnCommandline #blastn
from Bio.Blast.Applications import NcbiblastxCommandline #blastx
from Bio import SearchIO #Parse blast output
from itertools import chain
import numpy
from Bio import SeqIO
import operator #Sort

#-----------------------------------------------------------------------------#
######## Arguments/Help #######

parser = argparse.ArgumentParser(description="Select single copy sequences from data in fasta format")
parser.add_argument("input", help="Transcriptome(s) in fasta format. Must be in a directory")
parser.add_argument("database", help="File in fasta format used as a database to selected single copy")
parser.add_argument("dbtype", help="Type of database used. \'nucl\' if database's sequences are ADN or ARN. \'prot\' if database's sequences are protein.")
parser.add_argument("-NbT", type=int, nargs="?", default=2, 
	help="Minimum number of transcriptome where the selected sequence is found. Default = 2")
parser.add_argument("-np", nargs="?", default=1, help="Number of CPU used", type=int)
parser.add_argument("-l", "--length", type=int, nargs="?", 
	help="Sequences under this length will be removed")
parser.add_argument("-sim", type=int, default=75, nargs="?", help="Sequences under this percentage identity will be removed. Default = 75")
parser.parse_args()
ArgsDico = vars(parser.parse_args()) #Make a dictionnary to process arguments

#----------------------------------------------------------------------------#
######## Functions #########

#Database for Blast
def MakeDB(FastaFile):
	''' Commandline to make database for Blastn. The function take a
	file in fasta format as argument and return the database'''
	BlastDBCDline = "makeblastdb -in "+FastaFile+" -dbtype \"nucl\" -out "+FastaFile+".DB"
	os.system(BlastDBCDline)
	
#Blastn commandline
def MakeBlastn (Query, DBs, OutName, EVAL):
	'''Commandline for blastn. 
	Arguments: query sequences (fasta file with sequences), database, 
		name output, evalue
	Return: multi results output in xml format'''
	BlastNCline = NcbiblastnCommandline(query = Query, db=DBs, evalue=EVAL,
	 out=OutName, outfmt = 5) #Obtain blastn result in xml format
	stdout, stderr = BlastNCline()
	
#Blastx commandline
def MakeBlastx (Query, DBs, OutName, EVAL):
	'''Commandline for blastx. 
	Arguments: query sequence (fasta file with sequences, database, 
		name output, evalue)
	Return: multi results output in xml format'''
	BlastXCline = NcbiblastxCommandline(query = Query, db=DBs, evalue=EVAL,
	 out=OutName, outfmt = 5) #Obtain blastn result in xml format
	stdout, stderr = BlastXCline()
	
#Do blastn search
def DoBlastn(QueryDB):
	'''Parse list of query, database and evalue then do blastn search.
	Argument: List (ex: [query, database, evalue]
	Return: multi results output in xml format'''
	#Query's path
	QueryPath = QueryDB[0]
	#Database's path
	DB = QueryDB[1].replace(".nhr","")
	#To name output blast file, need the database species' name.
	DBStr2List = DB.split("/")
	NameOutBlast = QueryPath+"vsDB"+DBStr2List[-1].replace(".DB",".xml")
	#Blastn search commandline
	MakeBlastn(QueryPath, DB, NameOutBlast, QueryDB[2])
	
#Do blastn search
def DoBlastx(QueryDB):
	'''Parse list of query, database and evalue then do blastx search.
	Argument: List (ex: [query, database, evalue]
	Return: multi results output in xml format'''
	#Query's path
	QueryPath = QueryDB[0]
	#Database's path
	DB = QueryDB[1].replace(".phr","")
	#To name output blast file, need the database specie's name.
	DBStr2List = DB.split("/")
	NameOutBlast = QueryPath+"vsDB"+DBStr2List[-1].replace(".DB",".xml")
	#Blastn search commandline
	MakeBlastx(QueryPath, DB, NameOutBlast, QueryDB[2])

#Parse Blast output
def ParseBlast(BlastFile):
	''' Parsing blast output.
	Argument: File from blastn search (Transcriptome)
	Return: List. The first two argument are the name of the species used
		as query and database. Other are dictionnaries (keys: 
		Query species' name, Query sequence's accession number,
		Hit species' name, Hit sequence's accession number)'''
	ResultBlastList = [] #For list of dictionnaries
	for qresult in SearchIO.parse(BlastFile,"blast-xml"): 
		if len(qresult) != 0: #Keep all sequences with hits
			ResultBlast = {}
			BlastFileList = BlastFile.replace(".xml","").split("vsDB")
			NameQuery = BlastFileList[0].split("/")
			ResultBlast["Query species"] = NameQuery[-1]
			ResultBlast["Hit species"] = BlastFileList[-1]
			ResultBlast["Query id"] = qresult.id
			HitBlast = []
			for hit in qresult.hits:
				HitBlast.append(hit.id)
			ResultBlast["Hit id"] = HitBlast
			ResultBlastList.append(ResultBlast)
		else:
			continue
	return(ResultBlastList)
	
#In which transcriptome does this sequence have a homologue ? 
def FoundHomologueInTrans(ListIdBlastout):
	'''Count the number of transcriptome with homologuous sequences.
	Argument : List ([Id, Species, Result list of parsing blast output])
	Return : Dictionnary (Key: Query id, Query species, Species)'''
	DicoNew = {}
	IdQuery = ListIdBlastout[0]
	SpeciesQuery = ListIdBlastout[1]
	AllResultBlast = ListIdBlastout[2]
	DicoNew["Query id"]= IdQuery
	DicoNew["Query species"] = SpeciesQuery
	Species = []
	for Dico in AllResultBlast:
		if IdQuery == Dico["Query id"]:
			Species.append(Dico["Hit species"])
			Species.append(Dico["Query species"])
		if Dico["Hit id"].count(IdQuery) != 0:
			Species.append(Dico["Query species"])
			Species.append(Dico["Hit species"])
	#Do not repeat species => unique
	SpeciesArray = numpy.array(Species)
	SpeciesArrayUnique = numpy.unique(SpeciesArray)
	SpeciesOne = SpeciesArrayUnique.tolist()
	DicoNew["Species"] = SpeciesOne
	return(DicoNew)
	
# Parse Fasta File
def ExtractFromFa(FastaFile):
	'''Extract sequences informations from a fasta file
	Argument : file in fasta format
	Return : List of dictionnary with info about sequences'''
	DicoList = []
	#Fill dico with information coming from fasta file 
	for SeqRecord in SeqIO.parse(FastaFile, "fasta"):
		DicInput = {}
		DicInput["id"] = SeqRecord.id
		DicInput["Name"] = SeqRecord.description
		DicInput["Sequence"] = SeqRecord.seq
		DicoList.append(DicInput)
	return(DicoList)

#Given its id, find sequence and related information (length, complete name...)
def FindInfoInFa(ListId, ListDico):
	'''Given sequence's id, retrieve information coming from a parse
	fasta file (see function: ExtractFromFa)
	Arguments: List of id, Dictionnary coming from parsing
	Return: New list of dictionnaries'''
	ListSeq = []
	for Id in ListId:
		for DicoF in ListDico:
			if Id == DicoF["id"]:
				SeqInfo = {}
				SeqInfo["Name"] = DicoF["Name"]
				SeqInfo["Id"] = Id
				SeqInfo["Sequence"] = DicoF["Sequence"]
				ListSeq.append(SeqInfo)
	return(ListSeq)

	
#Write in file
def WriteTmpSelectedSequences(DicoTrans):
	'''Write selected sequences in a file
	Argument: Dictionnary of selected sequences from one species
	Return: file in fasta format'''
	Trans = DicoTrans["Species"]
	#Use function ExtractFromFa to obtain sequences of all sequences
	DicoFaTrans = ExtractFromFa(ArgsDico['input']+"/"+Trans)
	ListId = DicoTrans["Ids"]
	#Do not repeat id => unique
	ListIdArray = numpy.array(ListId)
	ListIdArrayUnique = numpy.unique(ListIdArray)
	ListIdOne = ListIdArrayUnique.tolist()
	#Use function FindInfoInFa to obtain sequence of selected sequences
	InfoSeq = FindInfoInFa(ListIdOne, DicoFaTrans)
	DicoTrans["Info Seq"] = InfoSeq
	#Remove list of id
	DicoTrans.pop("Ids", 0)
	NameSpe = Trans.replace(".fasta","")
	Output = open(ArgsDico['input']+"/TranscriptomeHomologues/"+NameSpe+"_SelectedSeq.fasta", "w")
	for SeqRecord in DicoTrans["Info Seq"]:
		Output.write(">"+str(SeqRecord["Name"])+"\n"+str(SeqRecord["Sequence"])+"\n")
	Output.close()

# Parse blast result against WGS or database => find single copy
def ParseBlastWGS(BlastFilePath):
	'''Parse blast output of blast selected sequences against whole genome 
	or know database
	Argument: blast output
	Return: List of dictionnaries'''
	ListIdDico = [] 
	Dico4Blast = [] #for dico
	# Parse blast output and make dictionnary with informations coming from parsing
	BlastFileList = BlastFilePath.replace(".xml","").split("/")
	NameTransQuery = BlastFileList[-1]
	ListIdDico.append(NameTransQuery)
	for qresult in SearchIO.parse(BlastFilePath,"blast-xml"):
		if len(qresult) == 1: #Keep single match
			Eval = qresult.hsps[0].evalue
			IdtResidu = qresult.hsps[0].ident_num
			LenAln = qresult.hsps[0][0].aln_span
			IdtPercentage = round(((float(IdtResidu)/LenAln)*100), 0)
			if IdtPercentage >= ArgsDico["sim"]:
				ListIdDico.append(qresult.id)
				#Dico to obtain percentage identity
				DicoB = {}
				DicoB["Query id"] = qresult.id
				DicoB["% identity"] = IdtPercentage
				Dico4Blast.append(DicoB)
	ListIdDico.append(Dico4Blast)
	return(ListIdDico)

#----------------------------------------------------------------------------#
######## Main #########

print(ArgsDico)
if __name__ == '__main__':
	p = Pool(ArgsDico['np']) #Parallel
	ListFile =  glob.glob(ArgsDico['input']+"/*" )
	print(ListFile)
	#Make Database
	ResultDB = p.map(MakeDB, ListFile)
	os.mkdir(ArgsDico['input']+"/Databases")
	#Moves files
	for File in glob.glob(ArgsDico['input']+"/*DB*"):
		shutil.move(File, ArgsDico['input']+"/Databases")
	DBList = glob.glob(ArgsDico['input']+"/Databases/*DB.nhr" )
	print("Making Blast")
	#Filter by length
	if ArgsDico['length'] is None: #No filter
		#Combine Trans and DB in a list of lists => parallel
		AllCombList = []
		for Trans in ListFile:
			TransPathList = Trans.split("/")
			NameTrans = TransPathList[-1]
			for DB in DBList:
				if DB.find(NameTrans) == -1:
					CombList = []
					CombList.append(Trans)
					CombList.append(DB)
					CombList.append(0.01) #Evalue for blastn search
					AllCombList.append(CombList)
	else: #Filter by length has been set in commandline
		os.mkdir(ArgsDico['input']+"/SeqAboveFilterLength")
		#Filter by length
		for Trans in ListFile:
			#Name's transcriptome
			TransPathList = Trans.split("/")
			NameTrans = TransPathList[-1]
			SelectByLenFile = open(ArgsDico['input']+"/SeqAboveFilterLength/"+NameTrans,"w")
			for SeqRecord in SeqIO.parse(Trans, "fasta"):
				if len(SeqRecord.seq) >= ArgsDico['length']:
					SelectByLenFile.write(">"+str(SeqRecord.name)+"\n"+str(SeqRecord.seq)+"\n")
		#Combine Trans and DB in a list of lists => parallel
		AllCombList = []
		NewTransFiles = glob.glob(ArgsDico['input']+"/SeqAboveFilterLength/*")
		for Trans in NewTransFiles:
			TransPathList = Trans.split("/")
			NameTrans = TransPathList[-1]
			for DB in DBList:
				if DB.find(NameTrans) == -1:
					CombList = []
					CombList.append(Trans)
					CombList.append(DB)
					CombList.append(0.01) #Evalue for blastn search
					AllCombList.append(CombList)
	#Do Blast 
	ResultBlast = p.map(DoBlastn, AllCombList)
	#Move Files
	os.mkdir(ArgsDico['input']+"/Blast")
	if ArgsDico['length'] is None: #No filter
		for File in glob.glob(ArgsDico['input']+"/*.xml"): 
			shutil.move(File, ArgsDico['input']+"/Blast")
	else: #Files in SeqAboveFilterLength folder
		for File in glob.glob(ArgsDico['input']+"/SeqAboveFilterLength/*.xml"): 
			shutil.move(File, ArgsDico['input']+"/Blast")
	#Parse Blast 
	print("Parse Blast files")
	ListBlastFile = glob.glob(ArgsDico['input']+"/Blast/*")
	ResultParse = p.map(ParseBlast, ListBlastFile)
	#Combine all result for all transcriptomes in one big list
	MergedResultParseBlast = list(chain(*ResultParse))
	print ("All results length: ")
	print(len(MergedResultParseBlast))
	# have sequences in the right format for next function and process in parallel
	QueryIdvsMergeDico = []
	for Dico in MergedResultParseBlast:
		NewList = []
		NewList.append(Dico['Query id'])
		NewList.append(Dico["Query species"])
		NewList.append(MergedResultParseBlast)
		QueryIdvsMergeDico.append(NewList)
	print("Merge")
	print(len(QueryIdvsMergeDico))
	# Find homologue(s) in transcriptomes
	AllIdSpe = p.map(FoundHomologueInTrans, QueryIdvsMergeDico)
	print("Homologues found: ")
	print(len(AllIdSpe))
	IdHomologues = []
	for DicoSp in AllIdSpe:
		#Keep homoloques which appear in X transcriptome(s)
		if len(DicoSp["Species"]) >= ArgsDico['NbT']:
			IdHomologuesDico = {}
			IdHomologuesDico["Query id"] = DicoSp["Query id"]
			IdHomologuesDico["Query species"] = DicoSp["Query species"]
			IdHomologues.append(IdHomologuesDico)
		else:
			continue
	print("Homologues in 2 or more transcriptomes: ") 
	print(len(IdHomologues))
	#Sort id by transcriptome
	SortIdHomologues = []
	for PathTrans in ListFile:
		TransHomologues = {}
		Query2List = PathTrans.split("/")
		Trans = Query2List [-1]
		TransHomologues["Species"] = Trans
		Id = []
		for IdDico in IdHomologues:
			if IdDico["Query species"] == Trans:
				Id.append(IdDico["Query id"])
			else:
				continue
		TransHomologues["Ids"] = Id
		SortIdHomologues.append(TransHomologues)
	print("Sort")
	print(len(SortIdHomologues))
	#Add sequences and Write selected sequences in file
	os.mkdir(ArgsDico['input']+"/TranscriptomeHomologues")
	p.map(WriteTmpSelectedSequences, SortIdHomologues)
	#Blast to select single copy sequences
	#Make Database for single copy sequences selection
	if ArgsDico['dbtype'] == "nucl":
		MakeDB(ArgsDico["database"])
	elif ArgsDico['dbtype'] == "prot":
		BlastDBCDline = "makeblastdb -in "+ArgsDico["database"]+" -dbtype \"prot\" -out "+ArgsDico["database"]+".DB"
		os.system(BlastDBCDline)
	#Need name of WGS file for blast
	if ArgsDico['dbtype'] == "nucl":
		NameSpeWGS = ArgsDico["database"]+".DB.nhr"
	elif ArgsDico['dbtype'] == "prot":
		NameSpeWGS = ArgsDico["database"]+".DB.phr"	
	# Files to use in blast
	SelectTranscriptomeFiles = glob.glob(ArgsDico['input']+"/TranscriptomeHomologues/*")
	#Combine Transcriptomes and DB in a list of list => parallel
	TransvsWGSList = []
	for Trans in SelectTranscriptomeFiles:
		CombList = []
		CombList.append(Trans)
		CombList.append(NameSpeWGS)
		CombList.append(0.001) #Evalue for blastn search
		TransvsWGSList.append(CombList)
	print(TransvsWGSList)
	#Do Blast
	if ArgsDico['dbtype'] == "nucl":
		WGSResultBlast = p.map(DoBlastn, TransvsWGSList)
	elif ArgsDico['dbtype'] == "prot":
		WGSResultBlast = p.map(DoBlastx, TransvsWGSList)	
	#Move database file
	DBFiles = glob.glob(ArgsDico["database"]+"*DB.*" )
	for DBFile in DBFiles:
		shutil.move(DBFile, ArgsDico['input']+"/Databases")
	#Move blast files
	os.mkdir(ArgsDico['input']+"/BlastWGS")
	for File in glob.glob(ArgsDico['input']+"/TranscriptomeHomologues/*.xml"):
		shutil.move(File, ArgsDico['input']+"/BlastWGS")
	#Parse blast files 
	BlastFileWGS = glob.glob(ArgsDico['input']+"/BlastWGS/*")
	ResultsBlastWGS = p.map(ParseBlastWGS, BlastFileWGS)
	#Need a new directory for final selected sequence dataset
	os.mkdir(ArgsDico['input']+"/Final_AllSelectedSequences")
	for BlastWGS in ResultsBlastWGS:
		TransNameSpeList = BlastWGS[0].split("_")
		TransNameSpe = TransNameSpeList[0]
		del BlastWGS[0]
		#To obtain % identity:
		IdentityDico = BlastWGS[-1]
		del BlastWGS[-1]
		DicoList = [] 
		print("Retrieve sequence in fasta file")
		DicoTransFa = ExtractFromFa(ArgsDico['input']+"/TranscriptomeHomologues/"+TransNameSpe+"_SelectedSeq.fasta")
		#Find selected sequence with id
		print("find selected sequence")
		SelectedSeqs = FindInfoInFa(BlastWGS, DicoTransFa)
		for SelectedSeq in SelectedSeqs:
			for Dico in IdentityDico:
				if SelectedSeq["Id"] == Dico["Query id"]:
					SelectedSeq["% identity"] = Dico["% identity"]
					SelectedSeq["Length Sequence"] = len(str(SelectedSeq["Sequence"]))
					break
		SortSelectedSeqs = sorted(SelectedSeqs, key=operator.itemgetter("Length Sequence", "% identity"))
		SortSelectedSeqs.reverse()
		print(len(SortSelectedSeqs))
		#Write selected sequences to file
		OutputSelected = open(ArgsDico['input']+"/Final_AllSelectedSequences/"+TransNameSpe, "w")
		for SeqDico in SortSelectedSeqs:
			OutputSelected.write(">"+str(SeqDico["Name"])+"\n"+str(SeqDico["Sequence"])+"\n")
		OutputSelected.close()

