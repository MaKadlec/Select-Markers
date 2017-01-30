#!/usr/bin/env python
# -*- coding: utf-8 -*-


##############################
##############################
## Kadlec Malvina
##############################

#-----------------------------------------------------------------------------#

#### BestMarkers.py version 1.0 #####

# Part II of pipeline used in paper "Targeted NGS for species level phylogenomics: 
#	“made to measure” or “one size fits all”? (Malvina Kadlec, Dirk U. Bellstedt, 
#	Nicholas C. Le Maitre, and Michael D. Pirie)
# BestMarkers select the best markers available in a pool a selected markers (file
#	in fasta format) according to length and/or similarity but also according
#	to the number, length and coverage of bapture library (or baits or probes)
#	to be used in Tageted capture.
# To select the best markers acording similarity, the script needs a multiresult 
#	blast file in xml format.
# Don't forget to make the script executable before using !


#----------------------------------------------------------------------------#

######## Module import #########

import argparse #For arguments
from Bio import SeqIO, SearchIO
from operator import itemgetter #To sort list of dico
import numpy

#-----------------------------------------------------------------------------#

######## Arguments/Help #######

parser = argparse.ArgumentParser(description="Select best markers from data in fasta format")
parser.add_argument("input", help="Sequences in fasta format")
parser.add_argument("-l", "--length", type=int, nargs="?", 
	help="Sequences with predicted length or 'real length' under setting will be removed")
parser.add_argument("-Blast_file", help="Blast output (xml format). For coding regions identification")
parser.add_argument("-intron?", choices=["Y", "N"], help="Do you want to take the presence of introns into account?")
parser.add_argument("-Length_intron", type=int, help="Set an arbitrary length for introns")
parser.add_argument("-only", choices=["length","similarity"], help="Sort sequences only by length OR similarity")
parser.add_argument("-priority", choices=["length","similarity"], help="Which parameters have the priority in the choice of marker?")
parser.add_argument("-Number_baits", type=int, help="Number of baits to reach", default=5770)
parser.add_argument("-Length_baits", type=int, help="Length of baits used", default=120)
parser.add_argument("-Coverage_baits", type=int, help="Length of baits used", default=4)
parser.parse_args()
ArgsDico = vars(parser.parse_args()) #Make a dictionnary to process arguments

#Add average length

#----------------------------------------------------------------------------#

######## Functions #########

# Parse blast result	
def ParseBlast(BlastFile):
	Dico4Blast = [] #for dico
	# Parse blast output and make dictionnary with informations coming from parsing
	for qresult in SearchIO.parse(BlastFile,"blast-xml"):
		if len(qresult) >= 1: #Keep single match
			DicoB = {}
			#Fill dico
			DicoB["Query id"] = qresult.id
			DicoB["Query name"] = qresult.description
			DicoB["Coding Region ?"] = len(qresult.hsps)
			if len(qresult.hsps) > 1: #Predict serveral coding regions
					PercentageIDCRs = [] #for percentage identity
					LenCRs = [] #22 v'la les flics
					CrRanges = [] #Where are the coding regions?
					LensIntron = [] #Intron length
					SeqCRs = [] #Codings regions' sequences
					for hsp in qresult.hsps:			
						IdtResiduH = hsp.ident_num
						for HspFragment in hsp:
							#Percentage identity
							LenAlnFrag = HspFragment.aln_span
							IdtFrag = round(((float(IdtResiduH)/LenAlnFrag)*100),0)
							LenCRs.append(LenAlnFrag)
							PercentageIDCRs.append(IdtFrag)
							#Coding regions' sequences
							Seq = HspFragment.aln[0].seq
							Seq.upper()
							SeqWithoutGap = str(Seq).replace("-","")
							SeqUp = SeqWithoutGap.upper()
							SeqCRs.append(SeqUp)
							# CR : starting and ending positions
							Range = []
							Range.append(HspFragment.hit_start)
							Range.append(HspFragment.hit_end)
							CrRanges.append(Range)
					#Introns length
					CrRanges.sort()
					for i in xrange(len(CrRanges)-1):
						if CrRanges[i][1] < CrRanges[i+1][0]:
							LenIntron = CrRanges[i+1][0] - CrRanges[i][1]
							LensIntron.append(LenIntron)
						else:
							continue
					#Fill dico
					DicoB["CR sequence"] = SeqCRs
					DicoB["Introns length"] = LensIntron
					DicoB["CR Length"] = LenCRs
					DicoB["Percentage identity"] = round(numpy.mean(PercentageIDCRs),0)
			else: #1 CR
					for hsp in qresult.hsps:			
						IdtResiduH = hsp.ident_num
						for HspFragment in hsp:
							#Identity
							LenAlnFrag = HspFragment.aln_span
							IdtFrag = round(((float(IdtResiduH)/LenAlnFrag)*100),0)
							#CR sequence
							Seq = HspFragment.aln[0].seq
							Seq.upper()
							SeqWithoutGap = str(Seq).replace("-","")
							SeqUp = SeqWithoutGap.upper()
					#Fill dico
					DicoB["Percentage identity"] = IdtFrag
					DicoB["CR Length"] = LenAlnFrag
					DicoB["Introns length"] = 0
					DicoB["CR sequence"] = SeqUp
			Dico4Blast.append(DicoB)
		else: #No match
			DicoB = {}
			DicoB["Query id"] = qresult.id
			DicoB["Query name"] = qresult.description
			DicoB["Coding Region ?"] = 1
			DicoB["Percentage identity"] = 0
			DicoB["CR Length"] = 0
			DicoB["Introns length"] = 0
			DicoB["CR sequence"] = "NA"
			Dico4Blast.append(DicoB)
	return Dico4Blast

#Fill a dictionnary with informations coming from file in fasta format
def DicoFasta(FastaFile):
	DicoList = [] 
	#Fill dico with information coming from fasta file (use for making database)
	for SeqRecord in SeqIO.parse(FastaFile, "fasta"):
		DicInput = {}
		DicInput["id"] = SeqRecord.id
		DicInput["Name"] = SeqRecord.description
		DicInput["Sequence"] = SeqRecord.seq
		DicInput["Length"] = len(str(SeqRecord.seq))
		DicoList.append(DicInput)
	return DicoList
	
#-----------------------------------------------------------------------------#

######## Main #########

for key, value in ArgsDico.iteritems() :
    print key, value
    
    
## All selected sequences : parse
DicoFastaFile = DicoFasta(ArgsDico['input'])

## Parse blast file if present
if ArgsDico["Blast_file"] is not None:
	DicoBlast = ParseBlast(ArgsDico["Blast_file"])
	#Add blast informations in Dictionnary of all selected sequences
	for DicoF in DicoFastaFile:
		for DicoBl in DicoBlast: 
			if DicoF["id"] == DicoBl["Query id"]:
				DicoF["Coding Region ?"] = DicoBl["Coding Region ?"]
				DicoF["CR sequence"] = DicoBl["CR sequence"]
				DicoF["Introns length"] = DicoBl["Introns length"]
				DicoF["CR Length"] = DicoBl["CR Length"]
				DicoF["Percentage identity"] = DicoBl["Percentage identity"]
				break
else:
	ArgsDico["intron?"] = "N"
	
#Taking introns into account ?
BestSeqs = []
#Baits coverage
LenCoverage = ArgsDico["Length_baits"]/ArgsDico["Coverage_baits"] 
LenBaits = ArgsDico["Length_baits"] #Bait's length
if ArgsDico["Blast_file"] is not None: #Take intron into account.
	#No arbitrary length is set for intron. Take mean length 
	if ArgsDico["Length_intron"] is None:
		#Introns : what is the average length ? (useful for next step)
		AllIntronsTotalLength = 0
		TotalNumberOfIntrons = 0
		CRNumber = 0
		CRsNumber = 0
		for DicoFF in DicoFastaFile: 
			IntronsLenList = DicoFF["Introns length"]
			if type(IntronsLenList) is list:
				CRsNumber = CRsNumber + 1
				for Len in IntronsLenList:
					AllIntronsTotalLength = AllIntronsTotalLength + Len
					TotalNumberOfIntrons = TotalNumberOfIntrons + len(IntronsLenList)
			else:
				CRNumber = CRNumber +1
				if DicoFF["CR Length"] == 0:
					DicoFF["CR Length"] = len(DicoFF["Sequence"])	
					DicoFF["CR sequence"] = DicoFF["Sequence"]
		#Average intron length		
		Mean = AllIntronsTotalLength/TotalNumberOfIntrons
		#Add average intron length in dico
		for DicoFF1 in DicoFastaFile:
			LenSeq = DicoFF1["Length"]
			if DicoFF1["Coding Region ?"] != 1: #Several CR
				#Add arbitrary length for introns
				Len2Add = (DicoFF1["Coding Region ?"]-1)*Mean
				TrueLenSeq = LenSeq + Len2Add
				DicoFF1["Predict Length"] = TrueLenSeq
			else:
				DicoFF1["Predict Length"] = LenSeq
		if ArgsDico["only"] is None: #Sort dico according to length and similarity.
			if ArgsDico["priority"] == "length": #First, sort by length, then, similarity
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter("Percentage identity"))
				DicoFastaSort = sorted(DicoFastaSort, key=itemgetter('Predict Length'), reverse=True)
			else: #sort by similarity, then length
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter('Predict Length'), reverse=True)
				DicoFastaSort = sorted(DicoFastaSort, key=itemgetter("Percentage identity"))
		else: #Sort by similarity OR length
			if ArgsDico["only"] == "length":
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter('Predict Length'), reverse=True)
			else: #Sort by similarity
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter("Percentage identity"))
		#Keep sequences until you reach X baits (X set by an option)
		TotalProbes = 0 #Total number of probes predicted
		for DicoSort in DicoFastaSort:
			if DicoSort["Predict Length"] >= ArgsDico["length"]:
				if DicoSort["Coding Region ?"] == 1: # 1 CR
					NumberProbe = 0
					Seq = DicoSort["Sequence"]					
					for i in xrange(0,len(Seq),LenCoverage):
					#Keep only part seq with a size of the baits
						if len(Seq[i:i+LenBaits]) == LenBaits:
							if TotalProbes < ArgsDico["Number_baits"]:
								NumberProbe = NumberProbe +1
								BestSeqs.append(DicoSort)
							else: 
								break
						else:
							continue
					TotalProbes = TotalProbes + NumberProbe	
				else: #several CRs
					CRsSeqs = DicoSort["CR sequence"]
					for CrSeq in CRsSeqs:
						Seq = str(CrSeq)
						NumberProbe = 0
						#Is sequence length above baits length ?
						if len(Seq) >= LenBaits: 
							for i in xrange(0,len(Seq),LenCoverage):
								#Keep only part seq with a size of the baits
								if len(Seq[i:i+LenBaits]) == LenBaits:
									if TotalProbes < ArgsDico["Number_baits"]:
										NumberProbe = NumberProbe +1
										BestSeqs.append(DicoSort)
									else: 
										break
								else:
									continue
							TotalProbes = TotalProbes + NumberProbe
						else:
							continue
			else:
				continue
	else: #Arbitrary length is set for intron	
		Mean = ArgsDico["Length_intron"]
		#Add average intron length in dico
		for DicoFF1 in DicoFastaFile:
			LenSeq = DicoFF1["Length"]
			if DicoFF1["Coding Region ?"] != 1: #Several CR
				#Add arbitrary length for introns
				Len2Add = (DicoFF1["Coding Region ?"]-1)*Mean
				TrueLenSeq = LenSeq + Len2Add
				DicoFF1["Predict Length"] = TrueLenSeq
			else:
				DicoFF1["Predict Length"] = LenSeq
		if ArgsDico["only"] is None: #Sort dico according to length and similarity.
			if ArgsDico["priority"] == "length": #First, sort by length, then, similarity
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter("Percentage identity"))
				DicoFastaSort = sorted(DicoFastaSort, key=itemgetter('Predict Length'), reverse=True)
			else: #sort by similarity, then length
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter('Predict Length'), reverse=True)
				DicoFastaSort = sorted(DicoFastaSort, key=itemgetter("Percentage identity"))
		else: #Sort by similarity OR length
			if ArgsDico["only"] == "length":
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter('Predict Length'), reverse=True)
			else: #Sort by similarity
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter("Percentage identity"))
		#Keep sequences until you reach X baits (X set by an option)
		TotalProbes = 0 #Total number of probes predicted
		for DicoSort in DicoFastaSort:
			if DicoSort["Predict Length"] >= ArgsDico["length"]:
				if DicoSort["Coding Region ?"] == 1: # 1 CR
					NumberProbe = 0
					Seq = DicoSort["Sequence"]
					for i in xrange(0,len(Seq),LenCoverage):
					#Keep only part seq with a size of the baits
						if len(Seq[i:i+LenBaits]) == LenBaits:
							if TotalProbes < ArgsDico["Number_baits"]:
								NumberProbe = NumberProbe +1
								BestSeqs.append(DicoSort)
							else: 
								break
						else:
							continue
					TotalProbes = TotalProbes + NumberProbe	
				else: #several CRs
					CRsSeqs = DicoSort["CR sequence"]
					for CrSeq in CRsSeqs:
						Seq = str(CrSeq)
						NumberProbe = 0
						#Is sequence length above baits length ?
						if len(Seq) >= LenBaits: 
							for i in xrange(0,len(Seq),LenCoverage):
								#Keep only part seq with a size of the baits
								if len(Seq[i:i+LenBaits]) == LenBaits:
									if TotalProbes < ArgsDico["Number_baits"]:
										NumberProbe = NumberProbe +1
										BestSeqs.append(DicoSort)
									else: 
										break
								else:
									continue
							TotalProbes = TotalProbes + NumberProbe
						else:
							continue
			else:
				continue
else:
	if ArgsDico["Blast_file"] is None:
		#Sort selected sequences by length
		DicoFastaSort = sorted(DicoFastaFile, key=itemgetter("Length"), reverse=True)
		TotalProbes = 0 #Total number of probes predicted
		for DicoSort in DicoFastaSort:
			if DicoSort["Length"] >= ArgsDico["length"]:
				NumberProbe = 0
				Seq = DicoSort["Sequence"]
				for i in xrange(0,len(Seq),LenCoverage):
					#Keep only part seq with a size of the baits
					if len(Seq[i:i+LenBaits]) == LenBaits:
						if TotalProbes < ArgsDico["Number_baits"]:
							NumberProbe = NumberProbe +1
							BestSeqs.append(DicoSort)
						else: 
							break
					else:
						continue
				TotalProbes = TotalProbes + NumberProbe
	else:
		#Introns : what is the average length (For recap files)
		for DicoFF in DicoFastaFile: 
			IntronsLenList = DicoFF["Introns length"]
			if type(IntronsLenList) is list:
				CRsNumber = CRsNumber + 1
				for Len in IntronsLenList:
					AllIntronsTotalLength = AllIntronsTotalLength + Len
					TotalNumberOfIntrons = TotalNumberOfIntrons + len(IntronsLenList)
			else:
				CRNumber = CRNumber +1
				if DicoFF["CR Length"] == 0:
					DicoFF["CR Length"] = len(DicoFF["Sequence"])	
					DicoFF["CR sequence"] = DicoFF["Sequence"]
		#Average intron length		
		Mean = AllIntronsTotalLength/TotalNumberOfIntrons
		#Add average intron length in dico
		for DicoFF1 in DicoFastaFile:
			LenSeq = DicoFF1["Length"]
			if DicoFF1["Coding Region ?"] != 1: #Several CR
				#Add arbitrary length for introns
				Len2Add = (DicoFF1["Coding Region ?"]-1)*Mean
				TrueLenSeq = LenSeq + Len2Add
				DicoFF1["Predict Length"] = TrueLenSeq
			else:
				DicoFF1["Predict Length"] = LenSeq
		if ArgsDico["only"] is None: # sort according to length and similarity
			if ArgsDico["priority"] == "length": #First, sort by length, then, similarity
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter("Percentage identity"))
				DicoFastaSort = sorted(DicoFastaSort, key=itemgetter('Length'), reverse=True)
			else: #sort by similarity, then length
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter('Length'), reverse=True)
				DicoFastaSort = sorted(DicoFastaSort, key=itemgetter("Percentage identity"))
		else: #sort only by length or similarity
			if ArgsDico["only"] == "length":
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter('Length'), reverse=True)
			else: #Sort by similarity
				DicoFastaSort = sorted(DicoFastaFile, key=itemgetter("Percentage identity"))
		#Keep sequences until you reach X baits (X set by an option)
		AllIntronsTotalLength = 0
		TotalNumberOfIntrons = 0
		CRNumber = 0
		CRsNumber = 0
		TotalProbes = 0 #Total number of probes predicted
		for DicoSort in DicoFastaSort:
			if DicoSort["Length"] >= ArgsDico["length"]:
				if DicoSort["Coding Region ?"] == 1: # 1 CR
					NumberProbe = 0
					Seq = DicoSort["Sequence"]
					for i in xrange(0,len(Seq),LenCoverage):
						#Keep only part seq with a size of the baits
						if len(Seq[i:i+LenBaits]) == LenBaits:
							if TotalProbes < ArgsDico["Number_baits"]:
								NumberProbe = NumberProbe +1
								BestSeqs.append(DicoSort)
							else: 
								break
						else:
							continue
					TotalProbes = TotalProbes + NumberProbe	
				else: #several CRs
					CRsSeqs = DicoSort["CR sequence"]
					for CrSeq in CRsSeqs:
						Seq = str(CrSeq)
						NumberProbe = 0
						#Is sequence length above baits length ?
						if len(Seq) >= LenBaits: 
							for i in xrange(0,len(Seq),LenCoverage):
								#Keep only part seq with a size of the baits
								if len(Seq[i:i+LenBaits]) == LenBaits:
									if TotalProbes < ArgsDico["Number_baits"]:
										NumberProbe = NumberProbe +1
										BestSeqs.append(DicoSort)
									else: 
										break
								else:
									continue
							TotalProbes = TotalProbes + NumberProbe
						else:
							continue
			else:
				continue
				

#Unique sequence
BestSeqsArray = numpy.array(BestSeqs)
BestSeqsArrayUnique = numpy.unique(BestSeqsArray)
BestSeqsListUnique = BestSeqsArrayUnique.tolist()

#For recap
TargetLength = []
PredictedLength = []
Identity = []
if ArgsDico["intron?"] == "Y":
	#Write output files
	# File 1 : candidates sequences in fasta format
	# File 2 : coding region of candidates sequences in fasta format
	# File 3 : Information about sequences : names, real length, predicted length, identity
	FastaCandidates = open("Candidates.fasta", "w") #File 1
	FastaCandidatesCR = open("CandidatesCodingRegionSeq.fasta", "w") #File 2
	TabInfoSeq = open("AllCandidatesInfo.txt", "w") #File 3
	TabInfoSeq.write("Marker id\tMarker name\tLength\tPredicted length\tPercentage identity\n")
	for DicoCandidate in BestSeqsListUnique:
		TargetLength.append(DicoCandidate["Length"])
		PredictedLength.append(DicoCandidate["Predict Length"])
		Identity.append(DicoCandidate["Percentage identity"])
		#Write File 1
		FastaCandidates.write(">"+str(DicoCandidate["Name"])+"\n"+str(DicoCandidate["Sequence"])+"\n")
		#Write File 2
		CRNumber = DicoCandidate["Coding Region ?"]
		if CRNumber == 1: # 1 CR
			SeqCr = DicoCandidate["CR sequence"]
			FastaCandidatesCR.write(">"+str(DicoCandidate["id"])+" "+str(DicoCandidate["Name"])+"\n"+str(SeqCr)+"\n")	
		else:
			CrSeq = DicoCandidate["CR sequence"]
			for i in xrange(CRNumber):
				FastaCandidatesCR.write(">"+str(DicoCandidate["id"])+" "+str(DicoCandidate["Name"])+"_Exon"+str(i)+"\n"+str(CrSeq[i])+"\n")
		#Write File 3
		PercentageStr = str(DicoCandidate["Percentage identity"])
		TabInfoSeq.write(str(DicoCandidate["id"])+"\t"+str(DicoCandidate["Name"])+"\t"+str(DicoCandidate["Length"])+"\t"+str(DicoCandidate["Predict Length"])+"\t"+PercentageStr+"\n")
	FastaCandidates.close()
	FastaCandidatesCR.close()
	TabInfoSeq.close()
else:
	#If blast data given
	if ArgsDico["Blast_file"] is not None:
		#Write output files
		# File 1 : candidates sequences in fasta format
		# File 2 : Information about sequences : names, real length, predicted length, identity
		FastaCandidates = open("Candidates.fasta", "w") #File 1
		TabInfoSeq = open("AllCandidatesInfo.txt", "w") #File 2
		TabInfoSeq.write("Marker id\tMarker name\tLength\tPredicted length\tPercentage identity\n")
		for DicoCandidate in BestSeqsListUnique:
			TargetLength.append(DicoCandidate["Length"])
			Identity.append(DicoCandidate["Percentage identity"])
			#Write File 1
			FastaCandidates.write(">"+str(DicoCandidate["Name"])+"\n"+str(DicoCandidate["Sequence"])+"\n")
			#Write File 2
			PercentageStr = str(DicoCandidate["Percentage identity"])
			TabInfoSeq.write(str(DicoCandidate["id"])+"\t"+str(DicoCandidate["Name"])+"\t"+str(DicoCandidate["Length"])+"\t"+str(DicoCandidate["Predict Length"])+"\t"+PercentageStr+"\n")
		FastaCandidates.close()
		TabInfoSeq.close()
	else:
		#Write output files
		# File 1 : candidates sequences in fasta format
		# File 2 : Information about sequences : names, real length
		FastaCandidates = open("Candidates.fasta", "w") #File 1
		FastaCandidatesCR = open("CandidatesCodingRegionSeq.fasta", "w") #File 2
		TabInfoSeq = open("AllCandidatesInfo.txt", "w")
		TabInfoSeq.write("Marker id\tMarker name\tMarker length\n")
		for DicoCandidate in BestSeqsListUnique:
			TargetLength.append(DicoCandidate["Length"])
			#Write File 1
			FastaCandidates.write(">"+str(DicoCandidate["Name"])+"\n"+str(DicoCandidate["Sequence"])+"\n")
			#Write File 2
			TabInfoSeq.write(str(DicoCandidate["id"])+"\t"+str(DicoCandidate["Name"])+"\t"+str(DicoCandidate["Length"])+"\n")
		FastaCandidates.close()
		TabInfoSeq.close()


Recap = open("Summary.txt", "w") #File 4
Recap.write("BestMarkers.py\n\n\n\nParameters:\nName input: "+str(ArgsDico["input"])+"\nNumber of baits: "+str(ArgsDico["Number_baits"])+"\nBaits length: "+str(ArgsDico["Length_baits"])+"\nBaits Coverage: "+str(ArgsDico["Coverage_baits"])+"\nMinimal length of sequences: "+str(ArgsDico["length"])+"\nTake introns presence into account for markers selection ? "+str(ArgsDico["intron?"]))
if ArgsDico["Length_intron"] is not None:
	Recap.write(" (Used arbitrary length: "+str(ArgsDico["Length_intron"])+")\n")
else:
	Recap.write(" (Used mean length)\n")
if ArgsDico["only"] is not None:
	Recap.write("Sequences have been sorted according to "+str(ArgsDico["only"])+"\n")
else:
	Recap.write("Sequences have been sorted according to length and similarity. Priority was given to "+str(ArgsDico["priority"])+"\n")
Recap.write("\n\n\nStat summary: \nNumber of sequences before markers selection: "+str(len(DicoFastaFile))+"\nNumber of sequences after markers selection: "+str(len(BestSeqsListUnique))+"\nTarget length: \n\tMax: "+str(max(TargetLength))+"\n\tMin: "+str(min(TargetLength))+"\n\tAverage: "+str(round(numpy.mean(TargetLength)))+"\n\tMedian: "+str(round(numpy.median(TargetLength)))+"\n\tStandart deviation: "+str(round(numpy.std(TargetLength)))+"\n")
if ArgsDico["Blast_file"] is not None:
	if ArgsDico["intron?"] == "Y":
		Recap.write("Predicted length: \n\tMax: "+str(max(PredictedLength))+"\n\tMin: "+str(min(PredictedLength))+"\n\tAverage: "+str(round(numpy.mean(PredictedLength)))+"\n\tMedian: "+str(round(numpy.median(PredictedLength)))+"\n\tStandart deviation: "+str(round(numpy.std(PredictedLength)))+"\nSimilarity: \n\tMax: "+str(max(Identity))+"\n\tMin: "+str(min(Identity))+"\n\tAverage: "+str(round(numpy.mean(Identity)))+"\n\tMedian: "+str(round(numpy.median(Identity)))+"\n\tStandart deviation: "+str(round(numpy.std(Identity))))
	else:
		Recap.write("Similarity: \n\tMax: "+str(max(Identity))+"\n\tMin: "+str(min(Identity))+"\n\tAverage: "+str(round(numpy.mean(Identity)))+"\n\tMedian: "+str(round(numpy.median(Identity)))+"\n\tStandart deviation: "+str(round(numpy.std(Identity))))
Recap.close()
