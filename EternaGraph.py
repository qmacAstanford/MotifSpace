from rdatkit.datahandlers import RDATFile
import numpy as np
import rdatkit.secondary_structure as ss
import sys

class Motifs:
     '''
     Motifs is a class for storing all the different types of motifs that appear in a problem. For example a tetraloop with sequence A(GAGA)U would be a motif.  Each motif is associated with a index number.  This index number is where it sequence is stored in Motifs.seq.  Likewise with Motifs.type, Motifs.msrmtsContaining, Motifs.motifInSet.  MotifTypes is a dictionary of dictionaries for look up a index number based on sequence.
     '''
     def __init__(self,motsConsid={'loops':True,'x2wayjunctions':True,'x3wayjunctions':True,'x4wayjunctions':True,'bulges':True,'pairs':False,'sstrand':False,'dangles':False}):
	self.dbns = {} # one list per measurement, multiple dbns in list
	self.sequences = {} #one per measurement
	self.msrmtsInSet = {} #True or false
	self.msrmtLabel = {} # one per sequence
	self.motifsConsidered=motsConsid
	self.seq = []  #an array of sequences
	self.type = [] #an array of words like "loop" and "helx"
	self.msrmtEdge = []
	self.motifInSet=[]
	self.loops = {} #a dictionary to go from sequence to index
	self.x2wayjunctions = {}
	self.x3wayjunctions = {}
	self.x4wayjunctions = {}
	self.bulges = {}
        self.pairs = {}
	self.sstrand = {}
	self.dangles = {}
	self.other = {}
	self.motifTypes = {'loops':self.loops,
			   'x2wayjunctions':self.x2wayjunctions,
			   'x3wayjunctions':self.x3wayjunctions,
			   'x4wayjunctions':self.x4wayjunctions,
			   'bulges':self.bulges,
			   'pairs':self.pairs,
			   'sstrand':self.sstrand,
			   'dangles':self.dangles,
			   'other':self.other}
	self.cumlativeNumR=[0]
     def add(self,seq,type,msrmtsNumber,length): #seq may have more charactors then length 
	self.msrmtsInSet[msrmtsNumber]=True
	mot =0
	if seq in self.motifTypes[type]:
	     exflag=1
	     mot = self.motifTypes[type][seq]
	     self.msrmtEdge[mot].append(msrmtsNumber)
	     return mot
	else:
	     self.motifTypes[type][seq]=len(self.seq) #add to dictionary
	     self.msrmtEdge.append([msrmtsNumber])
	     self.seq.append(seq)
	     self.motifInSet.append(1)
	     self.type.append(type)
	     self.cumlativeNumR.append(self.cumlativeNumR[-1]+length)
	     if len(self.msrmtEdge) != len(self.seq):
		print 'Warning, msrmtEdge is not the same lenght as seq'
		sys.exit()
	     if len(self.cumlativeNumR)-1 != len(self.seq):
		print 'Warning, wrong length of culativeNumR'
		sys.exit()
	     if type == 'loops' and len(seq)-4 != length:
		print 'tried to add loop', seq, 'with length', length
		sys.exit()
	     if type == 'bulges' and len(seq)-4 != length:
		print 'tried to add bulge', seq,'with length',length
		sys.exit()
	     if type == 'x2wayjunctions' and len(seq)-5 != length:
		print 'tried to add 2way', seq,'with length',length
	     return len(self.seq)-1
     def listMeasurementsInSet(self):
	#write now this lists all measurements in space
	msrmts = []
	for mot in range(0,len(self.seq)):
	     for msrmt in self.msrmtEdge[mot]:
		msrmts.append(msrmt)
	return sorted(set(msrmts))
     def motifLength(self,mot):
	return self.cumlativeNumR[mot+1]-self.cumlativeNumR[mot]
     def subSpace(self,chosenMotifs):
	subSpace = Motifs()
	for mot in range(0,len(self.seq)):
	     if self.motifInSet[mot] and (self.seq[mot] in chosenMotifs):
		for msrmt in self.msrmtEdge[mot]:
		    if self.msrmtsInSet[msrmt] and msrmt not in subSpace.msrmtsInSet:
			subSpace.addMeasurement(self.dbns[msrmt],self.sequences[msrmt],msrmt,Label=self.msrmtLabel[msrmt])
	subSpace.thin()
	return subSpace
     def includeKnownMots(self):
	outSpace = Motifs(motsConsid={'loops':True,'x2wayjunctions':True,'x3wayjunctions':True,'x4wayjunctions':True,'bulges':True,'pairs':True,'sstrand':True,'dangles':True})
	for mot in range(0,len(self.seq)):
	     if self.motifInSet[mot]:
		for msrmt in self.msrmtEdge[mot]:
		    if self.msrmtsInSet[msrmt] and msrmt not in outSpace.msrmtsInSet:
			outSpace.addMeasurement(self.dbns[msrmt],self.sequences[msrmt],msrmt,Label=self.msrmtLabel[msrmt])
	return outSpace
     def add_symmetric_difference(self,MotifA,MotifB,include_differnt_n_occur=True):
	if include_differnt_n_occur:
		#this assumes MotifA and MotifB came from a single mesurment
		for j in range(0,len(MotifA.seq)):
			motSeq=MotifA.seq[j]
			if motSeq is None:
				print 'Error, strange None detected in A'
				print 'j=',j,'len(MotifA.seq)=',len(MotifA.seq)
				MotifA.display()
				sys.exit()
			mot_in_A=j
			mot_in_B=MotifB.motifTypes[MotifA.type[j]].get(motSeq)
			if (mot_in_B is None) or len(MotifA.msrmtEdge[mot_in_A]) !=\
						 len(MotifB.msrmtEdge[mot_in_B]):
			    self.add(motSeq,MotifA.type[j],MotifA.msrmtEdge[j][0],
				 MotifA.motifLength(j))
		for j in range(0,len(MotifB.seq)):
			motSeq=MotifB.seq[j]
			if motSeq is None:
				print 'Error, strange None detected in B'
				sys.exit()
			mot_in_B=j
			mot_in_A=MotifA.motifTypes[MotifB.type[j]].get(motSeq)
			if (mot_in_A is None) or len(MotifA.msrmtEdge[mot_in_A]) !=\
						 len(MotifB.msrmtEdge[mot_in_B]):
			    self.add(motSeq,MotifB.type[j],MotifB.msrmtEdge[j][0],
				 MotifB.motifLength(j))
	else:	
		for j in range(0,len(MotifA.seq)):
		    if MotifA.seq[j] not in MotifB.motifTypes[MotifA.type[j]]:
			self.add(MotifA.seq[j],MotifA.type[j],MotifA.msrmtEdge[j][0],
				 MotifA.motifLength(j))
		for j in range(0,len(MotifB.seq)):
		    if MotifB.seq[j] not in MotifA.motifTypes[MotifB.type[j]]:
			self.add(MotifB.seq[j],MotifB.type[j],MotifB.msrmtEdge[j][0],
				 MotifB.motifLength(j))
     def structByMotif(self,dbn,seq):
	self.add('~','other',-1,1) #this is a blank motif for things we don't have data for
	struct_by_motif=[]
	first_base_numbers=[]
	astruct=ss.SecondaryStructure(dbn)
	mots = astruct.explode()
	if ('hairpins' in mots):
	    for aloop in mots['hairpins']:
		 loopseq=self.changeMotifFormat(seq,'hairpins',aloop)
		 if loopseq in self.loops:
			struct_by_motif.append(self.loops[loopseq])
			first_base_numbers.append(aloop[0])
		 else:
		 	for j in aloop:
			     struct_by_motif.append(self.other['~'])
			     first_base_numbers.append(j)
	if ('interiorloops' in mots):
	    for aloop in mots['interiorloops']:
		 loopseq=self.changeMotifFormat(seq,'interiorloops',aloop)
		 # ---  Warning ---False in below statment always neglects x3ways
		 if loopseq in self.x2wayjunctions and False:
			struct_by_motif.append(self.x2wayjunctions[loopseq])
			first_base_numbers.append(aloop[0])
		 else:
		 	for j in range(0,len(aloop)):
			     struct_by_motif.append(self.other['~'])
			     first_base_numbers.append(j)
	if ('3wayjunctions' in mots):
	    tuples_3way=[tuple(l) for l in mots['3wayjunctions']] #this line controlls for explode error
	    for aloop in set(tuples_3way):
		 loopseq=self.changeMotifFormat(seq,'3wayjunctions',aloop)
		 # ---  Warning ---False in below statment always neglects x3ways
		 if loopseq in self.x3wayjunctions and False:
			struct_by_motif.append(self.x3wyajunctions[loopseq])
			first_base_numbers.append(aloop[0])
		 else:
		 	for j in range(0,len(aloop)):
			     struct_by_motif.append(self.other['~'])
			     first_base_numbers.append(j)
	if ('4wayjunctions' in mots):
	    tuples_4way=[tuple(l) for l in mots['4wayjunctions']] #this line controlls for explode error
	    for aloop in set(tuples_4way):
		 loopseq=self.changeMotifFormat(seq,'4wayjunctions',aloop)
		 # ---  Warning ---False in below statment always neglects x3ways
		 if loopseq in self.x4wayjunctions and False:
			struct_by_motif.append(self.x4wyajunctions[loopseq])
			first_base_numbers.append(aloop[0])
		 else:
		 	for j in range(0,len(aloop)):
			     struct_by_motif.append(self.other['~'])
			     first_base_numbers.append(j)
	if ('5wayjunctions' in mots):
	    tuples_5way=[tuple(l) for l in mots['5wayjunctions']] #this line controlls for explode error
	    for aloop in set(tuples_5way):
		for j in range(0,len(aloop)):
		     struct_by_motif.append(self.other['~'])
		     first_base_numbers.append(j)

	if ('bulges' in mots):	
	    for aloop in mots['bulges']:
		 bulgseq=self.changeMotifFormat(seq,'bulges',aloop)
		 if loopseq in self.bulges:
			struct_by_motif.append(self.bulges[loopseq])
			first_base_numbers.append(aloop[0])
		 else:
		 	for j in range(0,len(aloop)):
			     struct_by_motif.append(self.other['~'])
			     first_base_numbers.append(j)
	if ('helices' in mots):
	    for ahelx in mots['helices']:
		for bace in range(0,len(ahelx)/2):
		    fivePrimePosition = ahelx[bace*2]
		    threePrimePosition = ahelx[bace*2+1]
		    fivePrimeBace = seq[fivePrimePosition]
		    threePrimeBace = seq[threePrimePosition]
		    loopseq = '*'+fivePrimeBace+threePrimeBace
		    if loopseq in self.pairs:
			struct_by_motif.append(self.pairs[loopseq])
			first_base_numbers.append(fivePrimePosition)
		    else:
			struct_by_motif.append(self.other['~'])
			first_base_numbers.append(fivePrimePosition)
		    loopseq = fivePrimeBace+threePrimeBace+'*'
		    if loopseq in self.pairs:
			struct_by_motif.append(self.pairs[loopseq])
			first_base_numbers.append(threePrimePosition)
		    else:
			struct_by_motif.append(self.other['~'])
			first_base_numbers.append(threePrimePosition)
	if ('dangles' in mots):
	    for aloop in mots['dangles']:
		loopseq=self.changeMotifFormat(seq,'dangles',aloop)
		if loopseq in self.dangles:
			struct_by_motif.append(self.dangles[loopseq])
			first_base_numbers.append(aloop[0])
		else:
		 	for j in aloop:
			     struct_by_motif.append(self.other['~'])
	if ('sstrand' in mots):
	    for aloop in mots['sstrand']:
		 loopseq=self.changeMotifFormat(seq,'sstrand',aloop)
		 if loopseq in self.sstrand:
			struct_by_motif.append(self.sstrand[loopseq])
			first_base_numbers.append(aloop[0])
		 else:
		 	for j in aloop:
			     struct_by_motif.append(self.other['~'])
			     first_base_numbers.append(j)
	TotalLength = 0
	for mot in struct_by_motif:
		TotalLength+=self.motifLength(mot)
	if TotalLength != len(dbn):
		print 'Error in structByMotif: some motif is missing'
		print 'dbn',dbn
		print [x for (y,x) in sorted(zip(first_base_numbers,struct_by_motif))]
		print ss.SecondaryStructure(dbn).explode()
		sys.exit()
	return [x for (y,x) in sorted(zip(first_base_numbers,struct_by_motif))]
	#make n way junctions blanks for the time being
     def thin(self,minimumNumberOfAppearances = 2):
	continue_thinning=True
	while continue_thinning:
	    continue_thinning=False
	    for mot in range(0,len(self.seq)):
		if self.motifInSet[mot]:
		    nEdges=0
		    for edge in self.msrmtEdge[mot]:
			if self.msrmtsInSet[edge]:
			     nEdges +=1
		    if nEdges < minimumNumberOfAppearances:
			for edge in self.msrmtEdge[mot]:
			     self.msrmtsInSet[edge] = False
			self.motifInSet[mot] = False
     			continue_thinning=True
     def changeMotifFormat(self,seq,type,baseList,closingPairs=1):
	if type == 'sstrand' or type == 'sstrand':
	     aloop=baseList
	     return seq[aloop[0]-closingPairs:aloop[0]]+'-'+\
		 		seq[aloop[0]:aloop[-1]+1]+'-'+\
				seq[aloop[-1]+1:aloop[-1]+1+closingPairs]
	elif type == 'hairpins' or type == 'loops':
	     aloop=baseList
	     return seq[aloop[0]-closingPairs:aloop[0]]+'('+\
		 		seq[aloop[0]:aloop[-1]+1]+')'+\
				seq[aloop[-1]+1:aloop[-1]+1+closingPairs]
	elif type == 'x2wayjunctions' or type == 'interiorloops':
	     breakpoint=0
	     aloop = baseList
	     for bace in range(0,len(aloop)-1):
		 if aloop[bace]+1 != aloop[bace+1]:
			 breakpoint=bace+1
			 break
	     return seq[aloop[0]-closingPairs:aloop[0]+\
			breakpoint+closingPairs]+'x'+\
			seq[aloop[breakpoint]-closingPairs:aloop[-1]+1+closingPairs]
	elif type == 'x3wayjunctions' or type =='3wayjunctions':
	    aloop=baseList
	    breakpoint=0
	    breakpoint2=0
	    for bace in range(0,len(aloop)-1):
	       if aloop[bace]+1 != aloop[bace+1]:
	   	 breakpoint=bace+1
	   	 break
	    for bace in range(breakpoint,len(aloop)-1):
	       if aloop[bace]+1 != aloop[bace+1]:
	   	 breakpoint2=bace+1
	   	 break
	    loopseq=(
	       seq[aloop[0]-closingPairs:aloop[0]+breakpoint+closingPairs]+
	       'x'+
	       seq[aloop[breakpoint]-closingPairs:
	   	    aloop[breakpoint]+breakpoint2-breakpoint+closingPairs]+
	       'x'+
	       seq[aloop[breakpoint2]-closingPairs:aloop[-1]+1+closingPairs]
	       )
	    return loopseq
	elif type == 'x4wayjunctions' or type =='4wayjunctions':
	    aloop=baseList
	    breakpoint=0
	    breakpoint2=0
	    breakpoint3=0
	    for bace in range(0,len(aloop)-1):
	       if aloop[bace]+1 != aloop[bace+1]:
	   	 breakpoint=bace+1
	   	 break
	    for bace in range(breakpoint,len(aloop)-1):
	       if aloop[bace]+1 != aloop[bace+1]:
	   	 breakpoint2=bace+1
	   	 break
	    for bace in range(breakpoint2,len(aloop)-1):
	       if aloop[bace]+1 != aloop[bace+1]:
	   	 breakpoint3=bace+1
	   	 break
	    loopseq=(
	       seq[aloop[0]-closingPairs:aloop[0]+breakpoint+closingPairs]+
	       'x'+
	       seq[aloop[breakpoint]-closingPairs:
	   	    aloop[breakpoint]+breakpoint2-breakpoint+closingPairs]+
	       'x'+
	       seq[aloop[breakpoint2]-closingPairs:
	   	    aloop[breakpoint2]+breakpoint3-breakpoint2+closingPairs]+
	       'x'+
	       seq[aloop[breakpoint3]-closingPairs:aloop[-1]+1+closingPairs]
	       )
	    return loopseq	
	elif type=='bulges':
	    bulge=baseList
	    loopseq = seq[bulge[0]-closingPairs:bulge[0]]+'['+\
		seq[bulge[0]:bulge[-1]+1]+']'+\
		seq[bulge[-1]+1:bulge[-1]+1+closingPairs]
	    return loopseq
	elif type =='helices':
	      print 'Error, helices not covered'
	      sys.exit()
	elif type == 'dangles':
	     adangle=baseList
	     if adangle[0]==0:
	        return "5'"+seq[adangle[0]:adangle[-1]+1]
	     else:
	        return seq[adangle[0]:adangle[-1]+1]+"3'"
	else:
	     raise ValueError('Type '+type+' not covered')
     def addMeasurement(self,dbns,seq,msrmtID,Label='-'):
	    self.dbns[msrmtID]=dbns
	    self.sequences[msrmtID]=seq
	    self.msrmtLabel[msrmtID]=Label
	    self.msrmtsInSet[msrmtID]=False
	    competingMotifs=[Motifs(),Motifs()]  #does this get deleted every time in garbage
	    for subStruct in range(0,len(dbns)):
		Motif = competingMotifs[subStruct]
		astruct=ss.SecondaryStructure(dbns[subStruct])
		mots = astruct.explode()
		if self.motifsConsidered['loops'] and ('hairpins' in mots):
		    for aloop in mots['hairpins']:
			 loopseq=self.changeMotifFormat(seq,'hairpins',aloop)
			 #print 'adding', loopseq, 'with length', len(aloop)
			 Motif.add(loopseq,'loops',msrmtID,len(aloop))
		if self.motifsConsidered['x2wayjunctions'] and ('interiorloops' in mots):
		    for aloop in mots['interiorloops']:
			 loopseq=self.changeMotifFormat(seq,'interiorloops',aloop)
			 Motif.add(loopseq,'x2wayjunctions',msrmtID,len(aloop))
		if self.motifsConsidered['x3wayjunctions'] and ('3wayjunctions' in mots):
		    for aloop in mots['3wayjunctions']:
			 loopseq=self.changeMotifFormat(seq,'3wayjunctions',aloop)
			 Motif.add(loopseq,'x3wayjunctions',msrmtID,len(aloop))
		if self.motifsConsidered['x4wayjunctions'] and ('4wayjunctions' in mots):
		    for aloop in mots['4wayjunctions']:
			 loopseq=self.changeMotifFormat(seq,'4wayjunctions',aloop)
			 Motif.add(loopseq,'x4wayjunctions',msrmtID,len(aloop))
		if self.motifsConsidered['bulges'] and ('bulges' in mots):	
		    for bulge in mots['bulges']:
			 bulgseq=self.changeMotifFormat(seq,'bulges',bulge)
			 Motif.add(bulgseq,'bulges',msrmtID,len(bulge))
		if self.motifsConsidered['pairs'] and ('helices' in mots):
		    for ahelx in mots['helices']:
			for bace in range(0,len(ahelx)/2):
			    bacePair=seq[ahelx[bace*2]]+seq[ahelx[bace*2+1]]
			    Motif.add('*'+bacePair,'pairs',msrmtID,1)
			    Motif.add(bacePair+'*','pairs',msrmtID,1)
		if self.motifsConsidered['dangles'] and ('dangles' in mots):
		    for adangle in mots['dangles']:
			loopseq=self.changeMotifFormat(seq,'dangles',adangle)
			Motif.add(loopseq,'dangles',msrmtID,len(adangle))
	    self.add_symmetric_difference(competingMotifs[0],competingMotifs[1])
     def display(self):			
	#print self.seq
	#print self.type
	#print self.msrmtEdge      
	#print 'msrmtsInSet', self.msrmtsInSet
        #print 'motifInSet', self.motifInSet
	if 0:
		for mot in range(0,len(self.seq)):
		   if self.motifInSet[mot]:
			 print self.type[mot],self.seq[mot],' in measuremnts [',
			 for mes in self.msrmtEdge[mot]:
				if self.msrmtsInSet[mes]:
				    print mes,',',
			 print ']'
	elif 1:
		print 'Source','Target'
		for mot in range(0,len(self.seq)):
		   if self.motifInSet[mot]:
			 for mes in self.msrmtEdge[mot]:
			    if self.msrmtsInSet[mes]:
			       print self.seq[mot],mes
	else:
		print 'ID','Label'
		for mes, Label in self.msrmtLabel.iteritems():
		   if self.msrmtsInSet[mes]:
			 print mes, Label
		for mot in range(0,len(self.seq)):
		   if self.motifInSet[mot]:
			 print self.seq[mot], self.seq[mot]	
'''
			       if self.type[mot]=='loops':
				  print self.seq[mot][0]+'('+self.seq[mot][1:-1]+')'+self.seq[mot][-1], mes
			       elif self.type[mot]=='bulges':
				  print self.seq[mot][0]+'['+self.seq[mot][1:-1]+']'+self.seq[mot][-1], mes
			       else:
				  print self.seq[mot],mes
'''
class Measurement:
     def __init__(self):
	self.StructsByMotif = []
	self.StructsDotBracket = []
	self.D_obs = []
	self.RNAstructure_E = []
