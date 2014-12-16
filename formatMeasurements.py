from EternaGraph import Motifs as motifSpace
from EternaGraph import Measurement as Measurement
import pickle
import sys
import pdb
from rdatkit.datahandlers import RDATFile
import numpy as np
'''
seq = 'AAACAAAAAUAAUCAAUAAAAAAAGAAAAAUAAA'
dbn1= '...(((...(((...)))..).))...(((.)))'
dbn2= '...(((...((((.))))..).))...(((.)))'
seq2= 'GAACAAAAAUAAUCAAUAAAAAAAGAAAAAUAAA'
dbn3= '...(((...(((...)))..).))...(((.)))'
dbn4= '...(((...((((.))))..).))...(((.)))'
space = motifSpace()

space.addMeasurement([dbn1,dbn2],seq,'m1')
space.addMeasurement([dbn3,dbn4],seq2,'m2')
print 'space:'
space.display()

subSpace = space.subSpace(['A(UCA)A','U(C)A'])
print 'subspace'
subSpace.display()


fullSpace = subSpace.includeKnownMots()
print 'fullSpace:'
fullSpace.display()

print 'structure by motif'
structbyMot = fullSpace.structByMotif(dbn2,seq)
for n in structbyMot:
	print fullSpace.seq[n]
'''

print 'loading Motif Space'
[AllMotifs,allMeasurements]=pickle.load(open('saveMotifSpace.p','rb'))
print 'done loading'


AllMotifs.display()
sys.exit()
Measurements=[]

chosenMots=[
'A(AUUCGU)U',
'A(UUCG)U',
'U(AUUCGU)A',
'C(AUUCGU)G',
'G(AUUCGU)C',
'G(AUUCGU)U'
]
subSpace = AllMotifs.subSpace(chosenMots)
subSpace.thin(minimumNumberOfAppearances = 4)
expandedSpace = subSpace.includeKnownMots()

'''
print 'allMotifs'
for mot in range(0,len(AllMotifs.seq)):
      print 'seq',AllMotifs.seq[mot],'length',AllMotifs.motifLength(mot)

print 'expandedSpace'
for mot in range(0,len(expandedSpace.seq)):
      print 'seq',expandedSpace.seq[mot],'length',expandedSpace.motifLength(mot)
'''


#print 'subspace:'
#subSpace.display()

#print 'Measurements in set'
#print expandedSpace.listMeasurementsInSet()

# loop over measurements and load them into Measurements
fileName = 'None'
for msrmtInSet in expandedSpace.listMeasurementsInSet():
     if fileName != allMeasurements[msrmtInSet][0]:
	fileName = allMeasurements[msrmtInSet][0]
	rdatFileName = 'ETERNA'+fileName[3:12]+'.rdat'
        rdat = RDATFile()
        rdat.load(open('/home/qmac/projects/testdir/'+rdatFileName))
        offset=0
     constructs = rdat.constructs.values()[0]
    # pdb.set_trace()
     dsection = constructs.data[allMeasurements[msrmtInSet][1]]
     if dsection.annotations['sequence'][0] != expandedSpace.sequences[msrmtInSet]:
	print 'Error, sequences not the same!'
	sys.exit()
     seq=dsection.annotations['sequence'][0] 

     countZeros = 0
     rdatLength = len(dsection.values)
     for j in range(1,rdatLength):
	if dsection.values[-j] != 0.0:
	    break
	else:
	    countZeros+=1
     signalLength = rdatLength -countZeros
     #print dsection.annotations['sequence'][0][signalLength-4:signalLength]
     if dsection.annotations['sequence'][0][signalLength-4:signalLength] != 'UUCG':
	if dsection.annotations['sequence'][0][signalLength-4:signalLength] == 'AUUC':
	    countZeros-=1
	    signalLength+=1
	else:
	    print 'Warning: Skipping sequence because it contains unidentified end'
	    continue 
     
	    
     
     aMeasurement = Measurement()
     aMeasurement.StructsDotBracket = expandedSpace.dbns[msrmtInSet]
     aMeasurement.StructsByMotif =\
           [expandedSpace.structByMotif(aMeasurement.StructsDotBracket[0],seq),\
            expandedSpace.structByMotif(aMeasurement.StructsDotBracket[1],seq)]
     aMeasurement.D_obs = np.clip(dsection.values[0:signalLength],a_min=0.0,a_max=1000)
     aMeasurement.RNAstructure_E = allMeasurements[msrmtInSet][2] 
     if len(dsection.values) != len(aMeasurement.StructsDotBracket[1]) and False:
	print 'Warning, there are not the same number of reactivity data points as bases'
	print 'You may need to set an offset'
	print 'number of reactivity values:', len(dsection.values)
	print 'length of dnb -------------:', len(aMeasurement.StructsDotBracket[1])
     Measurements.append(aMeasurement)

pickle.dump([Measurements,expandedSpace],open('chosenMeasurements.p','wb'))
     
    



