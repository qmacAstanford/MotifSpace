from rdatkit.datahandlers import RDATFile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import rdatkit.secondary_structure as ss
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infile",help="input file name, please end with .rdat")
parser.add_argument("outfile",help="output file name no extition needed")
args = parser.parse_args()


#import rdat data
rdat = RDATFile()
rdat.load(open('/home/qmac/projects/testdir/'+args.infile))
offset=0
constructs = rdat.constructs.values()[0]
competing_pairs = []
sequences_included=[]
msrmtsNumbers=[]
print 'lenth of constructs.data', len(constructs.data)

for count in range(0,len(constructs.data)):
     dsection = constructs.data[count]
     seq=dsection.annotations['sequence'][0]
     #structs=ss.fold(seq,nstructs=2)
     struct_energy_list =[(struct.dbn, energy) for struct, energy in zip(*ss.subopt(seq,nstructs=100,fraction=0.075,energies=True))]
     struct_energy_list_unique = list(set(struct_energy_list))
     struct_energy_list_unique = sorted(struct_energy_list_unique, key=lambda x: x[1])
     #print struct_energy_list_unique 
     #print 'length',len(struct_energy_list_unique)
     if len(struct_energy_list_unique) > 1:
	 competing_pairs.append(struct_energy_list_unique[0:2])
	 sequences_included.append(seq)
	 msrmtsNumbers.append(count)
     #print competing_pairs
     if count > 20000:
	break
pickle.dump([competing_pairs,sequences_included,msrmtsNumbers],open(args.outfile+".p","wb"))
