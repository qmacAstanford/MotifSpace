from rdatkit.datahandlers import RDATFile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import rdatkit.secondary_structure as ss
import reeffit.map_analysis_utils
from likelihood import *
import cPickle as pickle
import scipy as scipy
import sys
import time
import pickle
from EternaGraph import Motifs as MotifSpace
from EternaGraph import Measurement as Measurement
kT=0.616

#State of the system class
class State:
    def __init__(self,MotifSpace):
	nmot=len(MotifSpace.seq)
	self.E=np.zeros(nmot)
	self.R = []
	for mot in range(0,nmot):
	    self.R.append(zeros(MotifSpace.motifLength(mot)))
        self.post=0
	self.LogProb=-100000000 #This is the posterior
	self.sigma=0.1 #Error in data, this will become an array
	self.prior=0
    def clear(self):
	self.E *=0
	for motR in self.R:
	    motR*=0
	self.LogProb=-10000000
        self.post = 0
    def display(self,MotifSpace):
	for mot in range(0,len(MotifSpace.seq)):
	     print MotifSpace.seq[mot], 'Delta E:',self.E[mot],'R:',self.R[mot] 
    def clipReactivities(self):
	nmot=len(self.R)
	for mot in range(0,nmot):
	     self.R[mot]=np.clip(self.R[mot],a_min=0.00005,a_max=10)
    def reAportionR(self,Motif,DLratio=3):
	#this function needs to be updated

	#DLratio is the reactivity of average dangle over reactivity of average loop
	fl=1.0/(1.0+DLratio)#fraction of reactivity to give to loop 
	fd=1.0-fl
	for key,value in Motif.loops.iteritems():
		dangleNumber=Motif.dangles[key]
		totalR = self.R[value]+self.R[dangleNumber]
		self.R[value]=totalR*fl
		self.R[dangleNumber]=totalR*fd
    def reconstructR(self,R,MotifSpace):
	nmot = len(Motif.seq)
	if len(self.R) != nmot:
	    print 'Error: nmot', nmot,' != len(self.R)',len(self.R)
	for mot in range(0,nmot):
	     self.R[mot]=R[MotifSpace.cumlativeNumR[mot]:MotifSpace.cumlativeNumR[mot+1]]
	     #the above comand turns a list into a list of lists
    def MakeWmtrx(self,msrmts,Motif):
	nmeasurements=len(msrmts)
	nrct = 0
	for anR in self.R: #this for loop is used because not all motifs are the same length
	     nrct=nrct+len(anR)
	if nrct != Motif.cumlativeNumR[-1]:
	     print 'Error: nrct(',nrct ,') != cumlativeNumR,',Motif.cumlativeNumR[-1]
	nequ=0
	for mes in msrmts: #this for loop is used because sequances can have different lengths
	     nequ=nequ+len(mes.D_obs)
	Wmtrx=np.zeros((nequ,nrct))
	frim=0 #first row in measurement
	for mes in range(0,nmeasurements): #seq and mes are used interchanbeably
	     nstr=len(msrmts[mes].StructsByMotif) #this used to be 2
	     if nstr == 0:
		print 'msrmts[mes].StructsByMotif', msrmts[mes].StructsByMotif
		print 'Error: nstr = 0.  Occured on msrmt =',mes
	     Z=0
	     W=np.zeros(nstr)
	     for str in range(0,nstr):
		Energy=msrmts[mes].RNAstructure_E[str]
		for mot in range(0,len(msrmts[mes].StructsByMotif[str])):
		     Energy=Energy+self.E[msrmts[mes].StructsByMotif[str][mot]]
		W[str]=exp(-Energy/kT)
		Z=Z+W[str]
	     for str in range(0,nstr):
		W[str]=W[str]/Z
	     for str in range(0,nstr):
		posInSeq=0
		for mot in msrmts[mes].StructsByMotif[str]:
		     for posInMotif in range(0,len(self.R[mot])):
			Wmtrx[posInSeq+frim][Motif.cumlativeNumR[mot]+posInMotif]=(
			  Wmtrx[posInSeq+frim][Motif.cumlativeNumR[mot]+posInMotif]+W[str])
			posInSeq=posInSeq+1
	     frim=frim+posInSeq
	return Wmtrx
    def logPrior(self,doinglm=0):
	logprior=0
	#this code needs updated to handle more general input
	'''
	prior_mu_E=5.0  #This value needs to be set!
	prior_sigma_E=3.0  #This value needs to be set!
	for j in range(0,len(self.E)):
	     if Motif.type[j] == 'loop':
		     if doinglm:
			  logprior=logprior - (self.E[j] - prior_mu_E)**2/(2*prior_sigma_E**2)
		     else:
			  logprior=(logprior - 
				    (self.E[j] - prior_mu_E)**2/(2*prior_sigma_E**2)
				    -0.5*np.log(2*3.14159*prior_sigma_E**2))
	'''

	self.prior=np.exp(logprior)
	return logprior
    def evaluate(self,msrmts,MCMConly=1):
	self.LogProb=0
	errors=np.zeros(len(msrmts))
	for mes in range(0,len(msrmts)):  #Loops over measurements a.k.a sequences
	     nstr=len(msrmts[mes].StructsByMotif) #number of structures per sequence
	     Eg = np.empty(nstr)  #Energies
	     react = []
	     for struct in range(0,nstr):
		nmot = len(msrmts[mes].StructsByMotif[struct])
		Eg[struct]=msrmts[mes].RNAstructure_E[struct]
		for mot in range(0,nmot):
		     motifnumber=msrmts[mes].StructsByMotif[struct][mot]
		     yy=self.E[motifnumber]
		     Eg[struct] = Eg[struct]+self.E[motifnumber]
		     react.extend(self.R[motifnumber])
	     Z=0
	     W=np.zeros(nstr)
 	     for str in range(0,nstr):
	         W[str]=exp(-Eg[str]/kT)
	         Z=Z+W[str]
	     W=W/Z
	     npos=msrmts[mes].D_obs.shape[0]
	     predicted=zeros(npos)
 	     for str in range(0,nstr):
		 predicted=predicted+W[str]*react[str]
 	     for pos in range(0,npos):
		if MCMConly:
		     self.LogProb=self.LogProb-((predicted[pos]-msrmts[mes].D_obs[pos])**2)/(2*self.sigma**2)-0.5*np.log(2*3.14159*self.sigma**2)
		     #self.sigma can be made sequence and/or position dependent
		else:
		     errors[mes]=errors[mes]-((predicted[pos]-msrmts[mes].D_obs[pos])**2)/(2*self.sigma**2)  	
	if MCMConly:
	   self.LobProb=self.LogProb + self.logPrior()
           self.post=np.exp(self.LogProb)
	   print 'logLike', self.LogProb
	else:
	   print 'norm errors',np.linalg.norm(errors)
	   return errors
# ------------ end of State class---------

class dataForSearch: #this replaces 'stateAndObs' in the old version
	def __init__(self,state,Motif,measurements,energiesToVary):
		self.state=state
		self.msrmts=measurements
		self.energiesToVary=energiesToVary



def searchFunc(E_vec,searchData):
     print 'E_vec', E_vec
     for k in range(0,len(searchData.energiesToVary)):
	searchData.state.E[energiesToVary[k]]=E_vec[k]
     priorError = searchData.state.logPrior(doinglm=1)
     postError = searchData.state.evaluate(searchData.msrmts,MCMConly=0)
     fullError = np.append(priorError,postError)
     #print 'fullError', fullError
     #The sqrt abs comand is necessary to search the same distriution as MCMC
     return  np.sqrt(np.absolute(fullError))

#pickle.dump([Measurements,expandedSpace],open('chosenMeasurements.p','wb'))
[measurements,Motif]=pickle.load(open('chosenMeasurements.p','rb'))

D_obs_list = []
for msrmt in measurements:
      D_obs_list.extend(msrmt.D_obs)
D_obs_vec = np.array(D_obs_list)

energiesToVary=[] #to tell the search algorithmn which energies to vary
for key,value in Motif.loops.iteritems():# e.g. only vary loop energies
     energiesToVary.append(value)


dofullSearch=1
if dofullSearch:
     E_0=np.zeros(len(energiesToVary))
     #E_0=np.random.normal(loc=5.0,scale=1.0,size=len(Motif.loops))
     
     searchState=State(Motif)
     searchData = dataForSearch(searchState,Motif,measurements,energiesToVary)
     for iter in range(0,8):     
	print '------------------Starting linalg------------------'
	Wmtrx=searchData.state.MakeWmtrx(measurements,Motif)
	#(Rout,residuals,rank,singular)=np.linalg.lstsq(Wmtrx,D_obs_vec)
	(Rout,reactivityErrorNorm)=scipy.optimize.nnls(Wmtrx,D_obs_vec)
	print 'Reactivity Error:', reactivityErrorNorm
	searchState.reconstructR(Rout,Motif)
	searchState.display(Motif)
	searchState.clipReactivities()
	
	#if iter == 0:
	#     searchState.reAportionR(Motif,DLratio=2.0)
	#     searchState.display(Motif)

	
	print '-----------------Starting the search-----------------'
	(resultE,exitReason) = scipy.optimize.leastsq(searchFunc,E_0,args=searchData,maxfev=500,factor=0.5,xtol=0.0005)
	print 'Exit Reason:', exitReason
	rememberOptimalEnergies=0
	if rememberOptimalEnergies:
	     for k in range(0,len(searchData.energiesToVary)): 
		searchData.state.E[energiesToVary[k]]=resultE[k]
		E_0[k]=resultE[k]
	     E_0+=np.random.normal(loc=0.0,scale=0.5,size=E_0.shape)
	else:
	     for k in range(0,len(searchData.energiesToVary)): 
		searchData.state.E[energiesToVary[k]]=resultE[k]
	#searchState.display(Motif)
searchState.display(Motif)

doMCMC=0
if doMCMC:
	nmoves=8000
	#states=[State() for i in range(0,nmoves)]
	nextState=State(Motif)
	nextState.clear()
	nkeep=0
	
	def changeMtrx(maxstep=1.1,shape=1):
	    factors=np.random.uniform(low=1,high=maxstep,size=shape)
	    choiceMtrx=np.random.randint(2,size=shape)
	    factors = choiceMtrx*factors + (np.ones(shape)-choiceMtrx)/factors
	    return factors

	def MoveEnergies(currentState,nextState,stepSize):
	    for key, value in currentState.loopsE.iteritems():
		nextState.loopsE[key]=value+np.random.normal(loc=0,scale=stepSize,size=1)

	def MoveR(currentState,nextState,stepSize):
	    #This function isn't written yet
	    #Data set must have many of the same 
	    donothing=1

        def MoveEandR(currentState,nextState,Estep,Rstep):
	    #print 'currentState.R', currentState.R
	    #print 'currentState.E', currentState.E
	    nextState.clear()
	    for mot in range(0,len(currentState.R)):
		#nextState.R.append( np.random.normal(loc=0,scale=Rstep,size=currentState.R[mot].shape) + currentState.R[mot])
		nextState.R.append(currentState.R[mot]*changeMtrx(1.0+Rstep,currentState.R[mot].shape))
	    nextState.E = np.random.normal(loc=0,scale=Estep,size=currentState.E.shape) + currentState.E

	track = np.zeros([nmoves,3])
	#MCMC
	searchState.evaluate(measurements,MCMConly=1)
	for step in range(0,nmoves-1):
	     MoveEandR(searchState,nextState,0.025,0.075)
	     nextState.evaluate(measurements)
	     difference = nextState.LogProb - searchState.LogProb
	     if difference>0:
		print 'Keeping because more likely'
		nkeep+=1
		workptr=searchState
		searchState=nextState
	 	nextState=workptr
		nextState.clear()
	     else:
		ratio=exp(difference)
		randNumb = np.random.uniform(low=0,high=1)
		if (ratio<randNumb):
		   print 'reject by chance because ',ratio,'<',randNumb
		else:
		   print 'Keep by chance because ',ratio,'>',randNumb
		   nkeep+=1
		   workptr=searchState
		   searchState=nextState
	 	   nextState=workptr
		   nextState.clear()
	     track[step,:]=[searchState.R[2][1],searchState.R[3][2],searchState.R[1][3]]
        print 'kept',nkeep,'of',nmoves		
	searchState.display(Motif)
	
	
	plt.figure(1)
	plt.plot(track[:,0])
	plt.plot(track[:,1])
	plt.plot(track[:,2])
	plt.savefig('trackMCMC.png',dpi=256)

	sys.exit()
	plt.figure(1)
	H,xedges,yedges=np.histogram2d(Data[:,0],Data[:,1],bins=10)
	H=np.flipud(H)
	print H
	plt.xlabel('AAAA')
	plt.ylabel('GAAA')
	plt.imshow(H,vmin=0,vmax=H.max(),cmap=plt.get_cmap('coolwarm'),extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
	plt.savefig('post.png',dpi=100)
	#np.save('FinalD_fac',states[nmoves-1].D_fac)
	#np.save('FinalPsi',states[nmoves-1].Psi)

