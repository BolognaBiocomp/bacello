#!/usr/bin/env python
import os

BACELLO_HOME=os.environ.get("BACELLO_HOME")

try:
	import psyco
	psyco.full()
except ImportError:
	pass

import sys
import numpy
import os


#enviroment variables
#if CPU=='auto':
#	r,w=popen2.popen2("cat /proc/cpuinfo | grep '^processor'")
#	CPU=len(r.readlines())
	

def read1Fasta(fname):
    ''' read1Fasta(fname) return a sequence
    '''
    import sys
    try:
        lines=open(fname,'r').readlines()
    except:
        sys.stderr.write('cannot open file '+fname+'\n')
        sys.exit()
    seq=''
    for line in lines:
        if line[0]!='>': # assuming fasta file
           seq+=''.join(line.split())
    return seq


def readProf(fname,aaOrder='VLIMFWYGAPSTCHRKQEND'):
    ''' readProf(fname,aaOrder) returns a list with the profile elements and the residue order
        return aaOrder,prof[], otherfeatures[]
    '''
    try:
        lines=open(fname,'r').readlines()
    except:
        sys.stderr.write('cannot open file '+fname+'\n')
        sys.exit()
    prof=[]
    other_features=[]
    i=0
    while lines[i].find("#")>=0 and i <len(lines):
        pos=lines[i].find("ALPHABET=")
        if pos>0:
           v=lines[i][pos+len("ALPHABET="):].split()
           aaOrder=''.join(v[:20])
        i+=1
    if i >= len(lines):
        sys.stderr.write('Error in '+fname+'\n')
        sys.exit()
    pos=0
    while i<len(lines):
        v=lines[i].split()
        vect=[float(v[j]) for j in range(len(aaOrder))]
        prof.append(vect)
        if len(v) > len(vect):
            other_features.append(v[len(vect):])
        i+=1
    return aaOrder,prof,other_features



#----------- ANDREA PIERLEONI--------------#
def cal_freq_seq(sequenza, fin=0): 
	Laa = ['K','R','H',   'D','E',   'N','Q','S','T','Y',   'A','V','L','I','P','F','M','W','G','C']
	if fin > 0:
		sequenza=sequenza[:fin]
	elif fin < 0:
		sequenza=sequenza[fin:]
	else:
		pass
	Lrisultato=[]
	tot=float(len(sequenza))
	for aa in Laa:
		Lrisultato.append(sequenza.count(aa)/tot)
	return Lrisultato

def cal_freq_prof(profilo, fin=0,baco=0):
	'''calcola la frequenza (come vettore da 20 aminoacidi) di un profilo
	  modificare questa funzione se l'input del profilo e' in un formato diverso
	   FORMATO PROFILO ATTUALE:
	array numpy:
	primo livello: posizioni nella sequenza dall Nter al Cter
	secondo livello: frequenza aminocidica nella posizione, secondo l'ordine stabilito da Laa
	fin:	definisce la porzione di sequenza di cui calcolare la frequenza aminoacidica
		0---> sequenza intera
		numero positivo intero---> da N-ternimale a posizione specificata
		numero negativo intero--->da C-ternimale a posizione specificata
		
	baco= implementa il calcolo della frequenza usato nella prima versione di bacello'''
	tot=float(len(profilo))
	if fin > 0:
		profilo=profilo[:fin]
	elif fin < 0:
		profilo=profilo[fin:]
	else:
		pass
	Lrisultato=[]
	if baco:
		Lrisultato.extend(profilo.sum(0))
		for i in xrange(len(Lrisultato)):
			Lrisultato[i]=Lrisultato[i]/tot
	else:
		Lrisultato.extend(profilo.mean(0))
	return Lrisultato



def bacello(sequenza, profilo, classe='A'):
	''' versione light di BaCelLo:
		prende in input una sequenza aminoacidica il suo profilo e la classe di appartenza per il regno eucariotico:
		'A'===>Animali(Metazoa)
		'F'===>Funghi(Fungi)
		'P'===>Piante(Viridiplantae)
		
		da in output la localizzazione predetta'''
		
	if sequenza[0]=='M':#do not consider initial met
		sequenza=sequenza[1:]
	if classe=='P':
		Lmodelli=[os.path.join(BACELLO_HOME, 'MOD/M1'),
                          os.path.join(BACELLO_HOME, 'MOD/M2'),
                          os.path.join(BACELLO_HOME, 'MOD/M3'),
                          os.path.join(BACELLO_HOME, 'MOD/M4')]#nomi dei file dei 4 modelli di svm-light
		Lcodifiche=[[0,20,40,60],[0,20,40,60,-20,-50,-100],[0],[0]]#codice di codifica per i vari svm
		Lsoglia_med=[0.442, -0.201, 0.375, -0.327]# parametri di bilanciamento per svm
	elif (classe=='A' or classe=='F'):
		Lmodelli=[os.path.join(BACELLO_HOME, 'MOD/M1'),
                          os.path.join(BACELLO_HOME, 'MOD/M2'),
                          os.path.join(BACELLO_HOME, 'MOD/M3')]
		Lcodifiche=[[0,20,40,60],[0,20,40,60,-20,-50,-100],[0],[0]]
		if classe=='A':
			Lsoglia_med=[0.122, 0.467, -0.568]
		else:
			Lsoglia_med=[0.571, 0.344, -0.644]
	#Preparazione input (calcolo frequenze aminoacidiche di sequenza e profilo)
	input_completo=[]
	for codice in [0,20,40,60,-20,-50,-100]:
		input_completo.extend(cal_freq_prof(profilo,codice,baco=1))
		input_completo.extend(cal_freq_seq(sequenza,codice))
	Llimiti=[160,280,40,40]#lunghezza dei vettori di input per i vari modelli
	Linput=[]#contiene gli input dei vari modelli
	for k in xrange(len(Lmodelli)):
		Linput.append(numpy.array(input_completo[:Llimiti[k]]))
	#Predizione con Kernel RBF (Classe di piero)
	Lrisultati_svm=[]#contiene il risultato delle predizioni per ogni modello
	for k in xrange(len(Lmodelli)):
		svm=getSVMLight(Lmodelli[k]+'_'+classe)
		svmout=svm.predict(Linput[k])+Lsoglia_med[k]
		Lrisultati_svm.append(svmout)#calcolo con kernel RBF + aggiunta soglia di bilanciamento
	#Scelta localizzazione predetta
        # scores = numpy.array(Lrisultati_svm)/numpy.sum(numpy.array(Lrisultati_svm))
        x = numpy.array(Lrisultati_svm)
        scores = (x - numpy.min(x))/(numpy.max(x) - numpy.min(x))
        #print scores
        score = 0.0
        second = None
        
	if classe=='P':
		if Lrisultati_svm[0]>=0:
			prediction='Secretory'
			score = scores[0]
		else:
			if Lrisultati_svm[1]>=0:
				if Lrisultati_svm[3]>=0:
					prediction='Mitochondrion'
					score = scores[3]
				else:
					prediction='Chloroplast'
					score = 1 - scores[3]
			else:	
				if Lrisultati_svm[2]>=0:
					prediction='Nucleus'
					score = scores[2]
				else:
					prediction='Cytoplasm'
					score = 1 - scores[2]
	else:#animali e funghi
		if Lrisultati_svm[0]>=0:
			prediction='Secretory'
                        score = scores[0]
		else:
			if Lrisultati_svm[1]>=0:
				prediction='Mitochondrion'
                                score = scores[1]
			else:	
				if Lrisultati_svm[2]>=0:
					prediction='Nucleus'
                                        score = scores[2]
				else:
					prediction='Cytoplasm'
                                        score = 1 - scores[2]
	return prediction, score







#------------PIERO FARISELLI--------------#
class SVMLike:
    ''' This class implement the prediction phase of a SVM
        based on the model of svm light.
        Presently, there are several limitation including:
           - fixed dimension of the input vectors () 
           - !!! only rbf kernel implemente up to now !!!!
    '''

    def __init__(self,sv,ai,b,kernel,params):
        '''__init__(self,sv,ai,b,kerneltype,params)
               sv     -> support vectors
               ai     -> lagrange multipliers
               b      -> shifting threshold
               kernel -> 0: linear (default)
                         1: polynomial (s a*b+c)^d
                         2: radial basis function exp(-gamma ||a-b||^2)
                         3: sigmoid tanh(s a*b + c)
               params -> params['dim']: highest feature index
                         params['-d']:  degree polynomial
                         params['-g']:  gamma in rbf
                         params['-s']:  parameter s in sigmoid/poly kernel
                         params['-r']:  parameter c in sigmoid/poly kernel
        '''
        self.sv=sv
        self.numsv=len(sv)
        self.ai=ai
        self.b=b
        if kernel == 0:
           self.kernel=self._klin
        elif kernel == 1:
           self.kernel=self._kpoly
        elif kernel == 2:
           self.kernel=self._krbf
        elif kernel == 3:
           self.kernel=self._ksig
        else:
           self.kernel=_klin
        self.dim=params['dim']
        self.d=params['-d']
        self.g=params['-g']
        self.s=params['-s']
        self.t=params['-r']

    def _klin(self,v1,v2):
        ''' fake kernel '''
        return self._krbf(v1,v2)

    def _kpoly(self,v1,v2):
        ''' fake kernel '''
        return self._krbf(v1,v2)

    def _krbf(self,v1,v2):
        ''' _krbf(self,v1,v2)'''
        x=v1-v2
        return numpy.exp(-self.g*(numpy.dot(x,x)))

    def predict(self,x):
        kval=0.0
        for i in range(self.numsv):
            kval+=self.ai[i]*self.kernel(self.sv[i],x)
        return kval-self.b

    def svmSave(self,fname):
        ''' svmSave(self,fname) '''
        f=open(fname,'w')
        cPickle.dump(self,f)
        f.close()


def getsvmPickle(fname):
    ''' getsvmPickle(fname) '''
    f=open(fname)
    svm=cPickle.load(f)
    return svm

def unpacksvmVec(x,vecDim):
    ''' unpacksvmVec(x,vecDim) '''
    nv=numpy.zeros(vecDim,float)
    v=x.split('#')[0] # if contains comment
    v=v.split()
    first=float(v[0])
    for e in v[1:]:
        (idx,val)=e.split(':')
        idx,val=int(idx),float(val)
        nv[idx-1]=val
    return first,nv

def getSVMLight(fname):
    ''' getSVMMod(fname) '''
    lines=open(fname).readlines()
    # set kernel type
    i=0
    while lines[i].find('kernel type')<0:
        i+=1
    kernel=int(lines[i].split()[0])
    # parameters
    params={}
    while lines[i].find('threshold b')<0: # loop until threshold
        if lines[i].find('kernel parameter') >=0:
            v=lines[i].split()
            if v[-1] != '-u':
               params[v[-1]]=float(v[0])
        elif lines[i].find('highest feature index')>=0:
            params['dim']=int(lines[i].split()[0])
        i+=1
     
    b=float(lines[i].split()[0])
    # print "threshold found ",i
    i+=1
    ai=[]
    sv=[]
    while i < len(lines):
        a,v=unpacksvmVec(lines[i],params['dim'])
        ai.append(a)
        sv.append(v)
        i+=1
    return SVMLike(sv,ai,b,kernel,params)

#-----------MAIN--------------#


if __name__=='__main__':
        #testing
        if len(sys.argv) <4 :
                sys.stderr.write("usage: "+sys.argv[0]+''' fastafile profile type(A F P)
        A===>Animali(Metazoa)
        F===>Funghi(Fungi)
        P===>Piante(Viridiplantae)\n''')
                sys.stderr.write("="*30+'\n'+'!!! Warning !!!: the profile must be 20 letters\n')
                sys.stderr.write("="*30+'\n')
                sys.exit()
        seq=read1Fasta(sys.argv[1])
        aaOrder,prof,other_features=readProf(sys.argv[2],aaOrder='VLIMFWYGAPSTCHRKQEND')
        prof=numpy.array(prof)
        typePred=sys.argv[3]
        loc, sc = bacello(seq,prof,typePred)
        print loc, sc

