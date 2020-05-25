#------------PIERO FARISELLI--------------#
import numpy

class SVMLike:
    ''' This class implement the prediction phase of a SVM
    based on the model of svm light.
    Presently, there are several limitation including:
    - fixed dimension of the input vectors ()
    - !!! only rbf kernel implemente up to now !!!!'''

    def __init__(self,sv,ai,b,kernel,params):
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
