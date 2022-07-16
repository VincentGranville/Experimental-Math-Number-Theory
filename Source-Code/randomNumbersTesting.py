# Test rabdomness of binary sequences via the law of the iterated lograrithm
# By Vincent Granville, www.MLTechniques.com

import math
import random
import numpy as np
from primePy import primes

#--
def createRandomDigits(method,seed): 
  primeSign={}
  idx=0
  if method=='SQRT':
    y=seed[0]
    z=seed[1]
  elif method=='Python':
    random.seed(seed)
  else:
    x=seed
  start=2
  if method=='Dirichlet4':
    start=3
  for k in range(start,nterms):  
    if primes.check(k):
      primeSign[k]=1
      if method=='SQRT':
        if z<2*y:   
          y=4*y-2*z
          z=2*z+3
        else:     
          y=4*y
          z=2*z-1
          primeSign[k]=-1
      elif method=='Dirichlet4':
        if k%4==seed:
          primeSign[k]=-1
      elif method=='CounterExample':
        idx=idx+1
        if idx%2==seed:
          primeSign[k]=-1
      elif method=='Python':
        x=random.random()
      elif method=='Logistic':
        x=4*x*(1-x)
      elif method=='Base3':
        x=3*x-int(3*x)
      if method in ('Python','Logistic','Base3') and x>0.5:
        primeSign[k]=-1
  return(primeSign)

#--

def createSignHash2():
  signHash={}
  signHash[1]=1
  for p in primeSign:
    oldSignHash={}
    for k in signHash:
      oldSignHash[k]=signHash[k]
    for k in oldSignHash:
      pp=1
      power=0
      localProduct=oldSignHash[k]  ### signHash[k]
      while k*p*pp<nterms:
        pp=p*pp
        power=power+1
        new_k=k*pp
        localProduct=localProduct*primeSign[p]
        signHash[new_k]=localProduct  
  return(signHash)

#--

def createSignHash():
  # same as createSignHash() but for square-free integers only
  signHash={}
  signHash[1]=1
  for p in primeSign:
    oldSignHash={}
    for k in signHash:
      oldSignHash[k]=signHash[k]
    for k in oldSignHash:
      if k*p<nterms:
        new_k=k*p 
        signHash[new_k]=oldSignHash[k]*primeSign[p]  
  return(signHash)
 
#--

def testRandomness(category):
  signHash=createSignHash2()
  count=0
  count1=0
  sumL=0
  minL= 2*nterms
  maxL=-2*nterms
  argMin=-1
  argMax=-1
  for k in sorted(signHash):
    selected=False
    if category=='Prime' and k in primeSign:
      selected=True
    elif category=='nonPrime' and k not in primeSign:
      selected=True
    elif category=='All':
      selected=True
    if selected==True:
      if signHash[k]==1:
        count1=count1+1  
      count=count+1
      sumL=sumL+signHash[k]
      if sumL<minL:
        minL=sumL
        argMin=count
      if sumL>maxL:
        maxL=sumL
        argMax=count
  return(minL,argMin,maxL,argMax,count,count1)

#--

# Main Part. Requirements:
#   0 < seed < 1 for 'Base3' and 'Logistic'; rational numbers not random
#   seed=(y,z) with z>y, z!=2y, y!=2x and x,y>0 are integers for 'SQRT'
#   swapping -1/+1 for seed=(90,91) in 'SQRT' does well, the original does not

nterms=20000
seedMethod={}
seedMethod['Python']=(0,1,2,4,100,200,500)
seedMethod['Logistic']=(0.181517,0.72)
seedMethod['Base3']=(0.181517,0.72)
seedMethod['SQRT']=((2,5),(90,91))
seedMethod['Dirichlet4']=(1,3)
seedMethod['CounterExample']=(1,0)

OUT=open("prgTest.txt", "w")
for method in seedMethod:
  for seed in seedMethod[method]:

    primeSign=createRandomDigits(method,seed)
    
    [minL,argMin,maxL,argMax,count,count1]=testRandomness('Prime')

    print(method,'Prime',minL,maxL,"|",argMin,argMax,"|",count1,count,seed)
    line=method+'\tPrime\t'+str(minL)+'\t'+str(maxL)+'\t'+str(argMin)+'\t'
    line=line+str(argMax)+'\t'+str(count1)+'\t'+str(count)+'\t'+str(seed)+'\n'
    OUT.write(line)
    
    [minL,argMin,maxL,argMax,count,count1]=testRandomness('nonPrime')

    print(method,'nonPrime',minL,maxL,"|",argMin,argMax,"|",count1,count,seed)
    line=method+'\tnonPrime\t'+str(minL)+'\t'+str(maxL)+'\t'+str(argMin)+'\t'
    line=line+str(argMax)+'\t'+str(count1)+'\t'+str(count)+'\t'+str(seed)+'\n'
    OUT.write(line)


    [minL,argMin,maxL,argMax,count,count1]=testRandomness('All')

    print(method,'All',minL,maxL,"|",argMin,argMax,"|",count1,count,seed)
    line=method+'\tAll\t'+str(minL)+'\t'+str(maxL)+'\t'+str(argMin)+'\t'
    line=line+str(argMax)+'\t'+str(count1)+'\t'+str(count)+'\t'+str(seed)+'\n'
    OUT.write(line)

OUT.close()