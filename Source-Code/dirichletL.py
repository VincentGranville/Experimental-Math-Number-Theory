# DirichletL.py. Generate orbits of various Dirichlet-L and related functions
# By Vincent Granville, https://www.MLTechniques.com/resources/

import math
import random
from primePy import primes

nterms=2000 # increase to 10000 for sig = 0.5
method='Eta'
sig=0.9
Dirichlet=False 
x=1  # must have 0 < x <= 1; default is x=1
beta=0.5 # beta > 0.5 magnifies the hole of the orbit 

random.seed(1)
primeSign={}
start=2
if method=='Dirichlet4':
  start=3
idx=0
for k in range(start,nterms):
  if primes.check(k):
    idx=idx+1
    p=k
    xpow=x**(1/p)
    if method=='Beurling' and p==3:
      p=2+math.log(3)
    primeSign[p]=xpow
    if method=='Dirichlet4' and k%4==3:
      primeSign[p]=-xpow
    elif method=='Alternating' and idx%2==1:
      primeSign[p]=-xpow
    elif method=='Random' and random.random()>0.5:
      primeSign[p]=-xpow
    elif method=='Eta':
      Dirichlet=True

signHash={}
evenHash={}  
signHash[1]=1
evenHash[1]=0           # largest power of 2 dividing k
for p in primeSign:
  if p*math.pi %1 < 0.05: 
    print(p,"/",nterms) # show progress (where we are in the loop) 
  oldSignHash={}
  for k in signHash:
    oldSignHash[k]=signHash[k]
  for k in oldSignHash:
    pp=1
    power=0
    localProduct=oldSignHash[k]  
    while k*p*pp<nterms:
      pp=p*pp
      power=power+1 
      new_k=k*pp
      localProduct=localProduct*primeSign[p]
      signHash[new_k]=localProduct  
      if p==2:
        evenHash[new_k]=power
      else: 
        evenHash[new_k]=evenHash[k]

for k in sorted(evenHash):
  if Dirichlet and evenHash[k]>0: 
    signHash[k]=-signHash[k]

sumL=0
minL= 2*nterms
maxL=-2*nterms
argMin=-1
argMax=-1
denum={}
tlog={}
for k in sorted(signHash):
  denum[k]=signHash[k]/k**sig    
  tlog[k]=math.log(k)
  sumL=sumL+signHash[k]
  if sumL<minL:
    minL=sumL
    argMin=k
  if sumL>maxL:
    maxL=sumL
    argMax=k
denum[2]=signHash[2]/(1/beta)**sig 

def G(tau,sig,nterms):
  fetax=0
  fetay=0
  for j in sorted(signHash):
    fetax=fetax+math.cos(tau*tlog[j])*denum[j]
    fetay=fetay+math.sin(tau*tlog[j])*denum[j]
  return [fetax,fetay]

minT=0.0
maxT=2000.0 
increment=0.05

OUT  = open("dirichletL.txt", "w")
t=minT
loop=0
while t <maxT:
  if loop%100==0:
    print("t= %5.2f / %d" % (t,maxT))
  loop=loop+1
  (etax,etay)=G(t,sig,nterms)
  line=str(t)+"\t"+str(etax)+"\t"+str(etay)+"\n"
  OUT.write(line) 
  t=t+increment
OUT.close()

print("\n")
print(argMin,"-->",minL)
print(argMax,"-->",maxL)
