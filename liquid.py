import numpy as np 
from random import random 
import sys  

def get_minImageDist(dx,dy,dz,boxLen): 
    if(dx>0.5*boxLen):
        dx-=boxLen 
    elif(dx<-0.5*boxLen):
        dx+=boxLen 

    if(dy>0.5*boxLen):
        dy-=boxLen 
    elif(dy<-0.5*boxLen):
        dy+=boxLen 

    if(dz>0.5*boxLen):
        dz-=boxLen 
    elif(dz<-0.5*boxLen):
        dz+=boxLen 
    
    return(dx,dy,dz)

def get_LiquidCoord(boxLen,tol,max_try):
    
    x=[]
    y=[]
    z=[] 


    nInserted=0
    tol2=tol**2 

    iTry=0 
    while(iTry<max_try): 
        iTry+=1 
        xi=np.random.random()*boxLen 
        yi=np.random.random()*boxLen 
        zi=np.random.random()*boxLen 

        overlap=False 
        for i in range(nInserted):
            dx=xi-x[i]
            dy=yi-y[i]
            dz=zi-z[i]
            dx,dy,dz=get_minImageDist(dx,dy,dz,boxLen)
            rsq=dx**2 + dy**2 + dz**2 
            if(rsq<tol2):
                overlap=True 
                break 

        if(overlap == False):
            x.append(xi)
            y.append(yi)
            z.append(zi)
            nInserted+=1 
            print('Inserted: '+str(nInserted)+' atom/s;'+' nTrys: '+str(iTry))
            iTry=0 


    return(x,y,z) 

def write_coord(coord_file,x,y,z,symbol,nAtoms): 
    fw = open(coord_file, "w")
    fw.write('{0}\n'.format(nAtoms))
    fw.write('Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" Properties=species:S:1:pos:R:3 Zrme=0.0')
    fw.write('\n')
    for i in range(nAtoms): 
        fw.write('{0}   {1:.7f}   {2:.7f}   {3:.7f}\n'.format(symbol,x[i],y[i],z[i]))

max_try=1e6
boxLen=10.0 
tol=2.5 # minimum distance of atoms
symbol='Ni'
coord_file='Ni_liquid.xyz' 

x,y,z=get_LiquidCoord(boxLen,tol,max_try)
nAtoms=len(x)
X=np.asarray(x)/boxLen
Y=np.asarray(y)/boxLen
Z=np.asarray(z)/boxLen
#Y=y/boxLen
#Z=z/boxLen
write_coord(coord_file,x,y,z,symbol,nAtoms)
