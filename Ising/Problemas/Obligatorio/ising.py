# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:19:12 2024

@author: usuario
"""

def deltaE_T(i,j):
    ip=i
    im=i
    jp=j
    jm=j
    
    """
    Comprobamos la i y asignamos valores a ip (i+1) y a im (i-1)
    """
    if i==0:
        ip=i+1
        im=N-1
    elif i==N-1:
        ip=0
        im=j-1
    else:
        ip=i+1
        im=i-1
    
    """
    Comprobamos la j y asignamos valores a jp (j+1) y a jm (j-1)
    """
    if j==0:
        jp=j+1
        jm=N-1
    elif j==N-1:
        jp=0
        jm=j-1
    else:
        jp=j+1
        jm=j-1
        
    deentret=np.exp((-2*s[tempn,i,j]*(s[tempn,im,j]+s[tempn,ip,j]+s[tempn,i,jm]+s[tempn,i,jp]))/t)
    return deentret

#########################################################################################
#########################################################################################
#########################################################################################

import numpy as np

"""
Especificamos tamaño y generamos temperaturas entre 0 y 5 cada 0.5
"""
N=16
aT=np.arange(0,5,0.5)
aT[0]=0.1
print(aT)

"""
Abrimos ficheros para todo
"""
temp0=open('temp0.txt','w')
temp1=open('temp0.5.txt','w')
temp2=open('temp1.txt','w')
temp3=open('temp1.5.txt','w')
temp4=open('temp2.txt','w')
temp5=open('temp2.5.txt','w')
temp6=open('temp3.txt','w')
temp7=open('temp3.5.txt','w')
temp8=open('temp4.txt','w')
temp9=open('temp4.5.txt','w')
temp10=open('temp5.txt','w')
magne0=open('magnet0.txt','w')
magne1=open('magnet0.5.txt','w')
magne2=open('magnet1.txt','w')
magne3=open('magnet1.5.txt','w')
magne4=open('magnet2.txt','w')
magne5=open('magnet2.5.txt','w')
magne6=open('magnet3.txt','w')
magne7=open('magnet3.5.txt','w')
magne8=open('magnet4.txt','w')
magne9=open('magnet4.5.txt','w')
magne10=open('magnet5.txt','w')
magnem=open('magnetmedia.txt','w')

"""
Bucle para generar la configuración aleatoria de espines para las distintas temperaturas
"""
s=np.zeros((10,N,N))
smm=np.ones((10,N,N))           #Configuración ordenada de espín 1
prob=np.random.random((10,N,N))

for i in range(10):
    for j in range(N):
        for k in range(N):
            if prob[i,j,k] < 0.5:
                s[i,j,k]=-1
            else:
                s[i,j,k]=1
    #print(s[i,:,:])

tempn=0         #Valor que especificará en que posición de aT estamos

"""
Generamos un array con 2 valores, del cual sacaremos después el mínimo entre 1 y 
exp(-delta[E]/T)
"""
cond=np.ones(2)
condmm=np.ones(2)

"""
Creamos el bucle grande y empezamos los cálculos por cada temperatura
"""
vez=0
estacionario=np.array([False,False,False,False,False,False,False,False,False,False,False],dtype=bool)

magmm=(np.sum(smm[tempn,:,:]))/N**2
magnet=np.zeros(1e4)
for t in aT:
    #rand=np.random.randint(0,N,2)
    #print(rand)
    magnet[:]=0
    nmm=0
    desvest=0
    for i in np.arange(0,1e5):
        """
        Generamos la posición aleatoria que podrá cambiar el espin y dseta
        """
        rand=np.random.randint(0,N,2)
        dseta=np.random.random(1)
        
        """
        Comprobamos si cambia calculando el mínimo, comparamos con dseta y cambiamos si
        se cumple
        """
        cond[1]=deltaE_T(rand[0],rand[1])

        if np.min(cond)>dseta:
            s[tempn,rand[0],rand[1]]*=-1
            
        """
        Repetimos para la configuración ordenada con el objetivo de la magnetización media
        """
        rand=np.random.randint(0,N,2)
        dseta=np.random.random(1)
        condmm[1]=deltaE_T(rand[0],rand[1])
        
        if np.min(condmm)>dseta:
            smm[tempn,rand[0],rand[1]]*=-1
                    
        """
        Guardamos la magnetización cada pmc
        """
        mag=(np.sum(s[tempn,:,:]))/N**2
        
        magmm=(np.sum(smm[tempn,:,:]))/N**2
        
        if i%256==0:
            if tempn==0:
                magne0.write(f'{i/256}\t{mag}\n')
            elif tempn==1:
                magne1.write(f'{i/256}\t{mag}\n')
            elif tempn==2:
                magne2.write(f'{i/256}\t{mag}\n')
            elif tempn==3:
                magne3.write(f'{i/256}\t{mag}\n')
            elif tempn==4:
                magne4.write(f'{i/256}\t{mag}\n')
            elif tempn==5:
                magne5.write(f'{i/256}\t{mag}\n')
            elif tempn==6:
                magne6.write(f'{i/256}\t{mag}\n')
            elif tempn==7:
                magne7.write(f'{i/256}\t{mag}\n')
            elif tempn==8:
                magne8.write(f'{i/256}\t{mag}\n')
            elif tempn==9:
                magne9.write(f'{i/256}\t{mag}\n')
            elif tempn==10:
                magne10.write(f'{i/256}\t{mag}\n')
                
        """
        Guardamos los datos cada NxN pasos
        """
        if i%(N**2):
            if i>=9e4:
                nmm+=1
                
            if tempn==0:
                for j in range(N):
                    for k in range(N):
                        temp0.write(f'{s[tempn,j,k]}\t')
                temp0.write("\n")
            elif tempn==1:
                for j in range(N):
                    for k in range(N):
                        temp1.write(f'{s[tempn,j,k]}\t')
                temp1.write("\n")
            elif tempn==2:
                for j in range(N):
                    for k in range(N):
                        temp2.write(f'{s[tempn,j,k]}\t')
                temp2.write("\n")
            elif tempn==3:
                for j in range(N):
                    for k in range(N):
                        temp3.write(f'{s[tempn,j,k]}\t')
                temp3.write("\n")
            elif tempn==4:
                for j in range(N):
                    for k in range(N):
                        temp4.write(f'{s[tempn,j,k]}\t')
                temp4.write("\n")
            elif tempn==5:
                for j in range(N):
                    for k in range(N):
                        temp5.write(f'{s[tempn,j,k]}\t')
                temp5.write("\n")
            elif tempn==6:
                for j in range(N):
                    for k in range(N):
                        temp6.write(f'{s[tempn,j,k]}\t')
                temp6.write("\n")
            elif tempn==7:
                for j in range(N):
                    for k in range(N):
                        temp7.write(f'{s[tempn,j,k]}\t')
                temp7.write("\n")
            elif tempn==8:
                for j in range(N):
                    for k in range(N):
                        temp8.write(f'{s[tempn,j,k]}\t')
                temp8.write("\n")
            elif tempn==9:
                for j in range(N):
                    for k in range(N):
                        temp9.write(f'{s[tempn,j,k]}\t')
                temp9.write("\n")
            elif tempn==10:
                for j in range(N):
                    for k in range(N):
                        temp10.write(f'{s[tempn,j,k]}\t')
                temp10.write("\n")
            #print(f'pmc {pmc}')      #Vemos cuántos pmc tenemos por seguridad
            
            if nmm>=9e4:
                magnet[nmm]=magmm

    magnetmean=np.mean(magnet)
    desvest=np.std(magnet)
    magnem.write(f'{t}\t{magnetmean}\t{desvest}\n')
    tempn+=1    #Nos movemos a la siguiente temperatura
    
temp0.close()
temp1.close()
temp2.close()
temp3.close()
temp4.close()
temp5.close()
temp6.close()
temp7.close()
temp8.close()
temp9.close()
temp10.close()
magne0.close()
magne1.close()
magne2.close()
magne3.close()
magne4.close()
magne5.close()
magne6.close()
magne7.close()
magne8.close()
magne9.close()
magne10.close()
magnem.close()