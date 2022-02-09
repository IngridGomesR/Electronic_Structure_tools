import numpy as np
from kappa2 import coords

def gerafile(namefile, atomos, vetores):
    with open(namefile,'w') as f:
        f.write(str(np.shape(vetores)[0])+"\n"+"\n")
        for i in range(len(atomos)):
            f.write(str(atomos[i][0])+"\t"+str(vetores[i,0])+"\t"+str(vetores[i,1])+"\t"+str(vetores[i,2])+"\n")
    return 
    
def Rota(namefile, thetax, thetay, thetaz, newfile):
    atomos, vetores = coords(namefile)
    Rx = np.array([[1,0,0],[0,np.cos(thetax),-np.sin(thetax)],[0, np.sin(thetax), np.cos(thetax)]])
    Ry = np.array([[np.cos(thetay),0,np.sin(thetay)],[0,1,0],[-np.sin(thetay), 0, np.cos(thetay)]])
    Rz = np.array([[np.cos(thetaz), -np.sin(thetaz),0],[np.sin(thetaz), np.cos(thetaz),0],[0,0,1]])
    L =[]
    for vetor in vetores:
        rotax = np.dot(Rx,vetor)
        rotay = np.dot(Ry,rotax)
        rotaz = np.dot(Rz,rotay)
        L.append(rotaz)
    rota = np.array(L)
    gerafile(newfile, atomos, rota) 
    return rota

#rota = Rota('BenzenoTrans.xyz', np.pi/2, 0, np.pi/4, 'BenzenoRota.xyz')
#print(rota)



    
        