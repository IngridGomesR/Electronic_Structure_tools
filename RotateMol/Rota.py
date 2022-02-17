import numpy as np

def coords(fileDim):
    vetores = np.array([0,0,0])
    atomos  = np.array([0])
    with open(fileDim, 'r') as f:
        for line in f:
            line = line.split()
            if len(line) < 4:
                pass
            else:
                vetor   = np.array([float(line[1]),float(line[2]),float(line[3])])
                vetores = np.vstack((vetores,vetor))
                atomo   = np.array([line[0]])
                atomos  = np.vstack((atomos,atomo))
    atomos  = np.delete(atomos,0,0)
    vetores = np.delete(vetores,0,0)
    return atomos, vetores

def pega_geom(file):
    vetores = np.array([0,0,0])
    atomos  = np.array([0])
    arquivo = open(file, 'r')
    lista = arquivo.readlines()
    i=1
    for linha in lista:
        i += 1
        if 'Standard orientation:' in linha:
            j=i+3
        if "NAtoms=" in linha:
            aux=linha.split()
            num_a=int(aux[1])
    
    geom=lista[j:j+num_a]
    for linha in geom:
        aux=linha.split()
        vetor   = np.array([float(aux[3]),float(aux[4]),float(aux[5])])
        vetores = np.vstack((vetores,vetor))
        atomo   = np.array([aux[1]])
        atomos  = np.vstack((atomos,atomo))
    atomos  = np.delete(atomos,0,0)
    vetores = np.delete(vetores,0,0)    
    return atomos, vetores

def gerafile(namefile, atomos, vetores):
    with open(namefile,'w') as f:
        f.write(str(np.shape(vetores)[0])+"\n"+"\n")
        for i in range(len(atomos)):
            f.write(str(atomos[i][0])+"\t"+str(vetores[i,0])+"\t"+str(vetores[i,1])+"\t"+str(vetores[i,2])+"\n")
    return 
    
def Rota(namefile, thetax, thetay, thetaz, newfile):
    atomos, vetores = coords(namefile) #pega_geom(namefile) #
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

#rota = Rota('3bRota180.xyz', np.pi, 0, 0, '3bRota360.xyz')

def dif(fileS1, fileRota):
    atomosS1, vetoresS1 = pega_geom(fileS1)
    atomosRota, vetoresRota = coords(fileRota) 
    dif = abs(vetoresS1 - vetoresRota)
    print(dif)


#dif('3bS1.log', '3bRota360.xyz')

    
        