import numpy as np
import matplotlib.pyplot as plt
import random
import sys

sistema = int(float(sys.argv[1]))
n_grupos = float(sys.argv[2])
identidade = sys.argv[3]

nproc=30
mem= '28GB'
funcional = 'wb97xd/6-31G(d,p)'

lsize=1.41
w= np.sqrt(3)*lsize
h=2*lsize
x0=0
y0=0
z0=0

linhas=2*sistema #numeros de linhas
pontos_int= linhas + 1 #numero de pontos na primeira linha
pontos=pontos_int + 1
n = sistema + 1

coords = []

filename = 's_'+str(sistema)+'_g_'+str(n_grupos)+'_id_'+identidade+'_sistemas.xyz'

with open(filename, 'w') as f:
    f.write('%nproc='+str(nproc)+'\n')
    f.write('%mem='+str(mem)+'\n')
    f.write('%chk='+filename[:-3]+'chk\n')
    f.write('# '+funcional+' OPT \n')
    f.write('\n')
    f.write('title\n')
    f.write('\n')
    f.write('0 2\n')

hidros = []

#LINHAS DE BAIXO
with open(filename, 'a') as f:
    for i in range(1,n): #faz as linhas
        for j in range(1,pontos): #faz os pontos da linha
            k = j + 1
            xC = x0 + (j-1)*w*0.5
            yC = y0 - (k-int(k/2)*2)*h*0.25
            zC = z0
            coords.append([xC,yC,zC])
            f.write(str('C')+"\t"+str("{:.6f}".format(xC))+"\t"+str("{:.6f}".format(yC))+"\t"+str("{:.6f}".format(zC))+"\n")
            if i==1 and j==1:
                theta = np.arctan(0.25*h/(0.5*w))
                xH = xC - 0.96*np.cos(theta)
                yH = yC - 0.96*np.sin(theta)
                zH = zC
                hidros.append([xH,yH,zH])
                
            elif i==1 and j%2==0:
                xH = xC 
                yH = yC - 0.96
                zH = zC
                hidros.append([xH,yH,zH])
                
            elif i==1 and j==pontos_int:
                theta = np.arctan(0.25*h/(0.5*w))
                xH = xC + 0.96*np.cos(theta)
                yH = yC - 0.96*np.sin(theta)
                zH = zC
                hidros.append([xH,yH,zH])
                
            elif i!=1 and j==p:
                theta = np.arctan(0.25*h/(0.5*w))
                xH = xC + 0.96*np.cos(theta)
                yH = yC - 0.96*np.sin(theta)
                zH = zC
                hidros.append([xH,yH,zH])
                
            elif i!=1 and j==1:
                theta = np.arctan(0.25*h/(0.5*w))
                xH = xC - 0.96*np.cos(theta)
                yH = yC - 0.96*np.sin(theta)
                zH = zC
                hidros.append([xH,yH,zH])
                
        x0 = x0 - 0.5*w
        y0 = y0 + 0.75*h
        pontos = pontos + 2
        p = pontos - 1

#LINHAS DE CIMA
x0 = x0 + w*0.5
y0 = y0 - 0.75*h
z0 = 0
pontos = pontos - 2

y0 = y0 + lsize
p = pontos - 1

with open(filename, 'a') as f:
    for i in range(1,n):
        for j in range(1,pontos): 
            k = j + 1
            xC = x0 + (j-1)*w*0.5
            yC = y0 + (k-int(k/2)*2)*h*0.25
            zC = z0
            coords.append([xC,yC,zC])
            f.write(str('C')+"\t"+str("{:.6f}".format(xC))+"\t"+str("{:.6f}".format(yC))+"\t"+str("{:.6f}".format(zC))+"\n")
            a = n-1
            if i==a and j==1:
                theta = np.arctan(0.25*h/(0.5*w))
                xH = xC - 0.96*np.cos(theta)
                yH = yC + 0.96*np.sin(theta)
                zH = zC
                hidros.append([xH,yH,zH])
                
            elif i==a and j%2==0:
                xH = xC 
                yH = yC + 0.96
                zH = zC
                hidros.append([xH,yH,zH])
                
            elif i==a and j==pontos_int:
                theta = np.arctan(0.25*h/(0.5*w))
                xH = xC + 0.96*np.cos(theta)
                yH = yC + 0.96*np.sin(theta)
                zH = zC
                hidros.append([xH,yH,zH])
                
            elif i!=a and j==p:
                theta = np.arctan(0.25*h/(0.5*w))
                xH = xC + 0.96*np.cos(theta)
                yH = yC + 0.96*np.sin(theta)
                zH = zC
                hidros.append([xH,yH,zH])
                
            elif i!=a and j==1:
                theta = np.arctan(0.25*h/(0.5*w))
                xH = xC - 0.96*np.cos(theta)
                yH = yC + 0.96*np.sin(theta)
                zH = zC
                hidros.append([xH,yH,zH])
                
        x0 = x0 + 0.5*w
        y0 = y0 + 0.75*h
        pontos = pontos - 2
        p = pontos - 1

#grupo hidroxi
def bota_hidroxi(coords):
    xC = coords[0]
    yC = coords[1]
    zC = coords[2]
    with open(filename, 'a') as f:
        xO = xC
        yO = yC
        zO = zC + random.choice([-1,1])*1.47
        theta = 17.9*np.pi/180
        xH = xO + 0.98*np.cos(theta)
        yH = yO
        zH = zO + 0.98*np.sin(theta)
        f.write(str('H')+"\t"+str("{:.6f}".format(xH))+"\t"+str("{:.6f}".format(yH))+"\t"+str("{:.6f}".format(zH))+"\n")
        f.write(str('O')+"\t"+str("{:.6f}".format(xO))+"\t"+str("{:.6f}".format(yO))+"\t"+str("{:.6f}".format(zO))+"\n")


#grupo epoxy 
def bota_epoxy(coords):
    xC = coords[0]
    yC = coords[1]
    zC = coords[2]
    with open(filename, 'a') as f:
        yO = yC 
        xO = xC 
        zO = zC + random.choice([-1,1])*np.sqrt((1.44**2)-((h*0.125)**2))
        f.write(str('O')+"\t"+str("{:.6f}".format(xO))+"\t"+str("{:.6f}".format(yO))+"\t"+str("{:.6f}".format(zO))+"\n")

#grupo acido carboxilico
def bota_acido(coords):
    xC = coords[0]
    yC = coords[1]
    zC = coords[2]
    with open(filename, 'a') as f:
        theta = np.arctan(0.25*h/(0.5*w))
        xC1 = xC 
        yC1 = yC 
        nz = random.choice([-1,1])
        zC1 = zC + nz*lsize
        xO1 = xC1
        yO1 = yC1 - 1.47*np.cos(theta)
        zO1 = zC1 + nz*1.47*np.sin(theta)
        xO2 = xC1 
        yO2 = yC1 + 1.47*np.cos(theta)
        zO2 = zC1 + nz*1.47*np.sin(theta)
        xH = xO2 
        yH = yO2 
        zH = zO2 + nz*0.98
        f.write(str('C')+"\t"+str("{:.6f}".format(xC1))+"\t"+str("{:.6f}".format(yC1))+"\t"+str("{:.6f}".format(zC1))+"\n")
        f.write(str('O')+"\t"+str("{:.6f}".format(xO1))+"\t"+str("{:.6f}".format(yO1))+"\t"+str("{:.6f}".format(zO1))+"\n")
        f.write(str('O')+"\t"+str("{:.6f}".format(xO2))+"\t"+str("{:.6f}".format(yO2))+"\t"+str("{:.6f}".format(zO2))+"\n")
        f.write(str('H')+"\t"+str("{:.6f}".format(xH))+"\t"+str("{:.6f}".format(yH))+"\t"+str("{:.6f}".format(zH))+"\n")

def bota_hidros(hidros):
    with open(filename, 'a') as f:
        for hidro in hidros:
            xO = hidro[0] 
            yO = hidro[1]
            zO = hidro[2]
            f.write(str('H')+"\t"+str("{:.6f}".format(xO))+"\t"+str("{:.6f}".format(yO))+"\t"+str("{:.6f}".format(zO))+"\n")    


meiotas = []
xs = np.array([a[0] for a in coords])
ys = np.array([a[1] for a in coords])
zs = np.array([a[2] for a in coords])
for n in range(len(coords)):
    r = np.sqrt((xs - xs[n])**2 + (ys - ys[n])**2)
    for m in range(len(r)):
        if r[m] < 1.5 and r[m] > 0:
            x1 = max(xs[m],xs[n])
            x2 = min(xs[m],xs[n])
            y1 = max(ys[m],ys[n])
            y2 = min(ys[m],ys[n])
            X = np.round((x1 + x2)/2,2)
            Y = np.round((y1 + y2)/2,2)
            if [X,Y,0] not in meiotas:
                meiotas.append([X,Y,0])
            

funcao =  [bota_hidroxi,bota_acido,bota_epoxy]
total = 0

print(len(hidros))

n_dopant = int(n_grupos*(len(coords)+len(meiotas)+len(hidros)))

while total < n_dopant:
    func = random.choice(funcao)
    if func == bota_epoxy:
        try:
            coord = random.choice(meiotas)  
            del meiotas[meiotas.index(coord)]
            func(coord)
            total += 1
        except:
            pass
    elif func == bota_acido:
        try:
            coord = random.choice(hidros)  
            del hidros[hidros.index(coord)]
            func(coord)
            total += 1
        except:
            pass
    else:
        try:
            coord = random.choice(coords)
            del coords[coords.index(coord)]
            func(coord)
            total += 1
        except:
            pass    

print(len(hidros))
bota_hidros(hidros)



with open(filename, 'a') as f:
    f.write('\n')
    f.write('--Link1--\n')
    f.write('%nproc='+str(nproc)+'\n')
    f.write('%mem='+str(mem)+'\n')
    f.write('%chk='+filename[:-3]+'chk\n')
    f.write('# '+funcional+' TD=(Nstates=1) geom=allcheck \n')
    f.write('\n')
    
   

   











