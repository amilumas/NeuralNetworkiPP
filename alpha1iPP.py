import numpy as np
import makesemicrystallineLmpsfEMC as msc
import math


def linear2parabolaXZConnect(la, lb, pa1, ph1, pk1, pa2, ph2, pk2, start, end, mid, spacing):
    npoints = 1000
    xpos1 = np.linspace(start[0], mid[0], npoints)
    xpos2 = np.linspace(mid[0], end[0], npoints)
    ypos1 = la*xpos1 + lb*np.ones(npoints)
    ypos2 = la*xpos2 + lb*np.ones(npoints)
    zpos1 = pa1*(xpos1 - ph1)**2 + pk1
    zpos2 = pa2*(xpos2 - ph2)**2 + pk2
    #print("xpos1", xpos1)
    #print("zpos1", zpos1)
    dists1 = np.zeros(npoints)
    dists2 = np.zeros(npoints)
    conx = []
    cony = []
    conz = []
    sidex = []
    sidey = []
    sidez = []
    dists1[0] = 0
    pad = 1.5
    extra = 0
    nspacing = spacing
    for i in range(npoints-1):
        dists1[i+1] = dists1[i] + ((xpos1[i] - xpos1[i+1])**2 + (ypos1[i] - ypos1[i+1])**2 + (zpos1[i] - zpos1[i+1])**2)**0.5
        if dists1[i] >= nspacing:
            conx.append(xpos1[i+1])
            cony.append(ypos1[i+1])
            conz.append(zpos1[i+1])
            nspacing = nspacing + spacing
    dists2[0] = dists1[-1]
    N1 = len(conx)
    for i in range(N1):
        if i%2 == 0:
            sy = cony[i] - pad
            slope = -1*(pa1*2*(conx[i] - ph1))
            inter = conz[i] - slope*conx[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            #print("i", i,"normvec", normvec)
            sx = conx[i] + np.sign(normvec[0])*(spacing-pad + extra)*(normvec[0])
            sz = conz[i] + np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)

    for i in range(npoints-1):
        dists2[i+1] = dists2[i] + ((xpos2[i] - xpos2[i+1])**2 + (ypos2[i] - ypos2[i+1])**2 + (zpos2[i] - zpos2[i+1])**2)**0.5
        distend = ((xpos2[i+1] - end[0])**2 + (ypos2[i+1] - end[1])**2 + (zpos2[i+1] - end[2])**2)**2
        if dists2[i] >= nspacing and distend > spacing - 0.1:
            conx.append(xpos2[i+1])
            cony.append(ypos2[i+1])
            conz.append(zpos2[i+1])
            nspacing = nspacing + spacing
    N = len(conx)
    for i in range(N1, N):
        if i%2 == 0:
            sy = cony[i] + pad
            slope = -1*(pa2*2*(conx[i] - ph2))
            inter = conz[i] - slope*conx[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            sx = conx[i] -  np.sign(normvec[0])*(spacing-pad + extra)*normvec[0]
            sz = conz[i] -  np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)
    return conx, cony, conz, sidex, sidey, sidez

def linear2parabolaXZConnect2(la, lb, pa1, ph1, pk1, pa2, ph2, pk2, start, end, mid, spacing):
    npoints = 1000
    xpos1 = np.linspace(start[0], mid[0], npoints)
    xpos2 = np.linspace(mid[0], end[0], npoints)
    ypos1 = la*xpos1 + lb*np.ones(npoints)
    ypos2 = la*xpos2 + lb*np.ones(npoints)
    zpos1 = pa1*(xpos1 - ph1)**2 + pk1
    zpos2 = pa2*(xpos2 - ph2)**2 + pk2
    #print("xpos1", xpos1)
    #print("zpos1", zpos1)
    dists1 = np.zeros(npoints)
    dists2 = np.zeros(npoints)
    conx = []
    cony = []
    conz = []
    sidex = []
    sidey = []
    sidez = []
    dists1[0] = 0
    pad = 1.5
    extra = 0
    nspacing = spacing
    for i in range(npoints-1):
        dists1[i+1] = dists1[i] + ((xpos1[i] - xpos1[i+1])**2 + (ypos1[i] - ypos1[i+1])**2 + (zpos1[i] - zpos1[i+1])**2)**0.5
        if dists1[i] >= nspacing:
            conx.append(xpos1[i+1])
            cony.append(ypos1[i+1])
            conz.append(zpos1[i+1])
            nspacing = nspacing + spacing
    dists2[0] = dists1[-1]
    N1 = len(conx)
    for i in range(N1):
        if i%2 == 1:
            sy = cony[i] - pad
            slope = -1*(pa1*2*(conx[i] - ph1))
            inter = conz[i] - slope*conx[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            #print("i", i,"normvec", normvec)
            sx = conx[i] + np.sign(normvec[0])*(spacing-pad + extra)*(normvec[0])
            sz = conz[i] + np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)

    for i in range(npoints-1):
        dists2[i+1] = dists2[i] + ((xpos2[i] - xpos2[i+1])**2 + (ypos2[i] - ypos2[i+1])**2 + (zpos2[i] - zpos2[i+1])**2)**0.5
        distend = ((xpos2[i+1] - end[0])**2 + (ypos2[i+1] - end[1])**2 + (zpos2[i+1] - end[2])**2)**2
        if dists2[i] >= nspacing and distend > spacing - 0.1:
            conx.append(xpos2[i+1])
            cony.append(ypos2[i+1])
            conz.append(zpos2[i+1])
            nspacing = nspacing + spacing
    N = len(conx)
    for i in range(N1, N):
        if i%2 == 1:
            sy = cony[i] + pad
            slope = -1*(pa2*2*(conx[i] - ph2))
            inter = conz[i] - slope*conx[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            sx = conx[i] -  np.sign(normvec[0])*(spacing-pad + extra)*normvec[0]
            sz = conz[i] -  np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)
    return conx, cony, conz, sidex, sidey, sidez



def linear2parabolaConnect(la,lb,pa1, ph1, pk1, pa2, ph2, pk2, start, end,mid,spacing):
    npoints = 1000
    xpos1 = np.linspace(start[0], mid[0], npoints)
    xpos2 = np.linspace(mid[0], end[0], npoints)
    ypos1 = la*xpos1 + lb*np.ones(npoints)
    ypos2 = la*xpos2 + lb*np.ones(npoints)
    zpos1 = pa1*(ypos1 - ph1)**2 + pk1
    zpos2 = pa2*(ypos2 - ph2)**2 + pk2
    dists1 = np.zeros(npoints)
    dists2 = np.zeros(npoints)
    conx = []
    cony = []
    conz = []
    sidex = []
    sidey = []
    sidez = []
    dists1[0] = 0
    pad = 0
    extra = 0
    nspacing = spacing*0.8
    for i in range(npoints-1):
        dists1[i+1] = dists1[i] + ((xpos1[i] - xpos1[i+1])**2 + (ypos1[i] - ypos1[i+1])**2 + (zpos1[i] - zpos1[i+1])**2)**0.5
        if dists1[i] > nspacing:
            conx.append(xpos1[i+1])
            cony.append(ypos1[i+1])
            conz.append(zpos1[i+1])
            nspacing = nspacing + spacing
    dists2[0] = dists1[-1]
    N1 = len(conx)
    #put side groups for first parabola
    for i in range(N1):
        if i%2 == 0:

            sx = conx[i] 

            slope = -1*(pa1*2*(cony[i] - ph1))
            inter = conz[i] - slope*cony[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            #print("i", i,"normvec", normvec)
            

            sy = cony[i] - np.sign(normvec[0])*(spacing-pad + extra)*(normvec[0])
            sz = conz[i] - np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)

    for i in range(npoints-1):
        dists2[i+1] = dists2[i] + ((xpos2[i] - xpos2[i+1])**2 + (ypos2[i] - ypos2[i+1])**2 + (zpos2[i] - zpos2[i+1])**2)**0.5
        distend = ((xpos2[i+1] - end[0])**2 + (ypos2[i+1] - end[1])**2 + (zpos2[i+1] - end[2])**2)*2
        if dists2[i] > nspacing and distend > spacing-0.1:
            conx.append(xpos2[i+1])
            cony.append(ypos2[i+1])
            conz.append(zpos2[i+1])
            nspacing = nspacing + spacing
    N = len(conx)
    for i in range(N1, N):
        if i%2 == 0:
            sx = conx[i]
            slope = -1*(pa2*2*(cony[i] - ph2))
            inter = conz[i] - slope*cony[i]
            vec = np.array([slope, 1])
            distvec = sum(vec**2)**0.5
            normvec = vec/distvec
            sy = cony[i] +  np.sign(normvec[0])*(spacing-pad + extra)*normvec[0]
            sz = conz[i] + np.sign(normvec[0])*(spacing-pad + extra)*normvec[1]
            sidex.append(sx)
            sidey.append(sy)
            sidez.append(sz)

            
    return conx, cony, conz, sidex, sidey, sidez
    


def writexyz(xyzfile, atomsinfo):
    Natoms = len(atomsinfo)
    fout = open(xyzfile, 'w')
    fout.write(str(Natoms) + "\nAtoms\n")
    for i in range(Natoms):
        fout.write(str(int(atomsinfo[i,0])) + " " + str(int(atomsinfo[i,2])) + " " + str(atomsinfo[i,3]) + " " + str(atomsinfo[i,4]) + " " + str(atomsinfo[i,5]) + "\n")
    fout.close()


def setupOldInfiniteSystem(xyzfile, Rx, Ry, Rz):
    atomsinfo = []
    a = 6.63 
    b = 20.78
    c = 6.50
    atype = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    beta = 99.5 
    diffangle = 99.5 - 90
    #find unit vector of new c axis
    zoffset = np.cos(math.radians(diffangle))
    xoffset = -np.sin(math.radians(diffangle))
    print("xoffset", xoffset, "zoffset", zoffset)
    rcoords = np.zeros((9, 3))
    rcoords[:, 0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:, 1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:, 2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    widthy = max(rcoords[:, 1]) - min(rcoords[:, 1])
    widthx = max(rcoords[:, 0]) - min(rcoords[:, 0])
    spacing = widthy / 2

    xnegmin = min(-rcoords[:, 0])
    xposmin = min(rcoords[:, 0])
    diffxmin = xposmin - xnegmin
    symmetries = np.array([[1, 1, 1], [-1, -1, -1], [1, -1 + 0.5, 1 + 0.5], [-1, 1 + 0.5, -1 + 0.5]])
    symmetriesM = np.array([[1, 1, 1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array([[0.5, 0, 0], [0, 0.5, 1], [0, 0.5, 0], [0.5, 1, 1]])
    unitcell = np.zeros((len(symmetries) * len(rcoords), 3))
    ucount = 0
    sortedInds0 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds1 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds2 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds3 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    for s in range(len(symmetries)):
        #print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            if s == 0:
                unitcell[ucount, 1] = rcoords[sortedInds0[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds0[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds0[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1

            elif s ==1:
                unitcell[ucount, 1] = rcoords[sortedInds1[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds1[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds1[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
            elif s == 2:
                unitcell[ucount, 1] = rcoords[sortedInds2[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds2[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds2[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
            elif s == 3:
                unitcell[ucount, 1] = rcoords[sortedInds3[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds3[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds3[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
    
    #rcoords = np.zeros((8,3))

    #rcoords[:,0] = [0, 0, 0.5, 0.5, 0, 0, 0.5, 0.5]
    #rcoords[:,1] = [0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5]
    #rcoords[:,2] = [0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5 ]
    #rcoords = np.array([[0,0,0]])
    
    count = 1
    mol = 1
    for k in range(Rz):
        basez =k*c
        #basez = k*c
        for j in range(Ry):
            basey = j*b
            for i in range(Rx):
                basex = i*a 
                #basex  = i*a
                for s in range(len(symmetries)):
                    mol = 1 + s + k*len(symmetries)
                    for ci in range(len(rcoords)):
                        x = basex + unitcell[s*9 + ci,0]*a
                        y = basey + unitcell[s*9 + ci,1]*b
                        z = basez + unitcell[s*9 + ci,2]*c
                        tx = x + z*xoffset
                        ty = y 
                        tz = z*zoffset

                        #z = basez + (symmetries[s][2]*c + rcoords[ci][2]*c)
                        atomsinfo.append([count, mol, atype[ci], tx, ty, tz])
                        count = count + 1
                    

    atomsinfo = np.array(atomsinfo)
    writexyz(xyzfile, atomsinfo)
    #xlo = -(Rx+1)*a
    #ylo = -(Ry+1)*b
    #zlo = -(Rz+1)*c
    #xhi =  (Rx+1)*a
    #yhi = (Ry+1)*b*(1.5)
    #zhi = (Rz+1)*c*(1.5)
    xlo = min(atomsinfo[:,3]) - 5
    xhi = max(atomsinfo[:,3]) + 5
    ylo = min(atomsinfo[:,4]) - 5
    yhi = max(atomsinfo[:,4]) + 5
    zlo = min(atomsinfo[:,5]) - 5
    zhi = max(atomsinfo[:,5]) + 5
    return atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi

def readlammpsbondsPPctypes(filename,newfile):
    Cmasstypes = [13.019, 14.027, 15.035]
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata(filename)
    #change atomtypes in based on number of bonds
    #loop through bondsinfo
    numBonds = np.zeros(len(atomsinfo))
    for i in range(len(bondsinfo)):
        bond = bondsinfo[i,:]
        atom1 = bond[2]-1
        atom2 = bond[3]-1
        numBonds[int(atom1)] += 1
        numBonds[int(atom2)] += 1

    #update atomsinfo
    for i in range(len(atomsinfo)):
        atomsinfo[i,2] = numBonds[i]

    #write new lammps data with just atoms 
    msc.writelammpsdatajustatoms(newfile, boxcoords, Cmasstypes, len(atomsinfo), atomsinfo)



def setupInfiniteSystem(xyzfile, Rx, Ry, Rz):
    atomsinfo = []
    a = 6.63 
    b = 20.78
    c = 6.50
    atype = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    beta = 99.5
    #diffangle = 0
    diffangle = 99.5 - 90
    #find unit vector of new c axis
    zoffset = np.cos(math.radians(diffangle))
    xoffset = -np.sin(math.radians(diffangle))

    print("xoffset", xoffset, "zoffset", zoffset) 
    print("unit 001 in monoclinic", "x", xoffset, "y", 0, "z", zoffset) 
    rcoords = np.zeros((9,3))
    rcoords[:,0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:,1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:,2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    
    #rcoords = np.zeros((8,3))
    widthy = max(rcoords[:,1]) - min(rcoords[:,1])
    widthx = max(rcoords[:,0]) - min(rcoords[:,0])
    print("widthx", widthx)
    print("widthy",widthy)
    spacing = widthy/100
    print("widthy*4",widthy*4)

    #rcoords[:,0] = [0, 0, 0.5, 0.5, 0, 0, 0.5, 0.5]
    #rcoords[:,1] = [0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5]
    #rcoords[:,2] = [0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5 ]
    #rcoords = np.array([[0,0,0]])
    xnegmin = min(-rcoords[:,0])
    xposmin = min(rcoords[:,0])
    diffxmin = xposmin - xnegmin
    print("xnegmin", xnegmin, "xposmin", xposmin, "diffxmin", diffxmin)
    symmetries = np.array([[1,1,1], [-1,-1,-1], [1,-1+0.5, 1+0.5], [-1, 1+0.5, -1+0.5]])
    symmetriesM = np.array([[1,1,1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array([[0.5,0,0], [0, 0.5, 1], [0, 0.5, 0], [0.5, 1, 1]])
    #symmetriesNewM = np.array([[1,1,1],[-1,-1,-1],[1,-1,1],[-1,1,-1]])
    #symmetriesNewA = np.array([[]])
    unitcell = np.zeros((len(symmetries)*len(rcoords),3))
    ucount = 0
    for s in range(len(symmetries)):
        print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            #unitcell[ucount,0] = rcoords[ci,0]*np.sign(symmetries[s,0]) + np.sign(symmetries[s,0])*abs(abs(symmetries[s,0])-1) 
            #unitcell[ucount, 0] = rcoords[ci,0]*symmetriesM[s,0] + symmetriesA[s,0]
            unitcell[ucount, 1] = (rcoords[ci,1]*symmetriesM[s,1]+symmetriesA[s,1])
            unitcell[ucount, 2] = (rcoords[ci,2]*symmetriesM[s,2]+symmetriesA[s,2])
            unitcell[ucount, 0] = (rcoords[ci,0]*symmetriesM[s,0] +symmetriesA[s,0])
            #unitcell[ucount,1] = rcoords[ci,1]*np.sign(symmetries[s,1]) + np.sign(symmetries[s,1])*abs(abs(symmetries[s,1])-1)
            #unitcell[ucount,2] = rcoords[ci,2]*np.sign(symmetries[s,2]) + np.sign(symmetries[s,2])*abs(abs(symmetries[s,2])-1)
            print("unitcell[ucount]", unitcell[ucount])
            ucount +=1
    #calculate distances between ends:
    #pairs in unit cells
    
    #change unit cell order to mesh with order of beads
    first = list(unitcell[0])
    unitcell[0] = list(unitcell[1])
    unitcell[1] = first
    

    mid = list(unitcell[18])
    unitcell[18] = list(unitcell[19])
    unitcell[19] = mid
    print("unitcell", unitcell)




    ends1s = [[unitcell[8][0]*a, unitcell[8][1]*b, unitcell[8][2]*c], [unitcell[17][0]*a, unitcell[17][1]*b, unitcell[17][2]*c], [unitcell[26][0]*a, unitcell[26][1]*b, unitcell[26][2]*c+0.5], [unitcell[35][0]*a, unitcell[35][1]*b, unitcell[35][2]*c]]
    ends2s = [[unitcell[10][0]*a, unitcell[10][1]*b+0.8, unitcell[10][2]*c], [unitcell[19][0], unitcell[18][1]*b, unitcell[18][2]*c], [unitcell[28][0]*a, unitcell[28][1]*b+1, unitcell[28][2]*c], [unitcell[0,0]*a, (unitcell[0,1]+1)*b, unitcell[0,2]*c]]
    midpoints = []
    for i in range(len(ends1s)):
        #print("end1", ends1s[i], "end2", ends2s[i])
        xdist = (ends1s[i][0] - ends2s[i][0])*a
        ydist = (ends1s[i][1] - ends2s[i][1])*b
        zdist = (ends1s[i][2] - ends2s[i][2])*c
        print("xdist", xdist, "ydist", ydist, "zdist", zdist)
        dist = (xdist**2 + ydist**2 + zdist**2)**0.5
        print("connection", i+1, "dist", dist)
    for i in range(len(ends1s)):
        print("end1", ends1s[i][0], ends1s[i][1], ends1s[i][2])
        print("end2", ends2s[i][0], ends2s[i][1], ends2s[i][2])
        if i%2 == 0:
            print("midpoint", (ends1s[i][0]+ ends2s[i][0])/2, (ends1s[i][1] + ends2s[i][1])/2 , (ends1s[i][2] + ends2s[i][2])/2)
            midpoints.append([(ends1s[i][0] + ends2s[i][0])/2, (ends1s[i][1] + ends2s[i][1])/2 , (ends1s[i][2] + ends2s[i][2])/2 + c/2])
            print("midpoints", midpoints)
        else:
            print("midpoint", (ends1s[i][0] + ends2s[i][0])/2, (ends1s[i][1] + ends2s[i][1])/2 , (ends1s[i][2] + ends2s[i][2])/2)
            midpoints.append([(ends1s[i][0] + ends2s[i][0])/2, (ends1s[i][1] + ends2s[i][1])/2 , (ends1s[i][2] + ends2s[i][2])/2 - c/2])
            print("midpoints", midpoints)

    #connection 1 
    #parabola in z y
    # linear in x
    Aparabola1 = []
    Hparabola1 = []
    Kparabola1 = []
    Aparabola2 = []
    Hparabola2 = []
    Kparabola2 = []
    Aline = []
    Bline = []

    print('ends1s', ends1s, 'ends1s', ends2s, 'midpoints', midpoints)
    for i in range(len(ends1s)):
        print(ends1s[i][1], midpoints[i][1], ends2s[i][1], ends1s[i][2], midpoints[i][2], ends2s[i][2])
        ph1  = midpoints[i][1]
        pk1  = midpoints[i][2]
        pa1 = (ends1s[i][2] - pk1)/((ends1s[i][1] - ph1)**2)
        Aparabola1.append(pa1)
        Hparabola1.append(ph1)
        Kparabola1.append(pk1)
        ph2 = midpoints[i][1]
        pk2 = midpoints[i][2]
        pa2 = (ends2s[i][2] - pk2)/((ends2s[i][1] - ph2)**2)
        Aparabola2.append(pa2)
        Hparabola2.append(ph2)
        Kparabola2.append(pk2)
        la = (ends2s[i][1] - ends1s[i][1])/(ends2s[i][0] - ends1s[i][0])
        print("slope", la)
        lb = ends2s[i][1] - la*ends2s[i][0]
        #[la,lb] = np.polyfit([ends1s[i][1], ends2s[i][1]],[ends1s[i][0], ends2s[i][0]],1)
        Aline.append(la)
        Bline.append(lb) 
        
    connections1 = []
    connections2 = []
    connections3 = []
    connections4 = []
    connectionsR1 = []
    connectionsR2 = []
    connectionsR3 = []
    connectionsR4 = []
    npoints = [8, 6, 8, 6]
    spacing = [1.4, 1.4, 1.4, 1.4, 1.54]
    for i in range(4):
        conx, cony, conz, sidex, sidey, sidez = linear2parabolaConnect(Aline[i], Bline[i], Aparabola1[i], Hparabola1[i], Kparabola1[i], Aparabola2[i], Hparabola2[i], Kparabola2[i], ends1s[i], ends2s[i], midpoints[i], spacing[i])
        sidec = 0
        sidecR = len(sidex) -1
        for j in range(len(conx)):
            vars()["connections" + str(i + 1)].append([conx[j], cony[j], conz[j]])
            if j % 2 == 0:
                vars()["connections" + str(i + 1)].append([sidex[sidec], sidey[sidec], sidez[sidec]])
                sidec += 1
            jR = len(conx) -1 - j
            vars()["connectionsR" + str(i+1)].append([conx[jR], cony[jR], conz[jR]])
            if jR %2 == 0:
                vars()["connectionsR" + str(i+1)].append([sidex[sidecR], sidey[sidecR], sidez[sidecR]])
                sidecR = sidecR - 1
        #for j in range(len(conx)//2, len(conx)):
            #if j % 2 == 0:
                #vars()["connections" + str(i+1)].append([sidex[sidec], sidey[sidec], sidez[sidec]])
                #sidec +=1

            #vars()["connections" + str(i+1)].append([conx[j], cony[j], conz[j]])

    

    #planes connection
    
    pend1 = list(ends1s[3])
    pend2 = [ends1s[3][0]+a, ends1s[3][1], ends1s[3][2]]
    pla = (pend2[1] - pend1[1])/(pend2[0] - pend1[0])
    plb = pend2[1] - pla*pend2[0] 
    pmid = [(pend1[0] + pend2[0])/2, (pend1[1] + pend2[1])/2 , (pend1[2] + pend2[2])/2 - c/2]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]

    penda1 = (pend1[2] - pendk1)/((pend1[0] - pendh1)**2)
    penda2 = (pend2[2] - pendk2)/((pend2[0] - pendh2)**2)
    pconx, pcony, pconz, psidex, psidey, psidez = linear2parabolaXZConnect(pla,plb,penda1, pendh1, pendk1, penda2, pendh2, pendk2, pend1, pend2, pmid,1.6)
    print("pla", pla, "plb", plb)
    print("pend1", pend1)
    print("pend2", pend2)
    print("pmid", pmid)
    print("penda1", penda1, "penda2", penda2)


    p2end1 = [unitcell[0][0]*a , unitcell[0][1]*b, unitcell[0][2]*c]
    p2end2 = [unitcell[0][0]*a + a, p2end1[1]-0.5, p2end1[2]]
    pla = (p2end2[1] - p2end1[1])/(p2end2[0] - p2end1[0])
    plb = p2end2[1] - pla*p2end2[0]
    pmid = [(p2end1[0] + p2end2[0])/2, (p2end1[1] + p2end2[1])/2, (p2end1[2] + p2end2[2])/2 - c/2]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]
    
    penda1 = (p2end1[2] - pendk1)/((p2end1[0] - pendh1)**2)
    penda2 = (p2end2[2] - pendk2)/((p2end2[0] - pendh2)**2)
    print("p2end1", p2end1, "p2end2", p2end2)
    print("penda1", penda1, "penda2", penda2)
    p2conx, p2cony, p2conz,p2sidex, p2sidey, p2sidez = linear2parabolaXZConnect2(pla, plb, penda1, pendh1, pendk1, penda2, pendh2, pendk2, p2end1, p2end2, pmid, 1.6)
    print("p2sidex", p2sidex, "p2sidey", p2sidey, "p2sidez", p2sidez)
    print("p2conx", p2conx, "p2cony", p2cony, "p2conz", p2conz)
    pconnections1 = []
    cside = 0
    for i in range(len(pconx)):
        pconnections1.append([pconx[i], pcony[i], pconz[i]])
        if i % 2 == 0:
            pconnections1.append([psidex[cside], psidey[cside], psidez[cside]])
            cside += 1
    #for i in range(len(pconx)//2, len(pconx)):
        #if i % 2 == 0:
            #pconnections1.append([psidex[cside], psidey[cside], psidez[cside]])
            #cside += 1
        #pconnections1.append([pconx[i], pcony[i], pconz[i]])

    pconnections2 = []
    cside = 0
    for i in range(len(p2conx)):
        pconnections2.append([p2conx[i], p2cony[i], p2conz[i]])
        if i % 2 == 1:
            pconnections2.append([p2sidex[cside], p2sidey[cside], p2sidez[cside]])
            cside += 1
    #for i in range(len(p2conx)//2, len(p2conx)):
        #if i % 2 == 0:
            #pconnections2.append([p2sidex[cside], p2sidey[cside], p2sidez[cside]])
            #cside += 1
        #pconnections2.append([p2conx[i], p2cony[i], p2conz[i]])

    rcoords = np.zeros((9, 3))
    rcoords[:, 0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:, 1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:, 2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    widthy = max(rcoords[:, 1]) - min(rcoords[:, 1])
    widthx = max(rcoords[:, 0]) - min(rcoords[:, 0])
    spacing = widthy / 2

    xnegmin = min(-rcoords[:, 0])
    xposmin = min(rcoords[:, 0])
    diffxmin = xposmin - xnegmin
    symmetries = np.array([[1, 1, 1], [-1, -1, -1], [1, -1 + 0.5, 1 + 0.5], [-1, 1 + 0.5, -1 + 0.5]])
    symmetriesM = np.array([[1, 1, 1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array([[0.5, 0, 0], [0, 0.5, 1], [0, 0.5, 0], [0.5, 1, 1]])
    unitcell = np.zeros((len(symmetries) * len(rcoords), 3))
    ucount = 0
    sortedInds0 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds1 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds2 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds3 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    for s in range(len(symmetries)):
        #print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            if s == 0:
                unitcell[ucount, 1] = rcoords[sortedInds0[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds0[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds0[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1

            elif s ==1:
                unitcell[ucount, 1] = rcoords[sortedInds1[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds1[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds1[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
            elif s == 2:
                unitcell[ucount, 1] = rcoords[sortedInds2[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds2[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds2[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
            elif s == 3:
                unitcell[ucount, 1] = rcoords[sortedInds3[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds3[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds3[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1

    


    #Reverse
    unitcellReverse = np.copy(unitcell)
    started = True
    for i in range(36):
        unitcellReverse[i,:] = list(unitcell[35-i,:])
    sortedRInds = [0, 2, 1, 3, 5, 4, 6, 8, 7]
    unitcellRcopy = np.copy(unitcellReverse)
    for s in range(len(symmetries)):
        for ci in range(len(rcoords)): 
            unitcellReverse[s*9 + ci,:] = list(unitcellRcopy[s*9 + sortedRInds[ci],:])

  
    count = 1
    srclist = [2, 1, 0, 3]
    connect = True
    for i in range(Rx):
        basex = i*a
        if i%2 == 0:
            Rylist = list(range(Ry))
            reverse = False
            unitcellnew = np.copy(unitcell)
        else:
            Rylist = list(range(Ry-1, -1, -1))
            reverse = True
            unitcellnew = np.copy(unitcellReverse)
            #Rylist = list(range(Ry))
        for j in Rylist:
            basey = j*b
            for s in range(4):
                if (s%2 == 0):
                    Rzlist = list(range(Rz))
                    
                else:
                    Rzlist = list(range(Rz-1, -1,-1))
                    
                for k in Rzlist:
                    basez = k*c
                    
                    
                    for ci in range(9):
                            mol = ci//9 + 1
                            #x = basex + (symmetries[s][0]*a + rcoords[ci][0]*a) + xoffset*(symmetries[s][2]*c + rcoords[ci][2]*c)*a
                            x = basex + unitcellnew[s*9 + ci,0]*a 
                            #print("x", x,"basex", basex, "symmetries[s][0]*a", symmetries[s][0]*a, "symmetries[s][0]", symmetries[s][0])
                            y = basey + unitcellnew[s*9 + ci,1]*b
                            #print("y", y, "basey", basey, "symmetries[s][1]*b", symmetries[s][1]*b, "symmetries[s][1]", symmetries[s][1])
                            #z = basez + (symmetries[s][2]*c + rcoords[ci][2]*c)*zoffset 
                            z = basez + unitcellnew[s*9 + ci,2]*c
                            tx = x + z*xoffset
                            ty = y
                            tz = z*zoffset
                            atomsinfo.append([count, s+1, 1, tx, ty, tz])
                            count = count + 1

                    if reverse == False:
                        if k == Rzlist[-1] and s%2 == 0 and connect:
                            for ii in range(len(vars()["connections"+str(s+1)])):
                                x = basex + vars()["connections" + str(s+1)][ii][0]
                                y = basey + vars()["connections" + str(s+1)][ii][1]
                                z = basez + vars()["connections" + str(s+1)][ii][2]
                                tx = x + z*xoffset
                                ty = y 
                                tz = z*zoffset
                                atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                                count = count + 1
                        if k == Rzlist[-1] and s%2 ==1 and (s != 3 or j < Ry-1) and connect:
                            for ii in range(len(vars()["connections"+str(s+1)])):
                                x = basex + vars()["connections" + str(s+1)][ii][0]
                                y = basey + vars()["connections" + str(s+1)][ii][1]
                                z = vars()["connections" + str(s+1)][ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                                count = count + 1
                        if k == Rzlist[-1] and s==3 and j == Rylist[-1] and i%2 ==0 and connect:
                            for ii in range(len(pconnections1)):
                                x = basex + pconnections1[ii][0]
                                y = basey + pconnections1[ii][1]
                                z = pconnections1[ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset

                                atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                                count = count + 1

                        elif k == Rzlist[-1] and s==3 and j == Rylist[-1] and i%2 ==1 and connect:

                            for ii in range(len(pconnections2)):
                                x = basex  + pconnections2[ii][0]
                                y = pconnections2[ii][1]
                                z = pconnections2[ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                                count = count + 1

                    else: 
                        sr = srclist[s]
                        #print("sr", sr)
                        if k == Rzlist[-1] and s==3 and j == Rylist[-1] and i%2 ==0:
                            for ii in range(len(pconnections1)):
                                x = basex + pconnections1[ii][0]
                                y = basey + pconnections1[ii][1]
                                z = pconnections1[ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                
                                atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                                    
                                count = count + 1

                        elif k == Rzlist[-1] and s==3 and j == Rylist[-1] and i%2 ==1:
                            
                            for ii in range(len(pconnections2)):
                                x = basex  + pconnections2[ii][0]
                                y = pconnections2[ii][1]
                                z = pconnections2[ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                
                                atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                                    
                                
                                count = count + 1
                        elif k == Rzlist[-1] and sr%2 == 0:
                            
                            for ii in range(len(vars()["connectionsR"+str(sr+1)])):
                                x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                y = basey + vars()["connectionsR" + str(sr+1)][ii][1]
                                z = basez + vars()["connectionsR" + str(sr+1)][ii][2]
                                tx = x + z*xoffset
                                ty = y 
                                tz = z*zoffset
                                atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                                    
                                count = count + 1
                        elif k == Rzlist[-1] and sr%2 ==1 and (sr != 3 or j != Rylist[-1]):
                            
                            for ii in range(len(vars()["connectionsR"+str(sr+1)])):
                                if sr != 3:
                                    x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                    y = basey + vars()["connectionsR" + str(sr+1)][ii][1]
                                    z = vars()["connectionsR" + str(sr+1)][ii][2]
                                elif sr == 3:
                                    x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                    y = basey + vars()["connectionsR" + str(sr+1)][ii][1] - b
                                    z = vars()["connectionsR" + str(sr+1)][ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                atomsinfo.append([count, 5+s, 1, tx, ty, tz])
                                count = count + 1
                        
                    

    atomsinfo = np.array(atomsinfo)
    writexyz(xyzfile, atomsinfo)
    xlo = min(atomsinfo[:,3]) - 5
    xhi = max(atomsinfo[:,3]) + 5
    ylo = min(atomsinfo[:,4]) - 5
    yhi = max(atomsinfo[:,4]) + 5
    zlo = min(atomsinfo[:,5]) - 5
    zhi = max(atomsinfo[:,5]) + 5
    #xlo = -(Rx+1)*a
    #ylo = -(Ry+1)*b
    #zlo = -(Rz+1)*c
    #xhi =  (Rx+1)*a
    #yhi = (Ry+1)*b*(1.5)
    #zhi = (Rz+1)*c*(1.5)
    #xlo = -2*a + Rx*c*zoffset
    #xhi = Rx*a 
    #ylo = -2*b
    #yhi = Rx*b
    #zlo = -2*c*zoffset
    #zhi = Rz*c*zoffset
    return atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi

def readlammpsbondsPPctypes(filename,newfile):
    Cmasstypes = [15.035, 14.027, 13.019]
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata(filename)
    #change atomtypes in based on number of bonds
    #loop through bondsinfo
    numBonds = np.zeros(len(atomsinfo))
    print("bondsinfo", bondsinfo)
    for i in range(len(bondsinfo)):
        bond = bondsinfo[i,:]
        atom1 = int(bond[2])-1
        atom2 = int(bond[3])-1
        if atom1 < 5 or atom2 < 5:
            print("atom1", atom1, "atom2", atom2)
            print("bond", bond)
        numBonds[int(atom1)] = numBonds[int(atom1)] + 1
        numBonds[int(atom2)] = numBonds[int(atom2)] + 1

    #update atomsinfo
    for i in range(len(atomsinfo)):
        atomsinfo[i,2] = numBonds[i]

    #write new lammps data with just atoms 
    msc.writelammpsdatajustatoms(newfile, boxcoords, Cmasstypes, len(atomsinfo), atomsinfo)
    
    return atomsinfo 

def readxyz(filename):
    count = 0 
    with open(filename, 'r') as f:
        for i,line in enumerate(f):
            sstring = line.split()
            if i==0:
                N = int(sstring[0])
                atomsinfo = np.ones((N,6))
            elif i > 1:
                
                atomsinfo[count,0] = int(sstring[0])
                atomsinfo[count,1] = 1
                atomsinfo[count,2] = int(sstring[1])
                atomsinfo[count,3] = float(sstring[2])
                atomsinfo[count,4] = float(sstring[3])
                atomsinfo[count,5] = float(sstring[4])
                count = count + 1

    xlo = min(atomsinfo[:,3]) - 5
    xhi = max(atomsinfo[:,3]) + 5
    ylo = min(atomsinfo[:,4]) - 5
    yhi = max(atomsinfo[:,4]) + 5
    zlo = min(atomsinfo[:,5]) - 5
    zhi = max(atomsinfo[:,5]) + 5

    return atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi

def writeMoleculesinCrystalNoCons(xyzfile, moleculeClist, Rx, Ry, Rz):
    shownCs = []
    maxC = 0
    spaceC = [5 if m%2 == 0 else 5 for m in moleculeClist]
    nBeads = sum(moleculeClist)
    print("nBeads", nBeads)
    nMols  = len(moleculeClist)
    print("nMols", nMols)
    
    countRef = np.ones(nBeads + sum(spaceC[0:len(spaceC)-1]))

    print("len(countRef)", len(countRef))
    #from moleculeClist
    for mi,m in enumerate(moleculeClist):
        assert (m-2)%3 == 0
        clist = list(range(maxC, maxC + m))
        shownCs.extend(clist)
        for i in range(spaceC[mi]):
            if mi < nMols -1:
                countRef[maxC+m+i] = 0
        maxC = maxC + m + spaceC[mi]
    
    #print("shownCs", shownCs)


    atomsinfo = []
    a = 6.63
    b = 20.78
    c = 6.50
    atype = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    beta = 99.5
    diffangle = 99.5 - 90
    # find unit vector of new c axis
    zoffset = np.cos(math.radians(diffangle))
    xoffset = -np.sin(math.radians(diffangle))
    rcoords = np.zeros((9, 3))
    rcoords[:, 0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:, 1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:, 2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    widthy = max(rcoords[:, 1]) - min(rcoords[:, 1])
    widthx = max(rcoords[:, 0]) - min(rcoords[:, 0])
    spacing = widthy / 2

    xnegmin = min(-rcoords[:, 0])
    xposmin = min(rcoords[:, 0])
    diffxmin = xposmin - xnegmin
    symmetries = np.array([[1, 1, 1], [-1, -1, -1], [1, -1 + 0.5, 1 + 0.5], [-1, 1 + 0.5, -1 + 0.5]])
    symmetriesM = np.array([[1, 1, 1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array([[0.5,0,0], [0, 0.5, 1], [0, 0.5, 0], [0.5, 1, 1]])
    #symmetriesNewM = np.array([[1,1,1],[-1,-1,-1],[1,-1,1],[-1,1,-1]])
   
    unitcell = np.zeros((len(symmetries) * len(rcoords), 3))
    ucount = 0
    for s in range(len(symmetries)):
        #print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            unitcell[ucount, 1] = rcoords[ci, 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
            unitcell[ucount, 2] = (rcoords[ci, 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
            unitcell[ucount, 0] = rcoords[ci, 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
            ucount += 1
    # calculate distances between ends:
    # pairs in unit cells
    #print("unitcell before", unitcell)

    #change unit cell order to mesh with order of beads
    first = list(unitcell[0,:])
    #print("first", first)
    #print("second", unitcell[1])
    unitcell[0,:] = list(unitcell[1])
    unitcell[1,:] = first
    # swich another unit cell
    #switch  = list(unitcell[11])
    #unitcell[11] = list(unitcell[12])
    #unitcell[12] = switch

    mid = list(unitcell[18,:])
    unitcell[18,:] = list(unitcell[19,:])
    unitcell[19,:] = mid
    #print("unitcell", unitcell)
    ends1s = [[unitcell[8][0]*a, unitcell[8][1]*b, unitcell[8][2]*c], [unitcell[17][0]*a, unitcell[17][1]*b, unitcell[17][2]*c], [unitcell[26][0]*a, unitcell[26][1]*b, unitcell[26][2]*c+0.5], [unitcell[35][0]*a, unitcell[35][1]*b, unitcell[35][2]*c]]
    ends2s = [[unitcell[10][0]*a, unitcell[10][1]*b+0.8, unitcell[10][2]*c], [unitcell[19][0], unitcell[18][1]*b, unitcell[18][2]*c], [unitcell[28][0]*a, unitcell[28][1]*b+1, unitcell[28][2]*c], [unitcell[0,0]*a, (unitcell[0,1]+1)*b, unitcell[0,2]*c]]

    #ends1s = [[unitcell[8][0]*a, unitcell[8][1]*b, unitcell[8][2]*c], [unitcell[17][0]*a, unitcell[17][1]*b, unitcell[17][2]*c], [unitcell[26][0]*a, unitcell[26][1]*b, unitcell[26][2]*c+0.5], [unitcell[35][0]*a, unitcell[35][1]*b, unitcell[35][2]*c]]
    #ends2s = [[unitcell[10][0]*a, unitcell[10][1]*b+0.8, unitcell[10][2]*c], [unitcell[19][0], unitcell[18][1]*b, unitcell[18][2]*c], [unitcell[28][0]*a, unitcell[28][1]*b+1, unitcell[28][2]*c], [unitcell[0,0]*a, (unitcell[0,1]+1)*b, unitcell[0,2]*c]]
    midpoints = []

    
    for i in range(len(ends1s)):
        # print("end1", ends1s[i], "end2", ends2s[i])
        xdist = (ends1s[i][0] - ends2s[i][0]) * a
        ydist = (ends1s[i][1] - ends2s[i][1]) * b
        zdist = (ends1s[i][2] - ends2s[i][2]) * c
        dist = (xdist ** 2 + ydist ** 2 + zdist ** 2) ** 0.5
    for i in range(len(ends1s)):
        if i % 2 == 0:
            midpoints.append([(ends1s[i][0] + ends2s[i][0]) / 2, (ends1s[i][1] + ends2s[i][1]) / 2,(ends1s[i][2] + ends2s[i][2]) / 2 + c/2])
        else:
            midpoints.append([(ends1s[i][0] + ends2s[i][0]) / 2, (ends1s[i][1] + ends2s[i][1]) / 2,  (ends1s[i][2] + ends2s[i][2]) / 2 - c/2])

    # connection 1
    # parabola in z y
    # linear in x
    Aparabola1 = []
    Hparabola1 = []
    Kparabola1 = []
    Aparabola2 = []
    Hparabola2 = []
    Kparabola2 = []
    Aline = []
    Bline = []

    for i in range(len(ends1s)):
        #print(ends1s[i][1], midpoints[i][1], ends2s[i][1], ends1s[i][2], midpoints[i][2], ends2s[i][2])
        ph1 = midpoints[i][1]
        pk1 = midpoints[i][2]
        pa1 = (ends1s[i][2] - pk1) / ((ends1s[i][1] - ph1) ** 2)
        Aparabola1.append(pa1)
        Hparabola1.append(ph1)
        Kparabola1.append(pk1)
        ph2 = midpoints[i][1]
        pk2 = midpoints[i][2]
        pa2 = (ends2s[i][2] - pk2) / ((ends2s[i][1] - ph2) ** 2)
        Aparabola2.append(pa2)
        Hparabola2.append(ph2)
        Kparabola2.append(pk2)
        la = (ends2s[i][1] - ends1s[i][1]) / (ends2s[i][0] - ends1s[i][0])
        lb = ends2s[i][1] - la * ends2s[i][0]
        # [la,lb] = np.polyfit([ends1s[i][1], ends2s[i][1]],[ends1s[i][0], ends2s[i][0]],1)
        Aline.append(la)
        Bline.append(lb)

    connections1 = []
    connections2 = []
    connections3 = []
    connections4 = []
    connectionsR1 = []
    connectionsR2 = []
    connectionsR3 = []
    connectionsR4 = []
    npoints = [8, 6, 8, 6]
    spacing = [1.4, 1.4, 1.4, 1.4, 1.54]
    
    for i in range(4):
        conx, cony, conz, sidex, sidey, sidez = linear2parabolaConnect(Aline[i], Bline[i], Aparabola1[i], Hparabola1[i], Kparabola1[i], Aparabola2[i], Hparabola2[i], Kparabola2[i], ends1s[i], ends2s[i], midpoints[i], spacing[i])
        sidec = 0
        sidecR = len(sidex) -1
        for j in range(len(conx)):
            vars()["connections" + str(i + 1)].append([conx[j], cony[j], conz[j]])
            if j % 2 == 0:
                vars()["connections" + str(i + 1)].append([sidex[sidec], sidey[sidec], sidez[sidec]])
                sidec += 1
            jR = len(conx) -1 - j
            vars()["connectionsR" + str(i+1)].append([conx[jR], cony[jR], conz[jR]])
            if jR %2 == 0:
                vars()["connectionsR" + str(i+1)].append([sidex[sidecR], sidey[sidecR], sidez[sidecR]])
                sidecR = sidecR - 1
        #for j in range(len(conx)//2, len(conx)):
            #if j % 2 == 0:
                #vars()["connections" + str(i+1)].append([sidex[sidec], sidey[sidec], sidez[sidec]])
                #sidec +=1

            #vars()["connections" + str(i+1)].append([conx[j], cony[j], conz[j]])

    pend1 = list(ends1s[3])
    pend2 = [ends1s[3][0]+a, ends1s[3][1], ends1s[3][2]]

    
    pla = (pend2[1] - pend1[1]) / (pend2[0] - pend1[0])
    plb = pend2[1] - pla * pend2[0]
    pmid = [(pend1[0] + pend2[0]) / 2, (pend1[1] + pend2[1]) / 2, (pend1[2] + pend2[2]) / 2 - c/2]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]

    penda1 = (pend1[2] - pendk1) / ((pend1[0] - pendh1) ** 2)
    penda2 = (pend2[2] - pendk2) / ((pend2[0] - pendh2) ** 2)
    pconx, pcony, pconz, psidex, psidey, psidez = linear2parabolaXZConnect(pla, plb, penda1, pendh1, pendk1, penda2,
                                                                           pendh2, pendk2, pend1, pend2, pmid, 1.6)


    p2end1 = [unitcell[0][0]*a , unitcell[0][1]*b, unitcell[0][2]*c]
    p2end2 = [unitcell[0][0]*a + a, p2end1[1]-0.5, p2end1[2]]
    pla = (p2end2[1] - p2end1[1]) / (p2end2[0] - p2end1[0])
    plb = p2end2[1] - pla * p2end2[0]
    pmid = [(p2end1[0] + p2end2[0]) / 2, (p2end1[1] + p2end2[1]) / 2, (p2end1[2] + p2end2[2]) / 2 - c/2]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]

    penda1 = (p2end1[2] - pendk1) / ((p2end1[0] - pendh1) ** 2)
    penda2 = (p2end2[2] - pendk2) / ((p2end2[0] - pendh2) ** 2)
    p2conx, p2cony, p2conz, p2sidex, p2sidey, p2sidez = linear2parabolaXZConnect2(pla, plb, penda1, pendh1, pendk1,
                                                                                 penda2, pendh2, pendk2, p2end1, p2end2,
                                                                                 pmid, 1.6)
    pconnections1 = []
    cside = 0
    for i in range(len(pconx)):
        pconnections1.append([pconx[i], pcony[i], pconz[i]])
        if i % 2 == 0:
            pconnections1.append([psidex[cside], psidey[cside], psidez[cside]])
            cside += 1
    #for i in range(len(pconx)//2, len(pconx)):
        #if i % 2 == 0:
            #pconnections1.append([psidex[cside], psidey[cside], psidez[cside]])
            #cside += 1
        #pconnections1.append([pconx[i], pcony[i], pconz[i]])

    pconnections2 = []
    cside = 0
    for i in range(len(p2conx)):
        pconnections2.append([p2conx[i], p2cony[i], p2conz[i]])
        if i % 2 == 1:
            pconnections2.append([p2sidex[cside], p2sidey[cside], p2sidez[cside]])
            cside += 1
    #for i in range(len(p2conx)//2, len(p2conx)):
        #if i % 2 == 0:
            #pconnections2.append([p2sidex[cside], p2sidey[cside], p2sidez[cside]])
            #cside += 1
        #pconnections2.append([p2conx[i], p2cony[i], p2conz[i]])

    count = 0
    acount = 1
    mol = 1
    moli = 0
    madded = 0
    minspace = 5
    spacecount = minspace
    alldone = False
    nc = 0
    rcoords = np.zeros((9, 3))
    rcoords[:, 0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:, 1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:, 2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    widthy = max(rcoords[:, 1]) - min(rcoords[:, 1])
    widthx = max(rcoords[:, 0]) - min(rcoords[:, 0])
    spacing = widthy / 2

    xnegmin = min(-rcoords[:, 0])
    xposmin = min(rcoords[:, 0])
    diffxmin = xposmin - xnegmin
    symmetries = np.array([[1, 1, 1], [-1, -1, -1], [1, -1 + 0.5, 1 + 0.5], [-1, 1 + 0.5, -1 + 0.5]])
    symmetriesM = np.array([[1, 1, 1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array([[0.5,0,0], [0, 0.5, 1], [0, 0.5, 0], [0.5, 1, 1]]) 
    
    unitcell = np.zeros((len(symmetries) * len(rcoords), 3))
    ucount = 0
    sortedInds0 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds1 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds2 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds3 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    for s in range(len(symmetries)):
        #print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            if s == 0:
                unitcell[ucount, 1] = rcoords[sortedInds0[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds0[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds0[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1

            elif s ==1:
                unitcell[ucount, 1] = rcoords[sortedInds1[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds1[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds1[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
            elif s == 2:
                unitcell[ucount, 1] = rcoords[sortedInds2[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds2[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds2[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
            elif s == 3:
                unitcell[ucount, 1] = rcoords[sortedInds3[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds3[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds3[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1

    


    #Reverse
    unitcellReverse = np.copy(unitcell)
    started = True
    for i in range(36):
        unitcellReverse[i,:] = list(unitcell[35-i,:])
    sortedRInds = [0, 2, 1, 3, 5, 4, 6, 8, 7]
    unitcellRcopy = np.copy(unitcellReverse)
    for s in range(len(symmetries)):
        for ci in range(len(rcoords)): 
            unitcellReverse[s*9 + ci,:] = list(unitcellRcopy[s*9 + sortedRInds[ci],:])

    """
    for i in range(9):
        unitcellReverse[i,:] = list(unitcell[i+27,:])
        unitcellReverse[9+i,:] = list(unitcell[i+18,:])
        unitcellReverse[18+i,:] = list(unitcell[i+9,:])
        unitcellReverse[27+i,:] = list(unitcell[i,:])
    """
    srclist = [2, 1, 0, 3]
    avgxS0 = np.mean(unitcell[0:9,0])*a
    avgyS0 = np.mean(unitcell[0:9,1])*b
    avgzS0 = np.mean(unitcell[0:9,2])*c
    
    avgxS1 = np.mean(unitcell[9:18,0])*a
    avgyS1 = np.mean(unitcell[9:18,1])*b
    avgzS1 = np.mean(unitcell[9:18,2])*c
    
    avgxS2 = np.mean(unitcell[18:27,0])*a
    avgyS2 = np.mean(unitcell[18:27,1])*b
    avgzS2 = np.mean(unitcell[18:27,2])*c
    
    avgxS3 = np.mean(unitcell[27:36,0])*a
    avgyS3 = np.mean(unitcell[27:36,1])*b
    avgzS3 = np.mean(unitcell[27:36,2])*c

    dist01 = ((avgxS0 - avgxS1)**2 + (avgyS0 - avgyS1)**2 + (avgzS0 - avgzS1)**2)**0.5
    dist12 = ((avgxS2 - avgxS1)**2 + (avgyS2 - avgyS1)**2 + (avgzS2 - avgzS1)**2)**0.5
    dist23 = ((avgxS2 - avgxS3)**2 + (avgyS2 - avgyS3)**2 + (avgzS2 - avgzS3)**2)**0.5
    print(avgxS0, avgyS0, avgzS0)
    print(avgxS1, avgyS1, avgzS1)
    print(avgxS2, avgyS2, avgzS2)
    print(avgxS3, avgyS3, avgzS3)
    print("dist01", dist01, "dist12", dist12, "dist23", dist23)



    for i in range(Rx):
        basex = i*a
        if i%2 == 0:
            Rylist = list(range(Ry))
            reverse = False
            unitcellnew = np.copy(unitcell)
        else:
            Rylist = list(range(Ry-1, -1, -1))
            reverse = True
            unitcellnew = np.copy(unitcellReverse)
            #Rylist = list(range(Ry))
        for j in Rylist:
            basey = j*b
            for s in range(4):
                if (s%2 == 0):
                    Rzlist = list(range(Rz))
                    
                else:
                    Rzlist = list(range(Rz-1, -1,-1))
                    
                for k in Rzlist:
                    basez = k*c
                    
                    
                    for ci in range(9):
                            if not alldone and madded == 0 and spacecount >= minspace:
                                if not reverse:
                                    if s == 0 and ci == 0:
                                        started = True
                                    elif s == 1 and ci == 0:
                                        started = True
                                    elif s == 2 and ci == 0:
                                        started = True
                                    elif s == 3 and ci == 0:
                                        started == True
                                else:
                                    if s== 0 and ci == 0:
                                        started = True
                                    elif s == 1 and ci == 0:
                                        started = True
                                    elif s == 2 and ci == 0:
                                        started = True
                                    elif s == 3 and ci == 0:
                                        started = True

                            
                            #x = basex + (symmetries[s][0]*a + rcoords[ci][0]*a) + xoffset*(symmetries[s][2]*c + rcoords[ci][2]*c)*a
                            x = basex + unitcellnew[s*9 + ci,0]*a 
                            #print("x", x,"basex", basex, "symmetries[s][0]*a", symmetries[s][0]*a, "symmetries[s][0]", symmetries[s][0])
                            y = basey + unitcellnew[s*9 + ci,1]*b
                            #print("y", y, "basey", basey, "symmetries[s][1]*b", symmetries[s][1]*b, "symmetries[s][1]", symmetries[s][1])
                            #z = basez + (symmetries[s][2]*c + rcoords[ci][2]*c)*zoffset 
                            z = basez + unitcellnew[s*9 + ci,2]*c
                            tx = x + z*xoffset
                            ty = y
                            tz = z*zoffset
                            if not alldone and madded < moleculeClist[moli] and spacecount >= minspace and started: 
                                atomsinfo.append([acount, 1+s, 1, tx, ty, tz])
                                acount = acount + 1
                                madded  = madded + 1
                                #print("madded", madded)
                                #print("spacecount", spacecount)
                            elif (madded == 0 and spacecount < minspace):
                                if spacecount == 0:
                                    if ci > 0:
                                        spacecount = spacecount + 1
                                else:
                                    spacecount = spacecount + 1
                                #print("tx", tx, "ty", ty, "tz", tz)
                               
                            if not alldone and madded == moleculeClist[moli]:
                                moli = moli + 1
                                print("moli", moli)
                                madded = 0
                                spacecount = 0
                                if moli == len(moleculeClist):
                                    alldone = True
                                minspace = 5 +nc%2
                                started = False
                                #minspace = 5 
                                #print("minspace", minspace)
                            count = count + 1
                            

                    if reverse == False:
                        if k == Rzlist[-1] and s % 2 == 0:
                            nc = nc +1
                            for ii in range(len(vars()["connections" + str(s + 1)])):
                                x = basex + vars()["connections" + str(s + 1)][ii][0]
                                y = basey + vars()["connections" + str(s + 1)][ii][1]
                                z = basez + vars()["connections" + str(s + 1)][ii][2]
                                tx = x + z * xoffset
                                ty = y
                                tz = z * zoffset
                                #if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    #atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    #acount = acount + 1
                                    #madded = madded +1
                                #count = count + 1
                        elif k == Rzlist[-1] and s % 2 == 1 and (s != 3 or j < Ry - 1):
                            for ii in range(len(vars()["connections" + str(s + 1)])):
                                x = basex + vars()["connections" + str(s + 1)][ii][0]
                                y = basey + vars()["connections" + str(s + 1)][ii][1]
                                z = vars()["connections" + str(s + 1)][ii][2]
                                tx = x + z * xoffset
                                ty = y
                                tz = z * zoffset
                                #if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    #atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    #count = acount + 1
                                    #madded = madded + 1
                                #count = count + 1
                        elif k == Rzlist[-1] and s == 3 and j == Ry - 1 and i % 2 == 0:
                            nc = nc +1
                            for ii in range(len(pconnections1)):
                                x = basex + pconnections1[ii][0]
                                y = basey + pconnections1[ii][1]
                                z = pconnections1[ii][2]
                                tx = x + z * xoffset
                                ty = y
                                tz = z * zoffset
                                #if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    #atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    #acount = acount + 1
                                    #madded = madded + 1
                                #count = count + 1

                        elif k == Rzlist[-1] and s == 3 and j == Ry - 1 and i % 2 == 1:

                            for ii in range(len(pconnections2)):
                                x = basex + pconnections2[ii][0]
                                y = pconnections2[ii][1]
                                z = pconnections2[ii][2]
                                tx = x + z * xoffset
                                ty = y
                                tz = z * zoffset
                                #if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    #atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    #acount = acount + 1
                                    #madded = madded + 1
                                #count = count + 1
                        if not alldone and madded == moleculeClist[moli]:
                            madded = 0

                            moli = moli + 1
                            print("moli", moli)
                            if moli == len(moleculeClist):
                                alldone = True
                            #print("mol", mol)
                        #print("acount", acount)
                    
                    else: 
                        sr = srclist[s]
                        #print("sr", sr)
                        if k == Rzlist[-1] and s==3 and j == Rylist[-1] and i%2 ==0:
                            for ii in range(len(pconnections1)):
                                x = basex + pconnections1[ii][0]
                                y = basey + pconnections1[ii][1]
                                z = pconnections1[ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                #if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    #atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    #acount = acount + 1
                                    #madded = madded + 1
                                #count = count + 1

                        elif k == Rzlist[-1] and s==3 and j == Rylist[-1] and i%2 ==1:
                            nc = nc + 1
                            for ii in range(len(pconnections2)):
                                x = basex  + pconnections2[ii][0]
                                y = pconnections2[ii][1]
                                z = pconnections2[ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                #if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    #atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    #acount = acount + 1
                                    #madded = madded + 1
                                #count = count + 1
                        elif k == Rzlist[-1] and sr%2 == 0:
                            nc = nc + 1
                            for ii in range(len(vars()["connectionsR"+str(sr+1)])):
                                x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                y = basey + vars()["connectionsR" + str(sr+1)][ii][1]
                                z = basez + vars()["connectionsR" + str(sr+1)][ii][2]
                                tx = x + z*xoffset
                                ty = y 
                                tz = z*zoffset
                                #if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    #atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    #acount = acount + 1
                                    #madded = madded + 1
                                #count = count + 1
                        elif k == Rzlist[-1] and sr%2 ==1 and (sr != 3 or j != Rylist[-1]):
                            
                            for ii in range(len(vars()["connectionsR"+str(sr+1)])):
                                if sr != 3:
                                    x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                    y = basey + vars()["connectionsR" + str(sr+1)][ii][1]
                                    z = vars()["connectionsR" + str(sr+1)][ii][2]
                                elif sr == 3:
                                    x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                    y = basey + vars()["connectionsR" + str(sr+1)][ii][1] - b
                                    z = vars()["connectionsR" + str(sr+1)][ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                #if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    #atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    #acount = acount + 1
                                    #madded = madded + 1
                                #count = count + 1
                        if not alldone and madded == moleculeClist[moli]:
                            madded = 0
                            moli = moli + 1
                            print("moli", moli)
                            if moli == len(moleculeClist):
                                alldone = True
                            spacecount = 0
                            minspace = 5 +nc%2
                            started = False

    if count < len(countRef) or acount != nBeads+1:
        print("ERROR: NOT ALL WERE SHOWN")
        print("count", count, "len(shownCs)", len(shownCs))
        print("acount", acount, "nBeads+1", nBeads+1)
    atomsinfo = np.array(atomsinfo)
    writexyz(xyzfile, atomsinfo)
    xlo = min(atomsinfo[:, 3]) - 5
    xhi = max(atomsinfo[:, 3]) + 5
    ylo = min(atomsinfo[:, 4]) - 5
    yhi = max(atomsinfo[:, 4]) + 5
    zlo = min(atomsinfo[:, 5]) - 5
    zhi = max(atomsinfo[:, 5]) + 5
    return atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi


def writeMoleculesinCrystal(xyzfile, moleculeClist, Rx, Ry, Rz):
    shownCs = []
    maxC = 0
    spaceC = [5 if m%2 == 0 else 5 for m in moleculeClist]
    nBeads = sum(moleculeClist)
    print("nBeads", nBeads)
    nMols  = len(moleculeClist)
    print("nMols", nMols)
    
    countRef = np.ones(nBeads + sum(spaceC[0:len(spaceC)-1]))

    print("len(countRef)", len(countRef))
    #from moleculeClist
    for mi,m in enumerate(moleculeClist):
        assert (m-2)%3 == 0
        clist = list(range(maxC, maxC + m))
        shownCs.extend(clist)
        for i in range(spaceC[mi]):
            if mi < nMols -1:
                countRef[maxC+m+i] = 0
        maxC = maxC + m + spaceC[mi]
    
    #print("shownCs", shownCs)


    atomsinfo = []
    a = 6.63
    b = 20.78
    c = 6.50
    atype = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    beta = 99.5
    diffangle = 99.5 - 90
    # find unit vector of new c axis
    zoffset = np.cos(math.radians(diffangle))
    xoffset = -np.sin(math.radians(diffangle))
    rcoords = np.zeros((9, 3))
    rcoords[:, 0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:, 1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:, 2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    widthy = max(rcoords[:, 1]) - min(rcoords[:, 1])
    widthx = max(rcoords[:, 0]) - min(rcoords[:, 0])
    spacing = widthy / 2

    xnegmin = min(-rcoords[:, 0])
    xposmin = min(rcoords[:, 0])
    diffxmin = xposmin - xnegmin
    symmetries = np.array([[1, 1, 1], [-1, -1, -1], [1, -1 + 0.5, 1 + 0.5], [-1, 1 + 0.5, -1 + 0.5]])
    symmetriesM = np.array([[1, 1, 1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array([[0.5,0,0], [0, 0.5, 1], [0, 0.5, 0], [0.5, 1, 1]])
    #symmetriesNewM = np.array([[1,1,1],[-1,-1,-1],[1,-1,1],[-1,1,-1]])
   
    unitcell = np.zeros((len(symmetries) * len(rcoords), 3))
    ucount = 0
    for s in range(len(symmetries)):
        #print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            unitcell[ucount, 1] = rcoords[ci, 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
            unitcell[ucount, 2] = (rcoords[ci, 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
            unitcell[ucount, 0] = rcoords[ci, 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
            ucount += 1
    # calculate distances between ends:
    # pairs in unit cells
    #print("unitcell before", unitcell)

    #change unit cell order to mesh with order of beads
    first = list(unitcell[0,:])
    #print("first", first)
    #print("second", unitcell[1])
    unitcell[0,:] = list(unitcell[1])
    unitcell[1,:] = first
    # swich another unit cell
    #switch  = list(unitcell[11])
    #unitcell[11] = list(unitcell[12])
    #unitcell[12] = switch

    mid = list(unitcell[18,:])
    unitcell[18,:] = list(unitcell[19,:])
    unitcell[19,:] = mid
    #print("unitcell", unitcell)
    ends1s = [[unitcell[8][0]*a, unitcell[8][1]*b, unitcell[8][2]*c], [unitcell[17][0]*a, unitcell[17][1]*b, unitcell[17][2]*c], [unitcell[26][0]*a, unitcell[26][1]*b, unitcell[26][2]*c+0.5], [unitcell[35][0]*a, unitcell[35][1]*b, unitcell[35][2]*c]]
    ends2s = [[unitcell[10][0]*a, unitcell[10][1]*b+0.8, unitcell[10][2]*c], [unitcell[19][0], unitcell[18][1]*b, unitcell[18][2]*c], [unitcell[28][0]*a, unitcell[28][1]*b+1, unitcell[28][2]*c], [unitcell[0,0]*a, (unitcell[0,1]+1)*b, unitcell[0,2]*c]]

    #ends1s = [[unitcell[8][0]*a, unitcell[8][1]*b, unitcell[8][2]*c], [unitcell[17][0]*a, unitcell[17][1]*b, unitcell[17][2]*c], [unitcell[26][0]*a, unitcell[26][1]*b, unitcell[26][2]*c+0.5], [unitcell[35][0]*a, unitcell[35][1]*b, unitcell[35][2]*c]]
    #ends2s = [[unitcell[10][0]*a, unitcell[10][1]*b+0.8, unitcell[10][2]*c], [unitcell[19][0], unitcell[18][1]*b, unitcell[18][2]*c], [unitcell[28][0]*a, unitcell[28][1]*b+1, unitcell[28][2]*c], [unitcell[0,0]*a, (unitcell[0,1]+1)*b, unitcell[0,2]*c]]
    midpoints = []

    
    for i in range(len(ends1s)):
        # print("end1", ends1s[i], "end2", ends2s[i])
        xdist = (ends1s[i][0] - ends2s[i][0]) * a
        ydist = (ends1s[i][1] - ends2s[i][1]) * b
        zdist = (ends1s[i][2] - ends2s[i][2]) * c
        dist = (xdist ** 2 + ydist ** 2 + zdist ** 2) ** 0.5
    for i in range(len(ends1s)):
        if i % 2 == 0:
            midpoints.append([(ends1s[i][0] + ends2s[i][0]) / 2, (ends1s[i][1] + ends2s[i][1]) / 2,(ends1s[i][2] + ends2s[i][2]) / 2 + c/2])
        else:
            midpoints.append([(ends1s[i][0] + ends2s[i][0]) / 2, (ends1s[i][1] + ends2s[i][1]) / 2,  (ends1s[i][2] + ends2s[i][2]) / 2 - c/2])

    # connection 1
    # parabola in z y
    # linear in x
    Aparabola1 = []
    Hparabola1 = []
    Kparabola1 = []
    Aparabola2 = []
    Hparabola2 = []
    Kparabola2 = []
    Aline = []
    Bline = []

    for i in range(len(ends1s)):
        #print(ends1s[i][1], midpoints[i][1], ends2s[i][1], ends1s[i][2], midpoints[i][2], ends2s[i][2])
        ph1 = midpoints[i][1]
        pk1 = midpoints[i][2]
        pa1 = (ends1s[i][2] - pk1) / ((ends1s[i][1] - ph1) ** 2)
        Aparabola1.append(pa1)
        Hparabola1.append(ph1)
        Kparabola1.append(pk1)
        ph2 = midpoints[i][1]
        pk2 = midpoints[i][2]
        pa2 = (ends2s[i][2] - pk2) / ((ends2s[i][1] - ph2) ** 2)
        Aparabola2.append(pa2)
        Hparabola2.append(ph2)
        Kparabola2.append(pk2)
        la = (ends2s[i][1] - ends1s[i][1]) / (ends2s[i][0] - ends1s[i][0])
        lb = ends2s[i][1] - la * ends2s[i][0]
        # [la,lb] = np.polyfit([ends1s[i][1], ends2s[i][1]],[ends1s[i][0], ends2s[i][0]],1)
        Aline.append(la)
        Bline.append(lb)

    connections1 = []
    connections2 = []
    connections3 = []
    connections4 = []
    connectionsR1 = []
    connectionsR2 = []
    connectionsR3 = []
    connectionsR4 = []
    npoints = [8, 6, 8, 6]
    spacing = [1.4, 1.4, 1.4, 1.4, 1.54]
    
    for i in range(4):
        conx, cony, conz, sidex, sidey, sidez = linear2parabolaConnect(Aline[i], Bline[i], Aparabola1[i], Hparabola1[i], Kparabola1[i], Aparabola2[i], Hparabola2[i], Kparabola2[i], ends1s[i], ends2s[i], midpoints[i], spacing[i])
        sidec = 0
        sidecR = len(sidex) -1
        for j in range(len(conx)):
            vars()["connections" + str(i + 1)].append([conx[j], cony[j], conz[j]])
            if j % 2 == 0:
                vars()["connections" + str(i + 1)].append([sidex[sidec], sidey[sidec], sidez[sidec]])
                sidec += 1
            jR = len(conx) -1 - j
            vars()["connectionsR" + str(i+1)].append([conx[jR], cony[jR], conz[jR]])
            if jR %2 == 0:
                vars()["connectionsR" + str(i+1)].append([sidex[sidecR], sidey[sidecR], sidez[sidecR]])
                sidecR = sidecR - 1
        #for j in range(len(conx)//2, len(conx)):
            #if j % 2 == 0:
                #vars()["connections" + str(i+1)].append([sidex[sidec], sidey[sidec], sidez[sidec]])
                #sidec +=1

            #vars()["connections" + str(i+1)].append([conx[j], cony[j], conz[j]])

    pend1 = list(ends1s[3])
    pend2 = [ends1s[3][0]+a, ends1s[3][1], ends1s[3][2]]

    
    pla = (pend2[1] - pend1[1]) / (pend2[0] - pend1[0])
    plb = pend2[1] - pla * pend2[0]
    pmid = [(pend1[0] + pend2[0]) / 2, (pend1[1] + pend2[1]) / 2, (pend1[2] + pend2[2]) / 2 - c/2]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]

    penda1 = (pend1[2] - pendk1) / ((pend1[0] - pendh1) ** 2)
    penda2 = (pend2[2] - pendk2) / ((pend2[0] - pendh2) ** 2)
    pconx, pcony, pconz, psidex, psidey, psidez = linear2parabolaXZConnect(pla, plb, penda1, pendh1, pendk1, penda2,
                                                                           pendh2, pendk2, pend1, pend2, pmid, 1.6)


    p2end1 = [unitcell[0][0]*a , unitcell[0][1]*b, unitcell[0][2]*c]
    p2end2 = [unitcell[0][0]*a + a, p2end1[1]-0.5, p2end1[2]]
    pla = (p2end2[1] - p2end1[1]) / (p2end2[0] - p2end1[0])
    plb = p2end2[1] - pla * p2end2[0]
    pmid = [(p2end1[0] + p2end2[0]) / 2, (p2end1[1] + p2end2[1]) / 2, (p2end1[2] + p2end2[2]) / 2 - c/2]
    pendh1 = pmid[0]
    pendh2 = pmid[0]
    pendk1 = pmid[2]
    pendk2 = pmid[2]

    penda1 = (p2end1[2] - pendk1) / ((p2end1[0] - pendh1) ** 2)
    penda2 = (p2end2[2] - pendk2) / ((p2end2[0] - pendh2) ** 2)
    p2conx, p2cony, p2conz, p2sidex, p2sidey, p2sidez = linear2parabolaXZConnect2(pla, plb, penda1, pendh1, pendk1,
                                                                                 penda2, pendh2, pendk2, p2end1, p2end2,
                                                                                 pmid, 1.6)
    pconnections1 = []
    cside = 0
    for i in range(len(pconx)):
        pconnections1.append([pconx[i], pcony[i], pconz[i]])
        if i % 2 == 0:
            pconnections1.append([psidex[cside], psidey[cside], psidez[cside]])
            cside += 1
    #for i in range(len(pconx)//2, len(pconx)):
        #if i % 2 == 0:
            #pconnections1.append([psidex[cside], psidey[cside], psidez[cside]])
            #cside += 1
        #pconnections1.append([pconx[i], pcony[i], pconz[i]])

    pconnections2 = []
    cside = 0
    for i in range(len(p2conx)):
        pconnections2.append([p2conx[i], p2cony[i], p2conz[i]])
        if i % 2 == 1:
            pconnections2.append([p2sidex[cside], p2sidey[cside], p2sidez[cside]])
            cside += 1
    #for i in range(len(p2conx)//2, len(p2conx)):
        #if i % 2 == 0:
            #pconnections2.append([p2sidex[cside], p2sidey[cside], p2sidez[cside]])
            #cside += 1
        #pconnections2.append([p2conx[i], p2cony[i], p2conz[i]])

    count = 0
    acount = 1
    mol = 1
    moli = 0
    madded = 0
    minspace = 5
    spacecount = minspace
    alldone = False
    nc = 0
    rcoords = np.zeros((9, 3))
    rcoords[:, 0] = [float(i) for i in "-0.0727 -0.0765 -0.1021 -0.3087 -0.1146 -0.1044 0.2775 0.0872 0.1026".split()]
    rcoords[:, 1] = [float(i) for i in "0.2291 0.1592 0.1602 0.0589 0.0928 0.0854 0.0797 0.1156 0.1221".split()]
    rcoords[:, 2] = [float(i) for i in "0.2004 0.2788 0.5098 0.4941 0.6057 0.8428 0.9260 0.9730 1.2109".split()]
    widthy = max(rcoords[:, 1]) - min(rcoords[:, 1])
    widthx = max(rcoords[:, 0]) - min(rcoords[:, 0])
    spacing = widthy / 2

    xnegmin = min(-rcoords[:, 0])
    xposmin = min(rcoords[:, 0])
    diffxmin = xposmin - xnegmin
    symmetries = np.array([[1, 1, 1], [-1, -1, -1], [1, -1 + 0.5, 1 + 0.5], [-1, 1 + 0.5, -1 + 0.5]])
    symmetriesM = np.array([[1, 1, 1], [-1, -1, -1], [1, 1, 1], [-1, -1, -1]])
    symmetriesA = np.array([[0.5,0,0], [0, 0.5, 1], [0, 0.5, 0], [0.5, 1, 1]]) 
    
    unitcell = np.zeros((len(symmetries) * len(rcoords), 3))
    ucount = 0
    sortedInds0 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds1 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds2 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    sortedInds3 = [1, 0, 2, 4, 3, 5, 7, 6, 8]
    for s in range(len(symmetries)):
        #print("symmetries", symmetries[s])
        for ci in range(len(rcoords)):
            if s == 0:
                unitcell[ucount, 1] = rcoords[sortedInds0[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds0[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds0[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1

            elif s ==1:
                unitcell[ucount, 1] = rcoords[sortedInds1[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds1[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds1[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
            elif s == 2:
                unitcell[ucount, 1] = rcoords[sortedInds2[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds2[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds2[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1
            elif s == 3:
                unitcell[ucount, 1] = rcoords[sortedInds3[ci], 1] * symmetriesM[s, 1] + symmetriesA[s, 1]
                unitcell[ucount, 2] = (rcoords[sortedInds3[ci], 2] * symmetriesM[s, 2] + symmetriesA[s, 2])
                unitcell[ucount, 0] = rcoords[sortedInds3[ci], 0] * symmetriesM[s, 0] + symmetriesA[s, 0]
                ucount += 1

    


    #Reverse
    unitcellReverse = np.copy(unitcell)
    started = True
    for i in range(36):
        unitcellReverse[i,:] = list(unitcell[35-i,:])
    sortedRInds = [0, 2, 1, 3, 5, 4, 6, 8, 7]
    unitcellRcopy = np.copy(unitcellReverse)
    for s in range(len(symmetries)):
        for ci in range(len(rcoords)): 
            unitcellReverse[s*9 + ci,:] = list(unitcellRcopy[s*9 + sortedRInds[ci],:])

    """
    for i in range(9):
        unitcellReverse[i,:] = list(unitcell[i+27,:])
        unitcellReverse[9+i,:] = list(unitcell[i+18,:])
        unitcellReverse[18+i,:] = list(unitcell[i+9,:])
        unitcellReverse[27+i,:] = list(unitcell[i,:])
    """
    srclist = [2, 1, 0, 3]
    
    avgxS0 = np.mean(unitcell[0:9,0])
    avgyS0 = np.mean(unitcell[0:9,1])
    avgzS0 = np.mean(unitcell[0:9,2])
    
    avgxS1 = np.mean(unitcell[9:18,0])
    avgyS1 = np.mean(unitcell[9:18,1])
    avgzS1 = np.mean(unitcell[9:18,2])
    
    avgxS2 = np.mean(unitcell[18:27,0])
    avgyS2 = np.mean(unitcell[18:27,1])
    avgzS2 = np.mean(unitcell[18:27,2])
    
    avgxS3 = np.mean(unitcell[27:36,0])
    avgyS3 = np.mean(unitcell[27:36,1])
    avgzS3 = np.mean(unitcell[27:36,2])

    dist01 = ((avgxS0 - avgxS1)**2 + (avgyS0 - avgyS1)**2 + (avgzS0 - avgzS1)**2)**0.5
    dist12 = ((avgxS2 - avgxS1)**2 + (avgyS2 - avgyS1)**2 + (avgzS2 - avgzS1)**2)**0.5
    dist23 = ((avgxS2 - avgxS3)**2 + (avgyS2 - avgyS3)**2 + (avgzS2 - avgzS3)**2)**0.5

    print("dist01", dist01, "dist12", dist12, "dist23", dist23)
    





    for i in range(Rx):
        basex = i*a
        if i%2 == 0:
            Rylist = list(range(Ry))
            reverse = False
            unitcellnew = np.copy(unitcell)
        else:
            Rylist = list(range(Ry-1, -1, -1))
            reverse = True
            unitcellnew = np.copy(unitcellReverse)
            #Rylist = list(range(Ry))
        for j in Rylist:
            basey = j*b
            for s in range(4):
                if (s%2 == 0):
                    Rzlist = list(range(Rz))
                    
                else:
                    Rzlist = list(range(Rz-1, -1,-1))
                    
                for k in Rzlist:
                    basez = k*c
                    
                    
                    for ci in range(9):
                            if not alldone and madded == 0 and spacecount >= minspace:
                                if not reverse:
                                    if s == 0 and ci == 0:
                                        started = True
                                    elif s == 1 and ci == 0:
                                        started = True
                                    elif s == 2 and ci == 0:
                                        started = True
                                    elif s == 3 and ci == 0:
                                        started == True
                                else:
                                    if s== 0 and ci == 0:
                                        started = True
                                    elif s == 1 and ci == 0:
                                        started = True
                                    elif s == 2 and ci == 0:
                                        started = True
                                    elif s == 3 and ci == 0:
                                        started = True

                            
                            #x = basex + (symmetries[s][0]*a + rcoords[ci][0]*a) + xoffset*(symmetries[s][2]*c + rcoords[ci][2]*c)*a
                            x = basex + unitcellnew[s*9 + ci,0]*a 
                            #print("x", x,"basex", basex, "symmetries[s][0]*a", symmetries[s][0]*a, "symmetries[s][0]", symmetries[s][0])
                            y = basey + unitcellnew[s*9 + ci,1]*b
                            #print("y", y, "basey", basey, "symmetries[s][1]*b", symmetries[s][1]*b, "symmetries[s][1]", symmetries[s][1])
                            #z = basez + (symmetries[s][2]*c + rcoords[ci][2]*c)*zoffset 
                            z = basez + unitcellnew[s*9 + ci,2]*c
                            tx = x + z*xoffset
                            ty = y
                            tz = z*zoffset
                            if not alldone and madded < moleculeClist[moli] and spacecount >= minspace and started: 
                                atomsinfo.append([acount, 1+s, 1, tx, ty, tz])
                                acount = acount + 1
                                madded  = madded + 1
                                #print("madded", madded)
                                #print("spacecount", spacecount)
                            elif (madded == 0 and spacecount < minspace):
                                if spacecount == 0:
                                    if ci > 0:
                                        spacecount = spacecount + 1
                                else:
                                    spacecount = spacecount + 1
                                #print("tx", tx, "ty", ty, "tz", tz)
                               
                            if not alldone and madded == moleculeClist[moli]:
                                moli = moli + 1
                                print("moli", moli)
                                madded = 0
                                spacecount = 0
                                if moli == len(moleculeClist):
                                    alldone = True
                                minspace = 5 +nc%2
                                started = False
                                #minspace = 5 
                                #print("minspace", minspace)
                            count = count + 1
                            

                    if reverse == False:
                        if k == Rzlist[-1] and s % 2 == 0:
                            nc = nc +1
                            for ii in range(len(vars()["connections" + str(s + 1)])):
                                x = basex + vars()["connections" + str(s + 1)][ii][0]
                                y = basey + vars()["connections" + str(s + 1)][ii][1]
                                z = basez + vars()["connections" + str(s + 1)][ii][2]
                                tx = x + z * xoffset
                                ty = y
                                tz = z * zoffset
                                if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    acount = acount + 1
                                    madded = madded +1
                                count = count + 1
                        elif k == Rzlist[-1] and s % 2 == 1 and (s != 3 or j < Ry - 1):
                            for ii in range(len(vars()["connections" + str(s + 1)])):
                                x = basex + vars()["connections" + str(s + 1)][ii][0]
                                y = basey + vars()["connections" + str(s + 1)][ii][1]
                                z = vars()["connections" + str(s + 1)][ii][2]
                                tx = x + z * xoffset
                                ty = y
                                tz = z * zoffset
                                if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    acount = acount + 1
                                    madded = madded + 1
                                count = count + 1
                        elif k == Rzlist[-1] and s == 3 and j == Ry - 1 and i % 2 == 0:
                            nc = nc +1
                            for ii in range(len(pconnections1)):
                                x = basex + pconnections1[ii][0]
                                y = basey + pconnections1[ii][1]
                                z = pconnections1[ii][2]
                                tx = x + z * xoffset
                                ty = y
                                tz = z * zoffset
                                if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    acount = acount + 1
                                    madded = madded + 1
                                count = count + 1

                        elif k == Rzlist[-1] and s == 3 and j == Ry - 1 and i % 2 == 1:

                            for ii in range(len(pconnections2)):
                                x = basex + pconnections2[ii][0]
                                y = pconnections2[ii][1]
                                z = pconnections2[ii][2]
                                tx = x + z * xoffset
                                ty = y
                                tz = z * zoffset
                                if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    acount = acount + 1
                                    madded = madded + 1
                                count = count + 1
                        if not alldone and madded == moleculeClist[moli]:
                            madded = 0

                            moli = moli + 1
                            print("moli", moli)
                            if moli == len(moleculeClist):
                                alldone = True
                            #print("mol", mol)
                        #print("acount", acount)
                    
                    else: 
                        sr = srclist[s]
                        #print("sr", sr)
                        if k == Rzlist[-1] and s==3 and j == Rylist[-1] and i%2 ==0:
                            for ii in range(len(pconnections1)):
                                x = basex + pconnections1[ii][0]
                                y = basey + pconnections1[ii][1]
                                z = pconnections1[ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    acount = acount + 1
                                    madded = madded + 1
                                count = count + 1

                        elif k == Rzlist[-1] and s==3 and j == Rylist[-1] and i%2 ==1:
                            nc = nc + 1
                            for ii in range(len(pconnections2)):
                                x = basex  + pconnections2[ii][0]
                                y = pconnections2[ii][1]
                                z = pconnections2[ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    acount = acount + 1
                                    madded = madded + 1
                                count = count + 1
                        elif k == Rzlist[-1] and sr%2 == 0:
                            nc = nc + 1
                            for ii in range(len(vars()["connectionsR"+str(sr+1)])):
                                x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                y = basey + vars()["connectionsR" + str(sr+1)][ii][1]
                                z = basez + vars()["connectionsR" + str(sr+1)][ii][2]
                                tx = x + z*xoffset
                                ty = y 
                                tz = z*zoffset
                                if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    acount = acount + 1
                                    madded = madded + 1
                                count = count + 1
                        elif k == Rzlist[-1] and sr%2 ==1 and (sr != 3 or j != Rylist[-1]):
                            
                            for ii in range(len(vars()["connectionsR"+str(sr+1)])):
                                if sr != 3:
                                    x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                    y = basey + vars()["connectionsR" + str(sr+1)][ii][1]
                                    z = vars()["connectionsR" + str(sr+1)][ii][2]
                                elif sr == 3:
                                    x = basex + vars()["connectionsR" + str(sr+1)][ii][0]
                                    y = basey + vars()["connectionsR" + str(sr+1)][ii][1] - b
                                    z = vars()["connectionsR" + str(sr+1)][ii][2]
                                tx = x + z*xoffset
                                ty = y
                                tz = z*zoffset
                                if not alldone and madded < moleculeClist[moli] and madded > 0:
                                    atomsinfo.append([acount, 5+s, 1, tx, ty, tz])
                                    acount = acount + 1
                                    madded = madded + 1
                                count = count + 1
                        if not alldone and madded == moleculeClist[moli]:
                            madded = 0
                            moli = moli + 1
                            print("moli", moli)
                            if moli == len(moleculeClist):
                                alldone = True
                            spacecount = 0
                            minspace = 5 +nc%2
                            started = False

    if count < len(countRef) or acount != nBeads+1:
        print("ERROR: NOT ALL WERE SHOWN")
        print("count", count, "len(shownCs)", len(shownCs))
        print("acount", acount, "nBeads+1", nBeads+1)
    atomsinfo = np.array(atomsinfo)
    writexyz(xyzfile, atomsinfo)
    bspace =2
    xlo = min(atomsinfo[:, 3]) - bspace
    xhi = max(atomsinfo[:, 3]) + bspace
    ylo = min(atomsinfo[:, 4]) - bspace
    yhi = max(atomsinfo[:, 4]) + bspace
    zlo = min(atomsinfo[:, 5]) - bspace
    zhi = max(atomsinfo[:, 5]) + bspace
    return atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi

def makeMolList(Nx, Carbons):
    molList = []
    for i,x in enumerate(Carbons):
        reps = Nx[i]
        for ii in range(reps):
            molList.append(x)

    return molList 

def checkDistancesBefore(atomsinfo):
    N = len(atomsinfo)
    for i in range(2, N):
        x1 = atomsinfo[i-1][3]
        y1 = atomsinfo[i-1][4]
        z1 = atomsinfo[i-1][5]
        x2 = atomsinfo[i][3]
        y2 = atomsinfo[i][4]
        z2 = atomsinfo[i][5] 
        x3 = atomsinfo[i-2][3]
        y3 = atomsinfo[i-2][4]
        z3 = atomsinfo[i-2][5]
        dist1 = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
        dist2 = ((x2 - x3)**2 + (y2 - y3)**2 + (z2 - z3)**2)**0.5
        dist3 = ((x1 - x3)**2 + (y1 - y3)**2 + (z1 - z3)**2)**0.5
        if dist1 > 2:
            print("i", i, "dist1", dist1, "dist2", dist2, "dist3", dist3, "atom1", i, "atom2", i+1, "atom3", i-1)

def checkDistances(atomsinfo):
    N = len(atomsinfo)
    nbonds = 0
    for i in range(0,N):
        for j in range(i+1, min(i+5,N)):
            x1 = atomsinfo[i][3]
            y1 = atomsinfo[i][4]
            z1 = atomsinfo[i][5]
            x2 = atomsinfo[j][3]
            y2 = atomsinfo[j][4]
            z2 = atomsinfo[j][5]

            dist = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
            threshold = 1.863
            small = 1
            if dist < threshold:
                nbonds = nbonds + 1
                #print("nbonds", nbonds, "dist", vars()["dist" + str(ii+1)])
            if dist < small:
                print("too close", "i", i, "j", j, "dist", dist)
            
            
    print("nbonds", nbonds)



def checkNumberofPosBondsandAdd(atomsinfo):
    threshold = 1.863
    N = len(atomsinfo)
    print("natoms", N)
    nbonds = 0
    atoms1 = []
    atoms2 = []
    for i in range(0,N):
        for j in range(i+1, min(i+5, N)):
            x1 = atomsinfo[i,3]
            x2 = atomsinfo[j,3]
            y1 = atomsinfo[i,4]
            y2 = atomsinfo[j,4]
            z1 = atomsinfo[i,5]
            z2 = atomsinfo[j,5]
            dist = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
            if dist < threshold:
                nbonds = nbonds +1
                #print("nbonds", nbonds, "i", i, "j", j)
                atoms1.append(i+1)
                atoms2.append(j+1)

    print("added bonds", nbonds)
    bondsinfo = np.zeros((nbonds, 4))
    for i in range(nbonds):
        bondsinfo[i,0] = i+1
        bondsinfo[i,1] = 1
        bondsinfo[i,2] = atoms1[i]
        bondsinfo[i,3] = atoms2[i]
    return bondsinfo
def checkForTriangles(bondsinfo, natoms):
    print("looking for triangles")
    atomsBonds = [[] for i in range(natoms)]
    #print("atomsBonds", atomsBonds)
    for i in range(len(bondsinfo)):
        atom1 = int(bondsinfo[i,2])
        atom2 = int(bondsinfo[i,3])
        atomsBonds[atom1-1].append(atom2)
        atomsBonds[atom2-1].append(atom1)
    for i in range(natoms):
        bonds = atomsBonds[i]
        #print("bonds", bonds)
        for ii in range(len(bonds)-1): 
        
            batom1  = int(bonds[ii]) 
            batom2 = int(bonds[ii+1])
            if batom2 in atomsBonds[batom1-1] or batom1 in atomsBonds[batom2-1]:
                    print("there is triangle", "atom1", i+1, "atom2", batom1, "atom3", batom2)

        


def main():
    
    setupOldInfiniteSystem("test.xyz", 1, 1, 1)
    
    """
    
    atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = writeMoleculesinCrystalNoCons("MolinCryst.xyz", [12], 1, 1, 1)
    checkDistances(atomsinfo)
    msc.writelammpsdatajustatoms("MolinCryst.data", [xlo, xhi, ylo, yhi, zlo, zhi], [15], len(atomsinfo), atomsinfo)
    bondsinfo = checkNumberofPosBondsandAdd(atomsinfo)
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, oldbondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("MolinCryst.data")
    msc.writelammpsdataonebondtype("MolinCrystCustomBonds.data", boxcoords, masstypes, atoms, len(bondsinfo), angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    atomsinfo = readlammpsbondsPPctypes("MolinCrystCustomBonds.data", "MolinCrystCtype.data")
    msc.writelammpsdataonebondtype("MolinCrystCtypebonds.data", boxcoords, [15.035, 14.027, 13.019], atoms, len(bondsinfo), angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    """
    """
    p613MolList = makeMolList([9, 3, 1, 8, 13, 3, 1, 1], [926, 1436, 1946, 2456, 2966, 3476, 3986, 41216])
    print("p613MolList",p613MolList)
    atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = writeMoleculesinCrystal("p613MolinCryst.xyz", p613MolList, 20, 20, 20)
    #atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = writeMoleculesinCrystal("MolinCryst.xyz", p613MolList, 20, 20, 20)
    #msc.writelammpsdatajustatoms("MolinCryst.data",[xlo,xhi,ylo,yhi,zlo,zhi], [15], len(atomsinfo), atomsinfo)
    #checkDistances(atomsinfo)
    msc.writelammpsdatajustatoms("p613MolinCryst.data", [xlo, xhi, ylo, yhi, zlo, zhi], [15], len(atomsinfo), atomsinfo)
    bondsinfo = checkNumberofPosBondsandAdd(atomsinfo)
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, oldbondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("p613MolinCryst.data")
    msc.writelammpsdataonebondtype("p613MolinCrystCustomBonds.data", boxcoords, masstypes, atoms, len(bondsinfo), angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    
    
    atomsinfo = readlammpsbondsPPctypes("p613MolinCrystCustomBonds.data", "p613MolinCrystCtype.data")
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, oldatomsinfo, bondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("p613MolinCrystCustomBonds.data")
    
    msc.writelammpsdataonebondtype("p613MolinCrystCtypebonds.data", boxcoords, [15.035, 14.027, 13.019], atoms, bonds, angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    """
    """ 
    molList = makeMolList([9,7,1], [926,1436,3*1+2]) 
    print("molList", molList, "molecules", len(molList), "correct bonds", sum(molList) - len(molList))
    atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = writeMoleculesinCrystal("MolinCryst.xyz", molList, 20,20,20)
    checkDistances(atomsinfo)
    msc.writelammpsdatajustatoms("MolinCryst.data", [xlo,xhi,ylo,yhi,zlo,zhi], [15], len(atomsinfo), atomsinfo)
    bondsinfo = checkNumberofPosBondsandAdd(atomsinfo)
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, oldbondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("MolinCryst.data")
    msc.writelammpsdataonebondtype("MolinCrystCustomBonds.data", boxcoords, masstypes, atoms, len(bondsinfo), angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    checkForTriangles(bondsinfo, len(atomsinfo))
    checkDistances(atomsinfo)
    #boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("MolinCrystCustomBonds.data")
    #checkForTriangles(bondsinfo, len(atomsinfo))
    """
    """
    p813MolList = makeMolList([9,7,6,12,3,6,1,1,1,1], [926,1436,1946,2456,2966,3476,3986,11126,12656,32546])
    print("p813MolList",p813MolList)
    atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = writeMoleculesinCrystal("p813MolinCryst.xyz", p813MolList, 20, 20, 20)
    msc.writelammpsdatajustatoms("p813MolinCryst.data", [xlo, xhi, ylo, yhi, zlo, zhi], [15], len(atomsinfo), atomsinfo)
    bondsinfo = checkNumberofPosBondsandAdd(atomsinfo)
    checkDistances(atomsinfo)
    checkForTriangles(bondsinfo, len(atomsinfo))
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, oldbondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("p813MolinCryst.data")
    msc.writelammpsdataonebondtype("p813MolinCrystCustomBonds.data", boxcoords, masstypes, atoms, len(bondsinfo), angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    atomsinfo = readlammpsbondsPPctypes("p813MolinCrystCustomBonds.data", "p813MolinCrystCtype.data")
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, oldatomsinfo, bondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("p813MolinCrystCustomBonds.data")
    msc.writelammpsdataonebondtype("p813MolinCrystCtypebonds.data", boxcoords, [15.035, 14.027, 13.019], atoms, bonds, angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    
      
    
    pMPMolList = makeMolList([3, 4, 10, 11, 8, 4, 1, 1, 1, 1, 1], [926, 1526, 2126, 2726, 3326, 3926, 4526, 6326, 6926, 8126, 54326])
    print("pMPMolList",pMPMolList)
    atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = writeMoleculesinCrystal("pMPMolinCryst.xyz", pMPMolList, 20, 20, 20)
    msc.writelammpsdatajustatoms("pMPMolinCryst.data", [xlo, xhi, ylo, yhi, zlo, zhi], [15], len(atomsinfo), atomsinfo)
    bondsinfo = checkNumberofPosBondsandAdd(atomsinfo)
    checkDistances(atomsinfo)
    checkForTriangles(bondsinfo, len(atomsinfo))
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, oldbondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("pMPMolinCryst.data")
    msc.writelammpsdataonebondtype("pMPMolinCrystCustomBonds.data", boxcoords, masstypes, atoms, len(bondsinfo), angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    atomsinfo = readlammpsbondsPPctypes("pMPMolinCrystCustomBonds.data", "pMPMolinCrystCtype.data")
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, oldatomsinfo, bondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("pMPMolinCrystCustomBonds.data")
    msc.writelammpsdataonebondtype("pMPMolinCrystCtypebonds.data", boxcoords, [15.035, 14.027, 13.019], atoms, bonds, angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    """
    #readlammpsbondsPPctypes("pMPMolinCrystbonds.data", "pMPMolinCrystCtype.data")

    
    #atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = writeMoleculesinCrystal("MolinCryst.xyz",[3*10 + 2, 3*10 + 2, 3*10 + 2, 10*3 + 2], 5, 5, 5)
    #msc.writelammpsdatajustatoms("MolinCryst.data", [xlo, xhi, ylo, yhi, zlo, zhi], [15], len(atomsinfo), atomsinfo)
    #readlammpsbondsPPctypes("MolinCrystbonds.data", "MolinCrystCtype.data")
    #readlammpsbondsPPctypes("TrialInfa1iPPbonds.data", "TrialInfa1iPPC1type.data")
    """
    atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = setupInfiniteSystem("TrialInfa1iPP.xyz", 20, 20, 20)
    msc.writelammpsdatajustatoms("TrialInfa1iPP.data", [xlo, xhi, ylo, yhi, zlo, zhi], [15], len(atomsinfo), atomsinfo)
    checkDistances(atomsinfo)
    boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, oldbondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("TrialInfa1iPP.data")
    bondsinfo = checkNumberofPosBondsandAdd(atomsinfo) 
    msc.writelammpsdataonebondtype("TrialInfa1iPPbonds.data", boxcoords, masstypes, atoms, len(bondsinfo), angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo) 
    """
    #msc.writelammpsdatajustatoms("TrialInfa1iPP.data",[xlo,xhi,ylo,yhi,zlo,zhi], [15], len(atomsinfo), atomsinfo)
    #bondsinfo = checkNumberofPosBondsandAdd(atomsinfo)
    #boxcoords, masstypes, atoms, bonds, angles, dihedrals, atomsinfo, oldbondsinfo, anglesinfo, dihedralsinfo = msc.readlammpsdata("TrialInfa1iPP.data")
    #msc.writelammpsdataonebondtype("TrialInfa1iPPCustomBonds.data", boxcoords, masstypes, atoms, len(bondsinfo), angles, dihedrals, atomsinfo, bondsinfo, anglesinfo, dihedralsinfo)
    #checkDistances(atomsinfo)
    #atomsinfo, xlo, xhi, ylo, yhi, zlo, zhi = readxyz("custompp-iso.xyz")
    #msc.writelammpsdatajustatoms("custompp-iso.data", [xlo, xhi, ylo, yhi, zlo, zhi], [14, 1], len(atomsinfo), atomsinfo)
main()
