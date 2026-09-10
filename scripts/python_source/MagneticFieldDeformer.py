# Description
#------------
# Simple module to manipulate the (2D) magnetic field using simple
# and elementary transformations. The magnetic is assumed to be 
# provided in structured RZPsi format. 

import numpy as np
import scipy as sp
import copy

def Crop(R, Z, Psi, Rbox, Zbox):
    # Description
    #------------
    # Crop the magnetic field in a rather 'dumb' way by simply removing
    # all R and Z values that are not in the box defined by the minimal
    # and maximal value in Rbox, Zbox. 

    # Copy
    newR = copy.deepcopy(R)
    newZ = copy.deepcopy(Z)
    newPsi = copy.deepcopy(Psi)

    # Determine box
    Rmin = np.min(Rbox)
    Rmax = np.max(Rbox)
    Zmin = np.min(Zbox)
    Zmax = np.max(Zbox)

    # Determine cropping bounds
    keepR = (R >= Rmin) & (R <= Rmax)
    keepZ = (Z >= Zmin) & (Z <= Zmax)
    newR = R[keepR]
    newZ = Z[keepZ]
    newPsi = Psi[np.ix_(keepR, keepZ)]

    # Return 
    return newR, newZ, newPsi

def Translate(R, Z, Psi, dR, dZ):
    # Description
    #------------
    # Translate the magnetic field by dR in the R-direction and dZ in the
    # Z-direction. Only changes coordinates, not psi values

    newR = R + dR 
    newZ = Z + dZ 
    newPsi = Psi 

    return newR, newZ, newPsi 

def MirrorR(R, Z, Psi, R0, changesign): 
    # Description
    #------------
    # Mirror the field along the line described by R = R0. Note: R0 
    # should not be present in R except either at the end or at the 
    # start of R! Otherwise, an error will be thrown. 

    # Initialize
    nR = len(R)
    nZ = len(Z) 
    if changesign:
        mult = -1.0
    else:
        mult = 1.0

    # Checks
    R2L = True # Rotate from right to left?
    if all(R <= R0):
        R2L = False 
    elif all(R >= R0):
        R2L = True 
    else:
        # R0 is contained within R, throw error 
        raise ValueError("MirrorR: R0 is present in R vector, cannot mirror")  
    
    # Construct new R and Psi (Z remains the same)
    if R2L:
        # Check if we should duplicate the first value of R - don't if
        # it is equal to R0
        excludefirst = False 
        if R[0] == R0: 
            excludefirst = True

        # Initialize
        newZ = Z 
        if excludefirst: 
            newnR = nR*2-1
            startindR = 1
            newR = np.zeros(newnR, dtype=float)
            newPsi = np.zeros((newnR, nZ), dtype=float)
            newR[nR-1] = R[0]
            newPsi[nR-1, :] = Psi[0, :]
        else:
            newnR = 2*nR 
            startindR = 0
            newR = np.zeros(newnR, dtype=float)
            newPsi = np.zeros((newnR, nZ), dtype=float)
            

        # Extend
        i = startindR 
        while i < nR:
            # Build indices
            leftindnew = nR-1 - i
            rightindnew = i + nR - startindR

            # Copy old values at 'right' side of domain
            newR[rightindnew] = R[i]
            newPsi[rightindnew, :] = Psi[i, :]

            # Mirror values to 'left' side of domain
            newR[leftindnew] = -R[i] + 2*R0
            newPsi[leftindnew, :] = mult*Psi[i, :]  

            # Increment 
            i = i + 1
 
    else: 
        # Check if we should duplicate the last value of R - don't if
        # it is equal to R0
        excludelast = False 
        if R[nR-1] == R0: 
            excludelast = True

        # Initialize
        newZ = Z 
        if excludelast: 
            newnR = nR*2-1
            endindR = nR-1
            newR = np.zeros(newnR, dtype=float)
            newPsi = np.zeros((newnR, nZ), dtype=float)
            newR[nR-1] = R[nR-1]
            newPsi[nR-1, :] = Psi[nR-1, :]
        else:
            newnR = 2*nR 
            endindR = nR
            newR = np.zeros(newnR, dtype=float)
            newPsi = np.zeros((newnR, nZ), dtype=float)
            

        # Extend
        i = 0
        while i < endindR:
            # Build indices
            leftindnew = i
            rightindnew = i + nR

            # Copy old values at 'left' side of domain
            newR[leftindnew] = R[i]
            newPsi[leftindnew, :] = Psi[i, :]

            # Mirror values to 'right' side of domain
            newR[rightindnew] = -R[i] + 2*R0 
            newPsi[rightindnew, :] = mult*Psi[i, :]  
            
            # Increment 
            i = i + 1

    return newR, newZ, newPsi

def Perturb(R, Z, Psi, delta, diffmeth):
    # Description
    #------------
    # Perturb the magnetic field Psi values by adding an increment 
    # delta in some way defined by diffmeth:
    # - 'Psi'       just add to Psi values itself  
    # - 'R'         add delta*R (this means that you add a delta to dPsi/dR)
    # - 'Z'         add delta*Z (this means that you add a delta to dPsi/dZ)
    # - 'R2'        add 0.5*delta*R**2 (this means that you add a delta to 1/R dPsi/dR)    
    # - 'Z2'        add 0.5*delta*Z**2 (this means that you add a delta to 1/Z dPsi/dZ)      

    # Initialize
    newPsi = np.zeros((len(R), len(Z)), dtype=float) 
    newR = R
    newZ = Z

    # Perturb
    match diffmeth: 

        case 'Psi':

            newPsi = Psi + delta 

        case 'R':

            for i in np.arange(0, len(R), 1):
                newPsi[i, :] = Psi[i, :] + delta*R[i]

        case 'Z':

            for j in np.arange(0, len(Z), 1):
                newPsi[:, j] = Psi[:, j] + delta*Z[j]

        case 'R2':

            for i in np.arange(0, len(R), 1):
                newPsi[i, :] = Psi[i, :] + 0.5*delta*R[i]*R[i]

        case 'Z2':

            for j in np.arange(0, len(Z), 1):
                newPsi[:, j] = Psi[:, j] + 0.5*delta*Z[j]*Z[j]

        case _:
            raise ValueError('Perturb: unknown diffmeth')
            

    # Return 
    return newR, newZ, newPsi 
            

