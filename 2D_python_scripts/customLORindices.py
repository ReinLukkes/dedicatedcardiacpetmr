## Generate lookup table for the LOR-angle- and voxel-indices of the geometry
import numpy as np


def customLORindices(geom, m_nbRsectorsPerFullRing, m_nbModulesAxlFullRing, m_nbModulesTrs, m_nbCrystalsTrs, m_nbCrystalsAxl):

    m_nbRingsFullRing = m_nbModulesAxlFullRing * m_nbCrystalsAxl
    m_nbCrystalsPerFullRing = m_nbRsectorsPerFullRing * m_nbModulesTrs * m_nbCrystalsTrs

    # Minimum sector difference: make sure only LORs involving crystals that are at least this amount of rsectors away from each other are used
    m_minSectorDifference = 0
    m_minCrystalDifference = m_minSectorDifference * m_nbModulesTrs * m_nbCrystalsTrs

    # Sinogram dimensions (always defined for full ring system)
    m_nbSinogramBins = m_nbCrystalsPerFullRing - 2 * (m_minCrystalDifference - 1) - 1
    
    # m_nbSinogramBins is half the acceptation angle (angle / 2) but multiplied by a factor 2 because the LORs of
    # angles phi and phi+1 are both mapped to the same sinogram row (interleaved, to increase sampling)
    # see Bailey 2005, PET Basic Sciences, Figure 3.5
    m_nbSinogramAngles = m_nbCrystalsPerFullRing / 2 # only need to cover 180 degrees (other 180 are the same LORs)
    m_nbSinograms = m_nbRingsFullRing * m_nbRingsFullRing
    
    # determine transaxial ID of crystal (in its own ring) relative to the crystal at the center of the first rsector
    # this requires shifting all IDs by half the rsector size
    # the CASToR default is that the first rsector is at the top of the scanner (positive y-axis pointing towards the ceiling)
    # which implies that the top crystal's ID is 0 and all LORs having phi=0 are aligned with the positive y-axis
    distanceCrystalId0toFirstRsectorCenter = (m_nbModulesTrs * m_nbCrystalsTrs) / 2
    
    LOR_indices = np.zeros((m_nbCrystalsPerFullRing, m_nbCrystalsPerFullRing, 2)) - 1
    #for d1 = 1:m_nbCrystalsPerFullRing:
    for d1 in range(m_nbCrystalsPerFullRing):
        castorFullRingCrystalID1 = d1
        id1 = (castorFullRingCrystalID1 % m_nbCrystalsPerFullRing) - distanceCrystalId0toFirstRsectorCenter
        #for d2 = 1:m_nbCrystalsPerFullRing:
        for d2 in range(m_nbCrystalsPerFullRing):
            castorFullRingCrystalID2 = d2
            id2 = (castorFullRingCrystalID2 % m_nbCrystalsPerFullRing) - distanceCrystalId0toFirstRsectorCenter
            
            if (id1 < 0):
                id1 = id1 + m_nbCrystalsPerFullRing
            if (id2 < 0):
                id2 = id2 + m_nbCrystalsPerFullRing
            
            A = 0 
            B = 0
            if (id1 < id2):
                A = id1
                B = id2
                #ringIdA = castorFullRingCrystalID1 / m_nbCrystalsPerFullRing
                #ringIdB = castorFullRingCrystalID2 / m_nbCrystalsPerFullRing
            else:
                A = id2
                B = id1
                #ringIdA = castorFullRingCrystalID2 / m_nbCrystalsPerFullRing
                #ringIdB = castorFullRingCrystalID1 / m_nbCrystalsPerFullRing
            
            r = 0 
            phi = 0
            if (B - A < m_minCrystalDifference):
                continue
            else:
                if (A + B >= (3 * m_nbCrystalsPerFullRing) / 2 or A + B < m_nbCrystalsPerFullRing / 2):
                    if (A == B):
                        r = -m_nbCrystalsPerFullRing / 2
                    else:
                        r = ((B - A - 1) / 2) - ((m_nbCrystalsPerFullRing - (B - A + 1)) / 2)
                else:
                    if (A == B):
                        r = m_nbCrystalsPerFullRing / 2
                    else:
                        r = ((m_nbCrystalsPerFullRing - (B - A + 1)) / 2) - ((B - A - 1) / 2)
                
                r = int(r)
                
                if (A + B < m_nbCrystalsPerFullRing / 2):
                    phi = (2 * A + m_nbCrystalsPerFullRing + r) / 2
                else:
                    if (A + B >= (3 * m_nbCrystalsPerFullRing) / 2):
                        phi = (2 * A - m_nbCrystalsPerFullRing + r) / 2
                    else:
                        phi = (2 * A - r) / 2
                
                LOR_indices[d1, d2, 0] = int(phi)
                LOR_indices[d1, d2, 1] = int(r + m_nbSinogramBins / 2)
    
    return LOR_indices