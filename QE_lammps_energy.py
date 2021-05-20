import numpy as np
import os
import sys

bohr2ang = 0.529177249
ry2ev = 13.6056980659
rad2deg = 57.295779513


class QExpresso:
    def __init__(self, inFile):
        self.inFile = inFile

    def read(self):
        fileHandler = open(self.inFile, "r")
        self.lines = fileHandler.readlines()
        fileHandler.close()

        self.readLattice()
        self.readCellMat()
        self.readCoord()
        self.readEnergy()
        #self.readMagMoment()

    def readLattice(self):

        for line in self.lines:

            if "lattice parameter (alat)  =" in line:
                self.latParam = float(line.split("=")[1].split()[0])

    def readCellMat(self):

        self.cellMat = np.zeros((3, 3), dtype=np.double)

        for lineID, line in enumerate(self.lines):
            if "crystal axes: (cart. coord. in units of alat)" in line:
                line1 = self.lines[lineID + 1]
                line2 = self.lines[lineID + 2]
                line3 = self.lines[lineID + 3]

        rowID = 0
        for line in [line1, line2, line3]:
            ll = line.split("=")[1].replace("(", "").replace(")", "").split()
            vec = [float(Str) for Str in ll]
            self.cellMat[rowID, :] = vec
            rowID += 1

        self.cellMat *= self.latParam
        self.cellMat *= bohr2ang
        #print(self.cellMat)

    def readCoord(self):

        sIDx = None
        for lineID, line in enumerate(self.lines):
            if "site n.     atom                  positions (cryst. coord.)" in line:
                sIDx = lineID + 1

            if sIDx and not line.strip():
                eIDx = lineID
                break

        self.nAtoms = eIDx - sIDx

        self.crystalCoords = np.zeros((self.nAtoms, 3), dtype=np.double)
        self.symbols = []

        iatom = 0
        for line in self.lines[sIDx:eIDx]:

            ll = line.split("=")[1].replace("(", "").replace(")", "").split()
            xi, yi, zi = [float(Str) for Str in ll]

            self.crystalCoords[iatom, 0] = xi
            self.crystalCoords[iatom, 1] = yi
            self.crystalCoords[iatom, 2] = zi

            # get symbols
            isymbol = line.split()[1]
            for digit in range(10):
                digit = str(digit)
                isymbol = isymbol.replace(digit, "")
            self.symbols.append(isymbol)

            iatom += 1

        # type assiging
        unique_symbols = list(set(self.symbols))
        unique_symbols.sort()

        self.types = []

        for symbol in self.symbols:

            typeID = unique_symbols.index(symbol) + 1
            self.types.append(typeID)

    def readEnergy(self):

        for line in self.lines:

            if "!    total energy" in line:
                self.totEnr = float(line.split("=")[1].split()[0])
                self.totEnr *= ry2ev

#     def readMagMoment(self):
# 
#         sIDx = None
#         for lineID, line in enumerate(self.lines):
#             if "Magnetic moment per site" in line:
#                 sIDx = lineID + 1
# 
#             if sIDx and not line.strip():
#                 eIDx = lineID
#                 break
# 
#         if (eIDx - sIDx) != self.nAtoms:
#             raise RuntimeError("Parsing failed during Magnetic moment read")
# 
#         self.magMoment = []
# 
#         for line in self.lines[sIDx:eIDx]:
#             m = float(line.split("magn:")[1].split()[0])
# 
#             if m >= 0:
#                 m = 1.0
#             else:
#                 m = -1.0
# 
#             self.magMoment.append(m)

    def fixCellMat(self):

        a = self.cellMat[0, :]
        b = self.cellMat[1, :]
        c = self.cellMat[2, :]
        
        la = np.linalg.norm(a)
        lb = np.linalg.norm(b)
        lc = np.linalg.norm(c)
        alpha = self.vec2angle(b, c)
        beta = self.vec2angle(a, c)
        gamma = self.vec2angle(a, b)
        
        self.la = la
        self.lb = lb
        self.lc = lc
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        self.cellMat_fixed = np.zeros((3, 3), dtype=np.double)

        self.cellMat_fixed[0, 0] = la
        self.cellMat_fixed[0, 1] = 0.0
        self.cellMat_fixed[0, 2] = 0.0

        self.cellMat_fixed[1, 0] = lb * np.cos(gamma)
        self.cellMat_fixed[1, 1] = lb * np.sin(gamma)
        self.cellMat_fixed[1, 2] = 0.0

        self.cellMat_fixed[2, 0] = lc * np.cos(beta)
        self.cellMat_fixed[2, 1] = (
            lc * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
        )
        self.cellMat_fixed[2, 2] = (
            lc
            * np.sqrt(
                1
                + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
                - np.cos(alpha) ** 2
                - np.cos(beta) ** 2
                - np.cos(gamma) ** 2
            )
            / np.sin(gamma)
        )
        
        #print(self.cellMat_fixed)
        # convert t ocart coordinates
        self.cartCoords = np.matmul(self.crystalCoords, self.cellMat_fixed)

        # get bounds
        self.xy = self.cellMat_fixed[1, 0]
        self.yz = self.cellMat_fixed[2, 1]
        self.xz = self.cellMat_fixed[2, 0]
        
        self.aa = self.cellMat_fixed[0, 0]
        self.bb = self.cellMat_fixed[1, 1]
        self.cc = self.cellMat_fixed[2, 2]

        self.xlo_bound = 0.0+ min(0.0, self.xy, self.xz, self.xy + self.xz)
        self.xhi_bound = self.aa + max(0.0, self.xy, self.xz, self.xy + self.xz)
        self.ylo_bound = min(0.0, self.yz)
        self.yhi_bound = self.bb + max(0.0, self.yz)
        
        self.zlo_bound = 0.0
        self.zhi_bound = self.cc

    @staticmethod
    def vec2angle(vec1, vec2):

        angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        angle = np.arccos(angle)
        return angle

    def write(self, outFH, nFrames, iFrame):
        outFH.write("ITEM: TIMESTEP energy, energy_weight, force_weight, nsims\n")

        outFH.write("%-5d    %-.16f    1    1   %d\n" % (iFrame, self.totEnr, nFrames))
        outFH.write("ITEM: NUMBER OF ATOMS\n")
        outFH.write("%-10d\n" % self.nAtoms)

        outFH.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
        outFH.write(
            "%22.16f  %22.16f  %22.16f\n" % (self.xlo_bound, self.xhi_bound, self.xy)
        )
        outFH.write(
            "%22.16f  %22.16f  %22.16f\n" % (self.ylo_bound, self.yhi_bound, self.xz)
        )
        outFH.write(
            "%22.16f  %22.16f  %22.16f\n" % (self.zlo_bound, self.zhi_bound, self.yz)
        )

        outFH.write("ITEM: ATOMS id type x y z\n")
        for i in range(self.nAtoms):

            outFH.write(
                "%-5d  %-5d  %22.16f %22.16f %22.16f\n"
                % (
                    i + 1,
                    self.types[i],
                    self.cartCoords[i, 0],
                    self.cartCoords[i, 1],
                    self.cartCoords[i, 2],
                )
            )


def main():

    outFile = sys.argv[1]
    outFH = open(outFile, "w")

    files = []
    for file in os.listdir("./"):
        if file.endswith(".out"):
            files.append(file)

    files.sort()
    nFiles = len(files)
    
    for iFile, file in enumerate(files):
        qe = QExpresso(inFile=file)
        qe.read()
        qe.fixCellMat()
        qe.write(outFH, nFiles, iFile + 1)

    outFH.close()


main()
