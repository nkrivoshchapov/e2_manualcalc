import sys,os,glob,ntpath
import numpy as np

sep = ";"
HtoKC = 627.509474063
KCtoH = 1/HtoKC

class FockMatrix:
    def __init__(self, filename, nbasis):
        self.n = nbasis
        self.fm = np.empty([self.n, self.n])
        self.count = 0
        self.cur_row = 0
        self.cur_col = 0

        lines = open(filename, "r").readlines()
        for line in lines[3:]:
            parts = line.split()
            for part in parts:
                self.append_elem(float(part))

    def append_elem(self, newnum):
        if self.cur_col > self.cur_row:
            self.cur_row += 1
            self.cur_col = 0
        self.fm[self.cur_row][self.cur_col] = newnum
        self.fm[self.cur_col][self.cur_row] = newnum
        self.cur_col += 1

def gete2lines(lines):
    for i, line in enumerate(lines):
        if "Second Order Perturbation Theory Analysis of Fock Matrix in NBO Basis" in line:
            start = i + 8
    for i, line in enumerate(lines):
        if "Natural Bond Orbitals (Summary):" in line:
            end = i - 2
    if 'start' not in locals() or 'end' not in locals():
        raise Exception("Can't find second order analysis section")
    return start, end

def getocc(occ, lines):
    count = 0
    for line in lines:
        if "." in line:
            occ[count] = float(line[40:49])
            count += 1
    if count != nbasis:
        raise Exception("Error while parsing occupancy data")

nbopairs = []
pairlines = open("nbo_pairs.dat", "r").read().split("\n")
for line in pairlines:
    parts = line.split()
    nbopairs.append([int(parts[0]), int(parts[1])])

np.seterr(all='raise')
csvlines = [sep.join(["File", "NBO1", "NBO2", "E2(kcal/mol)"])]
for logfile in glob.glob("./*.log"):
    loglines = open(logfile, "r").readlines()

    for i, line in enumerate(loglines):
        if "NBasis=" in line:
            nbasis = int(line.split("RedAO")[0].split("=")[1])
        if "Natural Bond Orbitals (Summary):" in line:
            occstart = i + 5
        if "Total Lewis" in line:
            occend = i - 1

    e2start, e2end = gete2lines(loglines)
    myfm = FockMatrix(logfile.replace(".log", ".69").upper(), nbasis)
    occ = np.empty(nbasis)
    getocc(occ, loglines[occstart:occend])

    for pair in nbopairs:
        i = pair[0] - 1
        j = pair[1] - 1
        try:
            occsum = abs(occ[i] + occ[j])
            if occsum > 2:
                occsum = 4 - occsum
            fijsq = myfm.fm[i][j]**2
            enerdiff = abs(myfm.fm[i][i] - myfm.fm[j][j]) # in a.u.
            e2 = occsum*fijsq/enerdiff * HtoKC # in kcal/mol
        except:
            print("Error occurred:")
            print("OccSum = " + repr(occsum))
            print("Fij = " + repr(fijsq))
            print("dE = " + repr(enerdiff))
        csvlines.append(sep.join([ntpath.basename(logfile), str(pair[0]), str(pair[1]), str(e2)]))

wfile = open("res.csv", "w")
wfile.write("\n".join(csvlines))
wfile.close()