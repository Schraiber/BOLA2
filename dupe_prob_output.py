import argparse
import numpy as np
from subprocess import Popen, PIPE
import sys

parser = argparse.ArgumentParser(description="Compute the probability of an allele having a certain configuration in modern populations while being absent in archaics, using msms to simulate selection")
parser.add_argument("msdir",nargs=1,type=str)
parser.add_argument("-tNDYS",default=26000,type=float,help="Time of N/D divergence from modern humans")
parser.add_argument("-tND",default=21000,type=float,help="Time of N/D divergence")
parser.add_argument("-tYS",default=8000,type=float,help="Time of S/Y divergence")
parser.add_argument("-tYb",default=6300,type=float,help="Time of Y pop size increase") 
parser.add_argument("-tYs",default=1500,type=float,help="Time of Y pop size decline")
parser.add_argument("-nYs",default=10000,type=float,help="Size of Y currently")
parser.add_argument("-nYb",default=45000,type=float,help="Size of Y expansion")
parser.add_argument("-nS",default=10000,type=float,help="Size of S currently")
parser.add_argument("-nYS",default=24000,type=float,help="Size of Y/S ancestor (and part of Y)")
parser.add_argument("-nNDYS",default=21600,type=float,help="Size of N/D/Y/S ancestor")
parser.add_argument("-nND",default=500,type=float,help="Size of N/D")
parser.add_argument("-m",default=0,type=float,help="Mutation rate (nb: for now, this does nothing)")
parser.add_argument("-S",default=8,type=int,help="Number of S chromosomes sampled")
parser.add_argument("-Y",default=270,type=int,help="Number of Y chromosomes sampled")
parser.add_argument("-N",default=2,type=int,help="Number of N chromosomes sampled")
parser.add_argument("-D",default=2,type=int,help="Number of D chromosomes sampled")
parser.add_argument("-tSel",default=11280,type=float,help="Age of the selected allele")
parser.add_argument("-fi",default=1./20000.,type=float,help="Initial frequency of allele")
parser.add_argument("-sAA",default=0,type=float,help="Selection coefficient of AA genotype (not pop scaled!)")
parser.add_argument("-sAa",default=0,type=float,help="Selection coefficient of Aa genotype (not pop scaled!)")
parser.add_argument("-nreps",default=1,type=int,help="Number of replicates")

args = parser.parse_args()

tNDYS = float(args.tNDYS)/(4*args.nYs)
tND = float(args.tND)/(4*args.nYs)
tYS = float(args.tYS)/(4*args.nYS)
tYb = float(args.tYb)/(4*args.nYs)
tYs = float(args.tYs)/(4*args.nYs)
nYs = float(args.nYs)/args.nYs
nYb = float(args.nYb)/args.nYs
nS = float(args.nS)/args.nYs
nYS = float(args.nYS)/args.nYs
nNDYS = float(args.nNDYS)/args.nYs
nND = float(args.nND)/args.nYs
theta = 0
tSel = float(args.tSel)/(4*args.nYs)
sAA = 2*args.nYs*args.sAA
sAa = 2*args.nYs*args.sAa

outputFileName = "test_freq_age_%d.txt"%args.tSel

outputFile = open(outputFileName,"w")

myCMD = ["java","-jar",args.msdir[0],str(args.S+args.Y+args.N+args.D),str(args.nreps),"-t",str(theta),"-N", str(args.nYs), "-I","4",str(args.S),str(args.Y),str(args.N),str(args.D),"-n","1",str(nS),"-n","2",str(nYs),"-n","3",str(nND),"-n","4",str(nND),"-en",str(tYs),"2",str(nYb),"-en",str(tYb),"2",str(nYS),"-ej",str(tYS),"1","2","-ej",str(tND),"3","4","-ej",str(tNDYS),"4","2","-en",str(tNDYS),"2",str(nNDYS),"-SI",str(tSel),"4","0",str(args.fi),"0","0","-SAA",str(sAA),"-SaA",str(sAa), "-Smark"]

sys.stderr.write(' '.join(myCMD)+"\n")

unprocessed = Popen(myCMD,stdout=PIPE)
processed = Popen(["tail","-n","+4"],stdin=unprocessed.stdout,stdout=PIPE)
processed = Popen(["grep","-v","positions"],stdin=processed.stdout,stdout=PIPE)
ms_out = Popen(["grep","-v","//"],stdin=processed.stdout,stdout=PIPE)

good_sites = 0
sites = 0


outputFile.write('num_S\tnum_Y\tnum_N\tnum_D\n')
for line in ms_out.stdout:
	if line[0] == "s":
		sam = 0
		numSites = int(line.strip().split(" ")[1])
		sites += numSites
		cur_S = np.zeros((args.S,numSites))
		cur_Y = np.zeros((args.Y,numSites))
		cur_N = np.zeros((args.N,numSites))
		cur_D = np.zeros((args.D,numSites))
		continue
	if sam < args.S:
		cur_S[sam,:] = map(lambda x: int(x != "0"),list(line.strip()))
		sam += 1
		continue
	if sam < args.S+args.Y:
		cur_Y[sam-args.S,:] = map(lambda x: int(x != "0"),list(line.strip()))
		sam += 1
		continue
	if sam < args.S+args.Y+args.N:
		cur_N[sam-args.S-args.Y,:] = map(lambda x: int(x != "0"),list(line.strip()))
		sam += 1
		continue
	if sam < args.S+args.Y+args.N+args.D:
		cur_D[sam-args.S-args.Y-args.N,:] = map(lambda x: int(x != "0"),list(line.strip()))
		sam += 1
		continue
	if sam == args.S+args.Y+args.N+args.D:
		num_S = np.sum(cur_S,0)
		num_Y = np.sum(cur_Y,0)
		num_N = np.sum(cur_N,0)
		num_D = np.sum(cur_D,0)
		outputFile.write('%i\t%i\t%i\t%i\n'%(num_S,num_Y,num_N,num_D))
		#for site in range(numSites):
		#	if num_S[site] == args.S and num_Y[site] == args.Y and num_N[site] == 0 and num_D[site] == 0:
		#		good_sites += 1
		continue
				
num_S = np.sum(cur_S,0)
num_Y = np.sum(cur_Y,0)
num_N = np.sum(cur_N,0)
num_D = np.sum(cur_D,0)
outputFile.write('%i\t%i\t%i\t%i\n'%(num_S,num_Y,num_N,num_D))
#for site in range(numSites):
#	if num_S[site] == args.S and num_Y[site] == args.Y and num_N[site] == 0 and num_D[site] == 0:
#		good_sites += 1

