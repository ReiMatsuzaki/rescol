#/usr/bin/python
import argparse
from os.path import abspath, dirname, join
import sys
sys.path.append(join(dirname(dirname(__file__)), "src"))
sys.path.append("../src")
from rescol import *

parser = argparse.ArgumentParser(description = "compute matrix represented Spherical Harmonics Y")
#parser.set_defaults(l0=0, l1=2);
parser.add_argument("-l0", default=0, type=int)
parser.add_argument("-l1", default=2, type=int)
parser.add_argument("-m", default=0, type=int)
parser.add_argument("-t", "--target_dir", default=".")
args = parser.parse_args()

if (args.l0 - args.l1)%2 == 1:
    print("ArgumentError: l0-l1 must be even number ")
    print("l0: {0}".format(args.l0))
    print("l1: {0}".format(args.l1))
    print("calc_y1mat stop")
    sys.exit(1)

if args.l0>args.l1:
    print("ArgumentError: l0 must be greater or equal to l1")
    print("l0: {0}".format(args.l0))
    print("l1: {0}".format(args.l1))
    print("calc_y1mat stop")
    sys.exit(1)

if args.l0 < abs(args.m):
    print("ArgumentError: m must be in [-L,L] for all L");
    print("l0: {0}".format(args.l0))
    print("l1: {0}".format(args.l1))
    print("m:  {0}".format(args.m))
    print("calc_y1mat stop")
    sys.exit(1)

if not os.path.exists(args.target_dir):
    print("ArgumentError: target directory does not exist")
    print("target_dir:  {0}".format(args.target_dir))
    print("calc_y1mat stop")
    sys.exit(1)    

l0 = args.l0
l1 = args.l1
m = args.m

ys = get_coupledY_set(l1, m, True, True);
qmax = 2*max([y.L1 for y in ys])
qs = range(qmax+1)

cd(args.target_dir)
s_mat = [[ 1 if i==j else 0 
           for i in range(len(ys))]
         for j in range(len(ys))]
write_coo_mat(s_mat, "s_y2mat.dat")

l_1_mat = [[y.L1*(y.L1+1) if (i==j) else 0
           for (y, i) in with_index(ys)]
          for (yy, j) in with_index(ys)]
write_coo_mat(l_1_mat, "l_1_y2mat.dat")

l_2_mat = [[y.L2*(y.L2+1) if (i==j) else 0
           for (y, i) in with_index(ys)]
          for (yy, j) in with_index(ys)]
write_coo_mat(l_2_mat, "l_2_y2mat.dat")

for q in qs:
    pq_1A_mat = [[y2mat_Pq_r1A(y1, q, y2)
                  for y1 in ys] for y2 in ys]
    write_coo_mat(pq_1A_mat, "p{0}_A1_y2mat.dat".format(q), True)
    pq_2A_mat = [[y2mat_Pq_r2A(y1, q, y2)
                  for y1 in ys] for y2 in ys]
    write_coo_mat(pq_2A_mat, "p{0}_A2_y2mat.dat".format(q), True)
    pq_12_mat = [[y2mat_Pq_r12(y1, q, y2) 
                  for y1 in ys] for y2 in ys]
    write_coo_mat(pq_12_mat, "p{0}_12_y2mat.dat".format(q), True)
