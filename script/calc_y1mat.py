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

l_list = [L for L in range(l0, l1+2, 2)]
q_list = uniq(flatten([ls_non_zero_YYY(L1, L2)
                       for L1 in l_list for L2 in l_list]))

cd(args.target_dir)
for q in q_list:
    pq_mat = [[y1mat_Pq((L1,m), q, (L2,m))
               for L1 in l_list] for L2 in l_list]
    fn = "p{0}_y1mat.dat".format(q)
    write_coo_mat(pq_mat, fn)

s_mat = [[ 1 if L1==L2 else 0 
           for L1 in l_list]
         for L2 in l_list]
write_coo_mat(s_mat, "s_y1mat.dat")

l_mat = [[ L1*(L1+1) if L1==L2 else 0 
           for L1 in l_list] for L2 in l_list]
write_coo_mat(l_mat, "l_y1mat.dat")    
