import argparse
import sys
import re
import numpy as np


parser = argparse.ArgumentParser(description = "read file and extract data")

parser.add_argument("file_name", nargs='*')
parser.add_argument("-d", default=":")

parser.add_argument("-key", default=[], nargs='+')
args = parser.parse_args()

if len(args.file_name) == 0:
    print "not implemented yet."
    sys.exit(1)

re_list = [re.compile("^ *{0} *{1} *(.*)$".format(k, args.d, ))
           for k in args.key]

def maybe_group(res_search):
    if res_search:
        return res_search.group(1)
    else:
        return None

def remove_none(xs):
    return [x for x in xs if x!=None]

def connect_strs(str_list, sep):
    return reduce(lambda a,b:a+sep+b, str_list)

with open(args.file_name[0]) as f:
    lines = f.readlines()
    data_table = [remove_none([ maybe_group(re0.search(line))
                                for re0 in re_list ])
                  for line in lines]

    datas = np.array([remove_none([maybe_group(re0.search(line)) for line in lines])
                      for re0 in re_list])

    print connect_strs(args.key, " ")
    for ds in datas.T:
        print connect_strs(ds, " ")
        #print reduce(lambda a,b: a+" "+b, ds)

