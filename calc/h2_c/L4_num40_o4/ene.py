import os
import commands
x_list = range(25)

for x in x_list[1:]:
    bl = x/10.0
    e = commands.getoutput("grep eig0 {0}/h2mole.out.dat | cut -d: -f2".format(x))
    print "{0} {1} {2}".format(bl, e, float(e)+1.0/bl)
