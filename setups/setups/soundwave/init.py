import re
import os
from pylab import *

N = int(sys.argv[1])
# ======================================================
#### CREATE .opt file
# ======================================================                                                                                                  
opt = open("template.opt","r")
OPT = opt.readlines()
opt.close()

list_N = "0"
for i in range(1,N+1): list_N += " {:d}".format(i)

new_opt = open("soundwave.opt","w")
for line in OPT:
    if re.search("%FLUIDS",line):
        new_opt.write("FLUIDS := "+list_N+"\n")
        continue
    if re.search("%NFLUIDS",line):
        new_opt.write("NFLUIDS = {:d}".format(N+1)+"\n")
        continue
    new_opt.write(line)
new_opt.close()

# ======================================================
#### CREATE .bound file
# ======================================================                                                                                                 
for i in range(0,N+1): os.system("cp soundwave.bound soundwave.bound.{:d}".format(i))
