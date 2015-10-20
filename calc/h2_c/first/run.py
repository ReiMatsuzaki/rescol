import os

os.system('rm -r common')
os.system('mkdir common')
os.chdir('common')
os.system("""python ../../../../script/calc_y2mat.py
-l0 0 
-l1 2 
-l1guess 0 
-l2guess 0 | tee calc_y2mat.out""")
os.chdir('..')


