import os
import sys

src_dir =os.environ["RESCOL_DIR"] + "/src"
script_dir = os.environ["RESCOL_DIR"] + "/script"
orig_dir = os.path.abspath(".")

sys.path.append(src_dir)
from rescol import cd, with_dir, keyval_to_dict

bond_length_list = [(str(14), 1.4)]



