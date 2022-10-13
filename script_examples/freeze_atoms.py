
# this script is designed to be used from the command line
# by python /path/to/freeze_atoms.py [upper_threshold]
# every atom with z-coordinate greater than [upper_threshold] is allowed to relax
# other atoms are frozen  

import sys
import poscar_class as pc

poscar = pc.poscar("./POSCAR")
poscar.freeze(float(sys.argv[1]))
seld_false = sum(x == 'F F F' for x in poscar.coordinates_seldyn)
seld_true = sum(x == 'T T T' for x in poscar.coordinates_seldyn)
print("Not frozen:", seld_true, "Frozen:", seld_false)

poscar.save("POSCAR_f")



