# Python standard library
import numpy as np
import multiprocessing as mp
import os
import sys

# Local library:

# Third-party library:


# -------------------------------------------------------------------------
#                               xyz writer
# -------------------------------------------------------------------------


def xyz_writer(filename, mode, label, total_atoms, xyz, box):

    with open(filename, mode) as output:

        output.write("%d\n" % total_atoms)

        np.savetxt(output, np.c_[[box]])

        for i in range(total_atoms):

            output.write("%d %.15f %.15f %.15f\n"
                         % (label[i], xyz[i, 0], xyz[i, 1], xyz[i, 2]))

    return None


# -------------------------------------------------------------------------
#                          Fortran xyz writer
# -------------------------------------------------------------------------
