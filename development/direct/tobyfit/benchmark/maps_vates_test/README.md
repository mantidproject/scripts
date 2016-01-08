## TobyFit Single Q Coordinate Benchmark

This directory contains the files necessary to produce a simulation using
the QCoordinate model in Mantid to pick out a single Q coordinate as the
signal.

There are 3 lattice types:

- cubic
- hexagonal
- tetragonal

The data files are unfortunately too large to store here.

Information from email:

There are 18 files in total, which correspond to 3 different crystal coordinate systems (cubic, tetragonal and hexagonal), each with the signal given as h, k or l, calculated both by Tobyfit and Horace. The filename convention is fake_<crystal-system>_<Horace if calculated by Horace>_<coordinate>.sqw.

The lattice parameters I used were 5x5x5 Angstroms for the cubic, and a=b=3, c=10 Angstroms for the other two. For the cubic and tetragonal all lattice angles are 90 degrees, whereas for hexagonal it is alpha=beta=90, gamma=120 degrees. In all cases u=[0,0,1], v=[1,0,0], where u//ki, v defining the scattering plane together with u. The incident energy Ei=50, and all the goniometer offset angles were zero.
