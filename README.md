### INSTALL

Edit Makefile

make

### PERFORMANCE

Input: test4.inp

Solute: ProteinG (846 atoms)

Solvent: Water (Tip3p model)

Cell: 128^3 Angstrom^3 (256^3 grids)

<pre>
GPU      Time [s]
-----------------
K20C     43.8
P100     14.5
Titan V  12.5
V100      9.9
</pre>

### REFERENCE
1. Y. Maruyama, and N. Yoshida, "RISMiCal: A software package to perform fast RISM/3D-RISM calculations," J. Comput. Chem., (2024) 45, 1470-1482 (DOI: 10.1002/jcc.27340)

2. Y. Maruyama and F. Hirata, "Modified Anderson method for accelerating 3D-RISM calculations using graphics processing unit," J. Chem. Theory Comput., (2012) 8, 3015-3021 (DOI: 10.1021/ct300355r)
