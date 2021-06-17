# MISC
1. ` QE_lammps_alat.py` : The code will convert Quantum Espresso output file to lammps dump format. 
Default: 
- The Quantum espresso output files  extension is \*.out
- The python code and outputs should be in same directory. 
- The lammps dump format output should have extension .dump
# HOW to run
`
python3.9 QE_lammps_alat.py file.dump
`
2. `QE_lammps_spin_noncolin.py`: Conversion with non-collinear spin output in Quantum Espresso
3. `QE_lammps_spin.py`: Conversion with collinear spin output in Quantum Espresso
4. `liquid.py` : A python code to generate amorphous data, output will be in [extended xyz format](https://web.archive.org/web/20190811094343/https://libatoms.github.io/QUIP/io.html#extendedxyz)
