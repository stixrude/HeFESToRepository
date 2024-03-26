# HeFESTo
by L. Stixrude and C. Lithgow-Bertelloni, 2005-

### Citations
If you use HeFESTo results in a publication, please cite the following papers:
1. Stixrude, L. and C. Lithgow-Bertelloni, Thermodynamics of mantle minerals: 1. Physical properties, Geophysical Journal International, 162, 610-632, 2005.
2. Stixrude, L. and C. Lithgow-Bertelloni, Thermodynamics of mantle minerals II, Phase equilibria, Geophysical Journal International, 184, 1180-1213, 2011.
3. Stixrude, L. and C. Lithgow-Bertelloni, Thermodynamics of mantle minerals III, The role of iron, Geophysical Journal International, submitted, 2023.
4. The citation to the parameter file that you have chosen.

### Installation
> There is an involved installation doc that can be referenced for properly
> configuring HeFESTo (`docs/installation`).

## Quick Start
You will need a couple dependencies to run HeFESTo
- `gfortran`, `lapack`, `blas`, `make`, and `nlopt`
  - refer to [NLOPT](https://nlopt.readthedocs.io/en/latest/)'s website for the installation procedure

Copy the correct makefile from the architecture directory, `arch`
```bash
cp arch/makefile.xxx makefile
```

Compile the project
```bash
make all
```

You will need to gather a parameter set for HeFESTo to run on. These can be
found at
[https://github.com/stixrude?tab=repositories](https://github.com/stixrude?tab=repositories).
Where the parameters are repositories designated as HeFESTo_Parameters_xxxx.

In the `BENCHMARK` directory, the `control` file is how you modify your inputs.
Information about this can be found in `docs/hefesto.pdf`. You will need to modify
`line 13` so that it is pointed to your `parameter` directory.
> NOTE: a path that is too long will produce an error.

Copy the executable into `BENCHMARK` and run the program.
```bash
 cp main BENCHMARK/
 cd BENCHMARK
 ./main
```

If there is an issue, refer to `docs/installation` for properly configuring
HeFESTo


### __WARNING__
> This version is for testing purposes only, which means we grant the right to
> download, compile, and run the software, but do not grant the right to
> redistribute modified versions of the software either in compiled or in
> source-code form. If you have questions about the software, suggestions for
> improvements, or would like to ask for additional permissions or data files,
> please contact the authors:

###### Contact
Lars Stixrude: [lstixrude@epss.ucla.edu](mailto:lstixrude@epss.ucla.edu)

Carolina Lithgow-Bertelloni: [clb@epss.ucla.edu](mailto:clb@epss.ucla.edu)
