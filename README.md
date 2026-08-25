# Ternary-Valued Finite-Difference Time-Domain (FDTD) Method
Computer codes that create the figures for the `Ternary-Valued Finite-Difference Time-Domain Method` article

** Figure 1 **

```bash
rm -f convergence.txt
   for p in "10000 2000" "5000 4000" "2500 8000" \
            "1250 16000" "625 32000"; do
     set -- $p
     gcc -O2 -DMC=$1 -DNSTEP=$2 -o cv ternary_convergence_test.c -lm
     ./cv
   done
```
Then, use GNU Octave/MATLAB:

```bash
octave plot_convergence_V2.m
```
** Figure 2 **

First, generate `probe_fixed.txt`:

```bash
gcc -O3 -o aef ternary_fdtd_fixed_scattering.c -lm
./aef
```

Then, use GNU Octave/MATLAB:

```bash
octave compare.m
```

