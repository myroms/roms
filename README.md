# Regional Ocean Modeling System (ROMS)

The dynamical kernel of ROMS is comprised of four separate models, including the 
nonlinear (NLM), perturbation tangent linear (TLM), finite amplitude tangent 
linear (RPM), and adjoint (ADM). They are located in the Nonlinear, Tangent, 
Representer, and Adjoint sub-directories, respectively. The TLM and ADM were 
hand-coded from the discrete NLM code using the recipes of Giering and 
Kaminski (1998). The discrete adjoint is exact and is defined relative to 
the inner product that prescribes the L2-norm.

Detailed information about its numerical discretization, algorithms, usage, and 
tutorials can be found in the WikiROMS documentation portal at www.myroms.org/wiki

This version is mirrored from the community repositories at www.myroms.org for 
usage within the JCSDA Project and does not require registration on the ROMS 
website, and can be checkout from:
```
git clone https://github.com/JCSDA-internal/roms_src <source_dir>
```

Registered users of ROMS have access to:

- Official svn and git repositories
  ```
  svn checkout --username joe_roms https://www.myroms.org/svn/src/trunk <source_dir>

  git clone https://joe_roms@www.myroms.org/git/src  <source_dir>
  ```

- Idealized and Realistic test cases repository
  ```
  svn checkout https://www.myroms.org/svn/src/test  <test_dir>
  ```

- ROMS User's Forum:
  ```
  https://www.myroms.org/forum
  ```
