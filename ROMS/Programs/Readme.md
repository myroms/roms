*
* git $Id$
* svn $Id: Readme.md 1195 2023-08-31 01:52:55Z arango $
***********************************************************************
*  Copyright (c) 2002-2023 The ROMS/TOMS Group                        *
*    Licensed under a MIT/X style license                             *
*    See License_ROMS.txt                                             *
***********************************************************************
*

# Auxiliary ROMS Programs:

- **types.F**: It checks the precision and range of floating-point
  variables in a particular compiler architecture. To compile and
  link, use:

```
  ifort -o types.x types.F
  gfortran -o types.x types.F
```

- **yaml_parser_test.F**: It tests the yaml_parser.F available in
  ROMS. It requires the Fortran 2003 standard and a couple of features of
  the Fortran 2008 standard, which are available in modern compilers.
  We have reports from users having issues compiling with older versions of
  **gfortran**. To compile and link, use:

```
  ifort -o yaml_parser_test.x yaml_parser_test.F
  gfortran -o yaml_parser_test.x yaml_parser_test.F
```
It is easy to run by specifying the desired YAML files available in ROMS:

```
  > yaml_parser_test.x

  Enter YAML filename:  ../../ESM/coupling_esmf.yaml

  Enter YAML filename:  ../../ESM/roms_cmeps.yaml

  Enter YAML filename:  ../External/varinfo.yaml

  Enter YAML filename:
```
