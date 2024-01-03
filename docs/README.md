# Regional Ocean Modeling System (ROMS)

![ROMS_Picture](https://github.com/myroms/roms/assets/23062912/d72765ed-9d55-4109-84fc-c51b05832adb)

# License

**Copyright (c) 2002-2024 The ROMS/TOMS Group**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Overview

**ROMS** solves the free-surface, hydrostatic, flux form of the primitive
equations over variable bathymetry using stretched terrain following in the
vertical and orthogonal curvilinear coordinates in the horizontal. The finite
volume grid is discretized on a staggered Arakawa C-grid. Detailed information
about its governing equations, numerical discretization, algorithms, usage, and
tutorials is available in the **WikiROMS** documentation portal at
**`www.myroms.org/wiki`**.

The dynamical kernel of **ROMS** is comprised of four separate models, including
the nonlinear (**NLM**), perturbation tangent linear (**TLM**), finite amplitude
tangent linear (**RPM**), and adjoint (**ADM**). They are located in the
**Nonlinear**, **Tangent**, **Representer**, and **Adjoint** sub-directories,
respectively. The **TLM** and **ADM** were hand-coded from the discrete **NLM**
code using the recipes of Giering and Kaminski (1998). Therefore, any change to
its dynamical and numerical kernels will affect the symmetry of the **TLM** and
**ADM** operators. The discrete adjoint is exact and is defined relative to
the inner product that prescribes the L2-norm.

This official community version of **ROMS** is developed and maintained at Rutgers,
The State University of New Jersey, New Brunswick, New Jersey, USA.

# Instructions

The **ROMS**  framework is intended for users interested in ocean modeling. It
requires an extensive background in ocean dynamics, numerical modeling, and
computers to configure, run, and analyze the results to ensure you get the
correct solution for your application. Therefore, **we highly recommend** users
register at https://www.myroms.org and set up a **username** and **password** to
access the **ROMS** forum, email notifications for bugs/updates, technical
support from the community, **trac** code maintenance history, tutorials, workshops,
and publications. The user's **ROMS** forum has over 24,000 posts with helpful
information. Technical support is limited to registered users. We **do not**
provide user technical support, usage, or answers in **GitHub**.

The **GitHub** version becomes the official **git** repository for downloading,
updating, improving, and correcting defects/bugs to the **ROMS** source code.
Also, it is the version used in the **ROMS-JEDI** interface hosted at
https://github.com/JCSDA-internal, which is currently private. We will update
the git and svn repositories hosted at https://www.myroms.org for some time,
but encourage the User to transition to **GitHub**, which gives you access to the
feature branches. Use the following command to download the **ROMS** source code:
```
git clone https://github.com/myroms/roms.git
svn checkout --username joe_roms https://www.myroms.org/svn/src/trunk <source_dir>
```
The idealized and realistic **ROMS** Test Cases and the Matlab processing
software can be downloaded from:
```
git clone https://github.com/myroms/roms_test.git
git clone https://github.com/myroms/roms_matlab.git
```
We highly recommend that Users define the **`ROMS_ROOT_DIR`** variable in their
computer shell logging environment, specifying where the User cloned/downloaded
the **ROMS** source code, Test Cases, and Matlab processing software:
```
setenv ROMS_ROOT_DIR  MyDownlodLocationDirectory
```
The **build** scripts will use this environmental variable when compiling any of
the **ROMS Test Cases** without the need to customize the location of the
**ROMS** source code. Also, it is used for loading the path of Matlab scripts in
the **startup.m** configuration file.

The **doxygen** version of **ROMS** is available at:
```
https://www.myroms.org/doxygen
```
Registered users of **ROMS** have access to:

- **ROMS** User's Forum for technical support from the community:
  ```
  https://www.myroms.org/forum
  ```
- **ROMS Trac** source code maintenance and evolution:
  ```
  https://www.myroms.org/projects/src
  ```
- **WikiROMS** documentation and tutorials plus editing:
  ```
  https://www.myroms.org/wiki
  ```
