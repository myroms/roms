# **Regional Ocean Modeling System (ROMS)**

![roms-jedi logo](https://www.myroms.org/trac/roms_src_600px.png)

**ROMS** solves the free-surface, hydrostatic, flux form of the primitive
equations over variable bathymetry using stretched terrain-following in the
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
The State University of New Jersey, New Brunswick, New Jersey, USA. The scarlet
flag in the above logo indicates the location of our institution. Currently, this
**git** repository contains the following branches:

- **master**: Tagged versions and the latest stable release version of **ROMS**
- **develop**: Main developing branch of **ROMS**. Not recommended for public
 consumption but passes the internal tests. It is intended for ROMS superusers
 and beta testers.
- **feature branches**: Research and new development branches recommended to
superusers and beta testers.

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

This **GitHub** version becomes the official **git** repository for downloading,
updating, improving, and correcting defects/bugs to the **ROMS** source code.
Also, it is the version used in the **ROMS-JEDI** interface hosted at
https://github.com/JCSDA-internal, which is currently private. Use the following
command to download the **ROMS** source code:
```
git clone https://github.com/myroms/roms.git                 (default)
git clone https://github.com/myroms/roms.git <source_dir>
```
The idealized and realistic **ROMS** Test Cases and the Matlab processing
software can be downloaded from:
```
git clone https://github.com/myroms/roms_test.git
git clone https://github.com/myroms/roms_matlab.git
```
The **doxygen** version of **ROMS** is available at:
```
https://www.myroms.org/doxygen
```
Registered users of **ROMS** have access to:

- **ROMS** User's Forum for technical support from the community:
  ```
  https://www.myroms.org/forum
  ```
- **Trac** source code maintenance and evolution:
  ```
  https://www.myroms.org/projects/src
  ```
- **WikiROMS** documentation and tutorials plus editing:
  ```
  https://www.myroms.org/wiki
  ```
