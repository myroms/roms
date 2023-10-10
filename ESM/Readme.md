
<img width="500" alt="image" src="https://github.com/myroms/roms/assets/23062912/2fa815f4-df51-4671-b9ed-188995b87a7b">

## Earth System Model (ESM) Coupling:

This directory contains several files used for multi-model coupling
using the Earth System Modeling Framework (**ESMF**) with the National
Unified Operational Prediction Capability (**NUOPC**) layer.

The **NUOPC** layer is a simplified interface on top of the **ESMF** library
(version **8.0** or higher) that provides conventions and templates to facilitate
the smooth coupling between Earth System Model (**ESM**) components. **ROMS**
offers two distinct **NUOPC**-based modules for its native coupling framework
and a standalone interface for third-party coupling systems. That is, **ROMS**
coupling infrastructure allows both **DRIVER** and **COMPONENT** methods of operation.

<p align="center">
    <img src="https://www.myroms.org/trac/ROMS_Coupling.png" width=50% height=50%>
</p>

In the **`DRIVER`** method, it provides all the interfaces needed to couple
to other **ESM** components, including the executable driver, **NUOPC**-based generic
**ESM** component services, model gridded components or **NUOPC** _cap_ modules,
connectors between components for re-gridding source and destination
fields, input scripts, and coupling metadata management.

A **NUOPC** model _cap_ is a Fortran code layer that sits on top of the **ESM**
component, making calls to the numerical kernel via the _`initialize`_,
_`run`_, and _`finalize`_ computational phases.

Alternatively, in the **`COMPONENT`** method, the **NUOPC**-based **ROMS** _cap_
module **`Master/cmeps_roms.h`** is provided, and it can be adapted and
incorporated into other **NUOPC**-based coupling systems, like the
[UFS_coastal Framework](https://github.com/oceanmodeling/ufs-coastal)
using the Community Mediator for Earth Prediction Systems (**CMEPS**) and
the Community Data Models for Earth Predictive Systems (**CDEPS**) coupling
interfaces.

Currently, the following files are avaiabled for coupling with the **ESMF/NUOPC**
library:

| $\textcolor{blue}{\textsf{Coupling Files}}$ | $\textcolor{blue}{\textsf{Description}}$ |
|-----------------------|-------------|
| **cmeps_roms.h**      | **ROMS** standalone interface for the **UFS** using **CDEPS/CMEPS** |
| **coupler.F**         | **ESMF/NUOPC** or **MCT** native coupler module |
| **esmf_coupler.h**    | **ESMF** Regridding operators using connectors |
| **esmf_atm.F**        | Atmosphere gridded model (**ATM**) component module |
| **esmf_atm_coamps.h** | **COAMPS** gridded component **NUOPC** layer module |
| **esmf_atm_regcm.h**  | **RegCM** gridded component **NUOPC** layer module |
| **esmf_atm_wrf.h**    | **WRF** gridded component **NUOPC** layer module |
| **esmf_data.F**       | Coupling **DATA** component used in incongruent grids coupling |
| **emsf_driver.h**     | Configures, Creates, Initialize, Run, and Finalizes the **ESM** coupling system |
| **esmf_esm.F**        | Sets **ESM** gridded components services shared-objects and RunSequence |
| **esmf_ice.F**        | Seaice gridded model (**SEA-ICE**) component module |
| **esmf_ice_cice6.h**  | **CICE6** gridded component **NUOPC** layer module |
| **esmf_roms.h**       | **ROMS** native coupling framework **NUOPC** layer module |
| **esmf_roms.F**       | **ROMS** native or standalone Ocean component module |
| **esmf_wav.F**        | Waves gridded component (**WAVE**) component module |
| **esmf_wav_wam.h**    | **WAM** gridded component **NUOPC** layer module |
| **mod_esmf_esm.F**    | **ROMS** native coupling framework support module |

## Coupling Design:

The strategy is to couple to other **ESM** components with no or minimal
changes to the code distributed by developers or repositories.  The User is
responsible for subscribing to those repositories, downloading, and installing
the other **ESM** components. The coupling with such **ESM** components is expected
not to be affected by its version (previous, current, or future) since the
**NUOPC**-based layer is generic.

However, sometimes, we need to circumvent technical problems when coupling
to other **ESM** components and provide build scripts to facilitate
compiling, linking, and correcting interference to deprecated libraries.
Therefore, this directory also contains scripts and modified **ESM** component
files that substitute the ones distributed from source repositories
to solve such technical issues.

---

## WRF Coupling Notes:

* The coupling system has been tested with WRF versions 4.1 and up.

* To compile the coupled system with 'gfortran', you need to activate
  the environmental variable GFORTRAN_CONVERT_UNIT:

  export GFORTRAN_CONVERT_UNIT='big_endian'  or
  setenv GFORTRAN_CONVERT_UNIT big_endian

  to avoid end-of-file when reading binary files like 'RRTMG_LW_DATA'

* The 'gfortran' is more strict, and the execution fails because
  of unassigned INTENT(OUT) variables 'fname' and 'n2' in source file
  share/mediation_integrate.F

--- a/share/mediation_integrate.F
+++ b/share/mediation_integrate.F
@@ -2408,6 +2408,8 @@ SUBROUTINE open_hist_w ( grid , config_flags, stream, alarm_id, &
    ENDIF
    ierr = 0
+   fname = ""
+   n2 = ""
    ! Note that computation of fname and n2 are outside of the oid IF statement
    ! since they are OUT args and may be used by callers even if oid/=0.


### File Description:

* coupling_esmf.in:

  Standard input script for ROMS when coupling with the ESMF/NUOPC
  library. It is well documented and sets the coupling system. To
  submit a job, we could use for example:

  mpirun -np 8 romsM coupling.in > & log &

* coupling_esmf.dat:

  Coupling metadata defining import and export fields.

* build_cice.csh:

  A friendlier CSH script to compile CICE.

* build_wrf.csh, build_wrf.sh:

  CSH and BASH compiling scripts for WRF to facilitate easy compiling
  and linking.  It also corrects several technical issues with very old
  ESMF library interference and incorrect NetCDF4 library dependencies.
  Many of the sections of this script that do not require customization
  have been off-loaded into sub-scripts that are called from the WRF
  build script. They are: wrf_patch.*, wrf_links.* and wrf_move.*)

* wrf_restore.csh, wrf_restore.sh:

  CSH and BASH script to restore the WRF source code directory to its
  original checkout state by undoing the changes made by wrf_patch.csh
  or wrf_patch.sh.

* wrf_patch.csh, wrf_patch.sh:

  CSH and BASH to check whether the WRF souce code has been patched for
  NetCDF4 library dependencies, added configure options, creating clean
  f90 files for debugging, renaming modules to WRF_ESMF_*, and correcting
  optional argument from defaultCalendar or defaultCalKind in
  ESMF_Initialize call.

* wrf_move.csh, wrf_move.sh:

  If -move flag is used when executing build_wrf.*, this script moves
  the WRF objects and executables needed to run WRF in the coupled
  ESMF/NUOPC system.

* wrf_links.csh, wrf_links.sh:

  If -move flag is used when executeing build_wrf.* and ${WRF_CASE} is
  set to 'em_real', this script creates the data links for running the
  'em_real' executable.

* build_wps.csh, build_wps.sh

  CSH and BASH compiling scripts for WPS to facilitate easy compiling
  and linking. It patches the util/src/Makefile to allow compiling
  with parallel enabled NetCDF4/HDF5. If which_MPI is set to 'intel',
  it also corrects the MPI compiler names in configure.wps.

* *.runconfig

  The ESMF RunSequence configuration file sets how the ESM components
  are connected and coupled.

* roms_fields.yaml (currently not used)

  Right now this file is not used but we hope to use this file in the
  future to handle exchange fields for ESMF/NUOPC coupling.

The files below were adapted to work with WRF Versions 4.1.x - 4.3

  - wrf_configure:

       Replaces ${WRF_ROOT_DIR}/configure

       Reworking linking NetCDF4 library dependencies.

  - wrf_Makefile:

       Replaces ${WRF_ROOT_DIR}/Makefile

       Reworking linking NetCDF4 library dependencies.

  - wrf_postamble:

       Replaces ${WRF_ROOT_DIR}/arch/postamble

       Reworking linking NetCDF4 library dependencies.

  - wrf_configure.defaults:

       Replaces ${WRF_ROOT_DIR}/arch/configure.defaults

       Added Intel (ifort/icc) and OpenMPI for MacOS

   - wrf_Makefile.esmf:

       Replaces ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile

       Rename 'ESMF_*' modules to 'WRF_ESMF_*' to avoid conflicts with
       new versions of the ESMF/NUOPC library. Everything is done during
       C-preprocessing so original files are not modified.

    - wrf_Test1.F90:

       Replaces ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90

       Corrects bug in optional argument to 'ESMF_Initialize' call
       from 'defaultCalendar' to 'defaultCalKind'.
