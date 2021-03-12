! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: obs_def_reanalysis_bufr_mod.f90 11289 2017-03-10 21:56:06Z hendric@ucar.edu $

! BEGIN DART PREPROCESS KIND LIST
!MESONET_U_WIND_COMPONENT,  QTY_U_WIND_COMPONENT,     COMMON_CODE
!MESONET_V_WIND_COMPONENT,  QTY_V_WIND_COMPONENT,     COMMON_CODE
!MESONET_GEOPOTENTIAL_HGT,  QTY_GEOPOTENTIAL_HEIGHT,  COMMON_CODE
!MESONET_SURFACE_PRESSURE,  QTY_SURFACE_PRESSURE,     COMMON_CODE
!MESONET_TEMPERATURE,       QTY_TEMPERATURE,          COMMON_CODE
!MESONET_SPECIFIC_HUMIDITY, QTY_SPECIFIC_HUMIDITY,    COMMON_CODE
! END DART PREPROCESS KIND LIST

! !!! Note about Specific Humidity observations:
! !!! UNITS in original BUFR are g/kg; This is converted to kg/kg by
! !!! the BUFR to obs_sequence conversion programs making it unnecessary
! !!! to multiply by 1000 at assimilation time.
! !!! PLEASE pay attention to units for specific humidity in models.

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/observations/forward_operators/obs_def_reanalysis_bufr_mod.f90 $
! $Id: obs_def_reanalysis_bufr_mod.f90 11289 2017-03-10 21:56:06Z hendric@ucar.edu $
! $Revision: 11289 $
! $Date: 2017-03-10 14:56:06 -0700 (Fri, 10 Mar 2017) $
