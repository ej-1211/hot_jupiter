module mod_const

  implicit none

  integer, public, parameter :: r8 = SELECTED_REAL_KIND(12)   ! real r8

  !--- circular constant
  real(r8), parameter, public :: pai   = 4.0_r8 * datan(1.0_r8)
  real(r8), parameter, public :: hfpai = 0.5_r8 * pai
  real(r8), parameter, public :: dtr   = pai / 180.0_r8

  !--- Gas constant for dry air
  !real(r8), parameter, public :: Rd = 287.04_r8
  !--- Gas constant for water vapor
  !real(r8), parameter, public :: Rv = 461.52_r8
 ! !--- Gravity
 ! real(r8), parameter, public :: g = 9.80616_r8
 ! !--- Angular speed of the Earth
 ! real(r8), parameter, public :: omega = pai / 43200.0_r8
 ! !--- Radius of the Earth
 ! real(r8), parameter, public :: ae = 2.0d7 / pai
 ! !--- Surface Pressure
 ! real(r8), parameter, public :: ps = 1000.0_r8
  !--- Gravity
  real(r8), parameter, public :: g = 9.356_r8
  !--- Angular speed of the Earth
  real(r8), parameter, public :: omega = 1.2140388d-05
  !--- Radius of the Earth
  real(r8), parameter, public :: ae = 2.0d7 / pai * 15.6
  !--- Surface Pressure
  real(r8), parameter, public :: ps = 1000.0_r8 * 20
  !--- Boltzmann constant
  real(r8), parameter, public :: Boltzmann = 1.380658d-23
  !--- Avogadro constant
  real(r8), parameter, public :: Avogadro  = 6.0221367d+23
  !---
  real(r8), parameter, public :: UnGascon  = Boltzmann*Avogadro
  !!---
  !real(r8), parameter, public :: Md = 28.9644_r8
  !!---
  !real(r8), parameter, public :: Mv = 18.0153_r8
  !---
  real(r8), parameter, public :: Md = 2.36_r8
  !---
  real(r8), parameter, public :: Mv = 18.0153_r8
  !---
  real(r8), parameter, public :: Rd = 1000.0_r8*UnGascon/Md
  !---
  real(r8), parameter, public :: Rv = 1000.0_r8*UnGascon/Mv
  !--- Specific heat for dry air at constant pressure
  real(r8), parameter, public :: cpd = 3.5_r8*Rd
  !--- Specific heat for dry air at constant volume
  real(r8), parameter, public :: cvd = cpd*5/7
  !!--- Specific heat for dry air at constant volume
  !real(r8), parameter, public :: cvd = Cpd-Rd
  !--- cv/cp
  real(r8), parameter, public :: kap = cvd/cpd
  !--- Imaginary
  double complex, parameter, public :: img=(0.0_r8, 1.0_r8)

end module mod_const
