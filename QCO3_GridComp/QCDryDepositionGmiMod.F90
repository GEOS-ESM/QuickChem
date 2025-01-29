!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!          & Atmospheric Chemistry & Dynamics, Code 614                  !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  QCDryDepositionGmiMod
!
! !INTERFACE:
!

   module  QCDryDepositionGmiMod

! !USES:

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  qcO3_GMI_DryDeposition

!
! !DESCRIPTION:
!
!  This module implements dry deposition; the routines were
!   originally in GMI, under GMI_GridComp/GmiDeposition/
!
! !ORIGINAL COMMENTS:
!
!   CODE DEVELOPER
!     Dan Bergmann, LLNL
!     dbergmann@llnl.gov
!   FILE
!     drydep_update.F
!
! !REVISION HISTORY:
!
!  11Jul2014 - Michael Manyin, first crack
!              This is a modified version of the GMI Dry Dep routine
!              for use by TR
!
!EOP

!-------------------------------------------------------------------------
CONTAINS

!-----------------------------------------------------------------------------
!
! ROUTINE
!   qcO3_GMI_DryDeposition  (previously called Update_Drydep)
!
! DESCRIPTION
!   This routine updates const based on dry deposition using dry deposition
!   velocities calculated by the Harvard dry deposition model.
!
! ARGUMENTS
!   lwis_flags : array of flags that indicate land, water, ice, or snow
!   mcor       : area of grid box (m^2)  [only need if PRESENT(dry_depos)]
!   latdeg     : latitude  (deg)
!   londeg     : longitude (deg)
!   humidity   : specific humidity (g/kg)
!   max_cloud  : maximum overlap cloud fraction for LW
!   radswg     : net downward shortwave radiation at ground (W/m^2)
!   ran_cloud  : random  overlap cloud fraction for LW
!   surf_air_temp : surface air temperature (degK)
!   surf_rough    : surface roughness (m)
!   ustar      : friction velocity (m/s)
!   mass       : total mass of the atmosphere within each grid box (kg) [only need if PRESENT(dry_depos)]
!   dry_depos  : dry deposition accumulated since last output (kg/m^2)
!   diffaer    : aerosol diffusivity at bottom layer (m^2/s)
!   s_radius   : aerosol radius at bottom layer (m)
!   s_velocity : aerosol settling velocity at bottom layer (m/s)
!   one_species: species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine qcO3_GMI_DryDeposition  &
     &  (lwis_flags, mcor, cosSolarZenithAngle, &
     &   fracCloudCover, radswg, surf_air_temp,  &
     &   surf_rough, ustar, mass, &
     &   diffaer,  &
     &   s_radius, s_velocity, BoxHeightEdge, &
     &   ireg, iland, iuse, xlai, &
     &   aero_flag, hstar, f0, delH_298_over_R, &
     &   pr_diag, loc_proc, chem_opt, cdt,    &
     &   i1, i2, ju1, j2, k1, k2, &
     &   mw, one_species, dry_depos)


      implicit none


#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, chem_opt
      real*4 , intent(in) :: cdt
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      real   , intent(in) :: mw
      integer, intent(in) :: lwis_flags          (i1:i2, ju1:j2)
      real*4 , intent(in) :: mcor                (i1:i2, ju1:j2)
      real*8 , intent(in) :: fracCloudCover      (i1:i2, ju1:j2)
      real*8 , intent(in) :: cosSolarZenithAngle (i1:i2, ju1:j2)
      real*8 , intent(in) :: radswg              (i1:i2, ju1:j2)
      real*8 , intent(in) :: surf_air_temp       (i1:i2, ju1:j2)
      real*8 , intent(in) :: surf_rough          (i1:i2, ju1:j2)
      real*8 , intent(in) :: ustar               (i1:i2, ju1:j2)
      real*4 , intent(in) :: mass                (i1:i2, ju1:j2)
      real*4 , intent(in) :: BoxHeightEdge       (i1:i2, ju1:j2)
      real*8 , intent(in) :: diffaer             (i1:i2, ju1:j2)
      real*8 , intent(in) :: s_radius            (i1:i2, ju1:j2)
      real*8 , intent(in) :: s_velocity          (i1:i2, ju1:j2)
      integer, intent(in) :: ireg                (i1:i2, ju1:j2)
      integer, intent(in) :: iland               (i1:i2, ju1:j2, NTYPE)
      integer, intent(in) :: iuse                (i1:i2, ju1:j2, NTYPE)
      real*8 , intent(in) :: xlai                (i1:i2, ju1:j2, NTYPE)
      logical, intent(in) :: aero_flag
      real   , intent(in) :: hstar
      real   , intent(in) :: f0                        ! reactivity factor for oxidation of biological substances
      real   , intent(in) :: delH_298_over_R

      real*8 , intent(inout), OPTIONAL :: dry_depos   (i1:i2, ju1:j2)  !  optional
      real*4 , intent(inout), OPTIONAL :: one_species (i1:i2, ju1:j2)  !  optional

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: GOLDER_A(6) =  &
     &  (/ -0.096d0, -0.037d0, -0.002d0, 1.0d-5,   0.004d0,  0.035d0 /)

      real*8, parameter :: GOLDER_B(6) =  &
     &  (/  0.029d0,  0.029d0,  0.018d0, 0.000d0, -0.018d0, -0.036d0 /)

!     ----------------------
!     Variable declarations.
!     ----------------------

      logical, save :: first = .true.

      integer :: il, ij
      integer :: nsec_jan1

      integer :: pasquill(i1:i2)

      real*8  :: loc_depos(i1:i2)

      real*8  :: box_height(i1:i2)
      real*8  :: csza      (i1:i2)
      real*8  :: cz1       (i1:i2)
      real*8  :: obk       (i1:i2)
      real*8  :: windspeed (i1:i2)

      real*8  :: diffa(i1:i2)
      real*8  :: dvel (i1:i2)
      real*8  :: s_ra (i1:i2)
      real*8  :: s_vel(i1:i2)

      real*8  :: tdt


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'qcO3_GMI_DryDeposition called by ', loc_proc
      end if

      tdt = cdt

!     ==========
!     if (first) then
!     ==========

        first = .false.

        if (PRESENT(dry_depos)) dry_depos(i1:i2,ju1:j2) = 0.0d0

!     end if


      pasquill(:) = 0

      csza(:) = 0.0d0
      dvel(:) = 0.0d0


!cc!? For now just hardwire a few values for the radon/lead problem
!     (chem_opt = 1).  Setkin will eventually provide this information
!     for all species when doing a full chemistry calculation.

!     if (chem_opt == 1) then

!       aerosol(1) = 0
!       aerosol(2) = 1

!       hstar(1)   = 9.3d-3  (weakly soluble gas; since then we made it non-soluble)
!       hstar(2)   = 0.0d0

!       f0(1) = 0.0d0
!       f0(2) = 0.0d0

!     else if (chem_opt == 6) then

!       aerosol(1) = 15   (originally aerosol was index into arrays in gmi_aerosol.h)
!       aerosol(2) = 15

!       hstar(1)   = 0.0d0
!       hstar(2)   = 0.0d0

!       f0(1) = 0.0d0
!       f0(2) = 0.0d0

!     end if

!     ================
      do ij = ju1, j2
!     ================

!       -----------------------------------------------------------
!       Put reasonable values into met fields for testing purposes.
!       -----------------------------------------------------------

        box_height(i1:i2) =  BoxHeightEdge(i1:i2,ij)
        cz1       (i1:i2) =  BoxHeightEdge(i1:i2,ij)*0.50

!       -------------------------------------------------------
!       Calculate the Monin-Obukhov length by first calculating
!       the Pasquill stability class (A-F = 1-6) and then use
!       Golder's formula (see Seinfeld pp 509-511 1986).
!       -------------------------------------------------------

        windspeed(:) =  &
     &    (ustar(i1:i2,ij) / 0.4d0) *  &
     &    Log (10.0d0 / surf_rough(i1:i2,ij))

        do il = i1, i2

          if (cosSolarZenithAngle(il,ij) < -0.05d0) then
            if (fracCloudCover(il,ij) > 0.4375001d0) then
              if      (windspeed(il) < 2.0d0) then
                pasquill(il) = 6
              else if (windspeed(il) < 3.0d0)  then
                pasquill(il) = 5
              else
                pasquill(il) = 4
              end if
            else
              if      (windspeed(il) < 3.0d0) then
                pasquill(il) = 6
              else if (windspeed(il) < 5.0d0) then
                pasquill(il) = 5
              else
                pasquill(il) = 4
              end if
            end if
          else
            if (radswg(il,ij) > 700.0d0) then
              if      (windspeed(il) < 2.5d0) then
                pasquill(il) = 1
              else if (windspeed(il) < 5.0d0) then
                pasquill(il) = 2
              else
                pasquill(il) = 3
              end if
            else if (radswg(il,ij) < 350.0d0) then
              if      (windspeed(il) < 2.0d0) then
                pasquill(il) = 2
              else if (windspeed(il) < 5.0d0) then
                pasquill(il) = 3
              else
                pasquill(il) = 4
              end if
            else
              if      (windspeed(il) < 1.0d0) then
                pasquill(il) = 1
              else if (windspeed(il) < 4.0d0) then
                pasquill(il) = 2
              else if (windspeed(il) < 5.5d0) then
                pasquill(il) = 3
              else
                pasquill(il) = 4
              end if
            end if
          end if

        end do


        obk(:) =  &
     &    1.0d0 /  &
     &    (GOLDER_A(pasquill) +  &
     &     (GOLDER_B(pasquill) * Log (surf_rough(i1:i2,ij))))


!       --------------------------------------------------------
!       Calculate deposition velocities at this latitude for all
!       longitudes and all species.
!       --------------------------------------------------------

        diffa(i1:i2) = diffaer   (i1:i2,ij)
        s_ra (i1:i2) = s_radius  (i1:i2,ij)
        s_vel(i1:i2) = s_velocity(i1:i2,ij)

!       ===========
        call CalcDepositionVelocity  &
!       ===========
     &   (ij, radswg(:,ij), surf_air_temp(:,ij), cosSolarZenithAngle(:, ij), f0,  &
     &    hstar, mw, aero_flag, ustar(:,ij), cz1, obk, fracCloudCover(:,ij),  &
     &    lwis_flags(:,ij), dvel, ireg, iland, iuse, xlai, diffa,  &
     &    s_ra, s_vel, delh_298_over_r, i1, i2, ju1, j2)


!           --------------------------------------------------------
!           Subtract mixing ratio amount from existing mixing ratio
!           of const.
!           --------------------------------------------------------

          if (PRESENT(one_species)) then

            loc_depos(i1:i2)=  &
     &         one_species(i1:i2,ij) *  &
     &         (1.0d0 - Exp (-dvel(i1:i2) * tdt / box_height(i1:i2)))

            one_species(i1:i2,ij) = one_species(i1:i2,ij) - loc_depos(i1:i2)
          end if

          if (PRESENT(dry_depos)) then

!           ------------------------------------------------------
!           Accumulate the total dry deposited in units of kg/m^2.
!           ------------------------------------------------------

              dry_depos(i1:i2,ij) =  &
     &          dry_depos(i1:i2,ij) +  &
     &          loc_depos(i1:i2) * (mass(i1:i2,ij) / mcor(i1:i2,ij)) *  &
     &          mw / MWTAIR

          end if

!     ======
      end do
!     ======


      return

      end subroutine qcO3_GMI_DryDeposition


!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcDepositionVelocity  (previously called Depvel)
!
! DESCRIPTION
!   This routine computes the dry deposition velocities using a
!   resistance-in-series model.
!
!   Routine reads data which:
!     - converts land type id to deposition surface type id
!     - gives roughness heights for each land type id
!     - identifies water land type id's, for stability and z0 calculations
!     - reads surface resistance data for each deposition surface type id
!
!   Changes from version 3.1 to version 3.2:
!     * In unstable atmospheres with |zlmo| < zo, as can happen occasionally
!       under very low wind conditions with tall canopies, application of
!       Monin-Obukhov similarity yields negative values for ra.  This was a
!       problem in version 3.1.  In fact, Monin-Obukhov similarity does not
!       apply under such conditions, so we now set ra to zero and let the
!       boundary resistance rb define the overall aerodynamic resistance.
!       Since rb varies inversely with u*, it will impose a large aerodynamic
!       resistance under very low wind conditions.
!     * The range of applicability of stability correction functions to
!       Monin-Obukhov similarity has been extended to -2.5 < z/zmo < 1.5,
!       based on Figure 2 of Businger et al. [1971].  The range used to be
!       -1 < z/zmo < 1 in version 3.1.
!
!  Literature cited:
!     Baldocchi, D.D., B.B. Hicks, and P. Camara, A canopy stomatal
!       resistance model for gaseous deposition to vegetated surfaces,
!       Atmos. Environ. 21, 91-101, 1987.
!     Brutsaert, W., Evaporation into the Atmosphere, Reidel, 1982.
!     Businger, J.A., et al., Flux-profile relationships in the atmospheric
!       surface layer, J. Atmos. Sci., 28, 181-189, 1971.
!     Dwight, H.B., Tables of integrals and other mathematical data,
!       MacMillan, 1957.
!     Guenther, A., and 15 others, A global model of natural volatile
!       organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!     Hicks, B.B., and P.S. Liss, Transfer of SO2 and other reactive
!       gases across the air-sea interface, Tellus, 28, 348-354, 1976.
!     Jacob, D.J., and S.C. Wofsy, Budgets of reactive nitrogen,
!       hydrocarbons, and ozone over the Amazon forest during the wet season,
!       J.  Geophys. Res., 95, 16737-16754, 1990.
!     Jacob, D.J., and 9 others, Deposition of ozone to tundra,
!       J. Geophys. Res., 97, 16473-16479, 1992.
!     Levine, I.N., Physical Chemistry, 3rd ed., McGraw-Hill, New York, 1988.
!     Munger, J.W., and 8 others, Atmospheric deposition of reactive
!       nitrogen oxides and ozone in a temperate deciduous forest and a
!       sub-arctic woodland, J. Geophys. Res., in press, 1996.
!     Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, SO2, sulfate,
!       and HNO3 deposition velocities computed using regional landuse and
!       meteorological data, Atmos. Environ., 20, 949-964, 1986.
!     Wang, Y.H., paper in preparation, 1996.
!     Wesely, M.L, Improved parameterizations for surface resistance to
!       gaseous dry deposition in regional-scale numerical models,
!       Environmental Protection Agency Report EPA/600/3-88/025,
!       Research Triangle Park (NC), 1988.
!     Wesely, M.L., same title, Atmos. Environ., 23, 1293-1304, 1989.
!
! ARGUMENTS
!   ij      : tbd
!   radiat  : solar radiation (W*m^-2)
!   tempk   : surface air temperature (degK)
!   suncos  : cosine of solar zenith angle
!   f0      : reactivity factor for oxidation of biological substances
!   hstar   : Henry's Law constant
!   mw      : molecular weight; used to calculate molecular diffusivity
!             (kg/mol)
!   aerosol : 0 => gas-phase species, 1 => aerosol species
!   ustar   : friction velocity (m/s)
!   cz1     : altitude at which deposition velocity is computed (m)
!   obk     : Monin-Obukhov length, set to 1.0d5 under neutral conditions (m)
!   cfrac   : fractional cloud cover
!   lwis    : integer for snow and sea ice (0=>water, 1=>land, 2=>ice, 3=>snow)
!             used to be named lsnow
!   dvel    : deposition velocities (m/s)
!   ireg    : # of landtypes in grid square
!   iland   : land type id for elements ldt = 1, ireg;
!             could be from any source, mapped to deposition surface id
!   iuse    : fraction of gridbox area occupied by land type elements (mil^-1)
!   xlai    : leaf area index of land type elements
!   diffa   : aerosol diffusivity (m^2/s)
!   s_ra    : aerosol radius (m)
!   s_vel   : aerosol settling velocity (m/s)
!   delh_298_over_r : temperature scaling term for Henry's law
!
!-----------------------------------------------------------------------------

      subroutine CalcDepositionVelocity  &
     &  (ij, radiat, tempk, suncos, f0, hstar, mw, aero_flag,  &
     &   ustar, cz1, obk, cfrac, lwis, dvel, ireg, iland, iuse, xlai,  &
     &   diffa, s_ra, s_vel, delh_298_over_r, &
     &   i1, i2, ju1, j2)

      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"
#     include "gmi_drydep_data.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2
      integer :: ij
      real*8  :: radiat (i1:i2)
      real*8  :: tempk  (i1:i2)
      real*8  :: suncos (i1:i2)
      real    :: f0
      real    :: hstar
      real    :: mw
      logical :: aero_flag
      real*8  :: ustar  (i1:i2)
      real*8  :: cz1    (i1:i2)
      real*8  :: obk    (i1:i2)
      real*8  :: cfrac  (i1:i2)
      integer :: lwis   (i1:i2)
      real*8  :: dvel   (i1:i2)
      integer :: ireg   (i1:i2, ju1:j2)
      integer :: iland  (i1:i2, ju1:j2, 1:NTYPE)
      integer :: iuse   (i1:i2, ju1:j2, 1:NTYPE)
      real*8  :: xlai   (i1:i2, ju1:j2, 1:NTYPE)
      real*8  :: diffa  (i1:i2)
      real*8  :: s_ra   (i1:i2)
      real*8  :: s_vel  (i1:i2)
      real    :: delh_298_over_r


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: XCKMAN = 0.4d0


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Diffg


!     ----------------------
!     Variable declarations.
!     ----------------------

      logical :: is_found

      logical :: ldep

      logical :: lrgera(i1:i2)  ! T -> stable atmosphere; a high aerodynamic
                                ! resistance (ra = 1.0d4 m/s) is imposed; else
                                ! ra is calculated

      integer :: idep1
      integer :: iloop
      integer :: iolson
      integer :: iw
      integer :: ldt
      integer :: nwaterx

      real*8  :: c1
      real*8  :: ckustr
      real*8  :: corr1
      real*8  :: cz       ! altitude where deposition velocity is computed (m)
      real*8  :: dair
      real*8  :: dummy1, dummy2, dummy3, dummy4
      real*8  :: ra, rb
      real*8  :: rdc
      real*8  :: reyno
      real*8  :: schno
      real*8  :: stono
      real*8  :: ebr, eim, ein, r1
      real*8  :: rix
      real*8  :: rt
      real*8  :: tempc1   ! surface air temperatures in degC
      real*8  :: tempk1   ! surface air temperatures in degK
      real*8  :: xnu
      real*8  :: z0obk

      real*8, save :: press = 1.5d5

      real*8  :: c1x                ! total resistance to deposition for
                                    ! each species (s/m)
      real*8  :: vd
      real*8  :: vk

      real*8  :: henry_eff

      real*8  :: rac (NTYPE)
      real*8  :: rclo(NTYPE)
      real*8  :: rcls(NTYPE)
      real*8  :: rgso(NTYPE)
      real*8  :: rgss(NTYPE)
      real*8  :: ri  (NTYPE)
      real*8  :: rlu (NTYPE)

      real*8  :: zo  (NTYPE)  ! roughness height for specific surface types (m)

      real*8  :: rsurfc(NTYPE)  ! bulk surface resistance for
                                ! each species (s/m)

       character(len=30) :: err_msg

!     ----------------
!     Begin execution.
!     ----------------

!     --------------------------------------------------
!     If ldep==F, species does not deposit.
!     Deposition is applied only to species with ldep==T.
!     --------------------------------------------------

      ldep =  &
     &    ((hstar  > 0.0) .or.  &
     &     (f0     > 0.0) .or.  &
     &     (aero_flag .eqv. .TRUE.))


      dvel(i1:i2) = 0.0d0


!     ============================
      ILOOPX: do iloop = i1, i2
!     ============================

!       ------------------------------
!       Compute deposition velocities.
!       ------------------------------

        cz = cz1(iloop)

        tempk1 = tempk(iloop)
        tempc1 = tempk1 + ABS_ZERO


        vd = 0.0d0
        rsurfc(:) = 0.0d0


!       ---------------------------------------------------------------------
!       Calculate the kinematic viscosity xnu (m^2/s) of air as a function of
!       temperature.  The kinematic viscosity is used to calculate the
!       roughness heights over water surfaces and to diagnose whether such
!       surfaces are aerodynamically rough or smooth using a Reynolds number
!       criterion.  The expression for the temperature dependence of xnu is
!       from the Fortran code in Appendix II of Wesely [1988].
!       ---------------------------------------------------------------------

        c1  = -tempk1 / ABS_ZERO

        xnu = 0.151d0 * (c1**1.77d0) * 1.0d-04


!       -----------------------------------------------------------------
!       Compute bulk surface resistance for gases.
!
!       Adjust external surface resistances for temperature;
!       from Wesely [1989], expression given in text on p. 1296.
!
!       There is no evidence that the resistance continues to increase at
!       temperatures below -18 C, so at colder temperatures, hold the
!       resistance fixed.
!       -----------------------------------------------------------------

   !     rt = 1000.0d0 * Exp (-tempc1 - 4.0d0)

   !     if (tempc1 < -18.0d0) rt = 1.2d9
        ! Change to Zhang et al. 2003 resistance temperature dependence 
        rt = 1.0d0
        if (tempc1 < -1.0d0) rt = Exp (0.2d0*(-1.0d0 - tempc1))
        rt = Min (rt, 2.0d0)


!       ------------------------------------------------------------------
!       Get surface resistances - loop over land types ldt.
!
!       The land types within each grid square are defined using the Olson
!       land type database.  Each of the Olson land types is assigned a
!       corresponding "deposition land type" with characteristic values of
!       surface resistance components.  There are 74 Olson land-types, but
!       only 11 deposition land types (i.e., many of the Olson land types
!       share the same deposition characteristics).  Surface resistance
!       components for the "deposition land types" are from Wesely [1989],
!       except for tropical forests [Jacob and Wofsy, 1990] and for tundra
!       [Jacob et. al., 1992].  All surface resistance components are
!       normalized to a leaf area index of unity.
!
!       Olson land types, deposition land types, and surface resistance
!       components are read from the dry deposition file; check that file
!       for further details.
!       ------------------------------------------------------------------

!       ==========================================
        LDTLOOP1: do ldt = 1, ireg(iloop,ij)
!       ==========================================

!                                      ==============
          if (iuse(iloop,ij,ldt) == 0) cycle LDTLOOP1
!                                      ==============

          if (lwis(iloop) >= 2) then  !  ice=2 or snow=3
            idep1  = 1
          else
            iolson = iland(iloop,ij,ldt) + 1
            idep1  = IDEP(iolson)
          end if


!         =======================
          call Surface_Resistance  &
!         =======================
     &      (idep1, rac(ldt), rclo(ldt), rcls(ldt), rgso(ldt),  &
     &       rgss(ldt), ri(ldt), rlu(ldt), rt, tempc1,  &
     &       cfrac(iloop), radiat(iloop), suncos(iloop),  &
     &       xlai(iloop,ij,ldt), rix)


!         ----------------------------------------------------------------
!         Compute aerodynamic resistance to lower elements in lower part
!         of the canopy or structure, assuming level terrain; equation (5)
!         of Wesely [1989].
!         ----------------------------------------------------------------

          rdc =  &
     &      100.0d0 *  &
     &      (1.0d0 + (1000.0d0 / (radiat(iloop) + 10.0d0)))


!         =============================
          if (ldep) then
!         =============================

            if (aero_flag .eqv. .FALSE.) then

              henry_eff =  &
     &          DBLE(hstar)   *  &
     &          Exp (- DBLE(delh_298_over_r) *  &
     &               ((1.0d0 / tempk1) -  &
     &                (1.0d0 / 298.0d0)))

!             -----------------------------------------------------
!             Species-dependent corrections to resistances are from
!             equations (6)-(9) of Wesely [1989].
!             -----------------------------------------------------

!             ======================
              call Canopy_Resistance  &
!             ======================
     &          (rdc, rix, press, tempk1, DBLE(f0), henry_eff, DBLE(mw),  &
     &           rac(ldt), rclo(ldt), rcls(ldt), rgso(ldt), rgss(ldt),  &
     &           rlu(ldt), rsurfc(ldt))

            else

!             ------------------------------------
!             Get surface resistance for aerosols;
!             Eqs.(5)-(9) of Zhang et al. [2001]
!             ------------------------------------

              schno = xnu / diffa(iloop)

              if (A_VEG(idep1) > 0.0d0) then
                stono = s_vel(iloop) * ustar(iloop) /  &
     &                  (GMI_G * A_VEG(idep1))

                ein = 2.0d0 * (s_ra(iloop) / A_VEG(idep1)) ** 2.0d0
              else
                stono = s_vel(iloop) * ustar(iloop)**2.0d0 /  &
     &                  (GMI_G * xnu)
                ein = 0.0d0
              end if

              ebr = schno ** (-GAMA_BR(idep1))

              eim = (stono / (ALFA_IM(idep1) + stono)) ** 2.0d0

              if (IWET(idep1) > 0) then
                r1 = 1.0d0
              else
                r1 = Exp (-Sqrt (stono))
                if (r1 < 1.0d-10) r1 = 1.0d-10
              end if

              rsurfc(ldt) = 1.0d0 / (3.0d0 * ustar(iloop) *  &
     &                        (ebr + eim + ein) * r1)

            end if

!           ----------------------------------------------------
!           Set min and max values for bulk surface resistances.
!           ----------------------------------------------------

! Updated to Walmsley and Wesley, 1996
            rsurfc(ldt) =  &
     &        Max (1.0d0, Min (rsurfc(ldt), 1.0d25))


!         =============
          end if
!         =============

!       ===============
        end do LDTLOOP1
!       ===============


!       ==========================================
        LDTLOOP3: do ldt = 1, ireg(iloop,ij)
!       ==========================================

!         -------------------------------------------------------
!         Loop through the different landuse types present in the
!         grid square.
!         -------------------------------------------------------

!                                      ==============
          if (iuse(iloop,ij,ldt) == 0) cycle LDTLOOP3
!                                      ==============

          iolson = iland(iloop,ij,ldt) + 1


!         -------------------------------------------------------------
!         Get roughness heights; they are specified constants for each
!         surface type, except over water where zo = f(u*).  The latter
!         dependence is from equation (6) of Hicks and Liss [1976].
!         -------------------------------------------------------------

          nwaterx = NWATER

          is_found = .false.

          IWLOOP: do iw = 1, nwaterx

            if (iolson /= iwater(iw)) cycle IWLOOP

            zo(ldt) =  &
     &        1.4d-02 * ustar(iloop) * ustar(iloop) / 9.8d0 +  &
     &        1.1d-01 * xnu / ustar(iloop)

            is_found = .true.

            exit IWLOOP

          end do IWLOOP

          if (.not. is_found) then
            zo(ldt) = izo(iolson)
            zo(ldt) = zo(ldt) * 1.0d-4
          end if


!         ------------------------------------------------------------------
!         Get aerodynamic resistances ra and rb.
!
!         The aerodynamic resistance ra is integrated from altitude z0+d up
!         to the altitude z1, at which the dry deposition velocity is to be
!         referenced.  The integration corrects for stability using
!         Monin-Obukhov similarity formulas from Businger et. al. [1971],
!         which apply over the range -2.5 < z/zMO < 1.5  (see their
!         Figure 2).  Under very unstable conditions when z1 > -2.5 zMO, we
!         assume that there is no resistance to transfer in the convective
!         column between zmo and z1.  Under very stable conditions when
!         z1 > 1.5 zMO, we assume that vertical transfer in the column
!         between zmo and z1 is strongly suppressed so that the deposition
!         velocity at altitude z1 is very low.  Under these conditions, we
!         just specify a very large ra=1.0d4 s/m (lrgera = T).
!
!         The Reynolds number reyno diagnoses whether a surface is
!         aerodynamically rough (reyno > 10) or smooth.
!         NOTE: The criterion "reyno > 10" is now replaced by "reyno > 1".

!         Surface is rough in all cases except over water with low wind
!         speeds.  In the smooth case, vertical transport IN THE SUBLAYER
!         near the surface is limited by molecular diffusion and is
!         therefore very slow; we assign a large value of ra + rb to
!         account for this effect. (djj, hyl, bmy, 5/8/00)
!
!         In the aerodynamically rough case, the expression for ra is as
!         given in equation (5) of Jacob et. al. [1992]:
!           ra = (1/ku*)*int(from z0 to z1) (phi(x)/z)dz
!         where x = (z-d)/zmo, z is the height above ground, and d is the
!         displacement height, which is typically 70-80% of the canopy
!         height [Brutsaert, 1982].  We change the vertical coordinate so
!         that z=0 at the displacement height; that's OK since for all
!         practical applications, z1 >> d.  In this manner, we don't need to
!         assume any specific value for the displacement height.  Applying
!         the variable transformation z -> x = z/zmo, the equation above
!         becomes:
!           ra = (1/ku*)*int(from x0 to x1) (phi(x)/x)dx   with x=z/zmo
!         Here phi is a stability correction function originally formulated
!         by Businger et. al. [1971] and given in eqns 5a and 5b of Jacob
!         et. al. [1992].
!         For unstable conditions:
!           phi(x) = a/sqrt(1-bx)  where a=0.74, b = 9
!         The analytical solution to the integral is [Dwight, 1957, integral
!         192.11]:
!           int(dx/(x*sqrt(1-bx))) = log(abs((sqrt(1-bx)-1)/(sqrt(1-bx)+1)))
!         which yields the expression for ra used in the code for unstable
!         conditions.
!         For stable conditions,
!           phi(x) = a + bx        where a=0.74, b = 4.7
!         and the analytical solution to the integral is:
!           int((a/x)+b)dx = a*ln(x) + bx
!         which yields the expression of ra used in the code for stable
!         conditions.
!
!         The formulation of rb for gases is equation (12) of Walcek et. al.
!         [1986].
!         ------------------------------------------------------------------

          ckustr = XCKMAN * ustar(iloop)

          reyno = ustar(iloop) * zo(ldt) / xnu

          if (obk(iloop) == 0.0d0) then

            Write (6,910) obk(iloop), iloop, ldt

 910        format (1x, 'obk(iloop) = ', e11.2, 1x, 'iloop = ',  &
     &              i4, 1x, 'ldt = ', i3, /)

          end if

          corr1 = cz / obk(iloop)

          z0obk = zo(ldt) / obk(iloop)

          lrgera(iloop) = .false.

          if (corr1 > 0.0d0) then
            if (corr1 >   1.5d0) lrgera(iloop) = .true.
          else if (corr1 <= 0.0d0) then
            if (corr1 <= -2.5d0) corr1 = -2.5d0
          end if

          if (ckustr == 0.0d0) then

            Write (6,920) iloop, ckustr, XCKMAN, ustar(iloop)

 920        format (1x, 'iloop = ', i4, 1x, 'ckustr = ', e10.1, 1x,  &
     &              'XCKMAN = ', e12.4, 1x, 'ustar(iloop) = ', e12.4)

!c?
            Close (98)

            err_msg = ' Problem in Depvel '
            call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

          end if


!         ------------------------------------------------------------------
!         Aerodynamically rough or smooth surface.
!
!         In the classic study by Nikuradse (1933) the transition from
!         smooth to rough was examined in pipe flow.  He introduced a
!         roughness Reynolds number rr = u* z0 / nu, and found the flow to
!         be smooth for rr < 0.13, and rough for rr > 2.5, with a transition
!         regime in between (E.B. Kraus and J.A. Businger, Atmosphere-Ocean
!         Interaction, second edition, p.144-145, 1994).  Similar statements
!         can be found in the books:
!           Evaporation into the Atmosphere, by Wilfried Brutsaert, p.59,89,
!             1982; or
!           Seinfeld & Pandis, p.858, 1998.
!         Here we assume a sudden transition point rr = 1 from smooth to
!         rough, following L. Merlivat (1978, The dependence of bulk
!         evaporation coefficients on air-water interfacial conditions as
!         determined by the isotopic method, J. Geophys. Res.,
!         Oceans & Atmos., 83, C6, 2977-2980).  Also refer to Brutsaert's
!         book, p.125.  we used to use the criterion "reyno > 10" for
!         aerodynamically rough surface and now change to "reyno > 1".
!         (hyl, 10/15/99)
!         ------------------------------------------------------------------

!         ===================
          if (reyno >= 1.0d0) then
!         ===================

!           ------------------------------
!           Aerodynamically rough surface.
!           ------------------------------

            if ((corr1 <= 0.0d0) .and. (z0obk < -1.0d0)) then

!             ----------------------------------------
!             Unstable condition with Abs (zlmo) < Z0;
!             set ra to zero (version 3.2).
!             ----------------------------------------

              ra = 0.0d0

            else if ((corr1 <= 0.0d0) .and. (z0obk >= -1.0d0)) then

!             ---------------------------------------------------
!             Unstable conditions; compute ra as described above.
!             ---------------------------------------------------

              dummy1 = (1.0d0 - (9.0d0 * corr1))**0.5d0
              dummy2 = (1.0d0 - (9.0d0 * z0obk))**0.5d0
              dummy3 = Abs ((dummy1 - 1.0d0) / (dummy1 + 1.0d0))
              dummy4 = Abs ((dummy2 - 1.0d0) / (dummy2 + 1.0d0))

              ra = 0.74d0* (1.0d0 / ckustr) * Log (dummy3 / dummy4)

            else if ((corr1 > 0.0d0) .and.  &
     &               (.not. lrgera(iloop))) then

!             -----------------------------------------
!             Moderately stable conditions (z/zmo < 1);
!             compute ra as described above
!             -----------------------------------------

              ra =  &
     &          (1.0d0 / ckustr) *  &
     &          (0.74d0 * Log (corr1 / z0obk) +  &
     &           4.7d0  * (corr1 - z0obk))

            else if (lrgera(iloop)) then

!             -----------------------
!             Very stable conditions.
!             -----------------------

              ra = 1.0d+04

            end if


            if (ra < 0.0d0) then

!             ------------------------------------------------------------
!             If ra is negative (as occasionally happened in version 3.1),
!             send a warning message.
!             ------------------------------------------------------------

              Write (6,930) ra, cz, zo(ldt), obk(iloop)

 930          format ('WARNING:  ra < 0 in Depvel', 4(1x, e12.5))

            end if


!             ------------------------------------
!             Get total resistance for deposition.
!             ------------------------------------

            if (ldep) then

              if (aero_flag .eqv. .FALSE.) then

!               ------------------------------------------------------------
!               dair is the thermal diffusivity of air; value of
!               0.2*1.E-4 m^2/s, cited on p. 16,476 of Jacob et. al. [1992].
!               ------------------------------------------------------------

                dair = 0.2d0 * 1.0d-4

                rb = (2.0d0 / ckustr) *  &
     &            (dair / Diffg (tempk1, press, mw*KGPG))**0.667d0
                c1x = ra + rb + rsurfc(ldt)

              else

                c1x = ra + rsurfc(ldt)

              end if

            end if

!         ====
          else
!         ====

!           -------------------------------
!           Aerodynamically smooth surface.
!           -------------------------------

            if (ldep) then

!             --------------------------------------------------------------
!             Suppress drydep over smooth surfaces by setting ra to a
!             large value.  This prevents negative dry deposition
!             velocities when ustar is very small (djj, bmy, 5/8/00)
!             --------------------------------------------------------------

              ra = 1.0d4

              c1x = ra + rsurfc(ldt)

            end if

!         ======
          end if
!         ======


!         ----------------------------------------------------------------
!         iuse is the fraction of the grid square occupied by surface ldt
!         in units of mil^-1 (iuse=500 => 50% of the grid square).  Add
!         the contribution of surface type ldt to the deposition velocity;
!         this is a loop over all surface types in the gridbox.
!         ----------------------------------------------------------------

          if (ldep) then

            vk = vd
            vd = iuse(iloop,ij,ldt)

            if (aero_flag .eqv. .FALSE.) then
              vd = vk + (0.001d0 * vd / c1x)
            else
              vd = vk + (0.001d0 * vd * (1./c1x +  &
     &          s_vel(iloop)))
            end if

          end if


!       ===============
        end do LDTLOOP3
!       ===============


        if (ldep) then

          dvel(iloop) = vd

        end if


!     ==============
      end do ILOOPX
!     ==============


      return

      end subroutine CalcDepositionVelocity 

   end module QCDryDepositionGmiMod
