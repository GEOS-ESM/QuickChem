#include "MAPL_Generic.h"

module QC_EnvironmentMod

   use ESMF
   use MAPL

   implicit none
   private

   public :: QC_Environment

   type :: QC_Environment
!      type(Chem_Mie), dimension(2)    :: rad_MieTable, diag_MieTable
!      real, allocatable      :: radius(:)      ! particle effective radius [um]
!      real, allocatable      :: molwght(:)     ! molecular weight            !NOT UNIVERSAL ONLY FOR GASES, 
       integer                :: nbins
       integer                :: km             ! vertical grid dimension
       real                   :: CDT            ! chemistry timestep (secs)
       integer                :: instance       ! data or computational instance
!      real                   :: plid           ! pressure lid [hPa]
       integer                :: klid           ! vertical index of pressure lid
    contains
       procedure :: load_from_config
    end type QC_Environment


    !LOCALS
     integer :: status
     integer :: nbins
!    integer :: n_wavelengths_profile

 contains



    subroutine load_from_config(self, cfg, universal_cfg, rc)
       class(QC_Environment), intent(inout) :: self
       type(ESMF_Config), intent(inout) :: cfg
       type(ESMF_Config), intent(inout) :: universal_cfg
       integer, optional, intent(out) :: rc

       !   Get nbins from cfg
       call ESMF_ConfigGetAttribute (cfg, self%nbins, label='nbins:', __RC__)
       nbins = self%nbins

! Example from universal_cfg:
!      n_wavelengths_profile = ESMF_ConfigGetLen (universal_cfg, label='wavelengths_for_profile_aop_in_nm:', __RC__)  !! from QuickChem_GridComp.rc

! Example of vector allocation and reading into vector:
       !   Parse config file into private internal state
       !   ----------------------------------------------
!      allocate(self%radius(nbins), self%molwght(nbins), __STAT__)
       
!      call ESMF_ConfigGetAttribute (cfg, self%radius,     label='particle_radius_microns:', __RC__)  !! e.g. from XX_instance_XX.rc
!      call ESMF_ConfigGetAttribute (cfg, self%molwght,    label='molecular_weight:', __RC__)         !! e.g. from XX_instance_XX.rc

    end subroutine load_from_config

end module QC_EnvironmentMod
