#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: QuickChem_GridCompMod - The Quick Chemistry Grid Component

! !INTERFACE:

module QuickChem_GridCompMod

! !USES:

   use ESMF
   use MAPL
   use QuickChem_Generic

! !Establish the Childen's SetServices
 !-----------------------------------
   use OH_GridCompMod,    only   : OH_setServices  => SetServices
!  use XX_GridCompMod,    only   : XX_setServices  => SetServices

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices

  ! Private State
  type :: Instance
     integer :: id = -1
     logical :: is_active
     character(:), allocatable :: name
  end type Instance

  type Constituent
     type(Instance), allocatable :: instances(:)
     integer :: n_active
  end type Constituent

  type QuickChem_State
     private
     type(Constituent) :: OH
!    type(Constituent) :: XX
!    real, allocatable :: wavelengths_profile(:) ! wavelengths for profile aop [nm]
!    real, allocatable :: wavelengths_vertint(:) ! wavelengths for vertically integrated aop [nm]
  end type QuickChem_State

  type wrap_
     type (QuickChem_State), pointer     :: PTR => null()
  end type wrap_

! !DESCRIPTION:
!
!   {\tt QuickChem} is a gridded component that includes OH, etc

!
!
! !REVISION HISTORY:
!
!  2022.02.04  Manyin   First crack, adapted from GOCART2G.

!
!EOP
!============================================================================

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices (GC, RC)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                   :: RC  ! return code

! !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!   the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.

! !REVISION HISTORY: 

!  14oct2019  Sherman, da Silva, Darmenov, Clune - First attempt at refactoring for ESMF compatibility


!EOP
!============================================================================
!
!   Locals
    character (len=ESMF_MAXSTR)                   :: COMP_NAME 
    type (ESMF_Config)                            :: myCF      ! QuickChem_GridComp.rc
    type (ESMF_Config)                            :: cf        ! universal config
    type (QuickChem_State), pointer               :: self
    type (wrap_)                                  :: wrap

!   integer :: n_wavelengths_profile, n_wavelengths_vertint, n_wavelengths_diagmie
!   integer, allocatable, dimension(:) :: wavelengths_diagmie

    __Iam__('SetServices')

!****************************************************************************
! Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, name=comp_name, config=cf, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Set the Initialize, Run, Finalize entry points
!   ------------------------------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize,  __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run,  Run1, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run,  Run2, __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState (GC, 'QuickChem_State', wrap, STATUS)
    VERIFY_(STATUS)

    myCF = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (myCF, 'QuickChem_GridComp.rc', __RC__)

!!! OLD GOCART2G section, saved as an example:
!!!
!   Retrieve wavelengths from QuickChem_GridComp.rc
!   n_wavelengths_profile = ESMF_ConfigGetLen (myCF, label='wavelengths_for_profile_aop_in_nm:', __RC__)
!   n_wavelengths_vertint = ESMF_ConfigGetLen (myCF, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)
!   n_wavelengths_diagmie = ESMF_ConfigGetLen (myCF, label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', __RC__)

!   allocate(self%wavelengths_profile(n_wavelengths_profile), self%wavelengths_vertint(n_wavelengths_vertint), &
!            wavelengths_diagmie(n_wavelengths_diagmie), __STAT__)

!   call ESMF_ConfigGetAttribute (myCF, self%wavelengths_profile, label='wavelengths_for_profile_aop_in_nm:', __RC__)
!   call ESMF_ConfigGetAttribute (myCF, self%wavelengths_vertint, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)
!   call ESMF_ConfigGetAttribute (myCF, wavelengths_diagmie, label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', __RC__)

!   Set wavelengths in universal config

!   call MAPL_ConfigSetAttribute (cf, self%wavelengths_profile, label='wavelengths_for_profile_aop_in_nm:', __RC__)
!   call MAPL_ConfigSetAttribute (cf, self%wavelengths_vertint, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)
!   call MAPL_ConfigSetAttribute (cf, wavelengths_diagmie, label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', __RC__)
!!!

!   Get instances for each species
!   -----------------------------------------------------
    call getInstances_('OH', myCF, species=self%OH, __RC__)
!   call getInstances_('XX', myCF, species=self%XX, __RC__)

!   Nitrate currently only supports one active instance
!   if (self%NI%n_active > 1) then
!      if(mapl_am_i_root()) print*,'WARNING: GOCART can only support one active nitrate instance. Check the RC/QuickChem_GridComp.rc'
!   end if

    call ESMF_ConfigDestroy(myCF, __RC__)

!   Create children's gridded components and invoke their SetServices
!   Active instances are created first
!   -----------------------------------------------------------------
    call createInstances_(self, GC, __RC__)

!   Define EXPORT states

!   For an example of exporting a MAPL_StateItem and
!   a MAPL_BundleItem, see GOCART2G_GridCompMod.F90


#include "QuickChem_Export___.h"

!   Allow children of Chemistry to connect to these fields:
    if (size(self%OH%instances) > 0) call MAPL_AddExportSpec (GC, SHORT_NAME='OH', CHILD_ID=self%OH%instances(1)%id, __RC__)

!   For an example of adding connectivities between QuickChem species
!   see GOCART2G_GridCompMod.F90

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!============================================================================
!BOP

! !IROUTINE: Initialize -- Initialize method for the composite Gridded Component

! !INTERFACE:
  subroutine Initialize (GC, import, export, clock, RC)

! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION:  This initializes the QuickChem Grid Component. It primarily creates
!                its exports and establishes its children.

! !REVISION HISTORY: 
! 2022.02.04   Manyin  Adapted from GOCART2G

!EOP
!============================================================================

!   Locals 
    character (len=ESMF_MAXSTR)            :: COMP_NAME

    type (MAPL_MetaComp),       pointer    :: MAPL
    type (ESMF_GridComp),       pointer    :: gcs(:)
    type (ESMF_State),          pointer    :: gex(:)
    type (ESMF_Grid)                       :: grid
    type (ESMF_Config)                     :: CF

    type (QuickChem_State),     pointer    :: self
    type (wrap_)                           :: wrap

    integer                                :: dims(3)

    __Iam__('Initialize')

!****************************************************************************
! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

    if (mapl_am_i_root()) then
       print *, TRIM(Iam)//': Starting...'
       print *,' '
    end if

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

    call MAPL_GridGet ( grid, localCellCountPerDim=dims, __RC__ )

!   Call Generic Initialize
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState (GC, 'QuickChem_State', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    CF = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (CF, 'AGCM.rc', __RC__) ! should the rc file be changed?

!   Get children and their export states from my generic state
!   -----------------------------------------------------------
    call MAPL_Get (MAPL, gcs=gcs, gex=gex, __RC__ )


!   For an example of filling export states with the children's states
!   see GOCART2G_GridCompMod.F90


    RETURN_(ESMF_SUCCESS)


 end subroutine Initialize
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: RUN -- Run method for GOCART2G 


! !INTERFACE:

  subroutine Run1 (GC, import, export, clock, RC)

! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run method 

!EOP
!============================================================================

!   Locals
    character(len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_GridComp),      pointer  :: gcs(:)
    type (ESMF_State),         pointer  :: gim(:)
    type (ESMF_State),         pointer  :: gex(:)
    type (ESMF_State)                   :: internal

    integer                             :: i

    __Iam__('Run1')

!****************************************************************************
! Begin... 


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( MAPL, gcs=gcs, gim=gim, gex=gex, INTERNAL_ESMF_STATE=internal, __RC__ )

!   Run the children
!   -----------------
    do i = 1, size(gcs)
      call ESMF_GridCompRun (gcs(i), importState=gim(i), exportState=gex(i), phase=1, clock=clock, __RC__)
    end do


    RETURN_(ESMF_SUCCESS)

  end subroutine Run1

!============================================================================
!BOP
! !IROUTINE: RUN2 -- Run2 method for GOCART2G component

! !INTERFACE:

  subroutine Run2 (GC, import, export, clock, RC)

! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP
!============================================================================

!   Locals
    character(len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_GridComp),      pointer  :: gcs(:)
    type (ESMF_State),         pointer  :: gim(:)
    type (ESMF_State),         pointer  :: gex(:)
    type (ESMF_State)                   :: internal
    type (QuickChem_State),       pointer  :: self

    type (wrap_)                        :: wrap
    character(len=ESMF_MAXSTR)          :: child_name
    integer                             :: i, n, w
    real, pointer, dimension(:,:)       :: LATS
    real, pointer, dimension(:,:)       :: LONS

#include "QuickChem_DeclarePointer___.h"

    __Iam__('Run2')

!****************************************************************************
! Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( MAPL, gcs=gcs, gim=gim, gex=gex, INTERNAL_ESMF_STATE=internal, &
                    LONS=LONS, LATS=LATS, __RC__ )

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState (GC, 'QuickChem_State', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

#include "QuickChem_GetPointer___.h"


!   Run the children
!   -----------------
    do i = 1, size(gcs)
      call ESMF_GridCompGet (gcs(i), NAME=child_name, __RC__ )
      if ((index(child_name, 'data')) == 0) then ! only execute Run2 method if a computational instance
         call ESMF_GridCompRun (gcs(i), importState=gim(i), exportState=gex(i), phase=2, clock=clock, __RC__)
      end if
    end do

!   in G2G there was a section that computed total aerosol diagnostic values for export

    RETURN_(ESMF_SUCCESS)

  end subroutine Run2


!===============================================================================

  subroutine getInstances_ (aerosol, myCF, species, rc)

!   Description: Fills the QuickChem_State (aka, self%instance_XYZ) with user
!                defined instances from the QuickChem_GridComp.rc.

    implicit none

    character (len=*),                intent(in   )  :: aerosol
    type (ESMF_Config),               intent(inout)  :: myCF
    type(Constituent),                intent(inout)  :: species
    integer,                          intent(  out)  :: rc


!   locals
    integer                                          :: i
    integer                                          :: n_active
    integer                                          :: n_passive
    integer                                          :: n_instances
    character (len=ESMF_MAXSTR)                      :: inst_name

    __Iam__('QuickChem::getInstances_')

!--------------------------------------------------------------------------------------

!   Begin...
    n_active  = ESMF_ConfigGetLen (myCF, label='ACTIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    n_passive = ESMF_ConfigGetLen (myCF, label='PASSIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    n_instances = n_active + n_passive
    allocate (species%instances(n_instances), __STAT__)

!   !Fill the instances list with active instances first
    call ESMF_ConfigFindLabel (myCF, 'ACTIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    do i = 1, n_active
       call ESMF_ConfigGetAttribute (myCF, inst_name, __RC__)
       species%instances(i)%name = inst_name
       species%instances(i)%is_active = .true.
    end do
    species%n_active = n_active

!   !Now fill instances list with passive instances
    call ESMF_ConfigFindLabel (myCF, 'PASSIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    do i = n_active+1, n_active+n_passive
       call ESMF_ConfigGetAttribute (myCF, inst_name, __RC__)
       species%instances(i)%name = inst_name
       species%instances(i)%is_active = .false.
    end do


    RETURN_(ESMF_SUCCESS)

  end subroutine getInstances_


!====================================================================================
  subroutine createInstances_ (self, GC, rc)

!   Description: Creates QuickChem children. Active instances must be created first. If
!     additional QuickChem children are added, this subroutine will need to be updated.

    implicit none

    type (QuickChem_State), pointer,            intent(in   )     :: self
    type (ESMF_GridComp),                    intent(inout)     :: GC
    integer,                                 intent(  out)     :: rc

    ! locals
    integer                                                    :: i

    __Iam__('QuickChem::createInstances_')

!-----------------------------------------------------------------------------------
!   Begin...

!   Active instances must be created first! This ordering is necessary for
!   filing the QC_AERO states that are passed to radiation.
!   This is achieved by arranging the names of the active instances first.

    call addChildren__ (gc, self%OH, setServices=OH_setServices, __RC__)
!   call addChildren__ (gc, self%xx, SETsErvices=XX_setServices, __RC__)

    RETURN_(ESMF_SUCCESS)

    contains
    
        subroutine addChildren__ (gc, species, setServices, rc)
        
          type (ESMF_GridComp),            intent(inout)     :: gc
          type(Constituent),               intent(inout)     :: species
          external                                           :: setServices
          integer,                         intent(  out)     :: rc

          ! local
          integer  :: n

          __Iam__('QuickChem::createInstances_::addChildren__')

          n=size(species%instances)

          do i = 1, n
             species%instances(i)%id = MAPL_AddChild(gc, name=species%instances(i)%name, SS=SetServices, __RC__)
          end do

        RETURN_(ESMF_SUCCESS)

     end subroutine addChildren__

  end subroutine createInstances_

end module QuickChem_GridCompMod
