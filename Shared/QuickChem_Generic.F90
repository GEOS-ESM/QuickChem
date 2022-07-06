
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!

! !MODULE: QuickChem_Generic - Utilitarian subroutines used by QuickChem children. 
!                             
!
! !INTERFACE:
!
module  QuickChem_Generic

! !USES:
   use ESMF
   use MAPL

   implicit none
   private

!
! !PUBLIC MEMBER FUNCTIONS:

!!   public add_aero           -- example in G2G of adding a State for Export
!!   public append_to_bundle   -- example in G2G of adding a field to a field bundle
!                                 NOTE that there are SPECIES-SPECIFIC clauses in this
!!   public get_mixR           -- example in G2G of simple routine to sum over instances
!                                 NOTE that there are SPECIES-SPECIFIC clauses in this

   public determine_data_driven

!
! !DESCRIPTION:
!
!  These subroutines perform repetitive tasks needed by QuickChem children.
!
! !REVISION HISTORY:
!
!  March2020 Sherman, da Silva, Darmenov, Clune - created
!
!EOP
!-------------------------------------------------------------------------
contains

!=====================================================================================

  subroutine determine_data_driven(COMP_NAME, data_driven, RC)

    !ARGUMENTS:
    integer, optional,               intent(  out)   :: RC          ! Error code:
    character (len=ESMF_MAXSTR),     intent(in   )   :: COMP_NAME
    logical,                         intent(  out)   :: data_driven

    !Local
    integer                                          :: i

!   Description: Determines whether gridded component is data driven or not.

     __Iam__('determine_data_driven')

!   Begin... 

!   Is DU data driven?
!   ------------------
    data_driven = .false.

    i = index(COMP_NAME, 'data')
    if (i > 0) then
      data_driven = .true.
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine determine_data_driven

end module  QuickChem_Generic


