#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: OH_GridCompMod - Hydroxyl Radical gridded component 
! !                         (Member of QuickChem)

! !INTERFACE:
module OH_GridCompMod

!  !USES:
   use ESMF
   use MAPL
   use QuickChem_Generic
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr, c_int64_t, c_int, c_float

   use QC_EnvironmentMod
   use MAPL_StringTemplate

   use  xgb_fortran_api
   
   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

! Possible values for OH_data_source:
   INTEGER, parameter :: DS_UNDEFINED = 99
   INTEGER, parameter :: PRECOMPUTED  =  1
   INTEGER, parameter :: ONLINE_INST  =  2
   INTEGER, parameter :: ONLINE_AVG24 =  3

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices

! !DESCRIPTION: This module implements QuickChem's template (OH) Gridded Component.

! !REVISION HISTORY:
!  February 2022: M Manyin ported OH from Icarus version
!  July     2022: M Manyin made QuickChem a separate repository, and optimized OH

!EOP
!===========================================================================

!  !OH state
   TYPE, EXTENDS(QC_Environment) :: OH_GridComp
!      real, allocatable      :: rlow(:)        ! example of contents, was used in Run2
!      character(len=:), allocatable :: emission_scheme     ! See G2G DU for how this is used
!                                                           ! NOTE: it affects OH_Import___.h and OH_GetPointer___.h
       CHARACTER(LEN=ESMF_MAXPATHLEN) :: XGBoostFilePattern

! Some fields have optional sources
       INTEGER :: OH_data_source

!  If true,  OH imports both the instantaneous and 24hr average fields for T, P, etc
!            This allows for the use of instantaneous fields when tavg fields are all zero.
!  If false, OH imports only the 24hr average fields; NOTE - they must be in the OH_import_rst
       LOGICAL :: spinup_24hr_imports

!  If true, the model is within the first 24 hours of integration, daily means are not available
       LOGICAL :: use_inst_values

!  If true, only call Boost during the first timestep of the day
       LOGICAL :: compute_once_per_day

!  Get Scattering Coefficient values from GOCART2G
       INTEGER :: n_wavelengths_profile   ! Size of the 4th dimension for SCACOEF fields
       REAL    :: wavelength_for_scacoef  ! Wavelength value to be used
       INTEGER :: wavelength_index        ! Index of the desired value within the 4th dimension

!  Store the most recent version of OH from Boost in this variable,
!  which will persist from one timestep to the next:
       REAL, ALLOCATABLE         ::   OH_ML(:,:,:) ! top-down  (gather the answer in this)

   END TYPE OH_GridComp

   TYPE OH_BOOST_INPUT_DATA
!      Datasets used as input for the XGBoost model of OH, as specified by Dan Anderson
!        M2G    ==  from MERRA2-GMI
!        ONLINE ==  ONLINE_instantaneous, ONLINE_daily_mean, or PRECOMPUTED (M2G)
!        To save on unnecessary allocation, all entries are pointers
       REAL, POINTER, DIMENSION(:,:  ) ::  LAT                ! from the GEOS grid
       REAL, POINTER, DIMENSION(:,:,:) ::  PL                 ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  T                  ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  NO2                ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  O3                 ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  CH4                ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  CO                 ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  ISOP               ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  ACET               ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  C2H6               ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  C3H8               ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  PRPE               ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  ALK4               ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  MP                 ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  H2O2               ! M2G
       REAL, POINTER, DIMENSION(:,:,:) ::  TAUCLWDN           ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  TAUCLIDN           ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  TAUCLIUP           ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  TAUCLWUP           ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  CLOUD              ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  QV                 ! ONLINE
       REAL, POINTER, DIMENSION(:,:)   ::  GMISTRATO3         ! M2G
       REAL, POINTER, DIMENSION(:,:)   ::  ALBUV              ! from a climatology
       REAL, POINTER, DIMENSION(:,:,:) ::  AODUP              ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  AODDN              ! ONLINE
       REAL, POINTER, DIMENSION(:,:,:) ::  CH2O               ! M2G
       REAL, POINTER, DIMENSION(:,:)   ::  SZA                ! computed for local noon
   END TYPE OH_BOOST_INPUT_DATA

   type wrap_
      type (OH_GridComp), pointer     :: PTR => null()
   end type wrap_

contains

!-------------------------------------------------------------------------
   subroutine predict_OH_with_XGB( xgb_fname, icount, jcount, kcount,  &
                                   dynamic_k_range, tropp_min, pl, tropp, bb, OH_ML, rc)

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

   character(len=*),          intent(IN   ) :: xgb_fname
   integer,                   intent(IN   ) :: icount, jcount, kcount   ! assume all indices start at 1
   LOGICAL,                   intent(IN   ) :: dynamic_k_range  ! if TRUE, restrict slab using pl and tropp
   REAL,                      intent(IN   ) :: tropp_min    ! Pa - if NOT dynamic_k_range, use this min pressure
   REAL, ALLOCATABLE,         intent(IN   ) :: pl(:,:,:)    ! Pa - pressure at midlevel
   REAL, POINTER,             intent(IN   ) :: tropp(:,:)   ! Pa - tropopause pressure
   TYPE(OH_BOOST_INPUT_DATA), intent(IN   ) :: bb        ! XGBoost input data pointers
   REAL,                      intent(INOUT) :: OH_ML(:,:,:) ! mol/mol  top-down  (gather the answer in this)
   integer,                   intent(  OUT) :: rc

! !DESCRIPTION: Predict OH using XGBoost
! !             Based on original code from C. Keller
!
! !REVISION HISTORY:
! !             2021.10.07  Manyin - first crack
!
!
!EOP
!-------------------------------------------------------------------------


   CHARACTER(LEN=255) :: rcbasen = 'OH_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

! Calling XGBoost using the wrapper library of Christoph Keller that supports FORTRAN -> C

! 2022-07-22 - manyin - based on the sample code that Christoph provides with the library
! ----------------------------------------------------------------------------
!use iso_c_binding
!use xgb_fortran_api
!implicit none

!--- Local variables
 integer(c_int64_t)         :: xx_prediction_count  ! how many gridboxes to predict - varies
                                                    ! used in INIT and RUN
 ! for booster object 
 integer(c_int64_t)         :: xx_param_count       ! hard coded - number of params using in training
                                                    ! used in INIT and RUN
 ! for booster object 
 character(len=255)         :: xx_fname_bst         ! set to xgb_fname
                                                    ! used in INIT

 integer(c_int64_t)         :: xx_dmtrx_len         ! hard coded below
                                                    ! used in INIT

 type(c_ptr),          save :: xx_bst               ! INIT   - create and load
                                                    ! RUN    - use
                                                    ! FINAL  - free
 ! for XGDMatrix object
 type(c_ptr)                :: xx_dmtrx             ! INIT - create, use and free
                                                    ! RUN  - create, use and free

!integer(c_int64_t)         :: xx_dm_nrow, xx_dm_ncol   ! INIT - just check dims

 real(c_float), allocatable :: xx_carr(:,:)         ! allocated every time in RUN
 real(c_float), allocatable :: xx_carr_small(:,:)   ! allocated every time in INIT

 ! for prediction
 integer(c_int)             :: xx_option_mask, xx_ntree_limit, xx_training  ! RUN

 type(c_ptr)                :: xx_cpred  ! RUN

 integer(c_int64_t)         :: xx_pred_len   ! RUN

 real(c_float), pointer     :: xx_pred(:)  ! RUN

 integer :: i,j,k
 integer :: ksubcount   ! the vertical dimension of the slab sent to Boost
 integer :: k1,k2       ! the levels where predictions are needed

 integer :: m           ! index into xx_carr

 logical,              save :: first_time    = .TRUE.

!--- Parameters
 ! missing value 
 real(c_float), parameter :: xx_miss = -999.0   ! used in INIT and RUN

 ! debug flag
 logical, parameter       :: debug    = .false.

   Iam = "predict_OH_with_XGB"

!--- Settings

 ! Dan's OH model
 ! Dimension of input array. We assume it's 1x27 (1 prediction, 27 input variables)
 !xx_fname_bst = '/discover/nobackup/dcanders/QuickChem/Data/RTModels/xgboh_UpDwnALBUVSZAAll_NoGMIALB_NoScale_WithCH4_Depth18_eta1_100Trees_M07.joblib.dat'
 !xx_fname_bst = '/discover/nobackup/mmanyin/CCM/run/oh26/OH_model.bin'
 !xx_fname_bst = '/discover/nobackup/mmanyin/CCM/run/oh26/OH_model/OH_NoTune.bin'
 xx_fname_bst = TRIM(xgb_fname)
 xx_param_count      = 27     !  This entire file would need to change significantly, if the parameter count were to change

    ! Settings for prediction 
    xx_option_mask = 0
    xx_ntree_limit = 0    ! Per Dan: "The final version of the parameterization has 100 trees"
                          ! Tried value 100 and 0, apparently zero-diff
                          ! 100 seemed to take longer
    xx_training    = 0    ! actually FALSE or TRUE - new for v1.6.0


    if(debug) write(*,*) 'Starting predict_OH_with_XGB'

!!!!!!!!!!!!!!!!!!!!!!!!!   INIT

       ONE_TIME_SETUP: IF ( first_time ) THEN


         allocate(xx_carr_small(xx_param_count,1),__STAT__)   ! Just 1 prediction
         xx_carr_small(:,:) = 0.0

         !--- Load booster object in 3 steps

         ! Create (dummy) XGDMatrix
         rc = XGDMatrixCreateFromMat_f( xx_carr_small, 1, xx_param_count, xx_miss, xx_dmtrx )
         _ASSERT(rc==0,'Failed in XGDMatrixCreateFromMat_f')

         ! Create XGBooster object
         xx_dmtrx_len = 0    ! setting this to 27 results in a Seg Fault
         rc = XGBoosterCreate_f( xx_dmtrx, xx_dmtrx_len, xx_bst )
         _ASSERT(rc==0,'Failed in XGBoosterCreate_f')

         ! load XGBooster model from binary file
         if(debug) write(*,*) 'Reading '//trim(xx_fname_bst)
         rc = XGBoosterLoadModel_f( xx_bst, xx_fname_bst )
         _ASSERT(rc==0,'Failed in XGBoosterLoadModel_f')

         rc = XGDMatrixFree_f( xx_dmtrx )
         _ASSERT(rc==0,'Failed in XGDMatrixFree_f')

         if (allocated (xx_carr_small) ) deallocate(xx_carr_small)

         first_time  = .FALSE.

       END IF  ONE_TIME_SETUP

!!!!!!!!!!!!!!!!!!!!!!!!!   RUN

       IF ( dynamic_k_range ) THEN
         ! Determine the subset of K we need to work on:
         ! check the count in each column, and take the max
         ksubcount=0
         DO j=1,jcount
         DO i=1,icount
           k = COUNT( pl(i,j,:) .GT. tropp(i,j) )   ! Note the same comparison at @@@
           IF ( k > ksubcount ) ksubcount = k
         END DO
         END DO
       ELSE
         ! double check that the min pressure is sufficiently low
         k = COUNT( tropp .LE. tropp_min )
         _ASSERT(k==0,'OH Prediction: Minimum tropopause pressure is not low enough!')

         ! As above, but compare to a constant
         ksubcount=0
         DO j=1,jcount
         DO i=1,icount
           k = COUNT( pl(i,j,:) .GT. tropp_min )   ! Note the same comparison at @@@
           IF ( k > ksubcount ) ksubcount = k
         END DO
         END DO
       END IF

       k1 = kcount - ksubcount + 1
       k2 = kcount

       ! Array to hold input values
       ! (reverse the regular row/col order because we are calling C)
       xx_prediction_count  = icount * jcount * ksubcount
       allocate(xx_carr(xx_param_count,xx_prediction_count),__STAT__)

       m=1
       DO k=k1,k2
       DO j= 1,jcount
       DO i= 1,icount

         xx_carr( 1,m) = bb%LAT(       i,j  )          ! 2D
         xx_carr( 2,m) = bb%PL(        i,j,k) / 100.0  ! convert Pa to hPa
         xx_carr( 3,m) = bb%T(         i,j,k)
         xx_carr( 4,m) = bb%NO2(       i,j,k)
         xx_carr( 5,m) = bb%O3(        i,j,k)
         xx_carr( 6,m) = bb%CH4(       i,j,k)
         xx_carr( 7,m) = bb%CO(        i,j,k)
         xx_carr( 8,m) = bb%ISOP(      i,j,k)
         xx_carr( 9,m) = bb%ACET(      i,j,k)
         xx_carr(10,m) = bb%C2H6(      i,j,k)
         xx_carr(11,m) = bb%C3H8(      i,j,k)
         xx_carr(12,m) = bb%PRPE(      i,j,k)
         xx_carr(13,m) = bb%ALK4(      i,j,k)
         xx_carr(14,m) = bb%MP(        i,j,k)
         xx_carr(15,m) = bb%H2O2(      i,j,k)
         xx_carr(16,m) = bb%TAUCLWDN(  i,j,k)
         xx_carr(17,m) = bb%TAUCLIDN(  i,j,k)
         xx_carr(18,m) = bb%TAUCLIUP(  i,j,k)
         xx_carr(19,m) = bb%TAUCLWUP(  i,j,k)
         xx_carr(20,m) = bb%CLOUD(     i,j,k)
         xx_carr(21,m) = bb%QV(        i,j,k)
         xx_carr(22,m) = bb%GMISTRATO3(i,j  )          ! 2D
         xx_carr(23,m) = bb%ALBUV(     i,j  )          ! 2D
         xx_carr(24,m) = bb%AODUP(     i,j,k)
         xx_carr(25,m) = bb%AODDN(     i,j,k)
         xx_carr(26,m) = bb%CH2O(      i,j,k)
         xx_carr(27,m) = bb%SZA(       i,j  )          ! 2D

         m = m + 1

       END DO
       END DO
       END DO

       rc = XGDMatrixCreateFromMat_f(xx_carr, xx_prediction_count, xx_param_count, xx_miss, xx_dmtrx)
!      if(rc /= 0) then
!       DO m=1,xx_param_count
!         print*,'  Failing - xx_carr ', m, xx_carr(:,m)   REVISE THIS TO PRINT LESS
!       end do
!      endif
       _ASSERT(rc==0,'Failed in XGDMatrixCreateFromMat_f')

       ! Make prediction. The result will be stored in c pointer xx_cpred 
       rc = XGBoosterPredict_f(xx_bst,xx_dmtrx,xx_option_mask,xx_ntree_limit,xx_training,xx_pred_len,xx_cpred)
       !if (mapl_am_i_root()) print*,'BOOST PREDICTION LENGTH:',xx_pred_len   --->  72 when doing single columns
       _ASSERT(rc==0,'Failed in XGBoosterPredict_f')
       _ASSERT(xx_pred_len==xx_prediction_count,'Wrong value returned for xx_pred_len')

       ! Link to fortran pointer xx_pred 
       call c_f_pointer(xx_cpred, xx_pred, [xx_pred_len])

       m=1
       DO k=k1,k2
       DO j= 1,jcount
       DO i= 1,icount

         OH_ML(i,j,k) = 10.0 ** (xx_pred(m))   ! yields mol/mol
         m = m + 1

       END DO
       END DO
       END DO


       rc = XGDMatrixFree_f(xx_dmtrx)
       _ASSERT(rc==0,'Failed in XGDMatrixFree_f')


       if (associated(xx_pred) ) nullify(xx_pred)

       if (allocated (xx_carr) ) deallocate(xx_carr)

!!!!!!!!!!!!!!!!!!!!!!!!!   FINAL

!--- Cleanup

!!     Only free at the end of the program
!      rc = XGBoosterFree_f(xx_bst)
!      if(rc /= 0) write(*,*) __FILE__,__LINE__,'Failed in XGBoosterFree_f, Return code: ',rc
!      VERIFY_(rc)

      if(debug) write(*,*) 'Finished predict_OH_with_XGB'

      RETURN_(ESMF_SUCCESS)

   end subroutine predict_OH_with_XGB


   function computeSolarZenithAngle_LocalNoon (jday, latRad, lonRad, i1, i2, j1, j2) &
               result(this_)
!
! !INPUT PARAMETERS:
      integer :: i1, i2, j1, j2
      integer :: jday             ! day of year (1-366)
      real    :: latRad(i1:i2,j1:j2)    ! latitude (radians)
      real    :: lonRad(i1:i2,j1:j2)    ! longitude (radians)
!
! !RETURNED VALUE
      real    :: this_(i1:i2,j1:j2)
!
! !DESCRIPTION:
!  Computes the solar zenith angle for local noon
!  Based on the GMI function computeSolarZenithAngle_Photolysis
!
! !LOCAL VARIABLES:
    REAL    :: sindec, soldek, cosdec
    REAL    :: sinlat(i1:i2,j1:j2), sollat(i1:i2,j1:j2), coslat(i1:i2,j1:j2)
    REAL    :: cosz(i1:i2,j1:j2)
    REAL    :: tau, loct
    REAL    :: mylon
    integer :: i, j
!EOP
!------------------------------------------------------------------------------
!BOC
      sindec=0.3978*sin(0.9863*(jday-80.0)*MAPL_DEGREES_TO_RADIANS)
      soldek=asin(sindec)
      cosdec=cos(soldek)
      sinlat(i1:i2,j1:j2)=sin(latRad(i1:i2,j1:j2))
      sollat(i1:i2,j1:j2)=asin(sinlat(i1:i2,j1:j2))
      coslat(i1:i2,j1:j2)=cos(sollat(i1:i2,j1:j2))

      do j = j1, j2
       do i = i1, i2
! tau was originally "hours UTC" for model time
!     now we make it hours UTC for local noon
!     Greenwich uses 12 of course
         mylon = lonRad(i,j)*MAPL_RADIANS_TO_DEGREES
         ! We want  -180 < mylon < 180
         if ( mylon .GT.  180.0 ) mylon=mylon-360.0
         if ( mylon .LT. -180.0 ) mylon=mylon+360.0
         tau = 12.0 + (mylon/-180.0)*12.0

         loct       = ((tau*15.0)-180.0)*MAPL_DEGREES_TO_RADIANS + lonRad(i,j)
         cosz(i,j)  = cosdec*coslat(i,j)*cos(loct) + sindec*sinlat(i,j)

         IF ( cosz(i,j) < -1.0 ) THEN
           PRINT*,'COS LESS THAN -1 by:', cosz(i,j) + 1.0
           PRINT*,'<    DAY LAT LON:', jday, latRad(i,j)*MAPL_RADIANS_TO_DEGREES, lonRad(i,j)*MAPL_RADIANS_TO_DEGREES
         ENDIF

         IF ( cosz(i,j) >  1.0 ) THEN
           PRINT*,'COS GREATER THAN 1 by:', cosz(i,j) - 1.0
           PRINT*,'>    DAY LAT LON:', jday, latRad(i,j)*MAPL_RADIANS_TO_DEGREES, lonRad(i,j)*MAPL_RADIANS_TO_DEGREES
         ENDIF

! cosz was above 1 by 1e-7 on occasion
         cosz(i,j) = MIN(  1.0, cosz(i,j) )
         cosz(i,j) = MAX( -1.0, cosz(i,j) )

         this_(i,j) = acos(cosz(i,j))*MAPL_RADIANS_TO_DEGREES
       enddo
      enddo

   end function computeSolarZenithAngle_LocalNoon


!============================================================================
!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices (GC, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code

!   !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!     the Initialize and Finalize services to generic versions. It also
!     allocates our instance of a generic state and puts it in the 
!     gridded component (GC). Here we only set the two-stage run method
!     and declare the data services.

!   !REVISION HISTORY: 
!    Feb 2022 - M Manyin adapted from GOCART2G species

!EOP
!============================================================================

!
!   !Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (ESMF_Config)                 :: cfg
    type (ESMF_Config)                 :: universal_cfg
    type (wrap_)                       :: wrap
    type (OH_GridComp), pointer        :: self

    character (len=ESMF_MAXSTR)        :: field_name
    integer                            :: i
    logical                            :: data_driven = .true.

    type (ESMF_Config)                 :: g2g_cfg    ! We need values from GOCART2G
    real, allocatable                  :: wavelengths_profile(:) ! wavelengths for profile aop [nm]

    CHARACTER(LEN=255) :: token

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, config=universal_cfg, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self


! search for the desired wavelength in the array:


!
!   Load OH resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'OH_instance_'//trim(COMP_NAME)//'.rc', rc=status)

    if (status /= 0) then
        if (mapl_am_i_root()) print*,'OH_instance_'//trim(COMP_NAME)//'.rc does not exist! Loading OH_GridComp_OH.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'OH_instance_OH.rc', __RC__)
    end if

    ! process generic config items
    call self%QC_Environment%load_from_config(cfg, universal_cfg, __RC__)

!   allocate( self%rlow(self%nbins),  __STAT__)  !! example

    ! process DU-specific items
!   call ESMF_ConfigGetAttribute (cfg, self%rlow,       label='radius_lower:', __RC__)   !! example of reading a vector of values from OH_instance_OH.rc
    call ESMF_ConfigGetAttribute(cfg, self%XGBoostFilePattern, label='XGBoostFile:', __RC__)

    call ESMF_ConfigGetAttribute(cfg, token, label='OH_data_source:', __RC__)
    self%OH_data_source = DS_UNDEFINED
    IF ( TRIM(token) == "PRECOMPUTED"  ) self%OH_data_source = PRECOMPUTED
    IF ( TRIM(token) == "ONLINE_INST"  ) self%OH_data_source = ONLINE_INST
    IF ( TRIM(token) == "ONLINE_AVG24" ) self%OH_data_source = ONLINE_AVG24
    IF ( self%OH_data_source == DS_UNDEFINED ) THEN
      IF (mapl_am_i_root()) PRINT*,'Invalid OH_data_source: '//TRIM(token)
      VERIFY_(99)
    ELSE
      IF (mapl_am_i_root()) PRINT*,'OH_data source for ONLINE fields: '//TRIM(token)
    ENDIF

    call ESMF_ConfigGetAttribute(cfg, self%compute_once_per_day,   label='compute_once_per_day:',   __RC__)

    call ESMF_ConfigGetAttribute(cfg, self%spinup_24hr_imports,    label='spinup_24hr_imports:',    __RC__)

    call ESMF_ConfigGetAttribute(cfg, self%wavelength_for_scacoef, label='wavelength_for_scacoef:', __RC__)

    call ESMF_ConfigDestroy(cfg, __RC__)


!
!   Get information from GOCART2G
!   -----------------------------
    g2g_cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (g2g_cfg, 'GOCART2G_GridComp.rc', __RC__)

    self%n_wavelengths_profile = ESMF_ConfigGetLen (g2g_cfg,    label='wavelengths_for_profile_aop_in_nm:', __RC__)

    allocate(wavelengths_profile(   self%n_wavelengths_profile), __STAT__)

    call ESMF_ConfigGetAttribute (g2g_cfg, wavelengths_profile, label='wavelengths_for_profile_aop_in_nm:', __RC__)

!   Find the desired wavelength in the vector:
    self%wavelength_index = 0
    do i = 1, self%n_wavelengths_profile
      if ( wavelengths_profile(i) == self%wavelength_for_scacoef ) self%wavelength_index = i
    end do
    _ASSERT(self%wavelength_index > 0,'Did not find OH wavelength_for_scacoef in GOCART2G wavelengths_for_profile_aop_in_nm')

    deallocate( wavelengths_profile, __STAT__ )
    call ESMF_ConfigDestroy(g2g_cfg, __RC__)


!   Is OH data driven? -- not implemented yet
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run,         Run,        __RC__)
    if (data_driven .neqv. .true.) then
       call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)
    end if
!   call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Finalize,    Finalize,   __RC__)


!   IMPORT STATE
!   -------------
    if (data_driven) then

       call MAPL_AddInternalSpec(gc,&
          short_name='OH', &
          long_name='Hydroxyl Radical', &
          units='???', &
          dims=MAPL_DimsHorzVert, &
          vlocation=MAPL_VlocationCenter, &
          restart=MAPL_RestartOptional, &
          ungridded_dims=[self%nbins], &
!         friendlyto='DYNAMICS:TURBULENCE:MOIST', &   (we don't transport OH)
          add2export=.true., __RC__)

 
       do i = 1, self%nbins 
          write (field_name, '(A, I0.3)') '', i
          call MAPL_AddImportSpec(GC,                                        &
             short_name = 'climoh'//trim(field_name),                        &
             long_name  = 'OH Mixing Ratio (bin '//trim(field_name)//')',    &
             units      = 'kg kg-1',                                         &
             restart    = MAPL_RestartSkip,                                  &
             dims       = MAPL_DimsHorzVert,                                 &
             vlocation  = MAPL_VLocationCenter, __RC__)
        end do
    end if ! (data_driven)

!   Computational Instance
!   ----------------------
    if (.not. data_driven) then


#include "OH_Internal___.h"

#include "OH_Import___.h"


!!!
!!!  Some imports are optional, depending on the OH_data_source
!!!  Use macros to allow calls to MAPL_AddImportSpec using just a single line
!!!

#define __HO__       DIMS=MAPL_DimsHorzOnly
#define __HV__       DIMS=MAPL_DimsHorzVert
#define __VC__  VLOCATION=MAPL_VLocationCenter
#define __VE__  VLOCATION=MAPL_VLocationEdge
#define __SK__    RESTART=MAPL_RestartSkip
#define __24__  REFRESH_INTERVAL=60*60*24,AVERAGING_INTERVAL=60*60*24
#define __4D__  ungridded_dims=[self%n_wavelengths_profile]

!! INSTANTANEOUS
! center
#define ADD_IMPORT_3DC(sss,lll,uuu)          call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HV__, __VC__, __RC__)

! edge
#define ADD_IMPORT_3DE(sss,lll,uuu)          call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HV__, __VE__, __RC__)

! single level, skip rst
#define ADD_IMPORT_NORST_2D(sss,lll,uuu)     call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HO__, __SK__, __RC__)

! center, skip rst
#define ADD_IMPORT_NORST_3DC(sss,lll,uuu)    call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HV__, __VC__, __SK__, __RC__)
#define ADD_IMPORT_NORST_4DC(sss,lll,uuu)    call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HV__, __VC__, __SK__, __4D__, __RC__)

! edge, skip rst
#define ADD_IMPORT_NORST_3DE(sss,lll,uuu)    call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HV__, __VE__, __SK__, __RC__)


!! 24-HOUR MEANS  (always save in import_rst)
! single level
#define ADD_IMPORT_24_2D(sss,lll,uuu)        call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HO__, __24__, __RC__)

! center
#define ADD_IMPORT_24_3DC(sss,lll,uuu)       call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HV__, __VC__, __24__, __RC__)
#define ADD_IMPORT_24_4DC(sss,lll,uuu)       call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HV__, __VC__, __24__, __4D__, __RC__)

! edge
#define ADD_IMPORT_24_3DE(sss,lll,uuu)       call MAPL_AddImportSpec(GC, SHORT_NAME=sss, LONG_NAME=lll, UNITS=uuu, __HV__, __VE__, __24__, __RC__)

!!!
!!!  Always import (to do units conversion of the predicted OH, and
!!!                 to avoid setting values in the stratosphere):
!!!
     ADD_IMPORT_NORST_2D(  'TROPP', 'Tropopause Pressure',         'Pa'      )

     ADD_IMPORT_NORST_3DC( 'T',     'air_temperature',             'K'       )
     ADD_IMPORT_NORST_3DC( 'Q',     'specific_humidity',           'kg kg-1' )

     ADD_IMPORT_NORST_3DE( 'PLE',   'edge_pressures',              'Pa'      )

!!!
!!!  Optional imports:
!!!
! Importing INSTANTANEOUS fields:

!!!  (Also use the INST values as fallback values when the avg24 values are not available)

   IMPORT_INST: IF ( self%OH_data_source == ONLINE_INST .OR. &
                    (self%OH_data_source == ONLINE_AVG24 .AND. self%spinup_24hr_imports) ) THEN

     ! These are computed before OH is run, no restart needed:
     ADD_IMPORT_NORST_3DC( 'CH4',  'Methane',                      '?'       )
     ADD_IMPORT_NORST_3DC( 'CO',   'CO',                           '?'       )
!    ADD_IMPORT_NORST_3DC( 'T',    'air_temperature',              'K'       )     [always imported]
     ADD_IMPORT_NORST_3DC( 'FCLD', 'cloud_fraction_for_radiation', '1'       )
!    ADD_IMPORT_NORST_3DC( 'Q',    'specific_humidity',            'kg kg-1' )     [always imported]

!    ADD_IMPORT_NORST_3DE( 'PLE',  'edge_pressures',               'Pa'      )     [always imported]
     ADD_IMPORT_NORST_3DE( 'ZLE',  'edge_heights',                 'm'       )

     ADD_IMPORT_NORST_4DC( 'BCSCACOEF',   'Black Carbon Scattering Coefficient [550 nm]', 'm-1' )
     ADD_IMPORT_NORST_4DC( 'OCSCACOEF', 'Organic Carbon Scattering Coefficient [550 nm]', 'm-1' )
     ADD_IMPORT_NORST_4DC( 'DUSCACOEF', '          Dust Scattering Coefficient [550 nm]', 'm-1' )
     ADD_IMPORT_NORST_4DC( 'SUSCACOEF', '           SO4 Scattering Coefficient [550 nm]', 'm-1' )
     ADD_IMPORT_NORST_4DC( 'SSSCACOEF',       'Sea Salt Scattering Coefficient [550 nm]', 'm-1' )
     ADD_IMPORT_NORST_4DC( 'NISCACOEF',        'Nitrate Scattering Coefficient [550 nm]', 'm-1' )

     ! These are computed after OH is run:
     ADD_IMPORT_3DC( 'TAUCLW', 'in_cloud_optical_thickness_for_liquid_clouds', '1' )
     ADD_IMPORT_3DC( 'TAUCLI', 'in_cloud_optical_thickness_for_ice_clouds',    '1' )

   END IF IMPORT_INST


   IMPORT_24: IF ( self%OH_data_source == ONLINE_AVG24 ) THEN

     ! Save each of these in the import_rst
     ADD_IMPORT_24_3DC(       'CH4_avg24', 'daily_mean_methane',                                        '?'       )
     ADD_IMPORT_24_3DC(        'CO_avg24', 'daily_mean_CO',                                             '?'       )
     ADD_IMPORT_24_3DC(         'T_avg24', 'daily_mean_air_temperature',                                'K'       )
     ADD_IMPORT_24_3DC(      'FCLD_avg24', 'daily_mean_cloud_fraction_for_radiation',                   '1'       )
     ADD_IMPORT_24_3DC(         'Q_avg24', 'daily_mean_specific_humidity',                              'kg kg-1' )
     ADD_IMPORT_24_3DC(    'TAUCLW_avg24', 'daily_mean_in_cloud_optical_thickness_for_liquid_clouds',   '1'       )
     ADD_IMPORT_24_3DC(    'TAUCLI_avg24', 'daily_mean_in_cloud_optical_thickness_for_ice_clouds',      '1'       )

     ADD_IMPORT_24_3DE(       'PLE_avg24', 'daily_mean_edge_pressures',                                 'Pa'      )
     ADD_IMPORT_24_3DE(       'ZLE_avg24', 'daily_mean_edge_heights',                                   'm'       )

     ADD_IMPORT_24_4DC( 'BCSCACOEF_avg24', 'daily_mean_Black_Carbon_Scattering_Coefficient_[550_nm]',   'm-1'     )
     ADD_IMPORT_24_4DC( 'OCSCACOEF_avg24', 'daily_mean_Organic_Carbon_Scattering_Coefficient_[550_nm]', 'm-1'     )
     ADD_IMPORT_24_4DC( 'DUSCACOEF_avg24', 'daily_mean_Dust_Scattering_Coefficient_[550_nm]',           'm-1'     )
     ADD_IMPORT_24_4DC( 'SUSCACOEF_avg24', 'daily_mean_SO4_Scattering_Coefficient_[550_nm]',            'm-1'     )
     ADD_IMPORT_24_4DC( 'SSSCACOEF_avg24', 'daily_mean_Sea_Salt_Scattering_Coefficient_[550_nm]',       'm-1'     )
     ADD_IMPORT_24_4DC( 'NISCACOEF_avg24', 'daily_mean_Nitrate_Scattering_Coefficient_[550_nm]',        'm-1'     )


   END IF IMPORT_24


   IMPORT_PRECOMPUTED: IF ( self%OH_data_source == PRECOMPUTED ) THEN

     ADD_IMPORT_NORST_3DC( 'oh_CH4',    'Methane for OH algorithm',                          '1'       )
     ADD_IMPORT_NORST_3DC( 'oh_CO',     'CO for OH algorithm',                               '1'       )
     ADD_IMPORT_NORST_3DC( 'oh_T',      'air_temperature',                                   'K'       )
     ADD_IMPORT_NORST_3DC( 'oh_FCLD',   'cloud_area_fraction',                               '1'       )
     ADD_IMPORT_NORST_3DC( 'oh_Q',      'specific_humidity',                                 'kg kg-1' )
     ADD_IMPORT_NORST_3DC( 'oh_TAUCLW', 'in_cloud_optical_thickness_for_liquid_clouds',      '1'       )
     ADD_IMPORT_NORST_3DC( 'oh_TAUCLI', 'in_cloud_optical_thickness_for_ice_clouds',         '1'       )

     ADD_IMPORT_NORST_3DE( 'oh_PLE',    'edge_pressures',                                    '1'       )
     ADD_IMPORT_NORST_3DE( 'oh_ZLE',    'edge_heights',                                      'm'       )

     ! Archived SCACOEF fields are 3D
     ADD_IMPORT_NORST_3DC( 'oh_BCSCACOEF',   'Black Carbon Scattering Coefficient [550 nm]', 'm-1'     )
     ADD_IMPORT_NORST_3DC( 'oh_OCSCACOEF', 'Organic Carbon Scattering Coefficient [550 nm]', 'm-1'     )
     ADD_IMPORT_NORST_3DC( 'oh_DUSCACOEF',           'Dust Scattering Coefficient [550 nm]', 'm-1'     )
     ADD_IMPORT_NORST_3DC( 'oh_SUSCACOEF',            'SO4 Scattering Coefficient [550 nm]', 'm-1'     )
     ADD_IMPORT_NORST_3DC( 'oh_SSSCACOEF',       'Sea Salt Scattering Coefficient [550 nm]', 'm-1'     )
     ADD_IMPORT_NORST_3DC( 'oh_NISCACOEF',        'Nitrate Scattering Coefficient [550 nm]', 'm-1'     )

   END IF IMPORT_PRECOMPUTED


#include "OH_Export___.h"

    end if


!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'OH_GridComp', wrap, STATUS )
    VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!============================================================================
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:
  subroutine Initialize (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes OH's Grid Component.

! !REVISION HISTORY: 
! 16oct2019  E.Sherman, A.da Silva, T.Clune, A.Darmenov - First attempt at refactoring

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),       pointer  :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
!   type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg
    type (ESMF_Config)                   :: universal_cfg
    type (wrap_)                         :: wrap
    type (OH_GridComp), pointer          :: self

    integer                              :: i, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: bin_index, prefix
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)
    real, pointer, dimension(:,:,:)      :: int_ptr
    logical                              :: data_driven

    __Iam__('Initialize')

!****************************************************************************
!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, config=universal_cfg, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get my internal private state
!   -----------------------------
    call ESMF_UserCompGetInternalState(GC, 'OH_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    call MAPL_GridGet ( grid, globalCellCountPerDim=dims, __RC__ )

!   Dust emission tuning coefficient [kg s2 m-5]. NOT bin specific.
!   ---------------------------------------------------------------
!   self%Ch_DU = Chem_UtilResVal(dims(1), dims(2), self%Ch_DU_res(:), __RC__)  !! example of using res vec

!   Get dimensions
!   ---------------
    km = dims(3)
    self%km = km

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)                        
    call MAPL_GetResource(mapl, CDT, Label='QUICKCHEM_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'OH_instance_'//trim(COMP_NAME)//'.rc', RC=STATUS)
    if (status /= 0) then
        if (mapl_am_i_root()) print*,'OH_instance_'//trim(COMP_NAME)//'.rc does not exist! &
                                      loading OH_instance_OH.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'OH_instance_OH.rc', __RC__)
    end if

    allocate( self%OH_ML(1:dims(1),1:dims(2),1:dims(3)), __STAT__ )

!   Call Generic Initialize 
!   ------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is OH data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   If this is a data component, the data is provided in the import
!   state via ExtData instead of the actual QuickChem children
!   ---------------------------------------------------------------- 
!   if ( data_driven ) then
!       providerState = import
!       prefix = 'clim'
!   else
!       providerState = export
!       prefix = ''
!   end if

!   Add attribute information for OH export.
!   call ESMF_StateGet (export, 'OH', field, __RC__)
!   call ESMF_AttributeSet(field, NAME='radius', valueList=self%radius, itemCount=self%nbins, __RC__)
!   call ESMF_AttributeSet(field, NAME='fnum', valueList=self%fnum, itemCount=self%nbins, __RC__)

!   Add attribute information to internal state variables
!   -----------------------------------------------------

!   call ESMF_StateGet (internal, 'OH', field, __RC__)
!   call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', value=self%fscav(1), __RC__)

!   if (.not. data_driven) then
!
!      Set internal OH values to 0 where above klid
!      call MAPL_GetPointer (internal, int_ptr, 'OH', __RC__)
!      call setZeroKlid4d (self%km, self%klid, int_ptr)         ! see GOCART2G Generic F90
!
!   end if

    if (data_driven) then
       instance = instanceData
    else
       instance = instanceComputational
    end if

    self%instance = instance

    call ESMF_ConfigDestroy(cfg, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


!============================================================================
!BOP
! !IROUTINE: Run 

! !INTERFACE:
  subroutine Run (GC, import, export, clock, rc)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: rc     ! Error code:

!   !DESCRIPTION: Run method for the Dust Grid Component. Determines whether to
!                 run data or computational run method.

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal

    logical                           :: data_driven

    __Iam__('Run')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is this species data driven?
!   ----------------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Update INTERNAL state variables with ExtData
!   ---------------------------------------------
    if (data_driven) then
       call Run_data (GC, import, export, internal, __RC__)
    else
       call Run1 (GC, import, export, clock, __RC__)
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

!============================================================================
!BOP
! !IROUTINE: Run1 

! !INTERFACE:
  subroutine Run1 (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION:  Computes emissions/sources for this species

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (OH_GridComp), pointer       :: self
    type(ESMF_Time)                   :: time

    integer  :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
    integer  :: import_shape_3d(3), im, jm, km

!  Input fields from fvGCM
!  -----------------------

!  All 3D fields are top-down

      REAL, POINTER, DIMENSION(:,:,:) ::  TAUCLI     => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  TAUCLW     => null()

      REAL, POINTER, DIMENSION(:,:)   ::  gmito3     => null()
      REAL, POINTER, DIMENSION(:,:)   ::  gmitto3    => null()

      ! from original GOCART
      REAL, POINTER, DIMENSION(:,:,:)   ::  BCscacoef_3D  => null()
      REAL, POINTER, DIMENSION(:,:,:)   ::  OCscacoef_3D  => null()
      REAL, POINTER, DIMENSION(:,:,:)   ::  DUscacoef_3D  => null()
      REAL, POINTER, DIMENSION(:,:,:)   ::  SUscacoef_3D  => null()
      REAL, POINTER, DIMENSION(:,:,:)   ::  SSscacoef_3D  => null()
      REAL, POINTER, DIMENSION(:,:,:)   ::  NIscacoef_3D  => null()

      ! from GOCART2G
      REAL, POINTER, DIMENSION(:,:,:,:) ::  BCscacoef_4D  => null()
      REAL, POINTER, DIMENSION(:,:,:,:) ::  OCscacoef_4D  => null()
      REAL, POINTER, DIMENSION(:,:,:,:) ::  DUscacoef_4D  => null()
      REAL, POINTER, DIMENSION(:,:,:,:) ::  SUscacoef_4D  => null()
      REAL, POINTER, DIMENSION(:,:,:,:) ::  SSscacoef_4D  => null()
      REAL, POINTER, DIMENSION(:,:,:,:) ::  NIscacoef_4D  => null()

!     REAL, POINTER, DIMENSION(:,:,:,:) ::  DU_4d       => null()

      REAL, POINTER, DIMENSION(:,:,:) ::   ZLE_BST => null()   ! m   - for XGBoost computation
      REAL, POINTER, DIMENSION(:,:,:) ::   PLE_BST => null()   ! Pa  - for XGBoost computation
      REAL, POINTER, DIMENSION(:,:,:) ::   PLE_MOD => null()   ! Pa  - current value in the model

      REAL, POINTER, DIMENSION(:,:)   :: TROPP_MOD => null()   ! Pa  - current value in the model

      REAL, POINTER, DIMENSION(:,:)   ::  ptr2d    => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  ptr3d    => null()

      REAL, POINTER, DIMENSION(:,:,:) ::  default_OH    => null()

      REAL, POINTER, DIMENSION(:,:,:) ::  Q_MOD    => null()   ! top-down  - current value in the model
      REAL, POINTER, DIMENSION(:,:,:) ::  T_MOD    => null()   ! top-down  - current value in the model

      REAL, ALLOCATABLE         ::    TV_MOD(:,:,:) ! top-down  - current value in the model
      REAL, ALLOCATABLE         :: NDWET_MOD(:,:,:) ! top-down  - current value in the model
      REAL, ALLOCATABLE         ::    PL_MOD(:,:,:) ! top-down   ! Pa  - current value in the model
      REAL, ALLOCATABLE, TARGET ::    PL_BST(:,:,:) ! top-down   ! Pa  - value for XGBoost

      REAL, ALLOCATABLE, TARGET :: tauclwDN(:,:,:) ! top-down
      REAL, ALLOCATABLE, TARGET :: taucliDN(:,:,:) ! top-down
      REAL, ALLOCATABLE, TARGET :: taucliUP(:,:,:) ! top-down
      REAL, ALLOCATABLE, TARGET :: tauclwUP(:,:,:) ! top-down


      REAL, ALLOCATABLE, TARGET ::   latarr(:,:)   ! Array of latitudes
      REAL, ALLOCATABLE, TARGET ::  stratO3(:,:)   ! Total Strat Ozone
      REAL, ALLOCATABLE, TARGET :: sza_noon(:,:)   ! SZA for local noon

      REAL, ALLOCATABLE, TARGET ::   aodUP(:,:,:) ! top-down  Aero Optical Depth above
      REAL, ALLOCATABLE, TARGET ::   aodDN(:,:,:) ! top-down  Aero Optical Depth below
      REAL, ALLOCATABLE         ::     aod(:,:,:) ! top-down  (intermediate term)

      REAL, POINTER, DIMENSION(:,:)     :: LATS
      REAL, POINTER, DIMENSION(:,:)     :: LONS

      REAL*8,  ALLOCATABLE :: gridBoxThickness(:,:,:) !                       top-down

      LOGICAL :: dynamic_k_range  ! whether the prediction routine should limit work to the troposphere

      LOGICAL :: need_to_call_BOOST

      REAL    :: tropp_min

      INTEGER :: i1, i2, j1, j2
      INTEGER :: i, j, k

      INTEGER :: JDAY
!     REAL    :: radToDeg, degToRad, pi

!?    REAL, PARAMETER ::     epsilon = (MAPL_H2OMW/MAPL_AIRMW)

!!! If we need access to INTERNAL state:
!     type(MAPL_MetaComp), pointer :: genState    ! MAPL generic state
!     type(ESMF_State)             :: INTERNAL

      TYPE(OH_BOOST_INPUT_DATA)       :: bb
      CHARACTER(LEN=ESMF_MAXPATHLEN)  :: XGBoostFilename

      TYPE(ESMF_Alarm)  :: ALARM
      LOGICAL :: RunRightNow


#include "OH_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME)//'::'//'Run1'

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, &
                        LONS = LONS, &
                        LATS = LATS, __RC__ )

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'OH_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Extract nymd(yyyymmdd) from clock
!   ---------------------------------
    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)


!   NOTES on OH_DT
!   *  If OH_DT is not set, the alarm will ring at every timestep
!   *  OH_DT will not save much time if you are using compute_once_per_day
!   *  Set OH_REFERENCE_TIME to be the HEARTBEAT, to make sure this routine runs at
!      the BEGINNING of each OH_DT interval, e.g. in AGCM.rc:   OH_REFERENCE_TIME: 000730
!      otherwise the alarm will ring during the LAST timestep of each OH_DT interval
!   *  We silence the alarm at the end of Run2
!   ----------------------------------------------------------------------------------------
    CALL MAPL_Get(mapl, RUNALARM=ALARM, __RC__)
    RunRightNow = ESMF_AlarmIsRinging(ALARM, __RC__)
!   if (mapl_am_i_root()) PRINT*,'OH HMS1, alarm:',nhms,RunRightNow
    IF ( .NOT. RunRightNow ) THEN
      RETURN_(ESMF_SUCCESS)
    END IF

    call fill_grads_template(XGBoostFilename,self%XGBoostFilePattern,nymd=nymd,nhms=nhms,__RC__)

    IF ( self%compute_once_per_day .AND. nhms > 0 ) THEN
      need_to_call_BOOST = .FALSE.
    ELSE
      need_to_call_BOOST = .TRUE.
    END IF

#include "OH_GetPointer___.h"

!   Get dimensions
!   ---------------
    import_shape_3d = shape(oh_ISOP)   ! refer to a 3D field that is always imported
    im = import_shape_3d(1)
    jm = import_shape_3d(2)
    km = import_shape_3d(3)

!   Call the ML predictor:

!!!! At the end of this subroutine, set the INTERNAL field like this:
!!!! 
!!!! !   OH does not support nbins>1, so use 3 indices instead of 4:
!!!!     OH(:,:,self%km) = OH(:,:,self%km) + &
!!!!                       oh_src * (self%CDT * MAPL_GRAV / delp(:,:,self%km))

   !  Previously:
   !  Units for OH: w_c%qa(nbeg)%data3D(i1:i2,j1:j2,1:km) is always in terms of molec/cm3

!!!!  Take out oh_src



!     pi       = 4.00*ATAN(1.00)
!     degToRad = pi/180.00
!     radToDeg = 180.00/pi

!     radToDeg = 57.2957795

      !  Initialize local variables
      !  --------------------------
      i1 = 1
      i2 = im
      j1 = 1
      j2 = jm

! Get the current model values:
      CALL MAPL_GetPointer(import,     T_MOD,     'T',  __RC__ )
      CALL MAPL_GetPointer(import,     Q_MOD,     'Q',  __RC__ )
      CALL MAPL_GetPointer(import,   PLE_MOD,   'PLE',  __RC__ )
      CALL MAPL_GetPointer(import, TROPP_MOD, 'TROPP',  __RC__ )

!  !  Compute some MOD fields -- current values from MODEL:

      allocate(    PL_MOD(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(    TV_MOD(i1:i2,j1:j2,1:km), __STAT__ )
      allocate( NDWET_MOD(i1:i2,j1:j2,1:km), __STAT__ )

!  !  Layer mean pressures  (Pa)
!  !  --------------------
      _ASSERT(lbound(PLE_MOD,3)==0, "Error. Expecting PLE starting index 0")
      PL_MOD = (PLE_MOD(:,:,0:km-1)+PLE_MOD(:,:,1:km))*0.5

   !  Virtual Temperature (K)
      TV_MOD = T_MOD * (1.0 + Q_MOD/MAPL_EPSILON)/(1.0 + Q_MOD)

   !  Moist air number density  (molec/m3)
   !  NOTE: Odd units for MAPL_AVOGAD = molec / kmol   (usually expressed in terms of mol)
   !  NOTE: Odd units for MAPL_RUNIV  = J / (kmol K)   (usually expressed in terms of mol)
   !  Use online, instantaneous values to compute
   !  ------------------------
      NDWET_MOD = (MAPL_AVOGAD * PL_MOD) / (MAPL_RUNIV * TV_MOD)

  PREP_FOR_BOOST: IF ( need_to_call_BOOST ) THEN
      !  Initialize the BOOST pointers
      !  -----------------------------
      NULLIFY( bb%LAT        )   ! online
      NULLIFY( bb%PL         )   ! online (3 options)
      NULLIFY( bb%T          )   ! online (3 options)
      NULLIFY( bb%NO2        )   ! ExtData
      NULLIFY( bb%O3         )   ! ExtData
      NULLIFY( bb%CH4        )   ! online (3 options)
      NULLIFY( bb%CO         )   ! online (3 options)
      NULLIFY( bb%ISOP       )   ! ExtData
      NULLIFY( bb%ACET       )   ! ExtData
      NULLIFY( bb%C2H6       )   ! ExtData
      NULLIFY( bb%C3H8       )   ! ExtData
      NULLIFY( bb%PRPE       )   ! ExtData
      NULLIFY( bb%ALK4       )   ! ExtData
      NULLIFY( bb%MP         )   ! ExtData
      NULLIFY( bb%H2O2       )   ! ExtData
      NULLIFY( bb%TAUCLWDN   )   ! derived from online (3 options)
      NULLIFY( bb%TAUCLIDN   )   ! derived from online (3 options)
      NULLIFY( bb%TAUCLIUP   )   ! derived from online (3 options)
      NULLIFY( bb%TAUCLWUP   )   ! derived from online (3 options)
      NULLIFY( bb%CLOUD      )   ! online (3 options)
      NULLIFY( bb%QV         )   ! online (3 options)
      NULLIFY( bb%GMISTRATO3 )   ! computed from ExtData
      NULLIFY( bb%ALBUV      )   ! ExtData climatology
      NULLIFY( bb%AODUP      )   ! derived from 6 aerosols - (3 options)
      NULLIFY( bb%AODDN      )   ! derived from 6 aerosols - (3 options)
      NULLIFY( bb%CH2O       )   ! ExtData
      NULLIFY( bb%SZA        )   ! computed for local noon

      ! Reserve Some local work space
      !------------------------------
      allocate(          latarr(i1:i2,j1:j2),      __STAT__ )
      allocate(         stratO3(i1:i2,j1:j2),      __STAT__ )
      allocate(        sza_noon(i1:i2,j1:j2),      __STAT__ )

      allocate(gridBoxThickness(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(          PL_BST(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(        tauclwDN(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(        taucliDN(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(        taucliUP(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(        tauclwUP(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(             aod(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(           aodUP(i1:i2,j1:j2,1:km), __STAT__ )
      allocate(           aodDN(i1:i2,j1:j2,1:km), __STAT__ )

      !Get imports

      ! Initialize flag indicating that we are in the daily mean spinup period:
      self%use_inst_values = .FALSE.
      if ( self%OH_data_source == ONLINE_AVG24 ) then
        CALL MAPL_GetPointer(import, ptr3d, 'T_avg24', __RC__ )
        ! If we are within the first 24 hours of integration,
        ! and there was no field saved in the import_rst, then the Coupler returns
        ! a field of all zero, so we should use the Instantaneous values instead.
        ! We assume that if T is valid, all other daily mean imports are valid too.
        if ( ptr3d(1,1,1) == 0.0 ) self%use_inst_values = .TRUE.
      endif

      if (mapl_am_i_root()) then
       if (       self%use_inst_values ) PRINT*,'OH is in the SPINUP period for 24-hour averages'
       if ( .NOT. self%use_inst_values ) PRINT*,'OH is *NOT* in the SPINUP period for 24-hour averages'
      endif

! Get the value to use when calling XGBoost:
! (save directly into bb)
      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, bb%T,    'oh_T',           __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, bb%T,       'T',           __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
         if (       self%use_inst_values )       CALL MAPL_GetPointer(import, bb%T,       'T',           __RC__ )
         if ( .NOT. self%use_inst_values )       CALL MAPL_GetPointer(import, bb%T,       'T_avg24',     __RC__ )
      endif

! CALL MAPL_GetPointer(export, ptr3d, 'DIAG_T_in_OH',    __RC__)
! IF (ASSOCIATED(ptr3d))   ptr3d(:,:,:) = bb%T

! Get the value to use when calling XGBoost:
! (save directly into bb)
      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, bb%QV,   'oh_Q',           __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, bb%QV,      'Q',           __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
         if (       self%use_inst_values )       CALL MAPL_GetPointer(import, bb%QV,      'Q',           __RC__ )
         if ( .NOT. self%use_inst_values )       CALL MAPL_GetPointer(import, bb%QV,      'Q_avg24',     __RC__ )
      endif

! Get the value to use when calling XGBoost:
      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, PLE_BST, 'oh_PLE',         __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, PLE_BST,    'PLE',         __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
         if (       self%use_inst_values )       CALL MAPL_GetPointer(import, PLE_BST,    'PLE',         __RC__ )
         if ( .NOT. self%use_inst_values )       CALL MAPL_GetPointer(import, PLE_BST,    'PLE_avg24',   __RC__ )
      endif

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, ZLE_BST, 'oh_ZLE',         __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, ZLE_BST,    'ZLE',         __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, ZLE_BST,    'ZLE',         __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, ZLE_BST,    'ZLE_avg24',   __RC__ )
      endif

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, TAUCLW, 'oh_TAUCLW',       __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, TAUCLW,    'TAUCLW',       __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, TAUCLW,    'TAUCLW',       __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, TAUCLW,    'TAUCLW_avg24', __RC__ )
      endif

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, TAUCLI, 'oh_TAUCLI',       __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, TAUCLI,    'TAUCLI',       __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, TAUCLI,    'TAUCLI',       __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, TAUCLI,    'TAUCLI_avg24', __RC__ )
      endif



!   CALL MAPL_GetPointer(import, DU_4d, 'DUSCACOEF', __RC__)
!   du_shape = shape(DU_4d)
!   IF (MAPL_AM_I_ROOT()) print *, 'Shape of DU_4d:', du_shape(1), du_shape(2), du_shape(3), du_shape(4)


!!! If we need access to INTERNAL state:
!     CALL MAPL_GetObjectFromGC( this%gcESMF, genState, __RC__ )
!     CALL MAPL_Get( genState, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )


!     Note that the archived Scattering Coefficient fields are 3D, while the online fields are 4D:

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, BCscacoef_3D,  'oh_BCSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, BCscacoef_4D,     'BCSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, BCscacoef_4D,     'BCSCACOEF',       __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, BCscacoef_4D,     'BCSCACOEF_avg24', __RC__ )
      endif

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, OCscacoef_3D,  'oh_OCSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, OCscacoef_4D,     'OCSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, OCscacoef_4D,     'OCSCACOEF',       __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, OCscacoef_4D,     'OCSCACOEF_avg24', __RC__ )
      endif

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, DUscacoef_3D,  'oh_DUSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, DUscacoef_4D,     'DUSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, DUscacoef_4D,     'DUSCACOEF',       __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, DUscacoef_4D,     'DUSCACOEF_avg24', __RC__ )
      endif

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, SUscacoef_3D,  'oh_SUSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, SUscacoef_4D,     'SUSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, SUscacoef_4D,     'SUSCACOEF',       __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, SUscacoef_4D,     'SUSCACOEF_avg24', __RC__ )
      endif

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, SSscacoef_3D,  'oh_SSSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, SSscacoef_4D,     'SSSCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, SSscacoef_4D,     'SSSCACOEF',       __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, SSscacoef_4D,     'SSSCACOEF_avg24', __RC__ )
      endif

      if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, NIscacoef_3D,  'oh_NISCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, NIscacoef_4D,     'NISCACOEF',       __RC__ )
      if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )     CALL MAPL_GetPointer(import, NIscacoef_4D,     'NISCACOEF',       __RC__ )
           if ( .NOT. self%use_inst_values )     CALL MAPL_GetPointer(import, NIscacoef_4D,     'NISCACOEF_avg24', __RC__ )
      endif

!!! DOES the auto code generator provide for aliasing (pointer named differently from import)?


      gmito3  => oh_GMITO3
      gmitto3 => oh_GMITTO3

      latarr(:,:) = LATS(:,:)*MAPL_RADIANS_TO_DEGREES

      stratO3(:,:) = gmito3(:,:) - gmitto3(:,:)

   !  Cell depth (in meters)
   !  ---------------------
      _ASSERT(lbound(ZLE_BST,3)==0, "Error. Expecting ZLE starting index 0")
      gridBoxThickness(:,:,1:km) = ZLE_BST(:,:,0:km-1)-ZLE_BST(:,:,1:km)

   !  Aerosol Optical Depth
   !  ---------------------
      IF ( self%OH_data_source == PRECOMPUTED  ) THEN
        aod = gridBoxThickness * ( BCscacoef_3D + OCscacoef_3D + DUscacoef_3D + &
                                   SUscacoef_3D + SSscacoef_3D + NIscacoef_3D )
      ELSE
        aod = gridBoxThickness * ( BCscacoef_4D(:,:,:,self%wavelength_index) + &
                                   OCscacoef_4D(:,:,:,self%wavelength_index) + &
                                   DUscacoef_4D(:,:,:,self%wavelength_index) + &
                                   SUscacoef_4D(:,:,:,self%wavelength_index) + &
                                   SSscacoef_4D(:,:,:,self%wavelength_index) + &
                                   NIscacoef_4D(:,:,:,self%wavelength_index) )
      END IF

      DO k=1,km

        tauclwDN(:,:,k)  = SUM( TAUCLW(:,:,k:km), 3 )
        taucliDN(:,:,k)  = SUM( TAUCLI(:,:,k:km), 3 )
        taucliUP(:,:,k)  = SUM( TAUCLI(:,:,1:k ), 3 )
        tauclwUP(:,:,k)  = SUM( TAUCLW(:,:,1:k ), 3 )

           aodUP(:,:,k)  = SUM(    aod(:,:,1:k ), 3 )
           aodDN(:,:,k)  = SUM(    aod(:,:,k:km), 3 )

      END DO

      ! Compute SZA
      JDAY     = JulianDay(nymd)
      sza_noon = computeSolarZenithAngle_LocalNoon (JDAY, LATS, LONS, i1, i2, j1, j2)

!!! Assign values to all the input parameters held in bb:

        bb%LAT => latarr

        PL_BST = (PLE_BST(:,:,0:km-1)+PLE_BST(:,:,1:km))*0.5
        bb%PL => PL_BST

!       bb%T  set above

        CALL MAPL_GetPointer(import, bb%NO2,            'oh_NO2',  __RC__ )
        CALL MAPL_GetPointer(import, bb%O3,             'oh_O3',   __RC__ )

        if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, bb%CH4,  'oh_CH4',       __RC__ )
        if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, bb%CH4,     'CH4',       __RC__ )
        if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )       CALL MAPL_GetPointer(import, bb%CH4,     'CH4',       __RC__ )
           if ( .NOT. self%use_inst_values )       CALL MAPL_GetPointer(import, bb%CH4,     'CH4_avg24', __RC__ )
        endif

        if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, bb%CO,   'oh_CO',        __RC__ )
        if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, bb%CO,      'CO',        __RC__ )
        if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )       CALL MAPL_GetPointer(import, bb%CO,      'CO',        __RC__ )
           if ( .NOT. self%use_inst_values )       CALL MAPL_GetPointer(import, bb%CO,      'CO_avg24',  __RC__ )
        endif

        CALL MAPL_GetPointer(import, bb%ISOP,           'oh_ISOP', __RC__ )
        CALL MAPL_GetPointer(import, bb%ACET,           'oh_ACET', __RC__ )
        CALL MAPL_GetPointer(import, bb%C2H6,           'oh_C2H6', __RC__ )
        CALL MAPL_GetPointer(import, bb%C3H8,           'oh_C3H8', __RC__ )
        CALL MAPL_GetPointer(import, bb%PRPE,           'oh_PRPE', __RC__ )
        CALL MAPL_GetPointer(import, bb%ALK4,           'oh_ALK4', __RC__ )
        CALL MAPL_GetPointer(import, bb%MP,             'oh_MP',   __RC__ )
        CALL MAPL_GetPointer(import, bb%H2O2,           'oh_H2O2', __RC__ )

        bb%TAUCLWDN => tauclwDN
        bb%TAUCLIDN => taucliDN
        bb%TAUCLIUP => taucliUP
        bb%TAUCLWUP => tauclwUP

        if ( self%OH_data_source == PRECOMPUTED  ) CALL MAPL_GetPointer(import, bb%CLOUD,  'oh_FCLD',       __RC__ )
        if ( self%OH_data_source == ONLINE_INST  ) CALL MAPL_GetPointer(import, bb%CLOUD,     'FCLD',       __RC__ )
        if ( self%OH_data_source == ONLINE_AVG24 ) then
           if (       self%use_inst_values )       CALL MAPL_GetPointer(import, bb%CLOUD,     'FCLD',       __RC__ )
           if ( .NOT. self%use_inst_values )       CALL MAPL_GetPointer(import, bb%CLOUD,     'FCLD_avg24', __RC__ )
        endif

!       bb%QV  set above

        bb%GMISTRATO3 => stratO3   ! 2D  derived from M2G

        CALL MAPL_GetPointer(import, bb%ALBUV,         'oh_ALBUV', __RC__ )   ! single layer from climo (ExtData)

        bb%AODUP  =>  aodUP
        bb%AODDN  =>  aodDN

        CALL MAPL_GetPointer(import, bb%CH2O,           'oh_CH2O', __RC__ )

        bb%SZA    =>  sza_noon

  END IF PREP_FOR_BOOST


        ! default OH values; prediction will only be done in the troposphere
        CALL MAPL_GetPointer(import, default_OH,        'oh_OH',   __RC__ )

        CALL MAPL_MaxMin ( 'OH: OH From M2G ', default_OH )

!       ! capture MERRA2GMI version as diagnostic:
        CALL MAPL_GetPointer(export, ptr3d,     'DIAG_OH_M2G',     __RC__)
        IF (ASSOCIATED(ptr3d))   ptr3d(:,:,:) = default_OH(:,:,:)


  CALL_BOOST: IF ( need_to_call_BOOST ) THEN

        self%OH_ML(:,:,:) = 0.0

        dynamic_k_range = .NOT.  self%compute_once_per_day

        tropp_min = 40.0 * 100  ! convert hPa to Pa

        CALL predict_OH_with_XGB( XGBoostFilename, im, jm, km, dynamic_k_range, tropp_min, &
                                  PL_MOD, TROPP_MOD, bb, self%OH_ML, __RC__ )

        CALL MAPL_GetPointer(export, ptr3d,     'OH_boost',     __RC__)
        IF (ASSOCIATED(ptr3d))   ptr3d(:,:,:) = self%OH_ML(:,:,:)

  END IF CALL_BOOST


!       Only use Boost values in the troposphere
!       Note: we test using online values that change every timestep
        DO j= 1,jm
        DO i= 1,im
          WHERE ( PL_MOD(i,j,:) .GT. TROPP_MOD(i,j) )       ! Note the same comparison at @@@
            OH(i,j,:) = self%OH_ML(i,j,:)
          ELSEWHERE
            OH(i,j,:) = default_OH(i,j,:)
          END WHERE
        END DO
        END DO

!!
!!  March 2022:  Store the XGBoost version of OH in data3d; convert from mol/mol to molec/cm3
!!               This makes it available to CO and CH4
!!
      !   [mol_OH/mol_MoistAir] * nd_moist_air[molec_MoistAir/m3] = [molec_OH/m3]
      !   * 1.0e-6 = molec_OH/cm3
      OH = (OH * NDWET_MOD) * 1.0e-6

! debug - DIAG for moist air number density
      CALL MAPL_GetPointer(export, ptr3d,   'DIAG_NDWET',          __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = NDWET_MOD(:,:,:)


! Handle diagnostic and deallocation:
  AFTER_BOOST: IF ( need_to_call_BOOST ) THEN

      ! Diagnostics that are only meaningful when we call BOOST:

      CALL MAPL_GetPointer(export, ptr2d,      'DIAG_LAT',         __RC__)
      IF (ASSOCIATED(ptr2d))       ptr2d(:,:)   = bb%LAT

! IF ( self%OH_data_source == ONLINE_AVG24 ) THEN
!     CALL MAPL_GetPointer(export, ptr3d,   'DIAG_T_avg24',    __RC__)
!     IF (ASSOCIATED(ptr3d))  CALL MAPL_GetPointer(import, ptr3d, 'T_avg24', __RC__ )
! ENDIF

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_TAUCLWDN',    __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%TAUCLWDN

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_TAUCLIDN',    __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%TAUCLIDN

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_TAUCLIUP',    __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%TAUCLIUP

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_TAUCLWUP',    __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%TAUCLWUP

      CALL MAPL_GetPointer(export, ptr2d,      'DIAG_GMISTRATO3',  __RC__)
      IF (ASSOCIATED(ptr2d))       ptr2d(:,:)   = bb%GMISTRATO3

      CALL MAPL_GetPointer(export, ptr2d,      'DIAG_ALBUV',       __RC__)
      IF (ASSOCIATED(ptr2d))       ptr2d(:,:)   = bb%ALBUV

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_AODUP',       __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%AODUP

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_AODDN',       __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%AODDN

      CALL MAPL_GetPointer(export, ptr2d,      'DIAG_SZA',         __RC__)
      IF (ASSOCIATED(ptr2d))       ptr2d(:,:)   = bb%SZA

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_PL',          __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%PL

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_T',           __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%T

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_CH4',         __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%CH4

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_CO',          __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%CO

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_CLOUD',       __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%CLOUD

      CALL MAPL_GetPointer(export, ptr3d,      'DIAG_QV',          __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = bb%QV

      CALL MAPL_GetPointer(export, ptr3d,   'DIAG_ZLE',            __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = ZLE_BST

      CALL MAPL_GetPointer(export, ptr3d,   'DIAG_AOD',            __RC__)
      IF (ASSOCIATED(ptr3d))       ptr3d(:,:,:) = AOD


      CALL MAPL_GetPointer(export, ptr3d, 'DIAG_SC_BC',       __RC__)
      IF (ASSOCIATED(ptr3d)) THEN
        IF ( self%OH_data_source == PRECOMPUTED ) THEN
          ptr3d(:,:,:) =     BCscacoef_3D
        ELSE
          ptr3d(:,:,:) =     BCscacoef_4D(:,:,:,self%wavelength_index)
        ENDIF
      ENDIF
      CALL MAPL_GetPointer(export, ptr3d, 'DIAG_SC_OC',       __RC__)
      IF (ASSOCIATED(ptr3d)) THEN
        IF ( self%OH_data_source == PRECOMPUTED ) THEN
          ptr3d(:,:,:) =     OCscacoef_3D
        ELSE
          ptr3d(:,:,:) =     OCscacoef_4D(:,:,:,self%wavelength_index)
        ENDIF
      ENDIF
      CALL MAPL_GetPointer(export, ptr3d, 'DIAG_SC_DU',       __RC__)
      IF (ASSOCIATED(ptr3d)) THEN
        IF ( self%OH_data_source == PRECOMPUTED ) THEN
          ptr3d(:,:,:) =     DUscacoef_3D
        ELSE
          ptr3d(:,:,:) =     DUscacoef_4D(:,:,:,self%wavelength_index)
        ENDIF
      ENDIF
      CALL MAPL_GetPointer(export, ptr3d, 'DIAG_SC_SU',       __RC__)
      IF (ASSOCIATED(ptr3d)) THEN
        IF ( self%OH_data_source == PRECOMPUTED ) THEN
          ptr3d(:,:,:) =     SUscacoef_3D
        ELSE
          ptr3d(:,:,:) =     SUscacoef_4D(:,:,:,self%wavelength_index)
        ENDIF
      ENDIF
      CALL MAPL_GetPointer(export, ptr3d, 'DIAG_SC_SS',       __RC__)
      IF (ASSOCIATED(ptr3d)) THEN
        IF ( self%OH_data_source == PRECOMPUTED ) THEN
          ptr3d(:,:,:) =     SSscacoef_3D
        ELSE
          ptr3d(:,:,:) =     SSscacoef_4D(:,:,:,self%wavelength_index)
        ENDIF
      ENDIF
      CALL MAPL_GetPointer(export, ptr3d, 'DIAG_SC_NI',       __RC__)
      IF (ASSOCIATED(ptr3d)) THEN
        IF ( self%OH_data_source == PRECOMPUTED ) THEN
          ptr3d(:,:,:) =     NIscacoef_3D
        ELSE
          ptr3d(:,:,:) =     NIscacoef_4D(:,:,:,self%wavelength_index)
        ENDIF
      ENDIF


      DEALLOCATE(gridBoxThickness, PL_BST, latarr, stratO3, sza_noon,   &
                 tauclwDN, taucliDN, taucliUP, tauclwUP, &
                 aod, aodUP, aodDN, __STAT__ )

  END IF AFTER_BOOST

      DEALLOCATE( PL_MOD, TV_MOD, NDWET_MOD, __STAT__ )

    RETURN_(ESMF_SUCCESS)

  end subroutine Run1

!============================================================================
!BOP
! !IROUTINE: Run2 

! !INTERFACE:

  subroutine Run2 (GC, import, export, clock, RC)

    ! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run2 method for the Dust Grid Component.

!EOP
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal
    type (wrap_)                      :: wrap
    type (OH_GridComp), pointer     :: self

! #include "OH_DeclarePointer___.h"

    type(ESMF_Time)   :: time
    integer  :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
    TYPE(ESMF_Alarm)  :: ALARM
    LOGICAL :: RunRightNow

    __Iam__('Run2')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)


    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)

    CALL MAPL_Get(mapl, RUNALARM=ALARM, __RC__)
    RunRightNow = ESMF_AlarmIsRinging(ALARM, __RC__)
!   if (mapl_am_i_root()) PRINT*,'OH HMS2, alarm:',nhms,RunRightNow
    IF ( .NOT. RunRightNow ) THEN
      RETURN_(ESMF_SUCCESS)
    END IF


!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'OH_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!  #include "OH_GetPointer___.h"

!!  In G2G, here is where they call Chem_Settling, and (for each bin) DryDeposition

!!  In G2G, here is where they call WetRemovalGOCART2G for each bin, and Aero_Compute_Diags

    CALL ESMF_AlarmRingerOff(ALARM, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

!============================================================================
!BOP
! !IROUTINE: Run_data -- ExtData Dust Grid Component

! !INTERFACE:
  subroutine Run_data (GC, import, export, internal, RC)

    ! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC       ! Gridded component 
    type (ESMF_State),    intent(inout) :: import   ! Import state
    type (ESMF_State),    intent(inout) :: export   ! Export state
    type (ESMF_State),    intent(inout) :: internal ! Interal state
    integer, optional,    intent(  out) :: RC       ! Error code:

! !DESCRIPTION: Updates pointers in Internal state with fields from ExtData. 

! Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (wrap_)                       :: wrap
    type (OH_GridComp), pointer        :: self

    integer                            :: i
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:)    :: ptr3d_intern
    real, pointer, dimension(:,:,:)    :: ptr3d_import

    __Iam__('Run_data')

!EOP
!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'OH_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Update internal data pointers with ExtData
!   ------------------------------------------
    call MAPL_GetPointer (internal, NAME='OH', ptr=ptr3d_intern, __RC__)

    IF ( self%nbins /= 1 ) THEN
      PRINT*,'expecting only 1 OH bin'
      RETURN_(123)
    ENDIF

    do i = 1, self%nbins
!       write (field_name, '(A, I0.3)') 'oh', i
!       call MAPL_GetPointer (import,  NAME='clim'//trim(field_name), ptr=ptr3d_import, __RC__)

        CALL MAPL_GetPointer(import, ptr3d_import, 'oh_OH',  __RC__ )   ! top-down   FROM M2G (ExtData)

!       ptr4d_intern(:,:,:,i) = ptr3d_import
        ptr3d_intern(:,:,:  ) = ptr3d_import
    end do

    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data

!!
!!
!!  Finalize would require that QuickChem (parent GC) also has Finalize
!!
!!      If implemented, it should free  self%OH_ML
!!
!!


!EOC
!--------------------------------------------------------------------------
!BOP
      INTEGER FUNCTION JulianDay(nymd)
         INTEGER :: nymd
!
! !DESCRIPTION:
! Determines the Julian day: number between 1 and 365 (or 366).
         INTEGER :: ny, mm, dd
         INTEGER :: m, ds
         INTEGER :: days(12)
!EOP
!--------------------------------------------------------------------------
!BOC
         data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

         ny = nymd / 10000
         mm = mod(nymd, 10000) / 100
         dd = mod(nymd,   100)

         ds = dd

         if ( mm .ne. 1) then
            do m=1, mm-1
               if ( m.eq.2  .and. leap_year(ny) ) then
                  ds = ds + 29
               else
                  ds = ds + days(m)
               endif
            enddo
         endif

         JulianDay = ds
         
      END FUNCTION JulianDay
!EOC
!--------------------------------------------------------------------------
!BOP
      function leap_year(ny)
!
! Determine if year ny is a leap year
!
! Author: S.-J. Lin
      implicit none
      logical leap_year
      integer ny
      integer ny00
!EOP
!--------------------------------------------------------------------------
!BOC

!
! No leap years prior to 0000
!
      parameter ( ny00 = 0000 )   ! The threshold for starting leap-year

      if( ny >= ny00 ) then
         if( mod(ny,100) == 0. .and. mod(ny,400) == 0. ) then
             leap_year = .true.
         elseif( mod(ny,4) == 0. .and. mod(ny,100) /= 0.  ) then
             leap_year = .true.
         else
             leap_year = .false.
         endif
      else
          leap_year = .false.
      endif

      return
    end function leap_year

!EOC
!--------------------------------------------------------------------------

!-------------------------------------------------------------------------------------

! INTERESTING THINGS
!
!   _ASSERT(any(self%emission_scheme == [character(len=7) :: 'ginoux','other']), "Error. Unallowed emission scheme: "//trim(self%emission_scheme)//". Allowed: ginoux, other")

!     associate (scheme => self%emission_scheme)
! #include "OH_Import___.h"
!     end associate

!     associate (scheme => self%emission_scheme)
! #include "OH_GetPointer___.h"
    ! end associate

end module OH_GridCompMod
