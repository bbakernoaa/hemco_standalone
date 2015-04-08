!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: grid_mod
!
! !DESCRIPTION: Module GRID\_MOD contains variables and routines which are 
!  used to specify the parameters of a GEOS-Chem horizontal grid. Grid 
!  parameters are computed as 3D arrays, which are required for interfacing
!  with a GCM.
!\\  
!\\
! !INTERFACE: 
!
MODULE Grid_Mod
! 
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: DoGridComputation

  REAL(f8), PARAMETER :: PI     = 3.14159265358979323e+0_f8
  REAL(f8), PARAMETER :: PI_180 = PI / 180.0_f8
  REAL(f8), PARAMETER :: Re     = 6.375e+6_f8
!
! !REVISION HISTORY:
!  23 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  01 Mar 2012 - R. Yantosca - Validated for nested grids
!  03 Apr 2012 - M. Payer    - Added ySin for map_a2a regrid (M. Cooper)
!  04 Dec 2012 - R. Yantosca - Modified for GIGC running in ESMF environment
!  26 Feb 2013 - R. Yantosca - Fixed bug in computation of lons & lats when
!                              connecting GEOS-Chem to the GEOS-5 GCM
!  19 May 2013 - C. Keller   - Added wrapper routine DoGridComputation so that
!                              module can also be used by HEMCO.
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  26 Mar 2015 - R. Yantosca - Removed obsolete, commented-out code
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DoGridComputation 
!
! !DESCRIPTION: Subroutine DoGridComputation initializes the longitude, 
!  latitude and surface area arrays. This used to be subroutine COMPUTE\_GRID.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DoGridComputation( am_I_Root,                               &
                                I1,   I2,    J1,    J2,    JSP,   JNP,   &
                                L1,   L2,    DLON,  DLAT,  I_LO,  J_LO,  &
                                IOFF, JOFF,  XMD,   XDG,   YMD,   YDG,   & 
                                YSN,  YMDR,  YDGR,  YMDRW, YDGRW, AM2, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root                      ! Root CPU?

    ! Variables with local CPU indices
    INTEGER,  INTENT(IN)  :: I1,  I2                       ! Min lon index
    INTEGER,  INTENT(IN)  :: J1,  J2                       ! Local lat indices
    INTEGER,  INTENT(IN)  :: JSP, JNP                      ! Polar lat indices
    INTEGER,  INTENT(IN)  :: L1,  L2                       ! Local lev indices
    REAL(fp), INTENT(IN)  :: DLON(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lon [deg]
    REAL(fp), INTENT(IN)  :: DLAT(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lat [deg]

    ! Variables with global CPU indices
    INTEGER,  INTENT(IN)  :: I_LO                          ! Min global lon
    INTEGER,  INTENT(IN)  :: J_LO                          ! Min global lat

    ! Offsets (for nested grids) 
    INTEGER,  INTENT(IN)  :: IOFF
    INTEGER,  INTENT(IN)  :: JOFF
!
! !OUTPUT PARAMETERS:
! 
    REAL(fp), INTENT(OUT) :: XMD  (:,:,:)                  ! Lon centers [deg]
    REAL(fp), INTENT(OUT) :: XDG  (:,:,:)                  ! Lon edges [deg]
    REAL(fp), INTENT(OUT) :: YMD  (:,:,:)                  ! Lat centers [deg]
    REAL(fp), INTENT(OUT) :: YDG  (:,:,:)                  ! Lat edges [deg]
    REAL(fp), INTENT(OUT) :: YSN  (:,:,:)                  ! SIN( lat edges )
    REAL(fp), INTENT(OUT) :: YMDR (:,:,:)                  ! Lat centers [rad]
    REAL(fp), INTENT(OUT) :: YDGR (:,:,:)                  ! Lat edges [rad]
    REAL(fp), INTENT(OUT) :: YMDRW(:,:,:)                  ! window lat centers
    REAL(fp), INTENT(OUT) :: YDGRW(:,:,:)                  ! Window lat edes
    REAL(fp), INTENT(OUT) :: AM2  (:,:,:)                  ! Area [m2]
!
! !OUTPUT PARAMETERS:
!  
    INTEGER, INTENT(OUT) :: RC                             ! Success or failure?
!
! !REMARKS:
!  (1) Lon/lat loop indices IG, JG are global indices.
!  (2) Lon/lat loop indices I,  J  are local to each CPU.
!  (3) We do not need to have global loop indices for vertical levels,
!       because we will always decompose the grid for MPI parallelization
!       in longitude and/or latitude.  All vertical levels must be present
!       on each CPU for the grid-independent GEOS-Chem to function properly.
!
! !REVISION HISTORY:
!  23 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  03 Dec 2012 - R. Yantosca - Add RC to argument list
!  04 Dec 2012 - R. Yantosca - Now define arrays with local CPU lon/lat bounds
!  07 Dec 2012 - R. Yantosca - Bug fix: make sure the last longitude edge is
!                              computed properly.  Test for IG==I2, not I==I2.
!  07 Dec 2012 - R. Yantosca - Also do not apply half-polar boxes when running
!                              in ESMF environment
!  26 Feb 2013 - R. Yantosca - Bug fix: now compute IND_X and IND_Y properly
!                              when connecting GEOS-Chem to the GEOS-5 GCM
!  21 Mar 2013 - R. Yantosca - Add fix to prevent zero surface area at poles
!  21 Mar 2013 - R. Yantosca - Rename loop indices to prevent confusion
!  06 Jun 2013 - M. Payer    - Add fix to compute sine of last latitude edge
!                              for MAP_A2A regridding (C. Keller)
!  19 May 2014 - C. Keller   - Renamed from Compute_grid to DoGridComputation.
!  06 Nov 2014 - C. Keller   - Now use LBOUND to get leftmost index of YMDRW.
!  26 Mar 2015 - R. Yantosca - Fix apparent optimization error by using 
!                              scalars in call to the SIN function
!  26 Mar 2015 - R. Yantosca - Cosmetic changes; improve indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER  :: I,     J,       L,        IG,        JG,      LBND
    REAL(fp) :: SIN_N, SIN_S,   SIN_DIFF, YEDGE_VAL, YSIN_VAL

    ! Arrays
    REAL(fp) :: IND_X( I1:I2+1 )
    REAL(fp) :: IND_Y( J1:J2+1 )

    !======================================================================
    ! Initialization
    !======================================================================

    ! Index array for longitudes
    DO I = I1, I2+1
       IND_X(I) = ( ( I + IOFF - 1 ) * 1e+0_fp ) + ( I_LO - 1 )
    ENDDO

    ! Index array for latitudes
    DO J = J1, J2+1
       IND_Y(J) = ( ( J + JOFF - 1 ) * 1e+0_fp ) + ( J_LO - 1 )
    ENDDO

    ! Left bound of YMDRW. Since we now pass YMID_R_W as an argument to 
    ! this routine, we cannot know for sure that the 2nd subscript starts
    ! at index 0 (ckeller, 11/06/14).
    LBND = LBOUND(YMDRW,2)

    !======================================================================
    ! Compute longitude and latitude arrays
    !======================================================================
    
    ! We can set L=1 instead of looping over vertical levels
    L = 1
       
    !----------------------------------------------------------------------
    ! Longitude center and edge arrays
    !----------------------------------------------------------------------

    ! Loop over local latitudes
    DO J = J1, J2

       ! Loop over local longitudes
       DO I = I1, I2

          ! Longitude centers
          XMD(I,J,L)  = ( DLON(I,J,L) * IND_X(I) ) - 180e+0_fp
          
          ! Longitude edges
          XDG(I,J,L) = XMD(I,J,L) - ( DLON(I,J,L) * 0.5e+0_fp )

          ! Compute the last longitude edge
          IF ( I == I2 ) THEN
             XDG(I+1,J,L) = XDG(I,J,L) + DLON(I,J,L)
          ENDIF
             
       ENDDO
    ENDDO

    !----------------------------------------------------------------------
    ! Latitude center and edge arrays
    !-------------------------------------------------------------------===

    ! Loop over local latitudes
    DO J = J1, J2

       ! Global latitude index
       JG = J + ( J_LO - 1 )

       ! Loop over local longitudes
       DO I = I1, I2

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%   LATITUDE CENTERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !----------------------------------------------------------------
          !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
          !
          ! Do not define half-sized polar boxes (bmy, 3/21/13)
          !----------------------------------------------------------------
#else

# if defined( GCAP )

          !----------------------------------------------------------------
          !             %%%%% TRADITIONAL GEOS-Chem %%%%%
          !
          ! For the GCAP model, there are no half-size polar boxes.
          ! Compute the latitude centers accordingly.  (bmy, 7/2/13)
          !----------------------------------------------------------------

          ! Lat centers (degrees)
          YMD(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 88e+0_fp

#else

          !----------------------------------------------------------------
          !             %%%%% TRADITIONAL GEOS-Chem %%%%%
          !
          ! Current practice in the standard GEOS-Chem is to make
          ! the polar grid boxes (for non-GCAP global grids only) be
          ! half the size of other grid boxes.  This lets us make +90
          ! degrees and -90 degrees be the edges of the grid.
          ! (bmy, 7/2/13)
          !----------------------------------------------------------------

          ! Lat centers (degrees)
          YMD(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 90e+0_fp

          IF ( JG == JSP ) THEN
             YMD(I,J,L)  = -90e+0_fp + ( 0.5e+0_fp * DLAT(I,J,L) )   ! S pole
          ELSE IF ( JG == JNP ) THEN
             YMD(I,J,L)  = +90e+0_fp - ( 0.5e+0_fp * DLAT(I,J,L) )   ! N pole
          ENDIF

# endif
#endif
          ! Lat centers (radians)
          YMDR(I,J,L)   = ( PI_180 * YMD(I,J,L)  )

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA )

          !----------------------------------------------------------------
          !              %%%%% FOR NESTED GRIDS ONLY %%%%%
          !----------------------------------------------------------------

          ! Lat centers (radians), for nested grid window array
          YMDRW(I,LBND+J,L) = YMDR(I,J,L)

          ! Compute YMID_R_W at edges of nested region
          IF ( J == J1 ) THEN
             YMDRW(I,LBND+J-1,1) = YMDR(I,J,L) - ( DLAT(I,J,L) * PI_180 )
          ELSE IF ( J == J2 ) THEN
             YMDRW(I,LBND+J2+1,1) = YMDR(I,J,L) + ( DLAT(I,J,L) * PI_180 )
          ENDIF
                    
#endif

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%   LATITUDE EDGES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! Lat edges (degrees and radians)
          YDG(I,J,L)    = YMD(I,J,L) - ( DLAT(I,J,L) * 0.5e+0_fp )
            
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !----------------------------------------------------------------
          !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
          !
          ! Do not define half-sized polar boxes (bmy, 3/21/13)
          !----------------------------------------------------------------
#else
          !----------------------------------------------------------------
          !            %%%%% TRADITIONAL GEOS-Chem %%%%%
          !
          ! Current practice in the standard GEOS-Chem is to force
          ! the northern edge of grid boxes along the SOUTH POLE to
          ! be -90 degrees latitude. (bmy, 3/21/13)
          !----------------------------------------------------------------
          IF ( JG == JSP ) THEN
             YDG(I,J,L) = -90e+0_fp                             ! S pole
          ENDIF
#endif          

          ! Lat edges (radians)
          YDGR(I,J,L)   = ( PI_180  * YDG(I,J,L) )
          
          ! mjc - Compute sine of latitude edges (needed for map_a2a regrid)
!------------------------------------------------------------------------------
! Prior to 3/26/15:
! Fix apparent optimization error.  Now pass a scalar to the SIN function.
! Also save the result of SIN in a scalar before storing in the YSN array.
! Don't know why this happens but this seems to fix it. (bmy, 3/26/15)
!             YSN(I,J,L) = SIN ( YDGR(I,J,L) )
!------------------------------------------------------------------------------
          YEDGE_VAL     = YDGR(I,J,L)           ! Lat edge in radians
          YSIN_VAL      = SIN( YEDGE_VAL )      ! SIN( lat edge )
          YSN(I,J,L)    = YSIN_VAL              ! Store in YSN array

       ENDDO
    ENDDO

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%   LATITUDE EDGES (the last edge)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    ! Test for North Pole
    IF ( J2 == JNP ) THEN
          
       ! North pole case (global grids only)
       DO I = I1, I2
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !----------------------------------------------------------------
          !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
          !
          ! Do not define half-sized polar boxes (bmy, 3/21/13)
          !----------------------------------------------------------------
          YDG  (I,J2+1,L)  = YDG(I,J2,L)   + DLAT(I,J2,L)
#else
          !----------------------------------------------------------------
          !            %%%%% TRADITIONAL GEOS-Chem %%%%%
          !
          ! Current practice in the standard GEOS-Chem is to force
          ! the northern edge of grid boxes along the NORTH POLE to
          ! be +90 degrees latitude. (bmy, 3/21/13)
          !----------------------------------------------------------------
          YDG (I,J2+1,L)   = +90e+0_fp
#endif
          YDGR(I,J2+1,L)   = YDG(I,J2+1,L) * PI_180

             ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
!------------------------------------------------------------------------------
! Prior to 3/26/15:
! Fix apparent optimization error.  Now pass a scalar to the SIN function.
! Also save the result of SIN in a scalar before storing in the YSN array.
! Don't know why this happens but this seems to fix it. (bmy, 3/26/15)
!             YSN(I,J2+1,L) = SIN ( YDGR(I,J2+1,L) )
!------------------------------------------------------------------------------
          YEDGE_VAL        = YDGR(I,J2+1,L)     ! Lat edge in radians
          YSIN_VAL         = SIN( YEDGE_VAL )   ! SIN( lat edge )
          YSN(I,J2+1,L)    = YSIN_VAL           ! Store in YSN array

       ENDDO
          
    ELSE
          
       !-------------------------------------------------------------------
       !                %%%%% FOR NESTED GRIDS ONLY %%%%%
       !-------------------------------------------------------------------

       ! No north pole (nested grids only)
       DO I = I1, I2
          YDG (I,J2+1,L)  = YDG(I,J2,L  ) + DLAT(I,J2,L)
          YDGR(I,J2+1,L)  = YDG(I,J2+1,L) * PI_180

          ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
!------------------------------------------------------------------------------
! Prior to 3/26/15:
! Fix apparent optimization error.  Now pass a scalar to the SIN function.
! Also save the result of SIN in a scalar before storing in the YSN array.
! Don't know why this happens but this seems to fix it. (bmy, 3/26/15)
!             YSN(I,J2+1,L) = SIN ( YDGR(I,J2+1,L) )
!------------------------------------------------------------------------------
          YEDGE_VAL       = YDGR(I,J2+1,L)
          YSIN_VAL        = SIN( YEDGE_VAL )
          YSN(I,J2+1,L)   = YSIN_VAL

       ENDDO
    ENDIF

    !======================================================================
    ! Compute grid box surface areas (algorithm from old "input.f")
    !
    ! The surface area of a grid box is derived as follows:
    ! 
    !    Area = dx * dy
    !
    ! Where:
    !
    !    dx is the arc length of the box in longitude
    !    dy is the arc length of the box in latitude
    !  
    ! Which are computed as:
    !  
    !    dx = r * delta-longitude
    !       = ( Re * cos[ YMID[J] ] ) * ( 2 * PI / IIIPAR )
    !
    !    dy = r * delta-latitude
    !       = Re * ( YEDGE[J+1] - YEDGE[J] )
    !  
    ! Where:
    !    
    !    Re         is the radius of the earth
    !    YMID[J]    is the latitude at the center of box J
    !    YEDGE[J+1] is the latitude at the N. Edge of box J
    !    YEDGE[J]   is the latitude at the S. Edge of box J
    !
    ! So, the surface area is thus:
    ! 
    !    Area = ( Re * cos( YMID[J] ) * ( 2 * PI / IIIPAR ) *
    !             Re * ( YEDGE[J+1] - YEDGE[J] )
    !
    !    2*PI*Re^2    {                                            }      
    ! = ----------- * { cos( YMID[J] ) * ( YEDGE[J+1] - YEDGE[J] ) }
    !     IIIPAR      {                                            }
    !
    ! And, by using the trigonometric identity:
    !
    !    d sin(x) = cos x * dx
    !
    ! The following term:
    !
    !    cos( YMID[J] ) * ( YEDGE[J+1] - YEDGE[J] ) 
    !
    ! May also be written as a difference of sines:
    !
    !    sin( YEDGE[J+1] ) - sin( YEDGE[J] ) 
    ! 
    ! So the final formula for surface area of a grid box is:
    ! 
    !            2*PI*Re^2    {                                     }
    !    Area = ----------- * { sin( YEDGE[J+1] ) - sin( YEDGE[J] ) }
    !              IIIPAR     {                                     }
    !
    !
    ! NOTES:
    ! (1) The formula with sines is more numerically stable, and will 
    !      yield identical global total surface areas for all grids.
    ! (2) The units are determined by the radius of the earth Re.
    !      if you use Re [m], then surface area will be in [m2], or
    !      if you use Re [cm], then surface area will be in [cm2], etc.
    ! (3) The grid box surface areas only depend on latitude, as they
    !      are symmetric in longitude.  To compute the global surface
    !      area, multiply the surface area arrays below by the number
    !      of longitudes (e.g. IIIPAR).
    ! (4) At present, assumes that GEOS-Chem will work on a
    !      Cartesian grid.
    !
    ! (bmy, 4/20/06, 2/24/12)
    !====================================================================== 

    ! Loop over local latitudes
    DO J = J1, J2

       ! Global latitude index
       JG  = J + ( J_LO - 1 )
          
       ! Loop over local longitudes
       DO I = I1, I2

          ! Sine of latitudes at N and S edges of grid box (I,J,L)
          SIN_N       = SIN( YDGR(I,J+1,L) )
          SIN_S       = SIN( YDGR(I,J,  L) )

          ! Difference of sin(latitude) at N and S edges of grid box
          SIN_DIFF    = SIN_N - SIN_S

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !-------------------------------------------------------------
          !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
          !
          ! The GEOS-5 GCM is a grid-point model, with the polar
          ! boxes (+90 & -90 degrees latitude) being grid box
          ! centers.  But since GEOS-Chem also needs to know the
          ! latitudes at the north & south edges of each grid box,
          ! we have to make a kludge.
          !
          ! For all grid boxes along the NORTH POLE (+90 lat), we
          ! let the northern edge be greater than 90 degrees.  For
          ! example, if the grid spacing is 1 degree in latitude,
          ! then for all longitudes at the North Pole:
          !
          !   * The N. EDGE of each grid box is +90.5 degrees
          !   * The CENTER  of each grid box is +90   degrees
          !   * The S. EDGE of each grid box is +89.5 degrees.
          !
          ! Similarly, at for all grid boxes along the SOUTH POLE,
          ! we have this condition (also assuming a grid spacing
          ! of 1 degree in latitude):
          !
          !   * The N. EDGE of each grid box is -89.5 degrees
          !   * The CENTER  of each grid box is -90   degrees
          !   * The S. EDGE of each grid box is -90.5 degrees.
          !
          ! Therefore, at the poles, the latitudes at the northern
          ! and southern edges of each grid box are symmetric around
          ! either +90 degrees or -90 degrees.  When you take the
          ! difference of the sine of the latitudes at the north and
          ! south edges of a polar grid box, the terms will cancel
          ! each other out, resulting in a grid box surface area
          ! that is zero.
          !
          ! We can take advantage of this symmetry around +90 and
          ! -90 degrees to make a simple fix:
          !
          ! (1) AT THE SOUTH POLE: Subtract the sine of the latitude
          !     at the north edge of the grid box from the sine
          !     of -90 degrees (which is -1) and then multiply by 2.
          !
          ! (2) AT THE NORTH POLE: Subtract the sine of +90 degrees
          !     (which is 1) from the sine of the latitude at the
          !     south edge of the grid box, and then multiply by 2.
          !
          ! This fix avoids having polar grid boxes with zero area.
          !    -- Bob Yantosca (21 Mar 2013)
          !--------------------------------------------------------------

          ! South pole kludge
          IF ( JG == JSP ) THEN
             SIN_DIFF = 2e+0_fp * ( SIN_N - ( -1e+0_fp ) )
          ENDIF

          ! North pole kludge
          IF ( JG == JNP ) THEN
             SIN_DIFF = 2e+0_fp * ( 1e+0_fp - SIN_S )
          ENDIF
#endif

          ! Grid box surface areas [m2]
          AM2(I,J,L) = ( DLON(I,J,L) * PI_180 ) * ( Re**2 ) * SIN_DIFF

       ENDDO
    ENDDO

!    ! Return successfully
!    RC = GIGC_SUCCESS

    !======================================================================
    ! Echo info to stdout
    !======================================================================
    IF ( am_I_Root ) THEN
       WRITE( 6, '(''Nested-Grid X-offset [boxes]:'', i4 )' ) IOFF
       WRITE( 6, '(''Nested-Grid Y-offset [boxes]:'', i4 )' ) JOFF
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( XMD(I,1,1),  I=1,I2-I1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( XDG(I,1,1), I=1,I2+1-I1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( YMD(1,J,1),  J=1,J2-J1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( YDG(1,J,1), J=1,J2+1-J1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''SIN( grid box latitude edges )'')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( YSN(1,J,1), J=1,J2+1-J1+1 )
    ENDIF

  END SUBROUTINE DoGridComputation 
!EOC
END MODULE Grid_Mod
