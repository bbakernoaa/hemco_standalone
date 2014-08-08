!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: grid_mod.F90 
!
! !DESCRIPTION: Module grid\_mod contains routines to calculate the grid
! in the standalone HEMCO environment. This module is adapted from the
! corresponding module of GEOS-Chem.
!
! !INTERFACE:
!
MODULE GRID_MOD
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: DoGridComputation
!
! !REVISION HISTORY:
!  28 Jul 2014 - C. Keller   - Initial version (moved from hcoi_standalone_mod.F90) 
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DoGridComputation 
!
! !DESCRIPTION: Subroutine DoGridComputation initializes the longitude, 
!  latitude and surface area arrays. This used to be subroutine COMPUTE\_GRID.
!\\
!\\
!  NOTE: This routine was taken out of the GEOS-Chem module grid\_mod.F90
!  and moved into hcoi\_standalone\_mod.F90.  This prevents us from having 
!  to rely on GEOS-Chem modules in the HEMCO standalone build.
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
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root                     ! Root CPU?

    ! Variables with local CPU indices
    INTEGER, INTENT(IN)  :: I1,  I2                       ! Min lon index
    INTEGER, INTENT(IN)  :: J1,  J2                       ! Local lat indices
    INTEGER, INTENT(IN)  :: JSP, JNP                      ! Polar lat indices
    INTEGER, INTENT(IN)  :: L1,  L2                       ! Local lev indices
    REAL*8,  INTENT(IN)  :: DLON(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lon [deg]
    REAL*8,  INTENT(IN)  :: DLAT(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lat [deg]

    ! Variables with global CPU indices
    INTEGER, INTENT(IN)  :: I_LO                          ! Min global lon
    INTEGER, INTENT(IN)  :: J_LO                          ! Min global lat

    ! Offsets (for nested grids) 
    INTEGER, INTENT(IN)  :: IOFF
    INTEGER, INTENT(IN)  :: JOFF

    ! Arrays to be filled
    REAL*8,  INTENT(OUT) :: XMD  (:,:,:) 
    REAL*8,  INTENT(OUT) :: XDG  (:,:,:) 
    REAL*8,  INTENT(OUT) :: YMD  (:,:,:) 
    REAL*8,  INTENT(OUT) :: YDG  (:,:,:) 
    REAL*8,  INTENT(OUT) :: YSN  (:,:,:) 
    REAL*8,  INTENT(OUT) :: YMDR (:,:,:) 
    REAL*8,  INTENT(OUT) :: YDGR (:,:,:) 
    REAL*8,  INTENT(OUT) :: YMDRW(:,:,:) 
    REAL*8,  INTENT(OUT) :: YDGRW(:,:,:) 
    REAL*8,  INTENT(OUT) :: AM2  (:,:,:) 
!
! !OUTPUT PARAMETERS:
!  
    INTEGER, INTENT(OUT) :: RC                            ! Success or failure?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: I,     J,      L,        IG, JG
    REAL*8  :: SIN_N, SIN_S,  SIN_DIFF, TMP

    ! Arrays
    REAL*8  :: IND_X( I1:I2+1 )
    REAL*8  :: IND_Y( J1:J2+1 )

    ! Parameter
    REAL*8, PARAMETER :: PI     = 3.14159265358979323d0
    REAL*8, PARAMETER :: Re     = 6.375d6
    REAL*8, PARAMETER :: PI_180 = PI / 180d0

    !======================================================================
    ! Initialization
    !======================================================================

    ! Index array for longitudes
    DO I = I1, I2+1
       IND_X(I) = ( ( I + IOFF - 1 ) * 1d0 ) + ( I_LO - 1 )
    ENDDO

    ! Index array for latitudes
    DO J = J1, J2+1
       IND_Y(J) = ( ( J + JOFF - 1 ) * 1d0 ) + ( J_LO - 1 )
    ENDDO

    !======================================================================
    ! Compute longitude and latitude arrays
    !======================================================================
    
    ! Loop over levels
    DO L = L1, L2
       
       !-------------------------------------------------------------------
       ! Longitude center and edge arrays
       !-------------------------------------------------------------------

       ! Loop over local latitudes
       DO J = J1, J2

          ! Loop over local longitudes
          DO I = I1, I2

             ! Longitude centers
             XMD(I,J,L)  = ( DLON(I,J,L) * IND_X(I) ) - 180d0
          
             ! Longitude edges
             XDG(I,J,L) = XMD(I,J,L) - ( DLON(I,J,L) * 0.5d0 )

             ! Compute the last longitude edge
             IF ( I == I2 ) THEN
                XDG(I+1,J,L) = XDG(I,J,L) + DLON(I,J,L)
             ENDIF
             
          ENDDO
       ENDDO

       !-------------------------------------------------------------------
       ! Latitude center and edge arrays
       !-------------------------------------------------------------------

       ! Loop over local latitudes
       DO J = J1, J2

          ! Global latitude index
          JG = J + ( J_LO - 1 )

          ! Loop over local longitudes
          DO I = I1, I2

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !%%%%%%   LATITUDE CENTERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!             ! Lat centers (degrees)
!             YMD(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 90d0
          
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
             !-------------------------------------------------------------
             !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
             !
             ! Do not define half-sized polar boxes (bmy, 3/21/13)
             !-------------------------------------------------------------
#else

# if defined( GCAP )

             !-------------------------------------------------------------
             !             %%%%% TRADITIONAL GEOS-Chem %%%%%
             !
             ! For the GCAP model, there are no half-size polar boxes.
             ! Compute the latitude centers accordingly.  (bmy, 7/2/13)
             !-------------------------------------------------------------

             ! Lat centers (degrees)
             YMD(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 88d0

#else

             !-------------------------------------------------------------
             !             %%%%% TRADITIONAL GEOS-Chem %%%%%
             !
             ! Current practice in the standard GEOS-Chem is to make
             ! the polar grid boxes (for non-GCAP global grids only) be
             ! half the size of other grid boxes.  This lets us make +90
             ! degrees and -90 degrees be the edges of the grid.
             ! (bmy, 7/2/13)
             !-------------------------------------------------------------

             ! Lat centers (degrees)
             YMD(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 90d0

             IF ( JG == JSP ) THEN
                YMD(I,J,L)  = -90d0 + ( 0.5d0 * DLAT(I,J,L) )   ! S pole
             ELSE IF ( JG == JNP ) THEN
                YMD(I,J,L)  = +90d0 - ( 0.5d0 * DLAT(I,J,L) )   ! N pole
             ENDIF

# endif
#endif
             ! Lat centers (radians)
             YMDR(I,J,L)   = ( PI_180 * YMD(I,J,L)  )

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA )

             !-------------------------------------------------------------
             !              %%%%% FOR NESTED GRIDS ONLY %%%%%
             !-------------------------------------------------------------

             ! Lat centers (radians), for nested grid window array
             YMDRW(I,J,L) = YMDR(I,J,L)

             ! Compute YMID_R_W at edges of nested region
             IF ( J == J1 ) THEN
                !YMID_R_W(I,J1-1,1) = YMID_R(I,J1,L) - ( DLAT(I,J1,L) * PI_180 )
                YMDRW(I,J-1,1) = YMDR(I,J,L) - ( DLAT(I,J,L) * PI_180 )
             ELSE IF ( J == J2 ) THEN
                !YMID_R_W(I,J2+1,1) = YMID_R(I,J2,L) + ( DLAT(I,J2,L) * PI_180 )
                YMDRW(I,J+1,1) = YMDR(I,J,L) + ( DLAT(I,J,L) * PI_180 )
             ENDIF
                    
#endif

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !%%%%%%   LATITUDE EDGES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             ! Lat edges (degrees and radians)
             YDG(I,J,L)    = YMD(I,J,L) - ( DLAT(I,J,L) * 0.5d0 )
            
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
             !-----------------------------------------------------------
             !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
             !
             ! Do not define half-sized polar boxes (bmy, 3/21/13)
             !-----------------------------------------------------------
#else
             !-----------------------------------------------------------
             !            %%%%% TRADITIONAL GEOS-Chem %%%%%
             !
             ! Current practice in the standard GEOS-Chem is to force
             ! the northern edge of grid boxes along the SOUTH POLE to
             ! be -90 degrees latitude. (bmy, 3/21/13)
             !-----------------------------------------------------------
             IF ( JG == JSP ) THEN
                YDG(I,J,L) = -90d0                             ! S pole
             ENDIF
#endif          

             ! Lat edges (radians)
             YDGR(I,J,L)  = ( PI_180  * YDG(I,J,L) )

             ! mjc - Compute sine of latitude edges (needed for map_a2a regrid)
             YSN(I,J,L) = SIN ( YDGR(I,J,L) )
             
          ENDDO
       ENDDO

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%%%%   LATITUDE EDGES (the last edge)   %%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       ! Test for North Pole
       IF ( J2 == JNP ) THEN
          
          ! North pole case (global grids only)
          DO I = I1, I2
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
             !-----------------------------------------------------------
             !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
             !
             ! Do not define half-sized polar boxes (bmy, 3/21/13)
             !-----------------------------------------------------------
             YDG  (I,J2+1,L)   = YDG(I,J2,L)   + DLAT(I,J2,L)
#else
             !-----------------------------------------------------------
             !            %%%%% TRADITIONAL GEOS-Chem %%%%%
             !
             ! Current practice in the standard GEOS-Chem is to force
             ! the northern edge of grid boxes along the NORTH POLE to
             ! be +90 degrees latitude. (bmy, 3/21/13)
             !-----------------------------------------------------------
             YDG (I,J2+1,L)   = +90d0
#endif
             YDGR(I,J2+1,L)   = YDG(I,J2+1,L) * PI_180

             ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
             YSN(I,J2+1,L) = SIN ( YDGR(I,J2+1,L) )
          ENDDO
          
       ELSE
          
         !---------------------------------------------------------------
         !                %%%%% FOR NESTED GRIDS ONLY %%%%%
         !---------------------------------------------------------------

          ! No north pole (nested grids only)
          DO I = I1, I2
             YDG (I,J2+1,L)  = YDG(I,J2,L  ) + DLAT(I,J2,L)
             YDGR(I,J2+1,L)  = YDG(I,J2+1,L) * PI_180

             ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
             YSN(I,J2+1,L) = SIN ( YDGR(I,J2+1,L) )
          ENDDO
       ENDIF

       !=================================================================
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
       !=================================================================  

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
             !-----------------------------------------------------------
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
             !-----------------------------------------------------------

             ! South pole kludge
             IF ( JG == JSP ) THEN
                SIN_DIFF = 2d0 * ( SIN_N - ( -1d0 ) )
             ENDIF

             ! North pole kludge
             IF ( JG == JNP ) THEN
                SIN_DIFF = 2d0 * ( 1d0 - SIN_S )
             ENDIF
#endif

             ! Grid box surface areas [m2]
             AM2(I,J,L) = ( DLON(I,J,L) * PI_180 ) * ( Re**2 ) * SIN_DIFF

          ENDDO
       ENDDO

    ENDDO


    !=================================================================
    ! Echo info to stdout
    !=================================================================
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
    ENDIF

    ! Return w/ success

  END SUBROUTINE DoGridComputation 
!EOC
END MODULE GRID_MOD
