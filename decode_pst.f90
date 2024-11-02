      SUBROUTINE GTR_DECODE_PST(IDATA,PIXEL,HIT)
!-------------------------------------------------------------------------------
!     decode pattern trigger information 
!-------------------------------------------------------------------------------
!     Routine to take trigger words from 13 Hytec Pattern Selection
!     Trigger (PST) modules oriented to cover a 331 pixel hexagon and 
!     decipher which pixels caused a trigger. Each PST MODULE covers
!     59 pixels arranged in FIVE overlapping PATCHES of 19 PIXELS.
!-------------------------------------------------------------------------------

      IMPLICIT NONE
      
!---- arguments 
      LOGICAL, INTENT(IN)      :: IDATA(:)       ! data word array 
      INTEGER, POINTER         :: PIXEL(:,:)     ! pixel(u,v) look up table
      TYPE(GTR_HIT_T), POINTER :: HIT(:)         ! decoding result 

!---- 
      INTEGER  :: I,I1,I2                       ! loop counter abd limits 
      INTEGER  :: J,J1,J2                       ! loop counter abd limits 
      INTEGER  :: I32 

!---- declarations specifically for decoding position ------

      INTEGER N,M
      INTEGER IMODULE, IPATCH   ! module, patch
      INTEGER IPTCHU(19)        ! u-position within patch
      INTEGER IPTCHV(19)        ! v-position within patch
      INTEGER IMIDPTCHU(5)      ! u-position patch center
      INTEGER IMIDPTCHV(5)      ! v-position patch center
      INTEGER U,V
      INTEGER X,Y
      INTEGER NH
      INTEGER A,B,C,D           ! rotation matrix coeff 
      INTEGER NP                ! number of patch bits set 
 
 
!---- description of modules 
      TYPE :: MODULE_T  
        INTEGER U,V                    ! position of center 
        INTEGER A,B,C,D                ! rotation matrix coefficients 
      END TYPE 
      TYPE(MODULE_T) :: MODULE(13)  
  

!---- pixel position for each of 19 bits 
      DATA IPTCHU / 0,-1,-2,-2,-2,+1,0,-1,-1,-1,+2,+1,0,0,0,+2,1,+1,+2/
      DATA IPTCHV /-2,-1,0,+1,+2,-2,-1,0,+1,+2,-2,-1,0,+1,+2,-1,0,+1,0/

!---- patch position within module for each of 5 patches  
      DATA IMIDPTCHU / -4, -2, 0, +2, +4 /
      DATA IMIDPTCHV /  0,  0, 0,  0,  0 /

!---- position of module centre (channel 13 of patch 3) 
      DATA MODULE%U     /+4,+4,+4,+4,+4,-2, 0,+2,+4,-2,-4,-6,-8/
      DATA MODULE%V     / 0,-2, 0,+2,+4,-2,-4,-6,-8,+4,+4,+4,+4/
!     module rotation        <+60 deg >, < 180 deg >,  <-60 deg>   
      DATA MODULE%A    / 1, +1,+1,+1,+1, -1,-1,-1,-1,  0, 0, 0, 0 /  
      DATA MODULE%B    / 0, +1,+1,+1,+1,  0, 0, 0, 0, -1,-1,-1,-1 /  
      DATA MODULE%C    / 0, -1,-1,-1,-1,  0, 0, 0, 0, +1,+1,+1,+1 /  
      DATA MODULE%D    / 1,  0, 0, 0, 0, -1,-1,-1,-1, +1,+1,+1,+1 /  

!---- 
      INTEGER                      :: IERR 
      CHARACTER(LEN=132)           ::  MSG                      ! error message 
      CHARACTER(LEN=  *),PARAMETER ::  SRN = "GTR_DECODE_PAST" 
      LOGICAL                      :: TEST(5) 
!------------------------------------------------------------------------------
!     allocate memory  
!------------------------------------------------------------------------------
      NH = 0                                       ! reset bit set count 
      I1 = LBOUND(IDATA,1)                         ! lower index data words 
      I2 = UBOUND(IDATA,1)                         ! upper .... 
      DO I=I1,I2                                   ! all data words 
        DO J=1,19                                  ! all 19 channel bits
          I32 = IDATA(I)                           ! copy to integer 
          IF (BTEST(I32,J-1)) NH = NH + 1          ! count bits set = hits   
        ENDDO                                      ! next bit 
      ENDDO                                        ! next data word 
      IF (ASSOCIATED(HIT)) DEALLOCATE(HIT)         ! remove previous result  
      IF (NH.LE.0) RETURN                          ! any new hits or return?   
      ALLOCATE(HIT(NH))                            ! alloc max memory needed
!------------------------------------------------------------------------------
!--- NOW DECODE THE DATA WORDS INTO (u,v) PIXEL POSITIONS -------
!------------------------------------------------------------------------------
      NH = 0                                       ! reset hit count 
      DO I=I1,I2                                   ! all data words 

        I32     = IDATA(I)                         ! copy to integer 
        IMODULE = IAND(ISHFT(I32,-24),'7F'X)       ! decode module number
        IPATCH  = IAND(ISHFT(I32,-19),'1F'X)       ! get patch number bits  
        NP = 0                                     ! set bit count to zero
        DO J=0,4                                   ! all bit positions 
          IF (BTEST(IPATCH,J)) NP = NP+1           ! count bits 
        ENDDO  
        DO J=0,4                                   ! go through all five bits 
          IF (BTEST(IPATCH,J)) THEN                ! test if bit is on? 
            IPATCH = J+1                           ! if yes, set patch #  
            EXIT                                   ! found bit, exit loop 
          ENDIF                                    ! end test bit 
        ENDDO                                      ! try next bit 

        TEST(1) = IMODULE.LT.LBOUND(MODULE,1)      ! test module number 
        TEST(2) = IMODULE.GT.UBOUND(MODULE,1) 
        TEST(3) = IPATCH .LT.LBOUND(IPTCHU,1)      ! test patch number 
        TEST(4) = IPATCH .GT.UBOUND(IPTCHU,1) 
        TEST(5) = (NP.LT.1).OR.(NP.GT.5)           ! test single patch only 
        
        IF (ANY(TEST)) THEN 
          WRITE(MSG,*) "Illegal patch or module number",IPATCH,IMODULE,NP  
	  IERR = GOS_ERROR(SRN,TRIM(MSG),1)       ! send message 
          DEALLOCATE(HIT) 
          RETURN 
        ENDIF 

        DO J=1,19                                  ! all 19 channel bits
          I32     = IDATA(I)                       ! copy to integer 
          IF (BTEST(I32,J-1)) THEN                 ! test bit on 
            NH = NH + 1                            ! increment hits found
            HIT(NH)%MODULE = IMODULE               ! store module number
            HIT(NH)%PATCH  = IPATCH                ! store patch number
            X = IPTCHU(J)+IMIDPTCHU(IPATCH)        ! x-position within a module  
            Y = IPTCHV(J)+IMIDPTCHV(IPATCH)        ! y-position 
            A = MODULE(IMODULE)%A                  ! get rotation matrix
            B = MODULE(IMODULE)%B                  !   coefficients for 
            C = MODULE(IMODULE)%C                  !   this module 
            D = MODULE(IMODULE)%D                  !   ready to use 
            U = A*X + B*Y                          ! rotate, module 
            V = C*X + D*Y                          ! rotate, module 
            U = U + MODULE(IMODULE)%U              ! translate module in u  
            V = V + MODULE(IMODULE)%V              ! translate module in v 

! in November 1999 RWL/SMB inverted Granite display to match HV and CCD
! U and V co-ordinates of pixels changed as a result, so now trying to
! get PST display to match ADC display:

            U = -1*U
            V = -1*V

            HIT(NH)%UPOS  = U                      ! store u-position 
            HIT(NH)%VPOS  = V                      ! store v-position 
            HIT(NH)%PIXEL = PIXEL(U,V)             ! store pixel number
          ENDIF
        ENDDO                                      ! next bit

      ENDDO
         
      RETURN
      END SUBROUTINE GTR_DECODE_PST

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
