
    SUBROUTINE run_prosail ( N, Cab, Car, Cbrown, Cw, Cm, lai, LIDFa, LIDFb, &
                rsoil, psoil, hspot, tts, tto, psi, TypeLidf, retval, &
                soil_spectrum1, soil_spectrum2 )

    USE MOD_ANGLE               ! defines pi & rad conversion
    USE MOD_staticvar           ! static variables kept in memory for optimization
    USE MOD_flag_util           ! flags for optimization
    USE MOD_output_PROSPECT     ! output variables of PROSPECT
    USE MOD_SAIL                ! variables of SAIL
    USE MOD_dataSpec_P5B        
    IMPLICIT NONE

    ! LEAF BIOCHEMISTRY
    REAL*8, dimension(nw), intent(out) :: retval
    REAL*8, intent(in) :: N,Cab,Car,Cbrown,Cw,Cm
    ! CANOPY
    REAL*8, intent(in) :: lai,LIDFa,LIDFb,psoil,rsoil
    REAL*8, intent(in) :: hspot
    REAL*8, intent(in) :: tts,tto,psi
    INTEGER, intent(in) :: TypeLidf   
    REAL*8,ALLOCATABLE,SAVE :: resh(:),resv(:)
    REAL*8,ALLOCATABLE,SAVE :: rsoil0(:),PARdiro(:),PARdifo(:)
    REAL*8, dimension(nw), intent(in), optional :: soil_spectrum1
    REAL*8, dimension(nw), intent(in), optional :: soil_spectrum2
    INTEGER :: ii
    REAL*8 :: ihot, skyl
    ! ANGLE CONVERSION
    pi=3.151592d0
    rd=pi/180.d0

    ! PROSPECT output
    ALLOCATE (LRT(nw,2),rho(nw),tau(nw))
    ! SAIL
    ALLOCATE (sb(nw),sf(nw),vb(nw),vf(nw),w(nw))
    ALLOCATE (m(nw),m2(nw),att(nw),sigb(nw),rinf(nw))
    ALLOCATE (PARdiro(nw),PARdifo(nw))
    ALLOCATE(tsd(nw),tdd(nw),tdo(nw),rsd(nw),rdd(nw),rso(nw),rdo(nw))
    ALLOCATE(rddt(nw),rsdt(nw),rdot(nw),rsodt(nw),rsost(nw),rsot(nw),rsos(nw),rsod(nw))
    ALLOCATE(lidf(13))
    ! resh : hemispherical reflectance
    ! resv : directional reflectance
    ALLOCATE (resh(nw),resv(nw))
    ALLOCATE (rsoil_old(nw))

        !TypeLidf=1
        ! if 2-parameters LIDF: TypeLidf=1
        !!!IF (TypeLidf.EQ.1) THEN
            ! LIDFa LIDF parameter a, which controls the average leaf slope
            ! LIDFb LIDF parameter b, which controls the distribution's bimodality
            !   LIDF type       a        b
            !   Planophile      1        0
            !   Erectophile    -1        0
            !   Plagiophile     0       -1
            !   Extremophile    0        1
            !   Spherical      -0.35    -0.15
            !   Uniform 0 0
            !   requirement: |LIDFa| + |LIDFb| < 1  
           !!! LIDFa   =   -0.35
           !!! LIDFb   =   -0.15
        ! if ellipsoidal distribution: TypeLidf=2
        !!! ELSEIF (TypeLidf.EQ.2) THEN
            !   LIDFa   = average leaf angle (degrees) 0 = planophile   /   90 = erectophile
            !   LIDFb = 0
           !!! LIDFa   =   30
           !!! LIDFb   =   0
        !!!ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!LEAF CHEM & STR PROPERTIES!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! INITIAL PARAMETERS
!        Cab     =   40.     ! chlorophyll content (µg.cm-2) 
!        Car     =   8.      ! carotenoid content (µg.cm-2)
!        Cbrown  =   0.0     ! brown pigment content (arbitrary units)
!        Cw      =   0.01    ! EWT (cm)
!        Cm      =   0.009   ! LMA (g.cm-2)
!        N       =   1.5     ! structure coefficient

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  Soil Reflectance Properties !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! rsoil1 = dry soil
        ! rsoil2 = wet soil
        ALLOCATE (rsoil0(nw))
        !psoil   =   1.      ! soil factor (psoil=0: wet soil / psoil=1: dry soil)
        ! rsoil : soil brightness  term
        if ( present(soil_spectrum1) ) then
            Rsoil1 = soil_spectrum1
        endif
        if ( present ( soil_spectrum2 ) ) then
            Rsoil2 = soil_spectrum2
        endif
        rsoil0=rsoil*(psoil*Rsoil1+(1-psoil)*Rsoil2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  4SAIL canopy structure parm !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        LAI     =   3.      ! leaf area index (m^2/m^2)
!        hspot   =   0.01    ! hot spot
!        tts     =   30.     ! solar zenith angle (°)
!        tto     =   10.     ! observer zenith angle (°)
!        psi     =   0.      ! azimuth (°)

        init_completed=.false.  ! only at first call of PRO4SAIL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!        CALL PRO4SAIL         !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL PRO4SAIL(N,Cab,Car,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLIDF,LAI,hspot,tts,tto,psi,rsoil0)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  direct / diffuse light  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! the direct and diffuse light are taken into account as proposed by:
        ! Francois et al. (2002) Conversion of 4001100 nm vegetation albedo 
        ! measurements into total shortwave broadband albedo using a canopy 
        ! radiative transfer model, Agronomie
        !skyl    =   0.847- 1.61*sin((90-tts)*rd)+ 1.04*sin((90-tts)*rd)*sin((90-tts)*rd) ! % diffuse radiation
        ! Es = direct
        ! Ed = diffuse
        ! PAR direct
        !PARdiro =   (1-skyl)*Es
        ! PAR diffus
        !PARdifo =   (skyl)*Ed
        ! resv : directional reflectance
        
        !retval = (rdot*PARdifo+rsot*PARdiro)/(PARdiro+PARdifo)
        retval = rsot
        deallocate ( lrt )
        deallocate ( rho )
        deallocate ( tau )
        deallocate ( sb, sf, vb, vf, w )
        deallocate ( m, m2, att, sigb, rinf )
        deallocate ( pardiro, pardifo )
        deALLOCATE(tsd,tdd,tdo,rsd,rdd,rso,rdo)
        deALLOCATE(rddt,rsdt,rdot,rsodt,rsost,rsot,rsos,rsod)
        deALLOCATE(lidf)
        ! resh : hemispherical reflectance
        ! resv : directional reflectance
        deALLOCATE (resh,resv)
        deALLOCATE (rsoil_old)
        deALLOCATE (rsoil0)
        !WRITE( *,'(i4,f10.6)') (lambda(ii),run_prosail(ii), ii=1,nw)
    END subroutine run_prosail
    

    
    SUBROUTINE run_sail ( refl, trans, lai, LIDFa, LIDFb, &
                rsoil, psoil, hspot, tts, tto, psi, TypeLidf, retval )

    USE MOD_ANGLE               ! defines pi & rad conversion
    USE MOD_staticvar           ! static variables kept in memory for optimization
    USE MOD_flag_util           ! flags for optimization
    USE MOD_output_PROSPECT     ! output variables of PROSPECT
    USE MOD_SAIL                ! variables of SAIL
    USE MOD_dataSpec_P5B        
    IMPLICIT NONE

    ! LEAF BIOCHEMISTRY
    REAL*8, dimension(nw), intent(out) :: retval
    REAL*8, dimension(nw), intent(in) :: refl, trans
    ! CANOPY
    REAL*8, intent(in) :: lai,LIDFa,LIDFb,psoil,rsoil
    REAL*8, intent(in) :: hspot
    REAL*8, intent(in) :: tts,tto,psi
    INTEGER, intent(in) :: TypeLidf   
    REAL*8,ALLOCATABLE,SAVE :: resh(:),resv(:)
    REAL*8,ALLOCATABLE,SAVE :: rsoil0(:),PARdiro(:),PARdifo(:)
    INTEGER :: ii
    REAL*8 :: ihot, skyl
    ! ANGLE CONVERSION
    pi=3.151592d0
    rd=pi/180.d0

    ! PROSPECT output
    ALLOCATE (LRT(nw,2),rho(nw),tau(nw))
    ! SAIL
    ALLOCATE (sb(nw),sf(nw),vb(nw),vf(nw),w(nw))
    ALLOCATE (m(nw),m2(nw),att(nw),sigb(nw),rinf(nw))
    ALLOCATE (PARdiro(nw),PARdifo(nw))
    ALLOCATE(tsd(nw),tdd(nw),tdo(nw),rsd(nw),rdd(nw),rso(nw),rdo(nw))
    ALLOCATE(rddt(nw),rsdt(nw),rdot(nw),rsodt(nw),rsost(nw),rsot(nw),rsos(nw),rsod(nw))
    ALLOCATE(lidf(13))
    ! resh : hemispherical reflectance
    ! resv : directional reflectance
    ALLOCATE (resh(nw),resv(nw))
    ALLOCATE (rsoil_old(nw))

        !TypeLidf=1
        ! if 2-parameters LIDF: TypeLidf=1
        !!!IF (TypeLidf.EQ.1) THEN
            ! LIDFa LIDF parameter a, which controls the average leaf slope
            ! LIDFb LIDF parameter b, which controls the distribution's bimodality
            !   LIDF type       a        b
            !   Planophile      1        0
            !   Erectophile    -1        0
            !   Plagiophile     0       -1
            !   Extremophile    0        1
            !   Spherical      -0.35    -0.15
            !   Uniform 0 0
            !   requirement: |LIDFa| + |LIDFb| < 1  
           !!! LIDFa   =   -0.35
           !!! LIDFb   =   -0.15
        ! if ellipsoidal distribution: TypeLidf=2
        !!! ELSEIF (TypeLidf.EQ.2) THEN
            !   LIDFa   = average leaf angle (degrees) 0 = planophile   /   90 = erectophile
            !   LIDFb = 0
           !!! LIDFa   =   30
           !!! LIDFb   =   0
        !!!ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!LEAF CHEM & STR PROPERTIES!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! INITIAL PARAMETERS
!        Cab     =   40.     ! chlorophyll content (µg.cm-2) 
!        Car     =   8.      ! carotenoid content (µg.cm-2)
!        Cbrown  =   0.0     ! brown pigment content (arbitrary units)
!        Cw      =   0.01    ! EWT (cm)
!        Cm      =   0.009   ! LMA (g.cm-2)
!        N       =   1.5     ! structure coefficient

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  Soil Reflectance Properties !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! rsoil1 = dry soil
        ! rsoil2 = wet soil
        ALLOCATE (rsoil0(nw))
        !psoil   =   1.      ! soil factor (psoil=0: wet soil / psoil=1: dry soil)
        ! rsoil : soil brightness  term
        rsoil0=rsoil*(psoil*Rsoil1+(1-psoil)*Rsoil2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  4SAIL canopy structure parm !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        LAI     =   3.      ! leaf area index (m^2/m^2)
!        hspot   =   0.01    ! hot spot
!        tts     =   30.     ! solar zenith angle (°)
!        tto     =   10.     ! observer zenith angle (°)
!        psi     =   0.      ! azimuth (°)

        init_completed=.false.  ! only at first call of PRO4SAIL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!        CALL PRO4SAIL         !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL SAIL(refl, trans,LIDFa,LIDFb,TypeLIDF,LAI,hspot,tts,tto,psi,rsoil0)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  direct / diffuse light  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! the direct and diffuse light are taken into account as proposed by:
        ! Francois et al. (2002) Conversion of 4001100 nm vegetation albedo 
        ! measurements into total shortwave broadband albedo using a canopy 
        ! radiative transfer model, Agronomie
        !skyl    =   0.847- 1.61*sin((90-tts)*rd)+ 1.04*sin((90-tts)*rd)*sin((90-tts)*rd) ! % diffuse radiation
        ! Es = direct
        ! Ed = diffuse
        ! PAR direct
        !PARdiro =   (1-skyl)*Es
        ! PAR diffus
        !PARdifo =   (skyl)*Ed
        ! resv : directional reflectance
        
        !retval = (rdot*PARdifo+rsot*PARdiro)/(PARdiro+PARdifo)
        retval = rsot
        deallocate ( lrt )
        deallocate ( rho )
        deallocate ( tau )
        deallocate ( sb, sf, vb, vf, w )
        deallocate ( m, m2, att, sigb, rinf )
        deallocate ( pardiro, pardifo )
        deALLOCATE(tsd,tdd,tdo,rsd,rdd,rso,rdo)
        deALLOCATE(rddt,rsdt,rdot,rsodt,rsost,rsot,rsos,rsod)
        deALLOCATE(lidf)
        ! resh : hemispherical reflectance
        ! resv : directional reflectance
        deALLOCATE (resh,resv)
        deALLOCATE (rsoil_old)
        deALLOCATE (rsoil0)
        !WRITE( *,'(i4,f10.6)') (lambda(ii),run_prosail(ii), ii=1,nw)
    END subroutine run_sail
    
    