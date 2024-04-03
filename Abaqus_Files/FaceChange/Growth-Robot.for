C234567890123456789012345678901234567890123456789012345678901234567890
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C LOCAL ARRAYS
C ----------------------------------------------------------------
C AA - Elastic tensor at the end of the increment
C BAA - Left Cauchy-green elastic tensor at the end of the increment
C BAAB - Deviatoric left Cauchy-green elastic tensor at the end of the increment
C TRBAAB - Trace of BAAB
C Some geometry quantities:
C       Er, Fr, Gr, Lr, Mr, Nr are related to reference surface
C       EE, FF, GG, LL, MM, NN are related to target surface
C ----------------------------------------------------------------
C
        PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, 
     &    FOUR=4.D0, SIX=6.D0, NELEMENT=16030)

!      TODO: please modify the element number [NELEMENT]

        DOUBLE PRECISION :: AA(3,3), BAA(6), BAAB(6),
     &  Lambda1z0ImpSt1(NELEMENT), Lambda2z0ImpSt1(NELEMENT),
     &  LambdaSz0ImpSt1(NELEMENT), 
     &  Lambda1z1ImpSt1(NELEMENT), Lambda2z1ImpSt1(NELEMENT),
     &  LambdaSz1ImpSt1(NELEMENT), 
     &  a1ImpSt1(NELEMENT), a2ImpSt1(NELEMENT), a3ImpSt1(NELEMENT), 
     &  b1ImpSt1(NELEMENT), b2ImpSt1(NELEMENT), b3ImpSt1(NELEMENT), 
     &  Lambda1z0ImpSt2(NELEMENT), Lambda2z0ImpSt2(NELEMENT),
     &  LambdaSz0ImpSt2(NELEMENT), 
     &  Lambda1z1ImpSt2(NELEMENT), Lambda2z1ImpSt2(NELEMENT),
     &  LambdaSz1ImpSt2(NELEMENT), 
     &  a1ImpSt2(NELEMENT), a2ImpSt2(NELEMENT), a3ImpSt2(NELEMENT), 
     &  b1ImpSt2(NELEMENT), b2ImpSt2(NELEMENT), b3ImpSt2(NELEMENT)

        DOUBLE PRECISION :: Lambda1z0, Lambda1z1, Lambda2z0, Lambda2z1,
     &  LambdaSz0, LambdaSz1, Pi, TotalT,  X, Y, Z, 
     &  a1, a2, a3, b1, b2, b3, 
     &  FinalG11St1, FinalG12St1, FinalG13St1, 
     &  FinalG21St1, FinalG22St1, FinalG23St1, 
     &  FinalG31St1, FinalG32St1, FinalG33St1, 
     &  FinalG11St2, FinalG12St2, FinalG13St2, 
     &  FinalG21St2, FinalG22St2, FinalG23St2, 
     &  FinalG31St2, FinalG32St2, FinalG33St2, 
     &  Lam11, Lam12, Lam13, 
     &  Lam21, Lam22, Lam23, 
     &  Lam31, Lam32, Lam33, 
     &  G11, G12, G13, G21, G22, G23, G31, G32, G33, 
     &  ThickCoord, TotalTh

        SAVE Lambda1z0ImpSt1
        SAVE Lambda2z0ImpSt1
        SAVE LambdaSz0ImpSt1
        SAVE Lambda1z1ImpSt1
        SAVE Lambda2z1ImpSt1
        SAVE LambdaSz1ImpSt1

        SAVE Lambda1z0ImpSt2
        SAVE Lambda2z0ImpSt2
        SAVE LambdaSz0ImpSt2
        SAVE Lambda1z1ImpSt2
        SAVE Lambda2z1ImpSt2
        SAVE LambdaSz1ImpSt2

        SAVE a1ImpSt1
        SAVE a2ImpSt1
        SAVE a3ImpSt1
        SAVE b1ImpSt1
        SAVE b2ImpSt1
        SAVE b3ImpSt1

        SAVE a1ImpSt2
        SAVE a2ImpSt2
        SAVE a3ImpSt2
        SAVE b1ImpSt2
        SAVE b2ImpSt2
        SAVE b3ImpSt2

        data iread /1/
        SAVE iread
C
C ----------------------------------------------------------------
C UMAT FOR COMPRESSIBLE NEO-HOOKEAN HYPERELASTICITY
C CANNOT BE USED FOR PLANE STRESS
C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C STATEV(1) - Initial coordinate COORDS(1)
C STATEV(2) - Initial coordinate COORDS(2)
C STATEV(3) - Initial coordinate COORDS(3)

C STATEV(9) - Mod of growth tensor
C ----------------------------------------------------------------
C
C ELASTIC PROPERTIES
        EMOD=PROPS(1)
        ENU=PROPS(2)
        C10=EMOD/(FOUR*(ONE+ENU))
        D1=SIX*(ONE-TWO*ENU)/EMOD
        Pi=3.14159265359
C get coordinate from COORDS variable
      !   IF (STATEV(1) .EQ. 0) THEN
      !           STATEV(1)=COORDS(1)
      !   END IF
      !   IF (STATEV(2) .EQ. 0) THEN
      !           STATEV(2)=COORDS(2)
      !   END IF
      !   IF (STATEV(3) .EQ. 0) THEN
      !           STATEV(3)=COORDS(3)
      !   END IF
C store coordinate from STATEV
      !   X = STATEV(1)
      !   Y = STATEV(2)
      !   Z = STATEV(3)

!I cannot get total time from UMAT variable, please modify TotalT manually 
        TotalT=10.0

C calculate the thickness of the current Gauss point
! TODO: please modify the total thickness on Gauss points
        TotalTh = 0.0123653
        IF (NOEL .LE. NELEMENT) THEN
          IF (NPT .EQ. 1) THEN
            ThickCoord = -TotalTh/2.0*(0.5+0.5*Sqrt(3.0))
          END IF
          IF (NPT .EQ. 2) THEN
            ThickCoord = -TotalTh/2.0*(0.5-0.5*Sqrt(3.0))
          END IF
        ELSE IF (NOEL .GT. NELEMENT .AND. NOEL .LE. 2*NELEMENT) THEN
          IF (NPT .EQ. 1) THEN
            ThickCoord = TotalTh/2.0*(0.5-0.5*Sqrt(3.0))
          END IF
          IF (NPT .EQ. 2) THEN
            ThickCoord = TotalTh/2.0*(0.5+0.5*Sqrt(3.0))
          END IF
        END IF

        STATEV(8) = ThickCoord

C G is growth tensor
        call MutexInit( 1 )
        call MutexLock( 1 )
! TODO: please modify path
        IF (iread .EQ. 1) THEN
! --------------------------- stage 1 ---------------------------------
C read {Er,Fr,Gr,Lr,Mr,Nr} of reference surface
        open(301,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambda10St1.csv',status="old")
        read(301,*) Lambda1z0ImpSt1
        close(301)

        open(302,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambda20St1.csv',status="old")
        read(302,*) Lambda2z0ImpSt1
        close(302)

        open(303,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambdas0St1.csv',status="old")
        read(303,*) LambdaSz0ImpSt1
        close(303)

        open(304,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambda11St1.csv',status="old")
        read(304,*) Lambda1z1ImpSt1
        close(304)

        open(305,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambda21St1.csv',status="old")
        read(305,*) Lambda2z1ImpSt1
        close(305)

        open(306,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambdas1St1.csv',status="old")
        read(306,*) LambdaSz1ImpSt1
        close(306)

C read components of the two tengent vectors
        open(501,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'a1St1.csv',status="old")
        read(501,*) a1ImpSt1
        close(501)

        open(502,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'a2St1.csv',status="old")
        read(502,*) a2ImpSt1
        close(502)

        open(503,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'a3St1.csv',status="old")
        read(503,*) a3ImpSt1
        close(503)

        open(504,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'b1St1.csv',status="old")
        read(504,*) b1ImpSt1
        close(504)

        open(505,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'b2St1.csv',status="old")
        read(505,*) b2ImpSt1
        close(505)

        open(506,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'b3St1.csv',status="old")
        read(506,*) b3ImpSt1
        close(506)

! --------------------------- stage 2 ---------------------------------
C read {Er,Fr,Gr,Lr,Mr,Nr} of reference surface
        open(201,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambda10St2.csv',status="old")
        read(201,*) Lambda1z0ImpSt2
        close(201)

        open(202,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambda20St2.csv',status="old")
        read(202,*) Lambda2z0ImpSt2
        close(202)

        open(203,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambdas0St2.csv',status="old")
        read(203,*) LambdaSz0ImpSt2
        close(203)

        open(204,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambda11St2.csv',status="old")
        read(204,*) Lambda1z1ImpSt2
        close(204)

        open(205,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambda21St2.csv',status="old")
        read(205,*) Lambda2z1ImpSt2
        close(205)

        open(206,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'Lambdas1St2.csv',status="old")
        read(206,*) LambdaSz1ImpSt2
        close(206)

C read components of the two tengent vectors
        open(401,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'a1St2.csv',status="old")
        read(401,*) a1ImpSt2
        close(401)

        open(402,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'a2St2.csv',status="old")
        read(402,*) a2ImpSt2
        close(402)

        open(403,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'a3St2.csv',status="old")
        read(403,*) a3ImpSt2
        close(403)

        open(404,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'b1St2.csv',status="old")
        read(404,*) b1ImpSt2
        close(404)

        open(405,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'b2St2.csv',status="old")
        read(405,*) b2ImpSt2
        close(405)

        open(406,FILE='T:\Abaqus-Temp\20231212RobotToSophie\'//
     &  'b3St2.csv',status="old")
        read(406,*) b3ImpSt2
        close(406)

        iread=2
        END IF
        call MutexUnLock( 1 )


! G is growth tensor       
        TmpElem = NOEL
        IF (NOEL .GT. NELEMENT) THEN
          TmpElem = NOEL - NELEMENT
        END IF


! --------------------------- stage 1 ---------------------------------
C assign components of the two tengent vectors {g1,g2}
        a1 = a1ImpSt1(TmpElem)
        a2 = a2ImpSt1(TmpElem)
        a3 = a3ImpSt1(TmpElem)

        b1 = b1ImpSt1(TmpElem)
        b2 = b2ImpSt1(TmpElem)
        b3 = b3ImpSt1(TmpElem)

C assign growth functions Lambda
        Lambda1z0 = Lambda1z0ImpSt1(TmpElem)
        Lambda2z0 = Lambda2z0ImpSt1(TmpElem)
        LambdaSz0 = LambdaSz0ImpSt1(TmpElem)

        Lambda1z1 = Lambda1z1ImpSt1(TmpElem)
        Lambda2z1 = Lambda2z1ImpSt1(TmpElem)
        LambdaSz1 = LambdaSz1ImpSt1(TmpElem)

C growth functions in the curvilinear coordinates
        Lam11 = Lambda1z0 + ThickCoord*Lambda1z1
        Lam12 = Lambdasz0 + ThickCoord*Lambdasz1
        Lam21 = Lam12
        Lam22 = Lambda2z0 + ThickCoord*Lambda2z1
        Lam33 = 1.0

C FinalG is growth tensor at the end of the growth process
C tranform FinalG from curvilinear coordinates to Cartesian coordiantes
        FinalG11St1 = (a1*(b2**2+b3**2)*(a1*Lam11+b1*Lam21)-a2*b2
     &  *(a1**2*Lam12+b1**2*Lam21+a1*b1*(Lam11+Lam22))-a3*b3*(2*a2*b2
     &  +a1**2*Lam12+b1**2*Lam21+a1*b1*(Lam11+Lam22))+a3**2*(b2**2+b1
     &  *(a1*Lam12+b1*Lam22))+a2**2*(b3**2+b1*(a1*Lam12+b1*Lam22)))
     &  /(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2*(a1*b1+a3*b3)
     &  +a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG12St1 = (a1**3*b2*Lam12+a2*b1*(b1**2+b3**2)*Lam21+a1*(a2
     &  *b3**2*(-1+Lam11)-b1**2*b2*Lam21+a2*b1**2*(Lam11-Lam22))-a1**2
     &  *b1*(a2*Lam12+b2*(Lam11-Lam22))+a3**2*b2*(a1*Lam12+b1
     &  *(-1+Lam22))-a3*b3*(a1*b2*(-1+Lam11)+a1*a2*Lam12+b1*b2*Lam21
     &  +a2*b1*(-1+Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG13St1 = (a1**3*b3*Lam12+a3*b1*(b1**2+b2**2)*Lam21+a1*(a3
     &  *b2**2*(-1+Lam11)-b1**2*b3*Lam21+a3*b1**2*(Lam11-Lam22))-a1**
     &  2*b1*(a3*Lam12+b3*(Lam11-Lam22))+a2**2*b3*(a1*Lam12+b1*(-1
     &  +Lam22))-a2*b2*(a1*b3*(-1+Lam11)+a1*a3*Lam12+b1*b3*Lam21
     &  +a3*b1*(-1+Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))

        FinalG21St1 = (a1*(a2*b3**2*(-1+Lam11)-a2**2*b2*Lam12+b2*(b2**2
     &  +b3**2)*Lam21+a2*b2**2*(Lam11-Lam22))+a3**2*b1*(a2*Lam12+b2
     &  *(-1+Lam22))-a3*b3*(a2*b1*(-1+Lam11)+a1*a2*Lam12+b1*b2*Lam21
     &  +a1*b2*(-1+Lam22))+a2*b1*(a2**2*Lam12-b2**2*Lam21+a2*b2
     &  *(-Lam11+Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG22St1 = (a2*(b1**2+b3**2)*(a2*Lam11+b2*Lam21)-a1*b1*(a2
     &  **2*Lam12+b2**2*Lam21+a2*b2*(Lam11+Lam22))-a3*b3*(2*a1*b1+a2**2
     &  *Lam12+b2**2*Lam21+a2*b2*(Lam11+Lam22))+a3**2*(b1**2+b2*(a2
     &  *Lam12+b2*Lam22))+a1**2*(b3**2+b2*(a2*Lam12+b2*Lam22)))/(a3**2
     &  *(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2*(a1*b1+a3*b3)+a2**2
     &  *(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG23St1 = (a2**3*b3*Lam12+a2*(a1*b3*(b1-b1*Lam11+a1*Lam12)
     &  -b2**2*b3*Lam21+a3*(b1**2*(-1+Lam11)-a1*b1*Lam12+b2**2*(Lam11
     &  -Lam22)))-a2**2*b2*(a3*Lam12+b3*(Lam11-Lam22))+b2*(a3*(b1**2
     &  +b2**2)*Lam21+a1**2*b3*(-1+Lam22)+a1*b1*(a3-b3*Lam21
     &  -a3*Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))

        FinalG31St1 = (a1*(a3*b2**2*(-1+Lam11)-a3**2*b3*Lam12+b3*(b2**2
     &  +b3**2)*Lam21+a3*b3**2*(Lam11-Lam22))+a2**2*b1*(a3*Lam12+b3
     &  *(-1+Lam22))-a2*b2*(a3*b1*(-1+Lam11)+a1*a3*Lam12+b1*b3*Lam21
     &  +a1*b3*(-1+Lam22))+a3*b1*(a3**2*Lam12-b3**2*Lam21+a3*b3
     &  *(-Lam11+Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG32St1 = (a2*(a3*b1**2*(-1+Lam11)-a1*a3*b1*Lam12-a3**2*b3
     &  *Lam12+b3*(b1**2+b3**2)*Lam21+a3*b3**2*(Lam11-Lam22)+a1*b3
     &  *(b1-b1*Lam22))+b2*(a1*b1*(a3-a3*Lam11-b3*Lam21)+a1**2*(a3
     &  *Lam12+b3*(-1+Lam22))+a3*(a3**2*Lam12-b3**2*Lam21+a3*b3
     &  *(-Lam11+Lam22))))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG33St1 = (a3*(b1**2+b2**2)*(a3*Lam11+b3*Lam21)-a1*b1*(a3
     &  **2*Lam12+b3**2*Lam21+a3*b3*(Lam11+Lam22))-a2*b2*(2*a1*b1+a3**2
     &  *Lam12+b3**2*Lam21+a3*b3*(Lam11+Lam22))+a2**2*(b1**2+b3*(a3
     &  *Lam12+b3*Lam22))+a1**2*(b2**2+b3*(a3*Lam12+b3*Lam22)))/(a3**2
     &  *(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2*(a1*b1+a3*b3)+a2**2
     &  *(b1**2+b3**2)+a1**2*(b2**2+b3**2))

! --------------------------- stage 2 ---------------------------------
C assign components of the two tengent vectors {g1,g2}
        a1 = a1ImpSt2(TmpElem)
        a2 = a2ImpSt2(TmpElem)
        a3 = a3ImpSt2(TmpElem)

        b1 = b1ImpSt2(TmpElem)
        b2 = b2ImpSt2(TmpElem)
        b3 = b3ImpSt2(TmpElem)

C assign growth functions Lambda
        Lambda1z0 = Lambda1z0ImpSt2(TmpElem)
        Lambda2z0 = Lambda2z0ImpSt2(TmpElem)
        LambdaSz0 = LambdaSz0ImpSt2(TmpElem)

        Lambda1z1 = Lambda1z1ImpSt2(TmpElem)
        Lambda2z1 = Lambda2z1ImpSt2(TmpElem)
        LambdaSz1 = LambdaSz1ImpSt2(TmpElem)

C growth functions in the curvilinear coordinates
        Lam11 = Lambda1z0 + ThickCoord*Lambda1z1
        Lam12 = Lambdasz0 + ThickCoord*Lambdasz1
        Lam21 = Lam12
        Lam22 = Lambda2z0 + ThickCoord*Lambda2z1
        Lam33 = 1.0

C FinalG is growth tensor at the end of the growth process
C tranform FinalG from curvilinear coordinates to Cartesian coordiantes
        FinalG11St2 = (a1*(b2**2+b3**2)*(a1*Lam11+b1*Lam21)-a2*b2
     &  *(a1**2*Lam12+b1**2*Lam21+a1*b1*(Lam11+Lam22))-a3*b3*(2*a2*b2
     &  +a1**2*Lam12+b1**2*Lam21+a1*b1*(Lam11+Lam22))+a3**2*(b2**2+b1
     &  *(a1*Lam12+b1*Lam22))+a2**2*(b3**2+b1*(a1*Lam12+b1*Lam22)))
     &  /(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2*(a1*b1+a3*b3)
     &  +a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG12St2 = (a1**3*b2*Lam12+a2*b1*(b1**2+b3**2)*Lam21+a1*(a2
     &  *b3**2*(-1+Lam11)-b1**2*b2*Lam21+a2*b1**2*(Lam11-Lam22))-a1**2
     &  *b1*(a2*Lam12+b2*(Lam11-Lam22))+a3**2*b2*(a1*Lam12+b1
     &  *(-1+Lam22))-a3*b3*(a1*b2*(-1+Lam11)+a1*a2*Lam12+b1*b2*Lam21
     &  +a2*b1*(-1+Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG13St2 = (a1**3*b3*Lam12+a3*b1*(b1**2+b2**2)*Lam21+a1*(a3
     &  *b2**2*(-1+Lam11)-b1**2*b3*Lam21+a3*b1**2*(Lam11-Lam22))-a1**
     &  2*b1*(a3*Lam12+b3*(Lam11-Lam22))+a2**2*b3*(a1*Lam12+b1*(-1
     &  +Lam22))-a2*b2*(a1*b3*(-1+Lam11)+a1*a3*Lam12+b1*b3*Lam21
     &  +a3*b1*(-1+Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))

        FinalG21St2 = (a1*(a2*b3**2*(-1+Lam11)-a2**2*b2*Lam12+b2*(b2**2
     &  +b3**2)*Lam21+a2*b2**2*(Lam11-Lam22))+a3**2*b1*(a2*Lam12+b2
     &  *(-1+Lam22))-a3*b3*(a2*b1*(-1+Lam11)+a1*a2*Lam12+b1*b2*Lam21
     &  +a1*b2*(-1+Lam22))+a2*b1*(a2**2*Lam12-b2**2*Lam21+a2*b2
     &  *(-Lam11+Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG22St2 = (a2*(b1**2+b3**2)*(a2*Lam11+b2*Lam21)-a1*b1*(a2
     &  **2*Lam12+b2**2*Lam21+a2*b2*(Lam11+Lam22))-a3*b3*(2*a1*b1+a2**2
     &  *Lam12+b2**2*Lam21+a2*b2*(Lam11+Lam22))+a3**2*(b1**2+b2*(a2
     &  *Lam12+b2*Lam22))+a1**2*(b3**2+b2*(a2*Lam12+b2*Lam22)))/(a3**2
     &  *(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2*(a1*b1+a3*b3)+a2**2
     &  *(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG23St2 = (a2**3*b3*Lam12+a2*(a1*b3*(b1-b1*Lam11+a1*Lam12)
     &  -b2**2*b3*Lam21+a3*(b1**2*(-1+Lam11)-a1*b1*Lam12+b2**2*(Lam11
     &  -Lam22)))-a2**2*b2*(a3*Lam12+b3*(Lam11-Lam22))+b2*(a3*(b1**2
     &  +b2**2)*Lam21+a1**2*b3*(-1+Lam22)+a1*b1*(a3-b3*Lam21
     &  -a3*Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))

        FinalG31St2 = (a1*(a3*b2**2*(-1+Lam11)-a3**2*b3*Lam12+b3*(b2**2
     &  +b3**2)*Lam21+a3*b3**2*(Lam11-Lam22))+a2**2*b1*(a3*Lam12+b3
     &  *(-1+Lam22))-a2*b2*(a3*b1*(-1+Lam11)+a1*a3*Lam12+b1*b3*Lam21
     &  +a1*b3*(-1+Lam22))+a3*b1*(a3**2*Lam12-b3**2*Lam21+a3*b3
     &  *(-Lam11+Lam22)))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG32St2 = (a2*(a3*b1**2*(-1+Lam11)-a1*a3*b1*Lam12-a3**2*b3
     &  *Lam12+b3*(b1**2+b3**2)*Lam21+a3*b3**2*(Lam11-Lam22)+a1*b3
     &  *(b1-b1*Lam22))+b2*(a1*b1*(a3-a3*Lam11-b3*Lam21)+a1**2*(a3
     &  *Lam12+b3*(-1+Lam22))+a3*(a3**2*Lam12-b3**2*Lam21+a3*b3
     &  *(-Lam11+Lam22))))/(a3**2*(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2
     &  *(a1*b1+a3*b3)+a2**2*(b1**2+b3**2)+a1**2*(b2**2+b3**2))
        FinalG33St2 = (a3*(b1**2+b2**2)*(a3*Lam11+b3*Lam21)-a1*b1*(a3
     &  **2*Lam12+b3**2*Lam21+a3*b3*(Lam11+Lam22))-a2*b2*(2*a1*b1+a3**2
     &  *Lam12+b3**2*Lam21+a3*b3*(Lam11+Lam22))+a2**2*(b1**2+b3*(a3
     &  *Lam12+b3*Lam22))+a1**2*(b2**2+b3*(a3*Lam12+b3*Lam22)))/(a3**2
     &  *(b1**2+b2**2)-2*a1*a3*b1*b3-2*a2*b2*(a1*b1+a3*b3)+a2**2
     &  *(b1**2+b3**2)+a1**2*(b2**2+b3**2))

C G is growth tensor during the growth process
      IF ( (TIME(1)+DTIME) .LE. (TotalT/2.0)) THEN
        G11=1.0+(FinalG11St1-1.0)*2.0*(TIME(1)+DTIME)/TotalT
        G12=FinalG12St1*2.0*(TIME(1)+DTIME)/TotalT
        G13=FinalG13St1*2.0*(TIME(1)+DTIME)/TotalT
        G21=FinalG21St1*2.0*(TIME(1)+DTIME)/TotalT
        G22=1.0+(FinalG22St1-1.0)*2.0*(TIME(1)+DTIME)/TotalT
        G23=FinalG23St1*2.0*(TIME(1)+DTIME)/TotalT
        G31=FinalG31St1*2.0*(TIME(1)+DTIME)/TotalT
        G32=FinalG32St1*2.0*(TIME(1)+DTIME)/TotalT
        G33=1.0+(FinalG33St1-1.0)*2.0*(TIME(1)+DTIME)/TotalT
      ELSE
        G11=2.0*FinalG11St1-FinalG11St2
     &  +(FinalG11St2-FinalG11St1)*2.0*(TIME(1)+DTIME)/TotalT
        G12=2.0*FinalG12St1-FinalG12St2
     &  +(FinalG12St2-FinalG12St1)*2.0*(TIME(1)+DTIME)/TotalT
        G13=2.0*FinalG13St1-FinalG13St2
     &  +(FinalG13St2-FinalG13St1)*2.0*(TIME(1)+DTIME)/TotalT
        G21=2.0*FinalG21St1-FinalG21St2
     &  +(FinalG21St2-FinalG21St1)*2.0*(TIME(1)+DTIME)/TotalT
        G22=2.0*FinalG22St1-FinalG22St2
     &  +(FinalG22St2-FinalG22St1)*2.0*(TIME(1)+DTIME)/TotalT
        G23=2.0*FinalG23St1-FinalG23St2
     &  +(FinalG23St2-FinalG23St1)*2.0*(TIME(1)+DTIME)/TotalT
        G31=2.0*FinalG31St1-FinalG31St2
     &  +(FinalG31St2-FinalG31St1)*2.0*(TIME(1)+DTIME)/TotalT
        G32=2.0*FinalG32St1-FinalG32St2
     &  +(FinalG32St2-FinalG32St1)*2.0*(TIME(1)+DTIME)/TotalT
        G33=2.0*FinalG33St1-FinalG33St2
     &  +(FinalG33St2-FinalG33St1)*2.0*(TIME(1)+DTIME)/TotalT
      END IF

C Det(G)
        STATEV(9) = G11*G22*G33-G12*G21*G33+G12*G23*G31
     &      +G13*G32*G21-G13*G31*G22-G23*G32*G11

        STATEV(1) = FinalG11
        STATEV(2) = FinalG12
        STATEV(3) = FinalG13
        STATEV(4) = Lambda1z0
        STATEV(5) = Lambda1z1
        STATEV(6) = Lambda2z0
        STATEV(7) = Lambda2z1
C
C Reset components to be zero
C
        DO K1=1, 3
          DO K2=1, 3
            AA(K1, K2)=0.0
          END DO
        END DO
C
        DO K1=1, 6
            BAA(K1)=0.0
        END DO
C
        DO K1=1, 6
            BAAB(K1)=0.0
        END DO
C
C Elastic deformation tensor A=F.G**(-1)
        AA(1,1)=(DFGRD1(1,3)*G22*G31-DFGRD1(1,2)*G23*G31
     &  -DFGRD1(1,3)*G21*G32+DFGRD1(1,1)*G23*G32+DFGRD1(1,2)*G21*G33
     &  -DFGRD1(1,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        AA(1,2)=(DFGRD1(1,3)*G12*G31-DFGRD1(1,2)*G13*G31
     &  -DFGRD1(1,3)*G11*G32+DFGRD1(1,1)*G13*G32+DFGRD1(1,2)*G11*G33
     &  -DFGRD1(1,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
        AA(2,1)=(DFGRD1(2,3)*G22*G31-DFGRD1(2,2)*G23*G31
     &  -DFGRD1(2,3)*G21*G32+DFGRD1(2,1)*G23*G32+DFGRD1(2,2)*G21*G33
     &  -DFGRD1(2,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        AA(2,2)=(DFGRD1(2,3)*G12*G31-DFGRD1(2,2)*G13*G31
     &  -DFGRD1(2,3)*G11*G32+DFGRD1(2,1)*G13*G32+DFGRD1(2,2)*G11*G33
     &  -DFGRD1(2,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
        AA(3,3)=(DFGRD1(3,3)*G12*G21-DFGRD1(3,2)*G13*G21
     &  -DFGRD1(3,3)*G11*G22+DFGRD1(3,1)*G13*G22+DFGRD1(3,2)*G11*G23
     &  -DFGRD1(3,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
!       IF*(NTENS*.EQ.*3)*THEN
        AA(1,3)=(DFGRD1(1,3)*G12*G21-DFGRD1(1,2)*G13*G21
     &  -DFGRD1(1,3)*G11*G22+DFGRD1(1,1)*G13*G22+DFGRD1(1,2)*G11*G23
     &  -DFGRD1(1,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        AA(2,3)=(DFGRD1(2,3)*G12*G21-DFGRD1(2,2)*G13*G21
     &  -DFGRD1(2,3)*G11*G22+DFGRD1(2,1)*G13*G22+DFGRD1(2,2)*G11*G23
     &  -DFGRD1(2,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        AA(3,1)=(DFGRD1(3,3)*G22*G31-DFGRD1(3,2)*G23*G31
     &  -DFGRD1(3,3)*G21*G32+DFGRD1(3,1)*G23*G32+DFGRD1(3,2)*G21*G33
     &  -DFGRD1(3,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        AA(3,2)=(DFGRD1(3,3)*G12*G31-DFGRD1(3,2)*G13*G31
     &  -DFGRD1(3,3)*G11*G32+DFGRD1(3,1)*G13*G32+DFGRD1(3,2)*G11*G33
     &  -DFGRD1(3,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
        ! END IF
C
C Determinent of the elastic deformation tensor DETA=det(A)
C
        DETAA=AA(1, 1)*AA(2, 2)*AA(3, 3)
     &     -AA(1, 2)*AA(2, 1)*AA(3, 3)
C        IF(NSHR.EQ.3) THEN
        DETAA=DETAA+AA(1, 2)*AA(2, 3)*AA(3, 1)
     &     +AA(1, 3)*AA(3, 2)*AA(2, 1)
     &     -AA(1, 3)*AA(3,1)*AA(2, 2)
     &     -AA(2, 3)*AA(3, 2)*AA(1, 1)
C        END IF
C
C Left Cauchy-Green strain tensor BA=A*A**(T)
C
        BAA(1)=AA(1, 1)**2+AA(1, 2)**2
        BAA(2)=AA(2, 1)**2+AA(2, 2)**2
        BAA(3)=AA(3, 3)**2
        BAA(4)=AA(1, 1)*AA(2, 1)+AA(1, 2)*AA(2, 2)
C        IF(NSHR.EQ.3) THEN
          BAA(1)=BAA(1)+AA(1, 3)**2
          BAA(2)=BAA(2)+AA(2, 3)**2
          BAA(3)=BAA(3)+AA(3, 1)**2+AA(3, 2)**2
          BAA(4)=BAA(4)+AA(1, 3)*AA(2, 3)
          BAA(5)=AA(1, 1)*AA(3, 1)+AA(1, 2)*AA(3, 2)+AA(1, 3)*AA(3, 3)
          BAA(6)=AA(2, 1)*AA(3, 1)+AA(2, 2)*AA(3, 2)+AA(2, 3)*AA(3, 3)
C        END IF
C
C Deviatoric left Cauchy-Green strain tensor BAB=DetBA**(-2/3)*BA
        DETBAA=BAA(1)*BAA(2)*BAA(3)
     &     -BAA(4)*BAA(4)*BAA(3)
C        IF(NSHR.EQ.3) THEN
        DETBAA=DETBAA+BAA(4)*BAA(6)*BAA(5)
     &     +BAA(5)*BAA(6)*BAA(4)
     &     -BAA(5)*BAA(5)*BAA(2)
     &     -BAA(6)*BAA(6)*BAA(1)
C        END IF
        BAAB(1)=DETBAA**(-TWO/THREE)*BAA(1)
        BAAB(2)=DETBAA**(-TWO/THREE)*BAA(2)
        BAAB(3)=DETBAA**(-TWO/THREE)*BAA(3)
        BAAB(4)=DETBAA**(-TWO/THREE)*BAA(4)
C        IF(NSHR.EQ.3) THEN
        BAAB(5)=DETBAA**(-TWO/THREE)*BAA(5)
        BAAB(6)=DETBAA**(-TWO/THREE)*BAA(6)
C        END IF
C
C Calculate Cauchy stress
C
        TRBAAB=BAAB(1)+BAAB(2)+BAAB(3)
        EG=TWO*C10/DETAA
        PR=TWO/D1*(DETAA-ONE)
        DO K1=1,NDI
        STRESS(K1)=EG*(BAAB(K1)-TRBAAB/THREE)+PR
        END DO
        DO K1=NDI+1,NDI+NSHR
        STRESS(K1)=EG*BAAB(K1)
        END DO
C
C Calculate the consistent Jacobian
C
        EG23=EG*TWO/THREE
        EK=TWO/D1*(TWO*DETAA-ONE)
        DDSDDE(1, 1)= EG23*(BAAB(1)+TRBAAB/THREE)+EK
        DDSDDE(1, 2)=-EG23*(BAAB(1)+BAAB(2)-TRBAAB/THREE)+EK
        DDSDDE(1, 3)=-EG23*(BAAB(1)+BAAB(3)-TRBAAB/THREE)+EK
        DDSDDE(1, 4)= EG23*BAAB(4)/TWO
        DDSDDE(2, 2)= EG23*(BAAB(2)+TRBAAB/THREE)+EK
        DDSDDE(2, 3)=-EG23*(BAAB(2)+BAAB(3)-TRBAAB/THREE)+EK
        DDSDDE(2, 4)= EG23*BAAB(4)/TWO
        DDSDDE(3, 3)= EG23*(BAAB(3)+TRBAAB/THREE)+EK
        DDSDDE(3, 4)=-EG23*BAAB(4)
        DDSDDE(4, 4)= EG*(BAAB(1)+BAAB(2))/TWO
C        IF(NSHR.EQ.3) THEN
        DDSDDE(1, 5)= EG23*BAAB(5)/TWO
        DDSDDE(1, 6)=-EG23*BAAB(6)
        DDSDDE(2, 5)=-EG23*BAAB(5)
        DDSDDE(2, 6)= EG23*BAAB(6)/TWO
        DDSDDE(3, 5)= EG23*BAAB(5)/TWO
        DDSDDE(3, 6)= EG23*BAAB(6)/TWO
        DDSDDE(4, 5)= EG*BAAB(6)/TWO
        DDSDDE(4, 6)= EG*BAAB(5)/TWO
        DDSDDE(5, 5)= EG*(BAAB(1)+BAAB(3))/TWO
        DDSDDE(5, 6)= EG*BAAB(4)/TWO
        DDSDDE(6, 6)= EG*(BAAB(2)+BAAB(3))/TWO
C        END IF
        DO K1=1, NTENS
          DO K2=1, K1-1
            DDSDDE(K1, K2)=DDSDDE(K2, K1)
          END DO
        END DO
C
        RETURN
       END