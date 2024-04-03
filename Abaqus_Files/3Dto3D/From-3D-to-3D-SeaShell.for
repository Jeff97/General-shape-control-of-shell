C ----------------------------------------------------------------
C UMAT FOR COMPRESSIBLE NEO-HOOKEAN HYPERELASTICITY
C CANNOT BE USED FOR PLANE STRESS
C ----------------------------------------------------------------
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
C A1 - Elastic tensor at the end of the increment
C BA1 - Left Cauchy-green elastic tensor at the end of the increment
C BA1B - Deviatoric left Cauchy-green elastic tensor at the end of the increment
C TRBA1B - Trace of BA1B
C Some geometry quantities:
C       Er, Fr, Gr, Lr, Mr, Nr are related to reference surface
C       EE, FF, GG, LL, MM, NN are related to target surface
C ----------------------------------------------------------------
C
        DIMENSION A1(3,3), BA1(6), BA1B(6)
C
        PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0,
     &    SIX=6.D0)
        DOUBLE PRECISION :: Lambda1z0, Lambda1z1, Lambda2z0, Lambda2z1,
     &  LambdaSz0, LambdaSz1, Pi, TotalT, Delta, DeltaR, X, Y, Z, 
     &  Er, Fr, Gr, Lr, Mr, Nr, EE, FF, GG, LL, MM, NN,
     &  FinalG11, FinalG12, FinalG13, 
     &  FinalG21, FinalG22, FinalG23, 
     &  FinalG31, FinalG32, FinalG33, 
     &  LamrR, LamrT, LamrZ, 
     &  LamtR, LamtT, LamtZ,
     &  LamzR, LamzT, LamzZ,
     &  G11, G12, G13, G21, G22, G23, G31, G32, G33, 
     &  MMOD, Theta, SQMMOD, RadCylinder, Theta1, Theta2,
     &  FinalG11St1, FinalG12St1, FinalG13St1, 
     &  FinalG21St1, FinalG22St1, FinalG23St1, 
     &  FinalG31St1, FinalG32St1, FinalG33St1, 
     &  FinalG11St2, FinalG12St2, FinalG13St2, 
     &  FinalG21St2, FinalG22St2, FinalG23St2, 
     &  FinalG31St2, FinalG32St2, FinalG33St2

C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C STATEV(1) - Initial coordinate COORDS(1)
C STATEV(2) - Initial coordinate COORDS(2)
C STATEV(3) - Initial coordinate COORDS(3)
C STATEV(4) - 
C STATEV(5) - 
C STATEV(9) - Norm of growth tensor
C ----------------------------------------------------------------
C ELASTIC PROPERTIES
        EMOD=PROPS(1)
        ENU=PROPS(2)
        C10=EMOD/(FOUR*(ONE+ENU))
        D1=SIX*(ONE-TWO*ENU)/EMOD
        Pi=3.14159265359
C get coordinate from COORDS variable
        IF (STATEV(1) .EQ. 0) THEN
                STATEV(1)=COORDS(1)
        END IF
        IF (STATEV(2) .EQ. 0) THEN
                STATEV(2)=COORDS(2)
        END IF
        IF (STATEV(3) .EQ. 0) THEN
                STATEV(3)=COORDS(3)
        END IF

C store coordinate from STATEV
        X = STATEV(1)
        Y = STATEV(2)
        Z = STATEV(3)

C Calculate R and Theta
        MMOD = X*X+Y*Y
        SQMMOD = Sqrt(MMOD)
        IF (X .GE. 0 .AND. Y .GE. 0) THEN
            Theta = ASin((Y/SQMMOD))
        END IF
        IF (X .LE. 0 .AND. Y .GE. 0) THEN
            Theta = ASin((-X/SQMMOD))+0.5*Pi
        END IF
        IF (X .LE. 0 .AND. Y .LE. 0) THEN
            Theta = ASin((-Y/SQMMOD))+Pi
        END IF
        IF (X .GE. 0 .AND. Y .LE. 0) THEN
            Theta = ASin((X/SQMMOD))+1.5*Pi
        END IF
        STATEV(4) = Theta
        Theta1 = Theta
        Theta2 = Z

C the radius of cylinder in the reference configuration
        RadCylinder = 4.0

        ThickCoord = SQMMOD-RadCylinder

!I cannot get total time from UMAT variable, please modify TotalT manually 
        TotalT=10.0

C assign growth functions Lambda
        Er = 16.0
        Fr = 0.0
        Gr = 1.0
        Lr = -4.0
        Mr = 0.0
        Nr = 0.0

! --------------------------- stage 1 ---------------------------------

! half sea shell
        EE = 9*Sin((1+Theta2)/8.0)**2
        FF = -(3.0/4.0)*(6*Cos(Theta1)*Sin((1+Theta2)/8.0)**2
     &  +Sin(Theta1)*Sin((1+Theta2)/4.0))
        GG = 1.0/256*((4+3*Cos(Theta1))**2*(7*Cos((7*(1+Theta2))
     &  /8.0)-9*Cos((9*(1+Theta2))/8.0))**2+36*(Cos((1+Theta2)
     &  /8.0)*Sin(Theta1)-4*Sin((1+Theta2)/8.0))**2+(4
     &  +3*Cos(Theta1))**2*(7*Sin((7*(1+Theta2))/8.0)-9
     &  *Sin((9*(1+Theta2))/8.0))**2)
        LL = (24*Sqrt(2.0)*(4+3*Cos(Theta1))*Sin((1+Theta2)
     &  /8.0)**3)/(Sqrt((1401+1560*Cos(Theta1)+224*Cos(2*Theta1)
     &  -1367*Cos((1+Theta2)/4.0)-128*Cos(1.0/4*(1-8*Theta1+Theta2))
     &  -792*Cos(1.0/4*(1-4*Theta1+Theta2))-720*Cos(1.0/4
     &  *(1+4*Theta1+Theta2))-80*Cos(1.0/4*(1+8*Theta1+Theta2)))
     &  *Sin((1+Theta2)/8.0)**2))
        MM = -((3*Sqrt(2.0)*(3*Cos((1+Theta2)/8.0)-4*Cos(1.0/8
     &  *(1-8*Theta1+Theta2))+8*Cos(1.0/8*(1+8*Theta1+Theta2)))
     &  *Sin(Theta1)*Sin((1+Theta2)/8.0)**2)/(Sqrt((1401+1560
     &  *Cos(Theta1)+224*Cos(2*Theta1)-1367*Cos((1+Theta2)/4.0)
     &  -128*Cos(1.0/4*(1-8*Theta1+Theta2))-792*Cos(1.0/4*(1
     &  -4*Theta1+Theta2))-720*Cos(1.0/4*(1+4*Theta1+Theta2))
     &  -80*Cos(1.0/4*(1+8*Theta1+Theta2)))*Sin((1+Theta2)/8.0)**2)))
        NN = ((4+3*Cos(Theta1))*Sin((1+Theta2)/8.0)*(105+268
     &  *Cos(Theta1)+96*Cos(2*Theta1)-3*(31+84*Cos(Theta1)
     &  +32*Cos(2*Theta1))*Cos((1+Theta2)/4.0)-12*Sin(Theta1)
     &  *Sin((1+Theta2)/4.0)))/(8*Sqrt(2.0)*Sqrt((1401+1560
     &  *Cos(Theta1)+224*Cos(2*Theta1)-1367*Cos((1+Theta2)/4.0)
     &  -128*Cos(1.0/4*(1-8*Theta1+Theta2))-792*Cos(1.0/4
     &  *(1-4*Theta1+Theta2))-720*Cos(1.0/4*(1+4*Theta1+Theta2))
     &  -80*Cos(1.0/4*(1+8*Theta1+Theta2)))*Sin((1+Theta2)/8.0)**2))


        Delta = EE*GG - FF**TWO
        DeltaR = Er*Gr - Fr**TWO

        Lambda1z0 = (EE*Er+Sqrt(Er*Gr*Delta))/
     &              (Er*Sqrt(EE*Er+GG*Gr+TWO*Sqrt(Er*Gr*Delta)))
        Lambda2z0 = (GG*Gr+Sqrt(Er*Gr*Delta))/
     &              (Gr*Sqrt(EE*Er+GG*Gr+TWO*Sqrt(Er*Gr*Delta)))
        LambdaSz0 = FF/Sqrt(EE*Er+GG*Gr+TWO*Sqrt(Er*Gr*Delta))

        Lambda1z1 = (Lambda1z0*(((GG*Gr*Lr-Fr*GG*Mr-FF*Gr*Mr+FF*Fr*Nr)
     &  *Delta-GG*(GG*LL-TWO*FF*MM+EE*NN)*DeltaR+NN*Delta*DeltaR)
     &  *Lambda1z0-Delta*(FF*Fr*Lr-EE*Gr*Lr-Er*FF*Mr+EE*Fr*Mr
     &  +LL*DeltaR)*Lambda2z0)+(-Fr*(GG*Lr-TWO*FF*Mr+EE*Nr)
     &  *Delta+(Er*GG*Mr+EE*Gr*Mr-FF*(Gr*Lr+Er*Nr))*Delta+TWO
     &  *(FF*GG*LL-TWO*EE*GG*MM+EE*FF*NN+MM*Delta)*DeltaR)
     &  *Lambda1z0*LambdaSz0+((TWO*FF*Fr*Lr-EE*Gr*Lr-TWO*Er*FF*Mr
     &  +EE*Er*Nr)*Delta-EE*(GG*LL-TWO*FF*MM+EE*NN)*DeltaR+TWO
     &  *LL*Delta*DeltaR)*LambdaSz0**TWO)/
     &  (Delta*DeltaR*(GG*Lambda1z0+EE*Lambda2z0-TWO*FF*LambdaSz0))

        Lambda2z1 = (Lambda2z0*(-Delta*(Fr*GG*Mr-FF*Gr*Mr+FF*Fr*Nr
     &  -Er*GG*Nr+NN*DeltaR)*Lambda1z0+((FF*Fr*Lr-Er*FF*Mr-EE*Fr
     &  *Mr+EE*Er*Nr)*Delta-EE*(GG*LL-TWO*FF*MM+EE*NN)*DeltaR+LL
     &  *Delta*DeltaR)*Lambda2z0)+(-Fr*(GG*Lr-TWO*FF*Mr+EE*Nr)
     &  *Delta+(Er*GG*Mr+EE*Gr*Mr-FF*(Gr*Lr+Er*Nr))*Delta+TWO
     &  *(FF*GG*LL-TWO*EE*GG*MM+EE*FF*NN+MM*Delta)*DeltaR)
     &  *Lambda2z0*LambdaSz0-((-GG*Gr*Lr+TWO*FF*Gr*Mr-TWO*FF*Fr*Nr
     &  +Er*GG*Nr)*Delta+GG*(GG*LL-TWO*FF*MM+EE*NN)
     &  *DeltaR-TWO*NN*Delta*DeltaR)*LambdaSz0**TWO)
     &  /(Delta*DeltaR*(GG*Lambda1z0+EE*Lambda2z0-TWO*FF*LambdaSz0))

        LambdaSz1 = (Er*GG*Mr*Delta*Lambda1z0*Lambda2z0+EE*Gr*Mr*Delta
     &  *Lambda1z0*Lambda2z0+GG*Gr*Lr*Delta*Lambda1z0*LambdaSz0
     &  -GG**TWO*LL*DeltaR*Lambda1z0*LambdaSz0+EE*Er*Nr*Delta
     &  *Lambda2z0*LambdaSz0-EE**TWO*NN*DeltaR*Lambda2z0*LambdaSz0
     &  -FF*Gr*Delta*LambdaSz0*(Mr*Lambda1z0+Lr*LambdaSz0)
     &  +Delta*DeltaR*LambdaSz0*(NN*Lambda1z0+LL*Lambda2z0
     &  +TWO*MM*LambdaSz0)-EE*Fr*Delta*Lambda2z0*(Nr*Lambda1z0
     &  +Mr*LambdaSz0)-Fr*GG*Delta*Lambda1z0*(Lr*Lambda2z0
     &  +Mr*LambdaSz0)+FF*Fr*Delta*LambdaSz0*(Nr*Lambda1z0
     &  +Lr*Lambda2z0+TWO*Mr*LambdaSz0)-Er*FF*Delta*LambdaSz0
     &  *(Mr*Lambda2z0+Nr*LambdaSz0)+FF*GG*DeltaR*(TWO*MM*Lambda1z0
     &  *LambdaSz0+LL*(Lambda1z0*Lambda2z0+LambdaSz0**TWO))
     &  -EE*GG*DeltaR*((NN*Lambda1z0+LL*Lambda2z0)*LambdaSz0+TWO
     &  *MM*(Lambda1z0*Lambda2z0+LambdaSz0**TWO))+EE*FF*DeltaR
     &  *(TWO*MM*Lambda2z0*LambdaSz0
     &  +NN*(Lambda1z0*Lambda2z0+LambdaSz0**TWO)))
     &  /(Delta*DeltaR*(GG*Lambda1z0+EE*Lambda2z0-TWO*FF*LambdaSz0))


        LamrR = 1.0
        LamtT = Lambda1z0 + ThickCoord*Lambda1z1
        LamtZ = (Lambdasz0 + ThickCoord*Lambdasz1)*RadCylinder
        LamzT = (Lambdasz0 + ThickCoord*Lambdasz1)/RadCylinder
        LamzZ = Lambda2z0 + ThickCoord*Lambda2z1

C FinalG is growth tensor at the end of the growth process
C tranform FinalG from curvilinear coordinates to Cartesian coordiantes
        FinalG11St1 = (LamrR*X*X+LamtT*Y*Y)/MMOD
        FinalG22St1 = (LamrR*Y*Y+LamtT*X*X)/MMOD
        FinalG12St1 = (LamrR-LamtT)*X*Y/MMOD
        FinalG21St1 = FinalG12St1

        FinalG31St1 = -(LamzT*Y)/SQMMOD
        FinalG32St1 = (LamzT*X)/SQMMOD

        FinalG13St1 = -(LamtZ*Y)/SQMMOD
        FinalG23St1 = (LamtZ*X)/SQMMOD

        FinalG33St1 = LamzZ
! --------------------------- stage 2 ---------------------------------

! Full sea shell
        EE = 4*Sin((1+Theta2)/4.0)**2
        FF = 1.0/4*(-5*Cos(Theta1)+Cos(1.0/2*(1-2*Theta1+Theta2)))
     &  +Cos(1.0/2*(1+2*Theta1+Theta2))
        GG = 1.0/16*(371+390*Cos(Theta1)+64*Cos(2*Theta1)-2
     &  *(179+189*Cos(Theta1)+32*Cos(2*Theta1))*Cos((1+Theta2)/2.0)
     &  -10*Sin(Theta1)*Sin((1+Theta2)/2.0))
        LL = (16*Sqrt(2.0)*(3+2*Cos(Theta1))*Sin((1+Theta2)/4.0)**3)
     &  /(Sqrt((725+780*Cos(Theta1)+120*Cos(2*Theta1)-708
     &  *Cos((1+Theta2)/2.0)-63*Cos(1.0/2*(1-4*Theta1+Theta2))
     &  -388*Cos(1.0/2*(1-2*Theta1+Theta2))-368*Cos(1.0/2*(1
     &  +2*Theta1+Theta2))-48*Cos(1.0/2*(1+4*Theta1+Theta2)))
     &  *Sin((1+Theta2)/4.0)**2))
        MM = (4*Sqrt(2.0)*(-2*Cos((1+Theta2)/4.0)+Cos(1.0/4
     &  *(1-4*Theta1+Theta2))-4*Cos(1.0/4*(1+4*Theta1+Theta2)))
     &  *Sin(Theta1)*Sin((1+Theta2)/4.0)**2)/(Sqrt((725+780
     &  *Cos(Theta1)+120*Cos(2*Theta1)-708*Cos((1+Theta2)/2.0)
     &  -63*Cos(1.0/2*(1-4*Theta1+Theta2))-388*Cos(1.0/2*(1-2
     &  *Theta1+Theta2))-368*Cos(1.0/2*(1+2*Theta1+Theta2))
     &  -48*Cos(1.0/2*(1+4*Theta1+Theta2)))*Sin((1+Theta2)/4.0)**2))
        NN = ((3+2*Cos(Theta1))*Sin((1+Theta2)/4.0)*(70+201
     &  *Cos(Theta1)+64*Cos(2*Theta1)-(62+189*Cos(Theta1)
     &  +64*Cos(2*Theta1))*Cos((1+Theta2)/2.0)-5*Sin(Theta1)
     &  *Sin((1+Theta2)/2.0)))/(2*Sqrt(2.0)*Sqrt((725+780
     &  *Cos(Theta1)+120*Cos(2*Theta1)-708*Cos((1+Theta2)/2.0)
     &  -63*Cos(1.0/2*(1-4*Theta1+Theta2))-388*Cos(1.0/2
     &  *(1-2*Theta1+Theta2))-368*Cos(1.0/2*(1+2*Theta1+Theta2))
     &  -48*Cos(1.0/2*(1+4*Theta1+Theta2)))*Sin((1+Theta2)/4.0)**2))


        Delta = EE*GG - FF**TWO
        DeltaR = Er*Gr - Fr**TWO

        Lambda1z0 = (EE*Er+Sqrt(Er*Gr*Delta))/
     &              (Er*Sqrt(EE*Er+GG*Gr+TWO*Sqrt(Er*Gr*Delta)))
        Lambda2z0 = (GG*Gr+Sqrt(Er*Gr*Delta))/
     &              (Gr*Sqrt(EE*Er+GG*Gr+TWO*Sqrt(Er*Gr*Delta)))
        LambdaSz0 = FF/Sqrt(EE*Er+GG*Gr+TWO*Sqrt(Er*Gr*Delta))

        Lambda1z1 = (Lambda1z0*(((GG*Gr*Lr-Fr*GG*Mr-FF*Gr*Mr+FF*Fr*Nr)
     &  *Delta-GG*(GG*LL-TWO*FF*MM+EE*NN)*DeltaR+NN*Delta*DeltaR)
     &  *Lambda1z0-Delta*(FF*Fr*Lr-EE*Gr*Lr-Er*FF*Mr+EE*Fr*Mr
     &  +LL*DeltaR)*Lambda2z0)+(-Fr*(GG*Lr-TWO*FF*Mr+EE*Nr)
     &  *Delta+(Er*GG*Mr+EE*Gr*Mr-FF*(Gr*Lr+Er*Nr))*Delta+TWO
     &  *(FF*GG*LL-TWO*EE*GG*MM+EE*FF*NN+MM*Delta)*DeltaR)
     &  *Lambda1z0*LambdaSz0+((TWO*FF*Fr*Lr-EE*Gr*Lr-TWO*Er*FF*Mr
     &  +EE*Er*Nr)*Delta-EE*(GG*LL-TWO*FF*MM+EE*NN)*DeltaR+TWO
     &  *LL*Delta*DeltaR)*LambdaSz0**TWO)/
     &  (Delta*DeltaR*(GG*Lambda1z0+EE*Lambda2z0-TWO*FF*LambdaSz0))

        Lambda2z1 = (Lambda2z0*(-Delta*(Fr*GG*Mr-FF*Gr*Mr+FF*Fr*Nr
     &  -Er*GG*Nr+NN*DeltaR)*Lambda1z0+((FF*Fr*Lr-Er*FF*Mr-EE*Fr
     &  *Mr+EE*Er*Nr)*Delta-EE*(GG*LL-TWO*FF*MM+EE*NN)*DeltaR+LL
     &  *Delta*DeltaR)*Lambda2z0)+(-Fr*(GG*Lr-TWO*FF*Mr+EE*Nr)
     &  *Delta+(Er*GG*Mr+EE*Gr*Mr-FF*(Gr*Lr+Er*Nr))*Delta+TWO
     &  *(FF*GG*LL-TWO*EE*GG*MM+EE*FF*NN+MM*Delta)*DeltaR)
     &  *Lambda2z0*LambdaSz0-((-GG*Gr*Lr+TWO*FF*Gr*Mr-TWO*FF*Fr*Nr
     &  +Er*GG*Nr)*Delta+GG*(GG*LL-TWO*FF*MM+EE*NN)
     &  *DeltaR-TWO*NN*Delta*DeltaR)*LambdaSz0**TWO)
     &  /(Delta*DeltaR*(GG*Lambda1z0+EE*Lambda2z0-TWO*FF*LambdaSz0))

        LambdaSz1 = (Er*GG*Mr*Delta*Lambda1z0*Lambda2z0+EE*Gr*Mr*Delta
     &  *Lambda1z0*Lambda2z0+GG*Gr*Lr*Delta*Lambda1z0*LambdaSz0
     &  -GG**TWO*LL*DeltaR*Lambda1z0*LambdaSz0+EE*Er*Nr*Delta
     &  *Lambda2z0*LambdaSz0-EE**TWO*NN*DeltaR*Lambda2z0*LambdaSz0
     &  -FF*Gr*Delta*LambdaSz0*(Mr*Lambda1z0+Lr*LambdaSz0)
     &  +Delta*DeltaR*LambdaSz0*(NN*Lambda1z0+LL*Lambda2z0
     &  +TWO*MM*LambdaSz0)-EE*Fr*Delta*Lambda2z0*(Nr*Lambda1z0
     &  +Mr*LambdaSz0)-Fr*GG*Delta*Lambda1z0*(Lr*Lambda2z0
     &  +Mr*LambdaSz0)+FF*Fr*Delta*LambdaSz0*(Nr*Lambda1z0
     &  +Lr*Lambda2z0+TWO*Mr*LambdaSz0)-Er*FF*Delta*LambdaSz0
     &  *(Mr*Lambda2z0+Nr*LambdaSz0)+FF*GG*DeltaR*(TWO*MM*Lambda1z0
     &  *LambdaSz0+LL*(Lambda1z0*Lambda2z0+LambdaSz0**TWO))
     &  -EE*GG*DeltaR*((NN*Lambda1z0+LL*Lambda2z0)*LambdaSz0+TWO
     &  *MM*(Lambda1z0*Lambda2z0+LambdaSz0**TWO))+EE*FF*DeltaR
     &  *(TWO*MM*Lambda2z0*LambdaSz0
     &  +NN*(Lambda1z0*Lambda2z0+LambdaSz0**TWO)))
     &  /(Delta*DeltaR*(GG*Lambda1z0+EE*Lambda2z0-TWO*FF*LambdaSz0))


        LamrR = 1.0
        LamtT = Lambda1z0 + ThickCoord*Lambda1z1
        LamtZ = (Lambdasz0 + ThickCoord*Lambdasz1)*RadCylinder
        LamzT = (Lambdasz0 + ThickCoord*Lambdasz1)/RadCylinder
        LamzZ = Lambda2z0 + ThickCoord*Lambda2z1

C FinalG is growth tensor at the end of the growth process
C tranform FinalG from curvilinear coordinates to Cartesian coordiantes
        FinalG11St2 = (LamrR*X*X+LamtT*Y*Y)/MMOD
        FinalG22St2 = (LamrR*Y*Y+LamtT*X*X)/MMOD
        FinalG12St2 = (LamrR-LamtT)*X*Y/MMOD
        FinalG21St2 = FinalG12St2

        FinalG31St2 = -(LamzT*Y)/SQMMOD
        FinalG32St2 = (LamzT*X)/SQMMOD

        FinalG13St2 = -(LamtZ*Y)/SQMMOD
        FinalG23St2 = (LamtZ*X)/SQMMOD

        FinalG33St2 = LamzZ
C ---------------------------------------------------------------------

C G is growth tensor during the growth process
      G33 = 1.0
      IF ( (TIME(1)+DTIME) .LE. 5.0) THEN
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

C Reset components to be zero
        DO K1=1, 3
          DO K2=1, 3
            A1(K1, K2)=0.0
          END DO
        END DO

        DO K1=1, 6
            BA1(K1)=0.0
        END DO

        DO K1=1, 6
            BA1B(K1)=0.0
        END DO

C Elastic deformation tensor A=F.G**(-1)
        A1(1,1)=(DFGRD1(1,3)*G22*G31-DFGRD1(1,2)*G23*G31
     &  -DFGRD1(1,3)*G21*G32+DFGRD1(1,1)*G23*G32+DFGRD1(1,2)*G21*G33
     &  -DFGRD1(1,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(1,2)=(DFGRD1(1,3)*G12*G31-DFGRD1(1,2)*G13*G31
     &  -DFGRD1(1,3)*G11*G32+DFGRD1(1,1)*G13*G32+DFGRD1(1,2)*G11*G33
     &  -DFGRD1(1,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
        A1(2,1)=(DFGRD1(2,3)*G22*G31-DFGRD1(2,2)*G23*G31
     &  -DFGRD1(2,3)*G21*G32+DFGRD1(2,1)*G23*G32+DFGRD1(2,2)*G21*G33
     &  -DFGRD1(2,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(2,2)=(DFGRD1(2,3)*G12*G31-DFGRD1(2,2)*G13*G31
     &  -DFGRD1(2,3)*G11*G32+DFGRD1(2,1)*G13*G32+DFGRD1(2,2)*G11*G33
     &  -DFGRD1(2,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
        A1(3,3)=(DFGRD1(3,3)*G12*G21-DFGRD1(3,2)*G13*G21
     &  -DFGRD1(3,3)*G11*G22+DFGRD1(3,1)*G13*G22+DFGRD1(3,2)*G11*G23
     &  -DFGRD1(3,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
!       IF*(NTENS*.EQ.*3)*THEN
        A1(1,3)=(DFGRD1(1,3)*G12*G21-DFGRD1(1,2)*G13*G21
     &  -DFGRD1(1,3)*G11*G22+DFGRD1(1,1)*G13*G22+DFGRD1(1,2)*G11*G23
     &  -DFGRD1(1,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(2,3)=(DFGRD1(2,3)*G12*G21-DFGRD1(2,2)*G13*G21
     &  -DFGRD1(2,3)*G11*G22+DFGRD1(2,1)*G13*G22+DFGRD1(2,2)*G11*G23
     &  -DFGRD1(2,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(3,1)=(DFGRD1(3,3)*G22*G31-DFGRD1(3,2)*G23*G31
     &  -DFGRD1(3,3)*G21*G32+DFGRD1(3,1)*G23*G32+DFGRD1(3,2)*G21*G33
     &  -DFGRD1(3,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(3,2)=(DFGRD1(3,3)*G12*G31-DFGRD1(3,2)*G13*G31
     &  -DFGRD1(3,3)*G11*G32+DFGRD1(3,1)*G13*G32+DFGRD1(3,2)*G11*G33
     &  -DFGRD1(3,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
        ! END IF
C
C Determinent of the elastic deformation tensor DETA=det(A)
C
        DETA1=A1(1, 1)*A1(2, 2)*A1(3, 3)
     &     -A1(1, 2)*A1(2, 1)*A1(3, 3)
C        IF(NSHR.EQ.3) THEN
        DETA1=DETA1+A1(1, 2)*A1(2, 3)*A1(3, 1)
     &     +A1(1, 3)*A1(3, 2)*A1(2, 1)
     &     -A1(1, 3)*A1(3,1)*A1(2, 2)
     &     -A1(2, 3)*A1(3, 2)*A1(1, 1)
C        END IF
C
C Left Cauchy-Green strain tensor BA=A*A**(T)
C
        BA1(1)=A1(1, 1)**2+A1(1, 2)**2
        BA1(2)=A1(2, 1)**2+A1(2, 2)**2
        BA1(3)=A1(3, 3)**2
        BA1(4)=A1(1, 1)*A1(2, 1)+A1(1, 2)*A1(2, 2)
C        IF(NSHR.EQ.3) THEN
          BA1(1)=BA1(1)+A1(1, 3)**2
          BA1(2)=BA1(2)+A1(2, 3)**2
          BA1(3)=BA1(3)+A1(3, 1)**2+A1(3, 2)**2
          BA1(4)=BA1(4)+A1(1, 3)*A1(2, 3)
          BA1(5)=A1(1, 1)*A1(3, 1)+A1(1, 2)*A1(3, 2)+A1(1, 3)*A1(3, 3)
          BA1(6)=A1(2, 1)*A1(3, 1)+A1(2, 2)*A1(3, 2)+A1(2, 3)*A1(3, 3)
C        END IF
C
C Deviatoric left Cauchy-Green strain tensor BAB=DetBA**(-2/3)*BA
        DETBA1=BA1(1)*BA1(2)*BA1(3)
     &     -BA1(4)*BA1(4)*BA1(3)
C        IF(NSHR.EQ.3) THEN
        DETBA1=DETBA1+BA1(4)*BA1(6)*BA1(5)
     &     +BA1(5)*BA1(6)*BA1(4)
     &     -BA1(5)*BA1(5)*BA1(2)
     &     -BA1(6)*BA1(6)*BA1(1)
C        END IF
        BA1B(1)=DETBA1**(-TWO/THREE)*BA1(1)
        BA1B(2)=DETBA1**(-TWO/THREE)*BA1(2)
        BA1B(3)=DETBA1**(-TWO/THREE)*BA1(3)
        BA1B(4)=DETBA1**(-TWO/THREE)*BA1(4)
C        IF(NSHR.EQ.3) THEN
        BA1B(5)=DETBA1**(-TWO/THREE)*BA1(5)
        BA1B(6)=DETBA1**(-TWO/THREE)*BA1(6)
C        END IF
C
C Calculate Cauchy stress
C
        TRBA1B=BA1B(1)+BA1B(2)+BA1B(3)
        EG=TWO*C10/DETA1
        PR=TWO/D1*(DETA1-ONE)
        DO K1=1,NDI
        STRESS(K1)=EG*(BA1B(K1)-TRBA1B/THREE)+PR
        END DO
        DO K1=NDI+1,NDI+NSHR
        STRESS(K1)=EG*BA1B(K1)
        END DO
C
C Calculate the consistent Jacobian
C
        EG23=EG*TWO/THREE
        EK=TWO/D1*(TWO*DETA1-ONE)
        DDSDDE(1, 1)= EG23*(BA1B(1)+TRBA1B/THREE)+EK
        DDSDDE(1, 2)=-EG23*(BA1B(1)+BA1B(2)-TRBA1B/THREE)+EK
        DDSDDE(1, 3)=-EG23*(BA1B(1)+BA1B(3)-TRBA1B/THREE)+EK
        DDSDDE(1, 4)= EG23*BA1B(4)/TWO
        DDSDDE(2, 2)= EG23*(BA1B(2)+TRBA1B/THREE)+EK
        DDSDDE(2, 3)=-EG23*(BA1B(2)+BA1B(3)-TRBA1B/THREE)+EK
        DDSDDE(2, 4)= EG23*BA1B(4)/TWO
        DDSDDE(3, 3)= EG23*(BA1B(3)+TRBA1B/THREE)+EK
        DDSDDE(3, 4)=-EG23*BA1B(4)
        DDSDDE(4, 4)= EG*(BA1B(1)+BA1B(2))/TWO
C        IF(NSHR.EQ.3) THEN
        DDSDDE(1, 5)= EG23*BA1B(5)/TWO
        DDSDDE(1, 6)=-EG23*BA1B(6)
        DDSDDE(2, 5)=-EG23*BA1B(5)
        DDSDDE(2, 6)= EG23*BA1B(6)/TWO
        DDSDDE(3, 5)= EG23*BA1B(5)/TWO
        DDSDDE(3, 6)= EG23*BA1B(6)/TWO
        DDSDDE(4, 5)= EG*BA1B(6)/TWO
        DDSDDE(4, 6)= EG*BA1B(5)/TWO
        DDSDDE(5, 5)= EG*(BA1B(1)+BA1B(3))/TWO
        DDSDDE(5, 6)= EG*BA1B(4)/TWO
        DDSDDE(6, 6)= EG*(BA1B(2)+BA1B(3))/TWO
C        END IF
        DO K1=1, NTENS
          DO K2=1, K1-1
            DDSDDE(K1, K2)=DDSDDE(K2, K1)
          END DO
        END DO
C
        RETURN
       END