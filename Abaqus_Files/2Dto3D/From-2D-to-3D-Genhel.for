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
     &  G11St1, G12St1, G21St1, G22St1,
     &  G11St2, G12St2, G21St2, G22St2
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

!I cannot get total time from UMAT variable, please modify TotalT manually 
        TotalT=10.0

! --------------------------- stage 1 ---------------------------------

C assign growth functions Lambda
        Er = 1.0
        Fr = 0.0
        Gr = 1.0
        Lr = 0.0
        Mr = 0.0
        Nr = 0.0

        EE = 1.0/16*(1+4*Pi**2*(1-2*Y)**2)
        FF = 1.0/4*Pi*(-1+2*Y)
        GG = 1+Pi**2*(1-2*Y)**2
        LL = -((Pi**4*(-1+2*Y)**3)/Sqrt((1+2*Pi**2*(1-2*Y)**2)**2))
        MM = Pi/Sqrt((1+2*Pi**2*(1-2*Y)**2)**2)
        NN = (2*Pi**2*(2-4*Y))/Sqrt((1+2*Pi**2*(1-2*Y)**2)**2)

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

        G11St1 = Lambda1z0 + STATEV(3)*Lambda1z1
        G12St1 = LambdaSz0 + STATEV(3)*LambdaSz1
        G21St1 = G12St1
        G22St1 = Lambda2z0 + STATEV(3)*Lambda2z1

! --------------------------- stage 2 ---------------------------------

C assign growth functions Lambda
        Er = 1.0
        Fr = 0.0
        Gr = 1.0
        Lr = 0.0
        Mr = 0.0
        Nr = 0.0

        EE = 1.0/4+Pi**2*(1-2*Y)**2
        FF = Pi*(-(1.0/2)+Y)
        GG = 1+Pi**2*(1-2*Y)**2
        LL = -((4*Pi**4*(-1+2*Y)**3)/Sqrt((1+2*Pi**2*(1-2*Y)**2)**2))
        MM = (2*Pi)/Sqrt((1+2*Pi**2*(1-2*Y)**2)**2)
        NN = (2*Pi**2*(2-4*Y))/Sqrt((1+2*Pi**2*(1-2*Y)**2)**2)

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

        G11St2 = Lambda1z0 + STATEV(3)*Lambda1z1
        G12St2 = LambdaSz0 + STATEV(3)*LambdaSz1
        G21St2 = G12St2
        G22St2 = Lambda2z0 + STATEV(3)*Lambda2z1

C G is growth tensor
      G33 = 1.0
      IF ( (TIME(1)+DTIME) .LE. 5.0) THEN
        G11=1.0+(G11St1-1.0)*2.0*(TIME(1)+DTIME)/TotalT
        G22=1.0+(G22St1-1.0)*2.0*(TIME(1)+DTIME)/TotalT
        G12=G12St1*2.0*(TIME(1)+DTIME)/TotalT
        G21=G12
      ELSE
        G11=2.0*G11St1-G11St2
     &  +(G11St2-G11St1)*2.0*(TIME(1)+DTIME)/TotalT
        G22=2.0*G22St1-G22St2
     &  +(G22St2-G22St1)*2.0*(TIME(1)+DTIME)/TotalT
        G12=2.0*G12St1-G12St2
     &  +(G12St2-G12St1)*2.0*(TIME(1)+DTIME)/TotalT
        G21=G12
      END IF

C Det(G)
        STATEV(9) = G11*G22*G33-G12*G21*G33+G12*G23*G31
     &      +G13*G32*G21-G13*G31*G22-G23*G32*G11
C
! C Reset components to be zero
C
        DO K1=1, 3
          DO K2=1, 3
            A1(K1, K2)=0.0
          END DO
        END DO
C
        DO K1=1, 6
            BA1(K1)=0.0
        END DO
C
        DO K1=1, 6
            BA1B(K1)=0.0
        END DO
C
C Elastic deformation tensor A=F.G**(-1)
C
        A1(1, 1)=(DFGRD1(1, 2)*G21-DFGRD1(1, 1)*G22)
     & /(G12*G21-G11*G22)
        A1(1, 2)=(DFGRD1(1, 2)*G11-DFGRD1(1, 1)*G12)
     & /(-G12*G21+G11*G22)
        A1(2, 1)=(DFGRD1(2, 2)*G21-DFGRD1(2, 1)*G22)
     & /(G12*G21-G11*G22)
        A1(2, 2)=(DFGRD1(2, 2)*G11-DFGRD1(2, 1)*G12)
     & /(-G12*G21+G11*G22)
        A1(3, 3)=DFGRD1(3, 3)/G33
C        IF(NSHR.EQ.3) THEN
        A1(1, 3)=DFGRD1(1, 3)/G33
        A1(2, 3)=DFGRD1(2, 3)/G33
        A1(3, 1)=(DFGRD1(3, 2)*G21-DFGRD1(3, 1)*G22)
     & /(G12*G21-G11*G22)
        A1(3, 2)=(DFGRD1(3, 2)*G11-DFGRD1(3, 1)*G12)
     & /(-G12*G21+G11*G22)
C        END IF
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