      SUBROUTINE BLDPMAT
C
C     ------------------------------------------------------------------
C
C     DRIVER ROUTINE TO GENERATE DENSITY MATRIX.
C
C     BY ULISES REVELES, NOV. 2013
C
!     Modified - Raja July, 2015. 
!      Raja - Needs to be rewritten - It undoes all efforts to save memory. 
!
C     --- LOCAL DIMENSIONS ---
C
C     NBAS : Dimension of Slater type orbitals.
C     HFOMO: Highest fractional occupied MO.
C     LFOMO: Lowest fractional occupied MO.
C     NOMO : Number of occupied MOs.
C
C     --- LOCAL VARIABLES ---
C
C     KMA(B)    : MO coefficient matrix.
C     OCCUA(B)  : Vector with the occupation numbers of the MOs.
C     EIGVALA(B): MO energies.
C     P         : Density matrix.
C
C     ------------------------------------------------------------------
C 
C     --- COMMON VARIABLES ---
C   
      use common2,only : EGMAX,RIDT,LSYMMAX,N_CON,NFNCT,N_POS,
     %                   NSPN,NBO,PRTMOS,PRTPMAT,WFFILE
      use common5,only : PSI,EVLOCC,NWF,NWFS,EFERMI,OCCUPANCY
C
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:35 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I_POS, ICON, IDIM, IFAIL, IFNCT, ILB, ILOC, ISHDUMMY,
     & ISHELLA, ISIZE, IUB, IWF, J, J_POS, JLOC, JSPN, KLOC, L_NUC, LI,
     & LMX1, M_NUC, MCON, MU, NU, NOCUA, NOCUB, ISPN, MSPN, NWF_SAVE
       REAL*8 :: OCCP , PSIKH
C
      REAL*8 EPS
      PARAMETER (EPS=1.E-6)
      real*8 :: ffermi
      external ffermi
C
C     --- LOCAL WORK FIELDS ---
C
      DIMENSION ISIZE(3)
      DIMENSION KLOC(10,6,3),PSIKH(6)
C
C     --- LOCAL VARIABLES ---
C 
      CHARACTER*3 IATOM
      INTEGER BINDX,I,INUC,NBAS,KNT,NBASTOT
      INTEGER IMO,ISTO,JSTO,LFOMO,HFOMO,NOMO,NDIM,NPRT,OCCINT
      REAL*8 EFMIN,EFMAX,ERGMIN,ERGMAX,FACT,TEMP
      REAL*8 DIFF,OCCREA
      LOGICAL BETAMOS, FLOSIC, FLOWF
C
      DATA ISIZE/1,3,6/
C
C     --- LOCAL DYNAMICAL FIELDS ---
C
      INTEGER ALLOCATION
      REAL*8,ALLOCATABLE  :: RVECA(:,:)
      REAL*8,ALLOCATABLE  :: EIGVALA(:),EIGVALB(:),OCCUA(:),OCCUB(:),
     %                       KMA(:,:),KMB(:,:),P(:,:)
C
C     ------------------------------------------------------------------
C
C     --- FIRST BUILD MATRIX WITH MO COEFFICIENTS ---
C
C     --- CALCULATE ENERGY TRESHOLDS ---
C
      EFMIN=MIN(EFERMI(1),EFERMI(NSPN))
      EFMAX=MAX(EFERMI(1),EFERMI(NSPN))
C
C     --- GET ALL ENERGY LEVELS HERE ---
C 
      ERGMIN=EFMIN-(2000.)
      ERGMAX=EFMAX+(EGMAX)
      CALL WFWIND2(ERGMIN,ERGMAX,.TRUE.,.TRUE.,IFAIL,NBASTOT)
C
C     --- INITIALIZATION ---
C
      NPRT = NWFS(1)
      NBAS = NBASTOT


CJUR
      WRITE(*,*) 'ORBITAL ENERGIES:'
      DO I=1,NWF
      WRITE(*,*) I,EVLOCC(I)
      END DO
CJUR
      IF (NWF.GT.NWFS(1)) THEN
        IDIM = 2*NWFS(1)
      ELSE
        IDIM = NWFS(1)
      END IF
C
      BETAMOS = .FALSE.
      IF (NWF.GT.NWFS(1)) BETAMOS = .TRUE.
C
C     --- ALLOCATE LOCAL FIELDS ---
C
      ALLOCATE(RVECA(3,MX_GRP),EIGVALA(NBAS),OCCUA(NBAS),
     %         KMA(NBAS,NBAS),P(NBAS,NBAS),STAT=ALLOCATION)
      IF (BETAMOS) THEN
        ALLOCATE(EIGVALB(NBAS),OCCUB(NBAS),KMB(NBAS,NBAS),
     %  STAT=ALLOCATION)
      END IF
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR: ALLOCATION FAILED IN BLDPMAT'
      END IF
C
C     --- INITIALIZATION ---
C
      OCCUA(1:NBAS) = 0.0
      EIGVALA(1:NBAS) = 0.0
C
C     --- DO BETA MO'S IF NEEDED  ---
C
      IF (BETAMOS) THEN
        OCCUB(1:NBAS) = 0.0
        EIGVALB(1:NBAS) = 0.0
C
        ILB = NWFS(1) + 1
        IUB = 2*NWFS(1)
      END IF
C
      TEMP=0.0001D0
C
C     --- NOW RUN OVER ALL ALPHA AND BETA MOLECULAR ORBITALS WITH INDEX IWF ---
C
 
c Kamal Sharkas
C           NWF_SAVE=NWF

        INQUIRE(FILE='FRMORB',EXIST=FLOSIC)
        IF(FLOSIC)THEN
         OPEN(68,FILE='FRMORB')
         READ(68,*)NOCUA ,NOCUB
         CLOSE(68)
         ENDIF
C Kamal Sharkas

       DO 87 IWF=1,NWF
C
C     --- FIND THE NUMBER FOR THE BETA ORBITALS ---
C
        JSPN=1
        IF(IWF.GT.NWFS(1)) THEN
          JSPN=2
          BINDX = IWF - NWFS(1)
        END IF
C
C     --- INITIALIZE ATOMIC ORBITAL COUNTER ---
C   
        KNT = 0
C
C     --- OCCUPATION FOR ALPHA ORBITALS ---
C
        IF(JSPN.EQ.1) THEN
          EIGVALA(IWF) = EVLOCC(IWF)
          IF (ABS(EIGVALA(IWF)).LT.EPS) EIGVALA(IWF) = 0.0
          OCCP=FFERMI(EVLOCC(IWF),EFERMI(1),TEMP)
C          print*, 'PMAT test' ,IWF, EIGVALA(IWF), OCCP,EFERMI(1)
          OCCUA(IWF) = OCCP
C Kamal Sharkas
          IF(FLOSIC)   THEN
          IF (IWF.LE.NOCUA)  OCCUA(IWF)= 1.0
          ENDIF
C Kamal Sharkas

C     --- OCCUPATION FOR CLOSED SHELL SYSTEMS ---       
C
          IF (.NOT.BETAMOS) OCCUA(IWF) = 2.0*OCCUA(IWF)
C
C     --- NUMERICAL STABILITY TEST ---
C
          IF (ABS(OCCUA(IWF)).LT.EPS) OCCUA(IWF) = 0.0
C
C     --- BETA ORBITALS ---
C
        ELSE
          EIGVALB(BINDX) = EVLOCC(IWF)
          IF (ABS(EIGVALB(BINDX)).LT.EPS) EIGVALB(BINDX) = 0.0
          OCCP=FFERMI(EVLOCC(IWF),EFERMI(JSPN),TEMP)
          OCCUB(BINDX) = OCCP
C Kamal Sharkas
          IF(FLOSIC)  THEN
          IF(BINDX.LE.NOCUB)  OCCUB(BINDX)= 1.0
          ENDIF
C Kamal Sharkas

          IF (ABS(OCCUB(BINDX)).LT.EPS) OCCUB(BINDX) = 0.0
        END IF
C
        ISHELLA=0
        DO 86 IFNCT=1,NFNCT
          LMX1=LSYMMAX(IFNCT)+1
          DO 84 I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUMMY)
            DO 82 J_POS=1,M_NUC
              CALL UTRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     %                     RVECA,L_NUC,1)
              ILOC=0
              DO LI=1,LMX1
                DO MU=1,ISIZE(LI)
                  DO ICON=1,N_CON(LI,IFNCT)
                    ILOC=ILOC+1
                    KLOC(ICON,MU,LI)=ILOC
                  END DO
                END DO
              END DO
C
              DO LI=1,LMX1
                DO MCON=1,N_CON(LI,IFNCT)
                  DO MU=1,ISIZE(LI)
                    NU=0
                    DO ICON=1,N_CON(LI,IFNCT)
                      IF(ICON.EQ.MCON)THEN
                        JLOC=KLOC(ICON,MU,LI)
                        NU=MU   
                        PSIKH(NU)=PSI(JLOC,IWF,1)
                      END IF
                    END DO !ICON
                  END DO !MU
C
C    --- RUN OVER 'S', 'P', AND 'D' ORBITALS AND ---
C    --- SAVE MO COEFFICIENTS INTO KM MATRIX     ---
C
                  DO MU=1,ISIZE(LI)
                    KNT=KNT+1
C
                    IF(JSPN.EQ.1) THEN
                      KMA(KNT,IWF) = PSIKH(MU)
                    ELSE
                      KMB(KNT,BINDX) = PSIKH(MU)
                    END IF
C
                  END DO
C
                END DO
              END DO
   82       CONTINUE
   84     CONTINUE
   86   CONTINUE
   87 CONTINUE

C
C     --- PRINT ALPHA MOLECULAR ORBITALS ---
C
      IF (PRTMOS) THEN
        CALL PRIEIG(KMA,EIGVALA,OCCUA,NBAS,1,NBAS,1,NPRT,3,5,98,
     %              'BASIS',"ALPHA MOLECULAR ORBITALS")

C     --- PRINT BETA KM MATRIX IF NEEDED ---
C
        IF(BETAMOS) THEN
          CALL PRIEIG(KMB,EIGVALB,OCCUB,NBAS,1,NBAS,1,NPRT,3,5,98,
     %                'BASIS',"BETA MOLECULAR ORBITALS")
        END IF
      END IF
C
C     --- NOW CALCULATE DENSITY MATRIX ---
C
C     --- INITIALIZATION ---
C
      LFOMO = 0
      HFOMO = 0
      NOMO = 0
C
C     --- FIRST FIND LOWEST AND HIGHEST FRACTIONAL ---
C     --- OCCUPIED MOS: LFOMO AND HFOMO            ---  
C
      DO I=1,NBAS
        IF (OCCUA(I).GT.EPS) NOMO = NOMO + 1
        OCCINT = INT(OCCUA(I))
        OCCREA = FLOAT(OCCINT)
        DIFF = OCCREA - OCCUA(I)
C
C     --- THIS OCCUPATION IS FRACTIONAL ---
C
        IF (DIFF.GT.EPS) THEN
          IF (LFOMO.EQ.0) LFOMO = I
          IF (HFOMO.LT.LFOMO) HFOMO = I
        END IF
      END DO
C
      IF (LFOMO.EQ.0) LFOMO = NOMO + 1
      NDIM = LFOMO - 1
C
      P(1:NBAS,1:NBAS) = MATMUL(KMA(1:NBAS,1:NDIM),
     %                        TRANSPOSE(KMA(1:NBAS,1:NDIM)))
C
      FACT = 2.0
      IF (BETAMOS) FACT = 1.0
C
      DO ISTO=1,NBAS
        DO JSTO=1,NBAS
          P(ISTO,JSTO) = FACT*P(ISTO,JSTO)
        END DO
      END DO
C
C     --- ADD THE CONTRIBUTION FROM THE FRACTIONAL OCCUPIED MOS ---
C
      IF (HFOMO.GE.LFOMO) THEN
        WRITE(*,*) 'ADDING FRACTIONAL OCCUPATION'
C
        DO ISTO=1,NBAS
          DO IMO=LFOMO,HFOMO
            P(ISTO,ISTO) = P(ISTO,ISTO) + OCCUA(IMO)*KMA(ISTO,IMO)**2
          END DO
          DO JSTO=1,NBAS
            DO IMO=LFOMO,HFOMO
              P(ISTO,JSTO) = P(ISTO,JSTO) +
     %                     FACT*OCCUA(IMO)*KMA(ISTO,IMO)*KMA(JSTO,IMO)
            END DO
          END DO
        END DO
      END IF
C
C     --- PRINT DENSITY MATRIX ---
C
      IF (PRTPMAT) THEN
        CALL PRILMAT(P,NBAS,1,NBAS,1,NBAS,3,5,99,'BASIS','TOTAL',
     %               'ALPHA DENSITY MATRIX')
      END IF
C
C     --- WRITE ALPHA DENSITY MATRIX (LOWER TRIANGLE ONLY) ON TAPE ---
C
C     IF (NBO) THEN 
        OPEN(66,FILE='PAMAT',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(66)
        WRITE(66) NBAS
        DO I=1,NBAS
          DO J=1,I
            WRITE(66) P(I,J)
          END DO
        END DO
        CLOSE(66)
C     END IF      
C
C     --- NOW CALCULATE BETA DENSITY MATRIX IF NEEDED ---
C
      IF (BETAMOS) THEN
C
C     --- INITIALIZATION ---
C
        LFOMO = 0
        HFOMO = 0
        NOMO = 0
C
C     --- FIRST FIND LOWEST AND HIGHEST FRACTIONAL ---
C     --- OCCUPIED MOS: LFOMO AND HFOMO            ---  
C
        DO I=1,NBAS
          IF (OCCUB(I).GT.EPS) NOMO = NOMO + 1
          OCCINT = INT(OCCUB(I))
          OCCREA = FLOAT(OCCINT)
          DIFF = OCCREA - OCCUB(I)
C
C     --- THIS OCCUPATION IS FRACTIONAL ---
C
          IF (DIFF.GT.EPS) THEN
            IF (LFOMO.EQ.0) LFOMO = I
            IF (HFOMO.LT.LFOMO) HFOMO = I
          END IF
        END DO
C
        IF (LFOMO.EQ.0) LFOMO = NOMO + 1
        NDIM = LFOMO - 1
C
        P(1:NBAS,1:NBAS) = MATMUL(KMB(1:NBAS,1:NDIM),
     %                   TRANSPOSE(KMB(1:NBAS,1:NDIM)))
C
        DO ISTO=1,NBAS
          DO JSTO=1,NBAS
            P(ISTO,JSTO) = 1.0*P(ISTO,JSTO)
          END DO
        END DO
C
C     --- ADD THE CONTRIBUTION FROM THE FRACTIONAL OCCUPIED MOS ---
C
        IF (HFOMO.GE.LFOMO) THEN
          DO ISTO=1,NBAS
            DO IMO=LFOMO,HFOMO
              P(ISTO,ISTO) = P(ISTO,ISTO) + OCCUB(IMO)*KMB(ISTO,IMO)**2
            END DO
            DO JSTO=1,NBAS
              DO IMO=LFOMO,HFOMO
                P(ISTO,JSTO) = P(ISTO,JSTO) +
     %                     1.0*OCCUB(IMO)*KMB(ISTO,IMO)*KMB(JSTO,IMO)
              END DO
            END DO
          END DO
        END IF
C
C     --- PRINT DENSITY MATRIX ---
C
      IF (PRTPMAT) THEN
        CALL PRILMAT(P,NBAS,1,NBAS,1,NBAS,3,5,99,'BASIS','TOTAL',
     %               'BETA DENSITY MATRIX')
      END IF
C
C     --- WRITE BETA DENSITY MATRIX (LOWER TRIANGLE ONLY) ON TAPE ---
C 
C       IF (NBO) THEN 
          OPEN(67,FILE='PBMAT',FORM='UNFORMATTED',STATUS='UNKNOWN')
          REWIND(67)
          WRITE(67) NBAS
          DO I=1,NBAS
            DO J=1,I
              WRITE(67) P(I,J)
            END DO
          END DO
          CLOSE(67)
C       END IF      
C
      END IF
C
C     --- DEALLOCATE LOCAL FIELDS ---
C
      DEALLOCATE(RVECA,EIGVALA,OCCUA,KMA,P,STAT=ALLOCATION)
      IF(BETAMOS) DEALLOCATE(EIGVALB,OCCUB,KMB,STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR: DEALLOCATION FAILED IN BLDPMAT'
      END IF
C
      RETURN
C
C     ------------------------------------------------------------------
C
      END
