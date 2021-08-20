      SUBROUTINE NBODRV
C
C     ------------------------------------------------------------------
C
C     DRIVER ROUTINE TO PRINT DENSITY MATRIX AND TO WRITE NBO ANALYSIS
C     INPUT FILE.
C
C     BY ULISES REVELES, DEC. 2013
C
C     --- LOCAL DIMENSIONS ---
C
C     ------------------------------------------------------------------
C
      use xmol, only   : NUM_ATMS,XMOL_LIST,GET_LETTER,ATOMORB_LAB
      use common2,only : NBO,NFNCT,PRTSMAT
      use common5,only : NWF,NWFS
      Lind(I,J) = (max(I,J)*(max(I,J)-1))/2 + min(I,J) 
C
      INTEGER I,INUC,J,K,L,NBAS,NDUMMY,STREXT
C
      REAL*8 AU2ANG,GRAD, tmp,tmp2
      DATA AU2ANG/0.529177D0/
C
      CHARACTER*2 LETTER
      CHARACTER*10 LABNBAS
      CHARACTER*10 LABNATOMS
      CHARACTER*80 LINE,STRCOMP
C
      LOGICAL BETAMOS
C
C     --- LOCAL DYNAMICAL FIELDS ---
C
      INTEGER ALLOCATION
      INTEGER,ALLOCATABLE :: BASVEC(:,:),MOLPTR(:)
      REAL*8,ALLOCATABLE  :: OVERTOT(:,:),PMAT(:,:),WORK(:,:)
      REAL*8,ALLOCATABLE  :: OVERUNM(:,:)
      Real*8, Allocatable :: LTArray(:)
C
C     ------------------------------------------------------------------
C
C     --- READ NUMBER OF BASIS FUNCTIONS ---
C
      OPEN(66,FILE='OVLNRM',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND(66)
      READ(66) NDUMMY
C
C     --- INITIALIZATION ---
C
c      NBAS = NWFS(1)
       NBAS = NDUMMY 
      print*,'NBOdrv',NBAS
CJUR
!       NWF = 14  !just test
!       NWFS(1)=7 ! just test 
       print*, 'NFNCT in nbodrv is:', NFNCT
       print*, 'NWF in nbodrv is:', NWF
       print*, 'NWFS(1) in nbodrv is:', NWFS(1)
CJUR
C
C     --- CHECK IF BETA ORBITALS ARE NEEDED ---
C
      BETAMOS = .FALSE.
      IF (NWF.GT.NWFS(1)) BETAMOS = .TRUE.
C
C     --- ALLOCATE LOCAL FIELDS ---
C
      ALLOCATE(OVERTOT(NBAS,NBAS),STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR: ALLOCATION FAILED IN NBODRV'
      END IF
C
C     --- READ OVERLAP MATRIX (LOWER TRIANGLE ONLY) ---
C     --- AND COMPLETE UPPER TRIANGLE OF MATRIX     ---
C
      DO I=1,NBAS
        DO J=I, NBas
C           OVERTOT(J,I)  = HStor(Lind(I,J), 1)
          READ(66) OVERTOT(I,J)
          OVERTOT(J,I) = OVERTOT(I,J)
        END DO
      END DO
      CLOSE(66)
C
C      DO I=1,NBAS
C        DO J=I+1,NBAS
C          OVERTOT(I,J) = OVERTOT(J,I)
C        END DO
C      END DO
C
C     --- PRINT OVERLAP MATRIX ---
C
       print*,'Before PRILMAT'
      IF (PRTSMAT) THEN
        CALL PRILMAT(OVERTOT,NBAS,1,NBAS,1,NBAS,3,5,99,'BASIS','TOTAL',
     %              'OVERLAP MATRIX')
      END IF
       print*,'After PRILMAT'
C
C     --- QUICK RETURN IN NOT NBO ANALYSIS IS REQUESTED ---
C
      IF (.NOT.NBO) THEN
        DEALLOCATE(OVERTOT,STAT=ALLOCATION)
        DEALLOCATE(ATOMORB_LAB,STAT=ALLOCATION)
        IF (ALLOCATION.GT.0) THEN
          WRITE(*,*) 'ERROR: DEALLOCATION FAILED IN NBODRV'
        END IF
C
        RETURN
      END IF
C
C     --- ALLOCATE LOCAL FIELDS ---
C
c       call flush(6)
      ALLOCATE(BASVEC(NBAS,2),MOLPTR(NBAS),PMAT(NBAS,NBAS),
     %         WORK(NBAS,NBAS),STAT=ALLOCATION)
      ALLOCATE(OVERUNM(NBAS,NBAS),STAT=ALLOCATION)
      Allocate(LTArray(NBAS*(NBAS+1)/2) )
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR: ALLOCATION FAILED IN NBODRV'
      END IF
C
C     --- PREPARE POINTER FOR ORBITAL REORDERING                     ---
C     ---                                                            ---
C     --- NRLMOL                         NBO                         ---
C     --  DXX, DYY, DZZ, DXY, DXZ, DYZ   DXX, DXY, DXZ, DYY, DYZ, DZZ --
C
      
c       call flush(6)
C      IF(.not. ALLOCATED(ATOMORB_LAB)) 
C     %  ALLOCATE(ATOMORB_LAB(NDH),STAT=IERR)
      MOLPTR(1:NBAS) = 0  
      DO I=1,NBAS
        IF (ATOMORB_LAB(I)%ORBLAB(3:3).EQ.'S') THEN
          MOLPTR(I) = I
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:3).EQ.'P') THEN
          MOLPTR(I) = I
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:5).EQ.'DXX') THEN
          MOLPTR(I) = I
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:5).EQ.'DYY') THEN
          MOLPTR(I+2) = I
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:5).EQ.'DZZ') THEN
          MOLPTR(I+3) = I
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:5).EQ.'DXY') THEN
          MOLPTR(I-2) = I
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:5).EQ.'DXZ') THEN
          MOLPTR(I-2) = I
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:5).EQ.'DYZ') THEN
          MOLPTR(I-1) = I
        END IF
      END DO
c       call flush(6)
CJUR   
      WRITE(*,*) 'THIS IS MY POINTER'
      DO I=1,NBAS
       WRITE(*,*) I,MOLPTR(I)
      END DO
CJUR
C     --- NOW REORDER OVERLAP MATRIX ---
C
c       call flush(6)
      WORK(1:NBAS,1:NBAS) = OVERTOT(1:NBAS,1:NBAS)
CJEP make a copy of Overlap
      OVERUNM(1:NBAS,1:NBAS) = OVERTOT(1:NBAS,1:NBAS)
      OVERTOT(1:NBAS,1:NBAS) = 0.0
      DO I=1,NBAS
        DO J=1,NBAS
          OVERTOT(I,J) = WORK(MOLPTR(I),MOLPTR(J))
        END DO
      END DO
C
C     --- OPEN NBODRV FILE ---
C
c       call flush(6)
      OPEN(99,FILE='CLUSTER.NBO',FORM='FORMATTED',STATUS='UNKNOWN')
C
C     --- PREPARE HEADER FOR NBO FILE ---
C
      WRITE (LABNATOMS,'(I3)') NUM_ATMS
      LABNATOMS = STRCOMP(LABNATOMS)
      LENNUM = STREXT(LABNATOMS)
      LABNATOMS = 'NATOMS='//LABNATOMS(1:LENNUM)
C
      WRITE (LABNBAS,'(I6)') NBAS
      LABNBAS = STRCOMP(LABNBAS)
      LENNUM = STREXT(LABNBAS)
      LABNBAS = 'NBAS='//LABNBAS(1:LENNUM)
C
      IF (BETAMOS) THEN
          WRITE (LINE,5000) '$GENNBO  ',LABNATOMS,LABNBAS,
     $                      ' OPEN UPPER BODM $END'
          LINE = STRCOMP(LINE)
         ELSE
          WRITE (LINE,5000) '$GENNBO  ',LABNATOMS,LABNBAS,
     $                      ' UPPER BODM $END'
          LINE = STRCOMP(LINE)
      END IF
C
C     --- WRITE HEADER ---
C
      WRITE (99,'(A)') LINE
      WRITE (99,'(A)') '$NBO  NAOMO=PVAL NPA $END'
C
C     --- WRITE CARTESIAN COORDINATES INTO NBO FILE ---
C
      WRITE (99,'(A)') '$COORD'
      WRITE (99,'(A)') '***'
      DO INUC=1,NUM_ATMS
        CALL GET_LETTER(XMOL_LIST(INUC)%ANUM,LETTER)
        WRITE(99,5100) XMOL_LIST(INUC)%ANUM,XMOL_LIST(INUC)%ANUM,
     %             XMOL_LIST(INUC)%RX*AU2ANG,XMOL_LIST(INUC)%RY*AU2ANG, 
     %             XMOL_LIST(INUC)%RZ*AU2ANG
      END DO
      WRITE (99,'(A)') '$END'
C
C     --- WRITE BASIS SET ---
C
      WRITE (99,'(A)') '$BASIS'
C
      BASVEC(1:NBAS,1:2) = 0
      DO I=1,NBAS
        BASVEC(I,1) = ATOMORB_LAB(I)%ATNUM
C
        IF (ATOMORB_LAB(I)%ORBLAB(3:3).EQ.'S') THEN
          L = 0
          K = 0
          J = 1
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:3).EQ.'P') THEN
          L = 1
          K = 0
          IF (ATOMORB_LAB(I)%ORBLAB(4:4).EQ.'X') THEN
            J = 1
          ELSE IF (ATOMORB_LAB(I)%ORBLAB(4:4).EQ.'Y') THEN
            J = 2
          ELSE IF (ATOMORB_LAB(I)%ORBLAB(4:4).EQ.'Z') THEN
            J = 3
          END IF
        ELSE IF (ATOMORB_LAB(I)%ORBLAB(3:3).EQ.'D') THEN
          L = 2
          K = 0
          IF (ATOMORB_LAB(I)%ORBLAB(4:5).EQ.'XX') THEN
            J = 1
          ELSE IF (ATOMORB_LAB(I)%ORBLAB(4:5).EQ.'YY') THEN
            J = 2
          ELSE IF (ATOMORB_LAB(I)%ORBLAB(4:5).EQ.'ZZ') THEN
            J = 3
          ELSE IF (ATOMORB_LAB(I)%ORBLAB(4:5).EQ.'XY') THEN
            J = 4
          ELSE IF (ATOMORB_LAB(I)%ORBLAB(4:5).EQ.'XZ') THEN
            J = 5
          ELSE IF (ATOMORB_LAB(I)%ORBLAB(4:5).EQ.'YZ') THEN
            J = 6
          END IF
        END IF
C
        BASVEC(I,2) = L*100 + K + J
      END DO
C
      WRITE (99,5150) 'CENTER = ',(BASVEC(I,1),I=1,NBAS)
      WRITE (99,5150) ' LABEL = ',(BASVEC(I,2),I=1,NBAS)
      WRITE (99,'(A)') '$END'
C
C     --- WRITE OVERLAP MATRIX (TRIANGLE BLOCK) INTO NBODRV --- 
C
      LINE = '$OVERLAP ! Overlap matrix elements in the AO basis'
      WRITE (99,'(A)') LINE
CJEP 
      DO J=1,NBAS
         DO I=J, NBas
         LTArray(Lind(I,J)) = Overtot(I,J)
         ENDDO
      ENDDO
      WRITE (99,5200) (LTArray(I), I=1, NBAS*(NBAS+1)/2 )

C
C      DO J=1,NBAS
C        WRITE (99,5200) (OVERTOT(I,J),I=1,J)
C      END DO
      WRITE (99,'(A)') '$END'
C
C     --- WRITE DENSITY MATRIX INTO NBODRV ---
C
      LINE = '$DENSITY ! Density matrix elements in the AO basis'
      WRITE (99,'(A)') LINE
C
C     --- READ ALPHA DENSITY MATRIX ---
C
      OPEN(67,FILE='PAMAT',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND(67)
      READ(67) NDUMMY
C
      DO I=1,NBAS
        DO J=1,I
          READ(67) PMAT(I,J)
          PMAT(J,I) = PMAT(I,J)
        END DO
      END DO
      CLOSE(67)
C
C      DO I=1,NBAS
C        DO J=I+1,NBAS
C          PMAT(I,J) = PMAT(J,I)
C        END DO
C      END DO
C
C     --- REORDER ALPHA DENSITY MATRIX ---
C
      WORK(1:NBAS,1:NBAS) = PMAT(1:NBAS,1:NBAS)
      PMAT(1:NBAS,1:NBAS) = 0.0
      DO I=1,NBAS
        DO J=1,NBAS
          PMAT(I,J) = WORK(MOLPTR(I),MOLPTR(J))
        END DO
      END DO
C
CJEP 
      DO J=1,NBAS
         DO I=J, NBas
         LTArray(Lind(I,J)) = PMat(I,J)
         ENDDO
      ENDDO
      WRITE (99,5200) (LTArray(I), I=1, NBAS*(NBAS+1)/2 )

C      DO J=1,NBAS
C        WRITE (99,5200) (PMAT(I,J),I=1,J)
C      END DO
      IF (.NOT.BETAMOS) WRITE (99,'(A)') '$END'




CJEP check if the problem is writing 
c      call flush(6)
      tmp = 0.0D0
      tmp2 = 0.0D0
      DO J=1,NBAS
        DO I=1,NBAS
        tmp = tmp + PMAT(I,J)*OverTot(I,J)
        tmp2 = tmp2 + Work(I,J)*OVERUNM(I,J)
        END DO
      END DO
      Write(6,*) "Checking NBO Driver code."
      Print "(1X, A, F16.8)", "Tr(P.S) before re-order = ", Tmp2
      Print "(1X, A, F16.8)", "Tr(P.S) after  re-order = ", Tmp
C      Write(6,*)' # of Electrons from Tr(P.S) after  re-order = ', Tmp
C      Write(6,*)' Diagonal of the Overlap Matrix: '
C      Write(6,'(4e15.6 )') ( OVERUNM(i,i), i=1,NBas )

C
C     --- READ BETA DENSITY MATRIX ---
C
      IF (BETAMOS) THEN
        OPEN(68,FILE='PBMAT',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(68)
        READ(68) NDUMMY
C
        DO I=1,NBAS
          DO J=1,I
            READ(68) PMAT(I,J)
          END DO
        END DO
        CLOSE(68)
C
        DO I=1,NBAS
          DO J=I+1,NBAS
            PMAT(I,J) = PMAT(J,I)
          END DO
        END DO
C
C     --- REORDER BETA DENSITY MATRIX ---
C
        WORK(1:NBAS,1:NBAS) = PMAT(1:NBAS,1:NBAS)
        PMAT(1:NBAS,1:NBAS) = 0.0
        DO I=1,NBAS
          DO J=1,NBAS
            PMAT(I,J) = WORK(MOLPTR(I),MOLPTR(J))
          END DO
        END DO
C
CJEP 
      DO J=1,NBAS
         DO I=J, NBas
         LTArray(Lind(I,J)) = PMat(I,J)
         ENDDO
      ENDDO
      WRITE (99,5200) (LTArray(I), I=1, NBAS*(NBAS+1)/2 )

      tmp = 0.0D0
      tmp2 = 0.0D0
      DO J=1,NBAS
        DO I=1,NBAS
        tmp = tmp + PMAT(I,J)*OverTot(I,J)
        tmp2 = tmp2 + Work(I,J)*OVERUNM(I,J)
        END DO
      END DO
      Write(6,*) "Checking NBO Driver code."
      Print "(1X, A, F16.8)", "Tr(P.S) before re-order = ", Tmp2
      Print "(1X, A, F16.8)", "Tr(P.S) after  re-order = ", Tmp


C       DO J=1,NBAS
C          WRITE (99,5200) (PMAT(I,J),I=1,J)
C        END DO
        WRITE (99,'(A)') '$END'
      END IF
C
C     --- CLOSE NBO FILE ---
C
      CLOSE (99)
C
C     --- DEALLOCATE LOCAL FIELDS ---
C
      DEALLOCATE(BASVEC,MOLPTR,OVERTOT,PMAT,WORK,STAT=ALLOCATION)
      DEALLOCATE(OVERUNM,STAT=ALLOCATION)
      DeAllocate(LTArray)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR: DEALLOCATION FAILED IN NBODRV'
      END IF
C
C     --- DEALLOCATE ATOMORB_LAB VECTOR ---
C
      DEALLOCATE(ATOMORB_LAB,STAT=ALLOCATION)
C
C     --- DEFINITION OF FORMATS ---
C
 5000 FORMAT (A,A10,A10,A)
 5100 FORMAT ((2I4,2X,3(F10.6,4X)))
 5150 FORMAT (A,(T9,14(I3,2X)))
 5200 FORMAT ((2X,4(E14.7,2X)))
C
      RETURN
C
C     ------------------------------------------------------------------
C
      END
