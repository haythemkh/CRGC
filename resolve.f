C-------------------------------------------------------------------------
C-
C       --
C-   INSTITUT  DE RECHERCHE EN INFORMATIQUE DE TOULOUSE (I.R.I.T.)      --
C-       LABORATOIRE D'INFORMATIQUE et MATHEMATIQUE APPLIQUEE
C-                DE L'E.N.S.E.E.I.H.T.
C               --
C-
C-                   P.Amestoy, M. Dayde
CC      --
C-------------------------------------------------------------------------
       PROGRAM RESOLVE
       implicit none
C ------------------------
C Declarations des donnees
C ------------------------
C NMAX  : Dimension principale de la matrice (taille maximale)
C A    : Matrice reelle de taille maximale
C B    : Vecteur second membre de taille maximale
C R    : Vecteur residu
C N    : Dimension reelle de la matrice
C EPS  : Precision sur la norme du residu
C Xseq : Vecteur solution obtenu par le code sequentiel
C Xpar : Vecteur solution obtenu par le code parallele
       INTEGER NMAX
       PARAMETER(NMAX=20000)
       DOUBLE PRECISION EPS
       PARAMETER (EPS=1.0D-14)
       INTEGER N
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: B, Xseq, Xpar, R
C ----------------------------
C Variables locales de travail 
C ----------------------------
* nb_of_slaves    : Total number of PVM processes activated
* nb_of_processes : number of processes doing computation
*                   nb_of_processes=nb_of_slaves+1 (master works!)
* first_row(I)    : Index of the first row processed by task I
*                   The number of rows processed by task I is then
*                   first_row(I+1) - first_row(I)
*
       INTEGER slave_max
       DOUBLE PRECISION small
       PARAMETER (slave_max=32)
       PARAMETER (small = 1.0D-6)
       INTEGER I,J, first_row(slave_max+2),
     *         nb_of_slaves, nb_of_processes
       INTEGER IERR
       DOUBLE PRECISION t1,t2
* --------------
* pvm variables
* --------------
       INTEGER my_id, info
* -------------------
* Procedures externes
* -------------------
       EXTERNAL JACSEQ, JACPAR, secdeb, secfin, verif
C
       write(6,*) ' Version initiale du Jacobi '
       write(6,*) ' Parallelisation statique ',
     *    ' du produit matrice-vecteur dans ' 
       write(6,*) ' le calcul du residu '
*
* --------------------------
* Enroll this program in Pvm
* --------------------------
      call pvmfmytid(my_id)

      if (my_id .lt. 0) then
         write(*,*) 'failure in enrolling on host'
         stop
      endif

      write(6,*) ' master tid = ', my_id

*
* ---------------------------------------------
* BEGIN of data initialisation/preparation
* ---------------------------------------------
*
      N = 0
      DO WHILE ((N.LT.1) .OR. (N.GT.NMAX)) 
        WRITE(6,*) ' Taille de la matrice ? <= ',  NMAX
        READ(5,*) N
      END DO
*
      nb_of_slaves = 0
      DO WHILE ((nb_of_slaves.le.0).OR.(nb_of_slaves.gt.slave_max) )
        write(*,*) ' How many slave workstations will you used '
        read(*,*) nb_of_slaves
      END DO
      nb_of_processes = nb_of_slaves +1
*     ------------------------------------------------------------
*     compute the index of the first row processed by each task
*     ------------------------------------------------------------
      first_row(1)                 = 1
      first_row(nb_of_processes+1) = N+1
      j = (N / nb_of_processes)
      if (nb_of_slaves.gt.0) then
       do i=2,nb_of_processes
        first_row(i) = j*(i-1) +1
       enddo
      endif
*

      do i=1, nb_of_processes
         write(6,*) 'i, first_row(i) ',i,first_row(i)
      enddo

C     ---------------------------
C     Allocate matrix and vectors
C     ---------------------------
      ALLOCATE( A(N,N), B(N), Xseq(N), Xpar(N), R(N), stat=IERR)
      IF ( IERR .GT. 0 ) THEN
        WRITE(6,*) " Error in memory allocation of size ", N*N+4*N, 
     &             " reduce value of N =", N
        GOTO 1000  ! quit PVM and return
      ENDIF

C     --------------------------------
C     Initialisation de la matrice 
C     avec une diagonale principale
C     strictement dominante
C     --------------------------------
      DO I=1,N
           DO J=1,N
               A(I,J)=DFLOAT(1)
           END DO
**         A(I,I)=DFLOAT(6*N+1)
           A(I,I)=DFLOAT(3*N+1)
      END DO
C
C     Initialisation du second membre
      DO I=1,N
          B(I)=DFLOAT(1)
      END DO
C
C -----------------------------------------------------------------
C Appel de la version parallele de la resolution par Gauss-Seidel
C -----------------------------------------------------------------
       write (6,*) ' *** Debut Jacobi Parallele ****'
       call secdeb(t1)
       CALL JACPAR(A,N,B,Xpar,R,EPS,
     *             nb_of_slaves, first_row)
       call secfin(t1)
       write (6,*) ' *** Fin   Jacobi Parallele ****'
C -----------------------------------------------------------------
C Appel de la version sequentielle de la resolution par Gauss-Seidel
C -----------------------------------------------------------------
       write (6,*) ' *** Debut Jacobi Sequentiel ****'
       call secdeb(t2)
       CALL JACSEQ(A,N,B,Xseq,R,EPS)
       call secfin(t2)
       write (6,*) ' *** Fin   Jacobi Sequentiel ****'
*
* --------------------------------------------------
* Compare sequential and // solution
* Print max/relative error and performance
* --------------------------------------------------
       call verif(N,Xpar,Xseq)
*
*
      if ( (t1 .lt. small) .or. (t2 .lt. small) ) then
        write(6,9000)
        write(6,9001)
        goto 1000
      endif
      write (6,3000) t2,t1,t2/t1
*
1000   continue
C      Exit pvm and deallocate arrays
       call pvmfexit(info)
       IF (ALLOCATED(A)) DEALLOCATE(A) 
       IF (ALLOCATED(B)) DEALLOCATE(B) 
       IF (ALLOCATED(Xseq)) DEALLOCATE(Xseq) 
       IF (ALLOCATED(Xpar)) DEALLOCATE(Xpar) 
       IF (ALLOCATED(R)) DEALLOCATE(R) 
*
3000  format(4x,'Seq. Time',5x,'Para. Time',4x,'Speedup'/3(1pe14.5))
9000  format(1x,' Parallel timing was too small to do calculations')
9001  format(1x,' Use a larger value for problem size and rerun')
       STOP
       END
