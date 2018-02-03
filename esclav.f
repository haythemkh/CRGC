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
************************************************************************
C Esclave associe au calcul du produit matrice-vecteur dynamique
C --------------------------------------------------------------
      program esclav
      implicit none
* -------------------
* static parameters
* -------------------
      integer NMAX 
      parameter (NMAX=20000)
      double precision one, zero
      parameter ( one = 1.0D0, zero = 0.0D0)
* --------------------------
* matrix-vector related data
*     a matrix m*n
*     x vector of size n
*     y vector of size m
* --------------------------
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: a
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: x, y
      integer n,m,i,j
      logical fini
* --------------
* pvm variables
* --------------
      integer my_id, info, IERR

      include 'fpvm3.h'
*
* ----------------------------------
* variables used for message passing
* ----------------------------------
* no_slave Logical number of slave
*          It is received and then send back to the master 
*          The master need it to know
*          which part of the solution is computed
* bufidS buffer for sending data
* bufidR buffer for receiving data
* type = 0 to reveive data broadcasted by the master
*      = 1 to reveive local data from the master
*      = 2 to send results to the master
      integer no_slave, p_id, bufidS, bufidR,type
*
*
* ---------------------------
* Enroll this program in PVM_3
* Get the tid of the master's task id
* ---------------------------
      call pvmfmytid(my_id)
      call pvmfparent(p_id)
      if (p_id.le.0) stop
*
* ------------------------
* Receive broadcasted data
* ------------------------
*
* ------------------------------------------------------------
* Receive size of the block of row followed by sub-matrix
* ------------------------------------------------------------
      type = 1
      call pvmfrecv(p_id, type, bufidR)
      call pvmfunpack(INTEGER4, no_slave, 1, 1, info)
      call pvmfunpack(INTEGER4, m, 1, 1, info)
      call pvmfunpack(INTEGER4, n, 1, 1, info)
      call pvmfunpack(BYTE1, fini, 1, 1, info)
        write(6,*) ' slave: ',no_slave,' tid: ',my_id, ' m= ',m
        ALLOCATE( A(m,n), y(m), stat=IERR)
        IF (IERR .GT. 0 ) THEN
          WRITE(6,*) " Error in memory allocation of size ",m*n+n
          GOTO 1000  ! quit PVM and return
        ENDIF
        
      do i=1,n
        call pvmfunpack(REAL8, a(1,i), m, 1, info)
      enddo

      write(6,*) ' slave: ',my_id, ' end reception submatrix'
C       Allocate all Matrices and vectors of size n
        ALLOCATE( x(n), stat=IERR)
*
      fini = .false.
      DO WHILE(.not.fini)
      
      write(6,*) ' slave: ',my_id,' starts receiving'

      type = 0
      call pvmfrecv(p_id, type, bufidR)
      call pvmfunpack(INTEGER4, n, 1, 1, info)

        write(6,*) ' slave: ',my_id, ' n= ',n

        IF (IERR .GT. 0 ) THEN
          WRITE(6,*) " Error in memory allocation of size ", n, 
     &               " reduce value of N =", n
          GOTO 1000  ! quit PVM and return
        ENDIF

      call pvmfunpack(REAL8   , x, n, 1, info)


* ---------------------------
* compute his part of work
* --------------------------

      write(6,*) ' slave: ', my_id,' starts computation '

      do i=1,m
       y(i) = zero
       do j=1,n
        y(i) = y(i) + a(i,j)*x(j)
       enddo
      enddo
*
* --------------------------------------
* send back task_id and results to master
* --------------------------------------

      write(6,*) ' slave: ', no_slave, ' tid: ',my_id,
     *           ' sends results '

      call pvmfinitsend(PVMDATADEFAULT, bufidS)
      call pvmfpack(INTEGER4, no_slave, 1, 1, info)
      call pvmfpack(REAL8, y, m, 1, info)
      type = 2
      call pvmfsend(p_id, type, info)
      
      ENDDO

*
* ------------------------
* leave Pvm before exiting
* ------------------------

 1000 write(6,*) ' slave: ', my_id, ' leaves PVM '
*     -- free memory associated to a and y
      IF (ALLOCATED(a)) DEALLOCATE(a) 
      IF (ALLOCATED(y)) DEALLOCATE(y) 
      
      IF (ALLOCATED(x)) DEALLOCATE(x) 
      call pvmfexit(info)
*
      stop
      end 
      
    
