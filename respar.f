C
C
       SUBROUTINE RESPAR(A,N,X,B,R,
     *             nb_of_slaves, first_row, MACHI, FLAG, numt,
     *             inst, nb_of_processes, bufidS, no_slave, info, type,
     *             fini)
       implicit none
C ---------------------------------------------
C Calcul du residu R=B-A.X, A de dimension NxN
C Le produit A.x est effectue en parallele. 
C                            ------------
C   1/ Le resultat du produit // de A par X est range dans R
C   2/ R est alors mis a jour: R = B-R
C ---------------------------------------------
C Description des parametres:
C --------------------------
C   A             : Matrice initiale de taille NxN
C   B             : Vecteur de taille N
C   X             : Vecteur solution de taille N
C   R             : Vecteur residu (B-Ax)
C
C Parametres associes a la version parallele:
*   nb_of_slaves  : Total number of PVM processes activated
*                   nb_of_slaves > 0
*   first_row(I)  : Index of the first row processed by task I
*                   The number of rows processed by task I is then
*                   first_row(I+1) - first_row(I)

       INTEGER          N, I, J, fini
       DOUBLE PRECISION A(N,N),X(N),B(N),R(N)
       INTEGER          nb_of_slaves
       INTEGER          first_row(nb_of_slaves+2)
       DOUBLE PRECISION zero
       PARAMETER        (ZERO=0.0D0)
       INTEGER slave_max
       PARAMETER (slave_max=32)
       character*32 HFILE, MACHI(slave_max)
       INTEGER FLAG, numt, inst(slave_max), nb_of_processes
       INTEGER bufidS, no_slave, info, type
* --------------------
* PVM local variables
* --------------------
* inst : inst(I) is the PVM tid of slave process I 
* FLAG : specifies process spawning options
*        default is PVMTASKHOST
* type = 0 to broadcast initial informations
*      = 1 to distribute data to the slaves
*      = 2 to receive results from the slaves
* nb_of_processes : number of processes doing computation
*                    nb_of_processes=nb_of_slaves+1 (master works!)
*
       INTEGER bufidR

      include 'fpvm3.h'

        IF (fini.GT.0) THEN
            call pvmfpack(BYTE1, 1, 1, info)
            RETURN
        ENDIF

*
* -------------------------------
* broadcast n and x to the slaves
* -------------------------------
       write(6,*) ' broadcast n and x to all slaves '

       type =0
       call pvmfinitsend( PVMDATADEFAULT, bufidS)
       call pvmfpack(INTEGER4, n, 1, 1, info)
       call pvmfpack(REAL8,    x, n, 1, info)
       call pvmfmcast(nb_of_slaves, inst, type, info)
*
* ----------------------------
* Compute his part of work
* The  1st block of row of size 
* first_row(2)-first_row(1)
* is processed by the master
* ----------------------------
*
      do i=1, first_row(2)-first_row(1)
       R(i) = zero
       do j=1,N
        R(i) = R(i) + A(i,j) * X(j)
       enddo
      enddo
* 
* ---------------------
* Collect the results
* FIFO in reception
* ---------------------
*     nb_of_slaves.gt.0
      type = 2

      write(*,*) ' start to receive results from slaves '

      do i = 1, nb_of_slaves
*       from any slave
        call pvmfrecv(-1, type, bufidR)
        call pvmfunpack(INTEGER4, no_slave, 1, 1, info)
*       j = number of rows computed by slave no_slave
        j = first_row(no_slave+2) - first_row(no_slave+1)
        call pvmfunpack(REAL8,R(first_row(no_slave+1)),j,1,info)

        write(*,*) ' Reception results from slave  ', no_slave, 
     *    ' tid = ', inst(no_slave) , 
     *    ' nb Rows processed :', j

      enddo

* -------------------------
* End of // computation
* of R = A.X
* -------------------------
C
C
C prise en compte de B pour faire 
C le calcul final du residu.
       do i=1,N
            R(i)=B(i)-R(i)
       end do
C
9999   continue
       RETURN
       END
C
C

