       SUBROUTINE JACPAR(A,N,B,X,R,EPS, 
     *             nb_of_slaves, first_row )
       implicit none
C -----------------------------------------------------------
C Resolution du systeme A.X=B de dimension NxN, 
C precision EPS sur la norme du residu
C Version utilisant le produit matrice vacteur parallele
c                      ---------------------------------
C -----------------------------------------------------------
C Description des parametres:
C --------------------------
C   A             : Matrice initiale d'ordre N
C   B             : Vecteur second membre de taille N
C   X             : Vecteur solution de taille N
C   R             : Vecteur residu (b-Ax) de taille N
C
C   Parametres associes a la version parallele:
C   ------------------------------------------
*   nb_of_slaves  : Number of PVM processes activated
*   first_row     : INTEGER array of size nb_of_slaves+2
*                   first_row(I) holds the  index of the first row 
*                   processed by task I
*
       INTEGER N, j
       DOUBLE PRECISION  A(N,N),B(N),X(N),EPS
       INTEGER nb_of_slaves
       INTEGER first_row(nb_of_slaves+2)
       DOUBLE PRECISION R(N)
       INTEGER MAXIT
       PARAMETER (MAXIT=500)


       INTEGER slave_max
       PARAMETER (slave_max=32)
       character*32 HFILE, MACHI(slave_max)



      include 'fpvm3.h'

C 
C -----------------
C Variables locales
C -----------------
C   E    : Norme du residu
C   K    : Compteur d'iterations
       DOUBLE PRECISION E
       INTEGER K, I, INFO
       INTEGER FLAG, numt, inst(slave_max), nb_of_processes, bufidS
       INTEGER no_slave, type
C
C ------------------------------
C Procedures/fonctions utilisees
C --------------------------------
C   RESIDU  : Procedure de calcul du residu
C   NORM    : Fonction de calcul de la norme euclidienne
       EXTERNAL RESIDU
       DOUBLE PRECISION NORM
       EXTERNAL NORM
C
C Lecture dans MACHI(I) des entrees d un fichier de machines
C
       IF (nb_of_slaves.GT.slave_max)  THEN
        WRITE(6,*) ' Number of slaves should smaller than ', slave_max
        STOP
       ENDIF
       WRITE(6,*) 'Nom du fichier host:'
       read(5,*) HFILE
       open(unit =10, file =HFILE, STATUS='OLD',IOSTAT=info)
       IF (info.GT.0) THEN
        write(6,*) ' ** ERROR opening File ',HFILE,'does not exist??'
        STOP
       ENDIF
       DO I=1,nb_of_slaves
        call read_file(MACHI(I))
       ENDDO
       rewind(10)
       close (10)
C -------------
C DEBUT du code
C -------------
* --------------------------------------------------
* initialize nb_of_slaves instances of slave program
* --------------------------------------------------
      write(6,*) ' Create PVM slaves '
      nb_of_processes=nb_of_slaves+1
      if ((nb_of_slaves.LE.0).OR.(nb_of_slaves.GT.slave_max)) then
        write(6,*) ' ERROR with nb_of_slaves in RESPAR '
        RETURN
      endif
C
      if (nb_of_slaves.ge.1) then
       numt = 0
       FLAG=PVMTASKHOST
       DO i=1, nb_of_slaves
           call pvmfspawn('esclav',FLAG, MACHI(i), 1,
     *     inst(i),info)
           numt = numt + info
           if (info.ne.1) then
            write(6,*) ' failure spawning process no: ',i, 
     &     ' on machine ', MACHI(I)
            write(6,*) ' Error code is ', inst(i)
           endif
       ENDDO
       if (numt .ne. nb_of_slaves) then
         write(*,*) 'failure in spawning one slave '
         stop
       endif
      endif

*
* ------------------------------------
* send submatrix to each slave process
* ------------------------------------

        write(6,*) ' first block of size   : ',
     *       first_row(2)-first_row(1), ' processed by master'

       type = 1
       do no_slave = 1, nb_of_processes-1
*       number of components computed by the slave no_slave
        j = first_row(no_slave+2) - first_row(no_slave+1)
        call pvmfinitsend( PVMDATADEFAULT, bufidS)
        call pvmfpack(INTEGER4, no_slave, 1, 1, info)
        call pvmfpack(INTEGER4, j, 1, 1, info)
        call pvmfpack(INTEGER4, n, 1, 1, info)
        call pvmfpack(BYTE1, 1, 1, 1, info)
        do i=1, n
         call pvmfpack(REAL8, A(first_row(no_slave+1),i),j,1,info)
        enddo
        call pvmfsend(inst(no_slave), type, info)

        write(6,*) ' send submatrix of size: ',j,
     *      ' rows to slave: ', no_slave, ' tid= ',inst(no_slave)

       enddo
       
* --------------------------------------------------
C      Initialisation du compteur d'iterations
       K=0
C      Initialisation du vecteur iterant
       DO I=1,N
           X(I)=DFLOAT(1)
       END DO
      
C      Calcul du premier residu
       CALL RESPAR(A,N,X,B,R,
     *             nb_of_slaves, first_row, MACHI, FLAG, numt, 
     *             inst, nb_of_processes, bufidS, no_slave, info, type,
     *             0)
C      Calcul de la premiere norme du residu
       E=NORM(R,N)
C      Affichage de la premiere norme du residu
       WRITE(*,*)'***** Iteration ',K,' norme du residu ',E
       
C ---------------------
C Boucle de convergence
C ---------------------
1      DO WHILE ((E.GT.EPS).AND.(K.LT.MAXIT))
C                                       -1
C          mise a jour de X: X <-- X + D   * R
C          ou D est la diagonale de la matrice
           DO I = 1, N
            X(I) = X(I) + R(I)/A(I,I)
           ENDDO
C          Calcul du nouveau residu
           CALL RESPAR(A,N,X,B,R, 
     *             nb_of_slaves, first_row, MACHI, FLAG, numt, 
     *             inst, nb_of_processes, bufidS, no_slave, info, type,
     *             0)
C          Calcul de la nouvelle norme du residu
           E=NORM(R,N)
C          Incrementation de la norme du residu
           K=K+1
C          Affichage de la norme du residu
           WRITE(*,*)'***** Iteration ',K,' norme du residu ',E
       ENDDO

           CALL RESPAR(A,N,X,B,R, 
     *             nb_of_slaves, first_row, MACHI, FLAG, numt, 
     *             inst, nb_of_processes, bufidS, no_slave, info, type,
     *             1)
        

       RETURN
       END
C
C
