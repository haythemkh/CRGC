       SUBROUTINE JACSEQ(A,N,B,X,R,EPS)
C -----------------------------------------------------------
C Resolution du systeme A.X=B de dimension NxN, dimension
C principale de A = N, precision EPS sur la norme du residu
C Version sequentielle de reference
C         -------------------------
C -----------------------------------------------------------
       INTEGER N
       DOUBLE PRECISION A(N,N),B(N),X(N),EPS
C Vecteur residu
       DIMENSION R(N)
       DOUBLE PRECISION R
C Norme du residu
       DOUBLE PRECISION E
C Compteur d'iterations
       INTEGER K
C Procedure de calcul du residu
       EXTERNAL RESIDU
C Fonction de calcul de la norme euclidienne
       EXTERNAL NORM
       DOUBLE PRECISION NORM
       INTEGER I
       INTEGER MAXIT
       PARAMETER (MAXIT=500)
C Initialisation du compteur d'iterations
       K=0
C Initialisation du vecteur iterant
       DO I=1,N
           X(I)=DFLOAT(1)
       END DO
C Calcul du premier residu
       CALL RESIDU(A,N,X,B,R,N)
C Calcul de la premiere norme du residu
       E=NORM(R,N)
C Affichage de la premiere norme du residu
       WRITE(*,*)'***** Iteration ',K,' norme du residu ',E
C Boucle de convergence
       DO WHILE ((E.GT.EPS).AND.(K.LT.MAXIT)) 
C                                       -1
C          mise a jour de X: X <-- X + D   * R 
C          ou D est la diagonale de la matrice
           DO I = 1, N
            X(I) = X(I) +R(I)/A(I,I)
           ENDDO
C          Calcul du nouveau residu
           CALL RESIDU(A,N,X,B,R,N)
C          Calcul de la nouvelle norme du residu
           E=NORM(R,N)
C          Incrementation de la norme du residu
           K=K+1
C          Affichage de la norme du residu
           WRITE(*,*)'***** Iteration ',K,' norme du residu ',E
       ENDDO
C
       RETURN
       END
