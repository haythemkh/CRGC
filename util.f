C  --------------------------------
C  CALCUL de la norme d'un vecteur
C  --------------------------------
C
       DOUBLE PRECISION FUNCTION NORM(R,N)
C Calcul de la norme euclidienne du vecteur R, de longueur N
       DIMENSION R(N)
       DOUBLE PRECISION R
       INTEGER N, I
C Valeur provisoire de la norme
       DOUBLE PRECISION S
C Initialion
       S=DFLOAT(0)
C Calcul du produit scalaire
       DO I=1,N
           S=S+R(I)*R(I)
       END DO
C Calcul final de la norme
       NORM=DSQRT(S)
       RETURN
       END
C
C  ---------------------------------------
c  CALCUL du residu B-AX
C  ---------------------------------------
C
       SUBROUTINE RESIDU(A,LDA,X,B,R,N)
C Calcul du residu R=B-A.X, A de dimension NxN
C dimension principale de A = LDA
       INTEGER LDA, N, I, J
       DIMENSION A(LDA,*),X(N),B(N),R(N)
       DOUBLE PRECISION A,X,B,R
C Boucle sur les composantes du vecteur residu
       DO I=1,N
           R(I)=B(I)
           DO J=1,N
               R(I)=R(I)-A(I,J)*X(J)
           END DO
       END DO
       RETURN
       END
C
C  ---------------------------------
C  CALCUL l'erreur max entre deux 
C  vecteurs C et CT de taille M
C  ---------------------------------
C

        SUBROUTINE VERIF(M,CT,C)

        INTEGER  M
        double precision  C(M) , CT(M)
        double precision  ERMAX, MAXV, RELER,EPS, ONE
        INTEGER  I
        ONE  = 1.0D0
        EPS  = 1.0D-9
        ERMAX = ABS(CT(1)-C(1))
        MAXV  = ABS(C(1))
           DO  10 I=1,M
             ERMAX = MAX(ABS(CT(I)-C(I)), ERMAX)
             MAXV  = MAX(ABS(C(I)),MAXV)
 10        CONTINUE
        WRITE (6,*) 'Erreur maximum           =',ERMAX
        IF (MAXV.LE.EPS) MAXV = ONE
        RELER = ERMAX/MAXV
        WRITE (6,*) 'Erreur maximum  Relative =',RELER
        RETURN
        END

C
C  -------------------------------------------------
C  Lecture du nom des machines constituant la
C  machine virtuelle
C  -------------------------------------------------
       subroutine read_file(name)
        character*32 name
*
        character*32 machine
        integer i
 5      read(10,1000,END=100) machine
1000    format(a10)
        do 10, i=1,32
           name(i:i)=' '
10      continue
        i=1
20      continue
          if (machine(i:i).ne.' ') then
            name(i:i) = machine(i:i)
            i = i+1
            goto 20
          endif
cnot needed with gfortran        name(i:i)='\0'
        GOTO 500
******* error end of file ***************
*       we restart from beginning of file
100     rewind(10)
        GOTO 5
***********************************
 500    return
        end
C 
C*******************************************************
        subroutine secdeb(t)
C****************************************************
        real*8 t,t1
        real*8 ELAPSE8
        EXTERNAL ELAPSE8
*
        t=ELAPSE8(t1)
        return
        end
C*******************************************************
        subroutine secfin(t)
C****************************************************
        real*8 t
        real*8 ELAPSE8, ZERO
        EXTERNAL ELAPSE8
*
        ZERO = 0.0
        t=ELAPSE8(ZERO)-t
        return
        end
*
C*******************************************************
       real*8 function Max(n,x)
C*******************************************************
*
       integer n
       real*8 x(n)
*
       integer i
       Max = x(1)
       do 10, i=2,n
         if ( x(i).gt.Max ) then
           Max = x(i)
         endif
 10    continue
       return
       end
*
C*******************************************************
       real*8 function Min(n,x)
C*******************************************************
*
       integer n
       real*8 x(n)
*
       integer i
       Min = x(1)
       do 10, i=2,n
         if ( x(i).lt.Min ) then
           Min = x(i)
         endif
 10    continue
       return
       end

