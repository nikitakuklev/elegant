C+++
C     ..................................................................
C
C        SUBROUTINE QSF
C
C        PURPOSE
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN
C           EQUIDISTANT TABLE OF FUNCTION VALUES.
C
C        USAGE
C           CALL QSF (H,Y,Z,NDIM)
C
C        DESCRIPTION OF PARAMETERS
C           H      - THE INCREMENT OF ARGUMENT VALUES.
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE
C                    IDENTICAL WITH Y.
C           NDIM   - THE DIMENSION OF VECTORS Y AND Z.
C
C        REMARKS
C           NO ACTION IN CASE NDIM LESS THAN 3.
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY
C           MEANS OF SIMPSONS RULE TOGETHER WITH NEWTONS 3/8 RULE OR A
C           COMBINATION OF THESE TWO RULES. TRUNCATION ERROR IS OF
C           ORDER H**5 (I.E. FOURTH ORDER METHOD). ONLY IN CASE NDIM=3
C           TRUNCATION ERROR OF Z(2) IS OF ORDER H**4.
C           FOR REFERENCE, SEE
C           (1) F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,
C               MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.71-76.
C           (2) R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963,
C               PP.214-221.
C
C     ..................................................................
C---
      	SUBROUTINE QSF(H,Y,Z,NDIM)
	IMPLICIT REAL*8		(A-H,O-Z)
      	DIMENSION		Y(NDIM),Z(NDIM)
C
      	HT=1.0D0/3.0D0*H
	IF(NDIM-5)7,8,1
C
C  Ndim is greater than 5. preparations of integration loop
C
    1 	SUM1 = 4.0D0 * Y(2)
      	SUM1 = HT * (Y(1) + SUM1 + Y(3))
      	AUX1 = 4.0D0 * Y(4)
      	AUX1 = SUM1 + HT * (Y(3) + AUX1 + Y(5))
      	AUX2 = HT *(Y(1)+3.875D0*(Y(2)+Y(5))+2.625D0*(Y(3)+Y(4))+Y(6))
      	SUM2 = 4.0D0 * Y(5)
      	SUM2 = AUX2-HT*(Y(4)+SUM2+Y(6))
      	Z(1) = 0.0D0
      	AUX  = 4.0D0 * Y(3)
      	Z(2) = SUM2-HT*(Y(2)+AUX+Y(4))
      	Z(3) = SUM1
      	Z(4) = SUM2
      	IF(NDIM-6)5,5,2
C
C  Integration loop
C
    2 	DO 4 I=7,NDIM,2
      		SUM1   = AUX1
      		SUM2   = AUX2
      		AUX1   = 4.0D0*Y(I-1)
      		AUX1   = SUM1+HT*(Y(I-2)+AUX1+Y(I))
      		Z(I-2) = SUM1
      		IF(I-NDIM)3,6,6
C
    3 		AUX2 = 4.0D0*Y(I)
      		AUX2=SUM2+HT*(Y(I-1)+AUX2+Y(I+1))
    4 		Z(I-1) = SUM2
    5 		Z(NDIM-1) = AUX1
      		Z(NDIM) = AUX2
      	RETURN
    6 	Z(NDIM-1) = SUM2
      	Z(NDIM) = AUX1
      	RETURN
C
C  End of integration loop
C
    7 	IF(NDIM-3)12,11,8
C
C  Ndim is equal to 4 or 5
C
    8 	SUM2 = 1.125D0*HT*(Y(1) + 3.0D0*Y(2) + 3.0D0*Y(3) + Y(4))
      	SUM1 = 4.0D0*Y(2)
      	SUM1 = HT*(Y(1)+SUM1+Y(3))
      	Z(1) = 0.0D0
      	AUX1 = 4.0D0*Y(3)
      	Z(2) = SUM2-HT*(Y(2)+AUX1+Y(4))
      	IF(NDIM-5)10,9,9
    9 	AUX1 = 4.0D0*Y(4)
      	Z(5) = SUM1+HT*(Y(3)+AUX1+Y(5))
   10 	Z(3) = SUM1
      	Z(4) = SUM2
      	RETURN
C
C  Ndim is equal to 3
C
   11 	SUM1 = HT*(1.25D0*Y(1) + 2.0D0*Y(2) - .25D0*Y(3))
      	SUM2 = 4.0D0*Y(2)
      	Z(3) = HT*(Y(1)+SUM2+Y(3))
      	Z(1) = 0.0D0
      	Z(2) = SUM1
   12 	RETURN
      	END
