      SUBROUTINE CHKDERLPS(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)       
      INTEGER M,N,LDFJAC,MODE                                           
      DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),XP(N),FVECP(M),      
     *                 ERR(M)                                           

	double precision dercalc

C     **********                                                        
C                                                                       
C     SUBROUTINE CHKDER                                                 
C                                                                       
C     THIS SUBROUTINE CHECKS THE GRADIENTS OF M NONLINEAR FUNCTIONS     
C     IN N VARIABLES, EVALUATED AT A POINT X, FOR CONSISTENCY WITH      
C     THE FUNCTIONS THEMSELVES. THE USER MUST CALL CHKDER TWICE,        
C     FIRST WITH MODE = 1 AND THEN WITH MODE = 2.                       
C                                                                       
C     MODE = 1. ON INPUT, X MUST CONTAIN THE POINT OF EVALUATION.       
C               ON OUTPUT, XP IS SET TO A NEIGHBORING POINT.            
C                                                                       
C     MODE = 2. ON INPUT, FVEC MUST CONTAIN THE FUNCTIONS AND THE       
C                         ROWS OF FJAC MUST CONTAIN THE GRADIENTS       
C                         OF THE RESPECTIVE FUNCTIONS EACH EVALUATED    
C                         AT X, AND FVECP MUST CONTAIN THE FUNCTIONS    
C                         EVALUATED AT XP.                              
C               ON OUTPUT, ERR CONTAINS MEASURES OF CORRECTNESS OF      
C                          THE RESPECTIVE GRADIENTS.                    
C                                                                       
C     THE SUBROUTINE DOES NOT PERFORM RELIABLY IF CANCELLATION OR       
C     ROUNDING ERRORS CAUSE A SEVERE LOSS OF SIGNIFICANCE IN THE        
C     EVALUATION OF A FUNCTION. THEREFORE, NONE OF THE COMPONENTS       
C     OF X SHOULD BE UNUSUALLY SMALL (IN PARTICULAR, ZERO) OR ANY       
C     OTHER VALUE WHICH MAY CAUSE LOSS OF SIGNIFICANCE.                 
C                                                                       
C     THE SUBROUTINE STATEMENT IS                                       
C                                                                       
C       SUBROUTINE CHKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)     
C                                                                       
C     WHERE                                                             
C                                                                       
C       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER        
C         OF FUNCTIONS.                                                 
C                                                                       
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER        
C         OF VARIABLES.                                                 
C                                                                       
C       X IS AN INPUT ARRAY OF LENGTH N.                                
C                                                                       
C       FVEC IS AN ARRAY OF LENGTH M. ON INPUT WHEN MODE = 2,           
C         FVEC MUST CONTAIN THE FUNCTIONS EVALUATED AT X.               
C                                                                       
C       FJAC IS AN M BY N ARRAY. ON INPUT WHEN MODE = 2,                
C         THE ROWS OF FJAC MUST CONTAIN THE GRADIENTS OF                
C         THE RESPECTIVE FUNCTIONS EVALUATED AT X.                      
C                                                                       
C       LDFJAC IS A POSITIVE INTEGER INPUT PARAMETER NOT LESS THAN M    
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.      
C                                                                       
C       XP IS AN ARRAY OF LENGTH N. ON OUTPUT WHEN MODE = 1,            
C         XP IS SET TO A NEIGHBORING POINT OF X.                        
C                                                                       
C       FVECP IS AN ARRAY OF LENGTH M. ON INPUT WHEN MODE = 2,          
C         FVECP MUST CONTAIN THE FUNCTIONS EVALUATED AT XP.             
C                                                                       
C       MODE IS AN INTEGER INPUT VARIABLE SET TO 1 ON THE FIRST CALL    
C         AND 2 ON THE SECOND. OTHER VALUES OF MODE ARE EQUIVALENT      
C         TO MODE = 1.                                                  
C                                                                       
C       ERR IS AN ARRAY OF LENGTH M. ON OUTPUT WHEN MODE = 2,           
C         ERR CONTAINS MEASURES OF CORRECTNESS OF THE RESPECTIVE        
C         GRADIENTS. IF THERE IS NO SEVERE LOSS OF SIGNIFICANCE,        
C         THEN IF ERR(I) IS 1.0 THE I-TH GRADIENT IS CORRECT,           
C         WHILE IF ERR(I) IS 0.0 THE I-TH GRADIENT IS INCORRECT.        
C         FOR VALUES OF ERR BETWEEN 0.0 AND 1.0, THE CATEGORIZATION     
C         IS LESS CERTAIN. IN GENERAL, A VALUE OF ERR(I) GREATER        
C         THAN 0.5 INDICATES THAT THE I-TH GRADIENT IS PROBABLY         
C         CORRECT, WHILE A VALUE OF ERR(I) LESS THAN 0.5 INDICATES      
C         THAT THE I-TH GRADIENT IS PROBABLY INCORRECT.                 
C                                                                       
C     SUBPROGRAMS CALLED                                                
C                                                                       
C       MINPACK SUPPLIED ... DPMPAR                                     
C                                                                       
C       FORTRAN SUPPLIED ... DABS,DLOG10,DSQRT                          
C                                                                       
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.         
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             
c	Get EPSMCH from D1MACH instead of from DPMPAR.  Lars Stixrude, 2 September, 2018.
C                                                                       
C     **********                                                        
      INTEGER I,J                                                       
      DOUBLE PRECISION EPS,EPSF,EPSLOG,EPSMCH,FACTOR,ONE,TEMP,ZERO      
c      DOUBLE PRECISION DPMPAR
      DOUBLE PRECISION d1mach                                          
      DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/                          
C                                                                       
C     EPSMCH IS THE MACHINE PRECISION.                                  
C                                                                       
c      EPSMCH = DPMPAR(1)                                                
      EPSMCH = d1mach(1)                                                
C                                                                       
      EPS = DSQRT(EPSMCH)                                               
CLPS-->
	EPS = 1.e0*EPS
CLPS<--
C                                                                       
      IF (MODE .EQ. 2) GO TO 20                                         
C                                                                       
C        MODE = 1.                                                      
C                                                                       
         DO 10 J = 1, N                                                 
            TEMP = EPS*DABS(X(J))                                       
            IF (TEMP .EQ. ZERO) TEMP = EPS                              
            XP(J) = X(J) + TEMP                                         
   10       CONTINUE                                                    
         GO TO 70                                                       
   20 CONTINUE                                                          
C                                                                       
C        MODE = 2.                                                      
C                                                                       
         EPSF = FACTOR*EPSMCH                                           
         EPSLOG = DLOG10(EPS)                                           
         DO 30 I = 1, M                                                 
            ERR(I) = ZERO                                               
   30       CONTINUE                                                    
         DO 50 J = 1, N                                                 
            TEMP = DABS(X(J))                                           
            IF (TEMP .EQ. ZERO) TEMP = ONE                              
            DO 40 I = 1, M                                              
               ERR(I) = ERR(I) + TEMP*FJAC(I,J)                         
   40          CONTINUE                                                 
   50       CONTINUE                                                    
         DO 60 I = 1, M                                                 
            TEMP = ONE                                                  
            IF (FVEC(I) .NE. ZERO .AND. FVECP(I) .NE. ZERO              
     *          .AND. DABS(FVECP(I)-FVEC(I)) .GE. EPSF*DABS(FVEC(I)))   
     *         TEMP = EPS*DABS((FVECP(I)-FVEC(I))/EPS-ERR(I))           
     *                /(DABS(FVEC(I)) + DABS(FVECP(I)))                 
	dercalc = ERR(I)
            ERR(I) = ONE                                                
            IF (TEMP .GT. EPSMCH .AND. TEMP .LT. EPS)                   
     *         ERR(I) = (DLOG10(TEMP) - EPSLOG)/EPSLOG                  
            IF (TEMP .GE. EPS) ERR(I) = ZERO                            
            IF (TEMP .EQ. ONE) write(31,*) 'WARNING: chkder failed because change in f is too small'
     &                                     ,i,DABS(FVECP(I)-FVEC(I)),DABS(FVEC(I)),EPSF*DABS(FVEC(I))
            IF (ERR(I) .EQ. ZERO) write(31,*) 'WARNING: Error in derivative'
     &                                     ,i,DABS(FVECP(I)-FVEC(I))/EPS,dercalc,EPS

   60       CONTINUE                                                    
   70 CONTINUE                                                          
C                                                                       
      RETURN                                                            
C                                                                       
C     LAST CARD OF SUBROUTINE CHKDER.                                   
C                                                                       
      END                                                               
