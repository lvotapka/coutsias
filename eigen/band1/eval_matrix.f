C-----------------------------------------------------------------------
c---------------------<<<    eval_matrix.f   >>>-----------------------
C-----------------------------------------------------------------------
c
c....Depending on inputs from data file 'eval_in':
c     (1) compute Frobenius norm of matrices
c     (2) perturb matrices
c     (3) reverse roles of matrices
c     (4) transpose matrices
c
c....Matrices A and B are input to routine; 
c    eigenvalue problem is solved in the form: 
c
c              A * x = lambda * B * x
c
c----------------------------------------------------------------------
       subroutine eval_matrix(max_m,m,A,B,perturb,revers,trans,
     _                        A_norm,B_norm,debug)

       implicit none

       integer  i, j 
       integer  max_m,m
       integer  debug
       integer  trans, revers

c      double precision A(512,512), B(512,512)
       double precision A(max_m,*), B(max_m,*)

       double precision A_norm
       double precision B_norm
       double precision temp
       double precision perturb

c----------------------------------------------------------------------
c(1)....Compute Frobenius norms of matrices A and B
c----------------------------------------------------------------------
       if (debug .gt. 5) then
          call norm(max_m,m,A,A_norm)
          call norm(max_m,m,B,B_norm)
c         write(*,*) ' A_norm = ', A_norm, ' B_norm = ', B_norm
       end if

c----------------------------------------------------------------------
c(2)....Perturb matrices to calculate pseudospectrum
c----------------------------------------------------------------------
       if (perturb .ne. 0.d0) then
         call perturb_mat(max_m,m,A,perturb)
       endif

c----------------------------------------------------------------------
c(3)....Reverse roles of A and B matrices
c----------------------------------------------------------------------
       if (revers .eq. 1) then
          do 20 j = 1, m
             do 10 i = 1, m
                temp   = A(i,j)
                A(i,j) = B(i,j)
                B(i,j) = temp
 10          continue
 20       continue
       endif

c----------------------------------------------------------------------
c(4)....Transpose matrices
c----------------------------------------------------------------------
       if (trans .eq. 1)
     _   then
          do 40 j = 1, m
             do 30 i = j+1, m
C ----    Transpose A
                temp   = A(j,i)
                A(j,i) = A(i,j)
                A(i,j) = temp
C ----    Transpose B
                temp   = B(j,i)
                B(j,i) = B(i,j)
                B(i,j) = temp
 30          continue
 40       continue
       endif

       end
c----------------------------------------------------------------------
      subroutine perturb_mat(max_m,m,A,perturb)
c----------------------------------------------------------------------

      implicit none

      integer max_m,m,i,j
      double precision A(max_m,*),perturb,random,seed
      external random

      seed = 0.314159265d0
      do 20 i=1, m
        do 10 j=1, m
          A(i,j) = (1.d0 + perturb*random(seed))*A(i,j)
 10     continue
 20   continue

      end
c----------------------------------------------------------------------
      double precision FUNCTION RANDOM(IX)
c----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CALL UNIRAN (IX, R)
      RANDOM = R
      END
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c----------------------------------------------------------------------
      SUBROUTINE UNIRAN (IX, R)
c----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C SUPPLIES A UNIFORMLY DISTRIBUTED RANDOM NUMBER  R  IN  (0,1)
      DATA M/2796203/, K/125/
C
      IC = IX * K
      IY = IC / M
      IX = K * IX - M * IY
      X = IX
      R = X / M
C
      END
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------<<<      norm.f      >>>-----------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE NORM(maxm,MM,A,ANORM)
       INTEGER maxm
       INTEGER  I, J, MM
       DOUBLE PRECISION A(maxm,mm)
       DOUBLE PRECISION TEMP, ANORM, BNORM
C______________________________________________________________________
C   Compute Frobenius norms of matrix A
C______________________________________________________________________
       ANORM = 0.D0
       DO 80 I = 1, MM
          DO 85 J = 1, MM
          ANORM = ANORM + A(I,J)*A(I,J)
 85       CONTINUE
 80    CONTINUE
       ANORM = SQRT(ANORM)
       END
C-----------------------------------------------------------------------
C                  ---------<block_tau_2>---------
C copy a block matrix (2X2, later kXk) to one at lower truncation,
C with tau-constraints at lowest 2 (m) rows in each block.
C to build:  
C       subroutine block_tau_km(e_m_max,m_trunc,m_f,k,m,A,B,AT,BT)
C       e_m_max: leading column dimension of matrices
C       m_trunc: order of truncation (overall) (order of AT,BT)
C       m_f    : actual order of matrices A,B; max truncation
C       k      : number of blocks
C       m      ; number of tau constraints at bottom of each block
C-----------------------------------------------------------------------
	subroutine block_tau_2(e_m_max,m_trunc,m_f,A,B,AT,BT)
        implicit none
        integer e_m_max, m_trunc,m_f
        integer i,j
        double precision  A(e_m_max,*),  B(e_m_max,*)
        double precision AT(e_m_max,*), BT(e_m_max,*)
        do i = 1, m_trunc/2-4
          do j = 1, m_trunc/2
            AT(i          ,j          ) = A(i      ,j      )
            BT(i          ,j          ) = B(i      ,j      )
            AT(i          ,j+m_trunc/2) = A(i      ,j+m_f/2)
            BT(i          ,j+m_trunc/2) = B(i      ,j+m_f/2)
            AT(i+m_trunc/2,j          ) = A(i+m_f/2,j      )
            BT(i+m_trunc/2,j          ) = B(i+m_f/2,j      )
            AT(i+m_trunc/2,j+m_trunc/2) = A(i+m_f/2,j+m_f/2)
            BT(i+m_trunc/2,j+m_trunc/2) = B(i+m_f/2,j+m_f/2)
          end do
        end do
        do i = 0,3
          do j = 1, m_trunc/2
            AT(m_trunc/2-i,j          ) = A(m_f/2-i,j      )
            BT(m_trunc/2-i,j          ) = B(m_f/2-i,j      )
            AT(m_trunc/2-i,j+m_trunc/2) = A(m_f/2-i,j+m_f/2)
            BT(m_trunc/2-i,j+m_trunc/2) = B(m_f/2-i,j+m_f/2)
            AT(m_trunc-i  ,j          ) = A(m_f-i  ,j      )
            BT(m_trunc-i  ,j          ) = B(m_f-i  ,j      )
            AT(m_trunc-i  ,j+m_trunc/2) = A(m_f-i  ,j+m_f/2)
            BT(m_trunc-i  ,j+m_trunc/2) = B(m_f-i  ,j+m_f/2)
          end do
        end do
        return
        end
C-----------------------------------------------------------------------
