C-----------------------------------------------------------------------
         program jac
c        implements and test a general jacobi transform
C-----------------------------------------------------------------------
         implicit none
      integer     m_max, n_max
        parameter(m_max = 512, n_max = 0)
      DOUBLE PRECISION JACM(0:M_MAX,0:M_MAX)
      DOUBLE PRECISION   AA(0:M_MAX,0:M_MAX)
      DOUBLE PRECISION CHEBT(0:M_MAX,0:M_MAX)
      DOUBLE PRECISION T(0:M_MAX)
      DOUBLE PRECISION W(0:M_MAX)
      DOUBLE PRECISION G(0:M_MAX)
      DOUBLE PRECISION b(0:M_MAX)
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926535897932384626D0)
      double precision DATA(0:M_MAX,0:N_MAX+1), TRANSFORM(0:M_MAX)
      double precision DAT2(0:M_MAX)
      double precision al, be, x, nu, endpts(0:1)
      integer M, N, I, J, kk, kpts, kind, k
C-----------------------------------------------------------------------
        open(unit=13,file='in',status='old',form='formatted')
          read(13,*) x
          read(13,*) m
          read(13,*) al
          read(13,*) be
          read(13,*) nu
          read(13,*) kind
          read(13,*) kk
          read(13,*) kpts
          read(13,*) endpts(0)
          read(13,*) endpts(1)
        close(13)
c       write(*,*) 'input M = ? '
c       READ(*,*) M
        do i = 0,M
           do j = 0,N+1
              DATA(I,J) = 0.d0
              DAT2(I) = 0.d0
           end do
c          DATA(I,0) = 1.d0
c          DATA(I,0) = cos(PI*I/M)
           DATA(I,0) = 2.d0*cos(PI*I/M)**2-1.d0
           DAT2(I) = DATA(I,0)
        end do
c-----------------------------------------------------------------------
        CALL JACTINI(M_MAX,M, AL, BE, kind, kpts,endpts,JACM,T,W,G,b)
c       write(*,*) 'T'
c       write(*,*) (T(i), i= 0,m)
c       write(*,*) 'W'
c       write(*,*) (W(i), i= 0,m)
c       write(*,*) 'G'
c       write(*,*) (G(i), i= 0,m)
c       write(*,*) 
c       call mug(JACM,M_MAX,0,8,0,8)
c       write(*,*) 
C--- Test orthonormality
      Do i = 0,m
        Do j = 0,m
            AA(i,j) = 0.d0
          Do k= 0,m
            AA(i,j) = AA(i,j)+JACM(i,k)*G(j)*JACM(j,k)*W(k)
          end do
        end do
      end do
c     call mug(AA,M_MAX,0,8,0,8)
        write(*,*) 
        write(*,*) 
c       do i = 0,M
c          do j = 0,N+1
c             DATA(I,J) = 0.d0
c             DAT2(I) = 0.d0
c          end do
c          DATA(I,0) =  1.d0
c          DATA(I,0) = T(i)**3+T(i)
c          DATA(I,0) = 2.d0*cos(PI*I/M)**2-1.d0
c          DAT2(I) = DATA(I,0)
c       end do
c       write(*,*) 'data'
c       write(*,*) (DATA(I,0), I=0,M)
        CALL JACT(M_MAX,M,DAT2,TRANSFORM, 1,JACM,W,G,b)
        write(*,*) 'data'
        write(*,*) (DAT2(I), I=0,M)
        write(*,*) 'transform'
        write(*,*) (TRANSFORM(I), I=0,M)
        CALL JACT(M_MAX,M,DAT2,TRANSFORM,-1,JACM,W,G,b)
        write(*,*) 
        write(*,*) ' data, after xform (should be same)'
        write(*,*) 
        write(*,*) (DATA(I,0)-DAT2(I), I=0,M)
        write(*,*) (DAT2(I), I=0,M)
        write(*,*) 
        write(*,*)  ' original data'
        write(*,*) 
        write(*,*) (DATA(I,0), I=0,M)
C---- loop for fft
        write(*,*) 
        write(*,*) '---------------------------------' 
        write(*,*) 
        call cosftini(M_MAX,M,CHEBT)
c       call mug(CHEBT,M_MAX,0,8,0,8)
      CALL COST (M_MAX,M,DAT2,TRANSFORM, 1, CHEBT)
        write(*,*) ' Cheb. transform '
        write(*,*) (TRANSFORM(I), I=0,M)
      CALL COST (M_MAX,M,DAT2,TRANSFORM,-1, CHEBT)
        write(*,*) 
        write(*,*) ' differences CHEBT'
        write(*,*) 
c       write(*,*) (DATA(I,0)-DAT2(I), I=0,M)
        write(*,*) (DATA(I,0), I=0,M)
        write(*,*) 
        write(*,*)  'data'
        write(*,*) 
        write(*,*) (DAT2(I), I=0,M)
        END
C------------------------------------------------------
