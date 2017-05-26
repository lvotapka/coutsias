C-----------------------------------------------------------------------
         program test
c        implements and tests a general classical orthopoly transform
C-----------------------------------------------------------------------
         implicit none
      integer     m_max, n_max, buff, q_max
        parameter(m_max = 300, n_max = 0, buff = 12, q_max = 7)
      DOUBLE PRECISION JACM(0:M_MAX,0:M_MAX)
      DOUBLE PRECISION   AA(0:M_MAX,0:M_MAX)
        double precision     H5(0:m_max   ,  0:1       , 0:n_max/2)
        double precision      U(0:m_max   ,  0:1       , 0:n_max/2)
        double precision     U0(0:m_max   ,  0:1       , 0:n_max/2)
        double precision     O5(0:m_max+buff,  1:5       , 0:n_max/2)
        DOUBLE PRECISION     A1(0:M_MAX+buff,1:2*Q_MAX+1)
        DOUBLE PRECISION     A2(0:M_MAX+buff,1:2*Q_MAX+1)
        double precision     C5(0:m_max     ,  1:2       , 0:n_max/2)
        double precision     V5(1:2         ,  0:1       , 0:n_max/2)
        double precision    CH5(1:2         ,  0:1       , 0:n_max/2)
        double precision  alpha(0:2         ,  0:2                  )
      DOUBLE PRECISION T(0:M_MAX)
      DOUBLE PRECISION W(0:M_MAX)
      DOUBLE PRECISION G(0:M_MAX)
      DOUBLE PRECISION b(0:M_MAX)
      DOUBLE PRECISION AM(0:M_MAX+buff,-1:1)
      DOUBLE PRECISION BM(0:M_MAX+buff,-1:1)
        double precision   bb(0:m_max)
        double precision   cc(0:m_max)
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926535897932384626D0)
      double precision DATA(0:M_MAX,0:N_MAX+1), TRANSFORM(0:M_MAX)
      double precision DAT2(0:M_MAX)
      double precision al, be, x, endpts(0:1), parm(32)
      integer M, N, I, J, kk, kpts, kind, k, parno, debug, ii, problem
C-----------------------------------------------------------------------
        open(unit=13,file='in',status='old',form='formatted')
          read(13,*) x
          read(13,*) m
          read(13,*) al
          read(13,*) be
          read(13,*) kind
          read(13,*) kk
          read(13,*) kpts
          read(13,*) endpts(0)
          read(13,*) endpts(1)
          read(13,*) parno
          do i = 1, parno
            read(13,*) parm(i)
          end do
          read(13,*) debug
          read(13,*) problem
        close(13)
        write(*,*) 'kind = ', kind,m,al,be, kpts, endpts(0), endpts(1)
C-----------------------------------------------------------------------
C--- construction of X, B matrices: method specific.
        call rec_X(m_max+buff,AM,BM,kind,al,be)
C--- Initialize transform matrices: method specific.
        CALL JACTINI(M_MAX,M, AL, BE, kind, kpts,endpts,JACM,T,W,G,b)
C-----------------------------------------------------------------------
C -- Construct expansions of various functions by recurence/DE method
        if (problem .eq. 6) then
          do i = 0, m
            if (T(i) .le. parm(1)) then 
              DAT2(i) = 0.d0
            else
              DAT2(i) = 1.d0
            endif
          end do
          CALL JACT(M_MAX,M,DAT2,TRANSFORM, 1,JACM,W,G,b)
          CALL D1MAT(m_max, AA)
          do i = 0,m
            bb(i) = 0.d0
            do j = 0,m
              bb(i) = bb(i)+TRANSFORM(j)*AA(i,j)              
            end do
          end do
          CALL JACT(M_MAX,M,DAT2,bb,-1,JACM,W,G,b)
          do i = 0,m
            write(*,*) sngl(T(i)), sngl(dsqrt(1.d0-T(i)*T(i))*DAT2(i))
          end do
          STOP
        endif
        if (problem .eq. 0) then
          call cond1(m_max,1,1,A1,1.d0,debug,buff,AM,BM)
          call mug(A1,m_max+buff,0,6,0,8)
        STOP
        endif
        if (problem .eq. 1) then
        call cond1(m_max,1,1,A1,-parm(1),debug,buff,AM,BM)
        call cond1(m_max,0,0,A2,0,debug,buff,AM,BM)
        A2(0,0) =0.d0
        call bbadd(m_max,m,-2,0,-2,2,0,2,A1,A2,A1,1.d0,parm(2)/2,buff)
c       call cond1(m_max,1,0,A1,0.d0,0,buff,AM,BM)
c       call cond1(m_max,0,0,A2,0.d0,0,buff,AM,BM)
c       call bbadd(m_max,m,-1,0,-1,1,0,1,A1,A2,A1,parm(1),1.d0,buff)
C -- Define boundary conditions and rhs
        call jac_q(kind,0.d0,al,be,0,m,C5(0,1,0))
      V5(1,0,0) = parm(3)*dexp(-parm(1)**2/parm(2))
c     V5(1,0,0) = parm(2)
      call band_op(m_max,m,2,0,A1,O5,1,C5(0,1,0),V5,H5,CH5,debug,buff,
     _             bb,cc,U)
      call ul_sol(m_max,m,2,0,A1,O5,1,C5(0,1,0),V5,H5,CH5,U,0,0,1,debug,
     _             buff,U0,bb,cc)
c     call band_op(m_max,m,1,0,A1,O5,1,C5(0,1,0),V5,H5,CH5,0,buff,
c    _             bb,cc,U)
c     call  ul_sol(m_max,m,1,0,A1,O5,1,C5(0,1,0),V5,H5,CH5,U,0,0,1,0,
c    _             buff,U0,bb,cc)
      do i = 0,m
        write(*,*) i, U(i,0,0)
      end do
      ii = int(parm(4))
      do i = 0, 300
        x = dble(i)/60
        call jac_q(kind,x,al,be,0,ii, transform)
          DAT2(i) = 0.d0
        do j = 0,ii
          DAT2(i) = DAT2(i)+transform(j)*U(j,0,0)*dexp(-x)
        end do 
        write(*,*) sngl(x), sngl(DAT2(i))
      end do
      STOP
      end if
c-----------------------------------------------------------------------
C compare expansion found by integral operator inversion with direct
C construction from discrete Laguerre transform
      do i = ii+1,m
        U(i,0,0) = 0.d0
      end do
        CALL JACT(M_MAX,M,DAT2,U,-1,JACM,W,G,b)
        write(*,*) 'data from DE: pointvalues'
c     do i = 0,m
c       write(*,*) T(i), DAT2(i)*dexp(-T(i))
c     end do
        DO I = 0, int(parm(4))
          write(*,*) sngl(T(I)), sngl(DAT2(I)*dexp(-T(i))) 
        END DO
        write(*,*) 'direct evaluation and transform'
C--- Data
        do i = 0,M
           do j = 0,N+1
              DATA(I,J) = 0.d0
              DAT2(I) = 0.d0
           end do
C -- DATA FOR LAGUERRE AND HERMITE
        DATA(I,0) = parm(3)*dexp(-(t(i)-parm(1))**2/parm(2))
c          DATA(I,0) = dexp(-t(i)*parm(1))*parm(2)
           DAT2(I) = DATA(I,0)
           DATA(I,0) = DATA(I,0)*DEXP(-T(I))
        end do
        write(*,*) 'data'
        DO I = 0, M
          write(*,*) sngl(T(I)), sngl(DATA(I,0)) 
        END DO
c-----------------------------------------------------------------------
        CALL JACT(M_MAX,M,DAT2,TRANSFORM, 1,JACM,W,G,b)
c       write(*,*) ' data, original '
c       DO I = 0, M
c         write(*,*) sngl(T(I)), sngl(DAT2(I)) 
c       END DO
        write(*,*) 'transform'
        DO I = 0, M
          write(*,*) I, TRANSFORM(I) 
        END DO
        CALL JACT(M_MAX,M,DAT2,TRANSFORM,-1,JACM,W,G,b)
        write(*,*) 
        write(*,*) ' data-data_original, after xform (should be same)'
        write(*,*) 
        write(*,*) (DATA(I,0)-DAT2(I)*DEXP(-T(i)), I=0,M)
c       DO I = 0, M
c         write(*,*) sngl(T(I)), sngl(DAT2(I)) 
c       END DO
        write(*,*) 
        END
C------------------------------------------------------
