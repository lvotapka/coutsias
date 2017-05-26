c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c-----------------------<<< comp_eval.f       >>>-----------------------
c-----------------------------------------------------------------------
c....eigenvalue problem is solved in the form:
c              A * x = lambda * B * x
c....input
c      A:          matrix
c      B:          matrix
c....output
c      eval_beta(1):   Real part of largest e-val (i.e. "lambda")
c      eval_beta(2):   Imag part of largest e-val
c      i_max   :   index of largest real-part eigenpair.
c-----------------------------------------------------------------------
      subroutine comp_eval(e_m_max,m_0,m_f,A,B,step,window,matz,perturb,
     _        AT,BT,trans,revers,ddebug,evl_write,evlnam,evcnam,
     _        evlfil,evcfil,count,count_max,dist,
     _                  Z,alfr,alfi,beta,prevr,previ,eval_beta,tag,type)
      implicit none
      integer          e_m_max
      integer              tag(e_m_max),   type(e_m_max)
        integer ddebug
        integer matz
        integer evl_write
        double precision window
        integer m_0
        integer m_f
        integer step, count, count_max
        double precision perturb
        integer revers
        integer trans
        character*36 evlnam
        character*36 evcnam
      integer          i, ii, ierr, i_max
      integer          j
      integer          evlfil,   evcfil
      integer          m_eval_loop, m_trunc
      double precision  A(e_m_max,*),  B(e_m_max,*)
      double precision AT(e_m_max,*), BT(e_m_max,*)
      double precision Z(e_m_max,*)
      double precision  alfr(  e_m_max),   alfi(  e_m_max) 
      double precision  beta(  e_m_max),   dist(  e_m_max)
      double precision prevr(  e_m_max),  previ(  e_m_max)
      double precision eval_beta(  2)
      double precision A_norm, B_norm
      double precision dist1
c-----------------------------------------------------------------------
c....Initialize
c       Note:  eval_beta(1) keeps track of real part of largest e-val 
c              i_max keeps track of vector position of largest e-val
c-----------------------------------------------------------------------
      if (evl_write .eq. 1 .and. count .eq. 0) 
     _  open(unit=evlfil,file=evlnam,status='new')
c-----------------------------------------------------------------------
c....Main loop
c ____ Set up values; start computing over specified jumps (step)
c    from truncation m_0 up to truncation of order m_f
c-----------------------------------------------------------------------
      do 70 m_eval_loop = m_0, m_f, step
        m_trunc = m_eval_loop
c-----------------------------------------------------------------------
c........Initialize coefficient arrays
        call copym(e_m_max,e_m_max,m_trunc,m_trunc,A,AT)
        call copym(e_m_max,e_m_max,m_trunc,m_trunc,B,BT)
        call eval_matrix(e_m_max,m_trunc,AT,BT,perturb,revers,trans,
     _                   A_norm,B_norm,ddebug)
c-----------------------------------------------------------------------
c.......compute eigenvalues. 
        call rgg(e_m_max,m_trunc,AT,BT,alfr,alfi,beta,matz,Z,ierr)
        write(*,*) 'Actual Trunc = ',m_trunc,' Error code (1)  = ', ierr
c-----------------------------------------------------------------------
c.......convert eigenvalues to proper form
        do 10 j = 1, m_trunc
          tag(j)  = 0
          if(revers .eq. 1)
     _      beta(j) = (alfr(j)**2 + alfi(j)**2)/beta(j)
          alfr(j) = alfr(j)/beta(j)
          alfi(j) = alfi(j)/beta(j) 
 10     continue
c-----------------------------------------------------------------------
        eval_beta(1) =-999.d0
               i_max = 0
        do 40 i = 1, m_trunc
C---- Flag TYPE(1,2,3): alfi=><0 
          type(i) = 3
          if (alfi(i) .eq. 0) type(i) = 1
          if (alfi(i) .gt. 0) type(i) = 2
C---- Find stabilized eigenvalue with the largest real part
c    Note: the condition for alfr(i) .lt. 10 is arbitrary
          if ((eval_beta(1) .lt. alfr(i)) .and.
     _        (alfr(i) .lt. 10.d0)) then
            i_max = i
            eval_beta(1) = alfr(i)
            eval_beta(2) = alfi(i)
          end if
          if (m_trunc .gt. m_0) then
C---- find distance to closest eigenvalue from lower truncation set
            dist(i)  = 1953.d0
            do 30 ii = 1, m_trunc-step
              dist1 = sqrt((alfr(i)-prevr(ii))**2
     _                    +(alfi(i)-previ(ii))**2)
              if (dist1 .lt. dist(i)) dist(i) = dist1
 30         continue
C-- Flag TAG(1, 0)=(within,without) WINDOW to a lower trunc eval
            tag(i) = 0
            if (dist(i) .lt. window) tag(i) = 1
            if((ddebug .gt. 10).and.(tag(i) .eq. 1))
     _        write(*,5010) i, alfr(i), alfi(i), dist(i)
          endif
 40     continue
        if(m_trunc .gt. m_0)
     _    write(*,5010) i_max,alfr(i_max),alfi(i_max),dist(i_max)
        if(m_trunc .eq. m_0)
     _    write(*,5000) i_max,alfr(i_max),alfi(i_max)
c-----------------------------------------------------------------------
c........Save e-vals/vecs so can compare to next higher truncation.
        do 50 j = 1, m_trunc
          prevr(j) = alfr(j)
          previ(j) = alfi(j)
 50     continue 
c-----------------------------------------------------------------------
c....Write into e-val file the eigenvalues for this truncation (m_trunc)
        if (evl_write .ge. 1) then
          write(evlfil,5005) m_trunc
          do 60 i = 1, m_trunc
            write(evlfil,5000) tag(i), alfr(i), alfi(i)
 60       continue
        endif
c-----------------------------------------------------------------------
c....Write into e-vec file
c       if (matz .eq. 2 .and. count .eq. count_max) then
        if (matz .eq. 1) then
          if (i_max .eq. 0) return
          open(unit=evcfil,file=evcnam,status='new')
c         if (count .eq. 0) open(unit=evcfil,file=evcnam,status='new')
c         j = i_max
          write(evcfil,*)'count=',count,' trunc=',m_trunc
          do j = 1, m_trunc
c         write(evcfil,*) type(j), eval_beta(1), eval_beta(2)
          write(evcfil,*) type(j), alfr(j), alfi(j)
          if (type(j) .eq. 1) then
            write(evcfil,5015) (z(i,j), i = 1, m_trunc)
          else
            if(type(j) .eq. 2) j = j+1
          write(evcfil,*) type(j), alfr(j), alfi(j)
            write(evcfil,5016) (z(i,j-1),z(i,j), i = 1, m_trunc)
          endif
          end do
        end if
c-----------------------------------------------------------------------
c....End of Main loop
 70   continue 
c-----------------------------------------------------------------------
      if (count .eq. COUNT_MAX) then
        if (evl_write .eq. 1) close(unit=evlfil,status='keep')
        if (matz .eq. 1) close(unit=evcfil,status='keep')
c       if (matz .eq. 2) close(unit=evcfil,status='keep')
      endif
c-----------------------------------------------------------------------
 5000 format(i4,2(e22.14))
 5005 format(i6)
 5010 format(i4,2(e20.12),e13.5)
 5015 format(e22.14)
 5016 format(2(e22.14))
 9500 format(i3,1x,i3,1x,e22.14,1x,e22.14)
      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
