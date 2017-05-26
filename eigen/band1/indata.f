C----------------------------------------------------------------------|
      subroutine lin_indata(lu_in)
C----------------------------------------------------------------------|
      implicit none
        include 'limits.h'
        include 'indata.h'
        integer lu_in
C----------------------------------------------------------------------|
	open(unit=lu_in,file='in',status='old',form='formatted')
C----------------------------------------------------------------------|
c...read physical variables
C----------------------------------------------------------------------|
        read(lu_in,*) oper_case
        read(lu_in,*) nu
        read(lu_in,*) al
        read(lu_in,*) be
        read(lu_in,*) bc_case
        read(lu_in,*) mode
C----------------------------------------------------------------------|
        read(lu_in,*) debug
C----------------------------------------------------------------------|
        read(lu_in,*) m_0
        read(lu_in,*) m_f
        read(lu_in,*) step
        read(lu_in,*) window
        read(lu_in,*) perturb
        read(lu_in,*) matz
        read(lu_in,*) revers
        read(lu_in,*) trans
        read(lu_in,*) alpha_tol
        read(lu_in,*) count_max
        read(lu_in,*) evl_write
        read(lu_in,*) evlnam
        read(lu_in,*) evcnam
        read(lu_in,*) ddebug
        close(lu_in)
        write(*,*) 'op_case = ', oper_case
        write(*,*) 'nu = ',nu,' al,be = ',al,be
c       write(*,*) bc_case
c       write(*,*) mode
c       write(*,*) debug
c       write(*,*) m_0
c       write(*,*) m_f
c       write(*,*) step
c       write(*,*) window
c       write(*,*) perturb
c       write(*,*) matz
c       write(*,*) revers
c       write(*,*) trans
c       write(*,*) alpha_tol
c       write(*,*) evl_write
c       write(*,*) evlnam
c       write(*,*) evcnam
c       write(*,*) ddebug
c       stop
C----------------------------------------------------------------------|
c... check input data
C----------------------------------------------------------------------|
c.. calculate pi=3.1415926535897932384626.......
C----------------------------------------------------------------------|
        pi = 4.0d0*atan(1.0d0)
        write(*,*) pi*pi, pi*pi/4.d0
C----------------------------------------------------------------------|
        end
C----------------------------------------------------------------------|
        subroutine open_file(unit_file,file_file,status_file,
     #                       access_file,form_file)
C----------------------------------------------------------------------|
        implicit none
        integer       unit_file
        character*(*)  file_file
        character*(*)  status_file
        character*(*)  access_file
        character*(*)  form_file
c       open(unit=unit_file,file=file_file,status=status_file,
c    #       position=access_file,form=form_file)
        open(unit=unit_file,file=file_file,status=status_file,
     #       access=access_file,form=form_file)
        end
C----------------------------------------------------------------------|
