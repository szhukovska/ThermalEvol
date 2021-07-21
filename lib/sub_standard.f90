
        Function TrapInt(n,x,y)
        USE nrtype; 
	implicit none
	integer :: n,i
	real(dp) :: x(n), y(n), TrapInt
	TrapInt=0d0
	do i = 2, n
	  TrapInt = TrapInt + (y(i)+y(i-1))*0.5d0 * (x(i)-x(i-1))
	enddo
	return
End Function TrapInt

      SUBROUTINE FindInArray(GRM, DM, IDIMG, IO, IL)
        USE nrtype; 
        implicit none
        integer :: IM
      	integer, intent(in) :: IDIMG
	real(dp), intent(in) :: DM, GRM(IDIMG)
	integer, intent(out) :: io, il

        IO=IDIMG
        IM=( IDIMG + 1 ) / 2
        IL=1

100    CONTINUE
       IF(GRM(IM).le.DM) THEN
         IL = IM
       ELSE
         IO = IM
       END IF
       IF(IO-IL.ge.2) THEN
         IM = ( IO + IL ) / 2
         goto 100
       END IF
      RETURN
      END
!******************************************************************C
!	1d interpolation
!******************************************************************C
      Function Interpol1D(n,x,y,x0)
       use nrtype
       implicit none
        integer :: n, j1,j2
        real(dp):: x(n),y(n),x0, Interpol1D
 	call FindInArray(x,x0,n,j1,j2)
	Interpol1D = y(j1) + (y(j2) - y(j1))/(x(j2)-x(j1))*(x0-x(j1))
      END
!*************************************************************************!
       Subroutine FillArrayEquid(nmax,amin, amax, a)
	implicit none
	integer :: i, nmax
	real(16) :: a(nmax), amin, amax, q

	q = (amax - amin) /(nmax-1d0)

	Do i = 1, nmax
	  a(i) = amin + q * (i-1d0)
	Enddo
End Subroutine FillArrayEquid
!*************************************************************************!
      Subroutine FillArrayEquidLog(nmax,amin, amax, a)
	implicit none
	integer :: i, nmax
	real(16) :: a(nmax), amin, amax, q

	q = (amax/amin) **(1d0/(nmax-1d0))

	Do i = 1, nmax
	  a(i) = amin * q ** (i-1d0)
	Enddo
      Return
    End Subroutine FillArrayEquidLog
!*************************************************************************!
!	Interpolation with extrapolation = simply assigning last value of beyond interval  
!*************************************************************************!
      Subroutine NewGridArrayExtra(n1,n2,x1,x2,y1,y2)
	implicit none
	integer, intent(in) :: n1, n2
	real(16), intent(in) :: x1(n1), x2(n2), y1(n1)
	real(16), intent(out) :: y2(n2)
	integer :: i, j1, j2, max,min

        call FindInArray(x2,x1(1),n2,min,j1)
       call FindInArray(x2,x1(n1),n2,max,j1)
        y2(1:min)=0d0 
        y2(max:n2) = y1(n1)
	 Do i= min, max
	    call FindInArray(x1,x2(i),n1,j2,j1)
	    Y2(i)=y1(j1) + (y1(j2)-y1(j1))/(x1(j2)-x1(j1)) * (x2(i) - x1(j1))
	 End Do
      Return
      End

!*************************************************************************!
! Sorts an array a(1:n) into ascending numerical order by Shell’s
! method (diminishing increment sort). 
! n is input; a is replaced on output by its sorted rearrangement.
!*************************************************************************!
	SUBROUTINE SortArrayShell(n, a)
        use nrtype
	INTEGER :: n
	REAL(dp) :: a(n)
	INTEGER :: i,j,inc
	REAL(dp) :: v
	inc=1 !Determine the starting increment.
1 	inc=3*inc+1
	if(inc.le.n)goto 1
2 	continue !Loop over the partial sorts.
	inc=inc/3
	do i=inc+1, n !Outer loop of straight insertion.
	   v=a(i)
	   j=i
3 	   if(a(j-inc).gt.v)then !Inner loop of straight insertion.
	     a(j)=a(j-inc)
	     j=j-inc
             if(j.le.inc)goto 4
	     goto 3
	   endif
4 	   a(j)=v
	end do 
	if(inc.gt.1)goto 2
        END SUBROUTINE SortArrayShell
!*************************************************************************!
!	Sort Grid Points in array "a", "b" - array with function value on this grid
!*************************************************************************!
      SUBROUTINE Sort2ArrayShell(n, a, b)
        use nrtype
	INTEGER :: n
	REAL(dp) :: a(n), b(n)
	INTEGER :: i, j, inc
	REAL(dp):: v, vv
	inc = 1 !Determine the starting increment.
1 	inc = 3*inc + 1
	if(inc.le.n) goto 1
2 	continue !Loop over the partial sorts.
	inc=inc/3
	do i=inc+1, n !Outer loop of straight insertion.
	   v=a(i)
	   vv=b(i)
	   j=i
3 	   if(a(j-inc).gt.v)then !Inner loop of straight insertion.
	     a(j)=a(j-inc)
	     b(j)=b(j-inc)
	     j=j-inc
             if(j.le.inc)goto 4
	     goto 3
	   endif
4 	   a(j)=v
	   b(j)=vv
	end do 
	if(inc.gt.1)goto 2
       END SUBROUTINE Sort2ArrayShell

    FUNCTION zbrent_mod(func,x1,x2,tol,interval)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x1,x2,tol
	REAL(dp) :: zbrent_mod
        logical :: interval
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(DP), PARAMETER :: EPS=epsilon(x1)
	INTEGER(I4B) :: iter
	REAL(DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)
!       print*,'*** Start root finding in zbrent in the interval:'
!       print'(A,1p2e11.3)','a, fa', a, fa
!       print'(A,1p2e11.3)','b, fb', b, fb
 
	if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
		print*,'root must be bracketed for zbrent'
                interval = .false.
                zbrent_mod = 0d0
                return
        end if
        interval = .true.        
	c=b
	fc=fb
	do iter=1,ITMAX

		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
			fc=fa
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol1=2.0_dp*EPS*abs(b)+0.5_dp*tol
		xm=0.5_dp*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			zbrent_mod=b
                        write(*,'(A, 1p2e14.4, A, i2, A)') 'Found root x, f(x): ',b, fb, ' in ', iter, ' iterations'
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_dp*xm*s
				q=1.0_dp-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
				q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_dp*p  <  min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
		fb=func(b)
!!           write(*,'(A, i2,A, 1p2e14.4,A, 1pe14.4)') 'Current iteration i:', iter, '   x_i, f_i:',b, fb ,  '   error:',  abs(d)
	end do
	call nrerror('zbrent: exceeded maximum iterations')
	zbrent_mod=b
	END FUNCTION zbrent_mod
!------------------------------------

!*************************************************************************!
	Subroutine ShiftReadFile (nfile,nshift)
	implicit none
	integer nfile, nshift, i
	 do i=1,nshift
	  read(nfile,*,err=10)
	 enddo
	 return
10	write(5,*)'Error in ShiftReadFile, unit=',nfile
	stop
	End

