      subroutine dladgen(a,b,freq)
c
      real*8 a,b,freq(13),t
      real*8 dcum
c
      do i1=1,8
	 t=i1*10. 
	 freq(i1)=dcum(a,b,t)
      end do
c
      do i2=9,12
	 t=80.+(i2-8)*2.
	 freq(i2)=dcum(a,b,t)
      end do
c
      freq(13)=1.
c
      do i=13,2,-1
	 freq(i)=freq(i)-freq(i-1)
      end do
c
      return
      end


      function dcum(a,b,t)
c
	real*8 a,b,t,pi,rd,eps,delx,p,y,x,dx,dcum
      pi=atan(1.)*4.
      rd=pi/180.

	if (a.gt.1.) then
c
		dcum=1.-cos(rd*t)
c
	else

		eps=1.e-8
		delx=1.
c
		x=2*rd*t
		p=x
c
		do while (delx.gt.eps) 
			y = a*sin(x)+.5*b*sin(2.*x)
			dx=.5*(y-x+p) 
			x=x+dx
			delx=abs(dx)
		end do
c
		dcum=(2.*y+p)/pi
c
	endif
c
      return
      end
