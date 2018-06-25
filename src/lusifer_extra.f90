      subroutine lusifer_extra_max(maxex,maxgen) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      integer(kind=c_int) :: maxex,maxgen

      maxex = maxe
      maxgen = maxg
      end

      subroutine lusifer_extra_data(gen,nch) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer(kind=c_int) nchannel(maxg),nexternal(maxg),allbinary(maxg)
      common/lusifer_general/alphaisr,scale,meisr,s,p,mass,width, &
        nchannel,nexternal,allbinary
      bind(c) :: /lusifer_general/
      integer(kind=c_int), value :: gen
      integer(kind=c_int), intent(out) :: nch

      nch = nchannel(gen)
      end

      subroutine lusifer_extra_set(gen,nex,mw,gw,mz,gz,mh,gh,mt,gt) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer(kind=c_int) nchannel(maxg),nexternal(maxg),allbinary(maxg)
      common/lusifer_general/alphaisr,scale,meisr,s,p,mass,width, &
        nchannel,nexternal,allbinary
      bind(c) :: /lusifer_general/
      integer(kind=c_int), value :: gen,nex
      real(kind=c_double), value :: mw,gw,mz,gz,mh,gh,mt,gt

      nexternal(gen) = nex
      nchannel(gen) = 0
      mass(6) = mt
      width(6) = gt
      mass(23) = mz
      width(23) = gz
      mass(24) = mw
      width(24) = gw
      mass(25) = mh
      width(25) = gh
      end
