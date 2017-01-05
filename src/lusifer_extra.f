      subroutine lusifer_extra_max(maxex,maxgen)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      integer maxex,maxgen

      maxex = maxe
      maxgen = maxg
      end

      subroutine lusifer_extra_data(gen,nch)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
      common/lusifer_general/alphaisr,scale,meisr,s,p,mass,width,
     *  nchannel,nexternal,allbinary
      integer gen,nch

      nch = nchannel(gen)
      end

      subroutine lusifer_extra_set(gen,nex,mw,gw,mz,gz,mh,gh,mt,gt)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
      common/lusifer_general/alphaisr,scale,meisr,s,p,mass,width,
     *  nchannel,nexternal,allbinary
      integer gen,nex
      real*8 mw,gw,mz,gz,mh,gh,mt,gt

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
