      subroutine generatormax(maxex,maxgen)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=8,maxch=20000,maxg=1,maxv=40)
      integer maxex,maxgen

      maxex = maxe
      maxgen = maxg
      end

      subroutine generatordata(gen,nch)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=8,maxch=20000,maxg=1,maxv=40)
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
      common/general/alphaisr,scale,meisr,s,p,mass,width,nchannel,
     *  nexternal,allbinary
      integer gen,nch

      nch = nchannel(gen)
      end

      subroutine generatorset(gen,nex,mw,gw,mz,gz,mh,gh,mt,gt)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=8,maxch=20000,maxg=1,maxv=40)
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
      common/general/alphaisr,scale,meisr,s,p,mass,width,nchannel,
     *  nexternal,allbinary
      integer gen,nex
      real*8 mw,gw,mz,gz,mh,gh,mt,gt

      nexternal(gen) = nex
      nchannel(gen) = 0
      mass(6) = mt
      width(6) = gt
      mass(23) = mz
      width(23) = gz
      mass(24) = mh
      width(24) = gh
      mass(25) = mw
      width(25) = gw
      end
