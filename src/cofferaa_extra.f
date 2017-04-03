      subroutine cofferaa_extra_max(maxex,maxgen)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      integer maxex,maxgen

      maxex = maxe
      maxgen = maxg
      end

      subroutine cofferaa_extra_data(gen,nch)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
      common/mcgenerator/mass,width,power,nchannel
      integer gen,nch

      nch = nchannel(gen)
      end

      subroutine cofferaa_extra_set(gen,mw,gw,mz,gz,mh,gh,mt,gt)
      implicit none
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
      common/mcgenerator/mass,width,power,nchannel
      integer gen
      real*8 mw,gw,mz,gz,mh,gh,mt,gt

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
