!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     Monte Carlo generator                                       c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_generation(random,kbeam,k,g,channel, &
        generator,switch)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      real*8 random(3*maxe-10),p(0:3,2**maxe),s(2**maxe)
      real*8 k(maxe,0:3),kt(maxe,0:3),kbeam(2,0:3),cofferaa_h
      real*8 mmin,mmax,smin,smax,x,g,tcut1,tcut2
      integer i1,i2,i3,i4,ns,nt,channel,generator,switch,step
      integer ranstart,virtinv,inv1,inv2,em,sp,ga,cofferaa_pid
! mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
! mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
! mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
! mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
! mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
! mcoutput
      integer nout,numout,maxout
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv, &
        lmin,lmax
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process, &
        in2process,out1process,out2process,inprocess,virtprocess, &
        idhepprocess,nprocess
      common/mcoutput/nout,numout,maxout
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      if(switch.eq.0)return
!      step=2 ! calculation of local density
      step=1 ! calculation of momenta
      em=emmap(channel,generator)
      sp=spmap(channel,generator)
      ga=gamap(channel,generator)
      if(step.eq.1)then
! incoming momenta and definition of x
        ranstart=1+ninv(channel,generator) &
          +2*nprocess(channel,generator)+2*ndecay(channel,generator)
        if(em.eq.1.or.(sp.eq.1.and.em.ge.3))then
          x=1d0-cofferaa_h(random(ranstart),powermap,0d0,0d0,1d0,switch)
          do i1=0,3
            p(i1,1)=-x*kbeam(1,i1)
            p(i1,2)=-kbeam(2,i1)
          enddo
        elseif(em.eq.2.or.(sp.eq.2.and.em.ge.3))then
          x=1d0-cofferaa_h(random(ranstart),powermap,0d0,0d0,1d0,switch)
          do i1=0,3
            p(i1,1)=-kbeam(1,i1)
            p(i1,2)=-x*kbeam(2,i1)
          enddo
        else
          x=1d0
          do i1=0,3
            p(i1,1)=-kbeam(1,i1)
            p(i1,2)=-kbeam(2,i1)
          enddo
        endif
! remaining momenta
        do i1=0,3
          p(i1,3)=p(i1,1)+p(i1,2)
          p(i1,allbinary(ga,generator)-1)=-p(i1,1)
          p(i1,allbinary(ga,generator)-2)=-p(i1,2)
          p(i1,allbinary(ga,generator)-3)=-p(i1,3)
        enddo
! square of center-of-mass energy
        s(3)=p(0,3)**2-p(1,3)**2-p(2,3)**2-p(3,3)**2
        s(allbinary(ga,generator)-3)=s(3)
! event excluded by cuts after intial-state radiation
        if(s(3).le.mcutinv(allbinary(ga,generator)-3,generator)**2)then
          if(numout.lt.maxout)then
            write(nout,'(a)')' generation: s < smin'
            numout=numout+1
          endif
          switch=0
          return
        endif
      else
! calculation of local density
        g=1d0
! defining momenta
        do i1=1,nexternal(generator)
        do i2=0,3
          kt(i1,i2)=k(i1,i2)
        enddo
        enddo
        call cofferaa_mapin(kt,g,em,sp,ga,nexternal(generator),switch)
        do i1=1,allbinary(ga,generator)
          p(0,i1)=0d0
        enddo
        do i1=0,3
          p(i1,1)=-kt(1,i1)
          p(i1,2)=-kt(2,i1)
          do i2=3,nexternal(generator)
            p(i1,2**(i2-1))=kt(i2,i1)
          enddo
        enddo
! calculating invariants
        do i1=1,allbinary(ga,generator)-1
          if(p(0,allbinary(ga,generator)-i1).eq.0d0)then
            do i2=1,nexternal(generator)
              i3=i1/2**(i2-1)
              if(2*(i3/2).ne.i3.and.i1.ne.2**(i2-1))then
                do i4=0,3
                  p(i4,i1)=p(i4,i1-2**(i2-1))+p(i4,2**(i2-1))
                enddo
                goto 100
              endif
            enddo
 100        s(i1)=p(0,i1)**2-p(1,i1)**2-p(2,i1)**2-p(3,i1)**2
          else
            do i2=0,3
              p(i2,i1)=-p(i2,allbinary(ga,generator)-i1)
            enddo
            s(i1)=s(allbinary(ga,generator)-i1)
          endif
        enddo
      endif
! external masses
      do i1=1,nexternal(generator)
      if(i1.ne.ga)then
        s(2**(i1-1))=massext2(i1,generator)
        s(allbinary(ga,generator)-2**(i1-1))=massext2(i1,generator)
      endif
      enddo
! inv
      do ns=1,ninv(channel,generator)
        inv1=ininv(ns,channel,generator)
        inv2=allbinary(ga,generator)-3-ininv(ns,channel,generator)
        mmin=mcutinv(inv1,generator)
        mmax=mcutinv(inv2,generator)
        do i1=1,ns-1
          virtinv=ininv(i1,channel,generator)
          if(lmin(i1,ns,channel,generator).and. &
            dsqrt(s(virtinv)).gt.mcutinv(virtinv,generator))then
            mmin=mmin+dsqrt(s(virtinv))-mcutinv(inv1,generator) &
              +mcutinv(inv1-virtinv,generator)
            inv1=inv1-virtinv
          endif
          if(lmax(i1,ns,channel,generator).and. &
            dsqrt(s(virtinv)).gt.mcutinv(virtinv,generator))then
            mmax=mmax+dsqrt(s(virtinv))-mcutinv(inv2,generator) &
              +mcutinv(inv2-virtinv,generator)
            inv2=inv2-virtinv
          endif
        enddo
        smin=mmin**2
        smax=(dsqrt(s(3))-mmax)**2
        call cofferaa_inv( &
          random(ns), &                                  ! random number
          s(ininv(ns,channel,generator)), &              ! invariant mass
          g, &                                           ! local density
          mass(idhepinv(ns,channel,generator)), &        ! mass
          width(idhepinv(ns,channel,generator)), &       ! width
          powerinv(ns,channel,generator), &              ! mapping
          smin, &                                        ! lower bound
          smax, &                                        ! upper bound
          step,switch)
      enddo
! process
      ranstart=ninv(channel,generator)
      do nt=1,nprocess(channel,generator)
        tcut1=0d0
        tcut2=0d0
        if(nt.eq.1)then
          tcut1=tcutprocess(allbinary(ga,generator) &
            -virtprocess(nt,channel,generator),generator)
          tcut2=tcutprocess(virtprocess(nt,channel,generator), &
            generator)
        endif
        call cofferaa_process( &
          random(ranstart+2*nt-1), &                     ! random numbers
          p(0,in1process(nt,channel,generator)), &       ! incomming particle 1
          p(0,in2process(nt,channel,generator)), &       ! incomming particle 2
          p(0,out1process(nt,channel,generator)), &      ! outgoing particle 1
          p(0,out2process(nt,channel,generator)), &      ! outgoing particle 2
          p(0,virtprocess(nt,channel,generator)), &      ! virtual particle
          g, &                                           ! local density
          mass(idhepprocess(nt,channel,generator)), &    ! mass
          width(idhepprocess(nt,channel,generator)), &   ! width
          powerprocess(nt,channel,generator), &          ! mapping
          s(inprocess(nt,channel,generator)), &          ! incomming particles
          s(out1process(nt,channel,generator)), &        ! outgoing particle 1
          s(out2process(nt,channel,generator)), &        ! outgoing partcile 2
          s(virtprocess(nt,channel,generator)), &        ! virtual particle
          s(in1process(nt,channel,generator)), &         ! incomming particle 1
          s(in2process(nt,channel,generator)), &         ! incomming particle 2
          tcut1,tcut2, &                                 ! angular cuts
          step,switch)
      enddo
! decay
      ranstart=ninv(channel,generator)+2*nprocess(channel,generator)
      do ns=1,ndecay(channel,generator)
        call cofferaa_decay( &
          random(ranstart+2*ns-1), &                     ! random numbers
          p(0,indecay(ns,channel,generator)), &          ! incomming particle
          p(0,out1decay(ns,channel,generator)), &        ! outgoing particle 1
          p(0,out2decay(ns,channel,generator)), &        ! outgoing particle 2
          g, &                                           ! local density
          s(indecay(ns,channel,generator)), &            ! incomming particle
          s(out1decay(ns,channel,generator)), &          ! outgoing particle 1
          s(out2decay(ns,channel,generator)), &          ! outgoing particle 2
          step,switch)
      enddo
! defining external momenta
      do i1=0,3
        k(1,i1)=-p(i1,1)
        k(2,i1)=-p(i1,2)
        do i2=3,nexternal(generator)
          k(i2,i1)=p(i1,2**(i2-1))
        enddo
      enddo
! map
      if(step.eq.1)then
        ranstart=1+ninv(channel,generator) &
          +2*nprocess(channel,generator)+2*ndecay(channel,generator)
        call cofferaa_mapout(random(ranstart),k,x,g,em,sp,ga, &
          nexternal(generator),switch)
      endif
! restoring external masses
      do i1=1,nexternal(generator)
        k(i1,0)=dsqrt(massext2(i1,generator) &
          +k(i1,1)**2+k(i1,2)**2+k(i1,3)**2)
      enddo
!      call cofferaa_checkmom(k,nexternal(generator),switch)
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     calculating densities                                       c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_checkmom(k,next,switch)
      implicit none
! local variables
      integer maxe,maxg
      parameter(maxe=9,maxg=4)
      real*8 k(maxe,0:3),mom,acc
      integer i1,i2,next,switch
! mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
      common/mckinematics/massext2,allbinary,nexternal
      parameter(acc=1d-10)
      if(switch.eq.0)return
      do i1=0,3
        mom=k(1,i1)+k(2,i1)
        do i2=3,next
          mom=mom-k(i2,i1)
        enddo
        if(dabs(mom).gt.acc)then
          write(*,'(" momentum conservation violated:",i2,1d24.16)') &
            i1,mom
        endif
      enddo
      do i1=1,next
        mom=k(i1,0)**2
        do i2=1,3
          mom=mom-k(i1,i2)**2
        enddo
        if(dabs(mom).gt.acc)then
          write(*,'(" onshellness violated:",i2,1d24.16)')i1,mom
        endif
      enddo
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     calculating denities                                        c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_density(k,g,generator,switch)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      real*8 random(3*maxe-7),p(0:3,2**maxe,0:maxe,0:maxe,0:maxe)
      real*8 s(2**maxe,0:maxe,0:maxe,0:maxe),k(maxe,0:3),kt(maxe,0:3)
      real*8 gmap(0:maxe,0:maxe,0:maxe),ginv(maxch),gprocess(maxch)
      real*8 gdecay(maxch),g(maxch),mmin,mmax,smin,smax,tcut1,tcut2
      integer i1,i2,i3,i4,ns,nt,channel,generator,switch,step,pid
      integer virtinv,inv1,inv2,em,sp,ga
! mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
! mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
! mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
! mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
! mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
! mcdensity
      integer nsinv(maxch,maxg),chinv(maxch,maxg),maxinv(maxg)
      integer ntprocess(maxch,maxg),chprocess(maxch,maxg)
      integer maxprocess(maxg),nsdecay(maxch,maxg),chdecay(maxch,maxg)
      integer maxdecay(maxg),numinv(maxe,maxch,maxg)
      integer numprocess(maxe,maxch,maxg),numdecay(maxe,maxch,maxg)
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv, &
        lmin,lmax
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process, &
        in2process,out1process,out2process,inprocess,virtprocess, &
        idhepprocess,nprocess
      common/mcdensity/nsinv,chinv,maxinv,ntprocess,chprocess, &
        maxprocess,nsdecay,chdecay,maxdecay,numinv,numprocess,numdecay
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      if(switch.eq.0)return
      step=2
      do em=0,nexternal(generator)
      do sp=0,nexternal(generator)
      do ga=0,nexternal(generator)
        if(lmap(em,sp,ga))then
          do i1=1,nexternal(generator)
          do i2=0,3
            kt(i1,i2)=k(i1,i2)
          enddo
          enddo
!          call cofferaa_checkmom(k,nexternal(generator),switch)
! generating phasespace configurations
          call cofferaa_mapin(kt,gmap(em,sp,ga),em,sp,ga, &
            nexternal(generator),switch)
!          if(ga.gt.0)call cofferaa_checkmom(kt,nexternal(generator)-1,switch)
! calculating invariants
          do i1=1,allbinary(ga,generator)
            p(0,i1,em,sp,ga)=0d0
          enddo
          do i1=0,3
            p(i1,1,em,sp,ga)=-kt(1,i1)
            p(i1,2,em,sp,ga)=-kt(2,i1)
            do i2=3,nexternal(generator)
              p(i1,2**(i2-1),em,sp,ga)=kt(i2,i1)
            enddo
          enddo
! calculating invariants
          do i1=1,allbinary(ga,generator)-1
            if(p(0,allbinary(ga,generator)-i1,em,sp,ga).eq.0d0)then
              do i2=1,nexternal(generator)
                i3=i1/2**(i2-1)
                if(2*(i3/2).ne.i3.and.i1.ne.2**(i2-1))then
                  do i4=0,3
                    p(i4,i1,em,sp,ga)=p(i4,i1-2**(i2-1),em,sp,ga) &
                      +p(i4,2**(i2-1),em,sp,ga)
                  enddo
                  goto 100
                endif
              enddo
 100          s(i1,em,sp,ga)=p(0,i1,em,sp,ga)**2-p(1,i1,em,sp,ga)**2 &
                -p(2,i1,em,sp,ga)**2-p(3,i1,em,sp,ga)**2
            else
              do i2=0,3
                p(i2,i1,em,sp,ga)= &
                  -p(i2,allbinary(ga,generator)-i1,em,sp,ga)
              enddo
              s(i1,em,sp,ga)=s(allbinary(ga,generator)-i1,em,sp,ga)
            endif
          enddo
! external masses
          do i1=1,nexternal(generator)
          if(i1.ne.ga)then
            s(2**(i1-1),em,sp,ga)=massext2(i1,generator)
            s(allbinary(ga,generator)-2**(i1-1),em,sp,ga)= &
              massext2(i1,generator)
          endif
          enddo
        endif
      enddo
      enddo
      enddo
! inv
      do i1=1,maxinv(generator)
        ginv(i1)=1d0
        channel=chinv(i1,generator)
        em=emmap(channel,generator)
        sp=spmap(channel,generator)
        ga=gamap(channel,generator)
        ns=nsinv(i1,generator)
        inv1=ininv(ns,channel,generator)
        inv2=allbinary(ga,generator)-3-ininv(ns,channel,generator)
        mmin=mcutinv(inv1,generator)
        mmax=mcutinv(inv2,generator)
        do i2=1,ns-1
          virtinv=ininv(i2,channel,generator)
          if(lmin(i2,ns,channel,generator).and. &
            dsqrt(s(virtinv,em,sp,ga)).gt. &
            mcutinv(virtinv,generator))then
            mmin=mmin+dsqrt(s(virtinv,em,sp,ga)) &
              -mcutinv(inv1,generator) &
              +mcutinv(inv1-virtinv,generator)
            inv1=inv1-virtinv
          endif
          if(lmax(i2,ns,channel,generator).and. &
            dsqrt(s(virtinv,em,sp,ga)).gt. &
            mcutinv(virtinv,generator))then
            mmax=mmax+dsqrt(s(virtinv,em,sp,ga)) &
              -mcutinv(inv2,generator) &
              +mcutinv(inv2-virtinv,generator)
            inv2=inv2-virtinv
          endif
        enddo
        smin=mmin**2
        smax=(dsqrt(s(3,em,sp,ga))-mmax)**2
        call cofferaa_inv( &
          random(0), &                                      ! random number
          s(ininv(ns,channel,generator),em,sp,ga), &        ! invariant mass
          ginv(i1), &                                       ! local density
          mass(idhepinv(ns,channel,generator)), &           ! mass
          width(idhepinv(ns,channel,generator)), &          ! width
          powerinv(ns,channel,generator), &                 ! mapping
          smin, &                                           ! lower bound
          smax, &                                           ! upper bound
          step,switch)
      enddo
! process
      do i1=1,maxprocess(generator)
        gprocess(i1)=1d0
        channel=chprocess(i1,generator)
        em=emmap(channel,generator)
        sp=spmap(channel,generator)
        ga=gamap(channel,generator)
        nt=ntprocess(i1,generator)
        tcut1=0d0
        tcut2=0d0
        if(nt.eq.1)then
          tcut1=tcutprocess(allbinary(ga,generator) &
            -virtprocess(nt,channel,generator),generator)
          tcut2=tcutprocess(virtprocess(nt,channel,generator), &
            generator)
        endif
        call cofferaa_process( &
          random, &                                         ! random numbers
          p(0,in1process(nt,channel,generator),em,sp,ga), & ! incomming p. 1
          p(0,in2process(nt,channel,generator),em,sp,ga), & ! incomming p. 2
          p(0,out1process(nt,channel,generator),em,sp,ga), &! outgoing p. 1
          p(0,out2process(nt,channel,generator),em,sp,ga), &! outgoing p. 2
          p(0,virtprocess(nt,channel,generator),em,sp,ga), &! virtual p.
          gprocess(i1), &                                   ! local density
          mass(idhepprocess(nt,channel,generator)), &       ! mass
          width(idhepprocess(nt,channel,generator)), &      ! width
          powerprocess(nt,channel,generator), &             ! mapping
          s(inprocess(nt,channel,generator),em,sp,ga), &    ! incomming p.s
          s(out1process(nt,channel,generator),em,sp,ga), &  ! outgoing p. 1
          s(out2process(nt,channel,generator),em,sp,ga), &  ! outgoing p. 2
          s(virtprocess(nt,channel,generator),em,sp,ga), &  ! virtual p.
          s(in1process(nt,channel,generator),em,sp,ga), &   ! incomming p. 1
          s(in2process(nt,channel,generator),em,sp,ga), &   ! incomming p. 2
          tcut1,tcut2, &                                    ! angular cuts
          step,switch)
      enddo
! decay
      do i1=1,maxdecay(generator)
        gdecay(i1)=1d0
        channel=chdecay(i1,generator)
        em=emmap(channel,generator)
        sp=spmap(channel,generator)
        ga=gamap(channel,generator)
        ns=nsdecay(i1,generator)
        call cofferaa_decay( &
          random, &                                         ! random numbers
          p(0,indecay(ns,channel,generator),em,sp,ga), &    ! incomming p.
          p(0,out1decay(ns,channel,generator),em,sp,ga), &  ! outgoing p. 1
          p(0,out2decay(ns,channel,generator),em,sp,ga), &  ! outgoing p. 2
          gdecay(i1), &                                     ! local density
          s(indecay(ns,channel,generator),em,sp,ga), &      ! incomming p.
          s(out1decay(ns,channel,generator),em,sp,ga), &    ! outgoing p. 1
          s(out2decay(ns,channel,generator),em,sp,ga), &    ! outgoing p. 2
          step,switch)
       enddo
! calculating densities
      do channel=1,nchannel(generator)
        em=emmap(channel,generator)
        sp=spmap(channel,generator)
        ga=gamap(channel,generator)
        g(channel)=gmap(em,sp,ga)
! inv
        do ns=1,ninv(channel,generator)
          g(channel)=g(channel) &
            *ginv(numinv(ns,channel,generator))
        enddo
! process
        do nt=1,nprocess(channel,generator)
          g(channel)=g(channel) &
            *gprocess(numprocess(nt,channel,generator))
        enddo
! decay
        do ns=1,ndecay(channel,generator)
          g(channel)=g(channel) &
            *gdecay(numdecay(ns,channel,generator))
        enddo
      enddo
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialization generator                                    c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_initgenerator(energy,smin,hepnum,generator, &
        next,smself,sincludecuts,ssub,dipole_count,dipole_emitter, &
        dipole_unresolved,dipole_spectator)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv,maxo
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200,maxo=20)
      real*8 m2,m3,e2,e3,scutinv,energy,smin
      integer idhep(maxe,maxe),binary(maxe,maxe),idhep2,idhep3
      integer i1,i2,i3,i4,i5,i6,ns,nt,maxns,maxnt,binary1,binary2
      integer binary3,in1(maxe),in2(maxe),out1(maxe),out2(maxe), &
       cofferaa_pid
      integer virt(maxe),channel,generator,sincludecuts,prop2,prop3
      integer smself,smap,noutgen,hepnum(maxe),naux,nmap,next,ssub
      integer em,sp,ga,n
      logical cofferaa_vertex,cofferaa_included,exist, &
       cofferaa_comparechannel,cofferaa_compareinv
      logical cofferaa_comparedecay,cofferaa_compareprocess,schannel
      character*9 particle(maxv)
      character*80 vertices
      integer dipole_count,dipole_loop
      integer dipole_emitter(dipole_count)
      integer dipole_unresolved(dipole_count)
      integer dipole_spectator(dipole_count)
! mcparticle
      integer family(-maxv:maxv,6),light(-maxv:maxv)
      character*3 pname(-maxv:maxv),gname(-maxv:maxv)
! mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
! mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
! mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
! mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
! mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
! mcdensity
      integer nsinv(maxch,maxg),chinv(maxch,maxg),maxinv(maxg)
      integer ntprocess(maxch,maxg),chprocess(maxch,maxg)
      integer maxprocess(maxg),nsdecay(maxch,maxg),chdecay(maxch,maxg)
      integer maxdecay(maxg),numinv(maxe,maxch,maxg)
      integer numprocess(maxe,maxch,maxg),numdecay(maxe,maxch,maxg)
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
! mcoutput
      integer nout,numout,maxout
! mctechparam
      real*8 techparam(3)
! mcadaptopt
      real*8 alphaopt(maxch,maxo,maxg),betaopt(0:maxch,maxg)
      real*8 wi(maxch,maxg),alphamin
      integer nopt(0:maxo,maxg),opt(maxg)
! mccuts
      real*8 ecutp,ecutl,ecutq,scutqq,ccutpb,ccutpl,ccutpq,ccutll
      real*8 ccutqq,ccutlq,ccutlb,ccutqb,ecut(maxe,maxg)
      real*8 scut(maxe,maxe,maxg),ccut(maxe,maxe,maxg),xmin(maxg)
      common/mcparticle/family,light,pname,gname
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv, &
        lmin,lmax
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process, &
        in2process,out1process,out2process,inprocess,virtprocess, &
        idhepprocess,nprocess
      common/mcdensity/nsinv,chinv,maxinv,ntprocess,chprocess, &
        maxprocess,nsdecay,chdecay,maxdecay,numinv,numprocess,numdecay
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      common/mcoutput/nout,numout,maxout
      common/mctechparam/techparam
      common/mcadaptopt/alphaopt,betaopt,wi,alphamin,nopt,opt
      common/mccuts/ecutp,ecutl,ecutq,scutqq,ccutpb,ccutpl,ccutpq, &
        ccutll,ccutqq,ccutlq,ccutlb,ccutqb,ecut,scut,ccut,xmin
! options
      noutgen=1                     ! output for generator
      smap=1                        ! mapping
! technical parameter in h function and subroutine process
      techparam(1)=1d-4             ! small negative mass in h-function
      techparam(2)=1d-10*energy**2  ! allowed uncertainty tmax in process
      techparam(3)=1d-12*energy     ! allowed uncertainty q1,q2 in process
      powermap=0.8d0
! number of external particles
      nexternal(generator)=next
! initialize binary counting of particle number
      do i1=1,maxe
      do i2=1,maxe
        binary(i1,i2)=0
      enddo
      enddo
      do i1=1,nexternal(generator)
        binary(i1,1)=2**(i1-1)
      enddo
      allbinary(0,generator)=0
      do i1=1,nexternal(generator)
        allbinary(0,generator)=allbinary(0,generator)+2**(i1-1)
      enddo
      do i1=1,nexternal(generator)
        allbinary(i1,generator)=allbinary(0,generator)-2**(i1-1)
      enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialize particles                                        c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! generic particle names, see also subroutine vertexg
      do i1=1,maxv
        particle(i1)=' '
      enddo
! name, generic name, anti-field, family
      particle(1)='dq dq y 1'
      particle(2)='uq uq y 1'
      particle(3)='sq dq y 2'
      particle(4)='cq uq y 2'
      particle(5)='bq dq y 3'
      particle(6)='tq uq y 3'
      particle(11)='el el y 4'
      particle(12)='ne ne y 4'
      particle(13)='mu el y 5'
      particle(14)='nm ne y 5'
      particle(15)='ta el y 6'
      particle(16)='nt ne y 6'
      particle(22)='ph ph n  '
      particle(23)='Z0 Z0 n  '
      particle(24)='W- W- y  '
      particle(25)='h0 h0 n  '
      particle(26)='gl gl n  '
! supersymmetric particles
      particle(31)=''
      particle(32)=''
      particle(33)=''
      particle(34)=''
      particle(35)=''
      particle(36)=''
      particle(41)=''
      particle(42)=''
      particle(43)=''
      particle(44)=''
      particle(45)=''
      particle(46)=''
      particle(51)=''
      particle(52)=''
      particle(53)=''
      particle(54)=''
      particle(55)=''
      particle(56)=''
      particle(61)=''
      particle(63)=''
      particle(65)=''
      particle(71)=''
      particle(72)=''
      particle(73)=''
      particle(74)=''
      particle(70)=''
      particle(76)=''
      particle(81)=''
      particle(82)=''
      particle(83)=''
      particle(84)=''
! auxiliary particles
      naux=100
      particle(101)=' 1  1 y  '
      particle(102)=' 2  2 y  '
      particle(103)=' 3  3 y  '
      particle(104)=' 4  4 y  '
      particle(105)=' 5  5 y  '
      particle(106)=' 6  6 y  '
      particle(107)=' 7  7 y  '
      particle(108)=' 8  8 y  '
      particle(109)=''
      particle(110)=''
      particle(111)='11 11 y  '
      particle(112)='12 12 y  '
      particle(113)='13 13 y  '
      particle(114)='14 14 y  '
      particle(115)='15 15 y  '
      particle(116)=''  ! non-standard couplings
! additional particles for unstable particles
      nmap=150
      if(smap.eq.2)then
        i1=nmap
        do i2=1,naux-1
          if(width(i2).ne.0d0.and.particle(i2)(1:2).ne.'  ')then
            particle(i1)=particle(i2)
            i1=i1+1
            if(i1.gt.maxv)then
              write(*,'(a)')' initgenerator: maxv too small!'
              stop
            endif
          endif
        enddo
      endif
! setting particles
      gname(0)='   '
      do i1=1,maxv
        do i2=1,6
          family(i1,i2)=0
        enddo
        if(particle(i1)(9:9).eq.'1')family(i1,1)=1
        if(particle(i1)(9:9).eq.'2')family(i1,2)=1
        if(particle(i1)(9:9).eq.'3')family(i1,3)=1
        if(particle(i1)(9:9).eq.'4')family(i1,4)=1
        if(particle(i1)(9:9).eq.'5')family(i1,5)=1
        if(particle(i1)(9:9).eq.'6')family(i1,6)=1
        do i2=1,6
          family(-i1,i2)=-family(i1,i2)
        enddo
        pname(i1)=particle(i1)(1:3)
        pname(-i1)=pname(i1)
        gname(i1)=particle(i1)(4:6)
        gname(-i1)=gname(i1)
        if(particle(i1)(7:7).eq.'y')then
          pname(-i1)(3:3)='~'
          gname(-i1)(3:3)='~'
        endif
      enddo
! output of couplings
      if(nout.ne.0.and.noutgen.eq.2)then
        write(nout,'(a)')' '
        write(nout,'(a)')' Couplings:'
        i4=1
        vertices=' '
        do i1=-naux+1,naux-1
        do i2=i1+1,naux
        do i3=i2+1,naux
        if(cofferaa_vertex(i1,i2,i3,schannel,smself).and. &
          (gname(i1)(3:3).ne.' '.or.i1.gt.0).and. &
          (gname(i2)(3:3).ne.' '.or.i2.gt.0).and. &
          (gname(i3)(3:3).ne.' '.or.i3.gt.0))then
            vertices=vertices(1:i4)// &
              ' ('//pname(i1)//','//pname(i2)//','//pname(i3)//')'
            i4=i4+14
            if(i4.gt.60)then
              write(nout,'(a)')vertices
              vertices=' '
              i4=1
            endif
        endif
        enddo
        enddo
        enddo
        write(nout,'(a)')vertices
        i6=1
        vertices=' '
        do i1=naux,maxv
        do i2=-naux+1,naux-1
        do i3=i2+1,naux-1
        schannel=.true.
        if(cofferaa_vertex(i1,i2,i3,schannel,smself).and. &
          (gname(i2)(3:3).ne.' '.or.i2.gt.0).and. &
          (gname(i3)(3:3).ne.' '.or.i3.gt.0))then
          do i4=-naux+1,naux-1
          do i5=i4+1,naux-1
          if(cofferaa_vertex(-i1,i4,i5,schannel,smself).and. &
            (gname(i4)(3:3).ne.' '.or.i4.gt.0).and. &
            (gname(i5)(3:3).ne.' '.or.i5.gt.0))then
            vertices=vertices(1:i6)//' ('//pname(i2)//','// &
              pname(i3)//','//pname(i4)//','//pname(i5)//')'
            i6=i6+18
            if(i6.gt.60)then
              write(nout,'(a)')vertices
              vertices=' '
              i6=1
            endif
          endif
          enddo
          enddo
        endif
        enddo
        enddo
        enddo
        write(nout,'(a)')vertices
        write(nout,'(a)')' '
        write(nout,'(a)')' Channels:'
      endif
! mapping
      do i1=0,maxv
        if(mass(i1).eq.0d0)then
          power(i1)=0.8d0
        else
          power(i1)=2d0
        endif
        if(width(i1).gt.0d0)power(i1)=0d0
        if(i1.ge.naux)then
          mass(i1)=0d0
          width(i1)=0d0
          power(i1)=0d0
        endif
      enddo
! particle identity
      do i1=1,nexternal(generator)
        idhep(i1,1)=hepnum(i1)
        if(i1.ge.3)idhep(i1,1)=-idhep(i1,1)
        massext2(i1,generator)=mass(abs(hepnum(i1)))**2
      enddo
! uniform sampling for checks (after massext2 definition!)
      if(smap.eq.0)then
        do i1=1,naux-1
          if(mass(i1).ne.0d0)mass(i1)=1d-5*i1
          width(i1)=0d0
          power(i1)=1d-5*i1
        enddo
      endif
! small fermion-higgs couplings
      do i1=0,maxv
        light(i1)=0
        if(mass(i1).eq.0d0)light(i1)=1
        light(-i1)=light(i1)
      enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialize cuts                                             c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i1=0,allbinary(0,generator)
        mcutinv(i1,generator)=0d0
        if(i1.ge.1)tcutprocess(i1,generator)=0d0
      enddo
! angular cut
      if(sincludecuts.eq.1)then
        do i1=1,2
        do i2=3,nexternal(generator)
          if(ccut(i1,i2,generator).gt.0d0.and. &
            ccut(i1,i2,generator).lt.1d0)then
            tcutprocess(2**(i1-1)+2**(i2-1),generator)= &
              1d0/ccut(i1,i2,generator)**2-1d0
          endif
        enddo
        enddo
      endif
! invariant-mass cut
      do i1=1,allbinary(0,generator)
        scutinv=0d0
        do i2=1,nexternal(generator)
          if(cofferaa_included(i1,2**(i2-1),nexternal(generator)))then
             scutinv=scutinv+mass(abs(idhep(i2,1)))
           endif
        enddo
        scutinv=scutinv**2
        if(sincludecuts.eq.1)then
          do i2=1,nexternal(generator)
          do i3=i2+1,nexternal(generator)
            if(cofferaa_included(i1,2**(i2-1),nexternal(generator)).and. &
              cofferaa_included(i1,2**(i3-1),nexternal(generator)))then
              m2=mass(abs(idhep(i2,1)))
              m3=mass(abs(idhep(i3,1)))
              e2=max(ecut(i2,generator),m2)
              e3=max(ecut(i3,generator),m3)
              scutinv=scutinv-(m2+m3)**2 &
                +max(scut(i2,i3,generator),m2**2+m3**2+2d0*e2*e3 &
                    -2d0*dsqrt((e2**2-m2**2)*(e3**2-m3**2)) &
                      *ccut(i2,i3,generator))
            endif
          enddo
          enddo
        endif
        mcutinv(i1,generator)=dsqrt(scutinv)
        do i2=2,i1-1
          if(cofferaa_included(i1,i2,nexternal(generator)))then
            mcutinv(i1,generator)=max(mcutinv(i1,generator), &
              mcutinv(i2,generator)+mcutinv(i1-i2,generator))
          endif
        enddo
      enddo
      smin=mcutinv(allbinary(0,generator)-3,generator)**2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialize channels                                         c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      channel=nchannel(generator)+1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     loop over emitter and spectator                             c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do dipole_loop=0,dipole_count
      if(dipole_loop.eq.0)then
        em=0
        ga=0
        sp=0
      else
        em=dipole_emitter(dipole_loop)
        ga=dipole_unresolved(dipole_loop)
        sp=dipole_spectator(dipole_loop)
      endif
!      do em=0,nexternal(generator)
!      do sp=0,nexternal(generator)
!      do ga=0,nexternal(generator)
      n=0
      lmap(em,sp,ga)=.false.
      if(em.eq.0.and.sp.eq.0.and.ga.eq.0)then
        n=nexternal(generator)
      elseif(ssub.ne.0.and.em.ne.sp.and.em.ne.ga.and.sp.ne.ga.and. &
        em.gt.0.and.sp.gt.0.and.ga.gt.0)then
        if(pname(hepnum(ga)).eq.'gl '.and.ga.ge.3.and. &
          cofferaa_vertex(hepnum(em),-hepnum(em),hepnum(ga), &
          schannel,smself).and.cofferaa_vertex(hepnum(sp),-hepnum(sp), &
          hepnum(ga),schannel,smself))then
          n=nexternal(generator)-1
          idhep(ga,1)=0
        endif
      endif
      if(n.ne.0)then
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     s-channel propagators                                       c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! number of s-channel propagators
      do maxns=0,n-3
        ns=0
        if(maxns.eq.0)goto 500
! ns'th s-channel propagator
 100    ns=ns+1
        out1(ns)=2
! decay product 1 for ns'th decay
 200    out1(ns)=out1(ns)+1
        out2(ns)=out1(ns)
! decay product 2 for ns'th decay
 300    out2(ns)=out2(ns)+1
        virt(ns)=-maxv-1
! virtual particle for ns'th decay
 400    virt(ns)=virt(ns)+1
        if(gname(virt(ns)).eq.'   '.and.virt(ns).lt.maxv)goto 400
        if(gname(virt(ns)).eq.'   '.and.virt(ns).eq.maxv)goto 1500
        if(virt(ns).lt.0.and. &
          gname(virt(ns)).eq.gname(-virt(ns)))goto 400
! checking whether 3-particle vertex exists
        schannel=.true.
        if(.not.cofferaa_vertex(idhep(out1(ns),ns),idhep(out2(ns),ns), &
          virt(ns),schannel,smself))goto 1500
! checking last 3-particle vertex
        if(ns.eq.n-3.and..not.cofferaa_vertex(idhep(in1(ns),ns), &
          idhep(in2(ns),ns),-virt(ns),schannel,smself))goto 1500
! initializing next step
        do i2=1,n
          binary(i2,ns+1)=binary(i2,ns)
          idhep(i2,ns+1)=idhep(i2,ns)
        enddo
! combining particle out1 and out2 into new external particle out2
        binary(out1(ns),ns+1)=0
        idhep(out1(ns),ns+1)=0
        binary(out2(ns),ns+1)= &
          binary(out1(ns),ns)+binary(out2(ns),ns)
        idhep(out2(ns),ns+1)=-virt(ns)
! find n s-channel propagator for ns < maxns
        if(ns.lt.maxns)goto 100
 500    continue
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     t-channel propagators                                       c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! number of t-channel propagators
        maxnt=n-3
! loop of t-channel propagators
        if(maxns.eq.maxnt)goto 1000
        nt=maxns
! nt'th t-channel propagator
 600    nt=nt+1
        in1(nt)=0
! incoming particle 1 for nt's process
 700    in1(nt)=in1(nt)+1
! incoming particle 2 for nt's process
        in2(nt)=3-in1(nt)
        out1(nt)=2
! outgoing particle 1 for nt's process
 800    out1(nt)=out1(nt)+1
! virtual particle for nt's process
        virt(nt)=-nmap                  ! (-nmap <-> -maxv-1)
 900    virt(nt)=virt(nt)+1
        if(gname(virt(nt)).eq.'   '.and.virt(nt).lt.maxv)goto 900
        if(gname(virt(nt)).eq.'   '.and.virt(nt).eq.maxv)goto 1300
        if(virt(nt).lt.0.and. &
          gname(virt(nt)).eq.gname(-virt(nt)))goto 900
! avoid doube counting of diagrams
        if(nt.gt.maxns+1.and.in1(nt-1).eq.in1(nt))goto 1300
! checking whether 3-particle vertex exists
        schannel=.false.
        if(.not.cofferaa_vertex(idhep(in1(nt),nt),idhep(out1(nt),nt), &
          virt(nt),schannel,smself))goto 1300
! initializing n step
        do i2=1,n
          binary(i2,nt+1)=binary(i2,nt)
          idhep(i2,nt+1)=idhep(i2,nt)
        enddo
! combining particle in1 and out1 into new external particle in1
        binary(out1(nt),nt+1)=0
        idhep(out1(nt),nt+1)=0
        binary(in1(nt),nt+1)=binary(in1(nt),nt)+binary(out1(nt),nt)
        idhep(in1(nt),nt+1)=-virt(nt)
! find n t-channel propagator for nt < maxnt
        if(nt.lt.maxnt)goto 600
 1000   continue
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     checking last vertex                                        c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        i1=n-2
        do i2=3,n
          if(idhep(i2,i1).ne.0)i3=i2
        enddo
        schannel=.true.
        if(gname(idhep(i3,i1)).eq.'   '.or..not. &
          cofferaa_vertex(idhep(1,i1),idhep(2,i1),idhep(i3,i1), &
          schannel,smself)) goto 1300
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initializing channels                                       c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ninv(channel,generator)=n-4
        nprocess(channel,generator)=maxnt-maxns
        ndecay(channel,generator)=maxns
! initializing inv
        idhepfirst(channel,generator)=0
        if(nprocess(channel,generator).eq.0) &
          idhepfirst(channel,generator)=abs(idhep(i3,i1))
        do i1=1,maxns
          ininv(i1,channel,generator)= &
            binary(out1(i1),i1)+binary(out2(i1),i1)
          idhepinv(i1,channel,generator)=abs(virt(i1))
          powerinv(i1,channel,generator)=power(abs(virt(i1)))
        enddo
        do i1=maxns+1,maxnt-1
          i2=maxns+maxnt-i1
          ininv(i2,channel,generator)=allbinary(ga,generator)  &
            -binary(in1(i1),i1)-binary(in2(i1),i1)-binary(out1(i1),i1)
          idhepinv(i2,channel,generator)=0
          powerinv(i2,channel,generator)=0d0
        enddo
! minimal invariant-mass cuts for subroutine inv
        do i1=1,ninv(channel,generator)
        do i2=1,i1-1
          binary1=ininv(i1,channel,generator)
          binary2=ininv(i2,channel,generator)
          lmin(i2,i1,channel,generator)=.false.
          if(cofferaa_included(binary1,binary2,n))then
            lmin(i2,i1,channel,generator)=.true.
            do i3=1,i1-1
              binary3=ininv(i3,channel,generator)
              if(cofferaa_included(binary3,binary2,n) &
                .and.i2.ne.i3)then
                lmin(i2,i1,channel,generator)=.false.
              endif
            enddo
          endif
        enddo
        enddo
! maximal invariant-mass cuts for subroutine inv
        do i1=1,ninv(channel,generator)
        do i2=1,i1-1
          binary1=allbinary(ga,generator)-3-ininv(i1,channel,generator)
          binary2=ininv(i2,channel,generator)
          lmax(i2,i1,channel,generator)=.false.
          if(cofferaa_included(binary1,binary2,n))then
            lmax(i2,i1,channel,generator)=.true.
            do i3=1,i1-1
              binary3=ininv(i3,channel,generator)
              if(cofferaa_included(binary3,binary2,n) &
                .and.i2.ne.i3)then
                lmax(i2,i1,channel,generator)=.false.
              endif
            enddo
          endif
        enddo
        enddo
! initializing process
        do i1=maxns+1,maxnt
          i2=i1-maxns
          in1process(i2,channel,generator)=allbinary(ga,generator) &
            -binary(in1(i1),i1)
          in2process(i2,channel,generator)=allbinary(ga,generator) &
            -binary(in2(i1),i1)
          out1process(i2,channel,generator)=binary(out1(i1),i1)
          out2process(i2,channel,generator)=allbinary(ga,generator)  &
            -binary(in1(i1),i1)-binary(in2(i1),i1)-binary(out1(i1),i1)
          inprocess(i2,channel,generator)=allbinary(ga,generator) &
            -binary(in1(i1),i1)-binary(in2(i1),i1)
          virtprocess(i2,channel,generator)=allbinary(ga,generator) &
            -binary(in1(i1),i1)-binary(out1(i1),i1)
          idhepprocess(i2,channel,generator)=abs(virt(i1))
          powerprocess(i2,channel,generator)=power(abs(virt(i1)))
        enddo
! initializing decay
        do i1=1,maxns
          i2=maxns-i1+1
          indecay(i2,channel,generator)= &
            binary(out1(i1),i1)+binary(out2(i1),i1)
          out1decay(i2,channel,generator)=binary(out1(i1),i1)
          out2decay(i2,channel,generator)=binary(out2(i1),i1)
        enddo
! initializing map
         lmap(em,sp,ga)=.true.
         emmap(channel,generator)=em
         spmap(channel,generator)=sp
         gamap(channel,generator)=ga
! checking whether generator already exists
        do i1=1,channel-1
        if(cofferaa_comparechannel(i1,channel,smap,generator))then
          if(tcutprocess(allbinary(ga,generator) &
            -virtprocess(1,channel,generator),generator).ne.0d0.or. &
            tcutprocess(virtprocess(1,channel,generator),generator) &
            .ne.0d0)then
            call cofferaa_copychannel(channel,i1,generator)
          endif
          goto 1300
        endif
        enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     output                                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(nout.ne.0.and.noutgen.eq.2)then
          write(nout,'(a12,i6)')' channel   =',channel
          if(em.ne.0.or.sp.ne.0.or.ga.ne.0)then
            write(nout,'(a12,i6)')' emitter   =',em
            write(nout,'(a12,i6)')' spectator =',sp
            write(nout,'(a12,i6)')' photon    =',ga
          endif
          do i1=1,maxns
            write(nout,'(a30,i6," +",i6," ->",i6)') &
              '   s channel: ('//pname(-idhep(out1(i1),i1))//','// &
              pname(-idhep(out2(i1),i1))//','// &
              pname(-idhep(out2(i1),i1+1))//'), ', &
              cofferaa_pid(binary(out1(i1),i1)), &
              cofferaa_pid(binary(out2(i1),i1)), &
              cofferaa_pid(binary(out2(i1),i1+1))
          enddo
          do i1=maxns+1,maxnt
            write(nout,'(a30,i6," +",i6," ->",i6)') &
              '   t channel: ('//pname(idhep(in1(i1),i1))//','// &
              pname(-idhep(out1(i1),i1))//','// &
              pname(idhep(in1(i1),i1+1))//'), ', &
              cofferaa_pid(binary(in1(i1),i1)), &
              cofferaa_pid(binary(out1(i1),i1)), &
              cofferaa_pid(binary(in1(i1),i1+1))
          enddo
        endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     end of loops                                                c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        channel=channel+1
 1300   if(maxns.lt.maxnt)then
 1400     if(virt(nt).lt.nmap-1)goto 900    ! (nmap-1 <-> maxv)
          if(out1(nt).lt.n)goto 800
          if(in1(nt).lt.2)goto 700
          if(nt.gt.maxns+1)then
            nt=nt-1
            goto 1400
          endif
        endif
! end of t-channel propagators
 1500   if(maxns.gt.0)then
 1600     if(virt(ns).lt.maxv)goto 400
          if(out2(ns).lt.n)goto 300
          if(out1(ns).lt.n-1)goto 200
          if(ns.gt.1)then
            ns=ns-1
            goto 1600
          endif
        endif
! end of s-channel propagators
      enddo
! number of channels
      nchannel(generator)=channel-1
      if(nchannel(generator).gt.maxch)then
        write(*,'(a28,i6)')' initgenerator: reset maxch >', &
          nchannel(generator)
        stop
      endif
      if(nchannel(generator).eq.0)then
        write(*,'(a47,i2)') &
          ' initgenerator: no channel found for generator ',generator
        stop
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialize densities                                        c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! initializing inv
      maxinv(generator)=0
      do channel=1,nchannel(generator)
      do ns=1,ninv(channel,generator)
        numinv(ns,channel,generator)=0
        do i1=1,channel-1
        do i2=1,ninv(i1,generator)
          if(cofferaa_compareinv(ns,channel,i2,i1,generator)) &
            numinv(ns,channel,generator)=numinv(i2,i1,generator)
        enddo
        enddo
        if(numinv(ns,channel,generator).eq.0)then
          maxinv(generator)=maxinv(generator)+1
          numinv(ns,channel,generator)=maxinv(generator)
          nsinv(maxinv(generator),generator)=ns
          chinv(maxinv(generator),generator)=channel
        endif
      enddo
      enddo
! process
      maxprocess(generator)=0
      do channel=1,nchannel(generator)
      do nt=1,nprocess(channel,generator)
        numprocess(nt,channel,generator)=0
        do i1=1,channel-1
        do i2=1,nprocess(i1,generator)
          if(cofferaa_compareprocess(nt,channel,i2,i1,generator)) &
            numprocess(nt,channel,generator)=numprocess(i2,i1,generator)
        enddo
        enddo
        if(numprocess(nt,channel,generator).eq.0)then
          maxprocess(generator)=maxprocess(generator)+1
          numprocess(nt,channel,generator)=maxprocess(generator)
          ntprocess(maxprocess(generator),generator)=nt
          chprocess(maxprocess(generator),generator)=channel
        endif
      enddo
      enddo
! decay
      maxdecay(generator)=0
      do channel=1,nchannel(generator)
      do ns=1,ndecay(channel,generator)
        numdecay(ns,channel,generator)=0
        do i1=1,channel-1
        do i2=1,ndecay(i1,generator)
          if(cofferaa_comparedecay(ns,channel,i2,i1,generator)) &
            numdecay(ns,channel,generator)=numdecay(i2,i1,generator)
        enddo
        enddo
        if(numdecay(ns,channel,generator).eq.0)then
          maxdecay(generator)=maxdecay(generator)+1
          numdecay(ns,channel,generator)=maxdecay(generator)
          nsdecay(maxdecay(generator),generator)=ns
          chdecay(maxdecay(generator),generator)=channel
        endif
      enddo
      enddo
      endif
      enddo
!      enddo
!      enddo
      if(nout.ne.0.and.noutgen.ge.1)then
! output
        write(nout,'(a)')' '
        write(nout,'(a23,i2)')' Phase-space generator ',generator
        write(nout,'(" Number of channels            =",i6)') &
          nchannel(generator)
        write(nout,'(" Calculation of invariants     =",i6)') &
          maxinv(generator)
        write(nout,'(" Calculation of 2->2 processes =",i6)') &
          maxprocess(generator)
        write(nout,'(" Calculation of 1->2 decays    =",i6)') &
          maxdecay(generator)
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialize adaptive optimization                            c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i1=0,maxo
        nopt(i1,generator)=0
      enddo
      nopt(1,generator)=1000
      nopt(2,generator)=2000
      nopt(3,generator)=3000
      nopt(4,generator)=4000
      nopt(5,generator)=5000
      nopt(6,generator)=6000
      nopt(7,generator)=7000
      nopt(8,generator)=8000
      do i1=0,maxo
        nopt(i1,generator)=nopt(i1,generator)*nchannel(generator)
      enddo
      opt(generator)=1
      betaopt(0,generator)=0d0
      do i1=1,nchannel(generator)
        wi(i1,generator)=0d0
        alphaopt(i1,1,generator)=1d0/nchannel(generator)
        betaopt(i1,generator)=dble(i1)/nchannel(generator)
      enddo
      alphamin=0.01d0
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check whether generator already exists                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_comparechannel(ch1,ch2,smap,generator)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      integer i1,i2,ch1,ch2,generator,prop1,prop2,smap,idhep1,idhep2
      integer binary1,binary2,ga
      logical cofferaa_comparechannel,same,exist
! mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
! mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
! mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
! mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
! mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process, &
        in2process,out1process,out2process,inprocess,virtprocess, &
        idhepprocess,nprocess
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv, &
        lmin,lmax
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      cofferaa_comparechannel=.false.
! checking emitter, spectator, and photon
      if(emmap(ch1,generator).ne.emmap(ch2,generator))return
      if(spmap(ch1,generator).ne.spmap(ch2,generator))return
      if(gamap(ch1,generator).ne.gamap(ch2,generator))return
      ga=gamap(ch1,generator)
! checking number of s-channel propagators
      prop1=ndecay(ch1,generator)
      do i1=1,ndecay(ch1,generator)
        idhep1=idhepinv(i1,ch1,generator)
        if(mass(idhep1).eq.0d0.and.width(idhep1).eq.0d0.and. &
          power(idhep1).eq.0d0)prop1=prop1-1
      enddo
      prop2=ndecay(ch2,generator)
      do i1=1,ndecay(ch2,generator)
        idhep2=idhepinv(i1,ch2,generator)
        if(mass(idhep2).eq.0d0.and.width(idhep2).eq.0d0.and. &
          power(idhep2).eq.0d0)prop2=prop2-1
      enddo
! checking s-channel propagators
      if(prop1.eq.prop2)then
        do i1=1,ndecay(ch1,generator)
          binary1=ininv(i1,ch1,generator)
          idhep1=idhepinv(i1,ch1,generator)
          if(mass(idhep1).ne.0d0.or.width(idhep1).ne.0d0.or. &
            power(idhep1).ne.0d0)then
            exist=.false.
            do i2=1,ndecay(ch2,generator)
              binary2=ininv(i2,ch2,generator)
              idhep2=idhepinv(i2,ch2,generator)
              if(binary1.eq.binary2.and. &
                mass(idhep1).eq.mass(idhep2).and. &
                width(idhep1).eq.width(idhep2).and. &
                power(idhep1).eq.power(idhep2))exist=.true.
            enddo
            if(.not.exist)return
          endif
        enddo
      else
        return
      endif
! checking number of t-channel propagators
      prop1=nprocess(ch1,generator)
      do i1=1,nprocess(ch1,generator)
        idhep1=idhepprocess(i1,ch1,generator)
        if(mass(idhep1).eq.0d0.and.width(idhep1).eq.0d0.and. &
          power(idhep1).eq.0d0)prop1=prop1-1
      enddo
      prop2=nprocess(ch2,generator)
      do i1=1,nprocess(ch2,generator)
        idhep2=idhepprocess(i1,ch2,generator)
        if(mass(idhep2).eq.0d0.and.width(idhep2).eq.0d0.and. &
          power(idhep2).eq.0d0)prop2=prop2-1
      enddo
! checking t-channel propagators
      if(prop1.eq.prop2)then
        do i1=1,nprocess(ch1,generator)
          binary1=min(virtprocess(i1,ch1,generator), &
            allbinary(ga,generator)-virtprocess(i1,ch1,generator))
          idhep1=idhepprocess(i1,ch1,generator)
          if(mass(idhep1).ne.0d0.or.width(idhep1).ne.0d0.or. &
            power(idhep1).ne.0d0)then
            exist=.false.
            do i2=1,nprocess(ch2,generator)
              binary2=min(virtprocess(i2,ch2,generator), &
                allbinary(ga,generator)-virtprocess(i2,ch2,generator))
              idhep2=idhepprocess(i2,ch2,generator)
              if(binary1.eq.binary2.and. &
                mass(idhep1).eq.mass(idhep2).and. &
                width(idhep1).eq.width(idhep2).and. &
                power(idhep1).eq.power(idhep2))then
                  exist=.true.
              endif
            enddo
            if(.not.exist)return
          endif
        enddo
      else
        return
      endif
! checking first s-channel propagator
      if(smap.eq.0.and.nprocess(ch1,generator).eq.0.and. &
        nprocess(ch2,generator).eq.0)then
        idhep1=idhepfirst(ch1,generator)
        idhep2=idhepfirst(ch2,generator)
        if(mass(idhep1).ne.mass(idhep2).or. &
          width(idhep1).ne.width(idhep2).or. &
          power(idhep1).ne.power(idhep2))return
      endif
      cofferaa_comparechannel=.true.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     copy channel ch1 to channel ch2                             c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_copychannel(ch1,ch2,generator)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      integer i1,i2,ch1,ch2,generator
! mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
! mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
! mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
! mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mckinematics/massext2,allbinary,nexternal
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process, &
        in2process,out1process,out2process,inprocess,virtprocess, &
        idhepprocess,nprocess
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv, &
        lmin,lmax
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      do i1=1,maxe
! inv
        powerinv(i1,ch2,generator)=powerinv(i1,ch1,generator)
        ininv(i1,ch2,generator)=ininv(i1,ch1,generator)
        idhepinv(i1,ch2,generator)=idhepinv(i1,ch1,generator)
        idhepfirst(ch2,generator)=idhepfirst(ch1,generator)
        ninv(ch2,generator)=ninv(ch1,generator)
        do i2=1,maxe
          lmin(i1,i2,ch2,generator)=lmin(i1,i2,ch1,generator)
          lmax(i1,i2,ch2,generator)=lmax(i1,i2,ch1,generator)
        enddo
! decay
        indecay(i1,ch2,generator)=indecay(i1,ch1,generator)
        out1decay(i1,ch2,generator)=out1decay(i1,ch1,generator)
        out2decay(i1,ch2,generator)=out2decay(i1,ch1,generator)
        ndecay(ch2,generator)=ndecay(ch1,generator)
! process
        powerprocess(i1,ch2,generator)=powerprocess(i1,ch1,generator)
        in1process(i1,ch2,generator)=in1process(i1,ch1,generator)
        in2process(i1,ch2,generator)=in2process(i1,ch1,generator)
        out1process(i1,ch2,generator)=out1process(i1,ch1,generator)
        out2process(i1,ch2,generator)=out2process(i1,ch1,generator)
        inprocess(i1,ch2,generator)=inprocess(i1,ch1,generator)
        virtprocess(i1,ch2,generator)=virtprocess(i1,ch1,generator)
        idhepprocess(i1,ch2,generator)=idhepprocess(i1,ch1,generator)
        nprocess(ch2,generator)=nprocess(ch1,generator)
      enddo
! mapping
      emmap(ch2,generator)=emmap(ch1,generator)
      spmap(ch2,generator)=spmap(ch1,generator)
      gamap(ch2,generator)=gamap(ch1,generator)
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check whether two 1->2 decays are equal                     c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_comparedecay(ns1,ch1,ns2,ch2,generator)
      implicit none
! local variables
      integer maxe,maxch,maxg
      parameter(maxe=9,maxch=20000,maxg=4)
      integer ns1,ns2,ch1,ch2,generator
      logical cofferaa_comparedecay
! mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      cofferaa_comparedecay=.false.
      if(emmap(ch1,generator).ne.emmap(ch2,generator))return
      if(spmap(ch1,generator).ne.spmap(ch2,generator))return
      if(gamap(ch1,generator).ne.gamap(ch2,generator))return
      if(indecay(ns1,ch1,generator).ne.indecay(ns2,ch2,generator)) &
        return
      if(out1decay(ns1,ch1,generator).ne.out1decay(ns2,ch2,generator) &
        .and.out1decay(ns1,ch1,generator).ne. &
        out2decay(ns2,ch2,generator))return
      cofferaa_comparedecay=.true.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check whether two 2->2 processes are equal                  c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_compareprocess(ns1,ch1,ns2,ch2,generator)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      integer ns1,ns2,ch1,ch2,generator,idhep1,idhep2,ga
      logical cofferaa_compareprocess
! mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
! mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
! mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcprocess/powerprocess,tcutprocess,in1process, &
        in2process,out1process,out2process,inprocess,virtprocess, &
        idhepprocess,nprocess
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      cofferaa_compareprocess=.false.
      if(emmap(ch1,generator).ne.emmap(ch2,generator))return
      if(spmap(ch1,generator).ne.spmap(ch2,generator))return
      if(gamap(ch1,generator).ne.gamap(ch2,generator))return
      ga=gamap(ch1,generator)
      if(inprocess(ns1,ch1,generator).ne.inprocess(ns2,ch2,generator)) &
        return
      if(virtprocess(ns1,ch1,generator).ne. &
        virtprocess(ns2,ch2,generator).and. &
        virtprocess(ns1,ch1,generator).ne. &
        allbinary(ga,generator)-virtprocess(ns2,ch2,generator))return
      if(powerprocess(ns1,ch1,generator).ne. &
        powerprocess(ns2,ch2,generator))return
      idhep1=idhepprocess(ns1,ch1,generator)
      idhep2=idhepprocess(ns2,ch2,generator)
      if(mass(idhep1).ne.mass(idhep2).or.width(idhep1).ne.width(idhep2)) &
        return
      cofferaa_compareprocess=.true.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check whether calculation of invariants are equal           c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_compareinv(ns1,ch1,ns2,ch2,generator)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=4,maxv=200)
      integer i1,i2,ns1,ns2,ch1,ch2,generator,idhep1,idhep2
      logical cofferaa_compareinv,included
! mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
! mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv, &
        lmin,lmax
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      cofferaa_compareinv=.false.
      if(emmap(ch1,generator).ne.emmap(ch2,generator))return
      if(spmap(ch1,generator).ne.spmap(ch2,generator))return
      if(gamap(ch1,generator).ne.gamap(ch2,generator))return
      if(ininv(ns1,ch1,generator).ne.ininv(ns2,ch2,generator))return
      if(powerinv(ns1,ch1,generator).ne.powerinv(ns2,ch2,generator)) &
        return
      idhep1=idhepinv(ns1,ch1,generator)
      idhep2=idhepinv(ns2,ch2,generator)
      if(mass(idhep1).ne.mass(idhep2).or.width(idhep1).ne.width(idhep2)) &
        return
      do i1=1,ns1-1
        if(lmin(i1,ns1,ch1,generator))then
          included=.false.
          do i2=1,ns2-1
            if(lmin(i2,ns2,ch2,generator).and. &
              ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator)) &
              included=.true.
          enddo
          if(.not.included)return
        endif
        if(lmax(i1,ns1,ch1,generator))then
          included=.false.
          do i2=1,ns2-1
            if(lmax(i2,ns2,ch2,generator).and. &
              ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator)) &
              included=.true.
          enddo
          if(.not.included)return
        endif
      enddo
      do i2=1,ns2-1
        if(lmin(i2,ns2,ch2,generator))then
          included=.false.
          do i1=1,ns1-1
            if(lmin(i1,ns1,ch1,generator).and. &
              ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator)) &
              included=.true.
          enddo
          if(.not.included)return
        endif
        if(lmax(i2,ns2,ch2,generator))then
          included=.false.
          do i1=1,ns1-1
            if(lmax(i1,ns1,ch1,generator).and. &
              ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator)) &
              included=.true.
          enddo
          if(.not.included)return
        endif
      enddo
      cofferaa_compareinv=.true.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check if binary2 is included in binary1                     c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_included(binary1,binary2,nexternal)
      implicit none
! local variables
      integer binary1,binary2,i1,b1,b2,nexternal
      logical cofferaa_included
      cofferaa_included=.false.
      do i1=1,nexternal
        b1=binary1/2**(i1-1)
        b2=binary2/2**(i1-1)
        if(2*(b2/2).ne.b2.and.2*(b1/2).eq.b1)return
      enddo
      cofferaa_included=.true.
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                c
!     s-channel propagator                                       c
!                                                                c
!     written by Markus Roth                                     c
!                                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_inv(random,s,g,mass,width,power,smin,smax, &
       step,switch)
      implicit none
! local variable
      real*8 pi,random,s,g,mass,mass2,width,width2,power
      real*8 smax,smin,omax,omin,cofferaa_h,cofferaa_jacobian,denum
      integer step,switch
! mcoutput
      integer nout,numout,maxout
      common/mcoutput/nout,numout,maxout
      pi=4d0*datan(1d0)
      if(switch.eq.0)return
      mass2=mass*mass
      if(smax.le.smin)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' inv: smin >= smax'
          write(nout,'(6x,"smin=",d16.10," smax=",d16.10)')smin,smax
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(smin.lt.0d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' inv: smin < 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(width.gt.0d0)then
        width2=width*width
        omax=datan((smax-mass2)/mass/width)
        omin=datan((smin-mass2)/mass/width)
        if(step.eq.1)then
          s=mass2+mass*width*dtan(random*(omax-omin)+omin)
        elseif(step.eq.2)then
          denum=(omax-omin)*((s-mass2)**2+mass2*width2)
          if(denum.gt.0d0)then
            g=g*mass*width/denum
          else
            if(numout.lt.maxout)then
              write(nout,'(a)')' inv: 1/density <= 0 (width>0)'
              numout=numout+1
            endif
            switch=0
            return
          endif
        endif
      else
        if(step.eq.1)then
          s=cofferaa_h(random,power,mass2,smin,smax,switch)
        elseif(step.eq.2)then
          denum=cofferaa_jacobian(power,mass2,s,smin,smax,switch)
          if(denum.gt.0d0)then
            g=g/denum
          else
            if(nout.lt.maxout)then
              write(nout,'(a)')' inv: 1/density <= 0 (width=0)'
            endif
            numout=numout+1
            switch=0
            return
          endif
        endif
      endif
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                c
!     2->2 particle process with t-channel propagator            c
!                                                                c
!     written by Markus Roth                                     c
!                                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_process(random,q1,q2,p1,p2,qt,g,mass,width, &
        nu,s,s1,s2,t,t1,t2,tcut1,tcut2,step,switch)
      implicit none
! local variable
      real*8 pi,random(2),q1(0:3),q2(0:3),p1(0:3),p2(0:3),qt(0:3)
      real*8 g,mass,width,nu,s1,s2,t1,t2,tmin,tmax,tcut1,tcut2
      real*8 lambdas,lambdat,mass2,width2,t,phi,cost,s,roots
      real*8 gamma2,beta,gs,d,lambda,be,omega,omin,omax
      real*8 k1(0:3),q(0:3),cofferaa_h,cofferaa_jacobian,denum,cmax, &
       cmax1,cmax2
      integer i1,step,switch
! mcoutput
      integer nout,numout,maxout
! mctechparam
      real*8 techparam(3)
      common/mcoutput/nout,numout,maxout
      common/mctechparam/techparam
      if(switch.eq.0)return
      pi=4d0*datan(1d0)
      if(tcut1.ne.0d0.or.tcut2.ne.0d0)then
        beta=(q1(3)+q2(3))/(q1(0)+q2(0))
        gamma2=1d0/(1d0-beta*beta)
        lambda=(s-s1-s2)**2-4d0*s1*s2
        if(dabs(q1(1)+q2(1)).gt.techparam(3).and. &
          dabs(q1(2)+q2(2)).gt.techparam(3))then
          if(numout.lt.maxout)then
            write(nout,'(a)')' process: dabs(q1(1:2)+q2(1:2)) > acc '
            numout=numout+1
          endif
          switch=0
          return
        elseif(lambda.lt.0d0)then
          if(numout.lt.maxout)then
            write(nout,'(a)')' process: lambda < 0'
            numout=numout+1
          endif
          switch=0
          return
        else
          if(tcut1.ne.0d0)then    ! tcut1=tan**2(theta_1)
            be=beta*dsign(1d0,q1(3))
            gs=be*(s+s1-s2)/dsqrt(lambda)
            d=1d0+gamma2*(1d0-gs*gs)*tcut1
            if(d.lt.0d0)then
              if(numout.lt.maxout)then
                write(nout,'(a)')' process: d < 0'
                numout=numout+1
              endif
              switch=0
              return
            else
              cmax1=(-gs*gamma2*tcut1+dsqrt(d))/(1d0+gamma2*tcut1)
            endif
          else
            cmax1=1d0
          endif
          if(tcut2.ne.0d0)then    ! tcut2=tan**2(theta_2)
            be=beta*dsign(1d0,q2(3))
            gs=be*(s+s2-s1)/dsqrt(lambda)
            d=1d0+gamma2*(1d0-gs*gs)*tcut2
            if(d.lt.0d0)then
              if(numout.lt.maxout)then
                write(nout,'(a)')' process: d < 0'
                numout=numout+1
              endif
              switch=0
              return
            else
              cmax2=(-gs*gamma2*tcut2+dsqrt(d))/(1d0+gamma2*tcut2)
            endif
          else
            cmax2=1d0
          endif
        endif
        cmax=min(cmax1,cmax2)
      else
        cmax=1d0
      endif
      mass2=mass*mass
      lambdas=(s-s1-s2)*(s-s1-s2)-4d0*s1*s2
      if(lambdas.gt.0d0)then
        lambdas=dsqrt(lambdas)
      else
        if(numout.lt.maxout)then
          write(nout,'(a)')' process: lambdas <= 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      lambdat=(s-t1-t2)*(s-t1-t2)-4d0*t1*t2
      if(lambdat.gt.0d0)then
        lambdat=dsqrt(lambdat)
      else
        if(numout.lt.maxout)then
          write(nout,'(a)')' process: lambdat <= 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      tmin=s1+t1-((s+s1-s2)*(s+t1-t2)+lambdas*lambdat)/2d0/s
      tmax=s1+t1-((s+s1-s2)*(s+t1-t2)-lambdas*lambdat*cmax)/2d0/s
      if(tmax.le.tmin)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' process: tmin >= tmax'
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(tmax.gt.0d0)then
        if(tmax.lt.techparam(2))then
          tmax=0d0
        else
          if(numout.lt.maxout)then
            write(nout,'(a)')' process: tmax > 0'
            numout=numout+1
          endif
          switch=0
          return
        endif
      endif
      if(width.gt.0d0)then
        width2=width*width
        omin=datan((mass2-tmax)/mass/width)
        omax=datan((mass2-tmin)/mass/width)
        if(step.eq.1)then
          phi=2d0*pi*random(1)
          omega=omin+(omax-omin)*random(2)
          cost=((s+s1-s2)*(s+t1-t2)-2d0*s*(t1+s1)+2d0*mass2*s &
                -2d0*s*width*mass*dtan(omega))/lambdas/lambdat
        elseif(step.eq.2)then
          denum=pi*(omax-omin)*((mass2-t)**2+mass2*width2)
          if(denum.ne.0d0)then
            g=g*2d0*lambdat*mass*width/denum
          else
            if(numout.lt.maxout)then
              write(nout,'(a)')' process: 1/density = 0 (width>0)'
              numout=numout+1
            endif
            switch=0
            return
          endif
        endif
      else
        if(step.eq.1)then
          phi=2d0*pi*random(1)
          denum=lambdas*lambdat
          cost=((s+s1-s2)*(s+t1-t2)-2d0*s*(t1+s1) &
                -2d0*s*cofferaa_h(random(2),nu,-mass2,-tmax,-tmin, &
                switch))/lambdas/lambdat
        elseif(step.eq.2)then
          denum=pi*cofferaa_jacobian(nu,-mass2,-t,-tmax,-tmin,switch)
          if(denum.ne.0d0)then
            g=g*2d0*lambdat/denum
          else
            if(numout.lt.maxout)then
              write(nout,'(a)')' process: 1/density = 0 (width=0)'
            endif
            numout=numout+1
            switch=0
            return
          endif
        endif
      endif
      if(step.eq.1)then
        if(s.gt.0d0)then
          roots=dsqrt(s)
        else
          if(numout.lt.maxout)then
            write(nout,'(a)')' process: s<=0'
          endif
          numout=numout+1
          switch=0
          return
        endif
        p1(0)=(s+s1-s2)*0.5d0/roots
        p1(1)=0d0
        p1(2)=0d0
        p1(3)=lambdas*0.5d0/roots
        phi=dsign(phi,q1(3))
        call cofferaa_rotation(p1,phi,cost,switch)
        do i1=0,3
          q(i1)=q1(i1)+q2(i1)
          k1(i1)=q1(i1)
        enddo
        call cofferaa_boost(s,q,k1,1d0,switch)
        if(k1(1).eq.0d0)then
          phi=dsign(pi*0.5d0,k1(2))
        else
          phi=datan(k1(2)/k1(1))
        endif
        if(k1(1).lt.0d0)phi=phi+pi
        cost=k1(3)/dsqrt(k1(1)**2+k1(2)**2+k1(3)**2)
        call cofferaa_rotation(p1,-phi,cost,switch)
        call cofferaa_boost(s,q,p1,-1d0,switch)
        do i1=0,3
          p2(i1)=q(i1)-p1(i1)
          qt(i1)=q1(i1)-p1(i1)
        enddo
        t=qt(0)**2-qt(1)**2-qt(2)**2-qt(3)**2
      endif
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                c
!     1->2 particle decay without propagator                     c
!                                                                c
!     written by Markus Roth                                     c
!                                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_decay(random,q,p1,p2,g,s,s1,s2,step,switch)
      implicit none
! local variable
      real*8 pi,random(2),q(0:3),p1(0:3),p2(0:3)
      real*8 phi,cost,s,roots,lambda,denum,g,s1,s2
      integer i1,step,switch
! mcoutput
      integer nout,numout,maxout
      common/mcoutput/nout,numout,maxout
      pi=4d0*datan(1d0)
      if(switch.eq.0)return
      lambda=(s-s1-s2)*(s-s1-s2)-4d0*s1*s2
      if(lambda.gt.0d0)then
        lambda=dsqrt(lambda)
      else
        if(numout.lt.maxout)then
          write(nout,'(a)')' decay: lambda <= 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(step.eq.1)then
        if(s.gt.0d0)then
          roots=dsqrt(s)
        else
          if(numout.lt.maxout)then
            write(nout,'(a)')' decay: lambda <= 0'
            numout=numout+1
          endif
          switch=0
          return
        endif
        phi=2d0*pi*random(1)
        cost=2d0*random(2)-1d0
        p1(0)=(s+s1-s2)*0.5d0/roots
        p1(1)=0d0
        p1(2)=0d0
        p1(3)=lambda*0.5d0/roots
        call cofferaa_rotation(p1,phi,cost,switch)
        call cofferaa_boost(s,q,p1,-1d0,switch)
        do i1=0,3
          p2(i1)=q(i1)-p1(i1)
        enddo
      elseif(step.eq.2)then
        denum=pi*lambda
        if(denum.gt.0d0)then
          g=g*2d0*s/denum
        else
          if(numout.lt.maxout)then
            write(nout,'(a)')' dec: 1/density <= 0'
            numout=numout+1
          endif
          switch=0
          return
        endif
      endif
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                c
!     h function for importance sampling                         c
!                                                                c
!     written by Markus Roth                                     c
!                                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_h(random,nu,mass2,xmin,xmax,switch)
      implicit none
      real*8 cofferaa_h,random,nu,xmin,xmax,mass2,m2
      integer switch
! mctechparam
      real*8 techparam(3)
! mcoutput
      integer nout,numout,maxout
      common/mctechparam/techparam
      common/mcoutput/nout,numout,maxout
      cofferaa_h=0d0
      if(switch.eq.0)return
      m2=mass2-techparam(1)
      if(xmax-m2.lt.0d0.or.xmin-m2.lt.0d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' h: xmin-m2, xmax-m2 < 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(nu.eq.0d0)then
        cofferaa_h=random*xmax+(1d0-random)*xmin
      elseif(nu.eq.1d0)then
        cofferaa_h=dexp(random*dlog(xmax-m2)+(1d0-random)* &
                   dlog(xmin-m2))+m2
      else
        cofferaa_h=(random*(xmax-m2)**(1d0-nu) &
          +(1d0-random)*(xmin-m2)**(1d0-nu))**(1d0/(1d0-nu))+m2
      endif
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                c
!     jacobian of the h function                                 c
!                                                                c
!     written by Markus Roth                                     c
!                                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_jacobian(nu,mass2,x,xmin,xmax,switch)
      implicit none
      real*8 cofferaa_jacobian,nu,x,xmin,xmax,mass2,m2
      integer switch
! mctechparam
      real*8 techparam(3)
! mcoutput
      integer nout,numout,maxout
      common/mctechparam/techparam
      common/mcoutput/nout,numout,maxout
      cofferaa_jacobian=0d0
      if(switch.eq.0)return
      m2=mass2-techparam(1)
      if(xmax-m2.lt.0d0.or.xmin-m2.lt.0d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' jacobian: xmin-m2, xmax-m2 < 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(nu.eq.0)then
        cofferaa_jacobian=xmax-xmin
      elseif(nu.eq.1d0)then
        cofferaa_jacobian=(dlog(xmax-m2)-dlog(xmin-m2))*(x-m2)
      else
        cofferaa_jacobian=(((xmax-m2)**(1d0-nu)-(xmin-m2)**(1d0-nu)) &
          /(1d0-nu))*(x-m2)**nu
      endif
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                c
!     rotation                                                   c
!                                                                c
!     written by Markus Roth                                     c
!                                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_rotation(p,phi,cost,switch)
      implicit none
! local variable
      real*8 p(0:3),phi,cost,sint,cosp,sinp,px,py,pz
      integer switch
! mcoutput
      integer nout,numout,maxout
      common/mcoutput/nout,numout,maxout
      if(switch.eq.0)return
      cosp=dcos(phi)
      sinp=dsin(phi)
      if(dabs(cost).gt.1d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' rotation: |cost| > 1'
          numout=numout+1
        endif
        switch=0
        return
      endif
      sint=dsqrt((1d0-cost)*(1d0+cost))
      px=p(1)
      py=p(2)
      pz=p(3)
      p(1)=(pz*sint+px*cost)*cosp+py*sinp
      p(2)=-(pz*sint+px*cost)*sinp+py*cosp
      p(3)=pz*cost-px*sint
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                c
!     boost                                                      c
!                                                                c
!     written by Markus Roth                                     c
!                                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_boost(m2,q,p,dir,switch)
      implicit none
! local variable
      real*8 p(0:3),q(0:3),dir,m,m2,bx,by,bz,gamma
      real*8 p0,px,py,pz,bp,a
      integer switch
! mcoutput
      integer nout,numout,maxout
      common/mcoutput/nout,numout,maxout
      if(switch.eq.0)return
      if(m2.le.0d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' boost: m2 <= 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      m=dsqrt(m2)
      bx=-q(1)/m*dir
      by=-q(2)/m*dir
      bz=-q(3)/m*dir
      gamma=q(0)/m
      a=1d0/(1d0+gamma)
      p0=p(0)
      px=p(1)
      py=p(2)
      pz=p(3)
      bp=bx*px+by*py+bz*pz
      p(0)=gamma*p0+bp
      p(1)=px+bx*p0+a*bp*bx
      p(2)=py+by*p0+a*bp*by
      p(3)=pz+bz*p0+a*bp*bz
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     mapping from (n+1)- into n-particle phasespace              c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_mapin(kt,g,em,sp,ga,next,switch)
      implicit none
! local variables
      integer maxe,maxg,maxch
      parameter(maxe=9,maxg=4,maxch=20000)
      real*8 k(maxe,0:3),kt(maxe,0:3),g,h,denum,cofferaa_jacobian
      real*8 x,z,y,v,kk(maxe,maxe),pi,ck(0:3),ckt(0:3),st
      integer i1,i2,j,em,sp,ga,next,switch
! mcoutput
      integer nout,numout,maxout
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcoutput/nout,numout,maxout
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      if(switch.eq.0)return
      if(em.eq.0.or.sp.eq.0.or.ga.eq.0)then
        g=1d0
        return
      endif
      pi=4d0*datan(1d0)
      do i1=1,next
      do i2=0,3
        k(i1,i2)=kt(i1,i2)
      enddo
      enddo
      kk(em,sp)=k(em,0)*k(sp,0)-k(em,1)*k(sp,1)-k(em,2)*k(sp,2) &
        -k(em,3)*k(sp,3)
      kk(em,ga)=k(em,0)*k(ga,0)-k(em,1)*k(ga,1)-k(em,2)*k(ga,2) &
        -k(em,3)*k(ga,3)
      kk(sp,ga)=k(sp,0)*k(ga,0)-k(sp,1)*k(ga,1)-k(sp,2)*k(ga,2) &
        -k(sp,3)*k(ga,3)
! final-state emitter with final-state spectator
      if(em.ge.3.and.sp.ge.3)then
        x=1d0
        z=kk(em,sp)/(kk(em,sp)+kk(sp,ga))
        y=kk(em,ga)/(kk(em,sp)+kk(em,ga)+kk(sp,ga))
        do i1=0,3
          kt(em,i1)=k(em,i1)+k(ga,i1)-y/(1d0-y)*k(sp,i1)
          kt(sp,i1)=1d0/(1d0-y)*k(sp,i1)
          kt(ga,i1)=0d0
        enddo
        denum=pi*(1d0-y)*cofferaa_jacobian(powermap,0d0,1d0-y,0d0,1d0, &
          switch)*cofferaa_jacobian(powermap,0d0,1d0-z,0d0,1d0,switch)
! final-state emitter with initial-state spectator
      elseif(em.ge.3.and.sp.le.2)then
        x=(kk(em,sp)+kk(sp,ga)-kk(em,ga))/(kk(em,sp)+kk(sp,ga))
        z=kk(em,sp)/(kk(em,sp)+kk(sp,ga))
        do i1=0,3
          kt(em,i1)=k(em,i1)+k(ga,i1)-(1d0-x)*k(sp,i1)
          kt(sp,i1)=x*k(sp,i1)
          kt(ga,i1)=0d0
        enddo
        denum=pi*cofferaa_jacobian(powermap,0d0,1d0-x,0d0,1d0,switch) &
          *cofferaa_jacobian(powermap,0d0,1d0-z,0d0,1d0,switch)
! initial-state emitter with final-state spectator
      elseif(em.le.2.and.sp.ge.3)then
        x=(kk(em,sp)+kk(em,ga)-kk(sp,ga))/(kk(em,sp)+kk(em,ga))
        z=kk(em,sp)/(kk(em,sp)+kk(em,ga))
        do i1=0,3
          kt(sp,i1)=k(sp,i1)+k(ga,i1)-(1d0-x)*k(em,i1)
          kt(em,i1)=x*k(em,i1)
          kt(ga,i1)=0d0
        enddo
        denum=pi*cofferaa_jacobian(powermap,0d0,1d0-x,0d0,1d0,switch) &
          *cofferaa_jacobian(powermap,0d0,1d0-z,0d0,1d0,switch)
! initial-state emitter with initial-state spectator
      elseif(em.le.2.and.sp.le.2)then
        x=(kk(em,sp)-kk(em,ga)-kk(sp,ga))/kk(em,sp)
        v=kk(em,ga)/kk(em,sp)
        do i1=0,3
          ck(i1)=k(em,i1)+k(sp,i1)-k(ga,i1)
          ckt(i1)=x*k(em,i1)+k(sp,i1)
        enddo
        do j=3,6
          kk(em,j)=k(em,0)*k(j,0)-k(em,1)*k(j,1)-k(em,2)*k(j,2) &
            -k(em,3)*k(j,3)
          kk(sp,j)=k(sp,0)*k(j,0)-k(sp,1)*k(j,1)-k(sp,2)*k(j,2) &
            -k(sp,3)*k(j,3)
          kk(ga,j)=k(ga,0)*k(j,0)-k(ga,1)*k(j,1)-k(ga,2)*k(j,2) &
            -k(ga,3)*k(j,3)
          do i1=0,3
            kt(j,i1)=k(j,i1) &
              -((1d0+x)*kk(em,j)+2d0*kk(sp,j)-kk(ga,j)) &
                /(4d0*x+v-x*v)/kk(em,sp)*(ck(i1)+ckt(i1)) &
              +(kk(em,j)+kk(sp,j)-kk(ga,j))/x/kk(em,sp)*ckt(i1)
          enddo
        enddo
        do i1=0,3
          kt(em,i1)=x*k(em,i1)
          kt(sp,i1)=k(sp,i1)
          kt(ga,i1)=0d0
        enddo
        denum=pi*(1d0-x)*cofferaa_jacobian(powermap,0d0,1d0-x,0d0,1d0, &
          switch)
      endif
      st=kt(em,0)*kt(sp,0)-kt(em,1)*kt(sp,1)-kt(em,2)*kt(sp,2) &
        -kt(em,3)*kt(sp,3)
      if(denum.ge.0d0.and.st.ne.0d0)then
        g=x/denum/st
      else
        if(numout.lt.maxout)then
          write(nout,'(a21,i5,i9)')' mapin: 1/g <= 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     mapping from n- into (n+1)-particle phase space             c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_mapout(random,k,x,g,em,sp,ga,next,switch)
      implicit none
! local variables
      integer maxch,maxg,maxe
      parameter(maxch=20000,maxg=4,maxe=9)
      real*8 random(3),k(maxe,0:3),kt(maxe,0:3),g,m2,m,pi
      real*8 kem(0:3),ksp(0:3),ktr(0:3),ck(0:3),ckt(0:3)
      real*8 ksum(0:3),kj(maxe,0:3),y,cofferaa_h,z,f,v,c,x, &
       kk(maxe,maxe)
      integer i1,i2,j,em,sp,ga,next,switch
! mcoutput
      integer nout,numout,maxout
! mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcoutput/nout,numout,maxout
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      if(switch.eq.0)return
      if(em.eq.0.and.sp.eq.0.and.ga.eq.0)return
      pi=4d0*datan(1d0)
! mapping from n- into (n+1)-particle phasespace
      do i1=0,3
        kem(i1)=k(em,i1)
        ksp(i1)=k(sp,i1)
        ksum(i1)=kem(i1)+ksp(i1)
      enddo
      m2=ksum(0)**2-ksum(1)**2-ksum(2)**2-ksum(3)**2
      if(m2.le.0d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' mapout: m2 <= 0'
          numout=numout+1
        endif
        switch=0
        return
      else
        m=dsqrt(m2)
      endif
      if(em.ge.3.and.sp.ge.3)then
! final-state emitter with final-state spectator
        y=1d0-cofferaa_h(random(1),powermap,0d0,0d0,1d0,switch)
        z=1d0-cofferaa_h(random(2),powermap,0d0,0d0,1d0,switch)
        f=dsqrt(y*z/(1d0-z))
        ktr(0)=-f*f*m
        ktr(1)=f*m*dcos(2d0*pi*random(3))
        ktr(2)=f*m*dsin(2d0*pi*random(3))
        ktr(3)=-f*f*m
        call cofferaa_transverse(m2,ksum,ksp,ktr,switch)
        do i1=0,3
          k(em,i1)=y*(1d0+z)*ksp(i1)+z*kem(i1)+(1d0-z)*ktr(i1)
          k(ga,i1)=(1d0-z)*kem(i1)-y*z*ksp(i1)-(1d0-z)*ktr(i1)
          k(sp,i1)=(1d0-y)*ksp(i1)
        enddo
! final-state emitter with initial-state spectator
      elseif(em.ge.3.and.sp.le.2)then
        z=1d0-cofferaa_h(random(2),powermap,0d0,0d0,1d0,switch)
        f=dsqrt(z*(1d0-x)/(1d0-z)/x)
        ktr(0)=-f*f*m
        ktr(1)=f*m*dcos(2d0*pi*random(3))
        ktr(2)=f*m*dsin(2d0*pi*random(3))
        ktr(3)=-f*f*m
        call cofferaa_transverse(m2,ksum,ksp,ktr,switch)
        do i1=0,3
          k(ga,i1)=(1-z)*kem(i1)-(1d0-x)/x*z*ksp(i1)-(1d0-z)*ktr(i1)
          k(em,i1)=(1d0-x)/x*(1d0+z)*ksp(i1)+z*kem(i1) &
            +(1d0-z)*ktr(i1)
          k(sp,i1)=ksp(i1)/x
        enddo
! initial-state emitter with final-state spectator
      elseif(em.le.2.and.sp.ge.3)then
        z=1d0-cofferaa_h(random(2),powermap,0d0,0d0,1d0,switch)
        f=dsqrt((1d0-z)*(1d0-x)/z/x)
        ktr(0)=-f*f*m
        ktr(1)=f*m*dcos(2d0*pi*random(3))
        ktr(2)=f*m*dsin(2d0*pi*random(3))
        ktr(3)=-f*f*m
        call cofferaa_transverse(m2,ksum,kem,ktr,switch)
        do i1=0,3
          k(ga,i1)=(1d0-x)/x*(2d0-z)*kem(i1)+(1d0-z)*ksp(i1) &
            +z*ktr(i1)
          k(sp,i1)=z*ksp(i1)-(1d0-x)/x*(1d0-z)*kem(i1)-z*ktr(i1)
          k(em,i1)=kem(i1)/x
        enddo
! initial-state emitter with initial-state spectator
      elseif(em.le.2.and.sp.le.2)then
        do i1=3,next
        do i2=0,3
          kt(i1,i2)=k(i1,i2)
        enddo
        enddo
        v=random(2)*(1d0-x)
        f=dsqrt(v*(1d0-x-v)/x)
        ktr(0)=0d0
        ktr(1)=f*m*dcos(2d0*pi*random(3))
        ktr(2)=f*m*dsin(2d0*pi*random(3))
        ktr(3)=0d0
        call cofferaa_transverse(m2,ksum,kem,ktr,switch)
        do i1=0,3
          ck(i1)=(x+v)/x*kem(i1)+(1d0-v)*ksp(i1)-ktr(i1)
          ckt(i1)=kem(i1)+ksp(i1)
        enddo
        do i1=0,3
          k(ga,i1)=v*ksp(i1)+(1d0-x-v)/x*kem(i1)+ktr(i1)
          k(em,i1)=kem(i1)/x
        enddo
        do j=3,next
        if(j.ne.ga)then
          c=ktr(0)*kt(j,0)-ktr(1)*kt(j,1)-ktr(2)*kt(j,2) &
            -ktr(3)*kt(j,3)
          kk(em,sp)=kem(0)*ksp(0)-kem(1)*ksp(1)-kem(2)*ksp(2) &
            -kem(3)*ksp(3)
          kk(em,j)=kt(j,0)*kem(0)-kt(j,1)*kem(1)-kt(j,2)*kem(2) &
            -kt(j,3)*kem(3)
          kk(sp,j)=kt(j,0)*ksp(0)-kt(j,1)*ksp(1)-kt(j,2)*ksp(2) &
            -kt(j,3)*ksp(3)
          do i1=0,3
            k(j,i1)=kt(j,i1) &
              -((2d0*x+v)*kk(em,j)+(2d0-v)*x*kk(sp,j)-x*c) &
               /(4d0*x+v-x*v)/kk(em,sp)*(ck(i1)+ckt(i1)) &
              +(kk(em,j)+kk(sp,j))/kk(em,sp)*ck(i1)
          enddo
        endif
        enddo
      endif
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     determination of transverse momenta                         c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cofferaa_transverse(m2,ksum,k,ktr,switch)
      implicit none
! local variables
      real*8 k(0:3),ktr(0:3),ksum(0:3),p(0:3),pvec,pi,phi,cost,m2
      integer i1,switch
! output
      integer nout,numout,maxout
      common/output/nout,numout,maxout
      if(switch.eq.0)return
      pi=4d0*datan(1d0)
      do i1=0,3
        p(i1)=k(i1)
      enddo
      call cofferaa_boost(m2,ksum,p,1d0,switch)
      if(p(1).eq.0d0)then
        phi=dsign(pi*0.5d0,p(2))
      else
        phi=datan(p(2)/p(1))
      endif
      if(p(1).lt.0d0)phi=phi+pi
      if(p(1).gt.0d0.and.p(2).lt.0d0)phi=phi+2d0*pi
      pvec=dsqrt(p(1)**2+p(2)**2+p(3)**2)
      if(pvec.eq.0d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' transverse: pvec = 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      cost=p(3)/pvec
      call cofferaa_rotation(ktr,-phi,cost,switch)
      call cofferaa_boost(m2,ksum,ktr,-1d0,switch)
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     3-particle vertex                                           c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_vertex(i1,i2,i3,schannel,smself)
      implicit none
! local variables
      integer maxv
      parameter(maxv=200)
      integer i1,i2,i3,smself
      logical cofferaa_vertex,cofferaa_vertexg,schannel
! mcparticle
      integer family(-maxv:maxv,6),light(-maxv:maxv)
      character*3 pname(-maxv:maxv),gname(-maxv:maxv)
      common/mcparticle/family,light,pname,gname
      cofferaa_vertex=.false.
      if(gname(i1).eq.'   ')return
      if(gname(i2).eq.'   ')return
      if(gname(i3).eq.'   ')return
! family conservation
      if(family(i1,1)+family(i2,1)+family(i3,1).ne.0)return
      if(family(i1,2)+family(i2,2)+family(i3,2).ne.0)return
      if(family(i1,3)+family(i2,3)+family(i3,3).ne.0)return
      if(family(i1,4)+family(i2,4)+family(i3,4).ne.0)return
      if(family(i1,5)+family(i2,5)+family(i3,5).ne.0)return
      if(family(i1,6)+family(i2,6)+family(i3,6).ne.0)return
      cofferaa_vertex=.true.
! testing vertex
      if(cofferaa_vertexg(i1,i2,i3,schannel,smself))return
      if(cofferaa_vertexg(i1,i3,i2,schannel,smself))return
      if(cofferaa_vertexg(i2,i1,i3,schannel,smself))return
      if(cofferaa_vertexg(i2,i3,i1,schannel,smself))return
      if(cofferaa_vertexg(i3,i1,i2,schannel,smself))return
      if(cofferaa_vertexg(i3,i2,i1,schannel,smself))return
      cofferaa_vertex=.false.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     generic 3-particle vertex                                   c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_vertexg(i1,i2,i3,schannel,smself)
      implicit none
! local variables
      integer maxv
      parameter(maxv=200)
      integer i1,i2,i3,smself
      character*3 p1,p2,p3,g1,g2,g3
      logical cofferaa_vertexg,nonstandardcoup,schannel
! mcparticle
      integer family(-maxv:maxv,6),light(-maxv:maxv)
      character*3 pname(-maxv:maxv),gname(-maxv:maxv)
      common/mcparticle/family,light,pname,gname
      cofferaa_vertexg=.true.
      p1=pname(i1)
      p2=pname(i2)
      p3=pname(i3)
      g1=gname(i1)
      g2=gname(i2)
      g3=gname(i3)
! 2 leptons - higgs
      if(light(i2).eq.0)then
        if(g1.eq.'el '.and.g2.eq.'el~'.and.g3.eq.'h0 ')return
      endif
! 2 quarks - higgs
      if(light(i1).eq.0.and.light(i2).eq.0)then
        if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'h0 ')return
        if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'h0 ')return
      endif
! 2 leptons - gauge boson
      if(g1.eq.'ne '.and.g2.eq.'el~'.and.g3.eq.'W- ')return
      if(g1.eq.'ne~'.and.g2.eq.'el '.and.g3.eq.'W-~')return
! 2 quarks - gauge boson
      if(g1.eq.'uq '.and.g2.eq.'dq~'.and.g3.eq.'W- ')return
      if(g1.eq.'uq~'.and.g2.eq.'dq '.and.g3.eq.'W-~')return
! 2 quarks - gluon
      if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'gl ')return
      if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'gl ')return

      if(smself.eq.12)then
! 2 leptons - gauge boson
      if(g1.eq.'el '.and.g2.eq.'el~'.and.g3.eq.'ph ')return
      if(g1.eq.'ne '.and.g2.eq.'ne~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'el '.and.g2.eq.'el~'.and.g3.eq.'Z0 ')return
! 2 quarks - gauge boson
      if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'ph ')return
      if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'ph ')return
      if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'Z0 ')return
! 3 higgs
      if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.'h0 ')return
! higgs - 2 gauge bosons
      if(g1.eq.'h0 '.and.g2.eq.'Z0 '.and.g3.eq.'Z0 ')return
      if(g1.eq.'h0 '.and.g2.eq.'W- '.and.g3.eq.'W-~')return
! 3 gauge bosons
      if(g1.eq.'ph '.and.g2.eq.'W- '.and.g3.eq.'W-~')return
      if(g1.eq.'Z0 '.and.g2.eq.'W- '.and.g3.eq.'W-~')return
! 4 particle vertices
      if(g1.eq.'gl '.and.g2.eq.'ph '.and.g3.eq.' 1~')return
      if(g1.eq.'gl '.and.g2.eq.'Z0 '.and.g3.eq.' 1~')return
      if(g1.eq.'gl '.and.g2.eq.'W- '.and.g3.eq.' 2~')return
      if(g1.eq.'gl '.and.g2.eq.'W-~'.and.g3.eq.' 3 ')return
      if(schannel)then
        if(g1.eq.'gl '.and.g2.eq.'gl '.and.g3.eq.' 4 ')return
        if(g1.eq.'gl '.and.g2.eq.'gl '.and.g3.eq.' 4~')return
      endif
      if(g1.eq.'ph '.and.g2.eq.'W-~'.and.g3.eq.' 5 ')return
      if(g1.eq.'Z0 '.and.g2.eq.'W-~'.and.g3.eq.' 5 ')return
      if(g1.eq.'ph '.and.g2.eq.'W- '.and.g3.eq.' 5~')return
      if(g1.eq.'Z0 '.and.g2.eq.'W- '.and.g3.eq.' 5~')return
      if(g1.eq.'ph '.and.g2.eq.'W- '.and.g3.eq.' 6~')return
      if(g1.eq.'Z0 '.and.g2.eq.'W- '.and.g3.eq.' 6~')return
      if(g1.eq.'ph~'.and.g2.eq.'W-~'.and.g3.eq.' 6 ')return
      if(g1.eq.'Z0 '.and.g2.eq.'W-~'.and.g3.eq.' 6 ')return
      if(g1.eq.'W+ '.and.g2.eq.'W- '.and.g3.eq.' 7 ')return
      if(g1.eq.'W+ '.and.g2.eq.'W- '.and.g3.eq.' 7~')return
      if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.' 7~')return
      if(g1.eq.'ph '.and.g2.eq.'ph '.and.g3.eq.' 7~')return
      if(g1.eq.'ph '.and.g2.eq.'Z0 '.and.g3.eq.' 7~')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.' 7~')return
      if(schannel)then
        if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.' 8 ')return
        if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.' 8~')return
      endif
      if(g1.eq.'ph '.and.g2.eq.'ph '.and.g3.eq.'11 ')return
      if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.'11 ')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'12~')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'13~')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'14~')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'15~')return
      endif

! 3 gluons
      if(g1.eq.'gl '.and.g2.eq.'gl '.and.g3.eq.'gl ')return
! 4 particle vertices
      if(g1.eq.'gl '.and.g2.eq.'gl '.and.g3.eq.' 1~')return
! add non-standard couplings (with family conservation!)
!      if(nonstandardcoup(p1,p2,p3,g1,g2,g3,schannel))return
      cofferaa_vertexg=.false.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     identifying particles (usually not used)                    c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function cofferaa_pid(binary)
      implicit none
! local variables
      integer cofferaa_pid,binary,i1,i2,i3,nexternal
      i3=1
      nexternal=30
      cofferaa_pid=0
      do i1=nexternal,1,-1
        i2=binary/2**(i1-1)
        if((i2/2)*2.ne.i2)then
          cofferaa_pid=cofferaa_pid+i1*i3
          i3=10*i3
        endif
      enddo
      end
