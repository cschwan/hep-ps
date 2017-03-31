ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     Monte Carlo generator                                       c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine generation(random,kbeam,k,g,channel,generator,
     *  switch)
      implicit none
c local variables 
      integer maxe,maxch,maxg,maxv
      parameter(maxe=7,maxch=10000,maxg=4,maxv=200)
      real*8 random(3*maxe-7),p(0:3,2**maxe),s(2**maxe)
      real*8 k(maxe,0:3),kt(maxe,0:3),kbeam(2,0:3),h
      real*8 mmin,mmax,smin,smax,x,g,tcut1,tcut2
      integer i1,i2,i3,i4,ns,nt,channel,generator,switch,step
      integer ranstart,virtinv,inv1,inv2,em,sp,ga,pid
c mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
c mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
c mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
c mcoutput
      integer nout,numout,maxout
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv,
     *  lmin,lmax
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      common/mcoutput/nout,numout,maxout
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      if(switch.eq.0)return
c      step=2 ! calculation of local density
      step=1 ! calculation of momenta
      em=emmap(channel,generator)
      sp=spmap(channel,generator)
      ga=gamap(channel,generator)
      if(step.eq.1)then
c incoming momenta and definition of x        
        ranstart=1+ninv(channel,generator)
     *    +2*nprocess(channel,generator)+2*ndecay(channel,generator)
        if(em.eq.1.or.(sp.eq.1.and.em.ge.3))then
          x=1d0-h(random(ranstart),powermap,0d0,0d0,1d0,switch)
          do i1=0,3
            p(i1,1)=-x*kbeam(1,i1)
            p(i1,2)=-kbeam(2,i1)
          enddo
        elseif(em.eq.2.or.(sp.eq.2.and.em.ge.3))then
          x=1d0-h(random(ranstart),powermap,0d0,0d0,1d0,switch)
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
c remaining momenta 
        do i1=0,3
          p(i1,3)=p(i1,1)+p(i1,2)
          p(i1,allbinary(ga,generator)-1)=-p(i1,1)
          p(i1,allbinary(ga,generator)-2)=-p(i1,2)
          p(i1,allbinary(ga,generator)-3)=-p(i1,3)
        enddo
c square of center-of-mass energy
        s(3)=p(0,3)**2-p(1,3)**2-p(2,3)**2-p(3,3)**2
        s(allbinary(ga,generator)-3)=s(3)
c event excluded by cuts after intial-state radiation
        if(s(3).le.mcutinv(allbinary(ga,generator)-3,generator)**2)then
          if(numout.lt.maxout)then
            write(nout,'(a)')' generation: s < smin'
            numout=numout+1
          endif
          switch=0
          return
        endif
      else
c calculation of local density
        g=1d0
c defining momenta
        do i1=1,nexternal(generator)
        do i2=0,3
          kt(i1,i2)=k(i1,i2)
        enddo
        enddo
        call mapin(kt,g,em,sp,ga,nexternal(generator),switch)
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
c calculating invariants
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
c external masses
      do i1=1,nexternal(generator)
      if(i1.ne.ga)then
        s(2**(i1-1))=massext2(i1,generator)
        s(allbinary(ga,generator)-2**(i1-1))=massext2(i1,generator)
      endif
      enddo
c inv
      do ns=1,ninv(channel,generator)
        inv1=ininv(ns,channel,generator)
        inv2=allbinary(ga,generator)-3-ininv(ns,channel,generator)
        mmin=mcutinv(inv1,generator)
        mmax=mcutinv(inv2,generator)
        do i1=1,ns-1
          virtinv=ininv(i1,channel,generator)
          if(lmin(i1,ns,channel,generator).and.
     *      dsqrt(s(virtinv)).gt.mcutinv(virtinv,generator))then
            mmin=mmin+dsqrt(s(virtinv))-mcutinv(inv1,generator)
     *        +mcutinv(inv1-virtinv,generator)
            inv1=inv1-virtinv
          endif
          if(lmax(i1,ns,channel,generator).and.
     *      dsqrt(s(virtinv)).gt.mcutinv(virtinv,generator))then
            mmax=mmax+dsqrt(s(virtinv))-mcutinv(inv2,generator)
     *        +mcutinv(inv2-virtinv,generator)
            inv2=inv2-virtinv
          endif
        enddo        
        smin=mmin**2
        smax=(dsqrt(s(3))-mmax)**2
        call inv(
     *    random(ns),                                    ! random number
     *    s(ininv(ns,channel,generator)),                ! invariant mass
     *    g,                                             ! local density
     *    mass(idhepinv(ns,channel,generator)),          ! mass
     *    width(idhepinv(ns,channel,generator)),         ! width
     *    powerinv(ns,channel,generator),                ! mapping 
     *    smin,                                          ! lower bound 
     *    smax,                                          ! upper bound
     *    step,switch)
      enddo
c process
      ranstart=ninv(channel,generator)
      do nt=1,nprocess(channel,generator)
        tcut1=0d0
        tcut2=0d0
        if(nt.eq.1)then
          tcut1=tcutprocess(allbinary(ga,generator)
     *      -virtprocess(nt,channel,generator),generator)
          tcut2=tcutprocess(virtprocess(nt,channel,generator),
     *      generator)
        endif
        call process(
     *    random(ranstart+2*nt-1),                       ! random numbers
     *    p(0,in1process(nt,channel,generator)),         ! incomming particle 1
     *    p(0,in2process(nt,channel,generator)),         ! incomming particle 2
     *    p(0,out1process(nt,channel,generator)),        ! outgoing particle 1
     *    p(0,out2process(nt,channel,generator)),        ! outgoing particle 2
     *    p(0,virtprocess(nt,channel,generator)),        ! virtual particle
     *    g,                                             ! local density
     *    mass(idhepprocess(nt,channel,generator)),      ! mass
     *    width(idhepprocess(nt,channel,generator)),     ! width
     *    powerprocess(nt,channel,generator),            ! mapping 
     *    s(inprocess(nt,channel,generator)),            ! incomming particles 
     *    s(out1process(nt,channel,generator)),          ! outgoing particle 1
     *    s(out2process(nt,channel,generator)),          ! outgoing partcile 2
     *    s(virtprocess(nt,channel,generator)),          ! virtual particle
     *    s(in1process(nt,channel,generator)),           ! incomming particle 1
     *    s(in2process(nt,channel,generator)),           ! incomming particle 2
     *    tcut1,tcut2,                                   ! angular cuts
     *    step,switch)
      enddo
c decay
      ranstart=ninv(channel,generator)+2*nprocess(channel,generator)
      do ns=1,ndecay(channel,generator)
        call decay(
     *    random(ranstart+2*ns-1),                       ! random numbers
     *    p(0,indecay(ns,channel,generator)),            ! incomming particle 
     *    p(0,out1decay(ns,channel,generator)),          ! outgoing particle 1
     *    p(0,out2decay(ns,channel,generator)),          ! outgoing particle 2
     *    g,                                             ! local density
     *    s(indecay(ns,channel,generator)),              ! incomming particle
     *    s(out1decay(ns,channel,generator)),            ! outgoing particle 1
     *    s(out2decay(ns,channel,generator)),            ! outgoing particle 2
     *    step,switch)
      enddo
c defining external momenta
      do i1=0,3
        k(1,i1)=-p(i1,1)
        k(2,i1)=-p(i1,2)
        do i2=3,nexternal(generator)
          k(i2,i1)=p(i1,2**(i2-1))
        enddo
      enddo
c map
      if(step.eq.1)then
        ranstart=1+ninv(channel,generator)
     *    +2*nprocess(channel,generator)+2*ndecay(channel,generator)
        call mapout(random(ranstart),k,x,g,em,sp,ga,
     *    nexternal(generator),switch)
      endif
c restoring external masses
      do i1=1,nexternal(generator)
        k(i1,0)=dsqrt(massext2(i1,generator)
     *    +k(i1,1)**2+k(i1,2)**2+k(i1,3)**2)
      enddo
c      call checkmom(k,nexternal(generator),switch)
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     calculating densities                                       c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine checkmom(k,next,switch)
      implicit none
c local variables
      integer maxe,maxg
      parameter(maxe=7,maxg=4)
      real*8 k(maxe,0:3),mom,acc
      integer i1,i2,next,switch
c mckinemetics
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
          write(*,'(" momentum conservation violated:",i2,1d24.16)')
     *      i1,mom
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     calculating denities                                        c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine density(k,g,generator,switch)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=7,maxch=10000,maxg=4,maxv=200)
      real*8 random(3*maxe-7),p(0:3,2**maxe,0:maxe,0:maxe,0:maxe)
      real*8 s(2**maxe,0:maxe,0:maxe,0:maxe),k(maxe,0:3),kt(maxe,0:3)
      real*8 gmap(0:maxe,0:maxe,0:maxe),ginv(maxch),gprocess(maxch)
      real*8 gdecay(maxch),g(maxch),mmin,mmax,smin,smax,tcut1,tcut2
      integer i1,i2,i3,i4,ns,nt,channel,generator,switch,step,pid
      integer virtinv,inv1,inv2,em,sp,ga
c mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
c mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
c mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
c mcdensity
      integer nsinv(maxch,maxg),chinv(maxch,maxg),maxinv(maxg)
      integer ntprocess(maxch,maxg),chprocess(maxch,maxg)
      integer maxprocess(maxg),nsdecay(maxch,maxg),chdecay(maxch,maxg)
      integer maxdecay(maxg),numinv(maxe,maxch,maxg)
      integer numprocess(maxe,maxch,maxg),numdecay(maxe,maxch,maxg)
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv,
     *  lmin,lmax
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      common/mcdensity/nsinv,chinv,maxinv,ntprocess,chprocess,
     *  maxprocess,nsdecay,chdecay,maxdecay,numinv,numprocess,numdecay
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
c          call checkmom(k,nexternal(generator),switch)
c generating phasespace configurations
          call mapin(kt,gmap(em,sp,ga),em,sp,ga,nexternal(generator),
     *      switch)     
c          if(ga.gt.0)call checkmom(kt,nexternal(generator)-1,switch)
c calculating invariants
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
c calculating invariants
          do i1=1,allbinary(ga,generator)-1
            if(p(0,allbinary(ga,generator)-i1,em,sp,ga).eq.0d0)then
              do i2=1,nexternal(generator)
                i3=i1/2**(i2-1)
                if(2*(i3/2).ne.i3.and.i1.ne.2**(i2-1))then
                  do i4=0,3
                    p(i4,i1,em,sp,ga)=p(i4,i1-2**(i2-1),em,sp,ga)
     *                +p(i4,2**(i2-1),em,sp,ga)
                  enddo
                  goto 100
                endif
              enddo
 100          s(i1,em,sp,ga)=p(0,i1,em,sp,ga)**2-p(1,i1,em,sp,ga)**2
     *          -p(2,i1,em,sp,ga)**2-p(3,i1,em,sp,ga)**2
            else
              do i2=0,3
                p(i2,i1,em,sp,ga)=
     *            -p(i2,allbinary(ga,generator)-i1,em,sp,ga)
              enddo
              s(i1,em,sp,ga)=s(allbinary(ga,generator)-i1,em,sp,ga)
            endif
          enddo
c external masses
          do i1=1,nexternal(generator)
          if(i1.ne.ga)then
            s(2**(i1-1),em,sp,ga)=massext2(i1,generator)
            s(allbinary(ga,generator)-2**(i1-1),em,sp,ga)=
     *        massext2(i1,generator)
          endif
          enddo
        endif
      enddo
      enddo
      enddo
c inv
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
          if(lmin(i2,ns,channel,generator).and.
     *      dsqrt(s(virtinv,em,sp,ga)).gt.
     *      mcutinv(virtinv,generator))then
            mmin=mmin+dsqrt(s(virtinv,em,sp,ga))
     *        -mcutinv(inv1,generator)
     *        +mcutinv(inv1-virtinv,generator)
            inv1=inv1-virtinv
          endif
          if(lmax(i2,ns,channel,generator).and.
     *      dsqrt(s(virtinv,em,sp,ga)).gt.
     *      mcutinv(virtinv,generator))then
            mmax=mmax+dsqrt(s(virtinv,em,sp,ga))
     *        -mcutinv(inv2,generator)
     *        +mcutinv(inv2-virtinv,generator)
            inv2=inv2-virtinv
          endif
        enddo        
        smin=mmin**2
        smax=(dsqrt(s(3,em,sp,ga))-mmax)**2
        call inv(
     *    random,                                           ! random number
     *    s(ininv(ns,channel,generator),em,sp,ga),          ! invariant mass
     *    ginv(i1),                                         ! local density
     *    mass(idhepinv(ns,channel,generator)),             ! mass
     *    width(idhepinv(ns,channel,generator)),            ! width
     *    powerinv(ns,channel,generator),                   ! mapping 
     *    smin,                                             ! lower bound 
     *    smax,                                             ! upper bound
     *    step,switch)
      enddo
c process
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
          tcut1=tcutprocess(allbinary(ga,generator)
     *      -virtprocess(nt,channel,generator),generator)
          tcut2=tcutprocess(virtprocess(nt,channel,generator),
     *      generator)
        endif
        call process(
     *    random,                                           ! random numbers
     *    p(0,in1process(nt,channel,generator),em,sp,ga),   ! incomming p. 1
     *    p(0,in2process(nt,channel,generator),em,sp,ga),   ! incomming p. 2
     *    p(0,out1process(nt,channel,generator),em,sp,ga),  ! outgoing p. 1
     *    p(0,out2process(nt,channel,generator),em,sp,ga),  ! outgoing p. 2
     *    p(0,virtprocess(nt,channel,generator),em,sp,ga),  ! virtual p.
     *    gprocess(i1),                                     ! local density
     *    mass(idhepprocess(nt,channel,generator)),         ! mass
     *    width(idhepprocess(nt,channel,generator)),        ! width
     *    powerprocess(nt,channel,generator),               ! mapping 
     *    s(inprocess(nt,channel,generator),em,sp,ga),      ! incomming p.s 
     *    s(out1process(nt,channel,generator),em,sp,ga),    ! outgoing p. 1
     *    s(out2process(nt,channel,generator),em,sp,ga),    ! outgoing p. 2
     *    s(virtprocess(nt,channel,generator),em,sp,ga),    ! virtual p.
     *    s(in1process(nt,channel,generator),em,sp,ga),     ! incomming p. 1
     *    s(in2process(nt,channel,generator),em,sp,ga),     ! incomming p. 2
     *    tcut1,tcut2,                                      ! angular cuts
     *    step,switch)
      enddo
c decay
      do i1=1,maxdecay(generator)
        gdecay(i1)=1d0
        channel=chdecay(i1,generator)
        em=emmap(channel,generator)
        sp=spmap(channel,generator)
        ga=gamap(channel,generator)
        ns=nsdecay(i1,generator)
        call decay(
     *    random,                                           ! random numbers
     *    p(0,indecay(ns,channel,generator),em,sp,ga),      ! incomming p. 
     *    p(0,out1decay(ns,channel,generator),em,sp,ga),    ! outgoing p. 1
     *    p(0,out2decay(ns,channel,generator),em,sp,ga),    ! outgoing p. 2
     *    gdecay(i1),                                       ! local density
     *    s(indecay(ns,channel,generator),em,sp,ga),        ! incomming p.
     *    s(out1decay(ns,channel,generator),em,sp,ga),      ! outgoing p. 1
     *    s(out2decay(ns,channel,generator),em,sp,ga),      ! outgoing p. 2
     *    step,switch)
       enddo
c calculating densities        
      do channel=1,nchannel(generator)
        em=emmap(channel,generator)
        sp=spmap(channel,generator)
        ga=gamap(channel,generator)
        g(channel)=gmap(em,sp,ga)
c inv
        do ns=1,ninv(channel,generator)
          g(channel)=g(channel)
     *      *ginv(numinv(ns,channel,generator))
        enddo
c process
        do nt=1,nprocess(channel,generator)
          g(channel)=g(channel)
     *      *gprocess(numprocess(nt,channel,generator))
        enddo
c decay
        do ns=1,ndecay(channel,generator)
          g(channel)=g(channel)
     *      *gdecay(numdecay(ns,channel,generator))
        enddo
      enddo
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialization generator                                    c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initgenerator(energy,smin,hepnum,generator,next,
     *  smodel,sincludecuts,ssub)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv,maxo
      parameter(maxe=7,maxch=10000,maxg=4,maxv=200,maxo=20)
      real*8 m2,m3,e2,e3,scutinv,energy,smin
      integer idhep(maxe,maxe),binary(maxe,maxe),idhep2,idhep3
      integer i1,i2,i3,i4,i5,i6,ns,nt,maxns,maxnt,binary1,binary2
      integer binary3,in1(maxe),in2(maxe),out1(maxe),out2(maxe),pid
      integer virt(maxe),channel,generator,sincludecuts,prop2,prop3
      integer smodel,smap,noutgen,hepnum(maxe),naux,nmap,next,ssub
      integer em,sp,ga,n
      logical vertex,included,exist,comparechannel,compareinv
      logical comparedecay,compareprocess,schannel
      character*9 particle(maxv)
      character*80 vertices
c mcparticle
      integer family(-maxv:maxv,6),light(-maxv:maxv)
      character*3 pname(-maxv:maxv),gname(-maxv:maxv)
c mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
c mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
c mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
c mcdensity
      integer nsinv(maxch,maxg),chinv(maxch,maxg),maxinv(maxg)
      integer ntprocess(maxch,maxg),chprocess(maxch,maxg)
      integer maxprocess(maxg),nsdecay(maxch,maxg),chdecay(maxch,maxg)
      integer maxdecay(maxg),numinv(maxe,maxch,maxg)
      integer numprocess(maxe,maxch,maxg),numdecay(maxe,maxch,maxg)
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
c mcoutput
      integer nout,numout,maxout
c mctechparam
      real*8 techparam(3)
c mcadaptopt
      real*8 alphaopt(maxch,maxo,maxg),betaopt(0:maxch,maxg)
      real*8 wi(maxch,maxg),alphamin
      integer nopt(0:maxo,maxg),opt(maxg)
c mccuts
      real*8 ecutp,ecutl,ecutq,scutqq,ccutpb,ccutpl,ccutpq,ccutll
      real*8 ccutqq,ccutlq,ccutlb,ccutqb,ecut(maxe,maxg)
      real*8 scut(maxe,maxe,maxg),ccut(maxe,maxe,maxg),xmin(maxg)
      common/mcparticle/family,light,pname,gname
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv,
     *  lmin,lmax
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      common/mcdensity/nsinv,chinv,maxinv,ntprocess,chprocess,
     *  maxprocess,nsdecay,chdecay,maxdecay,numinv,numprocess,numdecay
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      common/mcoutput/nout,numout,maxout
      common/mctechparam/techparam
      common/mcadaptopt/alphaopt,betaopt,wi,alphamin,nopt,opt
      common/mccuts/ecutp,ecutl,ecutq,scutqq,ccutpb,ccutpl,ccutpq,
     *  ccutll,ccutqq,ccutlq,ccutlb,ccutqb,ecut,scut,ccut,xmin
c options
      noutgen=1                     ! output for generator
      smap=1                        ! mapping
c technical parameter in h function and subroutine process
      techparam(1)=1d-4             ! small negative mass in h-function
      techparam(2)=1d-10*energy**2  ! allowed uncertainty tmax in process
      techparam(3)=1d-12*energy     ! allowed uncertainty q1,q2 in process
      powermap=0.8d0
c number of external particles
      nexternal(generator)=next
c initialize binary counting of particle number
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize particles                                        c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c generic particle names, see also subroutine vertexg
      do i1=1,maxv
        particle(i1)=' '
      enddo
c name, generic name, anti-field, family 
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
c supersymmetric particles
      particle(31)='D1 DQ y 1'
      particle(32)='U1 UQ y 1'
      particle(33)='S1 DQ y 2'
      particle(34)='C1 UQ y 2'
      particle(35)='B1 DQ y 3'
      particle(36)='T1 UQ y 3'
      particle(41)='D2 DQ y 1'
      particle(42)='U2 UQ y 1'
      particle(43)='S2 DQ y 2'
      particle(44)='C2 UQ y 2'
      particle(45)='B2 DQ y 3'
      particle(46)='T2 UQ y 3'
      particle(51)='E1 EL y 4'
      particle(52)='NE NE y 4'
      particle(53)='M1 EL y 5'
      particle(54)='NM NE y 5'
      particle(55)='T1 EL y 6'
      particle(56)='NT NE y 6'
      particle(61)='E2 EL y 4'
      particle(63)='M2 EL y 5'
      particle(65)='T2 EL y 6'
      particle(71)='N1 C0 n  '
      particle(72)='N2 C0 n  '
      particle(73)='N3 C0 n  '
      particle(74)='N4 C0 n  '
      particle(70)='C1 C- y  '
      particle(76)='C2 C- y  '
      particle(81)='A0 A0 n  '
      particle(82)='H0 H0 n  '
      particle(83)='H- H- y  '   
      particle(84)='GL GL n  '
c auxiliary particles
      naux=100
      particle(101)=' 1  1 y  '
      particle(102)=' 2  2 y  '
      particle(103)=' 3  3 y  '
      particle(104)=' 4  4 y  '
      particle(105)=' 5  5 y  '
      particle(106)=' 6  6 y  '
      particle(107)=' 7  7 y  '
      particle(108)=' 8  8 y  '
      particle(109)=' 9  9 y  '
      particle(110)='10 10 y  '
      particle(111)='11 11 y  '
      particle(112)='12 12 y  '
      particle(113)='13 13 y  '
      particle(114)='14 14 y  '
      particle(115)='15 15 y  '
      particle(116)='16 16 y  '  ! non-standard couplings
c model settings
      if(smodel.eq.1)then        ! sm
        do i1=26,84
          particle(i1)=' '
        enddo
      elseif(smodel.eq.12)then   ! sm and qcd 
        do i1=31,84
          particle(i1)=' '
        enddo
      elseif(smodel.eq.2)then    ! qcd
        do i1=22,25  
          particle(i1)=' '
        enddo
        do i1=31,84
          particle(i1)=' '
        enddo
      elseif(smodel.eq.3)then    ! mssm
        particle(26)=' '
        particle(84)=' '
      elseif(smodel.eq.34)then   ! mssm and sqcd
      elseif(smodel.eq.4)then    ! sqcd
        do i1=22,25
          particle(i1)=' '
        enddo
        do i1=71,83
          particle(i1)=' '
        enddo
      endif
c additional particles for unstable particles
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
c setting particles
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
c output of couplings
      if(nout.ne.0.and.noutgen.eq.2)then
        write(nout,'(a)')' '
        write(nout,'(a)')' Couplings:'
        i4=1
        vertices=' '
        do i1=-naux+1,naux-1
        do i2=i1+1,naux
        do i3=i2+1,naux
        if(vertex(i1,i2,i3,schannel).and.
     *    (gname(i1)(3:3).ne.' '.or.i1.gt.0).and.
     *    (gname(i2)(3:3).ne.' '.or.i2.gt.0).and.
     *    (gname(i3)(3:3).ne.' '.or.i3.gt.0))then
            vertices=vertices(1:i4)//
     *        ' ('//pname(i1)//','//pname(i2)//','//pname(i3)//')'
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
        if(vertex(i1,i2,i3,schannel).and.
     *    (gname(i2)(3:3).ne.' '.or.i2.gt.0).and.
     *    (gname(i3)(3:3).ne.' '.or.i3.gt.0))then
          do i4=-naux+1,naux-1
          do i5=i4+1,naux-1
          if(vertex(-i1,i4,i5,schannel).and.
     *      (gname(i4)(3:3).ne.' '.or.i4.gt.0).and.
     *      (gname(i5)(3:3).ne.' '.or.i5.gt.0))then
            vertices=vertices(1:i6)//' ('//pname(i2)//','//
     *        pname(i3)//','//pname(i4)//','//pname(i5)//')'
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
c mapping
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
c particle identity
      do i1=1,nexternal(generator)
        idhep(i1,1)=hepnum(i1)
        if(i1.ge.3)idhep(i1,1)=-idhep(i1,1)
        massext2(i1,generator)=mass(abs(hepnum(i1)))**2
      enddo
c uniform sampling for checks (after massext2 definition!)
      if(smap.eq.0)then
        do i1=1,naux-1
          if(mass(i1).ne.0d0)mass(i1)=1d-5*i1
          width(i1)=0d0
          power(i1)=1d-5*i1
        enddo
      endif
c small fermion-higgs couplings
      do i1=0,maxv
        light(i1)=0
        if(mass(i1).eq.0d0)light(i1)=1
        light(-i1)=light(i1)
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize cuts                                             c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i1=0,allbinary(0,generator)
        mcutinv(i1,generator)=0d0
        if(i1.ge.1)tcutprocess(i1,generator)=0d0
      enddo
c angular cut
      if(sincludecuts.eq.1)then
        do i1=1,2
        do i2=3,nexternal(generator)
          if(ccut(i1,i2,generator).gt.0d0.and.
     *      ccut(i1,i2,generator).lt.1d0)then
            tcutprocess(2**(i1-1)+2**(i2-1),generator)=
     *        1d0/ccut(i1,i2,generator)**2-1d0
          endif
        enddo
        enddo
      endif
c invariant-mass cut
      do i1=1,allbinary(0,generator)
        scutinv=0d0
        do i2=1,nexternal(generator)
          if(included(i1,2**(i2-1),nexternal(generator)))then
             scutinv=scutinv+mass(abs(idhep(i2,1)))
           endif
        enddo
        scutinv=scutinv**2
        if(sincludecuts.eq.1)then
          do i2=1,nexternal(generator)
          do i3=i2+1,nexternal(generator)
            if(included(i1,2**(i2-1),nexternal(generator)).and.
     *        included(i1,2**(i3-1),nexternal(generator)))then
              m2=mass(abs(idhep(i2,1)))
              m3=mass(abs(idhep(i3,1)))
              e2=max(ecut(i2,generator),m2)
              e3=max(ecut(i3,generator),m3)
              scutinv=scutinv-(m2+m3)**2
     *          +max(scut(i2,i3,generator),m2**2+m3**2+2d0*e2*e3
     *              -2d0*dsqrt((e2**2-m2**2)*(e3**2-m3**2))
     *                *ccut(i2,i3,generator))
            endif
          enddo
          enddo
        endif
        mcutinv(i1,generator)=dsqrt(scutinv)
        do i2=2,i1-1
          if(included(i1,i2,nexternal(generator)))then
            mcutinv(i1,generator)=max(mcutinv(i1,generator),
     *        mcutinv(i2,generator)+mcutinv(i1-i2,generator))
          endif
        enddo
      enddo
      smin=mcutinv(allbinary(0,generator)-3,generator)**2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize channels                                         c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      channel=nchannel(generator)+1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     loop over emitter and spectator                             c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do em=0,nexternal(generator)
      do sp=0,nexternal(generator)
      do ga=0,nexternal(generator)
      n=0
      lmap(em,sp,ga)=.false.
      if(em.eq.0.and.sp.eq.0.and.ga.eq.0)then
        n=nexternal(generator)
      elseif(ssub.ne.0.and.em.ne.sp.and.em.ne.ga.and.sp.ne.ga.and.
     *  em.gt.0.and.sp.gt.0.and.ga.gt.0)then
        if(pname(hepnum(ga)).eq.'ph '.and.ga.ge.3.and.
     *    vertex(hepnum(em),-hepnum(em),hepnum(ga),schannel).and.
     *    vertex(hepnum(sp),-hepnum(sp),hepnum(ga),schannel))then
          n=nexternal(generator)-1
          idhep(ga,1)=0
        endif
      endif
      if(n.ne.0)then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     s-channel propagators                                       c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c number of s-channel propagators
      do maxns=0,n-3 
        ns=0
        if(maxns.eq.0)goto 500
c ns'th s-channel propagator
 100    ns=ns+1
        out1(ns)=2
c decay product 1 for ns'th decay
 200    out1(ns)=out1(ns)+1
        out2(ns)=out1(ns)
c decay product 2 for ns'th decay
 300    out2(ns)=out2(ns)+1
        virt(ns)=-maxv-1
c virtual particle for ns'th decay
 400    virt(ns)=virt(ns)+1        
        if(gname(virt(ns)).eq.'   '.and.virt(ns).lt.maxv)goto 400
        if(gname(virt(ns)).eq.'   '.and.virt(ns).eq.maxv)goto 1500
        if(virt(ns).lt.0.and.
     *    gname(virt(ns)).eq.gname(-virt(ns)))goto 400
c checking whether 3-particle vertex exists
        schannel=.true.
        if(.not.vertex(idhep(out1(ns),ns),idhep(out2(ns),ns),
     *    virt(ns),schannel))goto 1500
c checking last 3-particle vertex
        if(ns.eq.n-3.and..not.vertex(idhep(in1(ns),ns),
     *    idhep(in2(ns),ns),-virt(ns),schannel))goto 1500
c initializing next step
        do i2=1,n
          binary(i2,ns+1)=binary(i2,ns)
          idhep(i2,ns+1)=idhep(i2,ns)
        enddo
c combining particle out1 and out2 into new external particle out2 
        binary(out1(ns),ns+1)=0
        idhep(out1(ns),ns+1)=0
        binary(out2(ns),ns+1)=
     *    binary(out1(ns),ns)+binary(out2(ns),ns)
        idhep(out2(ns),ns+1)=-virt(ns)   
c find n s-channel propagator for ns < maxns
        if(ns.lt.maxns)goto 100
 500    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     t-channel propagators                                       c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c number of t-channel propagators
        maxnt=n-3
c loop of t-channel propagators
        if(maxns.eq.maxnt)goto 1000
        nt=maxns
c nt'th t-channel propagator
 600    nt=nt+1
        in1(nt)=0
c incoming particle 1 for nt's process
 700    in1(nt)=in1(nt)+1
c incoming particle 2 for nt's process
        in2(nt)=3-in1(nt)
        out1(nt)=2
c outgoing particle 1 for nt's process
 800    out1(nt)=out1(nt)+1
c virtual particle for nt's process 
        virt(nt)=-nmap                  ! (-nmap <-> -maxv-1)
 900    virt(nt)=virt(nt)+1
        if(gname(virt(nt)).eq.'   '.and.virt(nt).lt.maxv)goto 900
        if(gname(virt(nt)).eq.'   '.and.virt(nt).eq.maxv)goto 1300
        if(virt(nt).lt.0.and.
     *    gname(virt(nt)).eq.gname(-virt(nt)))goto 900
c avoid doube counting of diagrams
        if(nt.gt.maxns+1.and.in1(nt-1).eq.in1(nt))goto 1300
c checking whether 3-particle vertex exists
        schannel=.false.
        if(.not.vertex(idhep(in1(nt),nt),idhep(out1(nt),nt),
     *    virt(nt),schannel))goto 1300
c initializing n step
        do i2=1,n
          binary(i2,nt+1)=binary(i2,nt)
          idhep(i2,nt+1)=idhep(i2,nt)
        enddo
c combining particle in1 and out1 into new external particle in1 
        binary(out1(nt),nt+1)=0
        idhep(out1(nt),nt+1)=0
        binary(in1(nt),nt+1)=binary(in1(nt),nt)+binary(out1(nt),nt)
        idhep(in1(nt),nt+1)=-virt(nt)
c find n t-channel propagator for nt < maxnt
        if(nt.lt.maxnt)goto 600
 1000   continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     checking last vertex                                        c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        i1=n-2
        do i2=3,n
          if(idhep(i2,i1).ne.0)i3=i2
        enddo
        schannel=.true.
        if(gname(idhep(i3,i1)).eq.'   '.or..not.
     *    vertex(idhep(1,i1),idhep(2,i1),idhep(i3,i1),schannel))
     *    goto 1300
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initializing channels                                       c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ninv(channel,generator)=n-4
        nprocess(channel,generator)=maxnt-maxns
        ndecay(channel,generator)=maxns
c initializing inv
        idhepfirst(channel,generator)=0
        if(nprocess(channel,generator).eq.0)
     *    idhepfirst(channel,generator)=abs(idhep(i3,i1))
        do i1=1,maxns
          ininv(i1,channel,generator)=
     *      binary(out1(i1),i1)+binary(out2(i1),i1)
          idhepinv(i1,channel,generator)=abs(virt(i1))
          powerinv(i1,channel,generator)=power(abs(virt(i1)))
        enddo
        do i1=maxns+1,maxnt-1
          i2=maxns+maxnt-i1 
          ininv(i2,channel,generator)=allbinary(ga,generator) 
     *      -binary(in1(i1),i1)-binary(in2(i1),i1)-binary(out1(i1),i1)
          idhepinv(i2,channel,generator)=0
          powerinv(i2,channel,generator)=0d0
        enddo
c minimal invariant-mass cuts for subroutine inv
        do i1=1,ninv(channel,generator)
        do i2=1,i1-1
          binary1=ininv(i1,channel,generator)
          binary2=ininv(i2,channel,generator)
          lmin(i2,i1,channel,generator)=.false.
          if(included(binary1,binary2,n))then
            lmin(i2,i1,channel,generator)=.true.
            do i3=1,i1-1
              binary3=ininv(i3,channel,generator)
              if(included(binary3,binary2,n)
     *          .and.i2.ne.i3)then
                lmin(i2,i1,channel,generator)=.false.
              endif
            enddo  
          endif          
        enddo
        enddo
c maximal invariant-mass cuts for subroutine inv
        do i1=1,ninv(channel,generator)
        do i2=1,i1-1
          binary1=allbinary(ga,generator)-3-ininv(i1,channel,generator)
          binary2=ininv(i2,channel,generator)
          lmax(i2,i1,channel,generator)=.false.
          if(included(binary1,binary2,n))then
            lmax(i2,i1,channel,generator)=.true.
            do i3=1,i1-1
              binary3=ininv(i3,channel,generator)
              if(included(binary3,binary2,n)
     *          .and.i2.ne.i3)then
                lmax(i2,i1,channel,generator)=.false.
              endif
            enddo  
          endif
        enddo
        enddo
c initializing process
        do i1=maxns+1,maxnt
          i2=i1-maxns
          in1process(i2,channel,generator)=allbinary(ga,generator)
     *      -binary(in1(i1),i1)
          in2process(i2,channel,generator)=allbinary(ga,generator)
     *      -binary(in2(i1),i1)
          out1process(i2,channel,generator)=binary(out1(i1),i1)
          out2process(i2,channel,generator)=allbinary(ga,generator) 
     *      -binary(in1(i1),i1)-binary(in2(i1),i1)-binary(out1(i1),i1)
          inprocess(i2,channel,generator)=allbinary(ga,generator)
     *      -binary(in1(i1),i1)-binary(in2(i1),i1)
          virtprocess(i2,channel,generator)=allbinary(ga,generator)
     *      -binary(in1(i1),i1)-binary(out1(i1),i1)
          idhepprocess(i2,channel,generator)=abs(virt(i1))
          powerprocess(i2,channel,generator)=power(abs(virt(i1)))
        enddo
c initializing decay
        do i1=1,maxns
          i2=maxns-i1+1
          indecay(i2,channel,generator)=
     *      binary(out1(i1),i1)+binary(out2(i1),i1)
          out1decay(i2,channel,generator)=binary(out1(i1),i1)
          out2decay(i2,channel,generator)=binary(out2(i1),i1)
        enddo
c initializing map
         lmap(em,sp,ga)=.true.
         emmap(channel,generator)=em
         spmap(channel,generator)=sp
         gamap(channel,generator)=ga
c checking whether generator already exists
        do i1=1,channel-1
        if(comparechannel(i1,channel,smap,generator))then
          if(tcutprocess(allbinary(ga,generator)
     *      -virtprocess(1,channel,generator),generator).ne.0d0.or.
     *      tcutprocess(virtprocess(1,channel,generator),generator)
     *      .ne.0d0)then
            call copychannel(channel,i1,generator)
          endif
          goto 1300
        endif
        enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     output                                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(nout.ne.0.and.noutgen.eq.2)then
          write(nout,'(a12,i6)')' channel   =',channel
          if(em.ne.0.or.sp.ne.0.or.ga.ne.0)then
            write(nout,'(a12,i6)')' emitter   =',em
            write(nout,'(a12,i6)')' spectator =',sp
            write(nout,'(a12,i6)')' photon    =',ga
          endif
          do i1=1,maxns
            write(nout,'(a30,i6," +",i6," ->",i6)')
     *        '   s channel: ('//pname(-idhep(out1(i1),i1))//','//
     *        pname(-idhep(out2(i1),i1))//','//
     *        pname(-idhep(out2(i1),i1+1))//'), ',
     *        pid(binary(out1(i1),i1)),pid(binary(out2(i1),i1)),
     *        pid(binary(out2(i1),i1+1))
          enddo
          do i1=maxns+1,maxnt
            write(nout,'(a30,i6," +",i6," ->",i6)')
     *        '   t channel: ('//pname(idhep(in1(i1),i1))//','//
     *        pname(-idhep(out1(i1),i1))//','//
     *        pname(idhep(in1(i1),i1+1))//'), ',
     *        pid(binary(in1(i1),i1)),pid(binary(out1(i1),i1)),
     *        pid(binary(in1(i1),i1+1))
          enddo
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     end of loops                                                c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c end of t-channel propagators
 1500   if(maxns.gt.0)then
 1600     if(virt(ns).lt.maxv)goto 400
          if(out2(ns).lt.n)goto 300
          if(out1(ns).lt.n-1)goto 200
          if(ns.gt.1)then
            ns=ns-1
            goto 1600
          endif
        endif
c end of s-channel propagators
      enddo
c number of channels
      nchannel(generator)=channel-1
      if(nchannel(generator).gt.maxch)then
        write(*,'(a28,i6)')' initgenerator: reset maxch >',
     *    nchannel(generator)
        stop
      endif
      if(nchannel(generator).eq.0)then
        write(*,'(a47,i2)')
     *    ' initgenerator: no channel found for generator ',generator
        stop 
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize densities                                        c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initializing inv
      maxinv(generator)=0
      do channel=1,nchannel(generator)
      do ns=1,ninv(channel,generator)
        numinv(ns,channel,generator)=0
        do i1=1,channel-1
        do i2=1,ninv(i1,generator)
          if(compareinv(ns,channel,i2,i1,generator))
     *      numinv(ns,channel,generator)=numinv(i2,i1,generator)
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
c process
      maxprocess(generator)=0
      do channel=1,nchannel(generator)
      do nt=1,nprocess(channel,generator)
        numprocess(nt,channel,generator)=0
        do i1=1,channel-1
        do i2=1,nprocess(i1,generator)
          if(compareprocess(nt,channel,i2,i1,generator))
     *      numprocess(nt,channel,generator)=numprocess(i2,i1,generator)
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
c decay
      maxdecay(generator)=0
      do channel=1,nchannel(generator)
      do ns=1,ndecay(channel,generator)
        numdecay(ns,channel,generator)=0
        do i1=1,channel-1
        do i2=1,ndecay(i1,generator)
          if(comparedecay(ns,channel,i2,i1,generator))
     *      numdecay(ns,channel,generator)=numdecay(i2,i1,generator)
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
      enddo
      enddo
      if(nout.ne.0.and.noutgen.ge.1)then
c output
        write(nout,'(a)')' '
        write(nout,'(a23,i2)')' Phase-space generator ',generator
        write(nout,'(" Number of channels            =",i6)')
     *    nchannel(generator)  
        write(nout,'(" Calculation of invariants     =",i6)')
     *    maxinv(generator)
        write(nout,'(" Calculation of 2->2 processes =",i6)')
     *    maxprocess(generator)
        write(nout,'(" Calculation of 1->2 decays    =",i6)')
     *    maxdecay(generator)
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize adaptive optimization                            c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     check whether generator already exists                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function comparechannel(ch1,ch2,smap,generator)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=7,maxch=10000,maxg=4,maxv=200)
      integer i1,i2,ch1,ch2,generator,prop1,prop2,smap,idhep1,idhep2
      integer binary1,binary2,ga
      logical comparechannel,same,exist
c mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
c mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
c mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
c mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv,
     *  lmin,lmax
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      comparechannel=.false.
c checking emitter, spectator, and photon
      if(emmap(ch1,generator).ne.emmap(ch2,generator))return
      if(spmap(ch1,generator).ne.spmap(ch2,generator))return
      if(gamap(ch1,generator).ne.gamap(ch2,generator))return
      ga=gamap(ch1,generator)
c checking number of s-channel propagators
      prop1=ndecay(ch1,generator)
      do i1=1,ndecay(ch1,generator)
        idhep1=idhepinv(i1,ch1,generator)
        if(mass(idhep1).eq.0d0.and.width(idhep1).eq.0d0.and.
     *    power(idhep1).eq.0d0)prop1=prop1-1
      enddo
      prop2=ndecay(ch2,generator)
      do i1=1,ndecay(ch2,generator)
        idhep2=idhepinv(i1,ch2,generator)
        if(mass(idhep2).eq.0d0.and.width(idhep2).eq.0d0.and.
     *    power(idhep2).eq.0d0)prop2=prop2-1
      enddo
c checking s-channel propagators
      if(prop1.eq.prop2)then
        do i1=1,ndecay(ch1,generator)
          binary1=ininv(i1,ch1,generator)
          idhep1=idhepinv(i1,ch1,generator)
          if(mass(idhep1).ne.0d0.or.width(idhep1).ne.0d0.or.
     *      power(idhep1).ne.0d0)then
            exist=.false.
            do i2=1,ndecay(ch2,generator)
              binary2=ininv(i2,ch2,generator)
              idhep2=idhepinv(i2,ch2,generator)
              if(binary1.eq.binary2.and.
     *          mass(idhep1).eq.mass(idhep2).and.
     *          width(idhep1).eq.width(idhep2).and.
     *          power(idhep1).eq.power(idhep2))exist=.true.
            enddo
            if(.not.exist)return
          endif
        enddo 
      else 
        return
      endif
c checking number of t-channel propagators 
      prop1=nprocess(ch1,generator)
      do i1=1,nprocess(ch1,generator)
        idhep1=idhepprocess(i1,ch1,generator)
        if(mass(idhep1).eq.0d0.and.width(idhep1).eq.0d0.and.
     *    power(idhep1).eq.0d0)prop1=prop1-1
      enddo
      prop2=nprocess(ch2,generator)
      do i1=1,nprocess(ch2,generator)
        idhep2=idhepprocess(i1,ch2,generator)
        if(mass(idhep2).eq.0d0.and.width(idhep2).eq.0d0.and.
     *    power(idhep2).eq.0d0)prop2=prop2-1
      enddo
c checking t-channel propagators
      if(prop1.eq.prop2)then
        do i1=1,nprocess(ch1,generator)
          binary1=min(virtprocess(i1,ch1,generator),
     *      allbinary(ga,generator)-virtprocess(i1,ch1,generator))
          idhep1=idhepprocess(i1,ch1,generator)
          if(mass(idhep1).ne.0d0.or.width(idhep1).ne.0d0.or.
     *      power(idhep1).ne.0d0)then
            exist=.false.
            do i2=1,nprocess(ch2,generator)
              binary2=min(virtprocess(i2,ch2,generator),
     *          allbinary(ga,generator)-virtprocess(i2,ch2,generator))
              idhep2=idhepprocess(i2,ch2,generator)
              if(binary1.eq.binary2.and.
     *          mass(idhep1).eq.mass(idhep2).and.
     *          width(idhep1).eq.width(idhep2).and.
     *          power(idhep1).eq.power(idhep2))then
                  exist=.true.
              endif
            enddo 
            if(.not.exist)return
          endif
        enddo
      else 
        return
      endif  
c checking first s-channel propagator
      if(smap.eq.0.and.nprocess(ch1,generator).eq.0.and.
     *  nprocess(ch2,generator).eq.0)then
        idhep1=idhepfirst(ch1,generator)
        idhep2=idhepfirst(ch2,generator)
        if(mass(idhep1).ne.mass(idhep2).or.
     *    width(idhep1).ne.width(idhep2).or.
     *    power(idhep1).ne.power(idhep2))return
      endif
      comparechannel=.true.
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     copy channel ch1 to channel ch2                             c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copychannel(ch1,ch2,generator)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=7,maxch=10000,maxg=4,maxv=200)
      integer i1,i2,ch1,ch2,generator
c mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
c mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
c mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mckinematics/massext2,allbinary,nexternal
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcprocess/powerprocess,tcutprocess,in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv,
     *  lmin,lmax
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      do i1=1,maxe
c inv
        powerinv(i1,ch2,generator)=powerinv(i1,ch1,generator)
        ininv(i1,ch2,generator)=ininv(i1,ch1,generator)
        idhepinv(i1,ch2,generator)=idhepinv(i1,ch1,generator)
        idhepfirst(ch2,generator)=idhepfirst(ch1,generator)
        ninv(ch2,generator)=ninv(ch1,generator)
        do i2=1,maxe
          lmin(i1,i2,ch2,generator)=lmin(i1,i2,ch1,generator)
          lmax(i1,i2,ch2,generator)=lmax(i1,i2,ch1,generator)
        enddo
c decay
        indecay(i1,ch2,generator)=indecay(i1,ch1,generator)
        out1decay(i1,ch2,generator)=out1decay(i1,ch1,generator)
        out2decay(i1,ch2,generator)=out2decay(i1,ch1,generator)
        ndecay(ch2,generator)=ndecay(ch1,generator)
c process
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
c mapping
      emmap(ch2,generator)=emmap(ch1,generator)
      spmap(ch2,generator)=spmap(ch1,generator)
      gamap(ch2,generator)=gamap(ch1,generator)
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     check whether two 1->2 decays are equal                     c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function comparedecay(ns1,ch1,ns2,ch2,generator)
      implicit none
c local variables
      integer maxe,maxch,maxg
      parameter(maxe=7,maxch=10000,maxg=4)
      integer ns1,ns2,ch1,ch2,generator
      logical comparedecay
c mcdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcdecay/indecay,out1decay,out2decay,ndecay
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      comparedecay=.false.
      if(emmap(ch1,generator).ne.emmap(ch2,generator))return
      if(spmap(ch1,generator).ne.spmap(ch2,generator))return
      if(gamap(ch1,generator).ne.gamap(ch2,generator))return
      if(indecay(ns1,ch1,generator).ne.indecay(ns2,ch2,generator))
     *  return
      if(out1decay(ns1,ch1,generator).ne.out1decay(ns2,ch2,generator)
     *  .and.out1decay(ns1,ch1,generator).ne.
     *  out2decay(ns2,ch2,generator))return
      comparedecay=.true.
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     check whether two 2->2 processes are equal                  c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function compareprocess(ns1,ch1,ns2,ch2,generator)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=7,maxch=10000,maxg=4,maxv=200)
      integer ns1,ns2,ch1,ch2,generator,idhep1,idhep2,ga
      logical compareprocess
c mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
c mckinemetics
      real*8 massext2(maxe,maxg)
      integer allbinary(0:maxe,maxg),nexternal(maxg)
c mcprocess
      real*8 powerprocess(maxe,maxch,maxg),tcutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mckinematics/massext2,allbinary,nexternal
      common/mcprocess/powerprocess,tcutprocess,in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      compareprocess=.false.
      if(emmap(ch1,generator).ne.emmap(ch2,generator))return
      if(spmap(ch1,generator).ne.spmap(ch2,generator))return
      if(gamap(ch1,generator).ne.gamap(ch2,generator))return
      ga=gamap(ch1,generator)
      if(inprocess(ns1,ch1,generator).ne.inprocess(ns2,ch2,generator))
     *  return
      if(virtprocess(ns1,ch1,generator).ne.
     *  virtprocess(ns2,ch2,generator).and.
     *  virtprocess(ns1,ch1,generator).ne.
     *  allbinary(ga,generator)-virtprocess(ns2,ch2,generator))return
      if(powerprocess(ns1,ch1,generator).ne.
     *  powerprocess(ns2,ch2,generator))return
      idhep1=idhepprocess(ns1,ch1,generator)
      idhep2=idhepprocess(ns2,ch2,generator)
      if(mass(idhep1).ne.mass(idhep2).or.width(idhep1).ne.width(idhep2))
     *  return
      compareprocess=.true.
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     check whether calculation of invariants are equal           c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function compareinv(ns1,ch1,ns2,ch2,generator)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=7,maxch=10000,maxg=4,maxv=200)
      integer i1,i2,ns1,ns2,ch1,ch2,generator,idhep1,idhep2
      logical compareinv,included
c mcgenerator
      real*8 mass(0:maxv),width(0:maxv),power(0:maxv)
      integer nchannel(maxg)
c mcinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer idhepfirst(maxch,maxg),ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcgenerator/mass,width,power,nchannel
      common/mcinv/powerinv,mcutinv,ininv,idhepinv,idhepfirst,ninv,
     *  lmin,lmax
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      compareinv=.false.
      if(emmap(ch1,generator).ne.emmap(ch2,generator))return
      if(spmap(ch1,generator).ne.spmap(ch2,generator))return
      if(gamap(ch1,generator).ne.gamap(ch2,generator))return
      if(ininv(ns1,ch1,generator).ne.ininv(ns2,ch2,generator))return
      if(powerinv(ns1,ch1,generator).ne.powerinv(ns2,ch2,generator))
     *  return
      idhep1=idhepinv(ns1,ch1,generator)
      idhep2=idhepinv(ns2,ch2,generator)
      if(mass(idhep1).ne.mass(idhep2).or.width(idhep1).ne.width(idhep2))
     *  return
      do i1=1,ns1-1
        if(lmin(i1,ns1,ch1,generator))then
          included=.false.
          do i2=1,ns2-1
            if(lmin(i2,ns2,ch2,generator).and.
     *        ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator))
     *        included=.true.
          enddo
          if(.not.included)return
        endif
        if(lmax(i1,ns1,ch1,generator))then
          included=.false.
          do i2=1,ns2-1
            if(lmax(i2,ns2,ch2,generator).and.
     *        ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator))
     *        included=.true.
          enddo
          if(.not.included)return
        endif
      enddo    
      do i2=1,ns2-1
        if(lmin(i2,ns2,ch2,generator))then
          included=.false.
          do i1=1,ns1-1
            if(lmin(i1,ns1,ch1,generator).and.
     *        ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator))
     *        included=.true.
          enddo
          if(.not.included)return
        endif
        if(lmax(i2,ns2,ch2,generator))then
          included=.false.
          do i1=1,ns1-1
            if(lmax(i1,ns1,ch1,generator).and.
     *        ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator))
     *        included=.true.
          enddo
          if(.not.included)return
        endif
      enddo    
      compareinv=.true.
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     check if binary2 is included in binary1                     c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function included(binary1,binary2,nexternal)
      implicit none
c local variables
      integer binary1,binary2,i1,b1,b2,nexternal
      logical included
      included=.false.
      do i1=1,nexternal
        b1=binary1/2**(i1-1)
        b2=binary2/2**(i1-1)
        if(2*(b2/2).ne.b2.and.2*(b1/2).eq.b1)return
      enddo
      included=.true.
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     s-channel propagator                                       c
c                                                                c
c     written by Markus Roth                                     c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inv(random,s,g,mass,width,power,smin,smax,step,switch)
      implicit none
c local variable
      real*8 pi,random,s,g,mass,mass2,width,width2,power
      real*8 smax,smin,omax,omin,h,jacobian,denum
      integer step,switch
c mcoutput
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
          s=h(random,power,mass2,smin,smax,switch)
        elseif(step.eq.2)then
          denum=jacobian(power,mass2,s,smin,smax,switch)
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     2->2 particle process with t-channel propagator            c
c                                                                c
c     written by Markus Roth                                     c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine process(random,q1,q2,p1,p2,qt,g,mass,width,nu,
     *  s,s1,s2,t,t1,t2,tcut1,tcut2,step,switch)
      implicit none
c local variable
      real*8 pi,random(2),q1(0:3),q2(0:3),p1(0:3),p2(0:3),qt(0:3)
      real*8 g,mass,width,nu,s1,s2,t1,t2,tmin,tmax,tcut1,tcut2
      real*8 lambdas,lambdat,mass2,width2,t,phi,cost,s,roots
      real*8 gamma2,beta,gs,d,lambda,be,omega,omin,omax
      real*8 k1(0:3),q(0:3),h,jacobian,denum,cmax,cmax1,cmax2
      integer i1,step,switch
c mcoutput
      integer nout,numout,maxout
c mctechparam
      real*8 techparam(3)
      common/mcoutput/nout,numout,maxout
      common/mctechparam/techparam
      if(switch.eq.0)return
      pi=4d0*datan(1d0)
      if(tcut1.ne.0d0.or.tcut2.ne.0d0)then
        beta=(q1(3)+q2(3))/(q1(0)+q2(0))
        gamma2=1d0/(1d0-beta*beta)
        lambda=(s-s1-s2)**2-4d0*s1*s2
        if(dabs(q1(1)+q2(1)).gt.techparam(3).and.
     *    dabs(q1(2)+q2(2)).gt.techparam(3))then
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
          cost=((s+s1-s2)*(s+t1-t2)-2d0*s*(t1+s1)+2d0*mass2*s
     *          -2d0*s*width*mass*dtan(omega))/lambdas/lambdat
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
          cost=((s+s1-s2)*(s+t1-t2)-2d0*s*(t1+s1)
     *          -2d0*s*h(random(2),nu,-mass2,-tmax,-tmin,switch)
     *         )/lambdas/lambdat
        elseif(step.eq.2)then
          denum=pi*jacobian(nu,-mass2,-t,-tmax,-tmin,switch)
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
        call rotation(p1,phi,cost,switch)
        do i1=0,3
          q(i1)=q1(i1)+q2(i1) 
          k1(i1)=q1(i1)
        enddo
        call boost(q,k1,1d0,switch)
        if(k1(1).eq.0d0)then
          phi=dsign(pi*0.5d0,k1(2))
        else
          phi=datan(k1(2)/k1(1))
        endif
        if(k1(1).lt.0d0)phi=phi+pi
        cost=k1(3)/dsqrt(k1(1)**2+k1(2)**2+k1(3)**2)
        call rotation(p1,-phi,cost,switch)
        call boost(q,p1,-1d0,switch) 
        do i1=0,3
          p2(i1)=q(i1)-p1(i1)
          qt(i1)=q1(i1)-p1(i1)
        enddo
        t=qt(0)**2-qt(1)**2-qt(2)**2-qt(3)**2
      endif
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     1->2 particle decay without propagator                     c
c                                                                c
c     written by Markus Roth                                     c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine decay(random,q,p1,p2,g,s,s1,s2,step,switch)
      implicit none
c local variable
      real*8 pi,random(2),q(0:3),p1(0:3),p2(0:3)
      real*8 phi,cost,s,roots,lambda,denum,g,s1,s2
      integer i1,step,switch
c mcoutput
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
        call rotation(p1,phi,cost,switch)
        call boost(q,p1,-1d0,switch) 
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     h function for importance sampling                         c
c                                                                c
c     written by Markus Roth                                     c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function h(random,nu,mass2,xmin,xmax,switch)
      implicit none
      real*8 h,random,nu,xmin,xmax,mass2,m2
      integer switch
c mctechparam
      real*8 techparam(3)
c mcoutput
      integer nout,numout,maxout
      common/mctechparam/techparam
      common/mcoutput/nout,numout,maxout
      h=0d0
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
        h=random*xmax+(1d0-random)*xmin 
      elseif(nu.eq.1d0)then
        h=dexp(random*dlog(xmax-m2)+(1d0-random)*dlog(xmin-m2))+m2
      else
        h=(random*(xmax-m2)**(1d0-nu)
     *    +(1d0-random)*(xmin-m2)**(1d0-nu))**(1d0/(1d0-nu))+m2
      endif
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     jacobian of the h function                                 c
c                                                                c
c     written by Markus Roth                                     c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function jacobian(nu,mass2,x,xmin,xmax,switch)
      implicit none
      real*8 jacobian,nu,x,xmin,xmax,mass2,m2
      integer switch
c mctechparam
      real*8 techparam(3)
c mcoutput
      integer nout,numout,maxout
      common/mctechparam/techparam
      common/mcoutput/nout,numout,maxout
      jacobian=0d0
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
        jacobian=xmax-xmin
      elseif(nu.eq.1d0)then
        jacobian=(dlog(xmax-m2)-dlog(xmin-m2))*(x-m2)
      else
        jacobian=(((xmax-m2)**(1d0-nu)-(xmin-m2)**(1d0-nu))
     *    /(1d0-nu))*(x-m2)**nu
      endif
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     rotation                                                   c
c                                                                c
c     written by Markus Roth                                     c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rotation(p,phi,cost,switch)
      implicit none
c local variable
      real*8 p(0:3),phi,cost,sint,cosp,sinp,px,py,pz
      integer switch
c mcoutput
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     boost                                                      c
c                                                                c
c     written by Markus Roth                                     c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boost(q,p,dir,switch)
      implicit none
c local variable
      real*8 p(0:3),q(0:3),dir,m,m2,bx,by,bz,gamma
      real*8 p0,px,py,pz,bp,a
      integer switch
c mcoutput
      integer nout,numout,maxout
      common/mcoutput/nout,numout,maxout
      if(switch.eq.0)return
      m2=q(0)**2-q(1)**2-q(2)**2-q(3)**2
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     mapping from (n+1)- into n-particle phasespace              c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mapin(kt,g,em,sp,ga,next,switch)
      implicit none
c local variables
      integer maxe,maxg,maxch
      parameter(maxe=7,maxg=4,maxch=10000)
      real*8 k(maxe,0:3),kt(maxe,0:3),g,h,denum,jacobian
      real*8 x,z,y,v,kk(7,7),pi,ck(0:3),ckt(0:3),st
      integer i1,i2,j,em,sp,ga,next,switch
c mcoutput
      integer nout,numout,maxout
c mcmap
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
      kk(em,sp)=k(em,0)*k(sp,0)-k(em,1)*k(sp,1)-k(em,2)*k(sp,2)
     *  -k(em,3)*k(sp,3)
      kk(em,ga)=k(em,0)*k(ga,0)-k(em,1)*k(ga,1)-k(em,2)*k(ga,2)
     *  -k(em,3)*k(ga,3)
      kk(sp,ga)=k(sp,0)*k(ga,0)-k(sp,1)*k(ga,1)-k(sp,2)*k(ga,2)
     *  -k(sp,3)*k(ga,3)
c final-state emitter with final-state spectator
      if(em.ge.3.and.sp.ge.3)then
        x=1d0
        z=kk(em,sp)/(kk(em,sp)+kk(sp,ga))
        y=kk(em,ga)/(kk(em,sp)+kk(em,ga)+kk(sp,ga))
        do i1=0,3
          kt(em,i1)=k(em,i1)+k(ga,i1)-y/(1d0-y)*k(sp,i1)
          kt(sp,i1)=1d0/(1d0-y)*k(sp,i1)
          kt(ga,i1)=0d0
        enddo
        denum=pi*(1d0-y)*jacobian(powermap,0d0,1d0-y,0d0,1d0,switch)
     *    *jacobian(powermap,0d0,1d0-z,0d0,1d0,switch)
c final-state emitter with initial-state spectator
      elseif(em.ge.3.and.sp.le.2)then      
        x=(kk(em,sp)+kk(sp,ga)-kk(em,ga))/(kk(em,sp)+kk(sp,ga))
        z=kk(em,sp)/(kk(em,sp)+kk(sp,ga))
        do i1=0,3
          kt(em,i1)=k(em,i1)+k(ga,i1)-(1d0-x)*k(sp,i1)
          kt(sp,i1)=x*k(sp,i1)
          kt(ga,i1)=0d0
        enddo
        denum=pi*jacobian(powermap,0d0,1d0-x,0d0,1d0,switch)
     *    *jacobian(powermap,0d0,1d0-z,0d0,1d0,switch)            
c initial-state emitter with final-state spectator
      elseif(em.le.2.and.sp.ge.3)then      
        x=(kk(em,sp)+kk(em,ga)-kk(sp,ga))/(kk(em,sp)+kk(em,ga))
        z=kk(em,sp)/(kk(em,sp)+kk(em,ga))
        do i1=0,3
          kt(sp,i1)=k(sp,i1)+k(ga,i1)-(1d0-x)*k(em,i1)
          kt(em,i1)=x*k(em,i1)
          kt(ga,i1)=0d0
        enddo
        denum=pi*jacobian(powermap,0d0,1d0-x,0d0,1d0,switch)
     *    *jacobian(powermap,0d0,1d0-z,0d0,1d0,switch)            
c initial-state emitter with initial-state spectator
      elseif(em.le.2.and.sp.le.2)then      
        x=(kk(em,sp)-kk(em,ga)-kk(sp,ga))/kk(em,sp)
        v=kk(em,ga)/kk(em,sp)
        do i1=0,3
          ck(i1)=k(em,i1)+k(sp,i1)-k(ga,i1)
          ckt(i1)=x*k(em,i1)+k(sp,i1)
        enddo
        do j=3,6
          kk(em,j)=k(em,0)*k(j,0)-k(em,1)*k(j,1)-k(em,2)*k(j,2)
     *      -k(em,3)*k(j,3)
          kk(sp,j)=k(sp,0)*k(j,0)-k(sp,1)*k(j,1)-k(sp,2)*k(j,2)
     *      -k(sp,3)*k(j,3)
          kk(ga,j)=k(ga,0)*k(j,0)-k(ga,1)*k(j,1)-k(ga,2)*k(j,2)
     *      -k(ga,3)*k(j,3)
          do i1=0,3
            kt(j,i1)=k(j,i1)
     *        -((1d0+x)*kk(em,j)+2d0*kk(sp,j)-kk(ga,j))
     *          /(4d0*x+v-x*v)/kk(em,sp)*(ck(i1)+ckt(i1))
     *        +(kk(em,j)+kk(sp,j)-kk(ga,j))/x/kk(em,sp)*ckt(i1) 
          enddo
        enddo
        do i1=0,3
          kt(em,i1)=x*k(em,i1)
          kt(sp,i1)=k(sp,i1)
          kt(ga,i1)=0d0
        enddo 
        denum=pi*(1d0-x)*jacobian(powermap,0d0,1d0-x,0d0,1d0,switch)
      endif
      st=kt(em,0)*kt(sp,0)-kt(em,1)*kt(sp,1)-kt(em,2)*kt(sp,2)
     *  -kt(em,3)*kt(sp,3)
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     mapping from n- into (n+1)-particle phase space             c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mapout(random,k,x,g,em,sp,ga,next,switch)
      implicit none
c local variables
      integer maxch,maxg,maxe
      parameter(maxch=10000,maxg=4,maxe=7)
      real*8 random(3),k(maxe,0:3),kt(maxe,0:3),g,m2,m,pi
      real*8 kem(0:3),ksp(0:3),ktr(0:3),ck(0:3),ckt(0:3)
      real*8 ksum(0:3),kj(maxe,0:3),y,h,z,f,v,c,x,kk(maxe,maxe)
      integer i1,i2,j,em,sp,ga,next,switch
c mcoutput
      integer nout,numout,maxout
c mcmap
      real*8 powermap
      integer emmap(maxch,maxg),spmap(maxch,maxg),gamap(maxch,maxg)
      logical lmap(0:maxe,0:maxe,0:maxe)
      common/mcoutput/nout,numout,maxout
      common/mcmap/powermap,emmap,spmap,gamap,lmap
      if(switch.eq.0)return
      if(em.eq.0.and.sp.eq.0.and.ga.eq.0)return
      pi=4d0*datan(1d0)
c mapping from n- into (n+1)-particle phasespace
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
c final-state emitter with final-state spectator
        y=1d0-h(random(1),powermap,0d0,0d0,1d0,switch)
        z=1d0-h(random(2),powermap,0d0,0d0,1d0,switch)
        f=dsqrt(y*z/(1d0-z))
        ktr(0)=-f*f*m
        ktr(1)=f*m*dcos(2d0*pi*random(3))
        ktr(2)=f*m*dsin(2d0*pi*random(3))
        ktr(3)=-f*f*m
        call transverse(ksum,ksp,ktr,switch)
        do i1=0,3
          k(em,i1)=y*(1d0+z)*ksp(i1)+z*kem(i1)+(1d0-z)*ktr(i1)
          k(ga,i1)=(1d0-z)*kem(i1)-y*z*ksp(i1)-(1d0-z)*ktr(i1)
          k(sp,i1)=(1d0-y)*ksp(i1)
        enddo
c final-state emitter with initial-state spectator
      elseif(em.ge.3.and.sp.le.2)then      
        z=1d0-h(random(2),powermap,0d0,0d0,1d0,switch)
        f=dsqrt(z*(1d0-x)/(1d0-z)/x)
        ktr(0)=-f*f*m
        ktr(1)=f*m*dcos(2d0*pi*random(3))
        ktr(2)=f*m*dsin(2d0*pi*random(3))
        ktr(3)=-f*f*m
        call transverse(ksum,ksp,ktr,switch)
        do i1=0,3
          k(ga,i1)=(1-z)*kem(i1)-(1d0-x)/x*z*ksp(i1)-(1d0-z)*ktr(i1)
          k(em,i1)=(1d0-x)/x*(1d0+z)*ksp(i1)+z*kem(i1)
     *      +(1d0-z)*ktr(i1)
          k(sp,i1)=ksp(i1)/x
        enddo
c initial-state emitter with final-state spectator
      elseif(em.le.2.and.sp.ge.3)then
        z=1d0-h(random(2),powermap,0d0,0d0,1d0,switch)
        f=dsqrt((1d0-z)*(1d0-x)/z/x)
        ktr(0)=-f*f*m
        ktr(1)=f*m*dcos(2d0*pi*random(3))
        ktr(2)=f*m*dsin(2d0*pi*random(3))
        ktr(3)=-f*f*m
        call transverse(ksum,kem,ktr,switch)
        do i1=0,3
          k(ga,i1)=(1d0-x)/x*(2d0-z)*kem(i1)+(1d0-z)*ksp(i1)
     *      +z*ktr(i1)
          k(sp,i1)=z*ksp(i1)-(1d0-x)/x*(1d0-z)*kem(i1)-z*ktr(i1) 
          k(em,i1)=kem(i1)/x
        enddo       
c initial-state emitter with initial-state spectator
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
        call transverse(ksum,kem,ktr,switch)
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
          c=ktr(0)*kt(j,0)-ktr(1)*kt(j,1)-ktr(2)*kt(j,2)
     *      -ktr(3)*kt(j,3)
          kk(em,sp)=kem(0)*ksp(0)-kem(1)*ksp(1)-kem(2)*ksp(2)
     *      -kem(3)*ksp(3)
          kk(em,j)=kt(j,0)*kem(0)-kt(j,1)*kem(1)-kt(j,2)*kem(2)
     *      -kt(j,3)*kem(3)
          kk(sp,j)=kt(j,0)*ksp(0)-kt(j,1)*ksp(1)-kt(j,2)*ksp(2)
     *      -kt(j,3)*ksp(3)
          do i1=0,3
            k(j,i1)=kt(j,i1)
     *        -((2d0*x+v)*kk(em,j)+(2d0-v)*x*kk(sp,j)-x*c)
     *         /(4d0*x+v-x*v)/kk(em,sp)*(ck(i1)+ckt(i1))
     *        +(kk(em,j)+kk(sp,j))/kk(em,sp)*ck(i1) 
          enddo 
        endif 
        enddo
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     determination of transverse momenta                         c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine transverse(ksum,k,ktr,switch)
      implicit none
c local variables
      real*8 k(0:3),ktr(0:3),ksum(0:3),p(0:3),pvec,pi,phi,cost
      integer i1,switch
c output
      integer nout,numout,maxout
      common/output/nout,numout,maxout
      if(switch.eq.0)return
      pi=4d0*datan(1d0)
      do i1=0,3
        p(i1)=k(i1)
      enddo
      call boost(ksum,p,1d0,switch)
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
      call rotation(ktr,-phi,cost,switch)
      call boost(ksum,ktr,-1d0,switch)
      end       

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     3-particle vertex                                           c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function vertex(i1,i2,i3,schannel)
      implicit none
c local variables
      integer maxv
      parameter(maxv=200)
      integer i1,i2,i3
      logical vertex,vertexg,schannel
c mcparticle
      integer family(-maxv:maxv,6),light(-maxv:maxv)
      character*3 pname(-maxv:maxv),gname(-maxv:maxv)
      common/mcparticle/family,light,pname,gname
      vertex=.false.
      if(gname(i1).eq.'   ')return
      if(gname(i2).eq.'   ')return
      if(gname(i3).eq.'   ')return
c family conservation
      if(family(i1,1)+family(i2,1)+family(i3,1).ne.0)return
      if(family(i1,2)+family(i2,2)+family(i3,2).ne.0)return
      if(family(i1,3)+family(i2,3)+family(i3,3).ne.0)return
      if(family(i1,4)+family(i2,4)+family(i3,4).ne.0)return
      if(family(i1,5)+family(i2,5)+family(i3,5).ne.0)return
      if(family(i1,6)+family(i2,6)+family(i3,6).ne.0)return
      vertex=.true.
c testing vertex
      if(vertexg(i1,i2,i3,schannel))return
      if(vertexg(i1,i3,i2,schannel))return
      if(vertexg(i2,i1,i3,schannel))return
      if(vertexg(i2,i3,i1,schannel))return
      if(vertexg(i3,i1,i2,schannel))return
      if(vertexg(i3,i2,i1,schannel))return
      vertex=.false.
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     generic 3-particle vertex                                   c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function vertexg(i1,i2,i3,schannel)
      implicit none
c local variables
      integer maxv
      parameter(maxv=200)
      integer i1,i2,i3
      character*3 p1,p2,p3,g1,g2,g3
      logical vertexg,nonstandardcoup,schannel
c mcparticle
      integer family(-maxv:maxv,6),light(-maxv:maxv)
      character*3 pname(-maxv:maxv),gname(-maxv:maxv)
      common/mcparticle/family,light,pname,gname
      vertexg=.true.
      p1=pname(i1)
      p2=pname(i2)
      p3=pname(i3) 
      g1=gname(i1)
      g2=gname(i2)
      g3=gname(i3) 
c 2 charginos - higgs
      if(g1.eq.'C- '.and.g2.eq.'C-~'.and.g3.eq.'h0 ')return
      if(g1.eq.'C- '.and.g2.eq.'C-~'.and.g3.eq.'H0 ')return
      if(g1.eq.'C- '.and.g2.eq.'C-~'.and.g3.eq.'A0 ')return
c 2 leptons - higgs
      if(light(i2).eq.0)then
        if(g1.eq.'el '.and.g2.eq.'el~'.and.g3.eq.'h0 ')return
        if(g1.eq.'el '.and.g2.eq.'el~'.and.g3.eq.'H0 ')return
        if(g1.eq.'el '.and.g2.eq.'el~'.and.g3.eq.'A0 ')return
        if(g1.eq.'ne '.and.g2.eq.'el~'.and.g3.eq.'H- ')return
        if(g1.eq.'ne~'.and.g2.eq.'el '.and.g3.eq.'H-~')return
      endif
c 2 neutralinos - higgs
      if(g1.eq.'C0 '.and.g2.eq.'C0 '.and.g3.eq.'h0 ')return
      if(g1.eq.'C0 '.and.g2.eq.'C0 '.and.g3.eq.'H0 ')return
      if(g1.eq.'C0 '.and.g2.eq.'C0 '.and.g3.eq.'A0 ')return
c 2 quarks - higgs
      if(light(i1).eq.0.and.light(i2).eq.0)then
        if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'h0 ')return
        if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'h0 ')return
        if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'H0 ')return
        if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'H0 ')return
        if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'A0 ')return
        if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'A0 ')return
        if(g1.eq.'uq '.and.g2.eq.'dq~'.and.g3.eq.'H- ')return
        if(g1.eq.'uq~'.and.g2.eq.'dq '.and.g3.eq.'H-~')return
      endif
c chargino - lepton - slepton
      if(g1.eq.'ne '.and.g2.eq.'EL~'.and.g3.eq.'C- ')return
      if(g1.eq.'ne~'.and.g2.eq.'EL '.and.g3.eq.'C-~')return
      if(g1.eq.'NE '.and.g2.eq.'el~'.and.g3.eq.'C- ')return
      if(g1.eq.'NE~'.and.g2.eq.'el '.and.g3.eq.'C-~')return
c chargino - neutralino - higgs
      if(g1.eq.'C0 '.and.g2.eq.'C-~'.and.g3.eq.'H- ')return
      if(g1.eq.'C0 '.and.g2.eq.'C- '.and.g3.eq.'H-~')return
c chargino - quark - squark
      if(g1.eq.'uq '.and.g2.eq.'DQ~'.and.g3.eq.'C- ')return
      if(g1.eq.'uq~'.and.g2.eq.'DQ '.and.g3.eq.'C-~')return
      if(g1.eq.'UQ '.and.g2.eq.'dq~'.and.g3.eq.'C- ')return
      if(g1.eq.'UQ~'.and.g2.eq.'dq '.and.g3.eq.'C-~')return
c lepton - neutralino - slepton
      if(g1.eq.'ne '.and.g2.eq.'NE~'.and.g3.eq.'C0 ')return
      if(g1.eq.'ne~'.and.g2.eq.'NE '.and.g3.eq.'C0 ')return
      if(g1.eq.'el '.and.g2.eq.'EL~'.and.g3.eq.'C0 ')return
      if(g1.eq.'el~'.and.g2.eq.'EL '.and.g3.eq.'C0 ')return
c neutralino - quark - squark
      if(g1.eq.'uq '.and.g2.eq.'UQ~'.and.g3.eq.'C0 ')return
      if(g1.eq.'uq~'.and.g2.eq.'UQ '.and.g3.eq.'C0 ')return
      if(g1.eq.'dq '.and.g2.eq.'DQ~'.and.g3.eq.'C0 ')return
      if(g1.eq.'dq~'.and.g2.eq.'DQ '.and.g3.eq.'C0 ')return
c 2 charginos - gauge boson
      if(g1.eq.'C- '.and.g2.eq.'C-~'.and.g3.eq.'ph ')return
      if(g1.eq.'C- '.and.g2.eq.'C-~'.and.g3.eq.'Z0 ')return
c 2 leptons - gauge boson
      if(g1.eq.'el '.and.g2.eq.'el~'.and.g3.eq.'ph ')return
      if(g1.eq.'ne '.and.g2.eq.'ne~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'el '.and.g2.eq.'el~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'ne '.and.g2.eq.'el~'.and.g3.eq.'W- ')return
      if(g1.eq.'ne~'.and.g2.eq.'el '.and.g3.eq.'W-~')return
c 2 neutralinos - gauge bosons
      if(g1.eq.'C0 '.and.g2.eq.'C0 '.and.g3.eq.'Z0 ')return
c 2 quarks - gauge boson
      if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'ph ')return
      if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'ph ')return
      if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'uq '.and.g2.eq.'dq~'.and.g3.eq.'W- ')return
      if(g1.eq.'uq~'.and.g2.eq.'dq '.and.g3.eq.'W-~')return
c chargino - neutralino - gauge boson
      if(g1.eq.'C0 '.and.g2.eq.'C-~'.and.g3.eq.'W- ')return
      if(g1.eq.'C0 '.and.g2.eq.'C- '.and.g3.eq.'W-~')return
c 3 higgs
      if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.'h0 ')return
      if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.'H0 ')return
      if(g1.eq.'h0 '.and.g2.eq.'H0 '.and.g3.eq.'H0 ')return
      if(g1.eq.'H0 '.and.g2.eq.'H0 '.and.g3.eq.'H0 ')return
      if(g1.eq.'h0 '.and.g2.eq.'A0 '.and.g3.eq.'A0 ')return
      if(g1.eq.'H0 '.and.g2.eq.'A0 '.and.g3.eq.'A0 ')return
      if(g1.eq.'h0 '.and.g2.eq.'H-~'.and.g3.eq.'H- ')return
      if(g1.eq.'H0 '.and.g2.eq.'H-~'.and.g3.eq.'H- ')return
c higgs - 2 sleptons
      if(light(i2).eq.0)then
        if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'h0 ')return
        if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'H0 ')return
        if(g1.eq.'EL '.and.g2.eq.'EL~'.and.g3.eq.'h0 ')return
        if(g1.eq.'EL '.and.g2.eq.'EL~'.and.g3.eq.'H0 ')return
        if(g1.eq.'EL '.and.g2.eq.'EL~'.and.g3.eq.'A0 ')return
        if(g1.eq.'NE '.and.g2.eq.'EL~'.and.g3.eq.'H- ')return
        if(g1.eq.'NE~'.and.g2.eq.'EL '.and.g3.eq.'H-~')return
      endif
c higgs - 2 squarks
      if(light(i1).eq.0.and.light(i2).eq.0)then
        if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.'h0 ')return
        if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.'H0 ')return
        if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.'A0 ')return
        if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.'h0 ')return
        if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.'H0 ')return
        if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.'A0 ')return
        if(g1.eq.'UQ '.and.g2.eq.'DQ~'.and.g3.eq.'H- ')return
        if(g1.eq.'UQ~'.and.g2.eq.'DQ '.and.g3.eq.'H-~')return
      endif
c 2 higgs - gauge boson
      if(g1.eq.'h0 '.and.g2.eq.'A0 '.and.g3.eq.'Z0 ')return
      if(g1.eq.'H0 '.and.g2.eq.'A0 '.and.g3.eq.'Z0 ')return
      if(g1.eq.'H- '.and.g2.eq.'H-~'.and.g3.eq.'ph ')return
      if(g1.eq.'H- '.and.g2.eq.'H-~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'h0 '.and.g2.eq.'H-~'.and.g3.eq.'W- ')return
      if(g1.eq.'h0 '.and.g2.eq.'H- '.and.g3.eq.'W-~')return
      if(g1.eq.'H0 '.and.g2.eq.'H-~'.and.g3.eq.'W- ')return
      if(g1.eq.'H0 '.and.g2.eq.'H- '.and.g3.eq.'W-~')return
      if(g1.eq.'A0 '.and.g2.eq.'H-~'.and.g3.eq.'W- ')return
      if(g1.eq.'A0 '.and.g2.eq.'H- '.and.g3.eq.'W-~')return
c 2 sleptons - gauge boson
      if(g1.eq.'EL '.and.g2.eq.'EL~'.and.g3.eq.'ph ')return
      if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'EL '.and.g2.eq.'EL~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'NE '.and.g2.eq.'EL~'.and.g3.eq.'W- ')return
      if(g1.eq.'NE~'.and.g2.eq.'EL '.and.g3.eq.'W-~')return
c 2 squarks - gauge boson
      if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.'ph ')return
      if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.'ph ')return
      if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.'Z0 ')return
      if(g1.eq.'UQ '.and.g2.eq.'DQ~'.and.g3.eq.'W- ')return
      if(g1.eq.'UQ~'.and.g2.eq.'DQ '.and.g3.eq.'W-~')return
c higgs - 2 gauge bosons
      if(g1.eq.'h0 '.and.g2.eq.'Z0 '.and.g3.eq.'Z0 ')return
      if(g1.eq.'H0 '.and.g2.eq.'Z0 '.and.g3.eq.'Z0 ')return
      if(g1.eq.'h0 '.and.g2.eq.'W- '.and.g3.eq.'W-~')return
      if(g1.eq.'H0 '.and.g2.eq.'W- '.and.g3.eq.'W-~')return
c 3 gauge bosons
      if(g1.eq.'ph '.and.g2.eq.'W- '.and.g3.eq.'W-~')return
      if(g1.eq.'Z0 '.and.g2.eq.'W- '.and.g3.eq.'W-~')return
c gluino - quark - squark
      if(g1.eq.'uq '.and.g2.eq.'UQ~'.and.g3.eq.'GL ')return
      if(g1.eq.'uq~'.and.g2.eq.'UQ '.and.g3.eq.'GL ')return
      if(g1.eq.'dq '.and.g2.eq.'DQ~'.and.g3.eq.'GL ')return
      if(g1.eq.'dq~'.and.g2.eq.'DQ '.and.g3.eq.'GL ')return
c 2 gluinos - gluon
      if(g1.eq.'GL '.and.g2.eq.'GL '.and.g3.eq.'gl ')return
c 2 quarks - gluon
      if(g1.eq.'uq '.and.g2.eq.'uq~'.and.g3.eq.'gl ')return
      if(g1.eq.'dq '.and.g2.eq.'dq~'.and.g3.eq.'gl ')return
c 2 squarks - gluon
      if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.'gl ')return
      if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.'gl ')return
c 3 gluons
      if(g1.eq.'gl '.and.g2.eq.'gl '.and.g3.eq.'gl ')return
c 4 particle vertices
      if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.' 1 ')return
      if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.' 1 ')return
      if(g1.eq.'gl '.and.g2.eq.'gl '.and.g3.eq.' 1~')return
      if(g1.eq.'gl '.and.g2.eq.'ph '.and.g3.eq.' 1~')return
      if(g1.eq.'gl '.and.g2.eq.'Z0 '.and.g3.eq.' 1~')return
      if(g1.eq.'UQ '.and.g2.eq.'DQ~'.and.g3.eq.' 2 ')return
      if(g1.eq.'gl '.and.g2.eq.'W- '.and.g3.eq.' 2~')return
      if(g1.eq.'UQ~'.and.g2.eq.'DQ '.and.g3.eq.' 3~')return
      if(g1.eq.'gl '.and.g2.eq.'W-~'.and.g3.eq.' 3 ')return
      if(schannel)then
        if(g1.eq.'gl '.and.g2.eq.'gl '.and.g3.eq.' 4 ')return
        if(g1.eq.'gl '.and.g2.eq.'gl '.and.g3.eq.' 4~')return
      endif
      if(g1.eq.'NE '.and.g2.eq.'EL~'.and.g3.eq.' 5 ')return
      if(g1.eq.'UQ '.and.g2.eq.'DQ~'.and.g3.eq.' 5 ')return
      if(g1.eq.'ph '.and.g2.eq.'W-~'.and.g3.eq.' 5 ')return
      if(g1.eq.'Z0 '.and.g2.eq.'W-~'.and.g3.eq.' 5 ')return
      if(g1.eq.'ph '.and.g2.eq.'W- '.and.g3.eq.' 5~')return
      if(g1.eq.'Z0 '.and.g2.eq.'W- '.and.g3.eq.' 5~')return
      if(g1.eq.'NE~'.and.g2.eq.'EL '.and.g3.eq.' 6~')return
      if(g1.eq.'UQ~'.and.g2.eq.'DQ '.and.g3.eq.' 6~')return
      if(g1.eq.'ph '.and.g2.eq.'W- '.and.g3.eq.' 6~')return
      if(g1.eq.'Z0 '.and.g2.eq.'W- '.and.g3.eq.' 6~')return
      if(g1.eq.'ph~'.and.g2.eq.'W-~'.and.g3.eq.' 6 ')return
      if(g1.eq.'Z0 '.and.g2.eq.'W-~'.and.g3.eq.' 6 ')return
      if(g1.eq.'H+ '.and.g2.eq.'H- '.and.g3.eq.' 7 ')return
      if(g1.eq.'W+ '.and.g2.eq.'W- '.and.g3.eq.' 7 ')return
      if(g1.eq.'EL '.and.g2.eq.'EL~'.and.g3.eq.' 7 ')return
      if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.' 7 ')return
      if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.' 7 ')return
      if(g1.eq.'H+ '.and.g2.eq.'H- '.and.g3.eq.' 7~')return
      if(g1.eq.'W+ '.and.g2.eq.'W- '.and.g3.eq.' 7~')return
      if(g1.eq.'EL '.and.g2.eq.'EL~'.and.g3.eq.' 7~')return
      if(g1.eq.'UQ '.and.g2.eq.'UQ~'.and.g3.eq.' 7~')return
      if(g1.eq.'DQ '.and.g2.eq.'DQ~'.and.g3.eq.' 7~')return
      if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.' 7~')return
      if(g1.eq.'h0 '.and.g2.eq.'H0 '.and.g3.eq.' 7~')return
      if(g1.eq.'H0 '.and.g2.eq.'H0 '.and.g3.eq.' 7~')return
      if(g1.eq.'A0 '.and.g2.eq.'A0 '.and.g3.eq.' 7~')return
      if(g1.eq.'ph '.and.g2.eq.'ph '.and.g3.eq.' 7~')return
      if(g1.eq.'ph '.and.g2.eq.'Z0 '.and.g3.eq.' 7~')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.' 7~')return
      if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.' 7~')return
      if(schannel)then
        if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.' 8 ')return
        if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.' 8~')return
        if(g1.eq.'H0 '.and.g2.eq.'H0 '.and.g3.eq.' 9 ')return
        if(g1.eq.'H0 '.and.g2.eq.'H0 '.and.g3.eq.' 9~')return
        if(g1.eq.'A0 '.and.g2.eq.'A0 '.and.g3.eq.'10 ')return
        if(g1.eq.'A0 '.and.g2.eq.'A0 '.and.g3.eq.'10~')return
      endif
      if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'11 ')return
      if(g1.eq.'ph '.and.g2.eq.'ph '.and.g3.eq.'11 ')return
      if(g1.eq.'h0 '.and.g2.eq.'h0 '.and.g3.eq.'11 ')return
      if(g1.eq.'h0 '.and.g2.eq.'H0 '.and.g3.eq.'11~')return
      if(g1.eq.'H0 '.and.g2.eq.'H0 '.and.g3.eq.'11~')return
      if(g1.eq.'A0 '.and.g2.eq.'A0 '.and.g3.eq.'11~')return
      if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'11~')return
      if(g1.eq.'h0 '.and.g2.eq.'H0 '.and.g3.eq.'12 ')return
      if(g1.eq.'H0 '.and.g2.eq.'H0 '.and.g3.eq.'12~')return
      if(g1.eq.'A0 '.and.g2.eq.'A0 '.and.g3.eq.'12~')return
      if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'12~')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'12~')return
      if(g1.eq.'H0 '.and.g2.eq.'H0 '.and.g3.eq.'13 ')return
      if(g1.eq.'A0 '.and.g2.eq.'A0 '.and.g3.eq.'13~')return
      if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'13~')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'13~')return
      if(g1.eq.'A0 '.and.g2.eq.'A0 '.and.g3.eq.'14 ')return
      if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'14~')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'14~')return
      if(g1.eq.'NE '.and.g2.eq.'NE~'.and.g3.eq.'15 ')return
      if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'15~')return
c add non-standard couplings (with family conservation!)
c      if(nonstandardcoup(p1,p2,p3,g1,g2,g3,schannel))return
      vertexg=.false.
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     identifying particles (usually not used)                    c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function pid(binary)
      implicit none
c local variables
      integer pid,binary,i1,i2,i3,nexternal
      i3=1
      nexternal=30
      pid=0
      do i1=nexternal,1,-1
        i2=binary/2**(i1-1)
        if((i2/2)*2.ne.i2)then
          pid=pid+i1*i3
          i3=10*i3
        endif
      enddo
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     adaptive optimization                                       c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine optimization(weight,g,gsum,n,generator,nchannel)
      implicit none
c local variables
      integer maxch,maxo,maxg
      parameter(maxch=10000,maxo=20,maxg=4)
      real*8 weight,g(maxch),gsum
      integer n,nchannel,generator,i1
c mcadaptopt
      real*8 alphaopt(maxch,maxo,maxg),betaopt(0:maxch,maxg)
      real*8 wi(maxch,maxg),alphamin
      integer nopt(0:maxo,maxg),opt(maxg)
c mcoutput
      integer nout,numout,maxout
      common/mcoutput/nout,numout,maxout
      common/mcadaptopt/alphaopt,betaopt,wi,alphamin,nopt,opt
      if(nopt(opt(generator),generator).ne.0)then
        if(weight.ne.0d0.and.gsum.ne.0d0)then
          do i1=1,nchannel
            if(g(i1)*gsum.lt.0d0)then
              if(numout.lt.maxout)then
                write(nout,'(a)')' optimization: gi*gsum <= 0'
                numout=numout+1
              endif
              return
            else
              wi(i1,generator)=wi(i1,generator)
     *          +dsqrt(g(i1)*weight*weight/gsum)
            endif
          enddo
        endif
        if(n.ge.nopt(opt(generator),generator))then
          do i1=1,nchannel
            alphaopt(i1,opt(generator)+1,generator)=
     *        alphaopt(i1,opt(generator),generator)
     *          *dsqrt(wi(i1,generator))
            betaopt(i1,generator)=betaopt(i1-1,generator)
     *        +alphaopt(i1,opt(generator)+1,generator)
            wi(i1,generator)=0d0
          enddo
          opt(generator)=opt(generator)+1
          do i1=1,nchannel
            if(betaopt(nchannel,generator).ne.0d0)then
              alphaopt(i1,opt(generator),generator)=alphamin
     *            /dble(nchannel)
     *          +alphaopt(i1,opt(generator),generator)*(1d0-alphamin)
     *            /betaopt(nchannel,generator)
              betaopt(i1,generator)=alphamin*dble(i1)/dble(nchannel)
     *          +betaopt(i1,generator)*(1d0-alphamin)
     *            /betaopt(nchannel,generator)
            else
              alphaopt(i1,opt(generator),generator)=1d0/dble(nchannel)
              betaopt(i1,generator)=dble(i1)/dble(nchannel)
            endif
          enddo
        endif
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     pseudo-random number generator                              c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rans(random)
      implicit none 
c local variables
      real*8 s1,s2,s3,random
      integer init
      save init,s1,s2,s3
      data init/0/,s1/0d0/,s2/0d0/,s3/0d0/
      if (init.eq.0) then
        init=1
        s1=dsqrt(2d0)-1d0
        s2=dsqrt(3d0)-1d0
        s3=dsqrt(5d0)-2d0
      endif
      s1=dmod(s1+s2+s3,1d0)
      s2=dmod(s1+s2+s3,1d0)
      s3=dmod(s1+s2+s3,1d0)
      random=s1
      end  
