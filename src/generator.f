ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     General Monte Carlo generator                               c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phasespace(random,kbeam,k,x1,x2,g,channel,generator,
     *  switch)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      real*8 random(3*maxe-10),k(maxe,0:3),kbeam(2,0:3),g
      real*8 mmin,mmax,smin,smax,x1,x2
      integer i1,i2,i3,i4,ns,nt,channel,generator,switch
      integer ranstart,id,virtinv,inv1,inv2
c general
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
c cinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c cprocess
      real*8 powerprocess(maxe,maxch,maxg),ccutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c cdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
      common/general/alphaisr,scale,meisr,s,p,mass,width,nchannel,
     *  nexternal,allbinary
      common/cinv/powerinv,mcutinv,ininv,idhepinv,ninv,lmin,lmax
      common/cdecay/indecay,out1decay,out2decay,ndecay
      common/cprocess/powerprocess,ccutprocess, in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      if(switch.eq.0)return
      if(switch.eq.1)then
        do i1=1,allbinary(generator)
          p(0,i1)=0d0
        enddo
c incoming momenta
        do i1=0,3
          p(i1,1)=-x1*kbeam(1,i1)
          p(i1,2)=-x2*kbeam(2,i1)
          p(i1,3)=p(i1,1)+p(i1,2)
          p(i1,allbinary(generator)-1)=-p(i1,1)
          p(i1,allbinary(generator)-2)=-p(i1,2)
          p(i1,allbinary(generator)-3)=-p(i1,3)
        enddo
c square of center-of-mass energy
        s(3)=p(0,3)**2-p(1,3)**2-p(2,3)**2-p(3,3)**2
        s(allbinary(generator)-3)=s(3)
c event excluded by cuts after intial-state radiation
        if(s(3).le.mcutinv(allbinary(generator)-3,generator)**2)then
          switch=0
          return
        endif
      else
        g=1d0
      endif
c inv
      do ns=1,ninv(channel,generator)
        inv1=ininv(ns,channel,generator)
        inv2=allbinary(generator)-3-ininv(ns,channel,generator)
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
     *    switch)
      enddo
c process
      ranstart=ninv(channel,generator)
      do nt=1,nprocess(channel,generator)
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
     *    ccutprocess(allbinary(generator)               ! angular cuts
     *      -virtprocess(nt,channel,generator),generator),
     *    ccutprocess(virtprocess(nt,channel,generator),generator),  
     *    switch)
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
     *    switch)
      enddo
c calculating invariants
      if(switch.eq.1)then
c momenta
        do i1=0,3
          k(1,i1)=-p(i1,1)
          k(2,i1)=-p(i1,2)
          do i2=3,nexternal(generator)
            k(i2,i1)=p(i1,2**(i2-1))
          enddo
        enddo
c calculating invariants
        do i1=1,allbinary(generator)-1
          if(p(0,i1).eq.0d0)then
            if(p(0,allbinary(generator)-i1).eq.0d0)then
              do i2=1,nexternal(generator)
                i3=i1/2**(i2-1)
                if(2*(i3/2).ne.i3)then
                  do i4=0,3
                    p(i4,i1)=p(i4,i1-2**(i2-1))+p(i4,2**(i2-1))
                  enddo
                  goto 100
                endif
              enddo
 100          s(i1)=p(0,i1)**2-p(1,i1)**2-p(2,i1)**2-p(3,i1)**2
            else
              do i2=0,3
                p(i2,i1)=-p(i2,allbinary(generator)-i1)
              enddo
              s(i1)=s(allbinary(generator)-i1)
            endif
          endif
        enddo
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     calculating denities                                        c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine density(g,generator,switch)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      real*8 random(3*maxe-10),ginv(maxch),gprocess(maxch)
      real*8 gdecay(maxch),k(maxe,0:3),kbeam(2,0:3),x1,x2
      real*8 g(maxch),help,mmin,mmax,smin,smax,lambda,s1,s2
      real*8 t,t1,t2,d1,d2,d3,m(2**maxe)
      real*8 ccut1,ccut2,lambdas,lambdat,tmin,tmax
      real*8 mw,m2,p2,pi,omax,omin,jacobian,power
      real*8 q0,qvec,qvec1,qvec2,cmax1,cmax2,cmax
      integer i1,i2,i3,i4,i5,ns,nt,channel,generator,switch
      integer id,virtinv,inv1,inv2
c general
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
c cinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c cprocess
      real*8 powerprocess(maxe,maxch,maxg),ccutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c cdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
c cdensity
      integer nsinv(maxch,maxg),chinv(maxch,maxg),maxinv(maxg)
      integer ntprocess(maxch,maxg),chprocess(maxch,maxg)
      integer maxprocess(maxg),nsdecay(maxch,maxg),chdecay(maxch,maxg)
      integer maxdecay(maxg),numinv(maxe,maxch,maxg)
      integer numprocess(maxe,maxch,maxg),numdecay(maxe,maxch,maxg)
      common/general/alphaisr,scale,meisr,s,p,mass,width,nchannel,
     *  nexternal,allbinary
      common/cinv/powerinv,mcutinv,ininv,idhepinv,ninv,lmin,lmax
      common/cdecay/indecay,out1decay,out2decay,ndecay
      common/cprocess/powerprocess,ccutprocess, in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      common/cdensity/nsinv,chinv,maxinv,ntprocess,chprocess,
     *  maxprocess,nsdecay,chdecay,maxdecay,numinv,numprocess,numdecay
c inv
      do i1=1,maxinv(generator)
        ginv(i1)=1d0
        ns=nsinv(i1,generator)
        channel=chinv(i1,generator)
        inv1=ininv(ns,channel,generator)
        inv2=allbinary(generator)-3-ininv(ns,channel,generator)
        mmin=mcutinv(inv1,generator)
        mmax=mcutinv(inv2,generator)
        do i2=1,ns-1
          virtinv=ininv(i2,channel,generator)
          if(lmin(i2,ns,channel,generator).and.
     *      dsqrt(s(virtinv)).gt.mcutinv(virtinv,generator))then
            mmin=mmin+dsqrt(s(virtinv))-mcutinv(inv1,generator)
     *        +mcutinv(inv1-virtinv,generator)
            inv1=inv1-virtinv
          endif
          if(lmax(i2,ns,channel,generator).and.
     *      dsqrt(s(virtinv)).gt.mcutinv(virtinv,generator))then
            mmax=mmax+dsqrt(s(virtinv))-mcutinv(inv2,generator)
     *        +mcutinv(inv2-virtinv,generator)
            inv2=inv2-virtinv
          endif
        enddo        
        smin=mmin**2
        smax=(dsqrt(s(3))-mmax)**2
        call inv(
     *    random,                                        ! random number
     *    s(ininv(ns,channel,generator)),                ! invariant mass
     *    ginv(i1),                                      ! local density
     *    mass(idhepinv(ns,channel,generator)),          ! mass
     *    width(idhepinv(ns,channel,generator)),         ! width
     *    powerinv(ns,channel,generator),                ! mapping 
     *    smin,                                          ! lower bound 
     *    smax,                                          ! upper bound
     *    switch)
      enddo
c process
      do i1=1,maxprocess(generator)
        gprocess(i1)=1d0
        nt=ntprocess(i1,generator)
        channel=chprocess(i1,generator)
        call process(
     *    random,                                        ! random numbers
     *    p(0,in1process(nt,channel,generator)),         ! incomming particle 1
     *    p(0,in2process(nt,channel,generator)),         ! incomming particle 2
     *    p(0,out1process(nt,channel,generator)),        ! outgoing particle 1
     *    p(0,out2process(nt,channel,generator)),        ! outgoing particle 2
     *    p(0,virtprocess(nt,channel,generator)),        ! virtual particle
     *    gprocess(i1),                                  ! local density
     *    mass(idhepprocess(nt,channel,generator)),      ! mass
     *    width(idhepprocess(nt,channel,generator)),     ! width
     *    powerprocess(nt,channel,generator),            ! mapping 
     *    s(inprocess(nt,channel,generator)),            ! incomming particles 
     *    s(out1process(nt,channel,generator)),          ! outgoing particle 1
     *    s(out2process(nt,channel,generator)),          ! outgoing partcile 2
     *    s(virtprocess(nt,channel,generator)),          ! virtual particle
     *    s(in1process(nt,channel,generator)),           ! incomming particle 1
     *    s(in2process(nt,channel,generator)),           ! incomming particle 2
     *    ccutprocess(allbinary(generator)               ! angular cuts
     *      -virtprocess(nt,channel,generator),generator),
     *    ccutprocess(virtprocess(nt,channel,generator),generator),  
     *    switch)
      enddo
c decay
      do i1=1,maxdecay(generator)
        gdecay(i1)=1d0
        ns=nsdecay(i1,generator)
        channel=chdecay(i1,generator)
        call decay(
     *    random,                                        ! random numbers
     *    p(0,indecay(ns,channel,generator)),            ! incomming particle 
     *    p(0,out1decay(ns,channel,generator)),          ! outgoing particle 1
     *    p(0,out2decay(ns,channel,generator)),          ! outgoing particle 2
     *    gdecay(i1),                                    ! local density
     *    s(indecay(ns,channel,generator)),              ! incomming particle
     *    s(out1decay(ns,channel,generator)),            ! outgoing particle 1
     *    s(out2decay(ns,channel,generator)),            ! outgoing particle 2
     *    switch)
      enddo
c calculating densities        
      do channel=1,nchannel(generator)
        g(channel)=1d0
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
      subroutine initphasespace(name,generator,lightfermions,
     *  includecuts,sout)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv,maxo
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40,maxo=20)
      real*8 power(maxv),m2,m3,e2,e3,scutinv,pi,gamma
      integer idhep(maxe,maxe),binary(maxe,maxe),idhep2,idhep3
      integer i1,i2,i3,ns,nt,maxns,maxnt,binary1,binary2,binary3
      integer in1(maxe),in2(maxe),out1(maxe),out2(maxe),id,lightfermions
      integer virt(maxe),channel,generator,includecuts,prop2,prop3,sout
      character*3 gname(-maxv:maxv),name(maxe)
      logical vertex,included,exist,same,compareinv,comparedecay
      logical compareprocess
c general
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
c cinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
c cprocess
      real*8 powerprocess(maxe,maxch,maxg),ccutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
c cdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
c cdensity
      integer nsinv(maxch,maxg),chinv(maxch,maxg),maxinv(maxg)
      integer ntprocess(maxch,maxg),chprocess(maxch,maxg)
      integer maxprocess(maxg),nsdecay(maxch,maxg),chdecay(maxch,maxg)
      integer maxdecay(maxg),numinv(maxe,maxch,maxg)
      integer numprocess(maxe,maxch,maxg),numdecay(maxe,maxch,maxg)
c cisr
      real*8 gam,betae
c output
      integer nout,numout,maxout
c techparam
      real*8 a,techcut
c adaptopt
      real*8 alphaopt(maxch,maxo,maxg),betaopt(0:maxch,maxg)
      real*8 wi(maxch,maxg)
      integer nopt(0:maxo,maxg),opt(maxg)
c cuts
      real*8 ecut(maxe),scut(maxe,maxe),ccut(maxe,maxe)
      common/general/alphaisr,scale,meisr,s,p,mass,width,nchannel,
     *  nexternal,allbinary
      common/cinv/powerinv,mcutinv,ininv,idhepinv,ninv,lmin,lmax
      common/cdecay/indecay,out1decay,out2decay,ndecay
      common/cprocess/powerprocess,ccutprocess,in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      common/cdensity/nsinv,chinv,maxinv,ntprocess,chprocess,
     *  maxprocess,nsdecay,chdecay,maxdecay,numinv,numprocess,numdecay
      common/cisr/gam,betae
      common/output/nout,numout,maxout
      common/techparam/a,techcut
      common/adaptopt/alphaopt,betaopt,wi,nopt,opt
      common/cuts/ecut,scut,ccut
c technical parameter in h function and subroutine process
      a=1d-6
      techcut=1d-7
c generic particle names, see also subroutine vertexg
      do i1=-maxv,maxv
        gname(i1)='   '
      enddo
      gname(1)='dq '
      gname(2)='uq '
      gname(3)='sq '
      gname(4)='cq '
      gname(5)='bq '
      gname(6)='tq '
      gname(11)='el '
      gname(12)='ne '
      gname(13)='mu '
      gname(14)='nm '
      gname(15)='ta '
      gname(16)='nt '
      gname(22)='ga '
      gname(23)='Z0 '
      gname(24)='W+ '
      gname(25)='H0 '
c      gname(26)='gl '
      gname(30)='v1 ' ! auxiliary particle, 4-particle vertex
      gname(31)='v2 ' ! auxiliary particle, 4-particle vertex
      gname(32)='v3 ' ! auxiliary particle, 4-particle vertex
      do i1=1,maxv
        gname(-i1)=gname(i1)
        if(gname(i1).ne.'   '.and.gname(i1).ne.'ga '.and.
     *    gname(i1).ne.'Z0 '.and.gname(i1).ne.'H0 '.and.
     *    gname(i1).ne.'gl ')then
          gname(-i1)(3:3)='~'
        endif
      enddo
c particle identity
      do i1=1,nexternal(generator)
        idhep(i1,1)=0
        do i2=-maxv,maxv
          if(i1.le.2.and.name(i1).eq.gname(i2))idhep(i1,1)=i2
          if(i1.ge.3.and.name(i1).eq.gname(-i2))idhep(i1,1)=i2
        enddo
        if(idhep(i1,1).eq.0)then
          write(nout,'(a)')' Wrong particle name!'
          stop
        endif
      enddo
c initializing invariants
      do i1=1,nexternal(generator)
        s(2**(i1-1))=mass(abs(idhep(i1,1)))**2
      enddo
c mapping
      do i1=1,maxv
        power(i1)=0.9d0
        if(width(i1).gt.0d0)power(i1)=0d0
        if(gname(i1)(1:1).eq.'v')then
          mass(i1)=0d0
          width(i1)=0d0
          power(i1)=0d0
        endif
      enddo
c check with Madgraph
c      do i1=1,25
c        width(i1)=0d0
c        power(i1)=a+1d-5*i1
c      enddo
c initialize binary counting of particle number
      do i1=1,maxe
      do i2=1,maxe
        binary(i1,i2)=0
      enddo
      enddo
      allbinary(generator)=0
      do i1=1,nexternal(generator)
        binary(i1,1)=2**(i1-1)
        allbinary(generator)=allbinary(generator)+2**(i1-1)
      enddo
c initializing external masses
      do i1=1,nexternal(generator)
        s(allbinary(generator)-2**(i1-1))=s(2**(i1-1))
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize structure functions                              c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pi=4d0*datan(1d0)
      betae=2d0*alphaisr/pi*(2d0*dlog(scale/meisr)-1d0)
      gam=gamma(1d0+0.5d0*betae)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize cuts                                             c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i1=0,allbinary(generator)
        mcutinv(i1,generator)=0d0
        if(i1.ge.1)ccutprocess(i1,generator)=1d0
      enddo
c angular cut
      if(includecuts.eq.1)then
        do i1=1,2
        do i2=3,nexternal(generator)
          if(ccut(i1,i2).gt.0d0.and.ccut(i1,i2).lt.1d0)then
            ccutprocess(2**(i1-1)+2**(i2-1),generator)=ccut(i1,i2)
          endif
        enddo
        enddo
      endif
c invariant-mass cut
      do i1=1,allbinary(generator)
        scutinv=0d0
        do i2=1,nexternal(generator)
          if(included(i1,2**(i2-1),nexternal(generator)))then
             scutinv=scutinv+mass(abs(idhep(i2,1)))
           endif
        enddo
        scutinv=scutinv**2
        if(includecuts.eq.1)then
          do i2=1,nexternal(generator)
          do i3=i2+1,nexternal(generator)
            if(included(i1,2**(i2-1),nexternal(generator)).and.
     *        included(i1,2**(i3-1),nexternal(generator)))then
              m2=mass(abs(idhep(i2,1)))
              m3=mass(abs(idhep(i3,1)))
              e2=max(ecut(i2),m2)
              e3=max(ecut(i3),m3)
              scutinv=scutinv-(m2+m3)**2
     *          +max(scut(i2,i3),m2**2+m3**2+2d0*e2*e3
     *            -2d0*dsqrt((e2**2-m2**2)*(e3**2-m3**2))*ccut(i2,i3))
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize channels                                         c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      channel=nchannel(generator)+1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     s-channel propagators                                       c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c number of s-channel propagators
      do maxns=0,nexternal(generator)-3 
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
        if(.not.vertex(gname(idhep(out1(ns),ns)),
     *    gname(idhep(out2(ns),ns)),gname(virt(ns)),
     *    lightfermions))goto 1500
c checking last 3-particle vertex
        if(ns.eq.nexternal(generator)-3.and.
     *    .not.vertex(gname(idhep(in1(ns),ns)),
     *      gname(idhep(in2(ns),ns)),gname(-virt(ns)),
     *      lightfermions))goto 1500
c initializing next step
        do i2=1,nexternal(generator)
          binary(i2,ns+1)=binary(i2,ns)
          idhep(i2,ns+1)=idhep(i2,ns)
        enddo
c combining particle out1 and out2 into new external particle out2 
        binary(out1(ns),ns+1)=0
        idhep(out1(ns),ns+1)=0
        binary(out2(ns),ns+1)=
     *    binary(out1(ns),ns)+binary(out2(ns),ns)
        idhep(out2(ns),ns+1)=-virt(ns)   
c find next s-channel propagator for ns < maxns
        if(ns.lt.maxns)goto 100
 500    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     t-channel propagators                                       c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c number of t-channel propagators
        maxnt=nexternal(generator)-3
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
        virt(nt)=-maxv-1
 900    virt(nt)=virt(nt)+1
        if(gname(virt(nt)).eq.'   '.and.virt(nt).lt.maxv)goto 900
        if(gname(virt(nt)).eq.'   '.and.virt(nt).eq.maxv)goto 1300
        if(virt(nt).lt.0.and.
     *    gname(virt(nt)).eq.gname(-virt(nt)))goto 900
c avoid doube counting of diagrams
        if(nt.ge.2.and.in1(nt-1).eq.in1(nt))goto 1300
c checking whether 3-particle vertex exists
        if(.not.vertex(gname(idhep(in1(nt),nt)),
     *    gname(idhep(out1(nt),nt)),gname(virt(nt)),
     *    lightfermions))goto 1300
c initializing next step
        do i2=1,nexternal(generator)
          binary(i2,nt+1)=binary(i2,nt)
          idhep(i2,nt+1)=idhep(i2,nt)
        enddo
c combining particle in1 and out1 into new external particle in1 
        binary(out1(nt),nt+1)=0
        idhep(out1(nt),nt+1)=0
        binary(in1(nt),nt+1)=binary(in1(nt),nt)+binary(out1(nt),nt)
        idhep(in1(nt),nt+1)=-virt(nt)
c find next t-channel propagator for nt < maxnt
        if(nt.lt.maxnt)goto 600
 1000   continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     checking last vertex                                        c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        i1=nexternal(generator)-2
        do i2=3,nexternal(generator)
          if(idhep(i2,i1).ne.0)i3=i2
        enddo
        if(.not.vertex(gname(idhep(1,i1)),gname(idhep(2,i1)),
     *    gname(idhep(i3,i1)),lightfermions))goto 1300
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initializing channels                                       c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ninv(channel,generator)=nexternal(generator)-4
        nprocess(channel,generator)=maxnt-maxns
        ndecay(channel,generator)=maxns
c initializing inv
        do i1=1,maxns
          ininv(i1,channel,generator)=
     *      binary(out1(i1),i1)+binary(out2(i1),i1)
          idhepinv(i1,channel,generator)=abs(virt(i1))
          powerinv(i1,channel,generator)=power(abs(virt(i1)))
        enddo
        do i1=maxns+1,maxnt-1
          i2=maxns+maxnt-i1 
          ininv(i2,channel,generator)=allbinary(generator) 
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
          if(included(binary1,binary2,nexternal(generator)))then
            lmin(i2,i1,channel,generator)=.true.
            do i3=1,i1-1
              binary3=ininv(i3,channel,generator)
              if(included(binary3,binary2,nexternal(generator))
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
          binary1=allbinary(generator)-3-ininv(i1,channel,generator)
          binary2=ininv(i2,channel,generator)
          lmax(i2,i1,channel,generator)=.false.
          if(included(binary1,binary2,nexternal(generator)))then
            lmax(i2,i1,channel,generator)=.true.
            do i3=1,i1-1
              binary3=ininv(i3,channel,generator)
              if(included(binary3,binary2,nexternal(generator))
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
          in1process(i2,channel,generator)=allbinary(generator)
     *      -binary(in1(i1),i1)
          in2process(i2,channel,generator)=allbinary(generator)
     *      -binary(in2(i1),i1)
          out1process(i2,channel,generator)=binary(out1(i1),i1)
          out2process(i2,channel,generator)=allbinary(generator) 
     *      -binary(in1(i1),i1)-binary(in2(i1),i1)-binary(out1(i1),i1)
          inprocess(i2,channel,generator)=allbinary(generator)
     *      -binary(in1(i1),i1)-binary(in2(i1),i1)
          virtprocess(i2,channel,generator)=allbinary(generator)
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     check whether generator already exists                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i1=1,channel-1
          same=.true.
c checking number of s-channel propagators
          prop2=ndecay(channel,generator)
          do i2=1,ndecay(channel,generator)
            idhep2=idhepinv(i2,channel,generator)
            if(mass(idhep2).eq.0d0.and.width(idhep2).eq.0d0.and.
     *        power(idhep2).eq.0d0)prop2=prop2-1
          enddo
          prop3=ndecay(i1,generator)
          do i3=1,ndecay(i1,generator)
            idhep3=idhepinv(i3,i1,generator)
            if(mass(idhep3).eq.0d0.and.width(idhep3).eq.0d0.and.
     *        power(idhep3).eq.0d0)prop3=prop3-1
          enddo
c checking s-channel propagators
          if(prop2.eq.prop3)then
            do i2=1,ndecay(i1,generator)
              binary2=ininv(i2,i1,generator)
              idhep2=idhepinv(i2,i1,generator)
              if(mass(idhep2).ne.0d0.or.width(idhep2).ne.0d0.or.
     *          power(idhep2).ne.0d0)then
                exist=.false.
                do i3=1,ndecay(channel,generator)
                  binary3=ininv(i3,channel,generator)
                  idhep3=idhepinv(i3,channel,generator)
                  if(binary2.eq.binary3.and.
     *              mass(idhep2).eq.mass(idhep3).and.
     *              width(idhep2).eq.width(idhep3).and.
     *              power(idhep2).eq.power(idhep3))then
                      exist=.true.
                  endif
                enddo 
                if(.not.exist)same=.false.
              endif
            enddo 
          else 
            same=.false.
          endif
c checking number of t-channel propagators 
          prop2=nprocess(channel,generator)
          do i2=1,nprocess(channel,generator)
            idhep2=idhepprocess(i2,channel,generator)
            if(mass(idhep2).eq.0d0.and.width(idhep2).eq.0d0.and.
     *        power(idhep2).eq.0d0)prop2=prop2-1
          enddo
          prop3=nprocess(i1,generator)
          do i3=1,nprocess(i1,generator)
            idhep3=idhepprocess(i3,i1,generator)
            if(mass(idhep3).eq.0d0.and.width(idhep3).eq.0d0.and.
     *        power(idhep3).eq.0d0)prop3=prop3-1
          enddo
c checking t-channel propagators
          if(prop2.eq.prop3)then
            do i2=1,nprocess(i1,generator)
              binary2=min(virtprocess(i2,i1,generator),
     *          allbinary(generator)-virtprocess(i2,i1,generator))
              idhep2=idhepprocess(i2,i1,generator)
              if(mass(idhep2).ne.0d0.or.width(idhep2).ne.0d0.or.
     *          power(idhep2).ne.0d0)then
                exist=.false.
                do i3=1,nprocess(channel,generator)
                  binary3=min(virtprocess(i3,channel,generator),
     *              allbinary(generator)
     *              -virtprocess(i3,channel,generator))
                  idhep3=idhepprocess(i3,channel,generator)
                  if(binary2.eq.binary3.and.
     *              mass(idhep2).eq.mass(idhep3).and.
     *              width(idhep2).eq.width(idhep3).and.
     *              power(idhep2).eq.power(idhep3))then
                      exist=.true.
                  endif
                enddo 
                if(.not.exist)same=.false.
              endif
            enddo
          else 
            same=.false.
          endif  
          if(same)goto 1300
        enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     output                                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(sout.eq.2)then
          write(nout,*)'channel=',channel
          do i2=1,maxns
            write(nout,*)'s channel: ',
     *        id(binary(out1(i2),i2)),id(binary(out2(i2),i2)),' ->',
     *        id(binary(out2(i2),i2+1)),'   ',
     *        gname(-idhep(out1(i2),i2)),
     *        gname(-idhep(out2(i2),i2)),' -> ',
     *        gname(-idhep(out2(i2),i2+1))
          enddo
          do i2=maxns+1,maxnt
            write(nout,*)'t channel: ',
     *        id(binary(in1(i2),i2)),id(binary(out1(i2),i2)),' ->',
     *        id(binary(in1(i2),i2+1)),'   ',
     *        gname(idhep(in1(i2),i2)),
     *        gname(-idhep(out1(i2),i2)),' -> ',
     *        gname(idhep(in1(i2),i2+1))
          enddo
          write(nout,'(a)')'   '
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     end of loops                                                c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        channel=channel+1
 1300   if(maxns.lt.maxnt)then
 1400     if(virt(nt).lt.maxv)goto 900
          if(out1(nt).lt.nexternal(generator))goto 800
          if(in1(nt).lt.2)goto 700
          if(nt.gt.maxns+1)then
            nt=nt-1
            goto 1400
          endif
        endif
c end of t-channel propagators
 1500   if(maxns.gt.0)then
 1600     if(virt(ns).lt.maxv)goto 400
          if(out2(ns).lt.nexternal(generator))goto 300
          if(out1(ns).lt.nexternal(generator)-1)goto 200
          if(ns.gt.1)then
            ns=ns-1
            goto 1600
          endif
        endif
c end of s-channel propagators
      enddo
c number of channels
      nchannel(generator)=channel-1
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
      if(sout.ne.0)then
        write(nout,'(a)')' '
        write(nout,'(a)')' Phase-space generator:'
        write(nout,'(" number of channels            =",i6)')
     *    nchannel(generator)  
        write(nout,'(" calculation of invariants     =",i6)')
     *    maxinv(generator)
        write(nout,'(" calculation of 2->2 processes =",i6)')
     *    maxprocess(generator)
        write(nout,'(" calculation of 1->2 decays    =",i6)')
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
      nopt(1,generator)=100
      nopt(2,generator)=200
      nopt(3,generator)=300
      nopt(4,generator)=400
      nopt(5,generator)=500
      nopt(6,generator)=600
      nopt(7,generator)=700
      nopt(8,generator)=800
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
      integer maxe,maxch,maxg,maxv,maxo
      parameter(maxe=9,maxch=20000,maxg=1)
      integer ns1,ns2,ch1,ch2,generator
      logical comparedecay
c cdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
      common/cdecay/indecay,out1decay,out2decay,ndecay
      comparedecay=.false.
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
      integer maxe,maxch,maxg,maxv,maxo
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      integer ns1,ns2,ch1,ch2,generator,idhep1,idhep2
      logical compareprocess
c general
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
c cprocess
      real*8 powerprocess(maxe,maxch,maxg),ccutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
      common/general/alphaisr,scale,meisr,s,p,mass,width,nchannel,
     *  nexternal,allbinary
      common/cprocess/powerprocess,ccutprocess,in1process,
     *  in2process,out1process,out2process,inprocess,virtprocess,
     *  idhepprocess,nprocess
      compareprocess=.false.
      if(inprocess(ns1,ch1,generator).ne.inprocess(ns2,ch2,generator))
     *  return
      if(virtprocess(ns1,ch1,generator).ne.
     *  virtprocess(ns2,ch2,generator).and.
     *  virtprocess(ns1,ch1,generator).ne.
     *  allbinary(generator)-virtprocess(ns2,ch2,generator))return
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
c     check whether caculation of invariants are equal            c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function compareinv(ns1,ch1,ns2,ch2,generator)
      implicit none
c local variables
      integer maxe,maxch,maxg,maxv,maxo
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      integer i1,i2,ns1,ns2,ch1,ch2,generator,idhep1,idhep2
      logical compareinv,included
c general
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
c cinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
      common/general/alphaisr,scale,meisr,s,p,mass,width,nchannel,
     *  nexternal,allbinary
      common/cinv/powerinv,mcutinv,ininv,idhepinv,ninv,lmin,lmax
      compareinv=.false.
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
      subroutine inv(random,s,g,mass,width,power,smin,smax,switch)
      implicit none
c local variable
      real*8 pi,random,s,g,mass,mass2,width,width2,power
      real*8 smax,smin,omax,omin,h,jacobian,denum
      integer switch
c output
      integer nout,numout,maxout
      common/output/nout,numout,maxout
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
        if(switch.eq.1)then
          s=mass2+mass*width*dtan(random*(omax-omin)+omin)   
        elseif(switch.eq.2)then
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
        if(switch.eq.1)then
          s=h(random,power,mass2,smin,smax,switch)
        elseif(switch.eq.2)then
          denum=jacobian(power,mass2,s,smin,smax,switch)
          if(denum.gt.0d0)then
            g=g/denum
          else
            if(numout.lt.maxout)then
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
     *  s,s1,s2,t,t1,t2,ccut1,ccut2,switch)
      implicit none
c local variable
      real*8 pi,random(2),q1(0:3),q2(0:3),p1(0:3),p2(0:3),qt(0:3)
      real*8 g,mass,width,nu,s1,s2,t1,t2,tmin,tmax,ccut1,ccut2
      real*8 lambdas,lambdat,mass2,width2,t,phi,cost,s,roots
      real*8 qvec1,qvec2,q0,omega,omin,omax,k1(0:3),q(0:3),h
      real*8 jacobian,denum,cmax,cmax1,cmax2
      integer i1,switch
c output
      integer nout,numout,maxout
c techparam
      real*8 a,techcut
      common/output/nout,numout,maxout
      common/techparam/a,techcut
      pi=4d0*datan(1d0)
      if(switch.eq.0)return
      if((ccut1.ne.1d0.or.ccut2.ne.1d0).and.
     *  q1(1)+q2(1).eq.0d0.and.q1(2)+q2(2).eq.0d0)then
        q0=q1(0)+q2(0)
        if(s1.eq.0d0)then
          qvec1=(q1(3)+q2(3))*dsign(1d0,q1(3)) 
          cmax1=(q0*ccut1-qvec1)/(q0-qvec1*ccut1)
        else
          cmax1=1d0
        endif 
        if(s2.eq.0d0)then
          qvec2=(q1(3)+q2(3))*dsign(1d0,q2(3)) 
          cmax2=(q0*ccut2-qvec2)/(q0-qvec2*ccut2)
        else
          cmax2=1d0
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
      if(dabs(tmax).lt.techcut)tmax=0d0
      if(tmax.le.tmin)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' process: tmin >= tmax'
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(tmax.gt.0d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' process: tmax > 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(width.gt.0d0)then
        width2=width*width
        omin=datan((mass2-tmax)/mass/width)
        omax=datan((mass2-tmin)/mass/width)
        if(switch.eq.1)then
          phi=2d0*pi*random(1)
          omega=omin+(omax-omin)*random(2)
          cost=((s+s1-s2)*(s+t1-t2)-2d0*s*(t1+s1)+2d0*mass2*s
     *          -2d0*s*width*mass*dtan(omega))/lambdas/lambdat
        elseif(switch.eq.2)then
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
        if(switch.eq.1)then
          phi=2d0*pi*random(1)
          denum=lambdas*lambdat
          cost=((s+s1-s2)*(s+t1-t2)-2d0*s*(t1+s1)
     *          -2d0*s*h(random(2),nu,-mass2,-tmax,-tmin,switch)
     *         )/lambdas/lambdat
        elseif(switch.eq.2)then
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
      if(switch.eq.1)then
        roots=dsqrt(s)
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
        cost=k1(3)/dsqrt(k1(1)*k1(1)+k1(2)*k1(2)+k1(3)*k1(3))
        call rotation(p1,-phi,cost,switch)
        call boost(q,p1,-1d0,switch)        
        do i1=0,3
          p2(i1)=q(i1)-p1(i1)
          qt(i1)=q1(i1)-p1(i1)
        enddo
        t=qt(0)*qt(0)-qt(1)*qt(1)-qt(2)*qt(2)-qt(3)*qt(3)
      endif
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     1->2 particle decay without propagator                     c
c                                                                c
c     written by Markus Roth                                     c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine decay(random,q,p1,p2,g,s,s1,s2,switch)
      implicit none
c local variable
      real*8 pi,random(2),q(0:3),p1(0:3),p2(0:3)
      real*8 phi,cost,s,roots,lambda,denum,g,s1,s2
      integer i1,switch
c output
      integer nout,numout,maxout
      common/output/nout,numout,maxout
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
      if(switch.eq.1)then
        roots=dsqrt(s)
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
      elseif(switch.eq.2)then
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
c techparam
      real*8 a,techcut
c output
      integer nout,numout,maxout
      common/techparam/a,techcut
      common/output/nout,numout,maxout
      h=0d0
      if(switch.eq.0)return
      m2=mass2-a
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
c techparam
      real*8 a,techcut
c output
      integer nout,numout,maxout
      common/techparam/a,techcut
      common/output/nout,numout,maxout
      jacobian=0d0
      if(switch.eq.0)return
      m2=mass2-a
      if(xmax-m2.lt.0d0.or.xmin-m2.lt.0d0)then
        if(numout.lt.maxout)then
          write(nout,'(a)')' jacobian: xmin-m2, xmax-m2 < 0'
          numout=numout+1
        endif
        switch=0
        return
      endif
      if(nu.eq.0d0)then
        jacobian=xmax-xmin
      elseif(nu.eq.1d0)then
        jacobian=(dlog(xmax-m2)-dlog(xmin-m2))*(x-m2)
      else
        jacobian=(((xmax-m2)**(1d0-nu)-(xmin-m2)**(1d0-nu))
     *           /(1d0-nu))*(x-m2)**nu
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
c output
      integer nout,numout,maxout
      common/output/nout,numout,maxout
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
c output
      integer nout,numout,maxout
      common/output/nout,numout,maxout
      if(switch.eq.0)return
      m2=q(0)*q(0)-q(1)*q(1)-q(2)*q(2)-q(3)*q(3)
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
c     3-particle vertex                                           c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function vertex(p1,p2,p3,lightfermions)
      implicit none
c local variables
      character*3 p1,p2,p3
      integer lightfermions
      logical vertex,vertexg
      vertex=.true.
      if(vertexg(p1,p2,p3,lightfermions))return
      if(vertexg(p1,p3,p2,lightfermions))return
      if(vertexg(p2,p1,p3,lightfermions))return
      if(vertexg(p2,p3,p1,lightfermions))return
      if(vertexg(p3,p1,p2,lightfermions))return
      if(vertexg(p3,p2,p1,lightfermions))return
      vertex=.false.
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     generic 3-particle vertex                                   c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function vertexg(p1,p2,p3,lightfermions)
      implicit none
c local variables
      character*3 p1,p2,p3
      integer lightfermions
      logical vertexg
      vertexg=.true.
c ew three-particle vertices
      if(p1.eq.'el '.and.p2.eq.'ne~'.and.p3.eq.'W+ ')return
      if(p1.eq.'el~'.and.p2.eq.'ne '.and.p3.eq.'W+~')return
      if(p1.eq.'mu '.and.p2.eq.'nm~'.and.p3.eq.'W+ ')return
      if(p1.eq.'mu~'.and.p2.eq.'nm '.and.p3.eq.'W+~')return
      if(p1.eq.'ta '.and.p2.eq.'nt~'.and.p3.eq.'W+ ')return
      if(p1.eq.'ta~'.and.p2.eq.'nt '.and.p3.eq.'W+~')return
      if(p1.eq.'dq '.and.p2.eq.'uq~'.and.p3.eq.'W+ ')return
      if(p1.eq.'dq~'.and.p2.eq.'uq '.and.p3.eq.'W+~')return
      if(p1.eq.'sq '.and.p2.eq.'cq~'.and.p3.eq.'W+ ')return
      if(p1.eq.'sq~'.and.p2.eq.'cq '.and.p3.eq.'W+~')return
      if(p1.eq.'bq '.and.p2.eq.'tq~'.and.p3.eq.'W+ ')return
      if(p1.eq.'bq~'.and.p2.eq.'tq '.and.p3.eq.'W+~')return
      if(p1.eq.'ne '.and.p2.eq.'ne~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'el '.and.p2.eq.'el~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'el '.and.p2.eq.'el~'.and.p3.eq.'ga ')return
      if(p1.eq.'nm '.and.p2.eq.'nm~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'mu '.and.p2.eq.'mu~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'mu '.and.p2.eq.'mu~'.and.p3.eq.'ga ')return
      if(p1.eq.'nt '.and.p2.eq.'nt~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'ta '.and.p2.eq.'ta~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'ta '.and.p2.eq.'ta~'.and.p3.eq.'ga ')return
      if(p1.eq.'dq '.and.p2.eq.'dq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'dq '.and.p2.eq.'dq~'.and.p3.eq.'ga ')return
      if(p1.eq.'uq '.and.p2.eq.'uq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'uq '.and.p2.eq.'uq~'.and.p3.eq.'ga ')return
      if(p1.eq.'sq '.and.p2.eq.'sq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'sq '.and.p2.eq.'sq~'.and.p3.eq.'ga ')return
      if(p1.eq.'cq '.and.p2.eq.'cq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'cq '.and.p2.eq.'cq~'.and.p3.eq.'ga ')return
      if(p1.eq.'bq '.and.p2.eq.'bq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'bq '.and.p2.eq.'bq~'.and.p3.eq.'ga ')return
      if(p1.eq.'tq '.and.p2.eq.'tq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'tq '.and.p2.eq.'tq~'.and.p3.eq.'ga ')return
      if(p1.eq.'ga '.and.p2.eq.'W+ '.and.p3.eq.'W+~')return
      if(p1.eq.'Z0 '.and.p2.eq.'W+ '.and.p3.eq.'W+~')return
      if(p1.eq.'H0 '.and.p2.eq.'H0 '.and.p3.eq.'H0 ')return
      if(p1.eq.'H0 '.and.p2.eq.'W+ '.and.p3.eq.'W+~')return
      if(p1.eq.'H0 '.and.p2.eq.'Z0 '.and.p3.eq.'Z0 ')return
      if(lightfermions.eq.1)then
        if(p1.eq.'el '.and.p2.eq.'el~'.and.p3.eq.'H0 ')return
        if(p1.eq.'mu '.and.p2.eq.'mu~'.and.p3.eq.'H0 ')return
        if(p1.eq.'dq '.and.p2.eq.'dq~'.and.p3.eq.'H0 ')return
        if(p1.eq.'uq '.and.p2.eq.'uq~'.and.p3.eq.'H0 ')return
        if(p1.eq.'sq '.and.p2.eq.'sq~'.and.p3.eq.'H0 ')return
        if(p1.eq.'cq '.and.p2.eq.'cq~'.and.p3.eq.'H0 ')return
      endif
      if(lightfermions.eq.1.or.lightfermions.eq.2)then
        if(p1.eq.'ta '.and.p2.eq.'ta~'.and.p3.eq.'H0 ')return
        if(p1.eq.'bq '.and.p2.eq.'bq~'.and.p3.eq.'H0 ')return
      endif
      if(p1.eq.'tq '.and.p2.eq.'tq~'.and.p3.eq.'H0 ')return
c ew four-particle vertices (auxiliary particles v1,v2)
      if(p1.eq.'ga '.and.p2.eq.'ga '.and.p3.eq.'v1 ')return
      if(p1.eq.'ga '.and.p2.eq.'Z0 '.and.p3.eq.'v1 ')return
      if(p1.eq.'Z0 '.and.p2.eq.'Z0 '.and.p3.eq.'v1 ')return
      if(p1.eq.'W+ '.and.p2.eq.'W+~'.and.p3.eq.'v1 ')return
      if(p1.eq.'W+ '.and.p2.eq.'W+~'.and.p3.eq.'v1~')return
      if(p1.eq.'Z0 '.and.p2.eq.'Z0 '.and.p3.eq.'v2 ')return
      if(p1.eq.'W+ '.and.p2.eq.'W+~'.and.p3.eq.'v2 ')return
      if(p1.eq.'H0 '.and.p2.eq.'H0 '.and.p3.eq.'v2 ')return
      if(p1.eq.'H0 '.and.p2.eq.'H0 '.and.p3.eq.'v2~')return
c QCD vertices
      if(p1.eq.'dq '.and.p2.eq.'dq~'.and.p3.eq.'gl ')return
      if(p1.eq.'uq '.and.p2.eq.'uq~'.and.p3.eq.'gl ')return
      if(p1.eq.'sq '.and.p2.eq.'sq~'.and.p3.eq.'gl ')return
      if(p1.eq.'cq '.and.p2.eq.'cq~'.and.p3.eq.'gl ')return
      if(p1.eq.'bq '.and.p2.eq.'bq~'.and.p3.eq.'gl ')return
      if(p1.eq.'tq '.and.p2.eq.'tq~'.and.p3.eq.'gl ')return
      if(p1.eq.'gl '.and.p2.eq.'gl '.and.p3.eq.'gl ')return
c QCD four-gluon vertex (auxiliary particle v3)
      if(p1.eq.'gl '.and.p2.eq.'gl '.and.p3.eq.'v3 ')return
      if(p1.eq.'gl '.and.p2.eq.'gl '.and.p3.eq.'v3~')return
      vertexg=.false.
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     identifying particles                                       c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function id(binary)
      implicit none
c local variables
      integer id,binary,i1,i2,i3,nexternal
      i3=1
      nexternal=30
      id=0
      do i1=nexternal,1,-1
        i2=binary/2**(i1-1)
        if((i2/2)*2.ne.i2)then
          id=id+i1*i3
          i3=10*i3
        endif
      enddo
      end
