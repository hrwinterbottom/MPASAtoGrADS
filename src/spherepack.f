c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                       .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c ... file sphcom.f
c
c     this file must be loaded with all driver level files
c     in spherepack3.0.  it includes undocumented subroutines
c     called by some or all of the drivers
c
      subroutine dnlfk (m,n,cp)
c
c     cp requires n/2+1 double precision locations
c
      double precision cp,fnum,fden,fnmh,a1,b1,c1,cp2,fnnp1,fnmsq,fk,
     1       t1,t2,pm1,sc10,sc20,sc40
      dimension       cp(1)
      parameter (sc10=1024.d0)
      parameter (sc20=sc10*sc10)
      parameter (sc40=sc20*sc20)
c
      cp(1) = 0.
      ma = iabs(m)
      if(ma .gt. n) return
      if(n-1) 2,3,5
    2 cp(1) = dsqrt(2.d0)
      return
    3 if(ma .ne. 0) go to 4
      cp(1) = dsqrt(1.5d0)
      return
    4 cp(1) = dsqrt(.75d0)
      if(m .eq. -1) cp(1) = -cp(1)
      return
    5 if(mod(n+ma,2) .ne. 0) go to 10
      nmms2 = (n-ma)/2
      fnum = n+ma+1
      fnmh = n-ma+1
      pm1 = 1.d0
      go to 15
   10 nmms2 = (n-ma-1)/2
      fnum = n+ma+2
      fnmh = n-ma+2
      pm1 = -1.d0
c      t1 = 1.
c      t1 = 2.d0**(n-1)
c      t1 = 1.d0/t1
 15   t1 = 1.d0/sc20
      nex = 20
      fden = 2.d0
      if(nmms2 .lt. 1) go to 20
      do 18 i=1,nmms2
      t1 = fnum*t1/fden
      if(t1 .gt. sc20) then
      t1 = t1/sc40
      nex = nex+40
      end if
      fnum = fnum+2.
      fden = fden+2.
   18 continue
   20 t1 = t1/2.d0**(n-1-nex)
      if(mod(ma/2,2) .ne. 0) t1 = -t1
      t2 = 1. 
      if(ma .eq. 0) go to 26
      do 25 i=1,ma
      t2 = fnmh*t2/(fnmh+pm1)
      fnmh = fnmh+2.
   25 continue
   26 cp2 = t1*dsqrt((n+.5d0)*t2)
      fnnp1 = n*(n+1)
      fnmsq = fnnp1-2.d0*ma*ma
      l = (n+1)/2
      if(mod(n,2) .eq. 0 .and. mod(ma,2) .eq. 0) l = l+1
      cp(l) = cp2
      if(m .ge. 0) go to 29
      if(mod(ma,2) .ne. 0) cp(l) = -cp(l)
   29 if(l .le. 1) return
      fk = n
      a1 = (fk-2.)*(fk-1.)-fnnp1
      b1 = 2.*(fk*fk-fnmsq)
      cp(l-1) = b1*cp(l)/a1
   30 l = l-1
      if(l .le. 1) return
      fk = fk-2.
      a1 = (fk-2.)*(fk-1.)-fnnp1
      b1 = -2.*(fk*fk-fnmsq)
      c1 = (fk+1.)*(fk+2.)-fnnp1
      cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
      go to 30
      end
      subroutine dnlft (m,n,theta,cp,pb)
      double precision cp(*),pb,theta,cdt,sdt,cth,sth,chh
      cdt = dcos(theta+theta)
      sdt = dsin(theta+theta)
      nmod=mod(n,2)
      mmod=mod(m,2)
      if(nmod)1,1,2
1     if(mmod)3,3,4
c
c     n even, m even
c
3     kdo=n/2
      pb = .5*cp(1)
      if(n .eq. 0) return
      cth = cdt
      sth = sdt
      do 170 k=1,kdo
c     pb = pb+cp(k+1)*dcos(2*k*theta)
      pb = pb+cp(k+1)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  170 continue
      return
c
c     n even, m odd
c
    4 kdo = n/2
      pb = 0.
      cth = cdt
      sth = sdt
      do 180 k=1,kdo
c     pb = pb+cp(k)*dsin(2*k*theta)
      pb = pb+cp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  180 continue
      return
2     if(mmod)13,13,14
c
c     n odd, m even
c
13    kdo = (n+1)/2
      pb = 0.
      cth = dcos(theta)
      sth = dsin(theta)
      do 190 k=1,kdo
c     pb = pb+cp(k)*dcos((2*k-1)*theta)
      pb = pb+cp(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  190 continue
      return
c
c     n odd, m odd
c
   14 kdo = (n+1)/2
      pb = 0.
      cth = dcos(theta)
      sth = dsin(theta)
      do 200 k=1,kdo
c     pb = pb+cp(k)*dsin((2*k-1)*theta)
      pb = pb+cp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  200 continue
      return
      end
      subroutine dnlftd (m,n,theta,cp,pb)
c
c     computes the derivative of pmn(theta) with respect to theta
c
      dimension cp(1)
      double precision cp,pb,theta,cdt,sdt,cth,sth,chh
      cdt = dcos(theta+theta)
      sdt = dsin(theta+theta)
      nmod=mod(n,2)
      mmod=mod(abs(m),2)
      if(nmod)1,1,2
1     if(mmod)3,3,4
c
c     n even, m even
c
3     kdo=n/2
      pb = 0.d0
      if(n .eq. 0) return
      cth = cdt
      sth = sdt
      do 170 k=1,kdo
c     pb = pb+cp(k+1)*dcos(2*k*theta)
      pb = pb-2.d0*k*cp(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  170 continue
      return
c
c     n even, m odd
c
    4 kdo = n/2
      pb = 0.
      cth = cdt
      sth = sdt
      do 180 k=1,kdo
c     pb = pb+cp(k)*dsin(2*k*theta)
      pb = pb+2.d0*k*cp(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  180 continue
      return
2     if(mmod)13,13,14
c
c     n odd, m even
c
13    kdo = (n+1)/2
      pb = 0.
      cth = dcos(theta)
      sth = dsin(theta)
      do 190 k=1,kdo
c     pb = pb+cp(k)*dcos((2*k-1)*theta)
      pb = pb-(2.d0*k-1)*cp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  190 continue
      return
c
c     n odd, m odd
c
   14 kdo = (n+1)/2
      pb = 0.
      cth = dcos(theta)
      sth = dsin(theta)
      do 200 k=1,kdo
c     pb = pb+cp(k)*dsin((2*k-1)*theta)
      pb = pb+(2.d0*k-1)*cp(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  200 continue
      return
      end
      subroutine legin(mode,l,nlat,m,w,pmn,km)
c     this subroutine computes legendre polynomials for n=m,...,l-1
c     and  i=1,...,late (late=((nlat+mod(nlat,2))/2)gaussian grid
c     in pmn(n+1,i,km) using swarztrauber's recursion formula.
c     the vector w contains quantities precomputed in shigc.
c     legin must be called in the order m=0,1,...,l-1
c     (e.g., if m=10 is sought it must be preceded by calls with
c     m=0,1,2,...,9 in that order)
      dimension w(1),pmn(1)
c     set size of pole to equator gaussian grid
      late = (nlat+mod(nlat,2))/2
c     partition w (set pointers for p0n,p1n,abel,bbel,cbel,pmn)
      i1 = 1+nlat
      i2 = i1+nlat*late
      i3 = i2+nlat*late
      i4 = i3+(2*nlat-l)*(l-1)/2
      i5 = i4+(2*nlat-l)*(l-1)/2
      call legin1(mode,l,nlat,late,m,w(i1),w(i2),w(i3),w(i4),
     1            w(i5),pmn,km)
      return
      end
      subroutine legin1(mode,l,nlat,late,m,p0n,p1n,abel,bbel,cbel,
     1                  pmn,km)
      dimension p0n(nlat,late),p1n(nlat,late)
      dimension abel(1),bbel(1),cbel(1),pmn(nlat,late,3)
      data km0,km1,km2/ 1,2,3/
      save km0,km1,km2
c     define index function used in storing triangular
c     arrays for recursion coefficients (functions of (m,n))
c     for 2.le.m.le.n-1 and 2.le.n.le.l-1
      indx(m,n) = (n-1)*(n-2)/2+m-1
c     for l.le.n.le.nlat and 2.le.m.le.l
      imndx(m,n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1

c     set do loop indices for full or half sphere
      ms = m+1
      ninc = 1
      if (mode.eq.1) then
c     only compute pmn for n-m odd
      ms = m+2
      ninc = 2
      else if (mode.eq.2) then
c     only compute pmn for n-m even
      ms = m+1
      ninc = 2
      end if


      if (m.gt.1) then
      do 100 np1=ms,nlat,ninc
      n = np1-1
      imn = indx(m,n)
      if (n.ge.l) imn = imndx(m,n)
      do 100 i=1,late
      pmn(np1,i,km0) = abel(imn)*pmn(n-1,i,km2)
     1            +bbel(imn)*pmn(n-1,i,km0)
     2            -cbel(imn)*pmn(np1,i,km2)
  100 continue

      else if (m.eq.0) then
      do 101 np1=ms,nlat,ninc
      do 101 i=1,late
      pmn(np1,i,km0) = p0n(np1,i)
  101 continue

      else if (m.eq.1) then
      do 102 np1=ms,nlat,ninc
      do 102 i=1,late
      pmn(np1,i,km0) = p1n(np1,i)
  102 continue
      end if

c     permute column indices
c     km0,km1,km2 store m,m-1,m-2 columns
      kmt = km0
      km0 = km2
      km2 = km1
      km1 = kmt
c     set current m index in output param km
      km = kmt
      return
      end


      subroutine zfin (isym,nlat,nlon,m,z,i3,wzfin)
      dimension       z(1)        ,wzfin(1)
      imid = (nlat+1)/2
      lim = nlat*imid
      mmax = min0(nlat,nlon/2+1)
      labc = ((mmax-2)*(nlat+nlat-mmax-1))/2
      iw1 = lim+1
      iw2 = iw1+lim
      iw3 = iw2+labc
      iw4 = iw3+labc
c
c     the length of wzfin is 2*lim+3*labc
c
      call zfin1 (isym,nlat,m,z,imid,i3,wzfin,wzfin(iw1),wzfin(iw2),
     1            wzfin(iw3),wzfin(iw4))
      return
      end
      subroutine zfin1 (isym,nlat,m,z,imid,i3,zz,z1,a,b,c)
      dimension       z(imid,nlat,3),zz(imid,1),z1(imid,1),
     1                a(1),b(1),c(1)
      save i1,i2
      ihold = i1
      i1 = i2
      i2 = i3
      i3 = ihold
      if(m-1)25,30,35
   25 i1 = 1
      i2 = 2
      i3 = 3
      do 45 np1=1,nlat
      do 45 i=1,imid
      z(i,np1,i3) = zz(i,np1)
   45 continue
      return
   30 do 50 np1=2,nlat
      do 50 i=1,imid
      z(i,np1,i3) = z1(i,np1)
   50 continue
      return
   35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
      if(isym .eq. 1) go to 36
      do 85 i=1,imid
      z(i,m+1,i3) = a(ns)*z(i,m-1,i1)-c(ns)*z(i,m+1,i1)
   85 continue
   36 if(m .eq. nlat-1) return
      if(isym .eq. 2) go to 71
      ns = ns+1
      do 70 i=1,imid
      z(i,m+2,i3) = a(ns)*z(i,m,i1)-c(ns)*z(i,m+2,i1)
   70 continue
   71 nstrt = m+3
      if(isym .eq. 1) nstrt = m+4
      if(nstrt .gt. nlat) go to 80
      nstp = 2
      if(isym .eq. 0) nstp = 1
      do 75 np1=nstrt,nlat,nstp
      ns = ns+nstp
      do 75 i=1,imid
      z(i,np1,i3) = a(ns)*z(i,np1-2,i1)+b(ns)*z(i,np1-2,i3)
     1                              -c(ns)*z(i,np1,i1)
   75 continue
   80 return
      end
      subroutine zfinit (nlat,nlon,wzfin,dwork)
      dimension       wzfin(*)
      double precision dwork(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     the length of wzfin is 3*((l-3)*l+2)/2 + 2*l*imid
c     the length of dwork is nlat+2
c
      call zfini1 (nlat,nlon,imid,wzfin,wzfin(iw1),dwork,
     1                                       dwork(nlat/2+1))
      return
      end
      subroutine zfini1 (nlat,nlon,imid,z,abc,cz,work)
c
c     abc must have 3*((mmax-2)*(nlat+nlat-mmax-1))/2 locations
c     where mmax = min0(nlat,nlon/2+1)
c     cz and work must each have nlat+1 locations
c
      dimension z(imid,nlat,2),abc(1)
      double precision pi,dt,th,zh,cz(*),work(*)
      pi = 4.*datan(1.d0)
      dt = pi/(nlat-1)
      do 160 mp1=1,2
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dnzfk(nlat,m,n,cz,work)
      do 165 i=1,imid
      th = (i-1)*dt
      call dnzft(nlat,m,n,th,cz,zh)
      z(i,np1,mp1) = zh
  165 continue
      z(1,np1,mp1) = .5*z(1,np1,mp1)
  160 continue
      call rabcp(nlat,nlon,abc)
      return
      end
      subroutine dnzfk(nlat,m,n,cz,work)
c
c     dnzfk computes the coefficients in the trigonometric
c     expansion of the z functions that are used in spherical
c     harmonic analysis.
c
      dimension  cz(1),work(1)
c
c     cz and work must both have nlat/2+1 locations
c
      double precision sum,sc1,t1,t2,work,cz
      lc = (nlat+1)/2
      sc1 = 2.d0/float(nlat-1)
      call dnlfk(m,n,work)
      nmod = mod(n,2)
      mmod = mod(m,2)
      if(nmod)1,1,2
1     if(mmod)3,3,4
c
c     n even, m even
c
3     kdo = n/2+1
      do 5 idx=1,lc
      i = idx+idx-2
      sum = work(1)/(1.d0-i*i)
      if(kdo.lt.2) go to 29
      do 6 kp1=2,kdo
      k = kp1-1
      t1 = 1.d0-(k+k+i)**2
      t2 = 1.d0-(k+k-i)**2
8     sum = sum+work(kp1)*(t1+t2)/(t1*t2)
6     continue
29    cz(idx) = sc1*sum
5     continue
      return
c
c     n even, m odd
c
4     kdo = n/2
      do 9 idx=1,lc
      i = idx+idx-2
      sum = 0.
      do 101 k=1,kdo
      t1 = 1.d0-(k+k+i)**2
      t2 = 1.d0-(k+k-i)**2
12    sum=sum+work(k)*(t1-t2)/(t1*t2)
101   continue
      cz(idx) = sc1*sum
9     continue
      return
2     if(mmod)13,13,14
c
c     n odd, m even
c
13    kdo = (n+1)/2
      do 15 idx=1,lc
      i = idx+idx-1
      sum = 0.
      do 16 k=1,kdo
      t1 = 1.d0-(k+k-1+i)**2
      t2 = 1.d0-(k+k-1-i)**2
18    sum=sum+work(k)*(t1+t2)/(t1*t2)
16    continue
      cz(idx)=sc1*sum
15    continue
      return
c
c     n odd, m odd
c
14    kdo = (n+1)/2
      do 19 idx=1,lc
      i = idx+idx-3
      sum=0.
      do 20 k=1,kdo
      t1 = 1.d0-(k+k-1+i)**2
      t2 = 1.d0-(k+k-1-i)**2
22    sum=sum+work(k)*(t1-t2)/(t1*t2)
20    continue
      cz(idx)=sc1*sum
19    continue
      return
      end
      subroutine dnzft(nlat,m,n,th,cz,zh)
      dimension cz(1)
      double precision cz,zh,th,cdt,sdt,cth,sth,chh
      zh = 0.
      cdt = dcos(th+th)
      sdt = dsin(th+th)
      lmod = mod(nlat,2)
      mmod = mod(m,2)
      nmod = mod(n,2)
      if(lmod)20,20,10
   10 lc = (nlat+1)/2
      lq = lc-1
      ls = lc-2
      if(nmod)1,1,2
    1 if(mmod)3,3,4
c
c     nlat odd n even m even
c
    3 zh = .5*(cz(1)+cz(lc)*dcos(2*lq*th))
      cth = cdt
      sth = sdt
      do 201 k=2,lq
c     zh = zh+cz(k)*dcos(2*(k-1)*th)
      zh = zh+cz(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  201 continue
      return
c
c     nlat odd n even m odd
c
    4 cth = cdt
      sth = sdt
      do 202 k=1,ls
c     zh = zh+cz(k+1)*dsin(2*k*th)
      zh = zh+cz(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  202 continue
      return
c
c     nlat odd n odd, m even
c
    2 if(mmod)5,5,6
    5 cth = dcos(th)
      sth = dsin(th)
      do 203 k=1,lq
c     zh = zh+cz(k)*dcos((2*k-1)*th)
      zh = zh+cz(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  203 continue
      return
c
c     nlat odd n odd m odd
c
    6 cth = dcos(th)
      sth = dsin(th)
      do 204 k=1,lq
c     zh = zh+cz(k+1)*dsin((2*k-1)*th)
      zh = zh+cz(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  204 continue
      return
   20 lc = nlat/2
      lq = lc-1
      if(nmod)30,30,80
   30 if(mmod)40,40,60
c
c     nlat even n even m even
c
   40 zh = .5*cz(1)
      cth = cdt
      sth = sdt
      do 50 k=2,lc
c     zh = zh+cz(k)*dcos(2*(k-1)*th)
      zh = zh+cz(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   50 continue
      return
c
c     nlat even n even m odd
c
   60 cth = cdt
      sth = sdt
      do 70 k=1,lq
c     zh = zh+cz(k+1)*dsin(2*k*th)
      zh = zh+cz(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   70 continue
      return
c
c     nlat even n odd m even
c
   80 if(mmod)90,90,110
   90 zh = .5*cz(lc)*dcos((nlat-1)*th)
      cth = dcos(th)
      sth = dsin(th)
      do 100 k=1,lq
c     zh = zh+cz(k)*dcos((2*k-1)*th)
      zh = zh+cz(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  100 continue
      return
c
c     nlat even n odd m odd
c
  110 cth = dcos(th)
      sth = dsin(th)
      do 120 k=1,lq
c     zh = zh+cz(k+1)*dsin((2*k-1)*th)
      zh = zh+cz(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  120 continue
      return
      end
      subroutine alin (isym,nlat,nlon,m,p,i3,walin)
      dimension       p(1)        ,walin(1)
      imid = (nlat+1)/2
      lim = nlat*imid
      mmax = min0(nlat,nlon/2+1)
      labc = ((mmax-2)*(nlat+nlat-mmax-1))/2
      iw1 = lim+1
      iw2 = iw1+lim
      iw3 = iw2+labc
      iw4 = iw3+labc
c
c     the length of walin is ((5*l-7)*l+6)/2
c
      call alin1 (isym,nlat,m,p,imid,i3,walin,walin(iw1),walin(iw2),
     1            walin(iw3),walin(iw4))
      return
      end
      subroutine alin1 (isym,nlat,m,p,imid,i3,pz,p1,a,b,c)
      dimension       p(imid,nlat,3),pz(imid,1),p1(imid,1),
     1                a(1),b(1),c(1)
      save i1,i2
      ihold = i1
      i1 = i2
      i2 = i3
      i3 = ihold
      if(m-1)25,30,35
   25 i1 = 1
      i2 = 2
      i3 = 3
      do 45 np1=1,nlat
      do 45 i=1,imid
      p(i,np1,i3) = pz(i,np1)
   45 continue
      return
   30 do 50 np1=2,nlat
      do 50 i=1,imid
      p(i,np1,i3) = p1(i,np1)
   50 continue
      return
   35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
      if(isym .eq. 1) go to 36
      do 85 i=1,imid
      p(i,m+1,i3) = a(ns)*p(i,m-1,i1)-c(ns)*p(i,m+1,i1)
   85 continue
   36 if(m .eq. nlat-1) return
      if(isym .eq. 2) go to 71
      ns = ns+1
      do 70 i=1,imid
      p(i,m+2,i3) = a(ns)*p(i,m,i1)-c(ns)*p(i,m+2,i1)
   70 continue
   71 nstrt = m+3
      if(isym .eq. 1) nstrt = m+4
      if(nstrt .gt. nlat) go to 80
      nstp = 2
      if(isym .eq. 0) nstp = 1
      do 75 np1=nstrt,nlat,nstp
      ns = ns+nstp
      do 75 i=1,imid
      p(i,np1,i3) = a(ns)*p(i,np1-2,i1)+b(ns)*p(i,np1-2,i3)
     1                              -c(ns)*p(i,np1,i1)
   75 continue
   80 return
      end
      subroutine alinit (nlat,nlon,walin,dwork)
      dimension       walin(*)
      double precision dwork(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     the length of walin is 3*((l-3)*l+2)/2 + 2*l*imid
c     the length of work is nlat+1
c
      call alini1 (nlat,nlon,imid,walin,walin(iw1),dwork)
      return
      end
      subroutine alini1 (nlat,nlon,imid,p,abc,cp)
      dimension p(imid,nlat,2),abc(1),cp(1)
      double precision pi,dt,th,cp,ph
      pi = 4.*datan(1.d0)
      dt = pi/(nlat-1)
      do 160 mp1=1,2
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dnlfk (m,n,cp)
      do 160 i=1,imid
      th = (i-1)*dt
      call dnlft (m,n,th,cp,ph)
      p(i,np1,mp1) = ph
  160 continue
      call rabcp(nlat,nlon,abc)
      return
      end
      subroutine rabcp(nlat,nlon,abc)
c
c     subroutine rabcp computes the coefficients in the recurrence
c     relation for the associated legendre fuctions. array abc
c     must have 3*((mmax-2)*(nlat+nlat-mmax-1))/2 locations.
c
      dimension abc(1)
      mmax = min0(nlat,nlon/2+1)
      labc = ((mmax-2)*(nlat+nlat-mmax-1))/2
      iw1 = labc+1
      iw2 = iw1+labc
      call rabcp1(nlat,nlon,abc,abc(iw1),abc(iw2))
      return
      end
      subroutine rabcp1(nlat,nlon,a,b,c)
c
c     coefficients a, b, and c for computing pbar(m,n,theta) are
c     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
c
      dimension a(1),b(1),c(1)
      mmax = min0(nlat,nlon/2+1)
      do 215 mp1=3,mmax
      m = mp1-1
      ns = ((m-2)*(nlat+nlat-m-1))/2+1
      fm = float(m)
      tm = fm+fm
      temp = tm*(tm-1.)
      a(ns) = sqrt((tm+1.)*(tm-2.)/temp)
      c(ns) = sqrt(2./temp)
      if(m .eq. nlat-1) go to 215
      ns = ns+1
      temp = tm*(tm+1.)
      a(ns) = sqrt((tm+3.)*(tm-2.)/temp)
      c(ns) = sqrt(6./temp)
      mp3 = m+3
      if(mp3 .gt. nlat) go to 215
      do 210 np1=mp3,nlat
      n = np1-1
      ns = ns+1
      fn = float(n)
      tn = fn+fn
      cn = (tn+1.)/(tn-3.)
      fnpm = fn+fm
      fnmm = fn-fm
      temp = fnpm*(fnpm-1.)
      a(ns) = sqrt(cn*(fnpm-3.)*(fnpm-2.)/temp)
      b(ns) = sqrt(cn*fnmm*(fnmm-1.)/temp)
      c(ns) = sqrt((fnmm+1.)*(fnmm+2.)/temp)
  210 continue
  215 continue
      return
      end
      subroutine sea1(nlat,nlon,imid,z,idz,zin,wzfin,dwork)
      dimension z(idz,*),zin(imid,nlat,3),wzfin(*)
      double precision dwork(*)
      call zfinit(nlat,nlon,wzfin,dwork)
      mmax = min0(nlat,nlon/2+1)
      do 33 mp1=1,mmax
      m = mp1-1
      call zfin (0,nlat,nlon,m,zin,i3,wzfin)
      do 33 np1=mp1,nlat
      mn = m*(nlat-1)-(m*(m-1))/2+np1
      do 33 i=1,imid
      z(mn,i) = zin(i,np1,i3)
   33 continue
      return
      end
      subroutine ses1(nlat,nlon,imid,p,pin,walin,dwork)
      dimension p(imid,*),pin(imid,nlat,3),walin(*)
      double precision dwork(*)
      call alinit (nlat,nlon,walin,dwork)
      mmax = min0(nlat,nlon/2+1)
      do 10 mp1=1,mmax
      m = mp1-1
      call alin(0,nlat,nlon,m,pin,i3,walin)
      do 10 np1=mp1,nlat
      mn = m*(nlat-1)-(m*(m-1))/2+np1
      do 10 i=1,imid
      p(i,mn) = pin(i,np1,i3)
   10 continue
      return
      end
      subroutine zvinit (nlat,nlon,wzvin,dwork)
      dimension       wzvin(1)
      double precision dwork(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     the length of wzvin is 
c         2*nlat*imid +3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     the length of dwork is nlat+2
c
      call zvini1 (nlat,nlon,imid,wzvin,wzvin(iw1),dwork,
     1                                    dwork(nlat/2+2))
      return
      end
      subroutine zvini1 (nlat,nlon,imid,zv,abc,czv,work)
c
c     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     locations where mmax = min0(nlat,(nlon+1)/2)
c     czv and work must each have nlat/2+1  locations
c
      dimension zv(imid,nlat,2),abc(1)
      double precision pi,dt,czv(1),zvh,th,work(1)
      pi = 4.*datan(1.d0)
      dt = pi/(nlat-1)
      mdo = min0(2,nlat,(nlon+1)/2)
      do 160 mp1=1,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dzvk(nlat,m,n,czv,work)
      do 165 i=1,imid
      th = (i-1)*dt
      call dzvt(nlat,m,n,th,czv,zvh)
      zv(i,np1,mp1) = zvh
  165 continue
      zv(1,np1,mp1) = .5*zv(1,np1,mp1)
  160 continue
      call rabcv(nlat,nlon,abc)
      return
      end
      subroutine zwinit (nlat,nlon,wzwin,dwork)
      dimension       wzwin(1)
      double precision dwork(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     the length of wzvin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of dwork is nlat+2
c
      call zwini1 (nlat,nlon,imid,wzwin,wzwin(iw1),dwork,
     1                                        dwork(nlat/2+2))
      return
      end
      subroutine zwini1 (nlat,nlon,imid,zw,abc,czw,work)
c
c     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     locations where mmax = min0(nlat,(nlon+1)/2)
c     czw and work must each have nlat+1 locations
c
      dimension zw(imid,nlat,2),abc(1)
      double precision  pi,dt,czw(1),zwh,th,work(1)
      pi = 4.*datan(1.d0)
      dt = pi/(nlat-1)
      mdo = min0(3,nlat,(nlon+1)/2)
      if(mdo .lt. 2) return
      do 160 mp1=2,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dzwk(nlat,m,n,czw,work)
      do 165 i=1,imid
      th = (i-1)*dt
      call dzwt(nlat,m,n,th,czw,zwh)
      zw(i,np1,m) = zwh
  165 continue
      zw(1,np1,m) = .5*zw(1,np1,m)
  160 continue
      call rabcw(nlat,nlon,abc)
      return
      end
      subroutine zvin (ityp,nlat,nlon,m,zv,i3,wzvin)
      dimension       zv(1)        ,wzvin(1)
      imid = (nlat+1)/2
      lim = nlat*imid
      mmax = min0(nlat,(nlon+1)/2)
      labc = (max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      iw1 = lim+1
      iw2 = iw1+lim
      iw3 = iw2+labc
      iw4 = iw3+labc
c
c     the length of wzvin is 2*lim+3*labc
c
      call zvin1 (ityp,nlat,m,zv,imid,i3,wzvin,wzvin(iw1),wzvin(iw2),
     1            wzvin(iw3),wzvin(iw4))
      return
      end
      subroutine zvin1 (ityp,nlat,m,zv,imid,i3,zvz,zv1,a,b,c)
      dimension       zv(imid,nlat,3),zvz(imid,1),zv1(imid,1),
     1                a(1),b(1),c(1)
      save i1,i2
      ihold = i1
      i1 = i2
      i2 = i3
      i3 = ihold
      if(m-1)25,30,35
   25 i1 = 1
      i2 = 2
      i3 = 3
      do 45 np1=1,nlat
      do 45 i=1,imid
      zv(i,np1,i3) = zvz(i,np1)
   45 continue
      return
   30 do 50 np1=2,nlat
      do 50 i=1,imid
      zv(i,np1,i3) = zv1(i,np1)
   50 continue
      return
   35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
      if(ityp .eq. 1) go to 36
      do 85 i=1,imid
      zv(i,m+1,i3) = a(ns)*zv(i,m-1,i1)-c(ns)*zv(i,m+1,i1)
   85 continue
   36 if(m .eq. nlat-1) return
      if(ityp .eq. 2) go to 71
      ns = ns+1
      do 70 i=1,imid
      zv(i,m+2,i3) = a(ns)*zv(i,m,i1)-c(ns)*zv(i,m+2,i1)
   70 continue
   71 nstrt = m+3
      if(ityp .eq. 1) nstrt = m+4
      if(nstrt .gt. nlat) go to 80
      nstp = 2
      if(ityp .eq. 0) nstp = 1
      do 75 np1=nstrt,nlat,nstp
      ns = ns+nstp
      do 75 i=1,imid
      zv(i,np1,i3) = a(ns)*zv(i,np1-2,i1)+b(ns)*zv(i,np1-2,i3)
     1                              -c(ns)*zv(i,np1,i1)
   75 continue
   80 return
      end
      subroutine zwin (ityp,nlat,nlon,m,zw,i3,wzwin)
      dimension       zw(1)        ,wzwin(1)
      imid = (nlat+1)/2
      lim = nlat*imid
      mmax = min0(nlat,(nlon+1)/2)
      labc = (max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      iw1 = lim+1
      iw2 = iw1+lim
      iw3 = iw2+labc
      iw4 = iw3+labc
c
c     the length of wzwin is 2*lim+3*labc
c
      call zwin1 (ityp,nlat,m,zw,imid,i3,wzwin,wzwin(iw1),wzwin(iw2),
     1            wzwin(iw3),wzwin(iw4))
      return
      end
      subroutine zwin1 (ityp,nlat,m,zw,imid,i3,zw1,zw2,a,b,c)
      dimension       zw(imid,nlat,3),zw1(imid,1),zw2(imid,1),
     1                a(1),b(1),c(1)
      save i1,i2
      ihold = i1
      i1 = i2
      i2 = i3
      i3 = ihold
      if(m-2)25,30,35
   25 i1 = 1
      i2 = 2
      i3 = 3
      do 45 np1=2,nlat
      do 45 i=1,imid
      zw(i,np1,i3) = zw1(i,np1)
   45 continue
      return
   30 do 50 np1=3,nlat
      do 50 i=1,imid
      zw(i,np1,i3) = zw2(i,np1)
   50 continue
      return
   35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
      if(ityp .eq. 1) go to 36
      do 85 i=1,imid
      zw(i,m+1,i3) = a(ns)*zw(i,m-1,i1)-c(ns)*zw(i,m+1,i1)
   85 continue
   36 if(m .eq. nlat-1) return
      if(ityp .eq. 2) go to 71
      ns = ns+1
      do 70 i=1,imid
      zw(i,m+2,i3) = a(ns)*zw(i,m,i1)-c(ns)*zw(i,m+2,i1)
   70 continue
   71 nstrt = m+3
      if(ityp .eq. 1) nstrt = m+4
      if(nstrt .gt. nlat) go to 80
      nstp = 2
      if(ityp .eq. 0) nstp = 1
      do 75 np1=nstrt,nlat,nstp
      ns = ns+nstp
      do 75 i=1,imid
      zw(i,np1,i3) = a(ns)*zw(i,np1-2,i1)+b(ns)*zw(i,np1-2,i3)
     1                              -c(ns)*zw(i,np1,i1)
   75 continue
   80 return
      end
      subroutine vbinit (nlat,nlon,wvbin,dwork)
      dimension wvbin(1)
      double precision dwork(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of dwork is nlat+2
c
      call vbini1 (nlat,nlon,imid,wvbin,wvbin(iw1),dwork,
     1                                       dwork(nlat/2+2))
      return
      end
      subroutine vbini1 (nlat,nlon,imid,vb,abc,cvb,work)
c
c     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     locations where mmax = min0(nlat,(nlon+1)/2)
c     cvb and work must each have nlat+1 locations
c
      dimension vb(imid,nlat,2),abc(1)
      double precision pi,dt,cvb(1),th,vbh,work(1)
      pi = 4.*datan(1.d0)
      dt = pi/(nlat-1)
      mdo = min0(2,nlat,(nlon+1)/2)
      do 160 mp1=1,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dvbk(m,n,cvb,work)
      do 165 i=1,imid
      th = (i-1)*dt
      call dvbt(m,n,th,cvb,vbh)
      vb(i,np1,mp1) = vbh
  165 continue
  160 continue
      call rabcv(nlat,nlon,abc)
      return
      end
      subroutine wbinit (nlat,nlon,wwbin,dwork)
      dimension       wwbin(1)
      double precision dwork(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of dwork is nlat+2
c
      call wbini1 (nlat,nlon,imid,wwbin,wwbin(iw1),dwork,
     1                                        dwork(nlat/2+2))
      return
      end
      subroutine wbini1 (nlat,nlon,imid,wb,abc,cwb,work)
c
c     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     locations where mmax = min0(nlat,(nlon+1)/2)
c     cwb and work must each have nlat/2+1 locations
c
      dimension wb(imid,nlat,2),abc(1)
      double precision pi,dt,cwb(1),wbh,th,work(1)
      pi = 4.*datan(1.d0)
      dt = pi/(nlat-1)
      mdo = min0(3,nlat,(nlon+1)/2)
      if(mdo .lt. 2) return
      do 160 mp1=2,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dwbk(m,n,cwb,work)
      do 165 i=1,imid
      th = (i-1)*dt
      call dwbt(m,n,th,cwb,wbh)
      wb(i,np1,m) = wbh
  165 continue
  160 continue
      call rabcw(nlat,nlon,abc)
      return
      end
      subroutine vbin (ityp,nlat,nlon,m,vb,i3,wvbin)
      dimension       vb(1)        ,wvbin(1)
      imid = (nlat+1)/2
      lim = nlat*imid
      mmax = min0(nlat,(nlon+1)/2)
      labc = (max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      iw1 = lim+1
      iw2 = iw1+lim
      iw3 = iw2+labc
      iw4 = iw3+labc
c
c     the length of wvbin is 2*lim+3*labc
c
      call vbin1 (ityp,nlat,m,vb,imid,i3,wvbin,wvbin(iw1),wvbin(iw2),
     1            wvbin(iw3),wvbin(iw4))
      return
      end
      subroutine vbin1 (ityp,nlat,m,vb,imid,i3,vbz,vb1,a,b,c)
      dimension       vb(imid,nlat,3),vbz(imid,1),vb1(imid,1),
     1                a(1),b(1),c(1)
      save i1,i2
      ihold = i1
      i1 = i2
      i2 = i3
      i3 = ihold
      if(m-1)25,30,35
   25 i1 = 1
      i2 = 2
      i3 = 3
      do 45 np1=1,nlat
      do 45 i=1,imid
      vb(i,np1,i3) = vbz(i,np1)
   45 continue
      return
   30 do 50 np1=2,nlat
      do 50 i=1,imid
      vb(i,np1,i3) = vb1(i,np1)
   50 continue
      return
   35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
      if(ityp .eq. 1) go to 36
      do 85 i=1,imid
      vb(i,m+1,i3) = a(ns)*vb(i,m-1,i1)-c(ns)*vb(i,m+1,i1)
   85 continue
   36 if(m .eq. nlat-1) return
      if(ityp .eq. 2) go to 71
      ns = ns+1
      do 70 i=1,imid
      vb(i,m+2,i3) = a(ns)*vb(i,m,i1)-c(ns)*vb(i,m+2,i1)
   70 continue
   71 nstrt = m+3
      if(ityp .eq. 1) nstrt = m+4
      if(nstrt .gt. nlat) go to 80
      nstp = 2
      if(ityp .eq. 0) nstp = 1
      do 75 np1=nstrt,nlat,nstp
      ns = ns+nstp
      do 75 i=1,imid
      vb(i,np1,i3) = a(ns)*vb(i,np1-2,i1)+b(ns)*vb(i,np1-2,i3)
     1                              -c(ns)*vb(i,np1,i1)
   75 continue
   80 return
      end
      subroutine wbin (ityp,nlat,nlon,m,wb,i3,wwbin)
      dimension       wb(1)        ,wwbin(1)
      imid = (nlat+1)/2
      lim = nlat*imid
      mmax = min0(nlat,(nlon+1)/2)
      labc = (max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      iw1 = lim+1
      iw2 = iw1+lim
      iw3 = iw2+labc
      iw4 = iw3+labc
c
c     the length of wwbin is 2*lim+3*labc
c
      call wbin1 (ityp,nlat,m,wb,imid,i3,wwbin,wwbin(iw1),wwbin(iw2),
     1            wwbin(iw3),wwbin(iw4))
      return
      end
      subroutine wbin1 (ityp,nlat,m,wb,imid,i3,wb1,wb2,a,b,c)
      dimension       wb(imid,nlat,3),wb1(imid,1),wb2(imid,1),
     1                a(1),b(1),c(1)
      save i1,i2
      ihold = i1
      i1 = i2
      i2 = i3
      i3 = ihold
      if(m-2)25,30,35
   25 i1 = 1
      i2 = 2
      i3 = 3
      do 45 np1=2,nlat
      do 45 i=1,imid
      wb(i,np1,i3) = wb1(i,np1)
   45 continue
      return
   30 do 50 np1=3,nlat
      do 50 i=1,imid
      wb(i,np1,i3) = wb2(i,np1)
   50 continue
      return
   35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
      if(ityp .eq. 1) go to 36
      do 85 i=1,imid
      wb(i,m+1,i3) = a(ns)*wb(i,m-1,i1)-c(ns)*wb(i,m+1,i1)
   85 continue
   36 if(m .eq. nlat-1) return
      if(ityp .eq. 2) go to 71
      ns = ns+1
      do 70 i=1,imid
      wb(i,m+2,i3) = a(ns)*wb(i,m,i1)-c(ns)*wb(i,m+2,i1)
   70 continue
   71 nstrt = m+3
      if(ityp .eq. 1) nstrt = m+4
      if(nstrt .gt. nlat) go to 80
      nstp = 2
      if(ityp .eq. 0) nstp = 1
      do 75 np1=nstrt,nlat,nstp
      ns = ns+nstp
      do 75 i=1,imid
      wb(i,np1,i3) = a(ns)*wb(i,np1-2,i1)+b(ns)*wb(i,np1-2,i3)
     1                              -c(ns)*wb(i,np1,i1)
   75 continue
   80 return
      end
       subroutine dzvk(nlat,m,n,czv,work)
c
c     subroutine dzvk computes the coefficients in the trigonometric
c     expansion of the quadrature function zvbar(n,m,theta)
c
c     input parameters
c
c     nlat      the number of colatitudes including the poles.
c
c     n      the degree (subscript) of wbarv(n,m,theta)
c
c     m      the order (superscript) of wbarv(n,m,theta)
c
c     work   a work array with at least nlat/2+1 locations
c
c     output parameter
c
c     czv     the fourier coefficients of zvbar(n,m,theta).
c
      dimension czv(1),work(1)
      double precision czv,sc1,sum,work,t1,t2
      if(n .le. 0) return
      lc = (nlat+1)/2
      sc1 = 2.d0/float(nlat-1)
      call dvbk(m,n,work,czv)
      nmod = mod(n,2)
      mmod = mod(m,2)
      if(nmod .ne. 0) go to 1
      if(mmod .ne. 0) go to 2
c
c     n even, m even
c
      kdo = n/2
      do 9 id=1,lc
      i = id+id-2
      sum = 0.
      do 10 k=1,kdo
      t1 = 1.d0-(k+k+i)**2
      t2 = 1.d0-(k+k-i)**2
      sum = sum+work(k)*(t1-t2)/(t1*t2)
   10 continue
      czv(id) = sc1*sum
    9 continue
      return
c
c     n even, m odd
c
    2 kdo = n/2
      do 5 id=1,lc
      i = id+id-2
      sum = 0.
      do 6 k=1,kdo
      t1 = 1.d0-(k+k+i)**2
      t2 = 1.d0-(k+k-i)**2
      sum = sum+work(k)*(t1+t2)/(t1*t2)
    6 continue
      czv(id) = sc1*sum
    5 continue
      return
    1 if(mmod .ne. 0) go to 3
c
c     n odd, m even
c
      kdo = (n+1)/2
      do 19 id=1,lc
      i = id+id-3
      sum = 0.
      do 20 k=1,kdo
      t1 = 1.d0-(k+k-1+i)**2
      t2 = 1.d0-(k+k-1-i)**2
      sum = sum+work(k)*(t1-t2)/(t1*t2)
   20 continue
      czv(id) = sc1*sum
   19 continue
      return
c
c     n odd, m odd
c
    3 kdo = (n+1)/2
      do 15 id=1,lc
      i = id+id-1
      sum = 0.
      do 16 k=1,kdo
      t1 = 1.d0-(k+k-1+i)**2
      t2 = 1.d0-(k+k-1-i)**2
      sum = sum+work(k)*(t1+t2)/(t1*t2)
   16 continue
      czv(id) = sc1*sum
   15 continue
      return
      end
      subroutine dzvt(nlat,m,n,th,czv,zvh)
c
c     subroutine dzvt tabulates the function zvbar(n,m,theta)
c     at theta = th in double precision
c
c     input parameters
c
c     nlat      the number of colatitudes including the poles.
c
c     n      the degree (subscript) of zvbar(n,m,theta)
c
c     m      the order (superscript) of zvbar(n,m,theta)
c
c     czv     the fourier coefficients of zvbar(n,m,theta)
c             as computed by subroutine zwk.
c
c     output parameter
c
c     zvh     zvbar(m,n,theta) evaluated at theta = th
c
      dimension czv(1)
      double precision th,czv,zvh,cth,sth,cdt,sdt,chh
      zvh = 0.
      if(n .le. 0) return
      lc = (nlat+1)/2
      lq = lc-1
      ls = lc-2
      cth = dcos(th)
      sth = dsin(th)
      cdt = cth*cth-sth*sth
      sdt = 2.*sth*cth
      lmod = mod(nlat,2)
      mmod = mod(m,2)
      nmod = mod(n,2)
      if(lmod .eq. 0) go to 50
      if(nmod .ne. 0) go to 1
      cth = cdt
      sth = sdt
      if(mmod .ne. 0) go to 2
c
c     nlat odd  n even  m even
c
      do 10 k=1,ls
      zvh = zvh+czv(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   10 continue
      return
c
c     nlat odd  n even  m odd
c
    2 zvh = .5*czv(1)
      do 20 k=2,lq
      zvh = zvh+czv(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   20 continue
      zvh = zvh+.5*czv(lc)*dcos((nlat-1)*th)
      return
    1 if(mmod .ne. 0) go to 3
c
c     nlat odd  n odd  m even
c
      do 30 k=1,lq
      zvh = zvh+czv(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   30 continue
      return
c
c     nlat odd  n odd  m odd
c
    3 do 40 k=1,lq
      zvh = zvh+czv(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   40 continue
      return
   50 if(nmod .ne. 0) go to 51
      cth = cdt
      sth = sdt
      if(mmod .ne. 0) go to 52
c
c     nlat even  n even  m even
c
      do 55 k=1,lq
      zvh = zvh+czv(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   55 continue
      return
c
c     nlat even  n even  m odd
c
   52 zvh = .5*czv(1)
      do 57 k=2,lc
      zvh = zvh+czv(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   57 continue
      return
   51 if(mmod .ne. 0) go to 53
c
c     nlat even  n odd  m even
c
      do 58 k=1,lq
      zvh = zvh+czv(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   58 continue
      return
c
c     nlat even  n odd  m odd
c
   53 zvh = .5*czv(lc)*dcos((nlat-1)*th)
      do 60 k=1,lq
      zvh = zvh+czv(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   60 continue
      return
      end
      subroutine dzwk(nlat,m,n,czw,work)
c
c     subroutine dzwk computes the coefficients in the trigonometric
c     expansion of the quadrature function zwbar(n,m,theta)
c
c     input parameters
c
c     nlat      the number of colatitudes including the poles.
c
c     n      the degree (subscript) of zwbar(n,m,theta)
c
c     m      the order (superscript) of zwbar(n,m,theta)
c
c     work   a work array with at least nlat/2+1 locations
c
c     output parameter
c
c     czw     the fourier coefficients of zwbar(n,m,theta).
c
      dimension czw(1),work(1)
      double precision czw,work,sc1,sum,t1,t2
      if(n .le. 0) return
      lc = (nlat+1)/2
      sc1 = 2.d0/float(nlat-1)
      call dwbk(m,n,work,czw)
      nmod = mod(n,2)
      mmod = mod(m,2)
      if(nmod .ne. 0) go to 1
      if(mmod .ne. 0) go to 2
c
c     n even, m even
c
      kdo = n/2
      do 19 id=1,lc
      i = id+id-3
      sum = 0.
      do 20 k=1,kdo
      t1 = 1.d0-(k+k-1+i)**2
      t2 = 1.d0-(k+k-1-i)**2
      sum = sum+work(k)*(t1-t2)/(t1*t2)
   20 continue
      czw(id) = sc1*sum
   19 continue
      return
c
c     n even, m odd
c
    2 kdo = n/2
      do 15 id=1,lc
      i = id+id-1
      sum = 0.
      do 16 k=1,kdo
      t1 = 1.d0-(k+k-1+i)**2
      t2 = 1.d0-(k+k-1-i)**2
      sum = sum+work(k)*(t1+t2)/(t1*t2)
   16 continue
      czw(id) = sc1*sum
   15 continue
      return
    1 if(mmod .ne. 0) go to 3
c
c     n odd, m even
c
      kdo = (n-1)/2
      do 9 id=1,lc
      i = id+id-2
      sum = 0.
      do 10 k=1,kdo
      t1 = 1.d0-(k+k+i)**2
      t2 = 1.d0-(k+k-i)**2
      sum = sum+work(k)*(t1-t2)/(t1*t2)
   10 continue
      czw(id) = sc1*sum
    9 continue
      return
c
c     n odd, m odd
c
    3 kdo = (n+1)/2
      do 5 id=1,lc
      i = id+id-2
      sum = work(1)/(1.d0-i*i)
      if(kdo .lt. 2) go to 29
      do 6 kp1=2,kdo
      k = kp1-1
      t1 = 1.d0-(k+k+i)**2
      t2 = 1.d0-(k+k-i)**2
      sum = sum+work(kp1)*(t1+t2)/(t1*t2)
    6 continue
   29 czw(id) = sc1*sum
    5 continue
      return
      end
      subroutine dzwt(nlat,m,n,th,czw,zwh)
c
c     subroutine dzwt tabulates the function zwbar(n,m,theta)
c     at theta = th in double precision
c
c     input parameters
c
c     nlat      the number of colatitudes including the poles.
c            nlat must be an odd integer
c
c     n      the degree (subscript) of zwbar(n,m,theta)
c
c     m      the order (superscript) of zwbar(n,m,theta)
c
c     czw     the fourier coefficients of zwbar(n,m,theta)
c             as computed by subroutine zwk.
c
c     output parameter
c
c     zwh     zwbar(m,n,theta) evaluated at theta = th
c
      dimension czw(1)
      double precision czw,zwh,th,cth,sth,cdt,sdt,chh
      zwh = 0.
      if(n .le. 0) return
      lc = (nlat+1)/2
      lq = lc-1
      ls = lc-2
      cth = dcos(th)
      sth = dsin(th)
      cdt = cth*cth-sth*sth
      sdt = 2.*sth*cth
      lmod = mod(nlat,2)
      mmod = mod(m,2)
      nmod = mod(n,2)
      if(lmod .eq. 0) go to 50
      if(nmod .ne. 0) go to 1
      if(mmod .ne. 0) go to 2
c
c     nlat odd  n even  m even
c
      do 30 k=1,lq
      zwh = zwh+czw(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   30 continue
      return
c
c     nlat odd  n even  m odd
c
    2 do 40 k=1,lq
      zwh = zwh+czw(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   40 continue
      return
    1 cth = cdt
      sth = sdt
      if(mmod .ne. 0) go to 3
c
c     nlat odd  n odd  m even
c
      do 10 k=1,ls
      zwh = zwh+czw(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   10 continue
      return
c
c     nlat odd  n odd  m odd
c
    3 zwh = .5*czw(1)
      do 20 k=2,lq
      zwh = zwh+czw(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   20 continue
      zwh = zwh+.5*czw(lc)*dcos((nlat-1)*th)
      return
   50 if(nmod .ne. 0) go to 51
      if(mmod .ne. 0) go to 52
c
c     nlat even  n even  m even
c
      do 55 k=1,lq
      zwh = zwh+czw(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   55 continue
      return
c
c     nlat even  n even  m odd
c
   52 zwh = .5*czw(lc)*dcos((nlat-1)*th)
      do 60 k=1,lq
      zwh = zwh+czw(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   60 continue
      return
   51 cth = cdt
      sth = sdt
      if(mmod .ne. 0) go to 53
c
c     nlat even  n odd  m even
c
      do 65 k=1,lq
      zwh = zwh+czw(k+1)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   65 continue
      return
c
c     nlat even  n odd  m odd
c
   53 zwh = .5*czw(1)
      do 70 k=2,lc
      zwh = zwh+czw(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   70 continue
      return
      end
      subroutine dvbk(m,n,cv,work)
      double precision cv(1),work(1),fn,fk,cf
      cv(1) = 0.
      if(n .le. 0) return
      fn = n
      srnp1 = dsqrt(fn*(fn+1.))
      cf = 2.*m/srnp1
      modn = mod(n,2)
      modm = mod(m,2)
      call dnlfk(m,n,work)
      if(modn .ne. 0) go to 70
      ncv = n/2
      if(ncv .eq. 0) return
      fk = 0.
      if(modm .ne. 0) go to 60
c
c     n even m even
c
      do 55 l=1,ncv
      fk = fk+2.
      cv(l) = -fk*work(l+1)/srnp1
   55 continue
      return
c
c     n even m odd
c
   60 do 65 l=1,ncv
      fk = fk+2.
      cv(l) = fk*work(l)/srnp1
   65 continue
      return
   70 ncv = (n+1)/2
      fk = -1.
      if(modm .ne. 0) go to 80
c
c     n odd m even
c
      do 75 l=1,ncv
      fk = fk+2.
      cv(l) = -fk*work(l)/srnp1
   75 continue
      return
c
c     n odd m odd
c
   80 do 85 l=1,ncv
      fk = fk+2.
      cv(l) = fk*work(l)/srnp1
   85 continue
      return
      end
      subroutine dwbk(m,n,cw,work)
      double precision cw(1),work(1),fn,cf,srnp1
      cw(1) = 0.
      if(n.le.0 .or. m.le.0) return
      fn = n
      srnp1 = dsqrt(fn*(fn+1.))
      cf = 2.*m/srnp1
      modn = mod(n,2)
      modm = mod(m,2)
      call dnlfk(m,n,work)
      if(m .eq. 0) go to 50
      if(modn .ne. 0) go to 30
      l = n/2
      if(l .eq. 0) go to 50
      if(modm .ne. 0) go to 20
c
c     n even m even
c
      cw(l) = -cf*work(l+1)
   10 l = l-1
      if(l .le. 0) go to 50
      cw(l) = cw(l+1)-cf*work(l+1)
      go to 10
c
c     n even m odd
c
   20 cw(l) = cf*work(l)
   25 l = l-1
      if(l .le. 0) go to 50
      cw(l) = cw(l+1)+cf*work(l)
      go to 25
   30 if(modm .ne. 0) go to 40
      l = (n-1)/2
      if(l .eq. 0) go to 50
c
c     n odd m even
c
      cw(l) = -cf*work(l+1)
   35 l = l-1
      if(l .le. 0) go to 50
      cw(l) = cw(l+1)-cf*work(l+1)
      go to 35
c
c     n odd m odd
c
   40 l = (n+1)/2
      cw(l) = cf*work(l)
   45 l = l-1
      if(l .le. 0) go to 50
      cw(l) = cw(l+1)+cf*work(l)
      go to 45
   50 return
      end
      subroutine dvbt(m,n,theta,cv,vh)
      dimension cv(1)
      double precision cv,vh,theta,cth,sth,cdt,sdt,chh
      vh = 0.
      if(n.eq.0) return
      cth = dcos(theta)
      sth = dsin(theta)
      cdt = cth*cth-sth*sth
      sdt = 2.*sth*cth
      mmod = mod(m,2)
      nmod = mod(n,2)
      if(nmod .ne. 0) go to 1
      cth = cdt
      sth = sdt
      if(mmod .ne. 0) go to 2
c
c     n even  m even
c
      ncv = n/2
      do 10 k=1,ncv
      vh = vh+cv(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   10 continue
      return
c
c     n even  m odd
c
    2 ncv = n/2
      do 15 k=1,ncv
      vh = vh+cv(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   15 continue
      return
    1 if(mmod .ne. 0) go to 3
c
c     n odd m even
c
      ncv = (n+1)/2
      do 20 k=1,ncv
      vh = vh+cv(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   20 continue
      return
c
c case m odd and n odd
c
    3 ncv = (n+1)/2
      do 25 k=1,ncv
      vh = vh+cv(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   25 continue
      return
      end
      subroutine dwbt(m,n,theta,cw,wh)
      dimension cw(1)
      double precision theta,cw,wh,cth,sth,cdt,sdt,chh
      wh = 0.
      if(n.le.0 .or. m.le.0) return
      cth = dcos(theta)
      sth = dsin(theta)
      cdt = cth*cth-sth*sth
      sdt = 2.*sth*cth
      mmod=mod(m,2)
      nmod=mod(n,2)
      if(nmod .ne. 0) go to 1
      if(mmod .ne. 0) go to 2
c
c     n even  m even
c
      ncw = n/2
      do 10 k=1,ncw
      wh = wh+cw(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   10 continue
      return
c
c     n even  m odd
c
    2 ncw = n/2
      do 8 k=1,ncw
      wh = wh+cw(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
    8 continue
      return
    1 cth = cdt
      sth = sdt
      if(mmod .ne. 0) go to 3
c
c     n odd m even
c
      ncw = (n-1)/2
      do 20 k=1,ncw
      wh = wh+cw(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   20 continue
      return
c
c case m odd and n odd
c
    3 ncw = (n+1)/2
      wh = .5*cw(1)
      if(ncw.lt.2) return
      do 25 k=2,ncw
      wh = wh+cw(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   25 continue
      return
      end
      subroutine rabcv(nlat,nlon,abc)
c
c     subroutine rabcp computes the coefficients in the recurrence
c     relation for the functions vbar(m,n,theta). array abc
c     must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 locations.
c
      dimension abc(1)
      mmax = min0(nlat,(nlon+1)/2)
      labc = (max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      iw1 = labc+1
      iw2 = iw1+labc
      call rabcv1(nlat,nlon,abc,abc(iw1),abc(iw2))
      return
      end
      subroutine rabcv1(nlat,nlon,a,b,c)
c
c     coefficients a, b, and c for computing vbar(m,n,theta) are
c     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
c
      dimension a(1),b(1),c(1)
      mmax = min0(nlat,(nlon+1)/2)
      if(mmax .lt. 3) return
      do 215 mp1=3,mmax
      m = mp1-1
      ns = ((m-2)*(nlat+nlat-m-1))/2+1
      fm = float(m)
      tm = fm+fm
      temp = tm*(tm-1.)
      tpn = (fm-2.)*(fm-1.)/(fm*(fm+1.))
      a(ns) = sqrt(tpn*(tm+1.)*(tm-2.)/temp)
      c(ns) = sqrt(2./temp)
      if(m .eq. nlat-1) go to 215
      ns = ns+1
      temp = tm*(tm+1.)
      tpn = (fm-1.)*fm/((fm+1.)*(fm+2.))
      a(ns) = sqrt(tpn*(tm+3.)*(tm-2.)/temp)
      c(ns) = sqrt(6./temp)
      mp3 = m+3
      if(mp3 .gt. nlat) go to 215
      do 210 np1=mp3,nlat
      n = np1-1
      ns = ns+1
      fn = float(n)
      tn = fn+fn
      cn = (tn+1.)/(tn-3.)
      tpn = (fn-2.)*(fn-1.)/(fn*(fn+1.))
      fnpm = fn+fm
      fnmm = fn-fm
      temp = fnpm*(fnpm-1.)
      a(ns) = sqrt(tpn*cn*(fnpm-3.)*(fnpm-2.)/temp)
      b(ns) = sqrt(tpn*cn*fnmm*(fnmm-1.)/temp)
      c(ns) = sqrt((fnmm+1.)*(fnmm+2.)/temp)
  210 continue
  215 continue
      return
      end
      subroutine rabcw(nlat,nlon,abc)
c
c     subroutine rabcw computes the coefficients in the recurrence
c     relation for the functions wbar(m,n,theta). array abc
c     must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 locations.
c
      dimension abc(1)
      mmax = min0(nlat,(nlon+1)/2)
      labc = (max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      iw1 = labc+1
      iw2 = iw1+labc
      call rabcw1(nlat,nlon,abc,abc(iw1),abc(iw2))
      return
      end
      subroutine rabcw1(nlat,nlon,a,b,c)
c
c     coefficients a, b, and c for computing wbar(m,n,theta) are
c     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
c
      dimension a(1),b(1),c(1)
      mmax = min0(nlat,(nlon+1)/2)
      if(mmax .lt. 4) return
      do 215 mp1=4,mmax
      m = mp1-1
      ns = ((m-2)*(nlat+nlat-m-1))/2+1
      fm = float(m)
      tm = fm+fm
      temp = tm*(tm-1.)
      tpn = (fm-2.)*(fm-1.)/(fm*(fm+1.))
      tph = fm/(fm-2.)
      a(ns) = tph*sqrt(tpn*(tm+1.)*(tm-2.)/temp)
      c(ns) = tph*sqrt(2./temp)
      if(m .eq. nlat-1) go to 215
      ns = ns+1
      temp = tm*(tm+1.)
      tpn = (fm-1.)*fm/((fm+1.)*(fm+2.))
      tph = fm/(fm-2.)
      a(ns) = tph*sqrt(tpn*(tm+3.)*(tm-2.)/temp)
      c(ns) = tph*sqrt(6./temp)
      mp3 = m+3
      if(mp3 .gt. nlat) go to 215
      do 210 np1=mp3,nlat
      n = np1-1
      ns = ns+1
      fn = float(n)
      tn = fn+fn
      cn = (tn+1.)/(tn-3.)
      fnpm = fn+fm
      fnmm = fn-fm
      temp = fnpm*(fnpm-1.)
      tpn = (fn-2.)*(fn-1.)/(fn*(fn+1.))
      tph = fm/(fm-2.)
      a(ns) = tph*sqrt(tpn*cn*(fnpm-3.)*(fnpm-2.)/temp)
      b(ns) = sqrt(tpn*cn*fnmm*(fnmm-1.)/temp)
      c(ns) = tph*sqrt((fnmm+1.)*(fnmm+2.)/temp)
  210 continue
  215 continue
      return
      end
      subroutine vtinit (nlat,nlon,wvbin,dwork)
      dimension       wvbin(*)
      double precision dwork(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of dwork is nlat+2
c
      call vtini1 (nlat,nlon,imid,wvbin,wvbin(iw1),dwork,
     1                                       dwork(nlat/2+2))
      return
      end
      subroutine vtini1 (nlat,nlon,imid,vb,abc,cvb,work)
c
c     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     locations where mmax = min0(nlat,(nlon+1)/2)
c     cvb and work must each have nlat/2+1 locations
c
      dimension vb(imid,nlat,2),abc(1),cvb(1)
      double precision pi,dt,cvb,th,vbh,work(*)
      pi = 4.*datan(1.d0)
      dt = pi/(nlat-1)
      mdo = min0(2,nlat,(nlon+1)/2)
      do 160 mp1=1,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dvtk(m,n,cvb,work)
      do 165 i=1,imid
      th = (i-1)*dt
      call dvtt(m,n,th,cvb,vbh)
      vb(i,np1,mp1) = vbh
  165 continue
  160 continue
      call rabcv(nlat,nlon,abc)
      return
      end
      subroutine wtinit (nlat,nlon,wwbin,dwork)
      dimension       wwbin(1)
      double precision dwork(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of dwork is nlat+2
c
      call wtini1 (nlat,nlon,imid,wwbin,wwbin(iw1),dwork,
     1                                       dwork(nlat/2+2))
      return
      end
      subroutine wtini1 (nlat,nlon,imid,wb,abc,cwb,work)
c
c     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     locations where mmax = min0(nlat,(nlon+1)/2)
c     cwb and work must each have nlat/2+1 locations
c
      dimension wb(imid,nlat,2),abc(1)
      double precision pi,dt,cwb(*),wbh,th,work(*)
      pi = 4.*datan(1.d0)
      dt = pi/(nlat-1)
      mdo = min0(3,nlat,(nlon+1)/2)
      if(mdo .lt. 2) return
      do 160 mp1=2,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dwtk(m,n,cwb,work)
      do 165 i=1,imid
      th = (i-1)*dt
      call dwtt(m,n,th,cwb,wbh)
      wb(i,np1,m) = wbh
  165 continue
  160 continue
      call rabcw(nlat,nlon,abc)
      return
      end
      subroutine vtgint (nlat,nlon,theta,wvbin,work)
      dimension       wvbin(*)
      double precision theta(*), work(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     theta is a double precision array with (nlat+1)/2 locations
c     nlat is the maximum value of n+1
c     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of work is nlat+2
c
      call vtgit1 (nlat,nlon,imid,theta,wvbin,wvbin(iw1),
     +                        work,work(nlat/2+2))
      return
      end
      subroutine vtgit1 (nlat,nlon,imid,theta,vb,abc,cvb,work)
c
c     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     locations where mmax = min0(nlat,(nlon+1)/2)
c     cvb and work must each have nlat/2+1   locations
c
      dimension vb(imid,nlat,2),abc(*)
      double precision theta(*),cvb(*),work(*),vbh
      mdo = min0(2,nlat,(nlon+1)/2)
      do 160 mp1=1,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dvtk(m,n,cvb,work)
      do 165 i=1,imid
      call dvtt(m,n,theta(i),cvb,vbh)
      vb(i,np1,mp1) = vbh
  165 continue
  160 continue
      call rabcv(nlat,nlon,abc)
      return
      end
      subroutine wtgint (nlat,nlon,theta,wwbin,work)
      dimension       wwbin(*)
      double precision theta(*), work(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     theta is a double precision array with (nlat+1)/2 locations
c     nlat is the maximum value of n+1
c     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of work is nlat+2
c
      call wtgit1 (nlat,nlon,imid,theta,wwbin,wwbin(iw1),
     1                        work,work(nlat/2+2))
      return
      end
      subroutine wtgit1 (nlat,nlon,imid,theta,wb,abc,cwb,work)
c
c     abc must have 3*((nlat-3)*nlat+2)/2 locations
c     cwb and work must each have nlat/2+1 locations
c
      dimension wb(imid,nlat,2),abc(1)
      double precision theta(*), cwb(*), work(*), wbh
      mdo = min0(3,nlat,(nlon+1)/2)
      if(mdo .lt. 2) return
      do 160 mp1=2,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dwtk(m,n,cwb,work)
      do 165 i=1,imid
      call dwtt(m,n,theta(i),cwb,wbh)
      wb(i,np1,m) = wbh
  165 continue
  160 continue
      call rabcw(nlat,nlon,abc)
      return
      end
      subroutine dvtk(m,n,cv,work)
      double precision cv(*),work(*),fn,fk,cf,srnp1
      cv(1) = 0.
      if(n .le. 0) return
      fn = n
      srnp1 = dsqrt(fn*(fn+1.))
      cf = 2.*m/srnp1
      modn = mod(n,2)
      modm = mod(m,2)
      call dnlfk(m,n,work)
      if(modn .ne. 0) go to 70
      ncv = n/2
      if(ncv .eq. 0) return
      fk = 0.
      if(modm .ne. 0) go to 60
c
c     n even m even
c
      do 55 l=1,ncv
      fk = fk+2.
      cv(l) = -fk*fk*work(l+1)/srnp1
   55 continue
      return
c
c     n even m odd
c
   60 do 65 l=1,ncv
      fk = fk+2.
      cv(l) = -fk*fk*work(l)/srnp1
   65 continue
      return
   70 ncv = (n+1)/2
      fk = -1.
      if(modm .ne. 0) go to 80
c
c     n odd m even
c
      do 75 l=1,ncv
      fk = fk+2.
      cv(l) = -fk*fk*work(l)/srnp1
   75 continue
      return
c
c     n odd m odd
c
   80 do 85 l=1,ncv
      fk = fk+2.
      cv(l) = -fk*fk*work(l)/srnp1
   85 continue
      return
      end
      subroutine dwtk(m,n,cw,work)
      double precision cw(*),work(*),fn,cf,srnp1
      cw(1) = 0.
      if(n.le.0 .or. m.le.0) return
      fn = n
      srnp1 = dsqrt(fn*(fn+1.))
      cf = 2.*m/srnp1
      modn = mod(n,2)
      modm = mod(m,2)
      call dnlfk(m,n,work)
      if(m .eq. 0) go to 50
      if(modn .ne. 0) go to 30
      l = n/2
      if(l .eq. 0) go to 50
      if(modm .ne. 0) go to 20
c
c     n even m even
c
      cw(l) = -cf*work(l+1)
   10 l = l-1
      if(l .le. 0) go to 50
      cw(l) = cw(l+1)-cf*work(l+1)
      cw(l+1) = (l+l+1)*cw(l+1)
      go to 10
c
c     n even m odd
c
   20 cw(l) = cf*work(l)
   25 l = l-1
      if(l) 50,27,26
   26 cw(l) = cw(l+1)+cf*work(l)
   27 cw(l+1) = -(l+l+1)*cw(l+1)
      go to 25
   30 if(modm .ne. 0) go to 40
      l = (n-1)/2
      if(l .eq. 0) go to 50
c
c     n odd m even
c
      cw(l) = -cf*work(l+1)
   35 l = l-1
      if(l) 50,37,36
   36 cw(l) = cw(l+1)-cf*work(l+1)
   37 cw(l+1) = (l+l+2)*cw(l+1)
      go to 35
c
c     n odd m odd
c
   40 l = (n+1)/2
      cw(l) = cf*work(l)
   45 l = l-1
      if(l) 50,47,46
   46 cw(l) = cw(l+1)+cf*work(l)
   47 cw(l+1) = -(l+l)*cw(l+1)
      go to 45
   50 return
      end
      subroutine dvtt(m,n,theta,cv,vh)
      dimension cv(1)
      double precision cv,vh,theta,cth,sth,cdt,sdt,chh
      vh = 0.
      if(n.eq.0) return
      cth = dcos(theta)
      sth = dsin(theta)
      cdt = cth*cth-sth*sth
      sdt = 2.*sth*cth
      mmod = mod(m,2)
      nmod = mod(n,2)
      if(nmod .ne. 0) go to 1
      cth = cdt
      sth = sdt
      if(mmod .ne. 0) go to 2
c
c     n even  m even
c
      ncv = n/2
      do 10 k=1,ncv
      vh = vh+cv(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   10 continue
      return
c
c     n even  m odd
c
    2 ncv = n/2
      do 15 k=1,ncv
      vh = vh+cv(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   15 continue
      return
    1 if(mmod .ne. 0) go to 3
c
c     n odd m even
c
      ncv = (n+1)/2
      do 20 k=1,ncv
      vh = vh+cv(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   20 continue
      return
c
c case m odd and n odd
c
    3 ncv = (n+1)/2
      do 25 k=1,ncv
      vh = vh+cv(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   25 continue
      return
      end
      subroutine dwtt(m,n,theta,cw,wh)
      dimension cw(1)
      double precision theta,cw,wh,cth,sth,cdt,sdt,chh
      wh = 0.
      if(n.le.0 .or. m.le.0) return
      cth = dcos(theta)
      sth = dsin(theta)
      cdt = cth*cth-sth*sth
      sdt = 2.*sth*cth
      mmod=mod(m,2)
      nmod=mod(n,2)
      if(nmod .ne. 0) go to 1
      if(mmod .ne. 0) go to 2
c
c     n even  m even
c
      ncw = n/2
      do 10 k=1,ncw
      wh = wh+cw(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   10 continue
      return
c
c     n even  m odd
c
    2 ncw = n/2
      do 8 k=1,ncw
      wh = wh+cw(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
    8 continue
      return
    1 cth = cdt
      sth = sdt
      if(mmod .ne. 0) go to 3
c
c     n odd m even
c
      ncw = (n-1)/2
      do 20 k=1,ncw
      wh = wh+cw(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   20 continue
      return
c
c case m odd and n odd
c
    3 ncw = (n+1)/2
      wh = 0.
      if(ncw.lt.2) return
      do 25 k=2,ncw
      wh = wh+cw(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
   25 continue
      return
      end
      subroutine vbgint (nlat,nlon,theta,wvbin,work)
      dimension       wvbin(1)
      double precision theta(*),work(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     theta is a double precision array with (nlat+1)/2 locations
c     nlat is the maximum value of n+1
c     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of work is nlat+2
c
      call vbgit1 (nlat,nlon,imid,theta,wvbin,wvbin(iw1),
     +                        work,work(nlat/2+2))
      return
      end
      subroutine vbgit1 (nlat,nlon,imid,theta,vb,abc,cvb,work)
c
c     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c     locations where mmax = min0(nlat,(nlon+1)/2)
c     cvb and work must each have nlat/2+1 locations
c
      dimension vb(imid,nlat,2),abc(1)
      double precision cvb(1),theta(1),vbh,work(1)
      mdo = min0(2,nlat,(nlon+1)/2)
      do 160 mp1=1,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dvbk(m,n,cvb,work)
      do 165 i=1,imid
      call dvbt(m,n,theta(i),cvb,vbh)
      vb(i,np1,mp1) = vbh
  165 continue
  160 continue
      call rabcv(nlat,nlon,abc)
      return
      end
      subroutine wbgint (nlat,nlon,theta,wwbin,work)
      dimension       wwbin(1)
      double precision work(*),theta(*)
      imid = (nlat+1)/2
      iw1 = 2*nlat*imid+1
c
c     theta is a double precision array with (nlat+1)/2 locations
c     nlat is the maximum value of n+1
c     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
c     the length of work is nlat+2
c
      call wbgit1 (nlat,nlon,imid,theta,wwbin,wwbin(iw1),
     +                        work,work(nlat/2+2))
      return
      end
      subroutine wbgit1 (nlat,nlon,imid,theta,wb,abc,cwb,work)
c
c     abc must have 3*((nlat-3)*nlat+2)/2 locations
c     cwb and work must each have nlat/2+1 locations
c
      dimension wb(imid,nlat,2),abc(1)
      double precision cwb(1),theta(1),wbh,work(1)
      mdo = min0(3,nlat,(nlon+1)/2)
      if(mdo .lt. 2) return
      do 160 mp1=2,mdo
      m = mp1-1
      do 160 np1=mp1,nlat
      n = np1-1
      call dwbk(m,n,cwb,work)
      do 165 i=1,imid
      call dwbt(m,n,theta(i),cwb,wbh)
      wb(i,np1,m) = wbh
  165 continue
  160 continue
      call rabcw(nlat,nlon,abc)
      return
      end

 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c ... file hrfft.f
c
c     this file contains a multiple fft package for spherepack3.0.
c     it includes code and documentation for performing fast fourier
c     transforms (see subroutines hrffti,hrfftf and hrfftb)
c
c **********************************************************************
c
c     subroutine hrffti(n,wsave)
c
c     subroutine hrffti initializes the array wsave which is used in
c     both hrfftf and hrfftb. the prime factorization of n together 
c     with a tabulation of the trigonometric functions are computed and
c     stored in wsave.
c
c     input parameter
c
c     n       the length of the sequence to be transformed.
c
c     output parameter
c
c     wsave   a work array which must be dimensioned at least 2*n+15.
c             the same work array can be used for both hrfftf and 
c             hrfftb as long as n remains unchanged. different wsave 
c             arrays are required for different values of n. the 
c             contents of wsave must not be changed between calls 
c             of hrfftf or hrfftb.
c
c **********************************************************************
c
c     subroutine hrfftf(m,n,r,mdimr,wsave,work)
c
c     subroutine hrfftf computes the fourier coefficients of m real
c     perodic sequences (fourier analysis); i.e. hrfftf computes the
c     real fft of m sequences each with length n. the transform is 
c     defined below at output parameter r.
c
c     input parameters
c
c     m       the number of sequences.
c
c     n       the length of all m sequences.  the method is most
c             efficient when n is a product of small primes. n may
c             change as long as different work arrays are provided
c
c     r       r(m,n) is a two dimensional real array that contains m
c             sequences each with length n.
c
c     mdimr   the first dimension of the r array as it appears
c             in the program that calls hrfftf. mdimr must be
c             greater than or equal to m.
c
c
c     wsave   a work array with at least least 2*n+15 locations
c             in the program that calls hrfftf. the wsave array must be
c             initialized by calling subroutine hrffti(n,wsave) and a
c             different wsave array must be used for each different
c             value of n. this initialization does not have to be
c             repeated so long as n remains unchanged thus subsequent
c             transforms can be obtained faster than the first.
c             the same wsave array can be used by hrfftf and hrfftb.
c
c     work    a real work array with m*n locations.
c
c
c     output parameters
c
c     r      for all j=1,...,m
c  
c             r(j,1) = the sum from i=1 to i=n of r(j,i)
c
c             if n is even set l =n/2   , if n is odd set l = (n+1)/2
c
c               then for k = 2,...,l
c
c                  r(j,2*k-2) = the sum from i = 1 to i = n of
c
c                       r(j,i)*cos((k-1)*(i-1)*2*pi/n)
c
c                  r(j,2*k-1) = the sum from i = 1 to i = n of
c
c                      -r(j,i)*sin((k-1)*(i-1)*2*pi/n)
c
c             if n is even
c
c                  r(j,n) = the sum from i = 1 to i = n of
c
c                       (-1)**(i-1)*r(j,i)
c
c      *****  note
c                  this transform is unnormalized since a call of hrfftf
c                  followed by a call of hrfftb will multiply the input
c                  sequence by n.
c
c     wsave   contains results which must not be destroyed between
c             calls of hrfftf or hrfftb.
c
c     work    a real work array with m*n locations that does
c             not have to be saved.
c
c **********************************************************************
c
c     subroutine hrfftb(m,n,r,mdimr,wsave,work)
c
c     subroutine hrfftb computes the real perodic sequence of m
c     sequences from their fourier coefficients (fourier synthesis). 
c     the transform is defined below at output parameter r.
c
c     input parameters
c
c     m       the number of sequences.
c
c     n       the length of all m sequences.  the method is most
c             efficient when n is a product of small primes. n may
c             change as long as different work arrays are provided
c
c     r       r(m,n) is a two dimensional real array that contains
c             the fourier coefficients of m sequences each with 
c             length n.
c
c     mdimr   the first dimension of the r array as it appears
c             in the program that calls hrfftb. mdimr must be
c             greater than or equal to m.
c
c     wsave   a work array which must be dimensioned at least 2*n+15.
c             in the program that calls hrfftb. the wsave array must be
c             initialized by calling subroutine hrffti(n,wsave) and a
c             different wsave array must be used for each different
c             value of n. this initialization does not have to be
c             repeated so long as n remains unchanged thus subsequent
c             transforms can be obtained faster than the first.
c             the same wsave array can be used by hrfftf and hrfftb.
c
c     work    a real work array with m*n locations.
c
c
c     output parameters
c
c     r      for all j=1,...,m
c  
c             for n even and for i = 1,...,n
c
c                  r(j,i) = r(j,1)+(-1)**(i-1)*r(j,n)
c
c                       plus the sum from k=2 to k=n/2 of
c
c                        2.*r(j,2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                       -2.*r(j,2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c
c             for n odd and for i = 1,...,n
c
c                  r(j,i) = r(j,1) plus the sum from k=2 to k=(n+1)/2 of
c
c                       2.*r(j,2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                      -2.*r(j,2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c
c      *****  note
c                  this transform is unnormalized since a call of hrfftf
c                  followed by a call of hrfftb will multiply the input
c                  sequence by n.
c
c     wsave   contains results which must not be destroyed between
c             calls of hrfftb or hrfftf.
c
c     work    a real work array with m*n locations that does not
c             have to be saved
c
c **********************************************************************
c
c
c
      subroutine hrffti (n,wsave)
      dimension       wsave(n+15)                                              
      common /hrf/ tfft
      tfft = 0.
      if (n .eq. 1) return                                                     
      call hrfti1 (n,wsave(1),wsave(n+1))
      return                                                                   
      end                                                                      
      subroutine hrfti1 (n,wa,fac)
c                                                                              
c     a multiple fft package for spherepack
c                                                                              
      dimension       wa(n)      ,fac(15)    ,ntryh(4)                         
      double precision tpi,argh,argld,arg
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/                        
      nl = n                                                                   
      nf = 0           
      j = 0            
  101 j = j+1          
      if (j-4) 102,102,103                
  102 ntry = ntryh(j)                     
      go to 104        
  103 ntry = ntry+2    
  104 nq = nl/ntry     
      nr = nl-ntry*nq                     
      if (nr) 101,105,101                 
  105 nf = nf+1        
      fac(nf+2) = ntry                    
      nl = nq          
      if (ntry .ne. 2) go to 107          
      if (nf .eq. 1) go to 107            
      do 106 i=2,nf    
         ib = nf-i+2   
         fac(ib+2) = fac(ib+1)            
  106 continue         
      fac(3) = 2       
  107 if (nl .ne. 1) go to 104            
      fac(1) = n       
      fac(2) = nf      
      tpi = 8.d0*datan(1.d0)
      argh = tpi/float(n)                 
      is = 0           
      nfm1 = nf-1      
      l1 = 1           
      if (nfm1 .eq. 0) return             
      do 110 k1=1,nfm1                    
         ip = fac(k1+2)                   
         ld = 0        
         l2 = l1*ip    
         ido = n/l2    
         ipm = ip-1    
         do 109 j=1,ipm                   
            ld = ld+l1                    
            i = is     
            argld = float(ld)*argh        
            fi = 0.    
            do 108 ii=3,ido,2             
               i = i+2                    
               fi = fi+1.                 
               arg = fi*argld             
	       wa(i-1) = dcos(arg)
	       wa(i) = dsin(arg)
  108       continue   
            is = is+ido                   
  109    continue      
         l1 = l2       
  110 continue         
      return           
      end              
      subroutine hrfftf (m,n,r,mdimr,whrfft,work)
c                      
c     a multiple fft package for spherepack
c                      
      dimension       r(mdimr,n)  ,work(1)    ,whrfft(n+15)
      common /hrf/ tfft
      if (n .eq. 1) return                
c     tstart = second(dum)
      call hrftf1 (m,n,r,mdimr,work,whrfft,whrfft(n+1))
c     tfft = tfft+second(dum)-tstart
      return           
      end              
      subroutine hrftf1 (m,n,c,mdimc,ch,wa,fac)
c                      
c     a multiple fft package for spherepack
c                      
      dimension       ch(m,n) ,c(mdimc,n)  ,wa(n)   ,fac(15)
      nf = fac(2)      
      na = 1           
      l2 = n           
      iw = n           
      do 111 k1=1,nf   
         kh = nf-k1    
         ip = fac(kh+3)                   
         l1 = l2/ip    
         ido = n/l2    
         idl1 = ido*l1                    
         iw = iw-(ip-1)*ido               
         na = 1-na     
         if (ip .ne. 4) go to 102         
         ix2 = iw+ido                     
         ix3 = ix2+ido                    
         if (na .ne. 0) go to 101         
	 call hradf4 (m,ido,l1,c,mdimc,ch,m,wa(iw),wa(ix2),wa(ix3))
         go to 110     
  101    call hradf4 (m,ido,l1,ch,m,c,mdimc,wa(iw),wa(ix2),wa(ix3))
         go to 110     
  102    if (ip .ne. 2) go to 104         
         if (na .ne. 0) go to 103         
	 call hradf2 (m,ido,l1,c,mdimc,ch,m,wa(iw))
         go to 110     
  103    call hradf2 (m,ido,l1,ch,m,c,mdimc,wa(iw))
         go to 110     
  104    if (ip .ne. 3) go to 106         
         ix2 = iw+ido                     
         if (na .ne. 0) go to 105         
	 call hradf3 (m,ido,l1,c,mdimc,ch,m,wa(iw),wa(ix2))
         go to 110     
  105    call hradf3 (m,ido,l1,ch,m,c,mdimc,wa(iw),wa(ix2))
         go to 110     
  106    if (ip .ne. 5) go to 108         
         ix2 = iw+ido                     
         ix3 = ix2+ido                    
         ix4 = ix3+ido                    
         if (na .ne. 0) go to 107         
      call hradf5(m,ido,l1,c,mdimc,ch,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110     
  107 call hradf5(m,ido,l1,ch,m,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110     
  108    if (ido .eq. 1) na = 1-na        
         if (na .ne. 0) go to 109         
	 call hradfg (m,ido,ip,l1,idl1,c,c,c,mdimc,ch,ch,m,wa(iw))
         na = 1        
         go to 110     
  109    call hradfg (m,ido,ip,l1,idl1,ch,ch,ch,m,c,c,mdimc,wa(iw))
         na = 0        
  110    l2 = l1       
  111 continue         
      if (na .eq. 1) return
      do 112 j=1,n     
      do 112 i=1,m     
	 c(i,j) = ch(i,j)
  112 continue         
      return           
      end              
      subroutine hradf4 (mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2,wa3)
c                      
c     a multiple fft package for spherepack
c                      
      dimension    cc(mdimcc,ido,l1,4)   ,ch(mdimch,ido,4,l1)     ,
     1             wa1(ido)     ,wa2(ido)     ,wa3(ido)
      hsqt2=sqrt(2.)/2.                   
      do 101 k=1,l1    
         do 1001 m=1,mp                   
         ch(m,1,1,k) = (cc(m,1,k,2)+cc(m,1,k,4))             
     1      +(cc(m,1,k,1)+cc(m,1,k,3))    
         ch(m,ido,4,k) = (cc(m,1,k,1)+cc(m,1,k,3))           
     1      -(cc(m,1,k,2)+cc(m,1,k,4))    
         ch(m,ido,2,k) = cc(m,1,k,1)-cc(m,1,k,3)             
         ch(m,1,3,k) = cc(m,1,k,4)-cc(m,1,k,2)               
 1001    continue      
  101 continue         
      if (ido-2) 107,105,102              
  102 idp2 = ido+2     
      do 104 k=1,l1    
         do 103 i=3,ido,2                 
            ic = idp2-i                   
            do 1003 m=1,mp                
            ch(m,i-1,1,k) = ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*                 
     1       cc(m,i,k,2))+(wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*  
     1       cc(m,i,k,4)))+(cc(m,i-1,k,1)+(wa2(i-2)*cc(m,i-1,k,3)+             
     1       wa2(i-1)*cc(m,i,k,3)))       
            ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+(wa2(i-2)*cc(m,i-1,k,3)+           
     1       wa2(i-1)*cc(m,i,k,3)))-((wa1(i-2)*cc(m,i-1,k,2)+                  
     1       wa1(i-1)*cc(m,i,k,2))+(wa3(i-2)*cc(m,i-1,k,4)+  
     1       wa3(i-1)*cc(m,i,k,4)))       
            ch(m,i,1,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*   
     1       cc(m,i-1,k,2))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4)))+(cc(m,i,k,1)+(wa2(i-2)*cc(m,i,k,3)-               
     1       wa2(i-1)*cc(m,i-1,k,3)))     
            ch(m,ic,4,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*  
     1       cc(m,i-1,k,2))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4)))-(cc(m,i,k,1)+(wa2(i-2)*cc(m,i,k,3)-               
     1       wa2(i-1)*cc(m,i-1,k,3)))     
            ch(m,i-1,3,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)* 
     1       cc(m,i-1,k,2))-(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4)))+(cc(m,i-1,k,1)-(wa2(i-2)*cc(m,i-1,k,3)+           
     1       wa2(i-1)*cc(m,i,k,3)))       
            ch(m,ic-1,2,k) = (cc(m,i-1,k,1)-(wa2(i-2)*cc(m,i-1,k,3)+           
     1       wa2(i-1)*cc(m,i,k,3)))-((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*           
     1       cc(m,i-1,k,2))-(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4)))              
            ch(m,i,3,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)* 
     1       cc(m,i,k,4))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*  
     1       cc(m,i,k,2)))+(cc(m,i,k,1)-(wa2(i-2)*cc(m,i,k,3)-                 
     1       wa2(i-1)*cc(m,i-1,k,3)))     
            ch(m,ic,2,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*                  
     1       cc(m,i,k,4))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*  
     1       cc(m,i,k,2)))-(cc(m,i,k,1)-(wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*        
     1       cc(m,i-1,k,3)))              
 1003       continue   
  103    continue      
  104 continue         
      if (mod(ido,2) .eq. 1) return       
  105 continue         
      do 106 k=1,l1    
         do 1006 m=1,mp                   
            ch(m,ido,1,k) = (hsqt2*(cc(m,ido,k,2)-cc(m,ido,k,4)))+              
     1       cc(m,ido,k,1)                
            ch(m,ido,3,k) = cc(m,ido,k,1)-(hsqt2*(cc(m,ido,k,2)-                
     1       cc(m,ido,k,4)))              
            ch(m,1,2,k) = (-hsqt2*(cc(m,ido,k,2)+cc(m,ido,k,4)))-               
     1       cc(m,ido,k,3)                
            ch(m,1,4,k) = (-hsqt2*(cc(m,ido,k,2)+cc(m,ido,k,4)))+               
     1       cc(m,ido,k,3)                
 1006    continue      
  106 continue         
  107 return           
      end              
      subroutine hradf2 (mp,ido,l1,cc,mdimcc,ch,mdimch,wa1)
c                      
c     a multiple fft package for spherepack
c                      
      dimension   ch(mdimch,ido,2,l1)  ,cc(mdimcc,ido,l1,2)     ,
     1                wa1(ido)            
      do 101 k=1,l1    
         do 1001 m=1,mp                   
         ch(m,1,1,k) = cc(m,1,k,1)+cc(m,1,k,2)               
         ch(m,ido,2,k) = cc(m,1,k,1)-cc(m,1,k,2)             
 1001    continue      
  101 continue         
      if (ido-2) 107,105,102              
  102 idp2 = ido+2     
      do 104 k=1,l1    
         do 103 i=3,ido,2                 
            ic = idp2-i                   
            do 1003 m=1,mp                
            ch(m,i,1,k) = cc(m,i,k,1)+(wa1(i-2)*cc(m,i,k,2)- 
     1       wa1(i-1)*cc(m,i-1,k,2))      
            ch(m,ic,2,k) = (wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*   
     1       cc(m,i-1,k,2))-cc(m,i,k,1)   
            ch(m,i-1,1,k) = cc(m,i-1,k,1)+(wa1(i-2)*cc(m,i-1,k,2)+              
     1       wa1(i-1)*cc(m,i,k,2))        
            ch(m,ic-1,2,k) = cc(m,i-1,k,1)-(wa1(i-2)*cc(m,i-1,k,2)+             
     1       wa1(i-1)*cc(m,i,k,2))        
 1003       continue   
  103    continue      
  104 continue         
      if (mod(ido,2) .eq. 1) return       
  105 do 106 k=1,l1    
         do 1006 m=1,mp                   
         ch(m,1,2,k) = -cc(m,ido,k,2)     
         ch(m,ido,1,k) = cc(m,ido,k,1)    
 1006    continue      
  106 continue         
  107 return           
      end              
      subroutine hradf3 (mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2)
c                      
c     a multiple fft package for spherepack
c                      
      dimension   ch(mdimch,ido,3,l1)  ,cc(mdimcc,ido,l1,3)     ,
     1                wa1(ido)     ,wa2(ido)                 
      arg=2.*pimach()/3.                  
      taur=cos(arg)    
      taui=sin(arg)    
      do 101 k=1,l1    
         do 1001 m=1,mp                   
         ch(m,1,1,k) = cc(m,1,k,1)+(cc(m,1,k,2)+cc(m,1,k,3)) 
         ch(m,1,3,k) = taui*(cc(m,1,k,3)-cc(m,1,k,2))        
         ch(m,ido,2,k) = cc(m,1,k,1)+taur*                   
     1      (cc(m,1,k,2)+cc(m,1,k,3))     
 1001    continue      
  101 continue         
      if (ido .eq. 1) return              
      idp2 = ido+2     
      do 103 k=1,l1    
         do 102 i=3,ido,2                 
            ic = idp2-i                   
            do 1002 m=1,mp                
            ch(m,i-1,1,k) = cc(m,i-1,k,1)+((wa1(i-2)*cc(m,i-1,k,2)+             
     1       wa1(i-1)*cc(m,i,k,2))+(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*            
     1       cc(m,i,k,3)))                
            ch(m,i,1,k) = cc(m,i,k,1)+((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*          
     1       cc(m,i-1,k,2))+(wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*  
     1       cc(m,i-1,k,3)))              
            ch(m,i-1,3,k) = (cc(m,i-1,k,1)+taur*((wa1(i-2)*  
     1       cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2))+(wa2(i-2)*  
     1       cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))+(taui*((wa1(i-2)*            
     1       cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2))-(wa2(i-2)*  
     1       cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,3))))           
            ch(m,ic-1,2,k) = (cc(m,i-1,k,1)+taur*((wa1(i-2)* 
     1       cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2))+(wa2(i-2)*  
     1       cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))-(taui*((wa1(i-2)*            
     1       cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2))-(wa2(i-2)*  
     1       cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,3))))           
            ch(m,i,3,k) = (cc(m,i,k,1)+taur*((wa1(i-2)*cc(m,i,k,2)-             
     1       wa1(i-1)*cc(m,i-1,k,2))+(wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*            
     1       cc(m,i-1,k,3))))+(taui*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*          
     1       cc(m,i,k,3))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*  
     1       cc(m,i,k,2))))               
            ch(m,ic,2,k) = (taui*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*             
     1       cc(m,i,k,3))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*  
     1       cc(m,i,k,2))))-(cc(m,i,k,1)+taur*((wa1(i-2)*cc(m,i,k,2)-           
     1       wa1(i-1)*cc(m,i-1,k,2))+(wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*            
     1       cc(m,i-1,k,3))))             
 1002       continue   
  102    continue      
  103 continue         
      return           
      end              
      subroutine hradf5 (mp,ido,l1,cc,mdimcc,ch,mdimch,
     1                   wa1,wa2,wa3,wa4)
c                      
c     a multiple fft package for spherepack
c                      
      dimension  cc(mdimcc,ido,l1,5)    ,ch(mdimch,ido,5,l1)     ,
     1           wa1(ido)     ,wa2(ido)     ,wa3(ido)     ,wa4(ido)             
      arg=2.*pimach()/5.                  
      tr11=cos(arg)    
      ti11=sin(arg)    
      tr12=cos(2.*arg)                    
      ti12=sin(2.*arg)                    
      do 101 k=1,l1    
         do 1001 m=1,mp                   
         ch(m,1,1,k) = cc(m,1,k,1)+(cc(m,1,k,5)+cc(m,1,k,2))+                   
     1    (cc(m,1,k,4)+cc(m,1,k,3))       
         ch(m,ido,2,k) = cc(m,1,k,1)+tr11*(cc(m,1,k,5)+cc(m,1,k,2))+            
     1    tr12*(cc(m,1,k,4)+cc(m,1,k,3))  
         ch(m,1,3,k) = ti11*(cc(m,1,k,5)-cc(m,1,k,2))+ti12*  
     1    (cc(m,1,k,4)-cc(m,1,k,3))       
         ch(m,ido,4,k) = cc(m,1,k,1)+tr12*(cc(m,1,k,5)+cc(m,1,k,2))+            
     1    tr11*(cc(m,1,k,4)+cc(m,1,k,3))  
         ch(m,1,5,k) = ti12*(cc(m,1,k,5)-cc(m,1,k,2))-ti11*  
     1    (cc(m,1,k,4)-cc(m,1,k,3))       
 1001    continue      
  101 continue         
      if (ido .eq. 1) return              
      idp2 = ido+2     
      do 103 k=1,l1    
         do 102 i=3,ido,2                 
            ic = idp2-i                   
            do 1002 m=1,mp                
            ch(m,i-1,1,k) = cc(m,i-1,k,1)+((wa1(i-2)*cc(m,i-1,k,2)+             
     1       wa1(i-1)*cc(m,i,k,2))+(wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*            
     1       cc(m,i,k,5)))+((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*                   
     1       cc(m,i,k,3))+(wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4)))        
            ch(m,i,1,k) = cc(m,i,k,1)+((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*          
     1       cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*  
     1       cc(m,i-1,k,5)))+((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*                   
     1       cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4)))              
            ch(m,i-1,3,k) = cc(m,i-1,k,1)+tr11*              
     1      ( wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2)    
     1       +wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i,k,5))+tr12*                
     1      ( wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3)    
     1       +wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4))+ti11*                
     1      ( wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2)    
     1       -(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ti12*              
     1      ( wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,3)    
     1       -(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,k,4))) 
            ch(m,ic-1,2,k) = cc(m,i-1,k,1)+tr11*             
     1      ( wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2)    
     1       +wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i,k,5))+tr12*                
     1     ( wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3)     
     1      +wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4))-(ti11*                
     1      ( wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2)    
     1       -(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ti12*              
     1      ( wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,3)    
     1       -(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,k,4))))                   
            ch(m,i,3,k) = (cc(m,i,k,1)+tr11*((wa1(i-2)*cc(m,i,k,2)-             
     1       wa1(i-1)*cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*            
     1       cc(m,i-1,k,5)))+tr12*((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*              
     1       cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4))))+(ti11*((wa4(i-2)*cc(m,i-1,k,5)+                   
     1       wa4(i-1)*cc(m,i,k,5))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*            
     1       cc(m,i,k,2)))+ti12*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*              
     1       cc(m,i,k,4))-(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*  
     1       cc(m,i,k,3))))               
            ch(m,ic,2,k) = (ti11*((wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*             
     1       cc(m,i,k,5))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*  
     1       cc(m,i,k,2)))+ti12*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*              
     1       cc(m,i,k,4))-(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*  
     1       cc(m,i,k,3))))-(cc(m,i,k,1)+tr11*((wa1(i-2)*cc(m,i,k,2)-           
     1       wa1(i-1)*cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*            
     1       cc(m,i-1,k,5)))+tr12*((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*              
     1       cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4))))             
            ch(m,i-1,5,k) = (cc(m,i-1,k,1)+tr12*((wa1(i-2)*  
     1       cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2))+(wa4(i-2)*  
     1       cc(m,i-1,k,5)+wa4(i-1)*cc(m,i,k,5)))+tr11*((wa2(i-2)*              
     1       cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3))+(wa3(i-2)*  
     1       cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4))))+(ti12*((wa1(i-2)*            
     1       cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2))-(wa4(i-2)*cc(m,i,k,5)-         
     1       wa4(i-1)*cc(m,i-1,k,5)))-ti11*((wa2(i-2)*cc(m,i,k,3)-              
     1       wa2(i-1)*cc(m,i-1,k,3))-(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*            
     1       cc(m,i-1,k,4))))             
            ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+tr12*((wa1(i-2)* 
     1       cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2))+(wa4(i-2)*  
     1       cc(m,i-1,k,5)+wa4(i-1)*cc(m,i,k,5)))+tr11*((wa2(i-2)*              
     1       cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3))+(wa3(i-2)*  
     1       cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4))))-(ti12*((wa1(i-2)*            
     1       cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2))-(wa4(i-2)*cc(m,i,k,5)-         
     1       wa4(i-1)*cc(m,i-1,k,5)))-ti11*((wa2(i-2)*cc(m,i,k,3)-              
     1       wa2(i-1)*cc(m,i-1,k,3))-(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*            
     1       cc(m,i-1,k,4))))             
            ch(m,i,5,k) = (cc(m,i,k,1)+tr12*((wa1(i-2)*cc(m,i,k,2)-             
     1       wa1(i-1)*cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*            
     1       cc(m,i-1,k,5)))+tr11*((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*              
     1       cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4))))+(ti12*((wa4(i-2)*cc(m,i-1,k,5)+                   
     1       wa4(i-1)*cc(m,i,k,5))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*            
     1       cc(m,i,k,2)))-ti11*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*              
     1       cc(m,i,k,4))-(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*  
     1       cc(m,i,k,3))))               
            ch(m,ic,4,k) = (ti12*((wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*             
     1       cc(m,i,k,5))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*  
     1       cc(m,i,k,2)))-ti11*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*              
     1       cc(m,i,k,4))-(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*  
     1       cc(m,i,k,3))))-(cc(m,i,k,1)+tr12*((wa1(i-2)*cc(m,i,k,2)-           
     1       wa1(i-1)*cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*            
     1       cc(m,i-1,k,5)))+tr11*((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*              
     1       cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*  
     1       cc(m,i-1,k,4))))             
 1002       continue   
  102    continue      
  103 continue         
      return           
      end              
      subroutine hradfg (mp,ido,ip,l1,idl1,cc,c1,c2,mdimcc,
     1              ch,ch2,mdimch,wa)
c                      
c     a multiple fft package for spherepack
c                      
      dimension     ch(mdimch,ido,l1,ip)   ,cc(mdimcc,ido,ip,l1)  ,
     1            c1(mdimcc,ido,l1,ip)    ,c2(mdimcc,idl1,ip),
     2                ch2(mdimch,idl1,ip)           ,wa(ido)
      tpi=2.*pimach()                     
      arg = tpi/float(ip)                 
      dcp = cos(arg)   
      dsp = sin(arg)   
      ipph = (ip+1)/2                     
      ipp2 = ip+2      
      idp2 = ido+2     
      nbd = (ido-1)/2                     
      if (ido .eq. 1) go to 119           
      do 101 ik=1,idl1                    
         do 1001 m=1,mp                   
         ch2(m,ik,1) = c2(m,ik,1)         
 1001    continue      
  101 continue         
      do 103 j=2,ip    
         do 102 k=1,l1                    
            do 1002 m=1,mp                
            ch(m,1,k,j) = c1(m,1,k,j)     
 1002       continue   
  102    continue      
  103 continue         
      if (nbd .gt. l1) go to 107          
      is = -ido        
      do 106 j=2,ip    
         is = is+ido   
         idij = is     
         do 105 i=3,ido,2                 
            idij = idij+2                 
            do 104 k=1,l1                 
               do 1004 m=1,mp             
               ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j)+wa(idij)                
     1           *c1(m,i,k,j)             
               ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j)-wa(idij) 
     1           *c1(m,i-1,k,j)           
 1004          continue                   
  104       continue   
  105    continue      
  106 continue         
      go to 111        
  107 is = -ido        
      do 110 j=2,ip    
         is = is+ido   
         do 109 k=1,l1                    
            idij = is                     
            do 108 i=3,ido,2              
               idij = idij+2              
               do 1008 m=1,mp             
               ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j)+wa(idij)                
     1           *c1(m,i,k,j)             
               ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j)-wa(idij) 
     1           *c1(m,i-1,k,j)           
 1008          continue                   
  108       continue   
  109    continue      
  110 continue         
  111 if (nbd .lt. l1) go to 115          
      do 114 j=2,ipph                     
         jc = ipp2-j   
         do 113 k=1,l1                    
            do 112 i=3,ido,2              
               do 1012 m=1,mp             
               c1(m,i-1,k,j) = ch(m,i-1,k,j)+ch(m,i-1,k,jc)  
               c1(m,i-1,k,jc) = ch(m,i,k,j)-ch(m,i,k,jc)     
               c1(m,i,k,j) = ch(m,i,k,j)+ch(m,i,k,jc)        
               c1(m,i,k,jc) = ch(m,i-1,k,jc)-ch(m,i-1,k,j)   
 1012          continue                   
  112       continue   
  113    continue      
  114 continue         
      go to 121        
  115 do 118 j=2,ipph                     
         jc = ipp2-j   
         do 117 i=3,ido,2                 
            do 116 k=1,l1                 
               do 1016 m=1,mp             
               c1(m,i-1,k,j) = ch(m,i-1,k,j)+ch(m,i-1,k,jc)  
               c1(m,i-1,k,jc) = ch(m,i,k,j)-ch(m,i,k,jc)     
               c1(m,i,k,j) = ch(m,i,k,j)+ch(m,i,k,jc)        
               c1(m,i,k,jc) = ch(m,i-1,k,jc)-ch(m,i-1,k,j)   
 1016          continue                   
  116       continue   
  117    continue      
  118 continue         
      go to 121        
  119 do 120 ik=1,idl1                    
         do 1020 m=1,mp                   
         c2(m,ik,1) = ch2(m,ik,1)         
 1020    continue      
  120 continue         
  121 do 123 j=2,ipph                     
         jc = ipp2-j   
         do 122 k=1,l1                    
            do 1022 m=1,mp                
            c1(m,1,k,j) = ch(m,1,k,j)+ch(m,1,k,jc)           
            c1(m,1,k,jc) = ch(m,1,k,jc)-ch(m,1,k,j)          
 1022       continue   
  122    continue      
  123 continue         
c                      
      ar1 = 1.         
      ai1 = 0.         
      do 127 l=2,ipph                     
         lc = ipp2-l   
         ar1h = dcp*ar1-dsp*ai1           
         ai1 = dcp*ai1+dsp*ar1            
         ar1 = ar1h    
         do 124 ik=1,idl1                 
            do 1024 m=1,mp                
            ch2(m,ik,l) = c2(m,ik,1)+ar1*c2(m,ik,2)          
            ch2(m,ik,lc) = ai1*c2(m,ik,ip)                   
 1024       continue   
  124    continue      
         dc2 = ar1     
         ds2 = ai1     
         ar2 = ar1     
         ai2 = ai1     
         do 126 j=3,ipph                  
            jc = ipp2-j                   
            ar2h = dc2*ar2-ds2*ai2        
            ai2 = dc2*ai2+ds2*ar2         
            ar2 = ar2h                    
            do 125 ik=1,idl1              
               do 1025 m=1,mp             
               ch2(m,ik,l) = ch2(m,ik,l)+ar2*c2(m,ik,j)      
               ch2(m,ik,lc) = ch2(m,ik,lc)+ai2*c2(m,ik,jc)   
 1025          continue                   
  125       continue   
  126    continue      
  127 continue         
      do 129 j=2,ipph                     
         do 128 ik=1,idl1                 
            do 1028 m=1,mp                
            ch2(m,ik,1) = ch2(m,ik,1)+c2(m,ik,j)             
 1028       continue   
  128    continue      
  129 continue         
c                      
      if (ido .lt. l1) go to 132          
      do 131 k=1,l1    
         do 130 i=1,ido                   
            do 1030 m=1,mp                
            cc(m,i,1,k) = ch(m,i,k,1)     
 1030       continue   
  130    continue      
  131 continue         
      go to 135        
  132 do 134 i=1,ido   
         do 133 k=1,l1                    
            do 1033 m=1,mp                
            cc(m,i,1,k) = ch(m,i,k,1)     
 1033       continue   
  133    continue      
  134 continue         
  135 do 137 j=2,ipph                     
         jc = ipp2-j   
         j2 = j+j      
         do 136 k=1,l1                    
            do 1036 m=1,mp                
            cc(m,ido,j2-2,k) = ch(m,1,k,j)                   
            cc(m,1,j2-1,k) = ch(m,1,k,jc)                    
 1036       continue   
  136    continue      
  137 continue         
      if (ido .eq. 1) return              
      if (nbd .lt. l1) go to 141          
      do 140 j=2,ipph                     
         jc = ipp2-j   
         j2 = j+j      
         do 139 k=1,l1                    
            do 138 i=3,ido,2              
               ic = idp2-i                
               do 1038 m=1,mp             
               cc(m,i-1,j2-1,k) = ch(m,i-1,k,j)+ch(m,i-1,k,jc)                  
               cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j)-ch(m,i-1,k,jc)                 
               cc(m,i,j2-1,k) = ch(m,i,k,j)+ch(m,i,k,jc)     
               cc(m,ic,j2-2,k) = ch(m,i,k,jc)-ch(m,i,k,j)    
 1038          continue                   
  138       continue   
  139    continue      
  140 continue         
      return           
  141 do 144 j=2,ipph                     
         jc = ipp2-j   
         j2 = j+j      
         do 143 i=3,ido,2                 
            ic = idp2-i                   
            do 142 k=1,l1                 
               do 1042 m=1,mp             
               cc(m,i-1,j2-1,k) = ch(m,i-1,k,j)+ch(m,i-1,k,jc)                  
               cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j)-ch(m,i-1,k,jc)                 
               cc(m,i,j2-1,k) = ch(m,i,k,j)+ch(m,i,k,jc)     
               cc(m,ic,j2-2,k) = ch(m,i,k,jc)-ch(m,i,k,j)    
 1042          continue                   
  142       continue   
  143    continue      
  144 continue         
      return           
      end              
c      function pimach()                   
c      pimach=3.14159265358979             
c      return           
c      end              
      subroutine hrfftb(m,n,r,mdimr,whrfft,work)
c                      
c     a multiple fft package for spherepack
c                      
      dimension     r(mdimr,n)  ,work(1)    ,whrfft(n+15)
      common /hrf/ tfft
      if (n .eq. 1) return                
c     tstart = second(dum)
      call hrftb1 (m,n,r,mdimr,work,whrfft,whrfft(n+1))
c     tfft = tfft+second(dum)-tstart
      return           
      end              
      subroutine hrftb1 (m,n,c,mdimc,ch,wa,fac)
c                      
c     a multiple fft package for spherepack
c                      
      dimension       ch(m,n), c(mdimc,n), wa(n) ,fac(15)
      nf = fac(2)      
      na = 0           
      l1 = 1           
      iw = 1           
      do 116 k1=1,nf   
         ip = fac(k1+2)                   
         l2 = ip*l1    
         ido = n/l2    
         idl1 = ido*l1                    
         if (ip .ne. 4) go to 103         
         ix2 = iw+ido                     
         ix3 = ix2+ido                    
         if (na .ne. 0) go to 101         
	 call hradb4 (m,ido,l1,c,mdimc,ch,m,wa(iw),wa(ix2),wa(ix3))
         go to 102     
  101    call hradb4 (m,ido,l1,ch,m,c,mdimc,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na     
         go to 115     
  103    if (ip .ne. 2) go to 106         
         if (na .ne. 0) go to 104         
	 call hradb2 (m,ido,l1,c,mdimc,ch,m,wa(iw))
         go to 105     
  104    call hradb2 (m,ido,l1,ch,m,c,mdimc,wa(iw))
  105    na = 1-na     
         go to 115     
  106    if (ip .ne. 3) go to 109         
         ix2 = iw+ido                     
         if (na .ne. 0) go to 107         
	 call hradb3 (m,ido,l1,c,mdimc,ch,m,wa(iw),wa(ix2))
         go to 108     
  107    call hradb3 (m,ido,l1,ch,m,c,mdimc,wa(iw),wa(ix2))
  108    na = 1-na     
         go to 115     
  109    if (ip .ne. 5) go to 112         
         ix2 = iw+ido                     
         ix3 = ix2+ido                    
         ix4 = ix3+ido                    
         if (na .ne. 0) go to 110         
      call hradb5 (m,ido,l1,c,mdimc,ch,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111     
  110 call hradb5 (m,ido,l1,ch,m,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na     
         go to 115     
  112    if (na .ne. 0) go to 113         
	 call hradbg (m,ido,ip,l1,idl1,c,c,c,mdimc,ch,ch,m,wa(iw))
         go to 114     
  113    call hradbg (m,ido,ip,l1,idl1,ch,ch,ch,m,c,c,mdimc,wa(iw))
  114    if (ido .eq. 1) na = 1-na        
  115    l1 = l2       
         iw = iw+(ip-1)*ido               
  116 continue         
      if (na .eq. 0) return
      do 117 j=1,n     
      do 117 i=1,m     
	 c(i,j) = ch(i,j)
  117 continue         
      return           
      end              
      subroutine hradbg (mp,ido,ip,l1,idl1,cc,c1,c2,mdimcc,
     1          ch,ch2,mdimch,wa)
c                      
c     a multiple fft package for spherepack
c                      
      dimension    ch(mdimch,ido,l1,ip)    ,cc(mdimcc,ido,ip,l1) ,
     1           c1(mdimcc,ido,l1,ip)     ,c2(mdimcc,idl1,ip),
     2                ch2(mdimch,idl1,ip)       ,wa(ido)
      tpi=2.*pimach()                     
      arg = tpi/float(ip)                 
      dcp = cos(arg)   
      dsp = sin(arg)   
      idp2 = ido+2     
      nbd = (ido-1)/2                     
      ipp2 = ip+2      
      ipph = (ip+1)/2                     
      if (ido .lt. l1) go to 103          
      do 102 k=1,l1    
         do 101 i=1,ido                   
            do 1001 m=1,mp                
            ch(m,i,k,1) = cc(m,i,1,k)     
 1001       continue   
  101    continue      
  102 continue         
      go to 106        
  103 do 105 i=1,ido   
         do 104 k=1,l1                    
            do 1004 m=1,mp                
            ch(m,i,k,1) = cc(m,i,1,k)     
 1004       continue   
  104    continue      
  105 continue         
  106 do 108 j=2,ipph                     
         jc = ipp2-j   
         j2 = j+j      
         do 107 k=1,l1                    
            do 1007 m=1,mp                
            ch(m,1,k,j) = cc(m,ido,j2-2,k)+cc(m,ido,j2-2,k)  
            ch(m,1,k,jc) = cc(m,1,j2-1,k)+cc(m,1,j2-1,k)     
 1007       continue   
  107    continue      
  108 continue         
      if (ido .eq. 1) go to 116           
      if (nbd .lt. l1) go to 112          
      do 111 j=2,ipph                     
         jc = ipp2-j   
         do 110 k=1,l1                    
            do 109 i=3,ido,2              
               ic = idp2-i                
               do 1009 m=1,mp             
               ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k)+cc(m,ic-1,2*j-2,k)             
               ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k)-cc(m,ic-1,2*j-2,k)            
               ch(m,i,k,j) = cc(m,i,2*j-1,k)-cc(m,ic,2*j-2,k)                   
               ch(m,i,k,jc) = cc(m,i,2*j-1,k)+cc(m,ic,2*j-2,k)                  
 1009          continue                   
  109       continue   
  110    continue      
  111 continue         
      go to 116        
  112 do 115 j=2,ipph                     
         jc = ipp2-j   
         do 114 i=3,ido,2                 
            ic = idp2-i                   
            do 113 k=1,l1                 
               do 1013 m=1,mp             
               ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k)+cc(m,ic-1,2*j-2,k)             
               ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k)-cc(m,ic-1,2*j-2,k)            
               ch(m,i,k,j) = cc(m,i,2*j-1,k)-cc(m,ic,2*j-2,k)                   
               ch(m,i,k,jc) = cc(m,i,2*j-1,k)+cc(m,ic,2*j-2,k)                  
 1013          continue                   
  113       continue   
  114    continue      
  115 continue         
  116 ar1 = 1.         
      ai1 = 0.         
      do 120 l=2,ipph                     
         lc = ipp2-l   
         ar1h = dcp*ar1-dsp*ai1           
         ai1 = dcp*ai1+dsp*ar1            
         ar1 = ar1h    
         do 117 ik=1,idl1                 
            do 1017 m=1,mp                
            c2(m,ik,l) = ch2(m,ik,1)+ar1*ch2(m,ik,2)         
            c2(m,ik,lc) = ai1*ch2(m,ik,ip)                   
 1017       continue   
  117    continue      
         dc2 = ar1     
         ds2 = ai1     
         ar2 = ar1     
         ai2 = ai1     
         do 119 j=3,ipph                  
            jc = ipp2-j                   
            ar2h = dc2*ar2-ds2*ai2        
            ai2 = dc2*ai2+ds2*ar2         
            ar2 = ar2h                    
            do 118 ik=1,idl1              
               do 1018 m=1,mp             
               c2(m,ik,l) = c2(m,ik,l)+ar2*ch2(m,ik,j)       
               c2(m,ik,lc) = c2(m,ik,lc)+ai2*ch2(m,ik,jc)    
 1018          continue                   
  118       continue   
  119    continue      
  120 continue         
      do 122 j=2,ipph                     
         do 121 ik=1,idl1                 
            do 1021 m=1,mp                
            ch2(m,ik,1) = ch2(m,ik,1)+ch2(m,ik,j)            
 1021       continue   
  121    continue      
  122 continue         
      do 124 j=2,ipph                     
         jc = ipp2-j   
         do 123 k=1,l1                    
            do 1023 m=1,mp                
            ch(m,1,k,j) = c1(m,1,k,j)-c1(m,1,k,jc)           
            ch(m,1,k,jc) = c1(m,1,k,j)+c1(m,1,k,jc)          
 1023       continue   
  123    continue      
  124 continue         
      if (ido .eq. 1) go to 132           
      if (nbd .lt. l1) go to 128          
      do 127 j=2,ipph                     
         jc = ipp2-j   
         do 126 k=1,l1                    
            do 125 i=3,ido,2              
               do 1025 m=1,mp             
               ch(m,i-1,k,j) = c1(m,i-1,k,j)-c1(m,i,k,jc)    
               ch(m,i-1,k,jc) = c1(m,i-1,k,j)+c1(m,i,k,jc)   
               ch(m,i,k,j) = c1(m,i,k,j)+c1(m,i-1,k,jc)      
               ch(m,i,k,jc) = c1(m,i,k,j)-c1(m,i-1,k,jc)     
 1025          continue                   
  125       continue   
  126    continue      
  127 continue         
      go to 132        
  128 do 131 j=2,ipph                     
         jc = ipp2-j   
         do 130 i=3,ido,2                 
            do 129 k=1,l1                 
               do 1029 m=1,mp             
               ch(m,i-1,k,j) = c1(m,i-1,k,j)-c1(m,i,k,jc)    
               ch(m,i-1,k,jc) = c1(m,i-1,k,j)+c1(m,i,k,jc)   
               ch(m,i,k,j) = c1(m,i,k,j)+c1(m,i-1,k,jc)      
               ch(m,i,k,jc) = c1(m,i,k,j)-c1(m,i-1,k,jc)     
 1029          continue                   
  129       continue   
  130    continue      
  131 continue         
  132 continue         
      if (ido .eq. 1) return              
      do 133 ik=1,idl1                    
         do 1033 m=1,mp                   
         c2(m,ik,1) = ch2(m,ik,1)         
 1033    continue      
  133 continue         
      do 135 j=2,ip    
         do 134 k=1,l1                    
            do 1034 m=1,mp                
            c1(m,1,k,j) = ch(m,1,k,j)     
 1034       continue   
  134    continue      
  135 continue         
      if (nbd .gt. l1) go to 139          
      is = -ido        
      do 138 j=2,ip    
         is = is+ido   
         idij = is     
         do 137 i=3,ido,2                 
            idij = idij+2                 
            do 136 k=1,l1                 
               do 1036 m=1,mp             
               c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j)-wa(idij)*               
     1          ch(m,i,k,j)               
               c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j)+wa(idij)*                   
     1          ch(m,i-1,k,j)             
 1036          continue                   
  136       continue   
  137    continue      
  138 continue         
      go to 143        
  139 is = -ido        
      do 142 j=2,ip    
         is = is+ido   
         do 141 k=1,l1                    
            idij = is                     
            do 140 i=3,ido,2              
               idij = idij+2              
               do 1040 m=1,mp             
               c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j)-wa(idij)*               
     1          ch(m,i,k,j)               
               c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j)+wa(idij)*                   
     1          ch(m,i-1,k,j)             
 1040          continue                   
  140       continue   
  141    continue      
  142 continue         
  143 return           
      end              
      subroutine hradb4 (mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2,wa3)
c                      
c     a multiple fft package for spherepack
c                      
      dimension  cc(mdimcc,ido,4,l1)  ,ch(mdimch,ido,l1,4)    ,
     1                wa1(ido)  ,wa2(ido)  ,wa3(ido)         
      sqrt2=sqrt(2.)   
      do 101 k=1,l1    
          do 1001 m=1,mp                  
         ch(m,1,k,3) = (cc(m,1,1,k)+cc(m,ido,4,k))           
     1   -(cc(m,ido,2,k)+cc(m,ido,2,k))   
         ch(m,1,k,1) = (cc(m,1,1,k)+cc(m,ido,4,k))           
     1   +(cc(m,ido,2,k)+cc(m,ido,2,k))   
         ch(m,1,k,4) = (cc(m,1,1,k)-cc(m,ido,4,k))           
     1   +(cc(m,1,3,k)+cc(m,1,3,k))       
         ch(m,1,k,2) = (cc(m,1,1,k)-cc(m,ido,4,k))           
     1   -(cc(m,1,3,k)+cc(m,1,3,k))       
 1001     continue     
  101 continue         
      if (ido-2) 107,105,102              
  102 idp2 = ido+2     
      do 104 k=1,l1    
         do 103 i=3,ido,2                 
            ic = idp2-i                   
               do 1002 m=1,mp             
            ch(m,i-1,k,1) = (cc(m,i-1,1,k)+cc(m,ic-1,4,k))   
     1      +(cc(m,i-1,3,k)+cc(m,ic-1,2,k))                  
            ch(m,i,k,1) = (cc(m,i,1,k)-cc(m,ic,4,k))         
     1      +(cc(m,i,3,k)-cc(m,ic,2,k))   
            ch(m,i-1,k,2)=wa1(i-2)*((cc(m,i-1,1,k)-cc(m,ic-1,4,k))              
     1      -(cc(m,i,3,k)+cc(m,ic,2,k)))-wa1(i-1)            
     1      *((cc(m,i,1,k)+cc(m,ic,4,k))+(cc(m,i-1,3,k)-cc(m,ic-1,2,k)))        
            ch(m,i,k,2)=wa1(i-2)*((cc(m,i,1,k)+cc(m,ic,4,k)) 
     1      +(cc(m,i-1,3,k)-cc(m,ic-1,2,k)))+wa1(i-1)        
     1      *((cc(m,i-1,1,k)-cc(m,ic-1,4,k))-(cc(m,i,3,k)+cc(m,ic,2,k)))        
            ch(m,i-1,k,3)=wa2(i-2)*((cc(m,i-1,1,k)+cc(m,ic-1,4,k))              
     1      -(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))-wa2(i-1)        
     1      *((cc(m,i,1,k)-cc(m,ic,4,k))-(cc(m,i,3,k)-cc(m,ic,2,k)))            
            ch(m,i,k,3)=wa2(i-2)*((cc(m,i,1,k)-cc(m,ic,4,k)) 
     1      -(cc(m,i,3,k)-cc(m,ic,2,k)))+wa2(i-1)            
     1      *((cc(m,i-1,1,k)+cc(m,ic-1,4,k))-(cc(m,i-1,3,k)  
     1      +cc(m,ic-1,2,k)))             
            ch(m,i-1,k,4)=wa3(i-2)*((cc(m,i-1,1,k)-cc(m,ic-1,4,k))              
     1      +(cc(m,i,3,k)+cc(m,ic,2,k)))-wa3(i-1)            
     1     *((cc(m,i,1,k)+cc(m,ic,4,k))-(cc(m,i-1,3,k)-cc(m,ic-1,2,k)))         
            ch(m,i,k,4)=wa3(i-2)*((cc(m,i,1,k)+cc(m,ic,4,k)) 
     1      -(cc(m,i-1,3,k)-cc(m,ic-1,2,k)))+wa3(i-1)        
     1      *((cc(m,i-1,1,k)-cc(m,ic-1,4,k))+(cc(m,i,3,k)+cc(m,ic,2,k)))        
 1002          continue                   
  103    continue      
  104 continue         
      if (mod(ido,2) .eq. 1) return       
  105 continue         
      do 106 k=1,l1    
               do 1003 m=1,mp             
         ch(m,ido,k,1) = (cc(m,ido,1,k)+cc(m,ido,3,k))       
     1   +(cc(m,ido,1,k)+cc(m,ido,3,k))   
         ch(m,ido,k,2) = sqrt2*((cc(m,ido,1,k)-cc(m,ido,3,k))                   
     1   -(cc(m,1,2,k)+cc(m,1,4,k)))      
         ch(m,ido,k,3) = (cc(m,1,4,k)-cc(m,1,2,k))           
     1   +(cc(m,1,4,k)-cc(m,1,2,k))       
         ch(m,ido,k,4) = -sqrt2*((cc(m,ido,1,k)-cc(m,ido,3,k))                  
     1   +(cc(m,1,2,k)+cc(m,1,4,k)))      
 1003          continue                   
  106 continue         
  107 return           
      end              
      subroutine hradb2 (mp,ido,l1,cc,mdimcc,ch,mdimch,wa1)
c                      
c     a multiple fft package for spherepack
c                      
      dimension  cc(mdimcc,ido,2,l1)    ,ch(mdimch,ido,l1,2),
     1                wa1(ido)            
      do 101 k=1,l1    
          do 1001 m=1,mp                  
         ch(m,1,k,1) = cc(m,1,1,k)+cc(m,ido,2,k)             
         ch(m,1,k,2) = cc(m,1,1,k)-cc(m,ido,2,k)             
 1001     continue     
  101 continue         
      if (ido-2) 107,105,102              
  102 idp2 = ido+2     
      do 104 k=1,l1    
         do 103 i=3,ido,2                 
            ic = idp2-i                   
               do 1002 m=1,mp             
            ch(m,i-1,k,1) = cc(m,i-1,1,k)+cc(m,ic-1,2,k)     
            ch(m,i,k,1) = cc(m,i,1,k)-cc(m,ic,2,k)           
            ch(m,i-1,k,2) = wa1(i-2)*(cc(m,i-1,1,k)-cc(m,ic-1,2,k))             
     1      -wa1(i-1)*(cc(m,i,1,k)+cc(m,ic,2,k))             
            ch(m,i,k,2) = wa1(i-2)*(cc(m,i,1,k)+cc(m,ic,2,k))+wa1(i-1)          
     1      *(cc(m,i-1,1,k)-cc(m,ic-1,2,k))                  
 1002          continue                   
  103    continue      
  104 continue         
      if (mod(ido,2) .eq. 1) return       
  105 do 106 k=1,l1    
          do 1003 m=1,mp                  
         ch(m,ido,k,1) = cc(m,ido,1,k)+cc(m,ido,1,k)         
         ch(m,ido,k,2) = -(cc(m,1,2,k)+cc(m,1,2,k))          
 1003     continue     
  106 continue         
  107 return           
      end              
      subroutine hradb3 (mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2)
c                      
c     a multiple fft package for spherepack
c                      
      dimension  cc(mdimcc,ido,3,l1)    ,ch(mdimch,ido,l1,3),
     1                wa1(ido)   ,wa2(ido)                   
      arg=2.*pimach()/3.                  
      taur=cos(arg)    
      taui=sin(arg)    
      do 101 k=1,l1    
          do 1001 m=1,mp                  
         ch(m,1,k,1) = cc(m,1,1,k)+2.*cc(m,ido,2,k)          
         ch(m,1,k,2) = cc(m,1,1,k)+(2.*taur)*cc(m,ido,2,k)   
     1   -(2.*taui)*cc(m,1,3,k)           
         ch(m,1,k,3) = cc(m,1,1,k)+(2.*taur)*cc(m,ido,2,k)   
     1   +2.*taui*cc(m,1,3,k)             
 1001     continue     
  101 continue         
      if (ido .eq. 1) return              
      idp2 = ido+2     
      do 103 k=1,l1    
         do 102 i=3,ido,2                 
            ic = idp2-i                   
               do 1002 m=1,mp             
            ch(m,i-1,k,1) = cc(m,i-1,1,k)+(cc(m,i-1,3,k)+cc(m,ic-1,2,k))        
            ch(m,i,k,1) = cc(m,i,1,k)+(cc(m,i,3,k)-cc(m,ic,2,k))                
            ch(m,i-1,k,2) = wa1(i-2)*     
     1 ((cc(m,i-1,1,k)+taur*(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))- 
     * (taui*(cc(m,i,3,k)+cc(m,ic,2,k))))                    
     2                   -wa1(i-1)*       
     3 ((cc(m,i,1,k)+taur*(cc(m,i,3,k)-cc(m,ic,2,k)))+       
     * (taui*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))))                
            ch(m,i,k,2) = wa1(i-2)*       
     4 ((cc(m,i,1,k)+taur*(cc(m,i,3,k)-cc(m,ic,2,k)))+       
     8 (taui*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))))                
     5                  +wa1(i-1)*        
     6 ((cc(m,i-1,1,k)+taur*(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))- 
     8 (taui*(cc(m,i,3,k)+cc(m,ic,2,k))))                    
              ch(m,i-1,k,3) = wa2(i-2)*   
     7 ((cc(m,i-1,1,k)+taur*(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))+ 
     8 (taui*(cc(m,i,3,k)+cc(m,ic,2,k))))                    
     8   -wa2(i-1)*    
     9 ((cc(m,i,1,k)+taur*(cc(m,i,3,k)-cc(m,ic,2,k)))-       
     8 (taui*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))))                
            ch(m,i,k,3) = wa2(i-2)*       
     1 ((cc(m,i,1,k)+taur*(cc(m,i,3,k)-cc(m,ic,2,k)))-       
     8 (taui*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))))                
     2                 +wa2(i-1)*         
     3 ((cc(m,i-1,1,k)+taur*(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))+ 
     8 (taui*(cc(m,i,3,k)+cc(m,ic,2,k))))                    
 1002          continue                   
  102    continue      
  103 continue         
      return           
      end              
      subroutine hradb5 (mp,ido,l1,cc,mdimcc,ch,mdimch,
     1       wa1,wa2,wa3,wa4)
c                      
c     a multiple fft package for spherepack
c                      
      dimension  cc(mdimcc,ido,5,l1)    ,ch(mdimch,ido,l1,5),
     1             wa1(ido)     ,wa2(ido)     ,wa3(ido)     ,wa4(ido)           
      arg=2.*pimach()/5.                  
      tr11=cos(arg)    
      ti11=sin(arg)    
      tr12=cos(2.*arg)                    
      ti12=sin(2.*arg)                    
      do 101 k=1,l1    
      do 1001 m=1,mp   
         ch(m,1,k,1) = cc(m,1,1,k)+2.*cc(m,ido,2,k)+2.*cc(m,ido,4,k)            
         ch(m,1,k,2) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)    
     1   +tr12*2.*cc(m,ido,4,k))-(ti11*2.*cc(m,1,3,k)        
     1   +ti12*2.*cc(m,1,5,k))            
         ch(m,1,k,3) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)    
     1   +tr11*2.*cc(m,ido,4,k))-(ti12*2.*cc(m,1,3,k)        
     1   -ti11*2.*cc(m,1,5,k))            
         ch(m,1,k,4) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)    
     1   +tr11*2.*cc(m,ido,4,k))+(ti12*2.*cc(m,1,3,k)        
     1   -ti11*2.*cc(m,1,5,k))            
         ch(m,1,k,5) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)    
     1   +tr12*2.*cc(m,ido,4,k))+(ti11*2.*cc(m,1,3,k)        
     1   +ti12*2.*cc(m,1,5,k))            
 1001          continue                   
  101 continue         
      if (ido .eq. 1) return              
      idp2 = ido+2     
      do 103 k=1,l1    
         do 102 i=3,ido,2                 
            ic = idp2-i                   
      do 1002 m=1,mp   
            ch(m,i-1,k,1) = cc(m,i-1,1,k)+(cc(m,i-1,3,k)+cc(m,ic-1,2,k))        
     1      +(cc(m,i-1,5,k)+cc(m,ic-1,4,k))                  
            ch(m,i,k,1) = cc(m,i,1,k)+(cc(m,i,3,k)-cc(m,ic,2,k))                
     1      +(cc(m,i,5,k)-cc(m,ic,4,k))   
            ch(m,i-1,k,2) = wa1(i-2)*((cc(m,i-1,1,k)+tr11*   
     1      (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+tr12              
     1      *(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-(ti11*(cc(m,i,3,k)                 
     1      +cc(m,ic,2,k))+ti12*(cc(m,i,5,k)+cc(m,ic,4,k)))) 
     1      -wa1(i-1)*((cc(m,i,1,k)+tr11*(cc(m,i,3,k)-cc(m,ic,2,k))             
     1      +tr12*(cc(m,i,5,k)-cc(m,ic,4,k)))+(ti11*(cc(m,i-1,3,k)              
     1      -cc(m,ic-1,2,k))+ti12*(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))              
            ch(m,i,k,2) = wa1(i-2)*((cc(m,i,1,k)+tr11*(cc(m,i,3,k)              
     1      -cc(m,ic,2,k))+tr12*(cc(m,i,5,k)-cc(m,ic,4,k)))  
     1      +(ti11*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))+ti12       
     1      *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))+wa1(i-1)       
     1      *((cc(m,i-1,1,k)+tr11*(cc(m,i-1,3,k)             
     1      +cc(m,ic-1,2,k))+tr12*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))               
     1      -(ti11*(cc(m,i,3,k)+cc(m,ic,2,k))+ti12           
     1      *(cc(m,i,5,k)+cc(m,ic,4,k))))                    
            ch(m,i-1,k,3) = wa2(i-2)      
     1      *((cc(m,i-1,1,k)+tr12*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))                
     1      +tr11*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-(ti12*(cc(m,i,3,k)            
     1      +cc(m,ic,2,k))-ti11*(cc(m,i,5,k)+cc(m,ic,4,k)))) 
     1     -wa2(i-1)   
     1     *((cc(m,i,1,k)+tr12*(cc(m,i,3,k)-                 
     1      cc(m,ic,2,k))+tr11*(cc(m,i,5,k)-cc(m,ic,4,k)))   
     1      +(ti12*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))-ti11       
     1      *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))                
            ch(m,i,k,3) = wa2(i-2)        
     1     *((cc(m,i,1,k)+tr12*(cc(m,i,3,k)-                 
     1      cc(m,ic,2,k))+tr11*(cc(m,i,5,k)-cc(m,ic,4,k)))   
     1      +(ti12*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))-ti11       
     1      *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))                
     1      +wa2(i-1)                     
     1      *((cc(m,i-1,1,k)+tr12*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))                
     1      +tr11*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-(ti12*(cc(m,i,3,k)            
     1      +cc(m,ic,2,k))-ti11*(cc(m,i,5,k)+cc(m,ic,4,k)))) 
            ch(m,i-1,k,4) = wa3(i-2)      
     1      *((cc(m,i-1,1,k)+tr12*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))                
     1      +tr11*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+(ti12*(cc(m,i,3,k)            
     1      +cc(m,ic,2,k))-ti11*(cc(m,i,5,k)+cc(m,ic,4,k)))) 
     1      -wa3(i-1)                     
     1     *((cc(m,i,1,k)+tr12*(cc(m,i,3,k)-                 
     1      cc(m,ic,2,k))+tr11*(cc(m,i,5,k)-cc(m,ic,4,k)))   
     1      -(ti12*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))-ti11       
     1      *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))                
            ch(m,i,k,4) = wa3(i-2)        
     1     *((cc(m,i,1,k)+tr12*(cc(m,i,3,k)-                 
     1      cc(m,ic,2,k))+tr11*(cc(m,i,5,k)-cc(m,ic,4,k)))   
     1      -(ti12*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))-ti11       
     1      *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))                
     1      +wa3(i-1)                     
     1      *((cc(m,i-1,1,k)+tr12*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))                
     1      +tr11*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+(ti12*(cc(m,i,3,k)            
     1      +cc(m,ic,2,k))-ti11*(cc(m,i,5,k)+cc(m,ic,4,k)))) 
            ch(m,i-1,k,5) = wa4(i-2)      
     1      *((cc(m,i-1,1,k)+tr11*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))                
     1      +tr12*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+(ti11*(cc(m,i,3,k)            
     1      +cc(m,ic,2,k))+ti12*(cc(m,i,5,k)+cc(m,ic,4,k)))) 
     1      -wa4(i-1)                     
     1      *((cc(m,i,1,k)+tr11*(cc(m,i,3,k)-cc(m,ic,2,k))   
     1      +tr12*(cc(m,i,5,k)-cc(m,ic,4,k)))-(ti11*(cc(m,i-1,3,k)              
     1      -cc(m,ic-1,2,k))+ti12*(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))              
            ch(m,i,k,5) = wa4(i-2)        
     1      *((cc(m,i,1,k)+tr11*(cc(m,i,3,k)-cc(m,ic,2,k))   
     1      +tr12*(cc(m,i,5,k)-cc(m,ic,4,k)))-(ti11*(cc(m,i-1,3,k)              
     1      -cc(m,ic-1,2,k))+ti12*(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))              
     1      +wa4(i-1)                     
     1      *((cc(m,i-1,1,k)+tr11*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))                
     1      +tr12*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+(ti11*(cc(m,i,3,k)            
     1      +cc(m,ic,2,k))+ti12*(cc(m,i,5,k)+cc(m,ic,4,k)))) 
 1002          continue                   
  102    continue      
  103 continue         
      return           
      end              
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file vhsgc.f
c
c     this file contains code and documentation for subroutines
c     vhsgc and vhsgci
c
c ... files which must be loaded with vhsgc.f
c
c     sphcom.f, hrfft.f, gaqd.f
c
c     subroutine vhsgc(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
c    +                 mdab,ndab,wvhsgc,lvhsgc,work,lwork,ierror)
c                                                                              
c     subroutine vhsgc performs the vector spherical harmonic synthesis
c     of the arrays br, bi, cr, and ci and stores the result in the
c     arrays v and w. v(i,j) and w(i,j) are the colatitudinal 
c     (measured from the north pole) and east longitudinal components
c     respectively, located at the gaussian colatitude point theta(i)
c     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
c     representation of (v,w) is given below at output parameters v,w.
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are computed
c            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid point
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     ityp   = 0  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon.   
c
c            = 1  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 2  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 3  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 4  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 5  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 6  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 7  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 8  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c
c     nt     the number of syntheses.  in the program that calls vhsgc,
c            the arrays v,w,br,bi,cr, and ci can be three dimensional
c            in which case multiple syntheses will be performed.
c            the third index is the synthesis index which assumes the 
c            values k=1,...,nt.  for a single synthesis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that all the arrays are two
c            dimensional.
c
c     idvw   the first dimension of the arrays v,w as it appears in
c            the program that calls vhsgc. if ityp .le. 2 then idvw
c            must be at least nlat.  if ityp .gt. 2 and nlat is
c            even then idvw must be at least nlat/2. if ityp .gt. 2
c            and nlat is odd then idvw must be at least (nlat+1)/2.
c
c     jdvw   the second dimension of the arrays v,w as it appears in
c            the program that calls vhsgc. jdvw must be at least nlon.
c
c     br,bi  two or three dimensional arrays (see input parameter nt)
c     cr,ci  that contain the vector spherical harmonic coefficients
c            in the spectral representation of v(i,j) and w(i,j) given
c            below at the discription of output parameters v and w.
c
c     mdab   the first dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhsgc. mdab must be at
c            least min0(nlat,nlon/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     ndab   the second dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhsgc. ndab must be at
c            least nlat.
c
c     wvhsgc an array which must be initialized by subroutine vhsgci.
c            once initialized, wvhsgc can be used repeatedly by vhsgc
c            as long as nlon and nlat remain unchanged.  wvhsgc must
c            not be altered between calls of vhsgc.
c
c     lvhsgc the dimension of the array wvhsgc as it appears in the
c            program that calls vhsgc. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhsgc must be at least
c
c               4*nlat*l2+3*max0(l1-2,0)*(2*nlat-l1-1)+nlon+15
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhsgc. define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            if ityp .le. 2 then lwork must be at least
c
c                    nlat*(2*nt*nlon+max0(6*l2,nlon))
c
c            if ityp .gt. 2 then lwork must be at least
c
c                    l2*(2*nt*nlon+max0(6*nlat,nlon))
c
c     **************************************************************
c
c     output parameters
c
c     v,w    two or three dimensional arrays (see input parameter nt)
c            in which the synthesis is stored. v is the colatitudinal
c            component and w is the east longitudinal component. 
c            v(i,j),w(i,j) contain the components at the gaussian
c            colatitude theta(i) and longitude phi(j) = (j-1)*2*pi/nlon.
c            the index ranges are defined above at the input parameter
c            ityp. v and w are computed from the formulas given below.
c
c     define
c
c     1.  theta is colatitude and phi is east longitude
c
c     2.  the normalized associated legendre funnctions
c
c         pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)
c                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
c                        factorial(n)) times the (n+m)th derivative
c                        of (x**2-1)**n with respect to x=cos(theta)
c
c     3.  vbar(m,n,theta) = the derivative of pbar(m,n,theta) with
c                           respect to theta divided by the square
c                           root of n(n+1).
c
c         vbar(m,n,theta) is more easily computed in the form
c
c         vbar(m,n,theta) = (sqrt((n+m)*(n-m+1))*pbar(m-1,n,theta)
c         -sqrt((n-m)*(n+m+1))*pbar(m+1,n,theta))/(2*sqrt(n*(n+1)))
c
c     4.  wbar(m,n,theta) = m/(sin(theta))*pbar(m,n,theta) divided
c                           by the square root of n(n+1).
c
c         wbar(m,n,theta) is more easily computed in the form
c
c         wbar(m,n,theta) = sqrt((2n+1)/(2n-1))*(sqrt((n+m)*(n+m-1))
c         *pbar(m-1,n-1,theta)+sqrt((n-m)*(n-m-1))*pbar(m+1,n-1,theta))
c         /(2*sqrt(n*(n+1)))
c
c
c    the colatitudnal dependence of the normalized surface vector
c                spherical harmonics are defined by
c
c     5.    bbar(m,n,theta) = (vbar(m,n,theta),i*wbar(m,n,theta))
c
c     6.    cbar(m,n,theta) = (i*wbar(m,n,theta),-vbar(m,n,theta))
c
c
c    the coordinate to index mappings 
c
c     7.   phi(j) = (j-1)*2*pi/nlon, theta(i) is the i(th) guassian
c          point (see nlat as an input parameter).
c    
c     the maximum (plus one) longitudinal wave number
c
c     8.     mmax = min0(nlat,nlon/2) if nlon is even or
c            mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c    if we further define the output vector as
c
c     9.    h(i,j) = (v(i,j),w(i,j))
c
c    and the complex coefficients
c
c     10.   b(m,n) = cmplx(br(m+1,n+1),bi(m+1,n+1))
c
c     11.   c(m,n) = cmplx(cr(m+1,n+1),ci(m+1,n+1))
c
c
c    then for i=1,...,nlat and  j=1,...,nlon
c
c        the expansion for real h(i,j) takes the form
c
c     h(i,j) = the sum from n=1 to n=nlat-1 of the real part of
c
c         .5*(b(0,n)*bbar(0,n,theta(i))+c(0,n)*cbar(0,n,theta(i)))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
c     n=nlat-1 of the real part of
c
c              b(m,n)*bbar(m,n,theta(i))*exp(i*m*phi(j))
c             +c(m,n)*cbar(m,n,theta(i))*exp(i*m*phi(j))
c
c   *************************************************************
c
c   in terms of real variables this expansion takes the form
c
c             for i=1,...,nlat and  j=1,...,nlon
c
c     v(i,j) = the sum from n=1 to n=nlat-1 of
c
c               .5*br(1,n+1)*vbar(0,n,theta(i))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
c     n=nlat-1 of the real part of
c
c       (br(m+1,n+1)*vbar(m,n,theta(i))-ci(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *cos(m*phi(j))
c      -(bi(m+1,n+1)*vbar(m,n,theta(i))+cr(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *sin(m*phi(j))
c
c    and for i=1,...,nlat and  j=1,...,nlon
c
c     w(i,j) = the sum from n=1 to n=nlat-1 of
c
c              -.5*cr(1,n+1)*vbar(0,n,theta(i))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
c     n=nlat-1 of the real part of
c
c      -(cr(m+1,n+1)*vbar(m,n,theta(i))+bi(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *cos(m*phi(j))
c      +(ci(m+1,n+1)*vbar(m,n,theta(i))-br(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *sin(m*phi(j))
c
c
c      br(m+1,nlat),bi(m+1,nlat),cr(m+1,nlat), and ci(m+1,nlat) are
c      assumed zero for m even.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of ityp
c            = 4  error in the specification of nt
c            = 5  error in the specification of idvw
c            = 6  error in the specification of jdvw
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lvhsgc
c            = 10 error in the specification of lwork
c
c*************************************************************
c
c     subroutine vhsgci(nlat,nlon,wvhsgc,lvhsgc,dwork,ldwork,ierror)
c
c     subroutine vhsgci initializes the array wvhsgc which can then be
c     used repeatedly by subroutine vhsgc until nlat or nlon is changed.
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are computed
c            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid point
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     lvhsgc the dimension of the array wvhsgc as it appears in the
c            program that calls vhsgc. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhsgc must be at least
c
c               4*nlat*l2+3*max0(l1-2,0)*(2*nlat-l1-1)+nlon+15
c
c     work  a double precision work space that does not need to be saved
c
c     ldwork the dimension of the array dwork as it appears in the
c            program that calls vhsgsi. ldwork must be at least
c
c               2*nlat*(nlat+1)+1
c
c     **************************************************************
c
c     output parameters
c
c     wvhsgc an array which is initialized for use by subroutine vhsgc.
c            once initialized, wvhsgc can be used repeatedly by vhsgc
c            as long as nlat and nlon remain unchanged.  wvhsgc must not
c            be altered between calls of vhsgc.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lvhsgc
c            = 4  error in the specification of ldwork
c
      subroutine vhsgc(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
     +                 mdab,ndab,wvhsgc,lvhsgc,work,lwork,ierror)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          work(1),wvhsgc(1)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      if(ityp.lt.0 .or. ityp.gt.8) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      imid = (nlat+1)/2
      if((ityp.le.2 .and. idvw.lt.nlat) .or.
     1   (ityp.gt.2 .and. idvw.lt.imid)) return
      ierror = 6
      if(jdvw .lt. nlon) return
      ierror = 7
      mmax = min0(nlat,(nlon+1)/2)
      if(mdab .lt. mmax) return
      ierror = 8
      if(ndab .lt. nlat) return
      ierror = 9
      lzz1 = 2*nlat*imid
      labc = 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
c
c     check save work space length
c
      l1 = min0(nlat,(nlon+1)/2)
      l2 = (nlat+1)/2
      lwmin =   4*nlat*l2+3*max0(l1-2,0)*(2*nlat-l1-1)+nlon+15
      if (lvhsgc .lt. lwmin) return


c     if(lvhsgc .lt. 2*(lzz1+labc)+nlon+15) return
      ierror = 10
      if(ityp .le. 2 .and. 
     1         lwork .lt. nlat*(2*nt*nlon+max0(6*imid,nlon))) return
      if(ityp .gt. 2 .and. 
     1         lwork .lt. imid*(2*nt*nlon+max0(6*nlat,nlon))) return
      ierror = 0
      idv = nlat
      if(ityp .gt. 2) idv = imid
      lnl = nt*idv*nlon
      ist = 0
      if(ityp .le. 2) ist = imid
      iw1 = ist+1
      iw2 = lnl+1
      iw3 = iw2+ist
      iw4 = iw2+lnl
      iw5 = iw4+3*imid*nlat
      lzz1 = 2*nlat*imid
      labc = 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      lwzvin = lzz1+labc
      jw1 = lwzvin+1
      jw2 = jw1+lwzvin
      call vhsgc1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,ndab,
     1     br,bi,cr,ci,idv,work,work(iw1),work(iw2),work(iw3),
     2     work(iw4),work(iw5),wvhsgc,wvhsgc(jw1),wvhsgc(jw2))
      return
      end
      subroutine vhsgc1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,
     1   ndab,br,bi,cr,ci,idv,ve,vo,we,wo,vb,wb,wvbin,wwbin,wrfft)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          ve(idv,nlon,1),vo(idv,nlon,1),we(idv,nlon,1),
     3          wo(idv,nlon,1),wvbin(1),wwbin(1),wrfft(1),
     4          vb(imid,nlat,3),wb(imid,nlat,3)
      nlp1 = nlat+1
      mlat = mod(nlat,2)
      mlon = mod(nlon,2)
      mmax = min0(nlat,(nlon+1)/2)
      imm1 = imid
      if(mlat .ne. 0) imm1 = imid-1
      do 10 k=1,nt
      do 10 j=1,nlon
      do 10 i=1,idv
      ve(i,j,k) = 0.
      we(i,j,k) = 0.
   10 continue
      ndo1 = nlat
      ndo2 = nlat
      if(mlat .ne. 0) ndo1 = nlat-1
      if(mlat .eq. 0) ndo2 = nlat-1
   18 itypp = ityp+1
      go to (1,100,200,300,400,500,600,700,800),itypp
c
c     case ityp=0   no symmetries
c
    1 call vbin(0,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 15 k=1,nt
      do 15 np1=2,ndo2,2
      do 15 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1,iv)
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1,iv)
   15 continue
      do 16 k=1,nt
      do 16 np1=3,ndo1,2
      do 16 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1,iv)
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1,iv)
   16 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 30 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(0,nlat,nlon,m,vb,iv,wvbin)
      call wbin(0,nlat,nlon,m,wb,iw,wwbin)
      if(mp1 .gt. ndo1) go to 26
      do 25 k=1,nt
      do 24 np1=mp1,ndo1,2
      do 23 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,np1,iv)
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,np1,iw)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,np1,iv)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,np1,iw)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,np1,iv)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,np1,iw)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,np1,iv)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,np1,iw)
   23 continue
      if(mlat .eq. 0) go to 24
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,np1,iw)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,np1,iw)
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,np1,iw) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,np1,iw)
   24 continue
   25 continue
   26 if(mp2 .gt. ndo2) go to 30
      do 29 k=1,nt
      do 28 np1=mp2,ndo2,2
      do 27 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,np1,iv)
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,np1,iw)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,np1,iv)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,np1,iw)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,np1,iv)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,np1,iw)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,np1,iv)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,np1,iw)
   27 continue
      if(mlat .eq. 0) go to 28
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,np1,iv) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,np1,iv)
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,np1,iv)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,np1,iv)
   28 continue
   29 continue
   30 continue
      go to 950
c
c     case ityp=1   no symmetries,  cr and ci equal zero
c
  100 call vbin(0,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 115 k=1,nt
      do 115 np1=2,ndo2,2
      do 115 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1,iv)
  115 continue
      do 116 k=1,nt
      do 116 np1=3,ndo1,2
      do 116 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1,iv)
  116 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 130 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(0,nlat,nlon,m,vb,iv,wvbin)
      call wbin(0,nlat,nlon,m,wb,iw,wwbin)
      if(mp1 .gt. ndo1) go to 126
      do 125 k=1,nt
      do 124 np1=mp1,ndo1,2
      do 123 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,np1,iv)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,np1,iv)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,np1,iw)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,np1,iw)
  123 continue
      if(mlat .eq. 0) go to 124
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,np1,iw) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,np1,iw)
  124 continue
  125 continue
  126 if(mp2 .gt. ndo2) go to 130
      do 129 k=1,nt
      do 128 np1=mp2,ndo2,2
      do 127 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,np1,iv)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,np1,iv)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,np1,iw)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,np1,iw)
  127 continue
      if(mlat .eq. 0) go to 128
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,np1,iv) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,np1,iv)
  128 continue
  129 continue
  130 continue
      go to 950
c
c     case ityp=2   no symmetries,  br and bi are equal to zero
c
  200 call vbin(0,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 215 k=1,nt
      do 215 np1=2,ndo2,2
      do 215 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1,iv)
  215 continue
      do 216 k=1,nt
      do 216 np1=3,ndo1,2
      do 216 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1,iv)
  216 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 230 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(0,nlat,nlon,m,vb,iv,wvbin)
      call wbin(0,nlat,nlon,m,wb,iw,wwbin)
      if(mp1 .gt. ndo1) go to 226
      do 225 k=1,nt
      do 224 np1=mp1,ndo1,2
      do 223 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,np1,iw)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,np1,iw)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,np1,iv)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,np1,iv)
  223 continue
      if(mlat .eq. 0) go to 224
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,np1,iw)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,np1,iw)
  224 continue
  225 continue
  226 if(mp2 .gt. ndo2) go to 230
      do 229 k=1,nt
      do 228 np1=mp2,ndo2,2
      do 227 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,np1,iw)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,np1,iw)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,np1,iv)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,np1,iv)
  227 continue
      if(mlat .eq. 0) go to 228
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,np1,iv)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,np1,iv)
  228 continue
  229 continue
  230 continue
      go to 950
c
c     case ityp=3   v even,  w odd 
c
  300 call vbin(0,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 315 k=1,nt
      do 315 np1=2,ndo2,2
      do 315 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1,iv)
  315 continue
      do 316 k=1,nt
      do 316 np1=3,ndo1,2
      do 316 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1,iv)
  316 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 330 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(0,nlat,nlon,m,vb,iv,wvbin)
      call wbin(0,nlat,nlon,m,wb,iw,wwbin)
      if(mp1 .gt. ndo1) go to 326
      do 325 k=1,nt
      do 324 np1=mp1,ndo1,2
      do 323 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,np1,iw)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,np1,iw)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,np1,iv)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,np1,iv)
  323 continue
      if(mlat .eq. 0) go to 324
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,np1,iw)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,np1,iw)
  324 continue
  325 continue
  326 if(mp2 .gt. ndo2) go to 330
      do 329 k=1,nt
      do 328 np1=mp2,ndo2,2
      do 327 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,np1,iv)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,np1,iv)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,np1,iw)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,np1,iw)
  327 continue
      if(mlat .eq. 0) go to 328
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,np1,iv) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,np1,iv)
  328 continue
  329 continue
  330 continue
      go to 950
c
c     case ityp=4   v even,  w odd, and both cr and ci equal zero 
c
  400 call vbin(1,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 415 k=1,nt
      do 415 np1=2,ndo2,2
      do 415 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1,iv)
  415 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 430 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(1,nlat,nlon,m,vb,iv,wvbin)
      call wbin(1,nlat,nlon,m,wb,iw,wwbin)
      if(mp2 .gt. ndo2) go to 430
      do 429 k=1,nt
      do 428 np1=mp2,ndo2,2
      do 427 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,np1,iv)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,np1,iv)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,np1,iw)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,np1,iw)
  427 continue
      if(mlat .eq. 0) go to 428
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,np1,iv) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,np1,iv)
  428 continue
  429 continue
  430 continue
      go to 950
c
c     case ityp=5   v even,  w odd,     br and bi equal zero 
c
  500 call vbin(2,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 516 k=1,nt
      do 516 np1=3,ndo1,2
      do 516 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1,iv)
  516 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 530 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(2,nlat,nlon,m,vb,iv,wvbin)
      call wbin(2,nlat,nlon,m,wb,iw,wwbin)
      if(mp1 .gt. ndo1) go to 530
      do 525 k=1,nt
      do 524 np1=mp1,ndo1,2
      do 523 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,np1,iw)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,np1,iw)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,np1,iv)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,np1,iv)
  523 continue
      if(mlat .eq. 0) go to 524
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,np1,iw)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,np1,iw)
  524 continue
  525 continue
  530 continue
      go to 950
c
c     case ityp=6   v odd  ,  w even
c
  600 call vbin(0,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 615 k=1,nt
      do 615 np1=2,ndo2,2
      do 615 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1,iv)
  615 continue
      do 616 k=1,nt
      do 616 np1=3,ndo1,2
      do 616 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1,iv)
  616 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 630 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(0,nlat,nlon,m,vb,iv,wvbin)
      call wbin(0,nlat,nlon,m,wb,iw,wwbin)
      if(mp1 .gt. ndo1) go to 626
      do 625 k=1,nt
      do 624 np1=mp1,ndo1,2
      do 623 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,np1,iv)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,np1,iv)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,np1,iw)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,np1,iw)
  623 continue
      if(mlat .eq. 0) go to 624
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,np1,iw) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,np1,iw)
  624 continue
  625 continue
  626 if(mp2 .gt. ndo2) go to 630
      do 629 k=1,nt
      do 628 np1=mp2,ndo2,2
      do 627 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,np1,iw)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,np1,iw)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,np1,iv)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,np1,iv)
  627 continue
      if(mlat .eq. 0) go to 628
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,np1,iv)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,np1,iv)
  628 continue
  629 continue
  630 continue
      go to 950
c
c     case ityp=7   v odd, w even   cr and ci equal zero
c
  700 call vbin(2,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 716 k=1,nt
      do 716 np1=3,ndo1,2
      do 716 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1,iv)
  716 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 730 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(2,nlat,nlon,m,vb,iv,wvbin)
      call wbin(2,nlat,nlon,m,wb,iw,wwbin)
      if(mp1 .gt. ndo1) go to 730
      do 725 k=1,nt
      do 724 np1=mp1,ndo1,2
      do 723 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,np1,iv)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,np1,iv)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,np1,iw)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,np1,iw)
  723 continue
      if(mlat .eq. 0) go to 724
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,np1,iw) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,np1,iw)
  724 continue
  725 continue
  730 continue
      go to 950
c
c     case ityp=8   v odd,  w even   br and bi equal zero
c
  800 call vbin(1,nlat,nlon,0,vb,iv,wvbin)
c
c     case m = 0
c
      do 815 k=1,nt
      do 815 np1=2,ndo2,2
      do 815 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1,iv)
  815 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 830 mp1=2,mmax
      m = mp1-1
      mp2 = mp1+1
      call vbin(1,nlat,nlon,m,vb,iv,wvbin)
      call wbin(1,nlat,nlon,m,wb,iw,wwbin)
      if(mp2 .gt. ndo2) go to 830
      do 829 k=1,nt
      do 828 np1=mp2,ndo2,2
      do 827 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,np1,iw)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,np1,iw)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,np1,iv)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,np1,iv)
  827 continue
      if(mlat .eq. 0) go to 828
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,np1,iv)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,np1,iv)
  828 continue
  829 continue
  830 continue
  950 do 14 k=1,nt
      call hrfftb(idv,nlon,ve(1,1,k),idv,wrfft,vb)
      call hrfftb(idv,nlon,we(1,1,k),idv,wrfft,vb)
   14 continue
      if(ityp .gt. 2) go to 12
      do 60 k=1,nt
      do 60 j=1,nlon
      do 60 i=1,imm1
      v(i,j,k) = .5*(ve(i,j,k)+vo(i,j,k))
      w(i,j,k) = .5*(we(i,j,k)+wo(i,j,k))
      v(nlp1-i,j,k) = .5*(ve(i,j,k)-vo(i,j,k))
      w(nlp1-i,j,k) = .5*(we(i,j,k)-wo(i,j,k))
   60 continue
      go to 13
   12 do 11 k=1,nt
      do 11 j=1,nlon
      do 11 i=1,imm1
      v(i,j,k) = .5*ve(i,j,k)
      w(i,j,k) = .5*we(i,j,k)
   11 continue
   13 if(mlat .eq. 0) return
      do 65 k=1,nt
      do 65 j=1,nlon
      v(imid,j,k) = .5*ve(imid,j,k)
      w(imid,j,k) = .5*we(imid,j,k)
   65 continue
      return
      end
      subroutine vhsgci(nlat,nlon,wvhsgc,lvhsgc,dwork,ldwork,ierror)
      dimension wvhsgc(lvhsgc)
      double precision dwork(ldwork)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      imid = (nlat+1)/2
      lzz1 = 2*nlat*imid
      mmax = min0(nlat,(nlon+1)/2)
      labc = 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      if(lvhsgc .lt. 2*(lzz1+labc)+nlon+15) return
      ierror = 4
      if (ldwork .lt. 2*nlat*(nlat+1)+1) return
      ierror = 0
c
c     compute gaussian points in first nlat+1 words of dwork
c     double precision
c
c     lwk = 2*nlat*(nlat+2)
      jw1 = 1
      jw2 = jw1+nlat
      jw3 = jw2+nlat
c     jw2 = jw1+nlat+nlat
c     jw3 = jw2+nlat+nlat
      call gaqd(nlat,dwork(jw1),dwork(jw2),dwork(jw3),ldwork,ierror)
c     iwrk = nlat+2
      iwrk = (nlat+1)/2 + 1
      call vbgint (nlat,nlon,dwork,wvhsgc,dwork(iwrk))
      lwvbin = lzz1+labc
      iw1 = lwvbin+1
      call wbgint (nlat,nlon,dwork,wvhsgc(iw1),dwork(iwrk))
      iw2 = iw1+lwvbin
      call hrffti(nlon,wvhsgc(iw2))
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file shsgs.f
c
c     this file contains code and documentation for subroutines
c     shsgs and shsgsi
c
c ... files which must be loaded with shsgs.f
c
c     sphcom.f, hrfft.f, gaqd.f
c
c     subroutine shsgs(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
c    1                    wshsgs,lshsgs,work,lwork,ierror)
c
c     subroutine shsgs performs the spherical harmonic synthesis
c     on the arrays a and b and stores the result in the array g.
c     the synthesis is performed on an equally spaced longitude grid
c     and a gaussian colatitude grid.  the associated legendre functions
c     are stored rather than recomputed as they are in subroutine
c     shsgc.  the synthesis is described below at output parameter
c     g.
c
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are compu
c            in radians in theta(1),...,theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid poi
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c     isym   = 0  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 array g(i,j) for i=1,...,nlat and j=1,...,nlon.
c                 (see description of g below)
c
c            = 1  g is antisymmetric about the equator. the synthesis
c                 is performed on the northern hemisphere only.  i.e.
c                 if nlat is odd the synthesis is performed on the
c                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
c                 if nlat is even the synthesis is performed on the
c                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c
c            = 2  g is symmetric about the equator. the synthesis is
c                 performed on the northern hemisphere only.  i.e.
c                 if nlat is odd the synthesis is performed on the
c                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
c                 if nlat is even the synthesis is performed on the
c                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c     nt     the number of syntheses.  in the program that calls shsgs,
c            the arrays g,a and b can be three dimensional in which
c            case multiple synthesis will be performed.  the third
c            index is the synthesis index which assumes the values
c            k=1,...,nt.  for a single synthesis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that the arrays g,a and b
c            have only two dimensions.
c
c     idg    the first dimension of the array g as it appears in the
c            program that calls shagc. if isym equals zero then idg
c            must be at least nlat.  if isym is nonzero then idg must
c            be at least nlat/2 if nlat is even or at least (nlat+1)/2
c            if nlat is odd.
c
c     jdg    the second dimension of the array g as it appears in the
c            program that calls shagc. jdg must be at least nlon.
c
c     a,b    two or three dimensional arrays (see the input parameter
c            nt) that contain the coefficients in the spherical harmonic
c            expansion of g(i,j) given below at the definition of the
c            output parameter g.  a(m,n) and b(m,n) are defined for
c            indices m=1,...,mmax and n=m,...,nlat where mmax is the
c            maximum (plus one) longitudinal wave number given by
c            mmax = min0(nlat,(nlon+2)/2) if nlon is even or
c            mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     mdab   the first dimension of the arrays a and b as it appears
c            in the program that calls shsgs. mdab must be at least
c            min0((nlon+2)/2,nlat) if nlon is even or at least
c            min0((nlon+1)/2,nlat) if nlon is odd.
c
c     ndab   the second dimension of the arrays a and b as it appears
c            in the program that calls shsgs. ndab must be at least nlat
c
c     wshsgs an array which must be initialized by subroutine shsgsi.
c            once initialized, wshsgs can be used repeatedly by shsgs
c            as long as nlat and nlon remain unchanged.  wshsgs must
c            not be altered between calls of shsgs.
c
c     lshsgs the dimension of the array wshsgs as it appears in the
c            program that calls shsgs. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshsgs must be at least
c
c            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
c
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls shsgs. define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c
c            if isym is zero then lwork must be at least
c
c                  nlat*nlon*(nt+1)
c
c            if isym is nonzero then lwork must be at least
c
c                  l2*nlon*(nt+1)
c
c
c     **************************************************************
c
c     output parameters
c
c     g      a two or three dimensional array (see input parameter nt)
c            that contains the discrete function which is synthesized.
c            g(i,j) contains the value of the function at the gaussian
c            colatitude point theta(i) and longitude point
c            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined
c            above at the input parameter isym.  for isym=0, g(i,j)
c            is given by the the equations listed below.  symmetric
c            versions are used when isym is greater than zero.
c
c     the normalized associated legendre functions are given by
c
c     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
c                       *sin(theta)**m/(2**n*factorial(n)) times the
c                       (n+m)th derivative of (x**2-1)**n with respect
c                       to x=cos(theta)
c
c     define the maximum (plus one) longitudinal wave number
c     as   mmax = min0(nlat,(nlon+2)/2) if nlon is even or
c          mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     then g(i,j) = the sum from n=0 to n=nlat-1 of
c
c                   .5*pbar(0,n,theta(i))*a(1,n+1)
c
c              plus the sum from m=1 to m=mmax-1 of
c
c                   the sum from n=m to n=nlat-1 of
c
c              pbar(m,n,theta(i))*(a(m+1,n+1)*cos(m*phi(j))
c                                    -b(m+1,n+1)*sin(m*phi(j)))
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of isym
c            = 4  error in the specification of nt
c            = 5  error in the specification of idg
c            = 6  error in the specification of jdg
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lshsgs
c            = 10 error in the specification of lwork
c
c
c ****************************************************************
c
c     subroutine shsgsi(nlat,nlon,wshsgs,lshsgs,work,lwork,dwork,ldwork,
c    +                  ierror)
c
c     subroutine shsgsi initializes the array wshsgs which can then
c     be used repeatedly by subroutines shsgs. it precomputes
c     and stores in wshsgs quantities such as gaussian weights,
c     legendre polynomial coefficients, and fft trigonometric tables.
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are compu
c            in radians in theta(1),...,theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid poi
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c     wshsgs an array which must be initialized by subroutine shsgsi.
c            once initialized, wshsgs can be used repeatedly by shsgs
c            as long as nlat and nlon remain unchanged.  wshsgs must
c            not be altered between calls of shsgs.
c
c     lshsgs the dimension of the array wshsgs as it appears in the
c            program that calls shsgs. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshsgs must be at least
c
c            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
c
c     work   a real work space which need not be saved
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls shsgsi. lwork must be at least
c            4*nlat*(nlat+2)+2 in the routine calling shsgsi
c
c     dwork   a double precision work array that does not have to be saved.
c
c     ldwork  the length of dwork in the calling routine.  ldwork must
c             be at least nlat*(nlat+4)
c
c     output parameter
c
c     wshsgs an array which must be initialized before calling shsgs or
c            once initialized, wshsgs can be used repeatedly by shsgs or
c            as long as nlat and nlon remain unchanged.  wshsgs must not
c            altered between calls of shsgs.
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lshsgs
c            = 4  error in the specification of lwork
c            = 5  error in the specification of ldwork
c            = 5  failure in gaqd to compute gaussian points
c                 (due to failure in eigenvalue routine)
c
c
      subroutine shsgs(nlat,nlon,mode,nt,g,idg,jdg,a,b,mdab,ndab,
     1                    wshsgs,lshsgs,work,lwork,ierror)
      dimension g(idg,jdg,1),a(mdab,ndab,1),b(mdab,ndab,1),
     1          wshsgs(lshsgs),work(lwork)
c     check input parameters
      ierror = 1
      if (nlat.lt.3) return
      ierror = 2
      if (nlon.lt.4) return
      ierror = 3
      if (mode.lt.0 .or.mode.gt.2) return
      ierror = 4
      if (nt.lt.1) return
c     set limit on m subscript
      l = min0((nlon+2)/2,nlat)
c     set gaussian point nearest equator pointer
      late = (nlat+mod(nlat,2))/2
c     set number of grid points for analysis/synthesis
      lat = nlat
      if (mode.ne.0) lat = late
      ierror = 5
      if (idg.lt.lat) return
      ierror = 6
      if (jdg.lt.nlon) return
      ierror = 7
      if(mdab .lt. l) return
      ierror = 8
      if(ndab .lt. nlat) return
      l1 = l
      l2 = late
      ierror = 9
c     check permanent work space length
      lp=nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
      if(lshsgs.lt.lp) return
c     check temporary work space length
      ierror = 10
      if (mode.eq.0 .and. lwork.lt.nlat*nlon*(nt+1)) return
      if (mode.ne.0 .and. lwork.lt.l2*nlon*(nt+1)) return
      ierror = 0
c     starting address for fft values and legendre polys in wshsgs
      ifft = nlat+2*nlat*late+3*(l*(l-1)/2+(nlat-l)*(l-1))+1
      ipmn = ifft+nlon+15
c     set pointer for internal storage of g
      iw = lat*nlon*nt+1
      call shsgs1(nlat,nlon,l,lat,mode,g,idg,jdg,nt,a,b,mdab,ndab,
     1          wshsgs(ifft),wshsgs(ipmn),late,work,work(iw))
      return
      end

      subroutine shsgs1(nlat,nlon,l,lat,mode,gs,idg,jdg,nt,a,b,mdab,
     1                  ndab,wfft,pmn,late,g,work)
      dimension gs(idg,jdg,nt),a(mdab,ndab,nt),b(mdab,ndab,nt)
      dimension wfft(1),pmn(late,1),g(lat,nlon,nt),work(1)

c     reconstruct fourier coefficients in g on gaussian grid
c     using coefficients in a,b

c     initialize to zero
      do 100 k=1,nt
      do 100 j=1,nlon
      do 100 i=1,lat
      g(i,j,k) = 0.0
  100 continue

      lm1 = l
      if (nlon .eq. l+l-2) lm1 = l-1
      if (mode.eq.0) then
c     set first column in g
      m = 0
      mml1 = m*(2*nlat-m-1)/2
      do 101 k=1,nt
c     n even
      do 102 np1=1,nlat,2
      mn = mml1+np1
      do 102 i=1,late
      g(i,1,k) = g(i,1,k)+a(1,np1,k)*pmn(i,mn)
  102 continue
c     n odd
      nl2 = nlat/2
      do 103 np1=2,nlat,2
      mn = mml1+np1
      do 103 i=1,nl2
      is = nlat-i+1
      g(is,1,k) = g(is,1,k)+a(1,np1,k)*pmn(i,mn)
  103 continue
  101 continue

c     restore m=0 coefficients from odd/even
      do 112 k=1,nt
      do 112 i=1,nl2
      is = nlat-i+1
      t1 = g(i,1,k)
      t3 = g(is,1,k)
      g(i,1,k) = t1+t3
      g(is,1,k) = t1-t3
  112 continue

c     sweep interior columns of g
      do 104 mp1=2,lm1
      m = mp1-1
      mml1 = m*(2*nlat-m-1)/2
      mp2 = m+2
      do 105 k=1,nt
c     for n-m even store (g(i,p,k)+g(nlat-i+1,p,k))/2 in g(i,p,k) p=2*m,2*m+1
c     for i=1,...,late
      do 106 np1=mp1,nlat,2
      mn = mml1+np1
      do 107 i=1,late
      g(i,2*m,k) = g(i,2*m,k)+a(mp1,np1,k)*pmn(i,mn)
      g(i,2*m+1,k) = g(i,2*m+1,k)+b(mp1,np1,k)*pmn(i,mn)
  107 continue
  106 continue

c     for n-m odd store g(i,p,k)-g(nlat-i+1,p,k) in g(nlat-i+1,p,k)
c     for i=1,...,nlat/2 (p=2*m,p=2*m+1)
      do 108 np1=mp2,nlat,2
      mn = mml1+np1
      do 109 i=1,nl2
      is = nlat-i+1
      g(is,2*m,k) = g(is,2*m,k)+a(mp1,np1,k)*pmn(i,mn)
      g(is,2*m+1,k) = g(is,2*m+1,k)+b(mp1,np1,k)*pmn(i,mn)
  109 continue
  108 continue

c     now set fourier coefficients using even-odd reduction above
      do 110 i=1,nl2
      is = nlat-i+1
      t1 = g(i,2*m,k)
      t2 = g(i,2*m+1,k)
      t3 = g(is,2*m,k)
      t4 = g(is,2*m+1,k)
      g(i,2*m,k) = t1+t3
      g(i,2*m+1,k) = t2+t4
      g(is,2*m,k) = t1-t3
      g(is,2*m+1,k) = t2-t4
  110 continue

  105 continue
  104 continue

c     set last column (using a only) if necessary
      if (nlon.eq. l+l-2) then
      m = l-1
      mml1 = m*(2*nlat-m-1)/2
      do 111 k=1,nt
c     n-m even
      do 131 np1=l,nlat,2
      mn = mml1+np1
      do 131 i=1,late
      g(i,nlon,k) = g(i,nlon,k)+2.0*a(l,np1,k)*pmn(i,mn)

  131 continue
      lp1 = l+1
c     n-m odd
      do 132 np1=lp1,nlat,2
      mn = mml1+np1
      do 132 i=1,nl2
      is = nlat-i+1
      g(is,nlon,k) = g(is,nlon,k)+2.0*a(l,np1,k)*pmn(i,mn)
  132 continue
      do 133 i=1,nl2
      is = nlat-i+1
      t1 = g(i,nlon,k)
      t3 = g(is,nlon,k)
      g(i,nlon,k)= t1+t3
      g(is,nlon,k)= t1-t3
  133 continue
  111 continue
      end if

      else
c     half sphere (mode.ne.0)
c     set first column in g
      m = 0
      mml1 = m*(2*nlat-m-1)/2
      meo = 1
      if (mode.eq.1) meo = 2
      ms = m+meo
      do 113 k=1,nt
      do 113 np1=ms,nlat,2
      mn = mml1+np1
      do 113 i=1,late
      g(i,1,k) = g(i,1,k)+a(1,np1,k)*pmn(i,mn)
  113 continue

c     sweep interior columns of g

      do 114 mp1=2,lm1
      m = mp1-1
      mml1 = m*(2*nlat-m-1)/2
      ms = m+meo
      do 115 k=1,nt
      do 115 np1=ms,nlat,2
      mn = mml1+np1
      do 115 i=1,late
      g(i,2*m,k) = g(i,2*m,k)+a(mp1,np1,k)*pmn(i,mn)
      g(i,2*m+1,k) = g(i,2*m+1,k)+b(mp1,np1,k)*pmn(i,mn)
  115 continue
  114 continue

      if (nlon.eq.l+l-2) then
c     set last column
      m = l-1
      mml1 = m*(2*nlat-m-1)/2
      ns = l
      if (mode.eq.1) ns = l+1
      do 116 k=1,nt
      do 116 np1=ns,nlat,2
      mn = mml1+np1
      do 116 i=1,late
      g(i,nlon,k) = g(i,nlon,k)+2.0*a(l,np1,k)*pmn(i,mn)
  116 continue
      end if

      end if


c     do inverse fourier transform
      do 120 k=1,nt
      call hrfftb(lat,nlon,g(1,1,k),lat,wfft,work)
  120 continue
c     scale output in gs
      do 122 k=1,nt
      do 122 j=1,nlon
      do 122 i=1,lat
      gs(i,j,k) = 0.5*g(i,j,k)
  122 continue

      return
      end
      subroutine shsgsi(nlat,nlon,wshsgs,lshsgs,work,lwork,dwork,ldwork,
     +                  ierror)
c
c     this subroutine must be called before calling shags or shsgs with
c     fixed nlat,nlon. it precomputes the gaussian weights, points
c     and all necessary legendre polys and stores them in wshsgs.
c     these quantities must be preserved when calling shsgs
c     repeatedly with fixed nlat,nlon.
c
      dimension wshsgs(lshsgs),work(lwork)
      double precision dwork(ldwork)
      ierror = 1
      if (nlat.lt.3) return
      ierror = 2
      if (nlon.lt.4) return
c     set triangular truncation limit for spherical harmonic basis
      l = min0((nlon+2)/2,nlat)
c     set equator or nearest point (if excluded) pointer
      late = (nlat+1)/2
      l1 = l
      l2 = late
c     check permanent work space length
      ierror = 3
      lp=nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
      if(lshsgs.lt.lp) return
      ierror = 4
c     check temporary work space
      if (lwork.lt.4*nlat*(nlat+2)+2) return
      ierror = 5
      if (ldwork .lt. nlat*(nlat+4)) return
      ierror = 0
c     set preliminary quantites needed to compute and store legendre polys
      ldw = nlat*(nlat+4)
      call shsgsp(nlat,nlon,wshsgs,lshsgs,dwork,ldwork,ierror)
      if (ierror.ne.0) return
c     set legendre poly pointer in wshsgs
      ipmnf = nlat+2*nlat*late+3*(l*(l-1)/2+(nlat-l)*(l-1))+nlon+16
      call shsgss1(nlat,l,late,wshsgs,work,wshsgs(ipmnf))
      return
      end
      subroutine shsgss1(nlat,l,late,w,pmn,pmnf)
      dimension w(1),pmn(nlat,late,3),pmnf(late,1)
c     compute and store legendre polys for i=1,...,late,m=0,...,l-1
c     and n=m,...,l-1
      do i=1,nlat
	do j=1,late
	  do k=1,3
	    pmn(i,j,k) = 0.0
	  end do
	 end do
      end do
      do 100 mp1=1,l
      m = mp1-1
      mml1 = m*(2*nlat-m-1)/2
c     compute pmn for n=m,...,nlat-1 and i=1,...,(l+1)/2
      mode = 0
      call legin(mode,l,nlat,m,w,pmn,km)
c     store above in pmnf
      do 101 np1=mp1,nlat
      mn = mml1+np1
      do 102 i=1,late
      pmnf(i,mn) = pmn(np1,i,km)
  102 continue
  101 continue
  100 continue
      return
      end
      subroutine shsgsp(nlat,nlon,wshsgs,lshsgs,dwork,ldwork,ierror)
      dimension wshsgs(lshsgs)
      double precision dwork(ldwork)
      ierror = 1
      if (nlat.lt.3) return
      ierror = 2
      if (nlon.lt.4) return
c     set triangular truncation limit for spherical harmonic basis
      l = min0((nlon+2)/2,nlat)
c     set equator or nearest point (if excluded) pointer
      late = (nlat+mod(nlat,2))/2
      l1 = l
      l2 = late
      ierror = 3
c     check permanent work space length
      if (lshsgs .lt. nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
      ierror = 4
c     if (lwork.lt.4*nlat*(nlat+2)+2) return
      if (ldwork .lt. nlat*(nlat+4)) return
      ierror = 0
c     set pointers
      i1 = 1
      i2 = i1+nlat
      i3 = i2+nlat*late
      i4 = i3+nlat*late
      i5 = i4+l*(l-1)/2 +(nlat-l)*(l-1)
      i6 = i5+l*(l-1)/2 +(nlat-l)*(l-1)
      i7 = i6+l*(l-1)/2 +(nlat-l)*(l-1)
c     set indices in temp work for double precision gaussian wts and pts
      idth = 1
c     idwts = idth+2*nlat
c     iw = idwts+2*nlat
      idwts = idth+nlat
      iw = idwts+nlat
      call shsgsp1(nlat,nlon,l,late,wshsgs(i1),wshsgs(i2),wshsgs(i3),
     1wshsgs(i4),wshsgs(i5),wshsgs(i6),wshsgs(i7),dwork(idth),
     2dwork(idwts),dwork(iw),ierror)
      if (ierror.ne.0) ierror = 6
      return
      end
      subroutine shsgsp1(nlat,nlon,l,late,wts,p0n,p1n,abel,bbel,cbel,
     +                   wfft,dtheta,dwts,work,ier)
      dimension wts(nlat),p0n(nlat,late),p1n(nlat,late),abel(1),bbel(1),
     1 cbel(1),wfft(1),dtheta(nlat),dwts(nlat)
      double precision pb,dtheta,dwts,work(*)
      indx(m,n) = (n-1)*(n-2)/2+m-1
      imndx(m,n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1
      call hrffti(nlon,wfft)
c
c     compute double precision gaussian points and weights
c
      lw = nlat*(nlat+2)
      call gaqd(nlat,dtheta,dwts,work,lw,ier)
      if (ier.ne.0) return

c     store gaussian weights single precision to save computation
c     in inner loops in analysis
      do 100 i=1,nlat
      wts(i) = dwts(i)
  100 continue

c     initialize p0n,p1n using double precision dnlfk,dnlft
      do 101 np1=1,nlat
      do 101 i=1,late
      p0n(np1,i) = 0.0
      p1n(np1,i) = 0.0
  101 continue
c     compute m=n=0 legendre polynomials for all theta(i)
      np1 = 1
      n = 0
      m = 0
      call dnlfk(m,n,work)
      do 103 i=1,late
      call dnlft(m,n,dtheta(i),work,pb)
      p0n(1,i) = pb
  103 continue
c     compute p0n,p1n for all theta(i) when n.gt.0
      do 104 np1=2,nlat
      n = np1-1
      m = 0
      call dnlfk(m,n,work)
      do 105 i=1,late
      call dnlft(m,n,dtheta(i),work,pb)
      p0n(np1,i) = pb
  105 continue
c     compute m=1 legendre polynomials for all n and theta(i)
      m = 1
      call dnlfk(m,n,work)
      do 106 i=1,late
      call dnlft(m,n,dtheta(i),work,pb)
      p1n(np1,i) = pb
  106 continue
  104 continue
c
c     compute and store swarztrauber recursion coefficients
c     for 2.le.m.le.n and 2.le.n.le.nlat in abel,bbel,cbel
      do 107 n=2,nlat
      mlim = min0(n,l)
      do 107 m=2,mlim
      imn = indx(m,n)
      if (n.ge.l) imn = imndx(m,n)
      abel(imn)=sqrt(float((2*n+1)*(m+n-2)*(m+n-3))/
     1               float(((2*n-3)*(m+n-1)*(m+n))))
      bbel(imn)=sqrt(float((2*n+1)*(n-m-1)*(n-m))/
     1               float(((2*n-3)*(m+n-1)*(m+n))))
      cbel(imn)=sqrt(float((n-m+1)*(n-m+2))/
     1               float(((n+m-1)*(n+m))))
  107 continue
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                          SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c                             August 2003
c
c
c     This version of gaqd implements the method presented in:
c     P. N. swarztrauber, Computing the points and weights for 
c     Gauss-Legendre quadrature, SIAM J. Sci. Comput.,
c     24(2002) pp. 945-954.
c     
c     It the version that is new to spherepack 3.1
c     The w and lwork arrays are dummy and included only to
c     permit a simple pluggable exchange with the
c     old gaqd in spherepack 3.0. 
c
c
c
      subroutine gaqd(nlat,theta,wts,w,lwork,ierror)
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 2001 by ucar                 .
c  .                                                             .
c  .       university corporation for atmospheric research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c                             February 2002
c                        
c     gauss points and weights are computed using the fourier-newton
c     described in "on computing the points and weights for 
c     gauss-legendre quadrature", paul n. swarztrauber, siam journal 
c     on scientific computing that has been accepted for publication.
c     This routine is faster and more accurate than older program 
c     with the same name.
c
c     subroutine gaqd computes the nlat gaussian colatitudes and weights
c     in double precision. the colatitudes are in radians and lie in the
c     in the interval (0,pi).
c
c     input parameters
c
c     nlat    the number of gaussian colatitudes in the interval (0,pi)
c             (between the two poles).  nlat must be greater than zero.
c
c     w       unused double precision variable that permits a simple 
c             exchange with the old routine with the same name 
c             in spherepack.
c
c     lwork   unused variable that permits a simple exchange with the
c             old routine with the same name in spherepack.
c
c     output parameters
c
c     theta   a double precision array with length nlat
c             containing the gaussian colatitudes in
c             increasing radians on the interval (0,pi).
c
c     wts     a double precision array with lenght nlat
c             containing the gaussian weights.
c
c     ierror = 0 no errors
c            = 1 if nlat.le.0
c
c  *****************************************************************
c
      double precision theta(nlat),wts(nlat),w,
     1 x,pi,pis2,dtheta,dthalf,cmax,zprev,zlast,zero,
     2 zhold,pb,dpb,dcor,sum,cz
c
c     check work space length
c
      ierror = 1
      if (nlat.le.0) return
      ierror = 0
c
c     compute weights and points analytically when nlat=1,2
c
      if (nlat.eq.1) then
      theta(1) = dacos(0.0d0)
      wts(1) = 2.0d0
      return
      end if
      if (nlat.eq.2) then
      x = dsqrt(1.0d0/3.0d0)
      theta(1) = dacos(x)
      theta(2) = dacos(-x)
      wts(1) = 1.0d0
      wts(2) = 1.0d0
      return
      end if
      eps = sqrt(dzeps(1.0d0))
      eps = eps*sqrt(eps)
      pis2 = 2.0d0*datan(1.0d0)
      pi = pis2+pis2 
      mnlat = mod(nlat,2)
      ns2 = nlat/2
      nhalf = (nlat+1)/2
      idx = ns2+2
c
      call cpdp (nlat,cz,theta(ns2+1),wts(ns2+1))
c
      dtheta = pis2/nhalf
      dthalf = dtheta/2.0d0
      cmax = .2d0*dtheta
c
c     estimate first point next to theta = pi/2
c
      if(mnlat.ne.0) then
      zero = pis2-dtheta
      zprev = pis2
      nix = nhalf-1
      else
      zero = pis2-dthalf
      nix = nhalf
      end if
 9    it = 0
 10   it = it+1
      zlast = zero
c
c     newton iterations
c
      call tpdp (nlat,zero,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
      dcor = pb/dpb
      sgnd = 1.0
      if(dcor .ne. 0.0d0) sgnd = dcor/dabs(dcor)
      dcor = sgnd*min(dabs(dcor),cmax)
      zero = zero-dcor
      if(dabs(zero-zlast).gt.eps*dabs(zero)) go to 10
      theta(nix) = zero
      zhold = zero
c      wts(nix) = (nlat+nlat+1)/(dpb*dpb)
c    
c     yakimiw's formula permits using old pb and dpb
c
      wts(nix) = (nlat+nlat+1)/(dpb+pb*dcos(zlast)/dsin(zlast))**2
      nix = nix-1
      if(nix.eq.0) go to 30
      if(nix.eq.nhalf-1)  zero = 3.0*zero-pi
      if(nix.lt.nhalf-1)  zero = zero+zero-zprev
      zprev = zhold
      go to 9
c
c     extend points and weights via symmetries
c
 30   if(mnlat.ne.0) then
      theta(nhalf) = pis2
      call tpdp (nlat,pis2,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
      wts(nhalf) = (nlat+nlat+1)/(dpb*dpb)
      end if
      do i=1,ns2
      wts(nlat-i+1) = wts(i)
      theta(nlat-i+1) = pi-theta(i)
      end do
      sum = 0.0d0
      do i=1,nlat
      sum = sum+wts(i)
      end do
      do i=1,nlat
      wts(i) = 2.0d0*wts(i)/sum
      end do
      return
      end
      subroutine cpdp(n,cz,cp,dcp)
c
c     computes the fourier coefficients of the legendre
c     polynomial p_n^0 and its derivative. 
c     n is the degree and n/2 or (n+1)/2
c     coefficients are returned in cp depending on whether
c     n is even or odd. The same number of coefficients
c     are returned in dcp. For n even the constant 
c     coefficient is returned in cz. 
c
      double precision cp(n/2+1),dcp(n/2+1),
     1 t1,t2,t3,t4,cz
      ncp = (n+1)/2
      t1 = -1.0d0
      t2 = n+1.0d0
      t3 = 0.0d0
      t4 = n+n+1.0d0
      if(mod(n,2).eq.0) then
      cp(ncp) = 1.0d0
      do j = ncp,2,-1
      t1 = t1+2.0d0
      t2 = t2-1.0d0
      t3 = t3+1.0d0
      t4 = t4-2.0d0
      cp(j-1) = (t1*t2)/(t3*t4)*cp(j)
      end do
      t1 = t1+2.0d0
      t2 = t2-1.0d0
      t3 = t3+1.0d0
      t4 = t4-2.0d0
      cz = (t1*t2)/(t3*t4)*cp(1)
      do j=1,ncp
      dcp(j) = (j+j)*cp(j)
      end do
      else
      cp(ncp) = 1.0d0
      do j = ncp-1,1,-1
      t1 = t1+2.0d0
      t2 = t2-1.0d0
      t3 = t3+1.0d0
      t4 = t4-2.0d0
      cp(j) = (t1*t2)/(t3*t4)*cp(j+1)
      end do
      do j=1,ncp
      dcp(j) = (j+j-1)*cp(j)
      end do
      end if
      return
      end
      subroutine tpdp (n,theta,cz,cp,dcp,pb,dpb)
c
c     computes pn(theta) and its derivative dpb(theta) with 
c     respect to theta
c
      double precision cp(n/2+1),dcp(n/2+1),cz,
     1  pb,dpb,fn,theta,cdt,sdt,cth,sth,chh
c
      fn = n
      cdt = dcos(theta+theta)
      sdt = dsin(theta+theta)
      if(mod(n,2) .eq.0) then
c
c     n even
c
      kdo = n/2
      pb = .5d0*cz
      dpb = 0.0d0
      if(n .gt. 0) then
      cth = cdt
      sth = sdt
      do 170 k=1,kdo
c      pb = pb+cp(k)*cos(2*k*theta)
      pb = pb+cp(k)*cth
c      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta)
      dpb = dpb-dcp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  170 continue
      end if
      else
c
c     n odd
c
      kdo = (n+1)/2
      pb = 0.0d0
      dpb = 0.0d0
      cth = dcos(theta)
      sth = dsin(theta)
      do 190 k=1,kdo
c      pb = pb+cp(k)*cos((2*k-1)*theta)
      pb = pb+cp(k)*cth
c      dpb = dpb-(k+k-1)*cp(k)*sin((2*k-1)*theta)
      dpb = dpb-dcp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  190 continue
      end if
      return
      end
      real function dzeps (x)
      double precision  x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = abs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      dzeps = eps*dabs(x)
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file vrtgs.f
c
c     this file includes documentation and code for
c     subroutine divgs          i
c
c ... files which must be loaded with vrtgs.f
c
c     sphcom.f, hrfft.f, vhgsc.f, shsgs.f, gaqd.f
c
c     subroutine vrtgs(nlat,nlon,isym,nt,vort,ivrt,jvrt,cr,ci,mdc,ndc,
c    +                 wshsgs,lshsgs,work,lwork,ierror)
c
c     given the vector spherical harmonic coefficients cr and ci, precomputed
c     by subroutine vhags for a vector field (v,w), subroutine vrtgs
c     computes the vorticity of the vector field in the scalar array
c     vort.  vort(i,j) is the vorticity at the gaussian colatitude
c     theta(i) (see nlat as input parameter) and longitude
c     lambda(j) = (j-1)*2*pi/nlon on the sphere.  i.e.,
c
c            vort(i,j) =  [-dv/dlambda + d(sint*w)/dtheta]/sint
c
c     where sint = sin(theta(i)).  w is the east longitudinal and v
c     is the colatitudinal component of the vector field from which
c     cr,ci were precomputed.  required associated legendre polynomials
c     are stored rather than recomputed as they are in subroutine vrtgc.
c
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are computed
c            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid point
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than 3. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c
c     isym   a parameter which determines whether the vorticity is
c            computed on the full or half sphere as follows:
c
c      = 0
c            the symmetries/antsymmetries described in isym=1,2 below
c            do not exist in (v,w) about the equator.  in this case the
c            vorticity is neither symmetric nor antisymmetric about
c            the equator.  the vorticity is computed on the entire
c            sphere.  i.e., in the array vort(i,j) for i=1,...,nlat and
c            j=1,...,nlon.
c
c      = 1
c            w is antisymmetric and v is symmetric about the equator.
c            in this case the vorticity is symmetyric about the
c            equator and is computed for the northern hemisphere
c            only.  i.e., if nlat is odd the vorticity is computed
c            in the array vort(i,j) for i=1,...,(nlat+1)/2 and for
c            j=1,...,nlon.  if nlat is even the vorticity is computed
c            in the array vort(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c      = 2
c            w is symmetric and v is antisymmetric about the equator
c            in this case the vorticity is antisymmetric about the
c            equator and is computed for the northern hemisphere
c            only.  i.e., if nlat is odd the vorticity is computed
c            in the array vort(i,j) for i=1,...,(nlat+1)/2 and for
c            j=1,...,nlon.  if nlat is even the vorticity is computed
c            in the array vort(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c
c      nt    nt is the number of scalar and vector fields.  some
c            computational efficiency is obtained for multiple fields.
c            in the program that calls vrtgs, the arrays cr,ci, and vort
c            can be three dimensional corresponding to an indexed multiple
c            vector field.  in this case multiple scalar synthesis will
c            be performed to compute the vorticity for each field.  the
c            third index is the synthesis index which assumes the values
c            k=1,...,nt.  for a single synthesis set nt = 1.  the
c            description of the remaining parameters is simplified by
c            assuming that nt=1 or that all the arrays are two dimensional.
c
c     ivrt   the first dimension of the array vort as it appears in
c            the program that calls vrtgs. if isym = 0 then ivrt
c            must be at least nlat.  if isym = 1 or 2 and nlat is
c            even then ivrt must be at least nlat/2. if isym = 1 or 2
c            and nlat is odd then ivrt must be at least (nlat+1)/2.
c
c     jvrt   the second dimension of the array vort as it appears in
c            the program that calls vrtgs. jvrt must be at least nlon.
c
c    cr,ci   two or three dimensional arrays (see input parameter nt)
c            that contain vector spherical harmonic coefficients
c            of the vector field (v,w) as computed by subroutine vhags.
c     ***    cr and ci must be computed by vhags prior to calling
c            vrtgs.
c
c      mdc   the first dimension of the arrays cr and ci as it
c            appears in the program that calls vrtgs. mdc must be at
c            least min0(nlat,nlon/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c      ndc   the second dimension of the arrays cr and ci as it
c            appears in the program that calls vrtgs. ndc must be at
c            least nlat.
c
c   wshsgs   an array which must be initialized by subroutine shsgsi.
c            once initialized,
c            wshsgs can be used repeatedly by vrtgs as long as nlon
c            and nlat remain unchanged.  wshsgs must not be altered
c            between calls of vrtgs
c
c   lshsgs   the dimension of the array wshsgs   as it appears in the
c            program that calls vrtgs. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshsgs must be at least
c
c            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
c
c     work   a work array that does not have to be saved.
c
c    lwork   the dimension of the array work as it appears in the
c            program that calls vrtgs. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd.
c
c            if isym = 0 then lwork must be at least
c
c               nlat*((nt+1)*nlon+2*nt*l1+1)
c
c            if isym > 0 then lwork must be at least
c
c               (nt+1)*l2*nlon+nlat*(2*nt*l1+1)
c
c
c     **************************************************************
c
c     output parameters
c
c
c     vort   a two or three dimensional array (see input parameter nt)
c            that contains the vorticity of the vector field (v,w)
c            whose coefficients cr,ci where computed by subroutine vhags.
c            vort(i,j) is the vorticity at the gaussian colatitude point
c            theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon.
c            the index ranges are defined above at the input parameter
c            isym.
c
c
c   ierror   an error parameter which indicates fatal errors with input
c            parameters when returned positive.
c          = 0  no errors
c          = 1  error in the specification of nlat
c          = 2  error in the specification of nlon
c          = 3  error in the specification of isym
c          = 4  error in the specification of nt
c          = 5  error in the specification of ivrt
c          = 6  error in the specification of jvrt
c          = 7  error in the specification of mdc
c          = 8  error in the specification of ndc
c          = 9  error in the specification of lshsgs
c          = 10 error in the specification of lwork
c **********************************************************************
c                                                                              
c   
      subroutine vrtgs(nlat,nlon,isym,nt,vort,ivrt,jvrt,cr,ci,mdc,ndc,
     +                 wshsgs,lshsgs,work,lwork,ierror)

      dimension vort(ivrt,jvrt,nt),cr(mdc,ndc,nt),ci(mdc,ndc,nt)
      dimension wshsgs(lshsgs),work(lwork)
c
c     check input parameters
c
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 4) return
      ierror = 3
      if (isym.lt.0 .or. isym.gt.2) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      imid = (nlat+1)/2
      if((isym.eq.0 .and. ivrt.lt.nlat) .or.
     1   (isym.gt.0 .and. ivrt.lt.imid)) return
      ierror = 6
      if(jvrt .lt. nlon) return
      ierror = 7
      if(mdc .lt. min0(nlat,(nlon+1)/2)) return
      mmax = min0(nlat,(nlon+2)/2)
      ierror = 8
      if(ndc .lt. nlat) return
      ierror = 9
      imid = (nlat+1)/2
      lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2
      l2 = (nlat+mod(nlat,2))/2
      l1 = min0((nlon+2)/2,nlat)
      lp=nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
      if(lshsgs.lt.lp) return
      ierror = 10
c
c     verify unsaved work space (add to what shses requires, file f3)
c
c
c     set first dimension for a,b (as requried by shses)
c
      mab = min0(nlat,nlon/2+1)
      mn = mab*nlat*nt
      ls = nlat
      if(isym .gt. 0) ls = imid
      nln = nt*ls*nlon
      if(lwork.lt. nln+ls*nlon+2*mn+nlat) return
      ierror = 0
c
c     set work space pointers
c
      ia = 1
      ib = ia+mn
      is = ib+mn
      iwk = is+nlat
      lwk = lwork-2*mn-nlat
      call vrtgs1(nlat,nlon,isym,nt,vort,ivrt,jvrt,cr,ci,mdc,ndc,
     +work(ia),work(ib),mab,work(is),wshsgs,lshsgs,work(iwk),lwk,
     +ierror)
      return
      end

      subroutine vrtgs1(nlat,nlon,isym,nt,vort,ivrt,jvrt,cr,ci,mdc,ndc,
     +                  a,b,mab,sqnn,wsav,lwsav,wk,lwk,ierror)
      dimension vort(ivrt,jvrt,nt),cr(mdc,ndc,nt),ci(mdc,ndc,nt)
      dimension a(mab,nlat,nt),b(mab,nlat,nt),sqnn(nlat)
      dimension wsav(lwsav),wk(lwk)
c
c     set coefficient multiplyers
c
      do 1 n=2,nlat
      fn = float(n-1)
      sqnn(n) = sqrt(fn*(fn+1.))
    1 continue
c
c     compute divergence scalar coefficients for each vector field
c
      do 2 k=1,nt
      do 3 n=1,nlat
      do 4 m=1,mab
      a(m,n,k) = 0.0
      b(m,n,k) = 0.0
    4 continue
    3 continue
c
c     compute m=0 coefficients
c
      do 5 n=2,nlat
      a(1,n,k) = sqnn(n)*cr(1,n,k)
      b(1,n,k) = sqnn(n)*ci(1,n,k)
    5 continue
c
c     compute m>0 coefficients
c
      mmax = min0(nlat,(nlon+1)/2)
      do 6 m=2,mmax
      do 7 n=m,nlat
      a(m,n,k) = sqnn(n)*cr(m,n,k)
      b(m,n,k) = sqnn(n)*ci(m,n,k)
    7 continue
    6 continue
    2 continue
c
c     synthesize a,b into vort
c
      call shsgs(nlat,nlon,isym,nt,vort,ivrt,jvrt,a,b,
     +           mab,nlat,wsav,lwsav,wk,lwk,ierror)
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file shaes.f
c
c     this file contains code and documentation for subroutines
c     shaes and shaesi
c
c ... files which must be loaded with shaes.f
c
c     sphcom.f, hrfft.f
c
c     subroutine shaes(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
c    +                 wshaes,lshaes,work,lwork,ierror)
c
c     subroutine shaes performs the spherical harmonic analysis
c     on the array g and stores the result in the arrays a and b.
c     the analysis is performed on an equally spaced grid.  the
c     associated legendre functions are stored rather than recomputed
c     as they are in subroutine shaec.  the analysis is described
c     below at output parameters a,b.
c
c     sphcom.f, hrfft.f
c
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c     isym   = 0  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 array g(i,j) for i=1,...,nlat and j=1,...,nlon.
c                 (see description of g below)
c
c            = 1  g is antisymmetric about the equator. the analysis
c                 is performed on the northern hemisphere only.  i.e.
c                 if nlat is odd the analysis is performed on the
c                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
c                 if nlat is even the analysis is performed on the
c                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c
c            = 2  g is symmetric about the equator. the analysis is
c                 performed on the northern hemisphere only.  i.e.
c                 if nlat is odd the analysis is performed on the
c                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
c                 if nlat is even the analysis is performed on the
c                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c     nt     the number of analyses.  in the program that calls shaes,
c            the arrays g,a and b can be three dimensional in which
c            case multiple analyses will be performed.  the third
c            index is the analysis index which assumes the values
c            k=1,...,nt.  for a single analysis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that the arrays g,a and b
c            have only two dimensions.
c
c     g      a two or three dimensional array (see input parameter
c            nt) that contains the discrete function to be analyzed.
c            g(i,j) contains the value of the function at the colatitude
c            point theta(i) = (i-1)*pi/(nlat-1) and longitude point
c            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined
c            above at the input parameter isym.
c
c
c     idg    the first dimension of the array g as it appears in the
c            program that calls shaes.  if isym equals zero then idg
c            must be at least nlat.  if isym is nonzero then idg
c            must be at least nlat/2 if nlat is even or at least
c            (nlat+1)/2 if nlat is odd.
c
c     jdg    the second dimension of the array g as it appears in the
c            program that calls shaes.  jdg must be at least nlon.
c
c     mdab   the first dimension of the arrays a and b as it appears
c            in the program that calls shaes. mdab must be at least
c            min0(nlat,(nlon+2)/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     ndab   the second dimension of the arrays a and b as it appears
c            in the program that calls shaes. ndab must be at least nlat
c
c     wshaes an array which must be initialized by subroutine shaesi.
c            once initialized, wshaes can be used repeatedly by shaes
c            as long as nlon and nlat remain unchanged.  wshaes must
c            not be altered between calls of shaes.
c
c     lshaes the dimension of the array wshaes as it appears in the
c            program that calls shaes. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshaes must be at least
c
c               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls shaes.  define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            if isym is zero then lwork must be at least
c            (nt+1)*nlat*nlon. if isym is not zero then
c            lwork must be at least (nt+1)*l2*nlon.
c
c
c     **************************************************************
c
c     output parameters
c
c     a,b    both a,b are two or three dimensional arrays (see input
c            parameter nt) that contain the spherical harmonic
c            coefficients in the representation of g(i,j) given in the
c            discription of subroutine shses. for isym=0, a(m,n) and
c            b(m,n) are given by the equations listed below. symmetric
c            versions are used when isym is greater than zero.
c
c
c
c     definitions
c
c     1. the normalized associated legendre functions
c
c     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
c                       *sin(theta)**m/(2**n*factorial(n)) times the
c                       (n+m)th derivative of (x**2-1)**n with respect
c                       to x=cos(theta)
c
c     2. the normalized z functions for m even
c
c     zbar(m,n,theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of
c                       the integral from tau = 0 to tau = pi of
c                       cos(k*theta)*cos(k*tau)*pbar(m,n,tau)*sin(tau)
c                       (first and last terms in this sum are divided
c                       by 2)
c
c     3. the normalized z functions for m odd
c
c     zbar(m,n,theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of
c                       of the integral from tau = 0 to tau = pi of
c                       sin(k*theta)*sin(k*tau)*pbar(m,n,tau)*sin(tau)
c
c     4. the fourier transform of g(i,j).
c
c     c(m,i)          = 2/nlon times the sum from j=1 to j=nlon
c                       of g(i,j)*cos((m-1)*(j-1)*2*pi/nlon)
c                       (the first and last terms in this sum
c                       are divided by 2)
c
c     s(m,i)          = 2/nlon times the sum from j=2 to j=nlon
c                       of g(i,j)*sin((m-1)*(j-1)*2*pi/nlon)
c
c     5. the maximum (plus one) longitudinal wave number
c
c            mmax = min0(nlat,(nlon+2)/2) if nlon is even or
c            mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     then for m=0,...,mmax-1 and  n=m,...,nlat-1  the arrays a,b are
c     given by
c
c     a(m+1,n+1)      = the sum from i=1 to i=nlat of
c                       c(m+1,i)*zbar(m,n,theta(i))
c                       (first and last terms in this sum are
c                       divided by 2)
c
c     b(m+1,n+1)      = the sum from i=1 to i=nlat of
c                       s(m+1,i)*zbar(m,n,theta(i))
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of isym
c            = 4  error in the specification of nt
c            = 5  error in the specification of idg
c            = 6  error in the specification of jdg
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lshaes
c            = 10 error in the specification of lwork
c
c
c ****************************************************************
c     subroutine shaesi(nlat,nlon,wshaes,lshaes,work,lwork,dwork,
c    +                  ldwork,ierror)
c
c     subroutine shaesi initializes the array wshaes which can then
c     be used repeatedly by subroutine shaes
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c     lshaes the dimension of the array wshaes as it appears in the
c            program that calls shaesi. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshaes must be at least
c
c               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
c
c     work   a real   work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls shaesi.  define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lwork must be at least
c
c               5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
c
c
c     dwork  a double precision work array that does not have to be saved.
c
c     ldwork the dimension of the array dwork as it appears in the
c            program that calls shaesi.  ldwork must be at least nlat+1
c
c
c     output parameters
c
c     wshaes an array which is initialized for use by subroutine shaes.
c            once initialized, wshaes can be used repeatedly by shaes
c            as long as nlon and nlat remain unchanged.  wshaes must
c            not be altered between calls of shaes.
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lshaes
c            = 4  error in the specification of lwork
c            = 5  error in the specification of ldwork
c
c
c ****************************************************************
      subroutine shaes(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
     1                    wshaes,lshaes,work,lwork,ierror)
      dimension g(idg,jdg,1),a(mdab,ndab,1),b(mdab,ndab,1),wshaes(1),
     1          work(1)
      ierror = 1
      if(nlat.lt.3) return
      ierror = 2
      if(nlon.lt.4) return
      ierror = 3
      if(isym.lt.0 .or. isym.gt.2) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      if((isym.eq.0 .and. idg.lt.nlat) .or.
     1   (isym.ne.0 .and. idg.lt.(nlat+1)/2)) return
      ierror = 6
      if(jdg .lt. nlon) return
      ierror = 7
      mmax = min0(nlat,nlon/2+1)
      if(mdab .lt. mmax) return
      ierror = 8
      if(ndab .lt. nlat) return
      ierror = 9
      imid = (nlat+1)/2
      idz = (mmax*(nlat+nlat-mmax+1))/2
      lzimn = idz*imid
      if(lshaes .lt. lzimn+nlon+15) return
      ierror = 10
      ls = nlat
      if(isym .gt. 0) ls = imid
      nln = nt*ls*nlon
      if(lwork .lt. nln+ls*nlon) return
      ierror = 0
      ist = 0
      if(isym .eq. 0) ist = imid
      call shaes1(nlat,isym,nt,g,idg,jdg,a,b,mdab,ndab,wshaes,idz,
     1         ls,nlon,work,work(ist+1),work(nln+1),wshaes(lzimn+1))
      return
      end
      subroutine shaes1(nlat,isym,nt,g,idgs,jdgs,a,b,mdab,ndab,z,idz,
     1                  idg,jdg,ge,go,work,whrfft)
      dimension g(idgs,jdgs,1),a(mdab,ndab,1),b(mdab,ndab,1),z(idz,1),
     1          ge(idg,jdg,1),go(idg,jdg,1),work(1),whrfft(1)
      ls = idg
      nlon = jdg
      mmax = min0(nlat,nlon/2+1)
      mdo = mmax
      if(mdo+mdo-1 .gt. nlon) mdo = mmax-1
      nlp1 = nlat+1
      tsn = 2./nlon
      fsn = 4./nlon
      imid = (nlat+1)/2
      modl = mod(nlat,2)
      imm1 = imid
      if(modl .ne. 0) imm1 = imid-1
      if(isym .ne. 0) go to 15
      do 5 k=1,nt
      do 5 i=1,imm1
      do 5 j=1,nlon
      ge(i,j,k) = tsn*(g(i,j,k)+g(nlp1-i,j,k))
      go(i,j,k) = tsn*(g(i,j,k)-g(nlp1-i,j,k))
    5 continue
      go to 30
   15 do 20 k=1,nt
      do 20 i=1,imm1
      do 20 j=1,nlon
      ge(i,j,k) = fsn*g(i,j,k)
   20 continue
      if(isym .eq. 1) go to 27
   30 if(modl .eq. 0) go to 27
      do 25 k=1,nt
      do 25 j=1,nlon
      ge(imid,j,k) = tsn*g(imid,j,k)
   25 continue
   27 do 35 k=1,nt
      call hrfftf(ls,nlon,ge(1,1,k),ls,whrfft,work)
      if(mod(nlon,2) .ne. 0) go to 35
      do 36 i=1,ls
      ge(i,nlon,k) = .5*ge(i,nlon,k)
   36 continue
   35 continue
      do 40 k=1,nt
      do 40 mp1=1,mmax
      do 40 np1=mp1,nlat
      a(mp1,np1,k) = 0.
      b(mp1,np1,k) = 0.
   40 continue
      if(isym .eq. 1) go to 145
      do 110 k=1,nt
      do 110 i=1,imid
      do 110 np1=1,nlat,2
      a(1,np1,k) = a(1,np1,k)+z(np1,i)*ge(i,1,k)
  110 continue
      ndo = nlat
      if(mod(nlat,2) .eq. 0) ndo = nlat-1
      do 120 mp1=2,mdo
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      do 120 k=1,nt
      do 120 i=1,imid
      do 120 np1=mp1,ndo,2
      a(mp1,np1,k) = a(mp1,np1,k)+z(np1+mb,i)*ge(i,2*mp1-2,k)
      b(mp1,np1,k) = b(mp1,np1,k)+z(np1+mb,i)*ge(i,2*mp1-1,k)
  120 continue
      if(mdo .eq. mmax .or. mmax .gt. ndo) go to 135
      mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
      do 130 k=1,nt
      do 130 i=1,imid
      do 130 np1=mmax,ndo,2
      a(mmax,np1,k) = a(mmax,np1,k)+z(np1+mb,i)*ge(i,2*mmax-2,k)
  130 continue
  135 if(isym .eq. 2) return
  145 do 150 k=1,nt
      do 150 i=1,imm1
      do 150 np1=2,nlat,2
      a(1,np1,k) = a(1,np1,k)+z(np1,i)*go(i,1,k)
  150 continue
      ndo = nlat
      if(mod(nlat,2) .ne. 0) ndo = nlat-1
      do 160 mp1=2,mdo
      m = mp1-1
      mp2 = mp1+1
      mb = m*(nlat-1)-(m*(m-1))/2
      do 160 k=1,nt
      do 160 i=1,imm1
      do 160 np1=mp2,ndo,2
      a(mp1,np1,k) = a(mp1,np1,k)+z(np1+mb,i)*go(i,2*mp1-2,k)
      b(mp1,np1,k) = b(mp1,np1,k)+z(np1+mb,i)*go(i,2*mp1-1,k)
  160 continue
      mp2 = mmax+1
      if(mdo .eq. mmax .or. mp2 .gt. ndo) return
      mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
      do 170 k=1,nt
      do 170 i=1,imm1
      do 170 np1=mp2,ndo,2
      a(mmax,np1,k) = a(mmax,np1,k)+z(np1+mb,i)*go(i,2*mmax-2,k)
  170 continue
      return
      end

      subroutine shaesi(nlat,nlon,wshaes,lshaes,work,lwork,dwork,
     +                  ldwork,ierror)
      dimension wshaes(*),work(*)
      double precision dwork(*)
c
c     length of wshaes is (l*(l+1)*imid)/2+nlon+15
c     length of work is 5*l*imid + 3*((l-3)*l+2)/2
c
      ierror = 1
      if(nlat.lt.3) return
      ierror = 2
      if(nlon.lt.4) return
      ierror = 3
      mmax = min0(nlat,nlon/2+1)
      imid = (nlat+1)/2
      lzimn = (imid*mmax*(nlat+nlat-mmax+1))/2
      if(lshaes .lt. lzimn+nlon+15) return
      ierror = 4
      labc = 3*((mmax-2)*(nlat+nlat-mmax-1))/2
      if(lwork .lt. 5*nlat*imid + labc) return
      ierror = 5
      if (ldwork .lt. nlat+1) return
      ierror = 0
      iw1 = 3*nlat*imid+1
      idz = (mmax*(nlat+nlat-mmax+1))/2
      call sea1(nlat,nlon,imid,wshaes,idz,work,work(iw1),dwork)
      call hrffti(nlon,wshaes(lzimn+1))
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                       .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file shses.f
c
c     this file contains code and documentation for subroutines
c     shses and shsesi
c
c ... files which must be loaded with shses.f
c
c     sphcom.f, hrfft.f
c
c     subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
c    +                 wshses,lshses,work,lwork,ierror)
c
c     subroutine shses performs the spherical harmonic synthesis
c     on the arrays a and b and stores the result in the array g.
c     the synthesis is performed on an equally spaced grid.  the
c     associated legendre functions are stored rather than recomputed
c     as they are in subroutine shsec.  the synthesis is described
c     below at output parameter g.
c
c *** required files from spherepack2
c
c     sphcom.f, hrfft.f
c
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c     isym   = 0  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 array g(i,j) for i=1,...,nlat and j=1,...,nlon.
c                 (see description of g below)
c
c            = 1  g is antisymmetric about the equator. the synthesis
c                 is performed on the northern hemisphere only.  i.e.
c                 if nlat is odd the synthesis is performed on the
c                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
c                 if nlat is even the synthesis is performed on the
c                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c
c            = 2  g is symmetric about the equator. the synthesis is
c                 performed on the northern hemisphere only.  i.e.
c                 if nlat is odd the synthesis is performed on the
c                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
c                 if nlat is even the synthesis is performed on the
c                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c     nt     the number of syntheses.  in the program that calls shses,
c            the arrays g,a and b can be three dimensional in which
c            case multiple syntheses will be performed.  the third
c            index is the synthesis index which assumes the values
c            k=1,...,nt.  for a single synthesis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that the arrays g,a and b
c            have only two dimensions.
c
c     idg    the first dimension of the array g as it appears in the
c            program that calls shses.  if isym equals zero then idg
c            must be at least nlat.  if isym is nonzero then idg
c            must be at least nlat/2 if nlat is even or at least
c            (nlat+1)/2 if nlat is odd.
c
c     jdg    the second dimension of the array g as it appears in the
c            program that calls shses.  jdg must be at least nlon.
c
c     a,b    two or three dimensional arrays (see the input parameter
c            nt) that contain the coefficients in the spherical harmonic
c            expansion of g(i,j) given below at the definition of the
c            output parameter g.  a(m,n) and b(m,n) are defined for
c            indices m=1,...,mmax and n=m,...,nlat where mmax is the
c            maximum (plus one) longitudinal wave number given by
c            mmax = min0(nlat,(nlon+2)/2) if nlon is even or
c            mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     mdab   the first dimension of the arrays a and b as it appears
c            in the program that calls shses. mdab must be at least
c            min0(nlat,(nlon+2)/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     ndab   the second dimension of the arrays a and b as it appears
c            in the program that calls shses. ndab must be at least nlat
c
c     wshses an array which must be initialized by subroutine shsesi.
c            once initialized, wshses can be used repeatedly by shses
c            as long as nlon and nlat remain unchanged.  wshses must
c            not be altered between calls of shses.
c
c     lshses the dimension of the array wshses as it appears in the
c            program that calls shses. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshses must be at least
c
c               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls shses.  define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            if isym is zero then lwork must be at least
c
c               (nt+1)*nlat*nlon
c
c            if isym is nonzero lwork must be at least
c
c               (nt+1)*l2*nlon.
c
c     **************************************************************
c
c     output parameters
c
c     g      a two or three dimensional array (see input parameter
c            nt) that contains the spherical harmonic synthesis of
c            the arrays a and b at the colatitude point theta(i) =
c            (i-1)*pi/(nlat-1) and longitude point phi(j) =
c            (j-1)*2*pi/nlon. the index ranges are defined above at
c            at the input parameter isym.  for isym=0, g(i,j) is
c            given by the the equations listed below.  symmetric
c            versions are used when isym is greater than zero.
c
c     the normalized associated legendre functions are given by
c
c     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
c                       *sin(theta)**m/(2**n*factorial(n)) times the
c                       (n+m)th derivative of (x**2-1)**n with respect
c                       to x=cos(theta)
c
c     define the maximum (plus one) longitudinal wave number
c     as   mmax = min0(nlat,(nlon+2)/2) if nlon is even or
c          mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     then g(i,j) = the sum from n=0 to n=nlat-1 of
c
c                   .5*pbar(0,n,theta(i))*a(1,n+1)
c
c              plus the sum from m=1 to m=mmax-1 of
c
c                   the sum from n=m to n=nlat-1 of
c
c              pbar(m,n,theta(i))*(a(m+1,n+1)*cos(m*phi(j))
c                                    -b(m+1,n+1)*sin(m*phi(j)))
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of isym
c            = 4  error in the specification of nt
c            = 5  error in the specification of idg
c            = 6  error in the specification of jdg
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lshses
c            = 10 error in the specification of lwork
c
c
c ****************************************************************
c     subroutine shsesi(nlat,nlon,wshses,lshses,work,lwork,dwork,
c    +                  ldwork,ierror)
c
c     subroutine shsesi initializes the array wshses which can then
c     be used repeatedly by subroutine shses.
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c     lshses the dimension of the array wshses as it appears in the
c            program that calls shsesi. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshses must be at least
c
c               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
c
c     work   a real   work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in
c            the program that calls shsesi.  define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lwork must be at least
c
c               5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
c
c
c     dwork  a double precision work array that does not have to be saved.
c
c     ldwork the dimension of the array dwork as it appears in the
c            program that calls shsesi.  ldwork must be at least nlat+1
c
c
c     output parameters
c
c     wshses an array which is initialized for use by subroutine shses.
c            once initialized, wshses can be used repeatedly by shses
c            as long as nlon and nlat remain unchanged.  wshses must
c            not be altered between calls of shses.
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lshses
c            = 4  error in the specification of lwork
c            = 5  error in the specification of ldwork
c
c ****************************************************************
      subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
     1                    wshses,lshses,work,lwork,ierror)
      dimension g(idg,jdg,1),a(mdab,ndab,1),b(mdab,ndab,1),wshses(1),
     1          work(1)
      ierror = 1
      if(nlat.lt.3) return
      ierror = 2
      if(nlon.lt.4) return
      ierror = 3
      if(isym.lt.0 .or. isym.gt.2) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      if((isym.eq.0 .and. idg.lt.nlat) .or.
     1   (isym.ne.0 .and. idg.lt.(nlat+1)/2)) return
      ierror = 6
      if(jdg .lt. nlon) return
      ierror = 7
      mmax = min0(nlat,nlon/2+1)
      if(mdab .lt. mmax) return
      ierror = 8
      if(ndab .lt. nlat) return
      ierror = 9
      imid = (nlat+1)/2
      lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2
      if(lshses .lt. lpimn+nlon+15) return
      ierror = 10
      ls = nlat
      if(isym .gt. 0) ls = imid
      nln = nt*ls*nlon
      if(lwork.lt. nln+ls*nlon) return
      ierror = 0
      ist = 0
      if(isym .eq. 0) ist = imid
      call shses1(nlat,isym,nt,g,idg,jdg,a,b,mdab,ndab,wshses,imid,
     1       ls,nlon,work,work(ist+1),work(nln+1),wshses(lpimn+1))
      return
      end
      subroutine shses1(nlat,isym,nt,g,idgs,jdgs,a,b,mdab,ndab,p,imid,
     1                  idg,jdg,ge,go,work,whrfft)
      dimension g(idgs,jdgs,1),a(mdab,ndab,1),b(mdab,ndab,1),p(imid,1),
     1          ge(idg,jdg,1),go(idg,jdg,1),work(1),whrfft(1)
      ls = idg
      nlon = jdg
      mmax = min0(nlat,nlon/2+1)
      mdo = mmax
      if(mdo+mdo-1 .gt. nlon) mdo = mmax-1
      nlp1 = nlat+1
      modl = mod(nlat,2)
      imm1 = imid
      if(modl .ne. 0) imm1 = imid-1
      do 80 k=1,nt
      do 80 j=1,nlon
      do 80 i=1,ls
      ge(i,j,k) = 0.
 8000 continue
  800 continue
   80 continue
      if(isym .eq. 1) go to 125
      do 100 k=1,nt
      do 100 np1=1,nlat,2
      do 100 i=1,imid
      ge(i,1,k)=ge(i,1,k)+a(1,np1,k)*p(i,np1)
  100 continue
      ndo = nlat
      if(mod(nlat,2) .eq. 0) ndo = nlat-1
      do 110 mp1=2,mdo
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      do 110 np1=mp1,ndo,2
      mn = mb+np1
      do 110 k=1,nt
      do 110 i=1,imid
      ge(i,2*mp1-2,k) = ge(i,2*mp1-2,k)+a(mp1,np1,k)*p(i,mn)
      ge(i,2*mp1-1,k) = ge(i,2*mp1-1,k)+b(mp1,np1,k)*p(i,mn)
  110 continue
      if(mdo .eq. mmax .or. mmax .gt. ndo) go to 122
      mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
      do 120 np1=mmax,ndo,2
      mn = mb+np1
      do 120 k=1,nt
      do 120 i=1,imid
      ge(i,2*mmax-2,k) = ge(i,2*mmax-2,k)+a(mmax,np1,k)*p(i,mn)
  120 continue
  122 if(isym .eq. 2) go to 155
  125 do 140 k=1,nt
      do 140 np1=2,nlat,2
      do 140 i=1,imm1
      go(i,1,k)=go(i,1,k)+a(1,np1,k)*p(i,np1)
  140 continue
      ndo = nlat
      if(mod(nlat,2) .ne. 0) ndo = nlat-1
      do 150 mp1=2,mdo
      mp2 = mp1+1
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      do 150 np1=mp2,ndo,2
      mn = mb+np1
      do 150 k=1,nt
      do 150 i=1,imm1
      go(i,2*mp1-2,k) = go(i,2*mp1-2,k)+a(mp1,np1,k)*p(i,mn)
      go(i,2*mp1-1,k) = go(i,2*mp1-1,k)+b(mp1,np1,k)*p(i,mn)
  150 continue
      mp2 = mmax+1
      if(mdo .eq. mmax .or. mp2 .gt. ndo) go to 155
      mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
      do 152 np1=mp2,ndo,2
      mn = mb+np1
      do 152 k=1,nt
      do 152 i=1,imm1
      go(i,2*mmax-2,k) = go(i,2*mmax-2,k)+a(mmax,np1,k)*p(i,mn)
  152 continue
  155 do 160 k=1,nt
      if(mod(nlon,2) .ne. 0) go to 157
      do 156 i=1,ls
      ge(i,nlon,k) = 2.*ge(i,nlon,k)
  156 continue
  157 call hrfftb(ls,nlon,ge(1,1,k),ls,whrfft,work)
  160 continue
      if(isym .ne. 0) go to 180
      do 170 k=1,nt
      do 170 j=1,nlon
      do 175 i=1,imm1
      g(i,j,k) = .5*(ge(i,j,k)+go(i,j,k))
      g(nlp1-i,j,k) = .5*(ge(i,j,k)-go(i,j,k))
  175 continue
      if(modl .eq. 0) go to 170
      g(imid,j,k) = .5*ge(imid,j,k)
  170 continue
      return
  180 do 185 k=1,nt
      do 185 i=1,imid
      do 185 j=1,nlon
      g(i,j,k) = .5*ge(i,j,k)
  185 continue
      return
      end

      subroutine shsesi(nlat,nlon,wshses,lshses,work,lwork,dwork,
     +                  ldwork,ierror)
      dimension wshses(*),work(*)
      double precision dwork(*)
      ierror = 1
      if(nlat.lt.3) return
      ierror = 2
      if(nlon.lt.4) return
      ierror = 3
      mmax = min0(nlat,nlon/2+1)
      imid = (nlat+1)/2
      lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2
      if(lshses .lt. lpimn+nlon+15) return
      ierror = 4
      labc = 3*((mmax-2)*(nlat+nlat-mmax-1))/2
      if(lwork .lt. 5*nlat*imid + labc) return
      ierror = 5
      if (ldwork .lt. nlat+1) return
      ierror = 0
      iw1 = 3*nlat*imid+1
      CALL SES1(NLAT,NLON,IMID,WSHSES,WORK,WORK(IW1),DWORK)
      call hrffti(nlon,wshses(lpimn+1))
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file vhaes.f
c
c     this file contains code and documentation for subroutines
c     vhaes and vhaesi
c
c ... files which must be loaded with vhaes.f
c
c     sphcom.f, hrfft.f
c
c                                                                              
c     subroutine vhaes(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
c    +                 mdab,ndab,wvhaes,lvhaes,work,lwork,ierror)
c
c     subroutine vhaes performs the vector spherical harmonic analysis
c     on the vector field (v,w) and stores the result in the arrays
c     br, bi, cr, and ci. v(i,j) and w(i,j) are the colatitudinal 
c     (measured from the north pole) and east longitudinal components
c     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
c     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
c     representation of (v,w) is given at output parameters v,w in 
c     subroutine vhses.  
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     ityp   = 0  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon.   
c
c            = 1  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 2  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 3  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 4  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 5  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 6  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 7  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 8  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c
c     nt     the number of analyses.  in the program that calls vhaes,
c            the arrays v,w,br,bi,cr, and ci can be three dimensional
c            in which case multiple analyses will be performed.
c            the third index is the analysis index which assumes the 
c            values k=1,...,nt.  for a single analysis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that all the arrays are two
c            dimensional.
c
c     v,w    two or three dimensional arrays (see input parameter nt)
c            that contain the vector function to be analyzed.
c            v is the colatitudnal component and w is the east 
c            longitudinal component. v(i,j),w(i,j) contain the
c            components at colatitude theta(i) = (i-1)*pi/(nlat-1)
c            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
c            are defined above at the input parameter ityp.
c
c     idvw   the first dimension of the arrays v,w as it appears in
c            the program that calls vhaes. if ityp .le. 2 then idvw
c            must be at least nlat.  if ityp .gt. 2 and nlat is
c            even then idvw must be at least nlat/2. if ityp .gt. 2
c            and nlat is odd then idvw must be at least (nlat+1)/2.
c
c     jdvw   the second dimension of the arrays v,w as it appears in
c            the program that calls vhaes. jdvw must be at least nlon.
c
c     mdab   the first dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhaes. mdab must be at
c            least min0(nlat,nlon/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     ndab   the second dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhaes. ndab must be at
c            least nlat.
c
c     lvhaes an array which must be initialized by subroutine vhaesi.
c            once initialized, wvhaes can be used repeatedly by vhaes
c            as long as nlon and nlat remain unchanged.  wvhaes must
c            not be altered between calls of vhaes.
c
c     lvhaes the dimension of the array wvhaes as it appears in the
c            program that calls vhaes. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhaes must be at least
c
c            l1*l2(nlat+nlat-l1+1)+nlon+15
c
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhaes. define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            if ityp .le. 2 then lwork must be at least
c
c                       (2*nt+1)*nlat*nlon
c
c            if ityp .gt. 2 then lwork must be at least
c
c                        (2*nt+1)*l2*nlon   
c
c     **************************************************************
c
c     output parameters
c
c     br,bi  two or three dimensional arrays (see input parameter nt)
c     cr,ci  that contain the vector spherical harmonic coefficients
c            in the spectral representation of v(i,j) and w(i,j) given 
c            in the discription of subroutine vhses. br(mp1,np1),
c            bi(mp1,np1),cr(mp1,np1), and ci(mp1,np1) are computed 
c            for mp1=1,...,mmax and np1=mp1,...,nlat except for np1=nlat
c            and odd mp1. mmax=min0(nlat,nlon/2) if nlon is even or 
c            mmax=min0(nlat,(nlon+1)/2) if nlon is odd. 
c      
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of ityp
c            = 4  error in the specification of nt
c            = 5  error in the specification of idvw
c            = 6  error in the specification of jdvw
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lvhaes
c            = 10 error in the specification of lwork
c
c ********************************************************
c
c     subroutine vhaesi(nlat,nlon,wvhaes,lvhaes,work,lwork,dwork,
c    +                  ldwork,ierror)
c
c     subroutine vhaesi initializes the array wvhaes which can then be
c     used repeatedly by subroutine vhaes until nlat or nlon is changed.
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     lvhaes the dimension of the array wvhaes as it appears in the
c            program that calls vhaes. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhaes must be at least
c
c               l1*l2*(nlat+nlat-l1+1)+nlon+15
c
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhaes. lwork must be at least
c
c              3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+5*l2*nlat
c
c     dwork  an unsaved double precision work space
c
c     ldwork the length of the array dwork as it appears in the
c            program that calls vhaesi.  ldwork must be at least
c            2*(nlat+1)
c
c
c     **************************************************************
c
c     output parameters
c
c     wvhaes an array which is initialized for use by subroutine vhaes.
c            once initialized, wvhaes can be used repeatedly by vhaes
c            as long as nlat or nlon remain unchanged.  wvhaes must not
c            be altered between calls of vhaes.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lvhaes
c            = 4  error in the specification of lwork
c            = 5  error in the specification of ldwork
c
c
      subroutine vhaes(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
     1           mdab,ndab,wvhaes,lvhaes,work,lwork,ierror)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          work(1),wvhaes(1)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      if(ityp.lt.0 .or. ityp.gt.8) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      imid = (nlat+1)/2
      if((ityp.le.2 .and. idvw.lt.nlat) .or.
     1   (ityp.gt.2 .and. idvw.lt.imid)) return
      ierror = 6
      if(jdvw .lt. nlon) return
      ierror = 7
      mmax = min0(nlat,(nlon+1)/2)
      if(mdab .lt. mmax) return
      ierror = 8
      if(ndab .lt. nlat) return
      ierror = 9
      idz = (mmax*(nlat+nlat-mmax+1))/2
      lzimn = idz*imid
      if(lvhaes .lt. lzimn+lzimn+nlon+15) return
      ierror = 10
      idv = nlat
      if(ityp .gt. 2) idv = imid
      lnl = nt*idv*nlon
      if(lwork .lt. lnl+lnl+idv*nlon) return
      ierror = 0
      ist = 0
      if(ityp .le. 2) ist = imid
      iw1 = ist+1
      iw2 = lnl+1
      iw3 = iw2+ist
      iw4 = iw2+lnl
      jw1 = lzimn+1
      jw2 = jw1+lzimn
      call vhaes1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,ndab,
     1     br,bi,cr,ci,idv,work,work(iw1),work(iw2),work(iw3),
     2     work(iw4),idz,wvhaes,wvhaes(jw1),wvhaes(jw2))
      return
      end

      subroutine vhaes1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,
     1   ndab,br,bi,cr,ci,idv,ve,vo,we,wo,work,idz,zv,zw,wrfft)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          ve(idv,nlon,1),vo(idv,nlon,1),we(idv,nlon,1),
     3          wo(idv,nlon,1),work(1),wrfft(1),
     4          zv(idz,1),zw(idz,1)
      nlp1 = nlat+1
      tsn = 2./nlon
      fsn = 4./nlon
      mlat = mod(nlat,2)
      mlon = mod(nlon,2)
      mmax = min0(nlat,(nlon+1)/2)
      imm1 = imid
      if(mlat .ne. 0) imm1 = imid-1
      if(ityp .gt. 2) go to 3  
      do 5 k=1,nt 
      do 5 i=1,imm1
      do 5 j=1,nlon
      ve(i,j,k) = tsn*(v(i,j,k)+v(nlp1-i,j,k))
      vo(i,j,k) = tsn*(v(i,j,k)-v(nlp1-i,j,k))
      we(i,j,k) = tsn*(w(i,j,k)+w(nlp1-i,j,k))
      wo(i,j,k) = tsn*(w(i,j,k)-w(nlp1-i,j,k))
    5 continue
      go to 2
    3 do 8 k=1,nt
      do 8 i=1,imm1 
      do 8 j=1,nlon
      ve(i,j,k) = fsn*v(i,j,k)
      vo(i,j,k) = fsn*v(i,j,k)
      we(i,j,k) = fsn*w(i,j,k)
      wo(i,j,k) = fsn*w(i,j,k)
    8 continue
    2 if(mlat .eq. 0) go to 7
      do 6 k=1,nt 
      do 6 j=1,nlon
      ve(imid,j,k) = tsn*v(imid,j,k)
      we(imid,j,k) = tsn*w(imid,j,k)
    6 continue
    7 do 9 k=1,nt
      call hrfftf(idv,nlon,ve(1,1,k),idv,wrfft,work)
      call hrfftf(idv,nlon,we(1,1,k),idv,wrfft,work)
    9 continue 
      ndo1 = nlat
      ndo2 = nlat
      if(mlat .ne. 0) ndo1 = nlat-1
      if(mlat .eq. 0) ndo2 = nlat-1
      if(ityp.eq.2 .or. ityp.eq.5 .or. ityp.eq.8) go to 11 
      do 10 k=1,nt
      do 10 mp1=1,mmax
      do 10 np1=mp1,nlat
      br(mp1,np1,k)=0.
      bi(mp1,np1,k)=0.
   10 continue
   11 if(ityp.eq.1 .or. ityp.eq.4 .or. ityp.eq.7) go to 13 
      do 12 k=1,nt
      do 12 mp1=1,mmax
      do 12 np1=mp1,nlat
      cr(mp1,np1,k)=0.
      ci(mp1,np1,k)=0.
   12 continue
   13 itypp = ityp+1
      go to (1,100,200,300,400,500,600,700,800),itypp
c
c     case ityp=0 ,  no symmetries
c
c     case m=0
c
    1 do 15 k=1,nt
      do 15 i=1,imid
      do 15 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*ve(i,1,k)
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*we(i,1,k)
   15 continue
      do 16 k=1,nt
      do 16 i=1,imm1
      do 16 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*vo(i,1,k)
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*wo(i,1,k)
   16 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 20 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 17
      do 23 k=1,nt
      do 23 i=1,imm1
      do 23 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*we(i,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*ve(i,2*mp1-2,k)
   23 continue
      if(mlat .eq. 0) go to 17
      do 24 k=1,nt
      do 24 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zw(np1+mb,imid)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-zw(np1+mb,imid)*we(imid,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)+zw(np1+mb,imid)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zw(np1+mb,imid)*ve(imid,2*mp1-2,k)
   24 continue
   17 if(mp2 .gt. ndo2) go to 20
      do 21 k=1,nt
      do 21 i=1,imm1
      do 21 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*wo(i,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*vo(i,2*mp1-2,k)
   21 continue
      if(mlat .eq. 0) go to 20
      do 22 k=1,nt
      do 22 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-1,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-1,k)
   22 continue
   20 continue
      return
c
c     case ityp=1 ,  no symmetries but cr and ci equal zero
c
c     case m=0
c
  100 do 115 k=1,nt
      do 115 i=1,imid
      do 115 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*ve(i,1,k)
  115 continue
      do 116 k=1,nt
      do 116 i=1,imm1
      do 116 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*vo(i,1,k)
  116 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 120 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 117
      do 123 k=1,nt
      do 123 i=1,imm1
      do 123 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*we(i,2*mp1-2,k)
  123 continue
      if(mlat .eq. 0) go to 117
      do 124 k=1,nt
      do 124 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zw(np1+mb,imid)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-zw(np1+mb,imid)*we(imid,2*mp1-2,k)
  124 continue
  117 if(mp2 .gt. ndo2) go to 120
      do 121 k=1,nt
      do 121 i=1,imm1
      do 121 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*wo(i,2*mp1-2,k)
  121 continue
      if(mlat .eq. 0) go to 120
      do 122 k=1,nt
      do 122 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-1,k)
  122 continue
  120 continue
      return
c
c     case ityp=2 ,  no symmetries but br and bi equal zero   
c
c     case m=0
c
  200 do 215 k=1,nt
      do 215 i=1,imid
      do 215 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*we(i,1,k)
  215 continue
      do 216 k=1,nt
      do 216 i=1,imm1
      do 216 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*wo(i,1,k)
  216 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 220 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 217
      do 223 k=1,nt
      do 223 i=1,imm1
      do 223 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*ve(i,2*mp1-2,k)
  223 continue
      if(mlat .eq. 0) go to 217
      do 224 k=1,nt
      do 224 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+zw(np1+mb,imid)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zw(np1+mb,imid)*ve(imid,2*mp1-2,k)
  224 continue
  217 if(mp2 .gt. ndo2) go to 220
      do 221 k=1,nt
      do 221 i=1,imm1
      do 221 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*vo(i,2*mp1-2,k)
  221 continue
      if(mlat .eq. 0) go to 220
      do 222 k=1,nt
      do 222 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-1,k)
  222 continue
  220 continue
      return
c
c     case ityp=3 ,  v even , w odd
c
c     case m=0
c
  300 do 315 k=1,nt
      do 315 i=1,imid
      do 315 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*ve(i,1,k)
  315 continue
      do 316 k=1,nt
      do 316 i=1,imm1
      do 316 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*wo(i,1,k)
  316 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 320 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 317
      do 323 k=1,nt
      do 323 i=1,imm1
      do 323 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*ve(i,2*mp1-2,k)
  323 continue
      if(mlat .eq. 0) go to 317
      do 324 k=1,nt
      do 324 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+zw(np1+mb,imid)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zw(np1+mb,imid)*ve(imid,2*mp1-2,k)
  324 continue
  317 if(mp2 .gt. ndo2) go to 320
      do 321 k=1,nt
      do 321 i=1,imm1
      do 321 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*wo(i,2*mp1-2,k)
  321 continue
      if(mlat .eq. 0) go to 320
      do 322 k=1,nt
      do 322 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-1,k)
  322 continue
  320 continue
      return
c
c     case ityp=4 ,  v even, w odd, and cr and ci equal 0. 
c
c     case m=0
c
  400 do 415 k=1,nt
      do 415 i=1,imid
      do 415 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*ve(i,1,k)
  415 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 420 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 420
      do 421 k=1,nt
      do 421 i=1,imm1
      do 421 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*wo(i,2*mp1-2,k)
  421 continue
      if(mlat .eq. 0) go to 420
      do 422 k=1,nt
      do 422 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-1,k)
  422 continue
  420 continue
      return
c
c     case ityp=5   v even, w odd, and br and bi equal zero
c
c     case m=0
c
  500 do 516 k=1,nt
      do 516 i=1,imm1
      do 516 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*wo(i,1,k)
  516 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 520 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 520
      do 523 k=1,nt
      do 523 i=1,imm1
      do 523 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*ve(i,2*mp1-2,k)
  523 continue
      if(mlat .eq. 0) go to 520
      do 524 k=1,nt
      do 524 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+zw(np1+mb,imid)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zw(np1+mb,imid)*ve(imid,2*mp1-2,k)
  524 continue
  520 continue
      return
c
c     case ityp=6 ,  v odd , w even
c
c     case m=0
c
  600 do 615 k=1,nt
      do 615 i=1,imid
      do 615 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*we(i,1,k)
  615 continue
      do 616 k=1,nt
      do 616 i=1,imm1
      do 616 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*vo(i,1,k)
  616 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 620 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 617
      do 623 k=1,nt
      do 623 i=1,imm1
      do 623 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*we(i,2*mp1-2,k)
  623 continue
      if(mlat .eq. 0) go to 617
      do 624 k=1,nt
      do 624 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zw(np1+mb,imid)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-zw(np1+mb,imid)*we(imid,2*mp1-2,k)
  624 continue
  617 if(mp2 .gt. ndo2) go to 620
      do 621 k=1,nt
      do 621 i=1,imm1
      do 621 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*vo(i,2*mp1-2,k)
  621 continue
      if(mlat .eq. 0) go to 620
      do 622 k=1,nt
      do 622 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-1,k)
  622 continue
  620 continue
      return
c
c     case ityp=7   v odd, w even, and cr and ci equal zero
c
c     case m=0
c
  700 do 716 k=1,nt
      do 716 i=1,imm1
      do 716 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*vo(i,1,k)
  716 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 720 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 720
      do 723 k=1,nt
      do 723 i=1,imm1
      do 723 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*we(i,2*mp1-2,k)
  723 continue
      if(mlat .eq. 0) go to 720
      do 724 k=1,nt
      do 724 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zw(np1+mb,imid)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-zw(np1+mb,imid)*we(imid,2*mp1-2,k)
  724 continue
  720 continue
      return
c
c     case ityp=8   v odd, w even, and both br and bi equal zero
c
c     case m=0
c
  800 do 815 k=1,nt
      do 815 i=1,imid
      do 815 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*we(i,1,k)
  815 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 820 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 820
      do 821 k=1,nt
      do 821 i=1,imm1
      do 821 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*vo(i,2*mp1-2,k)
  821 continue
      if(mlat .eq. 0) go to 820
      do 822 k=1,nt
      do 822 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-1,k)
  822 continue
  820 continue
      return
      end
c
c     dwork must be of length at least 2*(nlat+1)
c
      subroutine vhaesi(nlat,nlon,wvhaes,lvhaes,work,lwork,dwork,
     +                  ldwork,ierror)
      dimension wvhaes(lvhaes),work(lwork)
      double precision dwork(ldwork)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      mmax = min0(nlat,(nlon+1)/2)
      imid = (nlat+1)/2
      lzimn = (imid*mmax*(nlat+nlat-mmax+1))/2
      if(lvhaes .lt. lzimn+lzimn+nlon+15) return
      ierror = 4
      labc = 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      if(lwork .lt. 5*nlat*imid+labc) return
      ierror = 5
      if (ldwork .lt. 2*(nlat+1)) return
      ierror = 0
      iw1 = 3*nlat*imid+1
      idz = (mmax*(nlat+nlat-mmax+1))/2
      CALL VEA1(NLAT,NLON,IMID,WVHAES,WVHAES(LZIMN+1),IDZ,
     +          WORK,WORK(IW1),DWORK)
      call hrffti(nlon,wvhaes(2*lzimn+1))
      return
      end
      subroutine vea1(nlat,nlon,imid,zv,zw,idz,zin,wzvin,dwork)
      dimension zv(idz,1),zw(idz,1),zin(imid,nlat,3),wzvin(1)
      double precision dwork(*)
      mmax = min0(nlat,(nlon+1)/2)
      call zvinit (nlat,nlon,wzvin,dwork)
      do 33 mp1=1,mmax
      m = mp1-1
      call zvin (0,nlat,nlon,m,zin,i3,wzvin)
      do 33 np1=mp1,nlat
      mn = m*(nlat-1)-(m*(m-1))/2+np1
      do 33 i=1,imid
      zv(mn,i) = zin(i,np1,i3)
   33 continue
      call zwinit (nlat,nlon,wzvin,dwork)
      do 34 mp1=1,mmax
      m = mp1-1
      call zwin (0,nlat,nlon,m,zin,i3,wzvin)
      do 34 np1=mp1,nlat
      mn = m*(nlat-1)-(m*(m-1))/2+np1
      do 34 i=1,imid
      zw(mn,i) = zin(i,np1,i3)
   34 continue
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file vhses.f
c
c     this file contains code and documentation for subroutines
c     vhses and vhsesi
c
c ... files which must be loaded with vhses.f
c
c     sphcom.f, hrfft.f
c
c   
c     subroutine vhses(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
c    +                 mdab,ndab,wvhses,lvhses,work,lwork,ierror)
c
c     subroutine vhses performs the vector spherical harmonic synthesis
c     of the arrays br, bi, cr, and ci and stores the result in the
c     arrays v and w. v(i,j) and w(i,j) are the colatitudinal 
c     (measured from the north pole) and east longitudinal components
c     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
c     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
c     representation of (v,w) is given below at output parameters v,w.
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     ityp   = 0  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon.   
c
c            = 1  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 2  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 3  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 4  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 5  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 6  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 7  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 8  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c
c     nt     the number of syntheses.  in the program that calls vhses,
c            the arrays v,w,br,bi,cr, and ci can be three dimensional
c            in which case multiple syntheses will be performed.
c            the third index is the synthesis index which assumes the 
c            values k=1,...,nt.  for a single synthesis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that all the arrays are two
c            dimensional.
c
c     idvw   the first dimension of the arrays v,w as it appears in
c            the program that calls vhaes. if ityp .le. 2 then idvw
c            must be at least nlat.  if ityp .gt. 2 and nlat is
c            even then idvw must be at least nlat/2. if ityp .gt. 2
c            and nlat is odd then idvw must be at least (nlat+1)/2.
c
c     jdvw   the second dimension of the arrays v,w as it appears in
c            the program that calls vhses. jdvw must be at least nlon.
c
c     br,bi  two or three dimensional arrays (see input parameter nt)
c     cr,ci  that contain the vector spherical harmonic coefficients
c            in the spectral representation of v(i,j) and w(i,j) given
c            below at the discription of output parameters v and w.
c
c     mdab   the first dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhses. mdab must be at
c            least min0(nlat,nlon/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     ndab   the second dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhses. ndab must be at
c            least nlat.
c
c     wvhses an array which must be initialized by subroutine vhsesi.
c            once initialized, wvhses can be used repeatedly by vhses
c            as long as nlon and nlat remain unchanged.  wvhses must
c            not be altered between calls of vhses.
c
c     lvhses the dimension of the array wvhses as it appears in the
c            program that calls vhses. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhses must be at least
c
c                 l1*l2*(nlat+nlat-l1+1)+nlon+15
c
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhses. define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            if ityp .le. 2 then lwork must be at least
c
c                       (2*nt+1)*nlat*nlon
c
c            if ityp .gt. 2 then lwork must be at least
c
c                        (2*nt+1)*l2*nlon 
c
c     **************************************************************
c
c     output parameters
c
c     v,w    two or three dimensional arrays (see input parameter nt)
c            in which the synthesis is stored. v is the colatitudinal
c            component and w is the east longitudinal component. 
c            v(i,j),w(i,j) contain the components at colatitude
c            theta(i) = (i-1)*pi/(nlat-1) and longitude phi(j) =
c            (j-1)*2*pi/nlon. the index ranges are defined above at
c            the input parameter ityp. v and w are computed from the 
c            formulas given below
c
c
c     define
c
c     1.  theta is colatitude and phi is east longitude
c
c     2.  the normalized associated legendre funnctions
c
c         pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)
c                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
c                        factorial(n)) times the (n+m)th derivative
c                        of (x**2-1)**n with respect to x=cos(theta)
c
c     3.  vbar(m,n,theta) = the derivative of pbar(m,n,theta) with
c                           respect to theta divided by the square
c                           root of n(n+1).
c
c         vbar(m,n,theta) is more easily computed in the form
c
c         vbar(m,n,theta) = (sqrt((n+m)*(n-m+1))*pbar(m-1,n,theta)
c         -sqrt((n-m)*(n+m+1))*pbar(m+1,n,theta))/(2*sqrt(n*(n+1)))
c
c     4.  wbar(m,n,theta) = m/(sin(theta))*pbar(m,n,theta) divided
c                           by the square root of n(n+1).
c
c         wbar(m,n,theta) is more easily computed in the form
c
c         wbar(m,n,theta) = sqrt((2n+1)/(2n-1))*(sqrt((n+m)*(n+m-1))
c         *pbar(m-1,n-1,theta)+sqrt((n-m)*(n-m-1))*pbar(m+1,n-1,theta))
c         /(2*sqrt(n*(n+1)))
c
c
c    the colatitudnal dependence of the normalized surface vector
c                spherical harmonics are defined by
c
c     5.    bbar(m,n,theta) = (vbar(m,n,theta),i*wbar(m,n,theta))
c
c     6.    cbar(m,n,theta) = (i*wbar(m,n,theta),-vbar(m,n,theta))
c
c
c    the coordinate to index mappings 
c
c     7.   theta(i) = (i-1)*pi/(nlat-1) and phi(j) = (j-1)*2*pi/nlon
c
c    
c     the maximum (plus one) longitudinal wave number
c
c     8.     mmax = min0(nlat,nlon/2) if nlon is even or
c            mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c    if we further define the output vector as
c
c     9.    h(i,j) = (v(i,j),w(i,j))
c
c    and the complex coefficients
c
c     10.   b(m,n) = cmplx(br(m+1,n+1),bi(m+1,n+1))
c
c     11.   c(m,n) = cmplx(cr(m+1,n+1),ci(m+1,n+1))
c
c
c    then for i=1,...,nlat and  j=1,...,nlon
c
c        the expansion for real h(i,j) takes the form
c
c     h(i,j) = the sum from n=1 to n=nlat-1 of the real part of
c
c         .5*(b(0,n)*bbar(0,n,theta(i))+c(0,n)*cbar(0,n,theta(i)))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
c     n=nlat-1 of the real part of
c
c              b(m,n)*bbar(m,n,theta(i))*exp(i*m*phi(j))
c             +c(m,n)*cbar(m,n,theta(i))*exp(i*m*phi(j))
c
c   *************************************************************
c
c   in terms of real variables this expansion takes the form
c
c             for i=1,...,nlat and  j=1,...,nlon
c
c     v(i,j) = the sum from n=1 to n=nlat-1 of
c
c               .5*br(1,n+1)*vbar(0,n,theta(i))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
c     n=nlat-1 of the real part of
c
c       (br(m+1,n+1)*vbar(m,n,theta(i))-ci(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *cos(m*phi(j))
c      -(bi(m+1,n+1)*vbar(m,n,theta(i))+cr(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *sin(m*phi(j))
c
c    and for i=1,...,nlat and  j=1,...,nlon
c
c     w(i,j) = the sum from n=1 to n=nlat-1 of
c
c              -.5*cr(1,n+1)*vbar(0,n,theta(i))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
c     n=nlat-1 of the real part of
c
c      -(cr(m+1,n+1)*vbar(m,n,theta(i))+bi(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *cos(m*phi(j))
c      +(ci(m+1,n+1)*vbar(m,n,theta(i))-br(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *sin(m*phi(j))
c
c
c      br(m+1,nlat),bi(m+1,nlat),cr(m+1,nlat), and ci(m+1,nlat) are
c      assumed zero for m even.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of ityp
c            = 4  error in the specification of nt
c            = 5  error in the specification of idvw
c            = 6  error in the specification of jdvw
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lvhses
c            = 10 error in the specification of lwork
c
c ************************************************************
c
c     subroutine vhsesi(nlat,nlon,wvhses,lvhses,work,lwork,dwork,
c    +                  ldwork,ierror)
c
c     subroutine vhsesi initializes the array wvhses which can then be
c     used repeatedly by subroutine vhses until nlat or nlon is changed.
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     lvhses the dimension of the array wvhses as it appears in the
c            program that calls vhses. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhses must be at least
c
c                  l1*l2*(nlat+nlat-l1+1)+nlon+15
c
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhses. lwork must be at least
c
c              3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+5*l2*nlat
c
c     dwork  an unsaved double precision work space
c
c     ldwork the length of the array dwork as it appears in the
c            program that calls vhsesi.  ldwork must be at least
c            2*(nlat+1)
c
c
c     **************************************************************
c
c     output parameters
c
c     wvhses an array which is initialized for use by subroutine vhses.
c            once initialized, wvhses can be used repeatedly by vhses
c            as long as nlat or nlon remain unchanged.  wvhses must not
c            be altered between calls of vhses.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lvhses
c            = 4  error in the specification of lwork
c            = 5  error in the specification of ldwork
c
c *****************************************
      subroutine vhses(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
     +                 mdab,ndab,wvhses,lvhses,work,lwork,ierror)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          work(1),wvhses(1)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      if(ityp.lt.0 .or. ityp.gt.8) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      imid = (nlat+1)/2
      if((ityp.le.2 .and. idvw.lt.nlat) .or.
     1   (ityp.gt.2 .and. idvw.lt.imid)) return
      ierror = 6
      if(jdvw .lt. nlon) return
      ierror = 7
      mmax = min0(nlat,(nlon+1)/2)
      if(mdab .lt. mmax) return
      ierror = 8
      if(ndab .lt. nlat) return
      ierror = 9
      idz = (mmax*(nlat+nlat-mmax+1))/2
      lzimn = idz*imid
      if(lvhses .lt. lzimn+lzimn+nlon+15) return
      ierror = 10
      idv = nlat
      if(ityp .gt. 2) idv = imid
      lnl = nt*idv*nlon
      if(lwork .lt. lnl+lnl+idv*nlon) return
      ierror = 0
      ist = 0
      if(ityp .le. 2) ist = imid
      iw1 = ist+1
      iw2 = lnl+1
      iw3 = iw2+ist
      iw4 = iw2+lnl
      jw1 = lzimn+1
      jw2 = jw1+lzimn
      call vhses1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,ndab,
     1     br,bi,cr,ci,idv,work,work(iw1),work(iw2),work(iw3),
     2     work(iw4),idz,wvhses,wvhses(jw1),wvhses(jw2))
      return
      end

      subroutine vhses1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,
     1   ndab,br,bi,cr,ci,idv,ve,vo,we,wo,work,idz,vb,wb,wrfft)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          ve(idv,nlon,1),vo(idv,nlon,1),we(idv,nlon,1),
     3          wo(idv,nlon,1),work(1),wrfft(1),
     4          vb(imid,1),wb(imid,1)
      nlp1 = nlat+1
      mlat = mod(nlat,2)
      mlon = mod(nlon,2)
      mmax = min0(nlat,(nlon+1)/2)
      imm1 = imid
      if(mlat .ne. 0) imm1 = imid-1
      do 10 k=1,nt
      do 10 j=1,nlon
      do 10 i=1,idv
      ve(i,j,k) = 0.
      we(i,j,k) = 0.
   10 continue
      ndo1 = nlat
      ndo2 = nlat
      if(mlat .ne. 0) ndo1 = nlat-1
      if(mlat .eq. 0) ndo2 = nlat-1
   18 itypp = ityp+1
      go to (1,100,200,300,400,500,600,700,800),itypp
c
c     case ityp=0   no symmetries
c
c     case m = 0
c
    1 do 15 k=1,nt
      do 15 np1=2,ndo2,2
      do 15 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1)
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1)
   15 continue
      do 16 k=1,nt
      do 16 np1=3,ndo1,2
      do 16 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1)
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1)
   16 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 30 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 26
      do 25 k=1,nt
      do 24 np1=mp1,ndo1,2
      mn = mb+np1
      do 23 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
   23 continue
      if(mlat .eq. 0) go to 24
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,mn)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,mn)
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,mn) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,mn)
   24 continue
   25 continue
   26 if(mp2 .gt. ndo2) go to 30
      do 29 k=1,nt
      do 28 np1=mp2,ndo2,2
      mn = mb+np1
      do 27 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
   27 continue
      if(mlat .eq. 0) go to 28
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,mn) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,mn)
   28 continue
   29 continue
   30 continue
      go to 950
c
c     case ityp=1   no symmetries,  cr and ci equal zero
c
c     case m = 0
c
  100 continue
      do 115 k=1,nt
      do 115 np1=2,ndo2,2
      do 115 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1)
  115 continue
      do 116 k=1,nt
      do 116 np1=3,ndo1,2
      do 116 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1)
  116 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 130 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 126
      do 125 k=1,nt
      do 124 np1=mp1,ndo1,2
      mn = mb+np1
      do 123 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  123 continue
      if(mlat .eq. 0) go to 124
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,mn) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,mn)
  124 continue
  125 continue
  126 if(mp2 .gt. ndo2) go to 130
      do 129 k=1,nt
      do 128 np1=mp2,ndo2,2
      mn = mb+np1
      do 127 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  127 continue
      if(mlat .eq. 0) go to 128
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,mn) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,mn)
  128 continue
  129 continue
  130 continue
      go to 950
c
c     case ityp=2   no symmetries,  br and bi are equal to zero
c
c     case m = 0
c
  200 do 215 k=1,nt
      do 215 np1=2,ndo2,2
      do 215 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1)
  215 continue
      do 216 k=1,nt
      do 216 np1=3,ndo1,2
      do 216 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1)
  216 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 230 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 226
      do 225 k=1,nt
      do 224 np1=mp1,ndo1,2
      mn = mb+np1
      do 223 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  223 continue
      if(mlat .eq. 0) go to 224
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,mn)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,mn)
  224 continue
  225 continue
  226 if(mp2 .gt. ndo2) go to 230
      do 229 k=1,nt
      do 228 np1=mp2,ndo2,2
      mn = mb+np1
      do 227 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  227 continue
      if(mlat .eq. 0) go to 228
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,mn)
  228 continue
  229 continue
  230 continue
      go to 950
c
c     case ityp=3   v even,  w odd 
c
c     case m = 0
c
  300 do 315 k=1,nt
      do 315 np1=2,ndo2,2
      do 315 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1)
  315 continue
      do 316 k=1,nt
      do 316 np1=3,ndo1,2
      do 316 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1)
  316 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 330 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 326
      do 325 k=1,nt
      do 324 np1=mp1,ndo1,2
      mn = mb+np1
      do 323 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  323 continue
      if(mlat .eq. 0) go to 324
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,mn)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,mn)
  324 continue
  325 continue
  326 if(mp2 .gt. ndo2) go to 330
      do 329 k=1,nt
      do 328 np1=mp2,ndo2,2
      mn = mb+np1
      do 327 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  327 continue
      if(mlat .eq. 0) go to 328
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,mn) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,mn)
  328 continue
  329 continue
  330 continue
      go to 950
c
c     case ityp=4   v even,  w odd, and both cr and ci equal zero 
c
c     case m = 0
c
  400 do 415 k=1,nt
      do 415 np1=2,ndo2,2
      do 415 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1)
  415 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 430 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 430
      do 429 k=1,nt
      do 428 np1=mp2,ndo2,2
      mn = mb+np1
      do 427 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  427 continue
      if(mlat .eq. 0) go to 428
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,mn) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,mn)
  428 continue
  429 continue
  430 continue
      go to 950
c
c     case ityp=5   v even,  w odd,     br and bi equal zero 
c
c     case m = 0
c
  500 do 516 k=1,nt
      do 516 np1=3,ndo1,2
      do 516 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1)
  516 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 530 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 530
      do 525 k=1,nt
      do 524 np1=mp1,ndo1,2
      mn = mb+np1
      do 523 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  523 continue
      if(mlat .eq. 0) go to 524
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,mn)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,mn)
  524 continue
  525 continue
  530 continue
      go to 950
c
c     case ityp=6   v odd  ,  w even
c
c     case m = 0
c
  600 do 615 k=1,nt
      do 615 np1=2,ndo2,2
      do 615 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1)
  615 continue
      do 616 k=1,nt
      do 616 np1=3,ndo1,2
      do 616 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1)
  616 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 630 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 626
      do 625 k=1,nt
      do 624 np1=mp1,ndo1,2
      mn = mb+np1
      do 623 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  623 continue
      if(mlat .eq. 0) go to 624
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,mn) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,mn)
  624 continue
  625 continue
  626 if(mp2 .gt. ndo2) go to 630
      do 629 k=1,nt
      do 628 np1=mp2,ndo2,2
      mn = mb+np1
      do 627 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  627 continue
      if(mlat .eq. 0) go to 628
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,mn)
  628 continue
  629 continue
  630 continue
      go to 950
c
c     case ityp=7   v odd, w even   cr and ci equal zero
c
c     case m = 0
c
  700 do 716 k=1,nt
      do 716 np1=3,ndo1,2
      do 716 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1)
  716 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 730 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 730
      do 725 k=1,nt
      do 724 np1=mp1,ndo1,2
      mn = mb+np1
      do 723 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  723 continue
      if(mlat .eq. 0) go to 724
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,mn) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,mn)
  724 continue
  725 continue
  730 continue
      go to 950
c
c     case ityp=8   v odd,  w even   br and bi equal zero
c
c     case m = 0
c
  800 do 815 k=1,nt
      do 815 np1=2,ndo2,2
      do 815 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1)
  815 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 830 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 830
      do 829 k=1,nt
      do 828 np1=mp2,ndo2,2
      mn = mb+np1
      do 827 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  827 continue
      if(mlat .eq. 0) go to 828
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,mn)
  828 continue
  829 continue
  830 continue
  950 do 14 k=1,nt
      call hrfftb(idv,nlon,ve(1,1,k),idv,wrfft,work)
      call hrfftb(idv,nlon,we(1,1,k),idv,wrfft,work)
   14 continue
      if(ityp .gt. 2) go to 12
      do 60 k=1,nt
      do 60 j=1,nlon
      do 60 i=1,imm1
      v(i,j,k) = .5*(ve(i,j,k)+vo(i,j,k))
      w(i,j,k) = .5*(we(i,j,k)+wo(i,j,k))
      v(nlp1-i,j,k) = .5*(ve(i,j,k)-vo(i,j,k))
      w(nlp1-i,j,k) = .5*(we(i,j,k)-wo(i,j,k))
   60 continue
      go to 13
   12 do 11 k=1,nt
      do 11 j=1,nlon
      do 11 i=1,imm1
      v(i,j,k) = .5*ve(i,j,k)
      w(i,j,k) = .5*we(i,j,k)
   11 continue
   13 if(mlat .eq. 0) return
      do 65 k=1,nt
      do 65 j=1,nlon
      v(imid,j,k) = .5*ve(imid,j,k)
      w(imid,j,k) = .5*we(imid,j,k)
   65 continue
      return
      end

      subroutine vhsesi(nlat,nlon,wvhses,lvhses,work,lwork,dwork,
     +                  ldwork,ierror)
      dimension wvhses(lvhses),work(lwork)
      double precision dwork(ldwork)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      mmax = min0(nlat,(nlon+1)/2)
      imid = (nlat+1)/2
      lzimn = (imid*mmax*(nlat+nlat-mmax+1))/2
      if(lvhses .lt. lzimn+lzimn+nlon+15) return
      ierror = 4
      labc = 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      if(lwork .lt. 5*nlat*imid+labc) return
      ierror = 5
      if (ldwork .lt. 2*(nlat+1)) return
      ierror = 0
      iw1 = 3*nlat*imid+1
      idz = (mmax*(nlat+nlat-mmax+1))/2
      call ves1(nlat,nlon,imid,wvhses,wvhses(lzimn+1),idz,work,
     1                                         work(iw1),dwork)
      call hrffti(nlon,wvhses(2*lzimn+1))
      return
      end
      subroutine ves1(nlat,nlon,imid,vb,wb,idz,vin,wzvin,dwork)
      dimension vb(imid,*),wb(imid,*),vin(imid,nlat,3),wzvin(*)
      double precision dwork(*)
      mmax = min0(nlat,(nlon+1)/2)
      call vbinit (nlat,nlon,wzvin,dwork)
      do 33 mp1=1,mmax
      m = mp1-1
      call vbin (0,nlat,nlon,m,vin,i3,wzvin)
      do 33 np1=mp1,nlat
      mn = m*(nlat-1)-(m*(m-1))/2+np1
      do 33 i=1,imid
      vb(i,mn) = vin(i,np1,i3)
   33 continue
      call wbinit (nlat,nlon,wzvin,dwork)
      do 34 mp1=1,mmax
      m = mp1-1
      call wbin (0,nlat,nlon,m,vin,i3,wzvin)
      do 34 np1=mp1,nlat
      mn = m*(nlat-1)-(m*(m-1))/2+np1
      do 34 i=1,imid
      wb(i,mn) = vin(i,np1,i3)
   34 continue
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file shags.f
c
c     this file contains code and documentation for subroutines
c     shags and shagsi
c
c ... files which must be loaded with shags.f
c
c     sphcom.f, hrfft.f, gaqd.f
c
c     subroutine shags(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
c    1                    wshags,lshags,work,lwork,ierror)
c
c     subroutine shags performs the spherical harmonic analysis
c     on the array g and stores the result in the arrays a and b.
c     the analysis is performed on a gaussian grid in colatitude
c     and an equally spaced grid in longitude.  the associated
c     legendre functions are stored rather than recomputed as they
c     are in subroutine shagc.  the analysis is described below
c     at output parameters a,b.
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are compu
c            in radians in theta(1),...,theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid poi
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c     isym   = 0  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 array g(i,j) for i=1,...,nlat and j=1,...,nlon.
c                 (see description of g below)
c
c            = 1  g is antisymmetric about the equator. the analysis
c                 is performed on the northern hemisphere only.  i.e.
c                 if nlat is odd the analysis is performed on the
c                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
c                 if nlat is even the analysis is performed on the
c                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c
c            = 2  g is symmetric about the equator. the analysis is
c                 performed on the northern hemisphere only.  i.e.
c                 if nlat is odd the analysis is performed on the
c                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
c                 if nlat is even the analysis is performed on the
c                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c     nt     the number of analyses.  in the program that calls shags,
c            the arrays g,a and b can be three dimensional in which
c            case multiple analyses will be performed.  the third
c            index is the analysis index which assumes the values
c            k=1,...,nt.  for a single analysis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that the arrays g,a and b
c            have only two dimensions.
c
c     g      a two or three dimensional array (see input parameter
c            nt) that contains the discrete function to be analyzed.
c            g(i,j) contains the value of the function at the gaussian
c            point theta(i) and longitude point phi(j) = (j-1)*2*pi/nlon
c            the index ranges are defined above at the input parameter
c            isym.
c
c     idg    the first dimension of the array g as it appears in the
c            program that calls shags. if isym equals zero then idg
c            must be at least nlat.  if isym is nonzero then idg must
c            be at least nlat/2 if nlat is even or at least (nlat+1)/2
c            if nlat is odd.
c
c     jdg    the second dimension of the array g as it appears in the
c            program that calls shags. jdg must be at least nlon.
c
c     mdab   the first dimension of the arrays a and b as it appears
c            in the program that calls shags. mdab must be at least
c            min0((nlon+2)/2,nlat) if nlon is even or at least
c            min0((nlon+1)/2,nlat) if nlon is odd.
c
c     ndab   the second dimension of the arrays a and b as it appears
c            in the program that calls shags. ndab must be at least nlat
c
c     wshags an array which must be initialized by subroutine shagsi.
c            once initialized, wshags can be used repeatedly by shags
c            as long as nlat and nlon remain unchanged.  wshags must
c            not be altered between calls of shags.
c
c     lshags the dimension of the array wshags as it appears in the
c            program that calls shags. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshags must be at least
c
c            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
c
c     work   a real work space which need not be saved
c
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls shags. define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c
c            if isym is zero then lwork must be at least
c
c                  nlat*nlon*(nt+1)
c
c            if isym is nonzero then lwork must be at least
c
c                  l2*nlon*(nt+1)
c
c     **************************************************************
c
c     output parameters
c
c     a,b    both a,b are two or three dimensional arrays (see input
c            parameter nt) that contain the spherical harmonic
c            coefficients in the representation of g(i,j) given in the
c            discription of subroutine shags. for isym=0, a(m,n) and
c            b(m,n) are given by the equations listed below. symmetric
c            versions are used when isym is greater than zero.
c
c     definitions
c
c     1. the normalized associated legendre functions
c
c     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
c                       *sin(theta)**m/(2**n*factorial(n)) times the
c                       (n+m)th derivative of (x**2-1)**n with respect
c                       to x=cos(theta).
c
c     2. the fourier transform of g(i,j).
c
c     c(m,i)          = 2/nlon times the sum from j=1 to j=nlon of
c                       g(i,j)*cos((m-1)*(j-1)*2*pi/nlon)
c                       (the first and last terms in this sum
c                       are divided by 2)
c
c     s(m,i)          = 2/nlon times the sum from j=2 to j=nlon of
c                       g(i,j)*sin((m-1)*(j-1)*2*pi/nlon)
c
c
c     3. the gaussian points and weights on the sphere
c        (computed by subroutine gaqd).
c
c        theta(1),...,theta(nlat) (gaussian pts in radians)
c        wts(1),...,wts(nlat) (corresponding gaussian weights)
c
c
c     4. the maximum (plus one) longitudinal wave number
c
c            mmax = min0(nlat,(nlon+2)/2) if nlon is even or
c            mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c
c     then for m=0,...,mmax-1 and n=m,...,nlat-1 the arrays a,b
c     are given by
c
c     a(m+1,n+1)     =  the sum from i=1 to i=nlat of
c                       c(m+1,i)*wts(i)*pbar(m,n,theta(i))
c
c     b(m+1,n+1)      = the sum from i=1 to nlat of
c                       s(m+1,i)*wts(i)*pbar(m,n,theta(i))
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of isym
c            = 4  error in the specification of nt
c            = 5  error in the specification of idg
c            = 6  error in the specification of jdg
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lshags
c            = 10 error in the specification of lwork
c
c
c ****************************************************************
c
c     subroutine shagsi(nlat,nlon,wshags,lshags,work,lwork,dwork,ldwork,
c    +                  ierror)
c
c     subroutine shagsi initializes the array wshags which can then
c     be used repeatedly by subroutines shags. it precomputes
c     and stores in wshags quantities such as gaussian weights,
c     legendre polynomial coefficients, and fft trigonometric tables.
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are compu
c            in radians in theta(1),...,theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid poi
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than or equal to 4. the efficiency of the computation is
c            improved when nlon is a product of small prime numbers.
c
c     wshags an array which must be initialized by subroutine shagsi.
c            once initialized, wshags can be used repeatedly by shags
c            as long as nlat and nlon remain unchanged.  wshags must
c            not be altered between calls of shags.
c
c     lshags the dimension of the array wshags as it appears in the
c            program that calls shags. define
c
c               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lshags must be at least
c
c            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
c
c     work   a real work space which need not be saved
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls shagsi. lwork must be at least
c            4*nlat*(nlat+2)+2 in the routine calling shagsi
c
c     dwork   a double precision work array that does not have to be saved.
c
c     ldwork  the length of dwork in the calling routine.  ldwork must
c             be at least nlat*(nlat+4)
c
c     output parameter
c
c     wshags an array which must be initialized before calling shags or
c            once initialized, wshags can be used repeatedly by shags or
c            as long as nlat and nlon remain unchanged.  wshags must not
c            altered between calls of shasc.
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lshags
c            = 4  error in the specification of lwork
c            = 5  error in the specification of ldwork
c            = 6  failure in gaqd to compute gaussian points
c                 (due to failure in eigenvalue routine)
c
c
c ****************************************************************
      subroutine shags(nlat,nlon,mode,nt,g,idg,jdg,a,b,mdab,ndab,
     1                    wshags,lshags,work,lwork,ierror)
c     subroutine shags performs the spherical harmonic analysis on
c     a gaussian grid on the array(s) in g and returns the coefficients
c     in array(s) a,b. the necessary legendre polynomials are fully
c     stored in this version.
c
      dimension g(idg,jdg,1),a(mdab,ndab,1),b(mdab,ndab,1),
     1          wshags(lshags),work(lwork)
c     check input parameters
      ierror = 1
      if (nlat.lt.3) return
      ierror = 2
      if (nlon.lt.4) return
      ierror = 3
      if (mode.lt.0 .or.mode.gt.2) return
c     set m limit for pmn
      l = min0((nlon+2)/2,nlat)
c     set gaussian point nearest equator pointer
      late = (nlat+mod(nlat,2))/2
c     set number of grid points for analysis/synthesis
      lat = nlat
      if (mode.ne.0) lat = late
      ierror = 4
      if (nt.lt.1) return
      ierror = 5
      if (idg.lt.lat) return
      ierror = 6
      if (jdg.lt.nlon) return
      ierror = 7
      if(mdab .lt. l) return
      ierror = 8
      if(ndab .lt. nlat) return
      l1 = l
      l2 = late
      ierror = 9
c     check permanent work space length
c
      lp= nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
      if(lshags.lt.lp) return
      ierror = 10
c     check temporary work space length
      if (mode.eq.0 .and. lwork.lt.nlat*nlon*(nt+1)) return
      if (mode.ne.0 .and. lwork.lt.l2*nlon*(nt+1)) return
      ierror = 0
c     set starting address for gaussian wts ,fft values,
c     and fully stored legendre polys in wshags
      iwts = 1
      ifft = nlat+2*nlat*late+3*(l*(l-1)/2+(nlat-l)*(l-1))+1
      ipmn = ifft+nlon+15
c     set pointer for internal storage of g
      iw = lat*nlon*nt+1
      call shags1(nlat,nlon,l,lat,mode,g,idg,jdg,nt,a,b,mdab,ndab,
     1wshags(iwts),wshags(ifft),wshags(ipmn),late,work,work(iw))
      return
      end

      subroutine shags1(nlat,nlon,l,lat,mode,gs,idg,jdg,nt,a,b,mdab,
     1                  ndab,wts,wfft,pmn,late,g,work)
      dimension gs(idg,jdg,nt),a(mdab,ndab,nt),
     1          b(mdab,ndab,nt),g(lat,nlon,nt)
      dimension wfft(1),pmn(late,1),wts(nlat),work(1)
c     set gs array internally in shags1
      do 100 k=1,nt
      do 100 j=1,nlon
      do 100 i=1,lat
      g(i,j,k) = gs(i,j,k)
  100 continue

c     do fourier transform
      do 101 k=1,nt
      call hrfftf(lat,nlon,g(1,1,k),lat,wfft,work)
  101 continue

c     scale result
      sfn = 2.0/float(nlon)
      do 102 k=1,nt
      do 102 j=1,nlon
      do 102 i=1,lat
      g(i,j,k) = sfn*g(i,j,k)
  102 continue

c     compute using gaussian quadrature
c     a(n,m) = s (ga(theta,m)*pnm(theta)*sin(theta)*dtheta)
c     b(n,m) = s (gb(theta,m)*pnm(theta)*sin(theta)*dtheta)
c     here ga,gb are the cos(phi),sin(phi) coefficients of
c     the fourier expansion of g(theta,phi) in phi.  as a result
c     of the above fourier transform they are stored in array
c     g as follows:
c     for each theta(i) and k= l-1
c     ga(0),ga(1),gb(1),ga(2),gb(2),...,ga(k-1),gb(k-1),ga(k)
c     correspond to
c     g(i,1),g(i,2),g(i,3),g(i,4),g(i,5),...,g(i,2l-4),g(i,2l-3),g(i,2l-2)
c     whenever 2*l-2 = nlon exactly
c     initialize coefficients to zero
      do 103 k=1,nt
      do 103 np1=1,nlat
      do 103 mp1=1,l
      a(mp1,np1,k) = 0.0
      b(mp1,np1,k) = 0.0
  103 continue
c     set mp1 limit on b(mp1) calculation
      lm1 = l
      if (nlon .eq. l+l-2) lm1 = l-1

      if (mode.eq.0) then
c     for full sphere (mode=0) and even/odd reduction:
c     overwrite g(i) with (g(i)+g(nlat-i+1))*wts(i)
c     overwrite g(nlat-i+1) with (g(i)-g(nlat-i+1))*wts(i)
      nl2 = nlat/2
      do 104 k=1,nt
      do 104 j=1,nlon
      do 105 i=1,nl2
      is = nlat-i+1
      t1 = g(i,j,k)
      t2 = g(is,j,k)
      g(i,j,k) = wts(i)*(t1+t2)
      g(is,j,k) = wts(i)*(t1-t2)
  105 continue
c     adjust equator if necessary(nlat odd)
      if (mod(nlat,2).ne.0) g(late,j,k) = wts(late)*g(late,j,k)
  104 continue
c     set m = 0 coefficients first
      mp1 = 1
      m = 0
      mml1 = m*(2*nlat-m-1)/2
      do 106 k=1,nt
      do 106 i=1,late
      is = nlat-i+1
      do 107 np1=1,nlat,2
c     n even
      a(1,np1,k) = a(1,np1,k)+g(i,1,k)*pmn(i,mml1+np1)
  107 continue
      do 108 np1=2,nlat,2
c     n odd
      a(1,np1,k) = a(1,np1,k)+g(is,1,k)*pmn(i,mml1+np1)
  108 continue
  106 continue
c     compute m.ge.1  coefficients next
      do 109 mp1=2,lm1
      m = mp1-1
      mml1 = m*(2*nlat-m-1)/2
      mp2 = mp1+1
      do 110 k=1,nt
      do 111 i=1,late
      is = nlat-i+1
c     n-m even
      do 112 np1=mp1,nlat,2
      a(mp1,np1,k) = a(mp1,np1,k)+g(i,2*m,k)*pmn(i,mml1+np1)
      b(mp1,np1,k) = b(mp1,np1,k)+g(i,2*m+1,k)*pmn(i,mml1+np1)
  112 continue
c     n-m odd
      do 113 np1=mp2,nlat,2
      a(mp1,np1,k) = a(mp1,np1,k)+g(is,2*m,k)*pmn(i,mml1+np1)
      b(mp1,np1,k) = b(mp1,np1,k)+g(is,2*m+1,k)*pmn(i,mml1+np1)
  113 continue
  111 continue
  110 continue
  109 continue
      if (nlon .eq. l+l-2) then
c     compute m=l-1, n=l-1,l,...,nlat-1 coefficients
      m = l-1
      mml1 = m*(2*nlat-m-1)/2
      do 114 k=1,nt
      do 114 i=1,late
      is = nlat-i+1
      do 124 np1=l,nlat,2
      mn = mml1+np1
      a(l,np1,k) = a(l,np1,k)+0.5*g(i,nlon,k)*pmn(i,mn)
  124 continue
c     n-m  odd
      lp1 = l+1
      do 125 np1=lp1,nlat,2
      mn = mml1+np1
      a(l,np1,k) = a(l,np1,k)+0.5*g(is,nlon,k)*pmn(i,mn)
  125 continue
  114 continue
      end if

      else

c     half sphere
c     overwrite g(i) with wts(i)*(g(i)+g(i)) for i=1,...,nlate/2
      nl2 = nlat/2
      do 116  k=1,nt
      do 116 j=1,nlon
      do 115 i=1,nl2
      g(i,j,k) = wts(i)*(g(i,j,k)+g(i,j,k))
  115 continue
c     adjust equator separately if a grid point
      if (nl2.lt.late) g(late,j,k) = wts(late)*g(late,j,k)
  116 continue

c     set m = 0 coefficients first
      mp1 = 1
      m = 0
      mml1 = m*(2*nlat-m-1)/2
      ms = 1
      if (mode.eq.1) ms = 2
      do 117 k=1,nt
      do 117 i=1,late
      do 117 np1=ms,nlat,2
      a(1,np1,k) = a(1,np1,k)+g(i,1,k)*pmn(i,mml1+np1)
  117 continue

c     compute m.ge.1  coefficients next
      do 118 mp1=2,lm1
      m = mp1-1
      mml1 = m*(2*nlat-m-1)/2
      ms = mp1
      if (mode.eq.1) ms = mp1+1
      do 119 k=1,nt
      do 119 i=1,late
      do 119 np1=ms,nlat,2
      a(mp1,np1,k) = a(mp1,np1,k)+g(i,2*m,k)*pmn(i,mml1+np1)
      b(mp1,np1,k) = b(mp1,np1,k)+g(i,2*m+1,k)*pmn(i,mml1+np1)
  119 continue
  118 continue

      if (nlon.eq.l+l-2) then
c     compute n=m=l-1 coefficients last
      m = l-1
      mml1 = m*(2*nlat-m-1)/2
c     set starting n for mode even
      ns = l
c     set starting n for mode odd
      if (mode.eq.1) ns = l+1
      do 120 k=1,nt
      do 120 i=1,late
      do 120 np1=ns,nlat,2
      mn = mml1+np1
      a(l,np1,k) = a(l,np1,k)+0.5*g(i,nlon,k)*pmn(i,mn)
  120 continue
      end if

      end if

      return
      end
      subroutine shagsi(nlat,nlon,wshags,lshags,work,lwork,dwork,ldwork,
     +                  ierror)
c
c     this subroutine must be called before calling shags or shsgs with
c     fixed nlat,nlon. it precomputes the gaussian weights, points
c     and all necessary legendre polys and stores them in wshags.
c     these quantities must be preserved when calling shags or shsgs
c     repeatedly with fixed nlat,nlon.  dwork must be of length at
c     least nlat*(nlat+4) in the routine calling shagsi.  This is
c     not checked.  undetectable errors will result if dwork is
c     smaller than nlat*(nlat+4).
c
      dimension wshags(lshags),work(lwork)
      double precision dwork(ldwork)
      ierror = 1
      if (nlat.lt.3) return
      ierror = 2
      if (nlon.lt.4) return
c     set triangular truncation limit for spherical harmonic basis
      l = min0((nlon+2)/2,nlat)
c     set equator or nearest point (if excluded) pointer
      late = (nlat+1)/2
      l1 = l
      l2 = late
c     check permanent work space length
      ierror = 3
      lp=nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
      if(lshags.lt.lp) return
      ierror = 4
c     check temporary work space
      if (lwork.lt.4*nlat*(nlat+2)+2) return
      ierror = 5
c     check double precision temporary space
      if (ldwork .lt. nlat*(nlat+4)) return
      ierror = 0
c     set preliminary quantites needed to compute and store legendre polys
      ldw = nlat*(nlat+4)
      call shagsp(nlat,nlon,wshags,lshags,dwork,ldwork,ierror)
      if (ierror.ne.0) return
c     set legendre poly pointer in wshags
      ipmnf = nlat+2*nlat*late+3*(l*(l-1)/2+(nlat-l)*(l-1))+nlon+16
      call shagss1(nlat,l,late,wshags,work,wshags(ipmnf))
      return
      end
      subroutine shagss1(nlat,l,late,w,pmn,pmnf)
      dimension w(1),pmn(nlat,late,3),pmnf(late,1)
c     compute and store legendre polys for i=1,...,late,m=0,...,l-1
c     and n=m,...,l-1
      do i=1,nlat
	do j=1,late
	  do k=1,3
	   pmn(i,j,k) = 0.0
	  end do
	end do
      end do
      do 100 mp1=1,l
      m = mp1-1
      mml1 = m*(2*nlat-m-1)/2
c     compute pmn for n=m,...,nlat-1 and i=1,...,(l+1)/2
      mode = 0
      call legin(mode,l,nlat,m,w,pmn,km)
c     store above in pmnf
      do 101 np1=mp1,nlat
      mn = mml1+np1
      do 102 i=1,late
      pmnf(i,mn) = pmn(np1,i,km)
  102 continue
  101 continue
  100 continue
      return
      end
      subroutine shagsp(nlat,nlon,wshags,lshags,dwork,ldwork,ierror)
      dimension wshags(lshags)
      double precision dwork(ldwork)
      ierror = 1
      if (nlat.lt.3) return
      ierror = 2
      if (nlon.lt.4) return
c     set triangular truncation limit for spherical harmonic basis
      l = min0((nlon+2)/2,nlat)
c     set equator or nearest point (if excluded) pointer
      late = (nlat+mod(nlat,2))/2
      l1 = l
      l2 = late
      ierror = 3
c     check permanent work space length
      if (lshags .lt. nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
      ierror = 4
c     if (lwork.lt.4*nlat*(nlat+2)+2) return
      if (ldwork.lt.nlat*(nlat+4))return
      ierror = 0
c     set pointers
      i1 = 1
      i2 = i1+nlat
      i3 = i2+nlat*late
      i4 = i3+nlat*late
      i5 = i4+l*(l-1)/2 +(nlat-l)*(l-1)
      i6 = i5+l*(l-1)/2 +(nlat-l)*(l-1)
      i7 = i6+l*(l-1)/2 +(nlat-l)*(l-1)
c     set indices in temp work for double precision gaussian wts and pts
      idth = 1
c     idwts = idth+2*nlat
c     iw = idwts+2*nlat
      idwts = idth+nlat
      iw = idwts+nlat
      call shagsp1(nlat,nlon,l,late,wshags(i1),wshags(i2),wshags(i3),
     1wshags(i4),wshags(i5),wshags(i6),wshags(i7),dwork(idth),
     2dwork(idwts),dwork(iw),ierror)
      if (ierror.ne.0) ierror = 6
      return
      end

      subroutine shagsp1(nlat,nlon,l,late,wts,p0n,p1n,abel,bbel,cbel,
     +                   wfft,dtheta,dwts,work,ier)
      dimension wts(nlat),p0n(nlat,late),p1n(nlat,late),abel(1),bbel(1),
     1 cbel(1),wfft(1),dtheta(nlat),dwts(nlat)
      double precision pb,dtheta,dwts,work(*)
      indx(m,n) = (n-1)*(n-2)/2+m-1
      imndx(m,n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1
      call hrffti(nlon,wfft)

c     compute double precision gaussian points and weights
c     lw = 4*nlat*(nlat+2)
      lw = nlat*(nlat+2)
      call gaqd(nlat,dtheta,dwts,work,lw,ier)
      if (ier.ne.0) return

c     store gaussian weights single precision to save computation
c     in inner loops in analysis
      do 100 i=1,nlat
      wts(i) = dwts(i)
  100 continue
c     initialize p0n,p1n using double precision dnlfk,dnlft
      do 101 np1=1,nlat
      do 101 i=1,late
      p0n(np1,i) = 0.0
      p1n(np1,i) = 0.0
  101 continue
c     compute m=n=0 legendre polynomials for all theta(i)
      np1 = 1
      n = 0
      m = 0
      call dnlfk(m,n,work)
      do 103 i=1,late
      call dnlft(m,n,dtheta(i),work,pb)
      p0n(1,i) = pb
  103 continue
c     compute p0n,p1n for all theta(i) when n.gt.0
      do 104 np1=2,nlat
      n = np1-1
      m = 0
      call dnlfk(m,n,work)
      do 105 i=1,late
      call dnlft(m,n,dtheta(i),work,pb)
      p0n(np1,i) = pb
  105 continue
c     compute m=1 legendre polynomials for all n and theta(i)
      m = 1
      call dnlfk(m,n,work)
      do 106 i=1,late
      call dnlft(m,n,dtheta(i),work,pb)
      p1n(np1,i) = pb
  106 continue
  104 continue
c
c     compute and store swarztrauber recursion coefficients
c     for 2.le.m.le.n and 2.le.n.le.nlat in abel,bbel,cbel
      do 107 n=2,nlat
      mlim = min0(n,l)
      do 107 m=2,mlim
      imn = indx(m,n)
      if (n.ge.l) imn = imndx(m,n)
      abel(imn)=sqrt(float((2*n+1)*(m+n-2)*(m+n-3))/
     1               float(((2*n-3)*(m+n-1)*(m+n))))
      bbel(imn)=sqrt(float((2*n+1)*(n-m-1)*(n-m))/
     1               float(((2*n-3)*(m+n-1)*(m+n))))
      cbel(imn)=sqrt(float((n-m+1)*(n-m+2))/
     1               float(((n+m-1)*(n+m))))
  107 continue
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file vhags.f
c
c     this file contains code and documentation for subroutines
c     vhags and vhagsi
c
c ... files which must be loaded with vhags.f
c
c     sphcom.f, hrfft.f, gaqd.f
c
c     subroutine vhags(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
c    +                 mdab,ndab,wvhags,lvhags,work,lwork,ierror)
c
c     subroutine vhags performs the vector spherical harmonic analysis
c     on the vector field (v,w) and stores the result in the arrays
c     br, bi, cr, and ci. v(i,j) and w(i,j) are the colatitudinal 
c     (measured from the north pole) and east longitudinal components
c     respectively, located at the gaussian colatitude point theta(i)
c     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
c     representation of (v,w) is given at output parameters v,w in 
c     subroutine vhses.  
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are computed
c            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid point
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     ityp   = 0  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon.   
c
c            = 1  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 2  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 3  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 4  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 5  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 6  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 7  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 8  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c
c     nt     the number of analyses.  in the program that calls vhags,
c            the arrays v,w,br,bi,cr, and ci can be three dimensional
c            in which case multiple analyses will be performed.
c            the third index is the analysis index which assumes the 
c            values k=1,...,nt.  for a single analysis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that all the arrays are two
c            dimensional.
c
c     v,w    two or three dimensional arrays (see input parameter nt)
c            that contain the vector function to be analyzed.
c            v is the colatitudnal component and w is the east 
c            longitudinal component. v(i,j),w(i,j) contain the
c            components at the gaussian colatitude point theta(i)
c            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
c            are defined above at the input parameter ityp.
c
c     idvw   the first dimension of the arrays v,w as it appears in
c            the program that calls vhags. if ityp .le. 2 then idvw
c            must be at least nlat.  if ityp .gt. 2 and nlat is
c            even then idvw must be at least nlat/2. if ityp .gt. 2
c            and nlat is odd then idvw must be at least (nlat+1)/2.
c
c     jdvw   the second dimension of the arrays v,w as it appears in
c            the program that calls vhags. jdvw must be at least nlon.
c
c     mdab   the first dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhags. mdab must be at
c            least min0(nlat,nlon/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     ndab   the second dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhags. ndab must be at
c            least nlat.
c
c     wvhags an array which must be initialized by subroutine vhgsi.
c            once initialized, wvhags can be used repeatedly by vhags
c            as long as nlon and nlat remain unchanged.  wvhags must
c            not be altered between calls of vhags.
c
c     lvhags the dimension of the array wvhags as it appears in the
c            program that calls vhags. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhags must be at least
c
c            l1*l2(nlat+nlat-l1+1)+nlon+15
c
c        ??? (nlat+1)*(nlat+1)*nlat/2 + nlon + 15
c
c
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhags. define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            if ityp .le. 2 then lwork must be at least
c            the larger of the two quantities
c
c               3*nlat*(nlat+1)+2  (required by vhagsi)
c
c            and
c
c               (2*nt+1)*nlat*nlon
c
c            if ityp .gt. 2 then lwork must be at least
c            the larger of the two quantities
c
c               3*nlat*(nlat+1)+2  (required by vhagsi)
c
c            and
c
c              (2*nt+1)*l2*nlon
c
c
c     **************************************************************
c
c     output parameters
c
c     br,bi  two or three dimensional arrays (see input parameter nt)
c     cr,ci  that contain the vector spherical harmonic coefficients
c            in the spectral representation of v(i,j) and w(i,j) given 
c            in the discription of subroutine vhses. br(mp1,np1),
c            bi(mp1,np1),cr(mp1,np1), and ci(mp1,np1) are computed 
c            for mp1=1,...,mmax and np1=mp1,...,nlat except for np1=nlat
c            and odd mp1. mmax=min0(nlat,nlon/2) if nlon is even or 
c            mmax=min0(nlat,(nlon+1)/2) if nlon is odd. 
c      
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of ityp
c            = 4  error in the specification of nt
c            = 5  error in the specification of idvw
c            = 6  error in the specification of jdvw
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lvhags
c            = 10 error in the specification of lwork
c
c
c     subroutine vhagsi(nlat,nlon,wvhags,lvhags,work,lwork,ierror)
c
c     subroutine vhagsi initializes the array wvhags which can then be
c     used repeatedly by subroutine vhags until nlat or nlon is changed.
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are computed
c            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid point
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     lvhags the dimension of the array wvhags as it appears in the
c            program that calls vhagsi.  define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhags must be at least
c
c               3*nlat*(nlat+1)+2  (required by vhagsi)
c
c     dwork  a double precision work space that does not need to be saved
c
c     ldwork the dimension of the array dwork as it appears in the
c            program that calls vhagsi. ldwork must be at least
c
c                   (3*nlat*(nlat+3)+2)/2
c
c     **************************************************************
c
c     output parameters
c
c     wvhags an array which is initialized for use by subroutine vhags.
c            once initialized, wvhags can be used repeatedly by vhags
c            as long as nlat and nlon remain unchanged.  wvhags must not
c            be altered between calls of vhags.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lvhags
c            = 4  error in the specification of ldwork
c
      subroutine vhags(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
     1           mdab,ndab,wvhags,lvhags,work,lwork,ierror)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          work(1),wvhags(1)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      if(ityp.lt.0 .or. ityp.gt.8) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      imid = (nlat+1)/2
      if((ityp.le.2 .and. idvw.lt.nlat) .or.
     1   (ityp.gt.2 .and. idvw.lt.imid)) return
      ierror = 6
      if(jdvw .lt. nlon) return
      ierror = 7
      mmax = min0(nlat,(nlon+1)/2)
      if(mdab .lt. mmax) return
      ierror = 8
      if(ndab .lt. nlat) return
      ierror = 9
      idz = (mmax*(nlat+nlat-mmax+1))/2
      lzimn = idz*imid
      if(lvhags .lt. lzimn+lzimn+nlon+15) return
      ierror = 10
      idv = nlat
      if(ityp .gt. 2) idv = imid
      lnl = nt*idv*nlon
      if(lwork .lt. lnl+lnl+idv*nlon) return
      ierror = 0
      ist = 0
      if(ityp .le. 2) ist = imid
c
c     set wvhags pointers
c
      lmn = nlat*(nlat+1)/2
      jw1 = 1
      jw2 = jw1+imid*lmn
      jw3 = jw2+imid*lmn
c
c     set work pointers
c
      iw1 = ist+1
      iw2 = lnl+1
      iw3 = iw2+ist
      iw4 = iw2+lnl
      call vhags1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,ndab,
     +     br,bi,cr,ci,idv,work,work(iw1),work(iw2),work(iw3),
     +     work(iw4),idz,wvhags(jw1),wvhags(jw2),wvhags(jw3))
      return
      end

      subroutine vhags1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,
     +ndab,br,bi,cr,ci,idv,ve,vo,we,wo,work,idz,vb,wb,wrfft)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          ve(idv,nlon,1),vo(idv,nlon,1),we(idv,nlon,1),
     3          wo(idv,nlon,1),work(1),
     4          vb(imid,1),wb(imid,1),wrfft(1)
      nlp1 = nlat+1
      tsn = 2./nlon
      fsn = 4./nlon
      mlat = mod(nlat,2)
      mlon = mod(nlon,2)
      mmax = min0(nlat,(nlon+1)/2)
      imm1 = imid
      if(mlat .ne. 0) imm1 = imid-1
      if(ityp .gt. 2) go to 3  
      do 5 k=1,nt 
      do 5 i=1,imm1
      do 5 j=1,nlon
      ve(i,j,k) = tsn*(v(i,j,k)+v(nlp1-i,j,k))
      vo(i,j,k) = tsn*(v(i,j,k)-v(nlp1-i,j,k))
      we(i,j,k) = tsn*(w(i,j,k)+w(nlp1-i,j,k))
      wo(i,j,k) = tsn*(w(i,j,k)-w(nlp1-i,j,k))
    5 continue
      go to 2
    3 do 8 k=1,nt
      do 8 i=1,imm1 
      do 8 j=1,nlon
      ve(i,j,k) = fsn*v(i,j,k)
      vo(i,j,k) = fsn*v(i,j,k)
      we(i,j,k) = fsn*w(i,j,k)
      wo(i,j,k) = fsn*w(i,j,k)
    8 continue
    2 if(mlat .eq. 0) go to 7
      do 6 k=1,nt 
      do 6 j=1,nlon
      ve(imid,j,k) = tsn*v(imid,j,k)
      we(imid,j,k) = tsn*w(imid,j,k)
    6 continue
    7 do 9 k=1,nt
      call hrfftf(idv,nlon,ve(1,1,k),idv,wrfft,work)
      call hrfftf(idv,nlon,we(1,1,k),idv,wrfft,work)
    9 continue 
      ndo1 = nlat
      ndo2 = nlat
      if(mlat .ne. 0) ndo1 = nlat-1
      if(mlat .eq. 0) ndo2 = nlat-1
      if(ityp.eq.2 .or. ityp.eq.5 .or. ityp.eq.8) go to 11 
      do 10 k=1,nt
      do 10 mp1=1,mmax
      do 10 np1=mp1,nlat
      br(mp1,np1,k)=0.
      bi(mp1,np1,k)=0.
   10 continue
   11 if(ityp.eq.1 .or. ityp.eq.4 .or. ityp.eq.7) go to 13 
      do 12 k=1,nt
      do 12 mp1=1,mmax
      do 12 np1=mp1,nlat
      cr(mp1,np1,k)=0.
      ci(mp1,np1,k)=0.
   12 continue
   13 itypp = ityp+1
      go to (1,100,200,300,400,500,600,700,800),itypp
c
c     case ityp=0 ,  no symmetries
c
c     case m=0
c
    1 do 15 k=1,nt
      do 15 i=1,imid
      do 15 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+vb(i,np1)*ve(i,1,k)
      cr(1,np1,k) = cr(1,np1,k)-vb(i,np1)*we(i,1,k)
   15 continue
      do 16 k=1,nt
      do 16 i=1,imm1
      do 16 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+vb(i,np1)*vo(i,1,k)
      cr(1,np1,k) = cr(1,np1,k)-vb(i,np1)*wo(i,1,k)
   16 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 20 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 17
      do 23 k=1,nt
      do 23 i=1,imm1
      do 23 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(i,np1+mb)*vo(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(i,np1+mb)*vo(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*we(i,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(i,np1+mb)*wo(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(i,np1+mb)*wo(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*ve(i,2*mp1-2,k)
   23 continue
      if(mlat .eq. 0) go to 17
      do 24 k=1,nt
      do 24 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+wb(imid,np1+mb)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-wb(imid,np1+mb)*we(imid,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)+wb(imid,np1+mb)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-wb(imid,np1+mb)*ve(imid,2*mp1-2,k)
   24 continue
   17 if(mp2 .gt. ndo2) go to 20
      do 21 k=1,nt
      do 21 i=1,imm1
      do 21 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(i,np1+mb)*ve(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(i,np1+mb)*ve(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*wo(i,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(i,np1+mb)*we(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(i,np1+mb)*we(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*vo(i,2*mp1-2,k)
   21 continue
      if(mlat .eq. 0) go to 20
      do 22 k=1,nt
      do 22 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(imid,np1+mb)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(imid,np1+mb)*ve(imid,2*mp1-1,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(imid,np1+mb)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(imid,np1+mb)*we(imid,2*mp1-1,k)
   22 continue
   20 continue
      return
c
c     case ityp=1 ,  no symmetries but cr and ci equal zero
c
c     case m=0
c
  100 do 115 k=1,nt
      do 115 i=1,imid
      do 115 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+vb(i,np1)*ve(i,1,k)
  115 continue
      do 116 k=1,nt
      do 116 i=1,imm1
      do 116 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+vb(i,np1)*vo(i,1,k)
  116 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 120 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 117
      do 123 k=1,nt
      do 123 i=1,imm1
      do 123 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(i,np1+mb)*vo(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(i,np1+mb)*vo(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*we(i,2*mp1-2,k)
  123 continue
      if(mlat .eq. 0) go to 117
      do 124 k=1,nt
      do 124 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+wb(imid,np1+mb)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-wb(imid,np1+mb)*we(imid,2*mp1-2,k)
  124 continue
  117 if(mp2 .gt. ndo2) go to 120
      do 121 k=1,nt
      do 121 i=1,imm1
      do 121 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(i,np1+mb)*ve(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(i,np1+mb)*ve(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*wo(i,2*mp1-2,k)
  121 continue
      if(mlat .eq. 0) go to 120
      do 122 k=1,nt
      do 122 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(imid,np1+mb)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(imid,np1+mb)*ve(imid,2*mp1-1,k)
  122 continue
  120 continue
      return
c
c     case ityp=2 ,  no symmetries but br and bi equal zero   
c
c     case m=0
c
  200 do 215 k=1,nt
      do 215 i=1,imid
      do 215 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-vb(i,np1)*we(i,1,k)
  215 continue
      do 216 k=1,nt
      do 216 i=1,imm1
      do 216 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-vb(i,np1)*wo(i,1,k)
  216 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 220 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 217
      do 223 k=1,nt
      do 223 i=1,imm1
      do 223 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(i,np1+mb)*wo(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(i,np1+mb)*wo(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*ve(i,2*mp1-2,k)
  223 continue
      if(mlat .eq. 0) go to 217
      do 224 k=1,nt
      do 224 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+wb(imid,np1+mb)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-wb(imid,np1+mb)*ve(imid,2*mp1-2,k)
  224 continue
  217 if(mp2 .gt. ndo2) go to 220
      do 221 k=1,nt
      do 221 i=1,imm1
      do 221 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(i,np1+mb)*we(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(i,np1+mb)*we(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*vo(i,2*mp1-2,k)
  221 continue
      if(mlat .eq. 0) go to 220
      do 222 k=1,nt
      do 222 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(imid,np1+mb)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(imid,np1+mb)*we(imid,2*mp1-1,k)
  222 continue
  220 continue
      return
c
c     case ityp=3 ,  v even , w odd
c
c     case m=0
c
  300 do 315 k=1,nt
      do 315 i=1,imid
      do 315 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+vb(i,np1)*ve(i,1,k)
  315 continue
      do 316 k=1,nt
      do 316 i=1,imm1
      do 316 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-vb(i,np1)*wo(i,1,k)
  316 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 320 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 317
      do 323 k=1,nt
      do 323 i=1,imm1
      do 323 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(i,np1+mb)*wo(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(i,np1+mb)*wo(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*ve(i,2*mp1-2,k)
  323 continue
      if(mlat .eq. 0) go to 317
      do 324 k=1,nt
      do 324 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+wb(imid,np1+mb)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-wb(imid,np1+mb)*ve(imid,2*mp1-2,k)
  324 continue
  317 if(mp2 .gt. ndo2) go to 320
      do 321 k=1,nt
      do 321 i=1,imm1
      do 321 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(i,np1+mb)*ve(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(i,np1+mb)*ve(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*wo(i,2*mp1-2,k)
  321 continue
      if(mlat .eq. 0) go to 320
      do 322 k=1,nt
      do 322 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(imid,np1+mb)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(imid,np1+mb)*ve(imid,2*mp1-1,k)
  322 continue
  320 continue
      return
c
c     case ityp=4 ,  v even, w odd, and cr and ci equal 0. 
c
c     case m=0
c
  400 do 415 k=1,nt
      do 415 i=1,imid
      do 415 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+vb(i,np1)*ve(i,1,k)
  415 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 420 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 420
      do 421 k=1,nt
      do 421 i=1,imm1
      do 421 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(i,np1+mb)*ve(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(i,np1+mb)*ve(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*wo(i,2*mp1-2,k)
  421 continue
      if(mlat .eq. 0) go to 420
      do 422 k=1,nt
      do 422 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(imid,np1+mb)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(imid,np1+mb)*ve(imid,2*mp1-1,k)
  422 continue
  420 continue
      return
c
c     case ityp=5   v even, w odd, and br and bi equal zero
c
c     case m=0
c
  500 do 516 k=1,nt
      do 516 i=1,imm1
      do 516 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-vb(i,np1)*wo(i,1,k)
  516 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 520 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 520
      do 523 k=1,nt
      do 523 i=1,imm1
      do 523 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(i,np1+mb)*wo(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(i,np1+mb)*wo(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*ve(i,2*mp1-2,k)
  523 continue
      if(mlat .eq. 0) go to 520
      do 524 k=1,nt
      do 524 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+wb(imid,np1+mb)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-wb(imid,np1+mb)*ve(imid,2*mp1-2,k)
  524 continue
  520 continue
      return
c
c     case ityp=6 ,  v odd , w even
c
c     case m=0
c
  600 do 615 k=1,nt
      do 615 i=1,imid
      do 615 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-vb(i,np1)*we(i,1,k)
  615 continue
      do 616 k=1,nt
      do 616 i=1,imm1
      do 616 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+vb(i,np1)*vo(i,1,k)
  616 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 620 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 617
      do 623 k=1,nt
      do 623 i=1,imm1
      do 623 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(i,np1+mb)*vo(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(i,np1+mb)*vo(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*we(i,2*mp1-2,k)
  623 continue
      if(mlat .eq. 0) go to 617
      do 624 k=1,nt
      do 624 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+wb(imid,np1+mb)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-wb(imid,np1+mb)*we(imid,2*mp1-2,k)
  624 continue
  617 if(mp2 .gt. ndo2) go to 620
      do 621 k=1,nt
      do 621 i=1,imm1
      do 621 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(i,np1+mb)*we(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(i,np1+mb)*we(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*vo(i,2*mp1-2,k)
  621 continue
      if(mlat .eq. 0) go to 620
      do 622 k=1,nt
      do 622 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(imid,np1+mb)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(imid,np1+mb)*we(imid,2*mp1-1,k)
  622 continue
  620 continue
      return
c
c     case ityp=7   v odd, w even, and cr and ci equal zero
c
c     case m=0
c
  700 do 716 k=1,nt
      do 716 i=1,imm1
      do 716 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+vb(i,np1)*vo(i,1,k)
  716 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 720 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 720
      do 723 k=1,nt
      do 723 i=1,imm1
      do 723 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+vb(i,np1+mb)*vo(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+vb(i,np1+mb)*vo(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*we(i,2*mp1-2,k)
  723 continue
      if(mlat .eq. 0) go to 720
      do 724 k=1,nt
      do 724 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+wb(imid,np1+mb)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-wb(imid,np1+mb)*we(imid,2*mp1-2,k)
  724 continue
  720 continue
      return
c
c     case ityp=8   v odd, w even, and both br and bi equal zero
c
c     case m=0
c
  800 do 815 k=1,nt
      do 815 i=1,imid
      do 815 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-vb(i,np1)*we(i,1,k)
  815 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 820 mp1=2,mmax
      m = mp1-1
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 820
      do 821 k=1,nt
      do 821 i=1,imm1
      do 821 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(i,np1+mb)*we(i,2*mp1-2,k)
     1                             +wb(i,np1+mb)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(i,np1+mb)*we(i,2*mp1-1,k)
     1                             -wb(i,np1+mb)*vo(i,2*mp1-2,k)
  821 continue
      if(mlat .eq. 0) go to 820
      do 822 k=1,nt
      do 822 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-vb(imid,np1+mb)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-vb(imid,np1+mb)*we(imid,2*mp1-1,k)
  822 continue
  820 continue
      return
      end
      subroutine vhagsi(nlat,nlon,wvhags,lvhags,dwork,ldwork,ierror)
      dimension wvhags(lvhags)
      double precision dwork(ldwork)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      imid = (nlat+1)/2
      lmn = (nlat*(nlat+1))/2
      if(lvhags .lt. 2*(imid*lmn)+nlon+15) return
      ierror = 4
c     if (ldwork.lt.nlat*(3*nlat+9)+2) return
      if (ldwork.lt.(nlat*(3*nlat+9)+2)/2) return
      ierror = 0
      jw1 = 1
      jw2 = jw1+imid*lmn
      jw3 = jw2+imid*lmn
      iw1 = 1
      iw2 = iw1+nlat
      iw3 = iw2+nlat
      iw4 = iw3+3*imid*nlat
c     iw2 = iw1+nlat+nlat
c     iw3 = iw2+nlat+nlat
c     iw4 = iw3+6*imid*nlat
      call vhgai1(nlat,imid,wvhags(jw1),wvhags(jw2),
     +dwork(iw1),dwork(iw2),dwork(iw3),dwork(iw4))
      call hrffti(nlon,wvhags(jw3))
      return
      end
      subroutine vhgai1(nlat,imid,vb,wb,dthet,dwts,dpbar,work)
      dimension vb(imid,*),wb(imid,*)
      double precision abel,bbel,cbel,ssqr2,dcf
      double precision dpbar(imid,nlat,3), dthet(*),dwts(*),work(*)
c     lwk = 4*nlat*(nlat+2)
      lwk = nlat*(nlat+2)
      call gaqd(nlat,dthet,dwts,dpbar,lwk,ierror)
c
c     compute associated legendre functions
c
c     compute m=n=0 legendre polynomials for all theta(i)
c
      ssqr2 = 1./dsqrt(2.d0)
      do 90 i=1,imid
      dpbar(i,1,1) = ssqr2
      vb(i,1) = 0.
      wb(i,1) = 0.
   90 continue
c
c     main loop for remaining vb, and wb
c
      do 100 n=1,nlat-1
      nm = mod(n-2,3)+1
      nz = mod(n-1,3)+1
      np = mod(n,3)+1
c
c     compute dpbar for m=0
c
      call dnlfk(0,n,work)
      mn = indx(0,n,nlat)
      do 105 i=1,imid
      call dnlft(0,n,dthet(i),work,dpbar(i,1,np))
  105 continue
c
c     compute dpbar for m=1
c
      call dnlfk(1,n,work)
      mn = indx(1,n,nlat)
      do 106 i=1,imid
      call dnlft(1,n,dthet(i),work,dpbar(i,2,np))
c      pbar(i,mn) = dpbar(i,2,np)
  106 continue
  104 continue
c
c     compute and store dpbar for m=2,n
c
      if(n.lt.2) go to 108
      do 107 m=2,n
      abel = dsqrt(dble(float((2*n+1)*(m+n-2)*(m+n-3)))/
     1                dble(float((2*n-3)*(m+n-1)*(m+n))))
      bbel = dsqrt(dble(float((2*n+1)*(n-m-1)*(n-m)))/
     1                dble(float((2*n-3)*(m+n-1)*(m+n))))
      cbel = dsqrt(dble(float((n-m+1)*(n-m+2)))/
     1                dble(float((m+n-1)*(m+n))))
      id = indx(m,n,nlat)
      if (m.ge.n-1) go to 102
      do 103 i=1,imid
      dpbar(i,m+1,np) = abel*dpbar(i,m-1,nm)+bbel*dpbar(i,m+1,nm)
     1                                         -cbel*dpbar(i,m-1,np)
  103 continue
      go to 107
  102 do 101 i=1,imid
      dpbar(i,m+1,np) = abel*dpbar(i,m-1,nm)-cbel*dpbar(i,m-1,np)
  101 continue
  107 continue
c
c     compute the derivative of the functions
c
  108 continue
      ix = indx(0,n,nlat)
      iy = indx(n,n,nlat)
      do 125 i=1,imid
      vb(i,ix) = -dpbar(i,2,np)*dwts(i)
      vb(i,iy) = dpbar(i,n,np)/dsqrt(dble(float(2*(n+1))))*dwts(i)
125   continue
c
      if(n.eq.1) go to 131 
      dcf = dsqrt(dble(float(4*n*(n+1))))
      do 130 m=1,n-1
      ix = indx(m,n,nlat)
      abel = dsqrt(dble(float((n+m)*(n-m+1))))/dcf
      bbel = dsqrt(dble(float((n-m)*(n+m+1))))/dcf
      do 130 i=1,imid
      vb(i,ix) = (abel*dpbar(i,m,np)-bbel*dpbar(i,m+2,np))*dwts(i)
130   continue
c
c     compute the vector harmonic w(theta) = m*pbar/cos(theta)
c
c     set wb=0 for m=0
c
  131 continue
      ix = indx(0,n,nlat)
      do 220 i=1,imid
      wb(i,ix) = 0.d0
220   continue
c
c     compute wb for m=1,n
c
      dcf = dsqrt(dble(float(n+n+1))/dble(float(4*n*(n+1)*(n+n-1))))
      do 230 m=1,n
      ix = indx(m,n,nlat)
      abel = dcf*dsqrt(dble(float((n+m)*(n+m-1))))
      bbel = dcf*dsqrt(dble(float((n-m)*(n-m-1))))
      if(m.ge.n-1) go to 231
      do 229 i=1,imid
      wb(i,ix) = (abel*dpbar(i,m,nz) + bbel*dpbar(i,m+2,nz))*dwts(i)
  229 continue
      go to 230
  231 do 228 i=1,imid
      wb(i,ix) = abel*dpbar(i,m,nz)*dwts(i)
  228 continue
  230 continue
  100 continue
      return 
      end

      function indx(m,n,nlat)
      integer indx
      indx = m*nlat-(m*(m+1))/2+n+1
      return
      end
 
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                         SPHEREPACK                          .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... file vhsgs.f
c
c     this file contains code and documentation for subroutines
c     vhsgs and vhsgsi
c
c ... files which must be loaded with vhsgs.f
c
c     sphcom.f, hrfft.f, gaqd.f
c
c     subroutine vhsgs(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
c    +                 mdab,ndab,wvhsgs,lvhsgs,work,lwork,ierror)
c                                                                              
c   
c     subroutine vhsgs performs the vector spherical harmonic synthesis
c     of the arrays br, bi, cr, and ci and stores the result in the
c     arrays v and w.  the synthesis is performed on an equally spaced
c     longitude grid and a gaussian colatitude grid (measured from
c     the north pole). v(i,j) and w(i,j) are the colatitudinal and
c     east longitudinal components respectively, located at the i(th)
c     colatitude gaussian point (see nlat below) and longitude
c     phi(j) = (j-1)*2*pi/nlon.  the spectral respresentation of (v,w)
c     is given below at output parameters v,w.
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are computed
c            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid point
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     ityp   = 0  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon.   
c
c            = 1  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 2  no symmetries exist about the equator. the synthesis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 3  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 4  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 5  v is symmetric and w is antisymmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 6  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 7  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 8  v is antisymmetric and w is symmetric about the 
c                 equator. the synthesis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the synthesis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the synthesis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c
c     nt     the number of syntheses.  in the program that calls vhsgs,
c            the arrays v,w,br,bi,cr, and ci can be three dimensional
c            in which case multiple syntheses will be performed.
c            the third index is the synthesis index which assumes the 
c            values k=1,...,nt.  for a single synthesis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that all the arrays are two
c            dimensional.
c
c     idvw   the first dimension of the arrays v,w as it appears in
c            the program that calls vhags. if ityp .le. 2 then idvw
c            must be at least nlat.  if ityp .gt. 2 and nlat is
c            even then idvw must be at least nlat/2. if ityp .gt. 2
c            and nlat is odd then idvw must be at least (nlat+1)/2.
c
c     jdvw   the second dimension of the arrays v,w as it appears in
c            the program that calls vhsgs. jdvw must be at least nlon.
c
c     br,bi  two or three dimensional arrays (see input parameter nt)
c     cr,ci  that contain the vector spherical harmonic coefficients
c            in the spectral representation of v(i,j) and w(i,j) given
c            below at the discription of output parameters v and w.
c
c     mdab   the first dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhsgs. mdab must be at
c            least min0(nlat,nlon/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     ndab   the second dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhsgs. ndab must be at
c            least nlat.
c
c     wvhsgs an array which must be initialized by subroutine vhsgsi.
c            once initialized, wvhsgs can be used repeatedly by vhsgs
c            as long as nlon and nlat remain unchanged.  wvhsgs must
c            not be altered between calls of vhsgs.
c
c     lvhsgs the dimension of the array wvhsgs as it appears in the
c            program that calls vhsgs. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhsgs must be at least
c
c                 l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat
c
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhsgs. define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            if ityp .le. 2 then lwork must be at least
c
c                       (2*nt+1)*nlat*nlon
c
c            if ityp .gt. 2 then lwork must be at least
c
c                        (2*nt+1)*l2*nlon 
c
c     **************************************************************
c
c     output parameters
c
c     v,w    two or three dimensional arrays (see input parameter nt)
c            in which the synthesis is stored. v is the colatitudinal
c            component and w is the east longitudinal component. 
c            v(i,j),w(i,j) contain the components at the guassian colatitude
c            point theta(i) and longitude phi(j) = (j-1)*2*pi/nlon.
c            the index ranges are defined above at the input parameter
c            ityp. v and w are computed from the formulas given below.
c
c
c     define
c
c     1.  theta is colatitude and phi is east longitude
c
c     2.  the normalized associated legendre funnctions
c
c         pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)
c                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
c                        factorial(n)) times the (n+m)th derivative
c                        of (x**2-1)**n with respect to x=cos(theta)
c
c     3.  vbar(m,n,theta) = the derivative of pbar(m,n,theta) with
c                           respect to theta divided by the square
c                           root of n(n+1).
c
c         vbar(m,n,theta) is more easily computed in the form
c
c         vbar(m,n,theta) = (sqrt((n+m)*(n-m+1))*pbar(m-1,n,theta)
c         -sqrt((n-m)*(n+m+1))*pbar(m+1,n,theta))/(2*sqrt(n*(n+1)))
c
c     4.  wbar(m,n,theta) = m/(sin(theta))*pbar(m,n,theta) divided
c                           by the square root of n(n+1).
c
c         wbar(m,n,theta) is more easily computed in the form
c
c         wbar(m,n,theta) = sqrt((2n+1)/(2n-1))*(sqrt((n+m)*(n+m-1))
c         *pbar(m-1,n-1,theta)+sqrt((n-m)*(n-m-1))*pbar(m+1,n-1,theta))
c         /(2*sqrt(n*(n+1)))
c
c
c    the colatitudnal dependence of the normalized surface vector
c                spherical harmonics are defined by
c
c     5.    bbar(m,n,theta) = (vbar(m,n,theta),i*wbar(m,n,theta))
c
c     6.    cbar(m,n,theta) = (i*wbar(m,n,theta),-vbar(m,n,theta))
c
c
c    the coordinate to index mappings 
c
c     7.   theta(i) = i(th) gaussian grid point and phi(j) = (j-1)*2*pi/nlon
c
c    
c     the maximum (plus one) longitudinal wave number
c
c     8.     mmax = min0(nlat,nlon/2) if nlon is even or
c            mmax = min0(nlat,(nlon+1)/2) if nlon is odd.
c
c    if we further define the output vector as
c
c     9.    h(i,j) = (v(i,j),w(i,j))
c
c    and the complex coefficients
c
c     10.   b(m,n) = cmplx(br(m+1,n+1),bi(m+1,n+1))
c
c     11.   c(m,n) = cmplx(cr(m+1,n+1),ci(m+1,n+1))
c
c
c    then for i=1,...,nlat and  j=1,...,nlon
c
c        the expansion for real h(i,j) takes the form
c
c     h(i,j) = the sum from n=1 to n=nlat-1 of the real part of
c
c         .5*(b(0,n)*bbar(0,n,theta(i))+c(0,n)*cbar(0,n,theta(i)))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
c     n=nlat-1 of the real part of
c
c              b(m,n)*bbar(m,n,theta(i))*exp(i*m*phi(j))
c             +c(m,n)*cbar(m,n,theta(i))*exp(i*m*phi(j))
c
c   *************************************************************
c
c   in terms of real variables this expansion takes the form
c
c             for i=1,...,nlat and  j=1,...,nlon
c
c     v(i,j) = the sum from n=1 to n=nlat-1 of
c
c               .5*br(1,n+1)*vbar(0,n,theta(i))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
c     n=nlat-1 of the real part of
c
c       (br(m+1,n+1)*vbar(m,n,theta(i))-ci(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *cos(m*phi(j))
c      -(bi(m+1,n+1)*vbar(m,n,theta(i))+cr(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *sin(m*phi(j))
c
c    and for i=1,...,nlat and  j=1,...,nlon
c
c     w(i,j) = the sum from n=1 to n=nlat-1 of
c
c              -.5*cr(1,n+1)*vbar(0,n,theta(i))
c
c     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
c     n=nlat-1 of the real part of
c
c      -(cr(m+1,n+1)*vbar(m,n,theta(i))+bi(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *cos(m*phi(j))
c      +(ci(m+1,n+1)*vbar(m,n,theta(i))-br(m+1,n+1)*wbar(m,n,theta(i)))
c                                          *sin(m*phi(j))
c
c
c      br(m+1,nlat),bi(m+1,nlat),cr(m+1,nlat), and ci(m+1,nlat) are
c      assumed zero for m even.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of ityp
c            = 4  error in the specification of nt
c            = 5  error in the specification of idvw
c            = 6  error in the specification of jdvw
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lvhsgs
c            = 10 error in the specification of lwork
c
c
c     subroutine vhsgsi(nlat,nlon,wvhsgs,lvhsgs,dwork,ldwork,ierror)
c
c     subroutine vhsgsi initializes the array wvhsgs which can then be
c     used repeatedly by subroutine vhsgs until nlat or nlon is changed.
c
c     input parameters
c
c     nlat   the number of points in the gaussian colatitude grid on the
c            full sphere. these lie in the interval (0,pi) and are computed
c            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
c            if nlat is odd the equator will be included as the grid point
c            theta((nlat+1)/2).  if nlat is even the equator will be
c            excluded as a grid point and will lie half way between
c            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
c            note: on the half sphere, the number of grid points in the
c            colatitudinal direction is nlat/2 if nlat is even or
c            (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     lvhsgs the dimension of the array wvhsgs as it appears in the
c            program that calls vhsgs. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhsgs must be at least
c
c                 l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat
c
c     dwork a double precision work array that does not need to be saved
c
c     ldwork the dimension of the array dwork as it appears in the
c            program that calls vhsgsi. ldwork must be at least
c
c                 (3*nlat*(nlat+3)+2)/2

c
c     **************************************************************
c
c     output parameters
c
c     wvhsgs an array which is initialized for use by subroutine vhsgs.
c            once initialized, wvhsgs can be used repeatedly by vhsgs
c            as long as nlat and nlon remain unchanged.  wvhsgs must not
c            be altered between calls of vhsgs.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lvhsgs
c            = 4  error in the specification of lwork
c
      subroutine vhsgs(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
     +                 mdab,ndab,wvhsgs,lvhsgs,work,lwork,ierror)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          work(1),wvhsgs(1)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      if(ityp.lt.0 .or. ityp.gt.8) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      imid = (nlat+1)/2
      if((ityp.le.2 .and. idvw.lt.nlat) .or.
     1   (ityp.gt.2 .and. idvw.lt.imid)) return
      ierror = 6
      if(jdvw .lt. nlon) return
      ierror = 7
      mmax = min0(nlat,(nlon+1)/2)
      if(mdab .lt. mmax) return
      ierror = 8
      if(ndab .lt. nlat) return
      ierror = 9
      idz = (mmax*(nlat+nlat-mmax+1))/2
      lzimn = idz*imid
      if(lvhsgs .lt. lzimn+lzimn+nlon+15) return
      ierror = 10
      idv = nlat
      if(ityp .gt. 2) idv = imid
      lnl = nt*idv*nlon
      if(lwork .lt. lnl+lnl+idv*nlon) return
      ierror = 0
      ist = 0
      if(ityp .le. 2) ist = imid
c
c     set wvhsgs pointers
c
      lmn = nlat*(nlat+1)/2
      jw1 = 1
      jw2 = jw1+imid*lmn
      jw3 = jw2+imid*lmn
c
c     set work pointers
c
      iw1 = ist+1
      iw2 = lnl+1
      iw3 = iw2+ist
      iw4 = iw2+lnl

      call vhsgs1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,ndab,
     +           br,bi,cr,ci,idv,work,work(iw1),work(iw2),work(iw3),
     +          work(iw4),idz,wvhsgs(jw1),wvhsgs(jw2),wvhsgs(jw3))
      return
      end

      subroutine vhsgs1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,
     1   ndab,br,bi,cr,ci,idv,ve,vo,we,wo,work,idz,vb,wb,wrfft)
      dimension v(idvw,jdvw,1),w(idvw,jdvw,1),br(mdab,ndab,1),
     1          bi(mdab,ndab,1),cr(mdab,ndab,1),ci(mdab,ndab,1),
     2          ve(idv,nlon,1),vo(idv,nlon,1),we(idv,nlon,1),
     3          wo(idv,nlon,1),work(1),wrfft(1),
     4          vb(imid,1),wb(imid,1)
      nlp1 = nlat+1
      mlat = mod(nlat,2)
      mlon = mod(nlon,2)
      mmax = min0(nlat,(nlon+1)/2)
      imm1 = imid
      if(mlat .ne. 0) imm1 = imid-1
      do 10 k=1,nt
      do 10 j=1,nlon
      do 10 i=1,idv
      ve(i,j,k) = 0.
      we(i,j,k) = 0.
   10 continue
      ndo1 = nlat
      ndo2 = nlat
      if(mlat .ne. 0) ndo1 = nlat-1
      if(mlat .eq. 0) ndo2 = nlat-1
   18 itypp = ityp+1
      go to (1,100,200,300,400,500,600,700,800),itypp
c
c     case ityp=0   no symmetries
c
c     case m = 0
c
    1 continue
      do 15 k=1,nt
      do 15 np1=2,ndo2,2
      do 15 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1)
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1)
   15 continue
      do 16 k=1,nt
      do 16 np1=3,ndo1,2
      do 16 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1)
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1)
   16 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 30 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 26
      do 25 k=1,nt
      do 24 np1=mp1,ndo1,2
      mn = mb+np1
      do 23 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
   23 continue
      if(mlat .eq. 0) go to 24
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,mn)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,mn)
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,mn) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,mn)
   24 continue
   25 continue
   26 if(mp2 .gt. ndo2) go to 30
      do 29 k=1,nt
      do 28 np1=mp2,ndo2,2
      mn = mb+np1
      do 27 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
   27 continue
      if(mlat .eq. 0) go to 28
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,mn) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,mn)
   28 continue
   29 continue
   30 continue
      go to 950
c
c     case ityp=1   no symmetries,  cr and ci equal zero
c
c     case m = 0
c
  100 continue
      do 115 k=1,nt
      do 115 np1=2,ndo2,2
      do 115 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1)
  115 continue
      do 116 k=1,nt
      do 116 np1=3,ndo1,2
      do 116 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1)
  116 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 130 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 126
      do 125 k=1,nt
      do 124 np1=mp1,ndo1,2
      mn = mb+np1
      do 123 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  123 continue
      if(mlat .eq. 0) go to 124
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,mn) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,mn)
  124 continue
  125 continue
  126 if(mp2 .gt. ndo2) go to 130
      do 129 k=1,nt
      do 128 np1=mp2,ndo2,2
      mn = mb+np1
      do 127 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  127 continue
      if(mlat .eq. 0) go to 128
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,mn) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,mn)
  128 continue
  129 continue
  130 continue
      go to 950
c
c     case ityp=2   no symmetries,  br and bi are equal to zero
c
c     case m = 0
c
  200 do 215 k=1,nt
      do 215 np1=2,ndo2,2
      do 215 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1)
  215 continue
      do 216 k=1,nt
      do 216 np1=3,ndo1,2
      do 216 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1)
  216 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 230 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 226
      do 225 k=1,nt
      do 224 np1=mp1,ndo1,2
      mn = mb+np1
      do 223 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  223 continue
      if(mlat .eq. 0) go to 224
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,mn)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,mn)
  224 continue
  225 continue
  226 if(mp2 .gt. ndo2) go to 230
      do 229 k=1,nt
      do 228 np1=mp2,ndo2,2
      mn = mb+np1
      do 227 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  227 continue
      if(mlat .eq. 0) go to 228
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,mn)
  228 continue
  229 continue
  230 continue
      go to 950
c
c     case ityp=3   v even,  w odd 
c
c     case m = 0
c
  300 do 315 k=1,nt
      do 315 np1=2,ndo2,2
      do 315 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1)
  315 continue
      do 316 k=1,nt
      do 316 np1=3,ndo1,2
      do 316 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1)
  316 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 330 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 326
      do 325 k=1,nt
      do 324 np1=mp1,ndo1,2
      mn = mb+np1
      do 323 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  323 continue
      if(mlat .eq. 0) go to 324
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,mn)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,mn)
  324 continue
  325 continue
  326 if(mp2 .gt. ndo2) go to 330
      do 329 k=1,nt
      do 328 np1=mp2,ndo2,2
      mn = mb+np1
      do 327 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  327 continue
      if(mlat .eq. 0) go to 328
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,mn) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,mn)
  328 continue
  329 continue
  330 continue
      go to 950
c
c     case ityp=4   v even,  w odd, and both cr and ci equal zero 
c
c     case m = 0
c
  400 do 415 k=1,nt
      do 415 np1=2,ndo2,2
      do 415 i=1,imid
      ve(i,1,k)=ve(i,1,k)+br(1,np1,k)*vb(i,np1)
  415 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 430 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 430
      do 429 k=1,nt
      do 428 np1=mp2,ndo2,2
      mn = mb+np1
      do 427 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  427 continue
      if(mlat .eq. 0) go to 428
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     +br(mp1,np1,k)*vb(imid,mn) 
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +bi(mp1,np1,k)*vb(imid,mn)
  428 continue
  429 continue
  430 continue
      go to 950
c
c     case ityp=5   v even,  w odd,     br and bi equal zero 
c
c     case m = 0
c
  500 do 516 k=1,nt
      do 516 np1=3,ndo1,2
      do 516 i=1,imm1
      wo(i,1,k)=wo(i,1,k)-cr(1,np1,k)*vb(i,np1)
  516 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 530 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 530
      do 525 k=1,nt
      do 524 np1=mp1,ndo1,2
      mn = mb+np1
      do 523 i=1,imm1
      ve(i,2*mp1-2,k) = ve(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      ve(i,2*mp1-1,k) = ve(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      wo(i,2*mp1-2,k) = wo(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      wo(i,2*mp1-1,k) = wo(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  523 continue
      if(mlat .eq. 0) go to 524
      ve(imid,2*mp1-2,k) = ve(imid,2*mp1-2,k)
     1                     -ci(mp1,np1,k)*wb(imid,mn)
      ve(imid,2*mp1-1,k) = ve(imid,2*mp1-1,k)
     1                     +cr(mp1,np1,k)*wb(imid,mn)
  524 continue
  525 continue
  530 continue
      go to 950
c
c     case ityp=6   v odd  ,  w even
c
c     case m = 0
c
  600 do 615 k=1,nt
      do 615 np1=2,ndo2,2
      do 615 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1)
  615 continue
      do 616 k=1,nt
      do 616 np1=3,ndo1,2
      do 616 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1)
  616 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 630 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 626
      do 625 k=1,nt
      do 624 np1=mp1,ndo1,2
      mn = mb+np1
      do 623 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  623 continue
      if(mlat .eq. 0) go to 624
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,mn) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,mn)
  624 continue
  625 continue
  626 if(mp2 .gt. ndo2) go to 630
      do 629 k=1,nt
      do 628 np1=mp2,ndo2,2
      mn = mb+np1
      do 627 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  627 continue
      if(mlat .eq. 0) go to 628
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,mn)
  628 continue
  629 continue
  630 continue
      go to 950
c
c     case ityp=7   v odd, w even   cr and ci equal zero
c
c     case m = 0
c
  700 do 716 k=1,nt
      do 716 np1=3,ndo1,2
      do 716 i=1,imm1
      vo(i,1,k)=vo(i,1,k)+br(1,np1,k)*vb(i,np1)
  716 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 730 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 730
      do 725 k=1,nt
      do 724 np1=mp1,ndo1,2
      mn = mb+np1
      do 723 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)+br(mp1,np1,k)*vb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+bi(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-bi(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)+br(mp1,np1,k)*wb(i,mn)
  723 continue
      if(mlat .eq. 0) go to 724
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -bi(mp1,np1,k)*wb(imid,mn) 
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     +br(mp1,np1,k)*wb(imid,mn)
  724 continue
  725 continue
  730 continue
      go to 950
c
c     case ityp=8   v odd,  w even   br and bi equal zero
c
c     case m = 0
c
  800 do 815 k=1,nt
      do 815 np1=2,ndo2,2
      do 815 i=1,imid
      we(i,1,k)=we(i,1,k)-cr(1,np1,k)*vb(i,np1)
  815 continue
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) go to 950
      do 830 mp1=2,mmax
      m = mp1-1
c     mb = m*(nlat-1)-(m*(m-1))/2
      mb = m*nlat-(m*(m+1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 830
      do 829 k=1,nt
      do 828 np1=mp2,ndo2,2
      mn = mb+np1
      do 827 i=1,imm1
      vo(i,2*mp1-2,k) = vo(i,2*mp1-2,k)-ci(mp1,np1,k)*wb(i,mn)
      vo(i,2*mp1-1,k) = vo(i,2*mp1-1,k)+cr(mp1,np1,k)*wb(i,mn)
      we(i,2*mp1-2,k) = we(i,2*mp1-2,k)-cr(mp1,np1,k)*vb(i,mn)
      we(i,2*mp1-1,k) = we(i,2*mp1-1,k)-ci(mp1,np1,k)*vb(i,mn)
  827 continue
      if(mlat .eq. 0) go to 828
      we(imid,2*mp1-2,k) = we(imid,2*mp1-2,k)
     1                     -cr(mp1,np1,k)*vb(imid,mn)
      we(imid,2*mp1-1,k) = we(imid,2*mp1-1,k)
     1                     -ci(mp1,np1,k)*vb(imid,mn)
  828 continue
  829 continue
  830 continue
  950 continue
      do 14 k=1,nt
      call hrfftb(idv,nlon,ve(1,1,k),idv,wrfft,work)
      call hrfftb(idv,nlon,we(1,1,k),idv,wrfft,work)
   14 continue
      if(ityp .gt. 2) go to 12
      do 60 k=1,nt
      do 60 j=1,nlon
      do 60 i=1,imm1
      v(i,j,k) = .5*(ve(i,j,k)+vo(i,j,k))
      w(i,j,k) = .5*(we(i,j,k)+wo(i,j,k))
      v(nlp1-i,j,k) = .5*(ve(i,j,k)-vo(i,j,k))
      w(nlp1-i,j,k) = .5*(we(i,j,k)-wo(i,j,k))
   60 continue
      go to 13
   12 do 11 k=1,nt
      do 11 j=1,nlon
      do 11 i=1,imm1
      v(i,j,k) = .5*ve(i,j,k)
      w(i,j,k) = .5*we(i,j,k)
   11 continue
   13 if(mlat .eq. 0) return
      do 65 k=1,nt
      do 65 j=1,nlon
      v(imid,j,k) = .5*ve(imid,j,k)
      w(imid,j,k) = .5*we(imid,j,k)
   65 continue
      return
      end
      subroutine vhsgsi(nlat,nlon,wvhsgs,lvhsgs,dwork,ldwork,ierror)
c
c     subroutine vhsfsi computes the gaussian points theta, gauss 
c     weights wts, and the components vb and wb of the vector 
c     harmonics. all quantities are computed internally in double 
c     precision but returned in single precision and are therfore 
c     accurate to single precision.
c
c     set imid = (nlat+1)/2 and lmn=(nlat*(nlat+1))/2 then
c     wvhsgs must have 2*(imid*lmn+nlat)+nlon+15 locations
c
c     double precision array dwork must have
c       3*nlat*(nlat+1)+5*nlat+1 = nlat*(3*nlat+8)+1 
c     locations which is determined by the size of dthet, 
c     dwts, dwork, and dpbar in vhsgs1
c
      dimension wvhsgs(*)
      double precision dwork(*)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      imid = (nlat+1)/2
      lmn = (nlat*(nlat+1))/2
      if(lvhsgs .lt. 2*(imid*lmn)+nlon+15) return
      ierror = 4
      if (ldwork .lt. (nlat*3*(nlat+3)+2)/2) return
      ierror = 0
c
c     set saved work space pointers
c
      jw1 = 1
      jw2 = jw1+imid*lmn
      jw3 = jw2+imid*lmn
c
c     set unsaved work space pointers
c
      iw1 = 1
      iw2 = iw1+nlat
      iw3 = iw2+nlat
      iw4 = iw3+3*imid*nlat
c     iw2 = iw1+nlat+nlat
c     iw3 = iw2+nlat+nlat
c     iw4 = iw3+6*imid*nlat
      call vhgsi1(nlat,imid,wvhsgs(jw1),wvhsgs(jw2),
     +dwork(iw1),dwork(iw2),dwork(iw3),dwork(iw4))
      call hrffti(nlon,wvhsgs(jw3))
      return
      end
      subroutine vhgsi1(nlat,imid,vb,wb,dthet,dwts,dpbar,work)
      dimension vb(imid,*),wb(imid,*)
      double precision abel,bbel,cbel,ssqr2,dcf
      double precision dthet(*),dwts(*),dpbar(imid,nlat,3),work(*)
c
c     compute gauss points and weights
c     use dpbar (length 3*nnlat*(nnlat+1)) as work space for gaqd
c
      lwk = nlat*(nlat+2)
      call gaqd(nlat,dthet,dwts,dpbar,lwk,ierror)
c
c     compute associated legendre functions
c
c     compute m=n=0 legendre polynomials for all theta(i)
c
      ssqr2 = 1./dsqrt(2.d0)
      do 90 i=1,imid
      dpbar(i,1,1) = ssqr2
      vb(i,1) = 0.
      wb(i,1) = 0.
   90 continue
c
c     main loop for remaining vb, and wb
c
      do 100 n=1,nlat-1
      nm = mod(n-2,3)+1
      nz = mod(n-1,3)+1
      np = mod(n,3)+1
c
c     compute dpbar for m=0
c
      call dnlfk(0,n,work)
      mn = indx(0,n,nlat)
      do 105 i=1,imid
      call dnlft(0,n,dthet(i),work,dpbar(i,1,np))
c      pbar(i,mn) = dpbar(i,1,np)
  105 continue
c
c     compute dpbar for m=1
c
      call dnlfk(1,n,work)
      mn = indx(1,n,nlat)
      do 106 i=1,imid
      call dnlft(1,n,dthet(i),work,dpbar(i,2,np))
c      pbar(i,mn) = dpbar(i,2,np)
  106 continue
  104 continue
c
c     compute and store dpbar for m=2,n
c
      if(n.lt.2) go to 108
      do 107 m=2,n
      abel = dsqrt(dble(float((2*n+1)*(m+n-2)*(m+n-3)))/
     1                dble(float((2*n-3)*(m+n-1)*(m+n))))
      bbel = dsqrt(dble(float((2*n+1)*(n-m-1)*(n-m)))/
     1                dble(float((2*n-3)*(m+n-1)*(m+n))))
      cbel = dsqrt(dble(float((n-m+1)*(n-m+2)))/
     1                dble(float((m+n-1)*(m+n))))
      id = indx(m,n,nlat)
      if (m.ge.n-1) go to 102
      do 103 i=1,imid
      dpbar(i,m+1,np) = abel*dpbar(i,m-1,nm)+bbel*dpbar(i,m+1,nm)
     1                                         -cbel*dpbar(i,m-1,np)
c      pbar(i,id) = dpbar(i,m+1,np)
  103 continue
      go to 107
  102 do 101 i=1,imid
      dpbar(i,m+1,np) = abel*dpbar(i,m-1,nm)-cbel*dpbar(i,m-1,np)
c      pbar(i,id) = dpbar(i,m+1,np)
  101 continue
  107 continue
c
c     compute the derivative of the functions
c
  108 ix = indx(0,n,nlat)
      iy = indx(n,n,nlat)
      do 125 i=1,imid
      vb(i,ix) = -dpbar(i,2,np)
      vb(i,iy) = dpbar(i,n,np)/dsqrt(dble(float(2*(n+1))))
125   continue
c
      if(n.eq.1) go to 131 
      dcf = dsqrt(dble(float(4*n*(n+1))))
      do 130 m=1,n-1
      ix = indx(m,n,nlat)     
      abel = dsqrt(dble(float((n+m)*(n-m+1))))/dcf
      bbel = dsqrt(dble(float((n-m)*(n+m+1))))/dcf
      do 130 i=1,imid
      vb(i,ix) = abel*dpbar(i,m,np)-bbel*dpbar(i,m+2,np)
130   continue
c
c     compute the vector harmonic w(theta) = m*pbar/cos(theta)
c
c     set wb=0 for m=0
c
  131 ix = indx(0,n,nlat)
      do 220 i=1,imid
      wb(i,ix) = 0.d0
220   continue
c
c     compute wb for m=1,n
c
      dcf = dsqrt(dble(float(n+n+1))/dble(float(4*n*(n+1)*(n+n-1))))
      do 230 m=1,n
      ix = indx(m,n,nlat)     
      abel = dcf*dsqrt(dble(float((n+m)*(n+m-1))))
      bbel = dcf*dsqrt(dble(float((n-m)*(n-m-1))))
      if(m.ge.n-1) go to 231
      do 229 i=1,imid
      wb(i,ix) = abel*dpbar(i,m,nz) + bbel*dpbar(i,m+2,nz)
  229 continue
      go to 230
  231 do 228 i=1,imid
      wb(i,ix) = abel*dpbar(i,m,nz)
  228 continue
  230 continue
  100 continue
      return 
      end
 
