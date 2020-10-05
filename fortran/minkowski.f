c       A program for determining a convex polyhedron from 
c       its facet areas+normals a la Minkowski
c
c       CR Mikko Kaasalainen 2005
c
c       See Kaasalainen et al. 2001 for details. The procedure is provably
c       convergent, so if something goes wrong it is due to some adjustments
c       of numerical parameters, or somehow the given set of facets is bad
c       (e.g. not a convex set, or some areas are numerically too small).
c       NOTE: the compiler must be able to do full double precision --
c       e.g., GNU g77 is not always accurate enough.
c
c       INPUT FILE:
c       FACETS: first line gives the number of facets, then follow
c       facet area
c       outward unit normal x,y,z coordinates
c
c       OUTPUT FILE:
c       SHAPE (VERTICES): the first line gives the numbers of vertices 
c       and facets, then follow the vertex x,y,z coordinates, then for each 
c       facet the number of vertices (arbitrary, use standardtri.f to get 
c       them rendered as triangles if neccessary) and, on a separate line, 
c       the order numbers of the facet vertices (anticlockwise seen from 
c       outside the body).
c
c       INFO PARAMETERS (SCREEN):
c       tk: step legth, not very important
c       volume: not very important (this is the maximized quantity)
c       alphak: this indicates the many-dimensional "angle" away from a
c       perfect solution. Values below 0.01 are usually very close to
c       final solution, this version stops iterating at alphak=epsilon=0.001.
c
c       NOTE: user can stop iteration at any time (e.g., when alphak=0.005)
c       by simple control-C interrupt command, the latest version is 
c       automatically stored in the output file (here with an interval of 
c       nwrite=20 iterations).
c       
c	modified by Josef Durech: reads the standard input and puts results to the standard
c	output 


c 	Copyright (C) 2005  Mikko Kaasalainen
c
c	This program is free software; you can redistribute it and/or
c	modify it under the terms of the GNU General Public License
c	as published by the Free Software Foundation; either version 2
c	of the License, or (at your option) any later version.
c
c	This program is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c	GNU General Public License for more details.
c
c	You should have received a copy of the GNU General Public License
c	along with this program; if not, write to the Free Software
c	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


        program Minkowski
        implicit none
        integer nmax,nedmax,interm,nwrite
        parameter(nmax=6000,nedmax=200)
        double precision initarea(nmax),anorm(3),dhelp,dot,
     1  epsilon,initm,bnorm(3),m,prevolume,eta,initareasum,
     2  initarealen,gdot,d(nmax),x(0:nmax),y(0:nmax),z(0:nmax),
     3  p1(3),p2(3),p0(3),cross(3),lower,rx(3,nmax),sqlength,
     4  maxpxk,h(1200),normal(3),max,rv(3),pxk(nmax),sv(3),tiny,
     5  volume,addvect(3),ra(3),rb(3),area(nmax),cmass(3),scale,
     6  avect(3),bvect(3),av(3),bv(3),cv(3),xk(nmax),scoef,
     7  retemp,premaxpxk,prearea(nmax),arealen,alphak,tk,coeff,
     8  f,ek(nmax),newx(nmax),w(nmax),wmin,norm(3),nor(nmax,3),
     9  xf(nmax),yf(nmax),zf(nmax)
        integer i,j,k,ncoef,posroute,numfaces,numedges(nmax),
     1  edge(nmax,nedmax),vert,n,modlo,l,temp,r(nmax,nedmax),
     2  numvert(nmax),numvertices,maxind,first,iter,iread,numfini
        character*15 filnam,filini
        common/facevalues/numfaces/points/x,y,z
        tiny=1e-8       
        iter=0
        nwrite=20
c       The vector initarea (areas of faces) and the normals are read.
c        write(6,*) 'Give the name of the file to be read'
c        read(5,*) filnam
        read(5,*) numfaces
        do i=1,numfaces
         read(5,*) initarea(i)
         read(5,*) anorm(1),anorm(2),anorm(3)
         dhelp=sqrt(dot(anorm,anorm))
         do j=1,3
          nor(i,j)=anorm(j)/dhelp
         end do
        end do
c       The checkvalue for stopping the iteration.
        epsilon=.005
c       An auxiliary value initm is computed.
        initm=0.
        do i=1,numfaces
         do k=1,3
          anorm(k)=nor(i,k)
         end do
         do j=1,numfaces
          do k=1,3
           bnorm(k)=nor(j,k)
          end do
          dhelp=1-(dot(anorm,bnorm))**2
          if (dhelp.gt.0.0001) then
           m=1./sqrt(dhelp)
           if (m.gt.initm) initm=m
          end if
         end do
        end do
        interm=0
        ncoef=0
        posroute=1
        prevolume=0.
        eta=1.3
c       The vector of the distances of the faces from the origin
c       is initialized.
        initareasum=0.
        do i=1,numfaces
         initareasum=initareasum+initarea(i)
        end do
c       Rescale such that initareasum=1000.
        scale=1000./initareasum
        do i=1,numfaces
         initarea(i)=scale*initarea(i)
        end do
        initareasum=100.
        initarealen=sqrt(gdot(initarea,initarea))
        do i=1,numfaces
         d(i)=initarealen**2/initareasum
        end do
c       The origin is x(0),y(0),z(0)
        x(0)=0.
        y(0)=0.
        z(0)=0.
c       The polytope in R^3 is transformed into a polytope in dual space,
c       i.e., the faces are transformed into points x(i),y(i),z(i).
10      do i=1,numfaces
         x(i)=nor(i,1)/d(i)
         y(i)=nor(i,2)/d(i)
         z(i)=nor(i,3)/d(i)
        end do
        call convhull(numfaces,numedges,edge,nmax,nedmax)
c        write(6,*) 'CH done'
c       (convhull's arguments include x,y,z as input as well;
c       they are given in the common block points)
c       The convex hull of the dual polytope is transformed into a 
c       polytope in R^3, i.e., the faces surrounding a vertex in the dual
c       space become vertices surrounding a face in R^3.
        k=1
c       k counts the total number of vertices, vert the vertices for face i
        do i=1,numfaces
         vert=0
         do n=1,numedges(i)
          call vector(i,edge(i,n),p1)
          call vector(i,edge(i,modlo(n+1,numedges(i))),p2)
          call vector(0,i,p0)
          call crossproduct(p1,p2,cross)
          lower=dot(p0,cross)
c       Coordinates of a vertex
          do l=1,3
           rx(l,k)=cross(l)/lower
          end do
c       Check if this vertex has been computed earlier
          do j=1,k-1
           sqlength=0.
           do l=1,3
            sqlength=sqlength+(rx(l,j)-rx(l,k))**2
           end do
c       The test allows for some rounding errors etc.
           if (sqlength.lt.tiny) then
            temp=j
            if (n.gt.1) then
c       Ignore this vertex if met earlier at this face
             if ((temp.eq.r(i,vert)).or.(temp.eq.r(i,1))) goto 70
            end if
c       If met earlier but not at this face, increase vert but not k
            goto 60
           end if
c       End for j
          end do
          temp=k
          k=k+1
          if (k.gt.(nmax-10)) write(6,*) k
60        vert=vert+1
c       r holds the vertices of facet i
          r(i,vert)=temp
c       End for n
70       end do
         numvert(i)=vert
c       End for i
        end do
        numvertices=k-1
c       Now the points of the dual space can be erased, and x,y,z will
c       represent coordinates in R^3.
        do i=1,numvertices
         x(i)=rx(1,i)
         y(i)=rx(2,i)
         z(i)=rx(3,i)
        end do
        maxpxk=0.
        do i=1,numfaces
         h(i)=0.
c       Compute 'vertices' and auxiliary values h for faces
c       that have no vertices
         if (numvert(i).eq.0) then
          do j=1,3
           normal(j)=nor(i,j)
          end do
          max=0.
          do l=1,numvertices
           call vector(0,l,rv)
           if (dot(rv,normal).gt.max) then
            max=dot(rv,normal)
            maxind=l
           end if
          end do
          numvert(i)=1
          r(i,1)=maxind
          do l=1,numvertices
           call vector(0,l,rv)
           if ((dot(rv,normal).eq.max).and.(l.ne.maxind)) then
            numvert(i)=2
            r(i,2)=l
           end if
          end do
          h(i)=d(i)-max
         end if
c       The circuits of the edges of each face are computed
c       and the largest circuit is seeked.
         pxk(i)=0
         do j=1,numvert(i)
          call vector(r(i,j),r(i,modlo((j+1),numvert(i))),sv)
          pxk(i)=pxk(i)+sqrt(dot(sv,sv))
         end do
         if (maxpxk.lt.pxk(i)) maxpxk=pxk(i)
        end do
c       Some auxiliary values are computed
        volume=0.
        do i=1,numfaces
         do j=1,3
          addvect(j)=0.
         end do
         do n=1,numvert(i)
          call vector(0,r(i,n),ra)
          call vector(0,r(i,modlo(n+1,numvert(i))),rb)
          call crossproduct(ra,rb,cross)
          do j=1,3
           addvect(j)=addvect(j)+cross(j)
          end do
         end do
         area(i)=sqrt(dot(addvect,addvect))/2.
         volume=volume+d(i)*area(i)/3.
        end do
        iter=iter+1
        if (volume.ge.prevolume) then
         do j=1,3
          cmass(j)=0.
         end do
         do i=1,numfaces
          do n=1,numvert(i)-2
           call vector(r(i,1),r(i,n+1),avect)
           call vector(r(i,1),r(i,n+2),bvect)
           call crossproduct(avect,bvect,cross)
           call vector(0,r(i,1),av)
           call vector(0,r(i,n+1),bv)
           call vector(0,r(i,n+2),cv)
           do j=1,3
            cmass(j)=cmass(j)+d(i)*
     1      sqrt(dot(cross,cross))*(av(j)+bv(j)+cv(j))
           end do
          end do
         end do
         do j=1,3
          cmass(j)=cmass(j)/(24*volume)
         end do
         do i=1,numfaces
          do j=1,3
           norm(j)=nor(i,j)
          end do
          xk(i)=d(i)-dot(norm,cmass)-h(i)
         end do
         scoef=initarealen**2/gdot(xk,initarea)
         do i=1,numfaces
          xk(i)=scoef*xk(i)
          area(i)=scoef**2*area(i)
         end do
         maxpxk=scoef*maxpxk
         volume=scoef**3*volume
         if (prevolume.ne.0) then
          retemp=eta**ncoef*(arealen*sin(alphak))**2/(6*initm*
     1    premaxpxk)
         else
          retemp=0.
         end if
         if (((volume-prevolume).ge.retemp).and.(posroute.eq.1)) then
          ncoef=ncoef+1
         else
          ncoef=ncoef-1
         end if
c       else, i.e., if volume.lt.prevolume
        else
         ncoef=ncoef-1
c       No changes are made to xk:s.
         volume=prevolume
         maxpxk=premaxpxk
         do i=1,numfaces
          area(i)=prearea(i)
         end do
        end if
c        if((((volume-prevolume)/volume).lt.1.0d-7)
c     -  .and.(((volume-prevolume)/volume).
c     -     gt.1.0d-15)) goto 100
 80     prevolume=volume
        premaxpxk=maxpxk
        do i=1,numfaces
         prearea(i)=area(i)
        end do
c       The iteration step and the angle alphak are computed.
        arealen=sqrt(gdot(area,area))
        alphak=acos(gdot(area,initarea)/(arealen*initarealen))
c       Test whether to stop iterating
        if (alphak.lt.epsilon) goto 100
        tk=eta**ncoef*arealen*sin(alphak)/(3*initm*maxpxk)
        coeff=gdot(area,initarea)/(initarealen**2)
        do i=1,numfaces
         f=area(i)-coeff*initarea(i)
         ek(i)=f/(arealen*sin(alphak))
        end do
        do i=1,numfaces
         newx(i)=xk(i)+tk*ek(i)
         if (newx(i).le.0) goto 90
        end do
        do i=1,numfaces
         d(i)=newx(i)
        end do
        posroute=1
        goto 10
90      posroute=-1
        first=0
        do i=1,numfaces
         w(i)=xk(i)/(coeff*initarea(i)-area(i))
         if ((w(i).gt.0).and.(first.eq.0)) then
          wmin=w(i)
          first=1
         end if
         if ((w(i).gt.0).and.(w(i).lt.wmin)) wmin=w(i)
        end do
        tk=0.9*arealen*wmin*sin(alphak)
        do i=1,numfaces
         d(i)=xk(i)+tk*ek(i)
        end do
        goto 10
100     write(6,*) numvertices,numfaces
        scoef=scoef*sqrt(initarealen/arealen)
        scoef=scoef/sqrt(scale)
        do i=1,numvertices
         xf(i)=x(i)-cmass(1)
         yf(i)=y(i)-cmass(2)
         zf(i)=z(i)-cmass(3)
         write(6,110) scoef*xf(i),scoef*yf(i),scoef*zf(i)
        end do
        do i=1,numfaces
          write(6,*) numvert(i)
          write(6,*) (r(i,j),j=1,numvert(i))
        end do
 110    format(3(E23.16,1X))
        end

c       The dot product in gradient space
        double precision function gdot(avect,bvect)
        implicit none
        integer i,numfaces,nmax
        parameter(nmax=6000)
        double precision avect(nmax),bvect(nmax)
        common/facevalues/numfaces
        gdot=0.
        do i=1,numfaces
         gdot=gdot+avect(i)*bvect(i)
        end do
        end

        subroutine convhull(numpoints,numedges,edge,nedge1,nedge2)
        implicit double precision (a-h,o-z)
        parameter(nmax=6000,tiny=1e-16)
        double precision x(0:nmax),y(0:nmax),z(0:nmax),
     1  stick(3),tvec(3),base(3),vectn(3),crvec(3),cvec(3) 
        integer numedges(nedge1),edge(nedge1,nedge2),
     1  listch(nedge1),inlist(nedge1)
        common/points/x,y,z
c       origin is the zeroth point (if not inside, use cmass)
        x(0)=0.
        y(0)=0.
        z(0)=0.
c       First find the largest z-value (any min,max,coord would do)
c       to establish one vertex of the CH
        zmax=z(1)
        mind=1
        do i=2,numpoints
         if (z(i).gt.zmax) then
          zmax=z(i)
          mind=i
         end if
        end do
c       listch is the list of all CH vertices; inlist is either
c       1 or 0 depending on whether the point is in listch
        do i=1,numpoints
         inlist(i)=0
         numedges(i)=0
        end do 
        listch(1)=mind
        inlist(mind)=1
c       Then find the point with the largest 'umbrella angle'
c       from the first point (umbrella 'stick' towards origin)
c       to establish one edge of the CH
        call vector(mind,0,stick)
        sticklen=sqrt(dot(stick,stick))
        dmin=1.
        do i=1,numpoints
         if (i.ne.mind) then
          call vector(mind,i,tvec)
          tveclen=sqrt(dot(tvec,tvec))
          dtest=dot(stick,tvec)/(sticklen*tveclen)
          if (dtest.lt.dmin) then
           dmin=dtest
           mind2=i
          end if
         end if
        end do
        listch(2)=mind2
        inlist(mind2)=1
        edge(mind,1)=mind2
        edge(mind2,1)=mind
        numedges(mind)=1
        numedges(mind2)=1
c       nvert counts the vertices for which all edges have been 
c       determined; nconn counts the points that have been
c       connected to any other point (and are thus vertices of CH)
        nvert=0
        nconn=2
c       Begin nvert-loop
10      inda=listch(nvert+1)
        indb=edge(inda,1)
        call vector(inda,indb,base)
c       Begin numedges-loop: find a new edge for the vertex inda.
c       Go through all points; i is the trial point, indn the
c       reference point
20      do i=0,numpoints
         if ((i.ne.inda).and.(i.ne.indb)) then
          call vector(inda,i,cvec)
          clen=sqrt(dot(cvec,cvec))
          if ((i.eq.0).or.
     1        (dot(cvec,crvec).lt.(clen*crlen*(-tiny)))) then
c       the trial point becomes the new reference point (.lt.(-tiny) in
c       the above test for counterclockwise ordering of edges,
c       .gt.tiny for clockwise ordering)
           indn=i
           do j=1,3
            vectn(j)=cvec(j)
           end do
           call crossproduct(vectn,base,crvec)
           crlen=sqrt(dot(crvec,crvec))
          end if
         end if
        end do
c       check if we've not yet gone round inda (not found all its edges)
        if ((numedges(inda).eq.1).or.
     1  ((indn.ne.edge(inda,1)).and.(indn.ne.edge(inda,2)))) then
c       check if a new vertex of CH has been found
         if (inlist(indn).eq.0) then
          nconn=nconn+1
          listch(nconn)=indn
          inlist(indn)=1
          numedges(indn)=1
          edge(indn,1)=inda
         end if
         numedges(inda)=numedges(inda)+1
         edge(inda,numedges(inda))=indn
         indb=indn
         do j=1,3
          base(j)=vectn(j)
         end do
         goto 20
        end if
        nvert=nvert+1
c       check if we haven't yet found edges for all connected points
c       (otherwise we have determined the CH)
        if (nvert.lt.nconn) goto 10
        return
        end

c       The cross product of two vectors
        subroutine crossproduct(avect,bvect,crossvect)
        implicit none
        double precision avect(3),bvect(3),crossvect(3)
        crossvect(1)=avect(2)*bvect(3)-avect(3)*bvect(2)
        crossvect(2)=avect(3)*bvect(1)-avect(1)*bvect(3)
        crossvect(3)=avect(1)*bvect(2)-avect(2)*bvect(1)
        end

c       The dot product of two vectors
        double precision function dot(avect,bvect)
        implicit none
        double precision avect(3),bvect(3)
        integer n
        dot=0.
        do n=1,3
         dot=dot+avect(n)*bvect(n)
        end do
        end

c       The modulo function
        integer function modlo(index,modsize)
        implicit none
        integer index,modsize
        if (index.gt.modsize) then
         modlo=index-modsize
        else
         if (index.lt.1) then
          modlo=index+modsize
         else
          modlo=index
         end if
        end if
        end

c       A procedure for constructing a vector from two points
        subroutine vector(start,fin,vect)
        implicit none
        integer fin,start,nmax
        parameter(nmax=6000)
        double precision vect(3),x(0:nmax),y(0:nmax),z(0:nmax)
        common/points/x,y,z
        vect(1)=x(fin)-x(start)
        vect(2)=y(fin)-y(start)
        vect(3)=z(fin)-z(start)
        end



