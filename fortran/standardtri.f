c       Removes facets with less than three vertices, and splits
c       those with more than three into triangles
c
c       CR Mikko Kaasalainen 2005
c
c       INPUT FILE: General polyhedron vertex list (from minkowski.f)
c       OUTPUT FILE: The same polyhedron with triangle facets
c
c	modified by Josef Durech: reads the standard input and puts results to the standard 
c	output


c	Copyright (C) 2005  Mikko Kaasalainen
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



        program standardtri
        implicit double precision (a-h,o-z)
        double precision xver(10000,3)
        integer iver(10000,3),itmp(100)
        character*15 filnam
        read(5,*) nver,nfac
        do i=1,nver
         read(5,*) (xver(i,j),j=1,3)
        end do
        nf=0
        do i=1,nfac
         read(5,*) npts
         if (npts.lt.3) then
          read(5,*) idum
          goto 10
         end if
         read(5,*) (itmp(j),j=1,npts)
         do j=2,npts-1
          nf=nf+1
          iver(nf,1)=itmp(1)
          iver(nf,2)=itmp(j)
          iver(nf,3)=itmp(j+1)
         end do
 10     end do
        write(6,*) nver,nf
        do i=1,nver
         write(6,*) (xver(i,j),j=1,3)
        end do
        do i=1,nf
c         write(6,*) 3
         write(6,*) (iver(i,j),j=1,3)
        end do
        end
