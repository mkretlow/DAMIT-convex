/* cross product of two vectors c = a x b */

/*
Copyright (C) 2006  Josef Durech, Mikko Kaasalainen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

void cross_product(double a[], double b[], double c[])
{
  c[1] = a[2] * b[3] - a[3] * b[2];
  c[2] = -a[1] * b[3] + a[3] * b[1];
  c[3] = a[1] * b[2] - a[2] * b[1];
}

