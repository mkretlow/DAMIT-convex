/*
Copyright (C) 2006  Mikko Kaasalainen, Josef Durech

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

#include <stdlib.h>
#include <stdio.h>

double *vector_double(int length)
{
   double *p_x;

   if ((p_x = (double *) malloc((length + 1) * sizeof(double))) == NULL)
   {
      fprintf(stderr, "failure in 'vector_double()' \n");
      fflush(stderr);
   }
   return (p_x);
}
  
int *vector_int(int length)
{
   int *p_x;

   if ((p_x = (int *) malloc((length + 1) * sizeof(long int))) == NULL)
   {
      fprintf(stderr, "failure in 'vector_int()' \n");
      fflush(stderr);
    }
    return (p_x);
}

double **matrix_double(int rows, int columns)
{
   double **p_x;
   int i;
  
   p_x = (double **)malloc((rows + 1) * sizeof(double *));
   for (i = 0; (i <= rows) && (!i || p_x[i-1]); i++) 
      p_x[i] = (double *) malloc((columns + 1) * sizeof(double));
   if (i < rows) 
   {
      fprintf(stderr,"failure in 'matrix_double()' \n");
      fflush(stderr);
   }
   return (p_x);
}

int **matrix_int(int rows, int columns)
{
   int **p_x;
   int i;

   p_x = (int **) malloc((rows + 1) * sizeof(int *));
   for (i = 0; (i <= rows) && (!i || p_x[i-1]); i++) 
      p_x[i] = (int *) malloc((columns + 1) * sizeof(int));
   if (i < rows)
   {
      fprintf(stderr,"failure in 'matrix_int()' \n");
      fflush(stderr);
   }
   return (p_x);
}

double ***matrix_3_double(int n_1, int n_2, int n_3)
{
   int i, j;
   double ***p_x;
 
   p_x = (double ***) malloc((n_1 + 1) * sizeof(double **));  
   for (i = 0; i <= n_1; i++)   
   {
      p_x[i] = (double **) malloc((n_2 + 1) * sizeof(double *));
      for (j = 0; j <= n_2; j++) 
         p_x[i][j] = (double *)malloc((n_3 + 1) * sizeof(double));
      if (j < n_2) 
      {
         fprintf(stderr,"failure in 'matrix_3_double' \n");
         fflush(stderr);
      }
   }
   if (i < n_1) 
   {
      fprintf(stderr,"failure in 'matrix_3_double' \n");
      fflush(stderr);
   }

   return (p_x);
}
  
void deallocate_vector(void *p_x)
{
   free((void *) p_x);
   p_x = NULL;
}

void deallocate_matrix(void **p_x, int rows)
{
   int i;
    
   for (i = 0; i <= rows; i++)
   {
      free(p_x[i]);
      p_x[i] = NULL;
   }
   free(p_x);
}

void deallocate_matrix_3(void ***p_x, int n_1, int n_2)
{
   int i, j;
    
   for (i = 0; i <= n_1; i++)
   {  
      for (j = 1; j <= n_2; j++)
      {
         free(p_x[i][j]);
         p_x[i][j] = NULL;
      }
      free(p_x[i]);
      p_x[i] = NULL;
   }
   free(p_x);
}

