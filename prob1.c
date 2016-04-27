#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>



#define FROM_MASTER 1
#define FROM_WORKER 2
#define N_ROWS 256
#define N_COL  256
#define root 0

MPI_Status status;


/* Define globally the matrices */
                                     /*  BT will be B transpose
                                       C will be result */


int **
Matrix_Alloc (int m, int n)
{
     int **X;
     X = (int **) calloc (m , sizeof (int *));
     X[0] = (int *) calloc (sizeof (int) , m * n);
      int i;
      for (i = 0 ; i < m; i++)
        X[i] =  (*X + n * i);

     return X;
}


void
Matrix_init (int **X, int m , int n)
{
  int i, j;
  for (i= 0; i < m; i++)
    {
      for (j=0; j < n; j ++)
         {
            X[i][j] = rand () % 100;
         }
    }

}

int **
Transpose (int **X, int m, int n) 
{

   int **Y;
   Y = Matrix_Alloc (m,n);
   int i,j;
   for (i = 0; i < m; i++)
     for (j = 0; j < n; j++)
       {
          Y[i][j] = X[j][i];   
       }

    return Y;
}

void
print_Mat (int **X, int rows)
{
   int i,j;
   for (i =0 ; i < rows; i++)
     {
       for (j=0; j < rows; j++)
         {
           printf ("%d \t", X[i][j]);
         }
        printf ("\n");
     }
}

double
Verify (int **X, int **Y, int **Z, int m)
{
  int **Q;
  Q = Matrix_Alloc (m,m); 
 
  int i,j,k;

  double start, end;
  start = MPI_Wtime ();
   for (i=0; i<m; i++)
    {
     for (j = 0; j <m; j++)
       {
         for (k=0; k<m; k++)
           Q[i][j] = Q[i][j] + X[i][k]*Y[k][j];
       }

    }

   end = MPI_Wtime ();
//   print_Mat (Z, m);
  // print_Mat (Q,m);
   int sum = 0;
   for (j = 0; j <m; j++)
     {
       for (k=0; k<m; k++)
          sum  = Q[j][k] - Z[j][k];
     }

       printf ("error :%d\n", sum);
      return (end - start);
}





int
main (int argc, char **argv)
{

   int rc;
   int world_size;
   int m,n,i;
   int tag;
   int rows;                        /* An index for rows per node */

   int **A, **B, **BT, **C, **CT;      /* A,B will be input matrices,*/
   rc = MPI_Init (&argc, &argv);
   if (rc != MPI_SUCCESS)
     MPI_Abort (MPI_COMM_WORLD, rc);

   MPI_Comm_size (MPI_COMM_WORLD, &world_size);

   int world_rank;
   MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);

   m = N_ROWS;
   n = N_COL;
  
   A = Matrix_Alloc (m,n);
   B = Matrix_Alloc (m,n);
   Matrix_init (A,m,n);
   Matrix_init (B,m,n);
   BT = Transpose (B,m,n);   
   CT = Matrix_Alloc (m,n);
   int workers = world_size - 1;


   double start_time;
   start_time = MPI_Wtime ();

   int column_per_node = n/workers; /* No of columns each node gets */
   int excess = n % workers;        /* Excess if N_COL is not perfectly divisible by nodes */
   
   int offset = 0;
   int source;
   
   /* Broadcasting array A to all nodes */
   MPI_Bcast (&(A[0][0]), m * n , MPI_INT, root , MPI_COMM_WORLD);

   /* For root node */
   if (world_rank == 0)
   {
   int dest;
   tag = FROM_MASTER;
  
   /* Sending data for workers */
   for (dest = 1; dest <= workers;  dest ++)
     {
        /* "Sending rows of BT", because the rows of BT is actually the column of B */
        rows = (dest <= excess) ? column_per_node + 1 : column_per_node;
        MPI_Send (&offset, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
        MPI_Send (&rows, 1 , MPI_INT, dest, tag, MPI_COMM_WORLD);
        MPI_Send (&BT[offset][0], rows * n, MPI_INT, dest, tag, MPI_COMM_WORLD);
        
        offset = offset + rows;
     }

   /* wait for result from workers */
   tag = FROM_WORKER;
   for (i = 1; i <= workers; i++)
     {
        source = i;
        MPI_Recv (&offset, 1 , MPI_INT, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv (&rows, 1 , MPI_INT, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv (&CT[offset][0], rows * N_COL, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
     } 

   double end_time = MPI_Wtime ();
   
   C = Transpose (CT,m,n);
   double serial;
   serial = Verify (A, B, C, m);

   printf ("spedup: %f " , serial/(end_time - start_time));
     
   }

   /* Code to be execute by worker nodes */
   if (world_rank >= 1)
     {
       tag = FROM_MASTER;
       source = 0;                  /* Root (Master) is set as Source */                
       MPI_Recv (&offset, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
       MPI_Recv (&rows, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
       MPI_Recv (&BT[offset][0], N_COL * rows, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
       
       int j, k;
       /* TODO : Check logic and variable scope */
       for (k = 0; k < N_ROWS; k++)
         {
           for (i = offset; i < offset + rows; i++)
             {
               for (j = 0; j < N_COL; j++)
                 CT[i][k] = CT[i][k] + A[k][j] * BT [i][j];
             }
          }
        
        MPI_Send (&offset, 1, MPI_INT, 0, FROM_WORKER, MPI_COMM_WORLD);
        MPI_Send (&rows, 1, MPI_INT, 0, FROM_WORKER, MPI_COMM_WORLD);
        MPI_Send (&CT[offset][0], rows * N_COL, MPI_INT, 0, FROM_WORKER, MPI_COMM_WORLD);
       }

   MPI_Finalize ();


  return 0;
}
       


    
