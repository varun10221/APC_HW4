#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <stdlib.h>

#define NRA 4
#define NCA 4
#define NRB 4
#define NCB 4



double A[NRA][NRA],
       B[NRA][NRA],
       C[NRA][NRA];



void print_Mat (int m , double X[][m]);
int coordinates [2], dims[2], period[2], send_coordinates [2];

void
Init_zero_Mat (int m, double X[][m])
{
    int i,j;
    for (i = 0; i < m ; i++)
      {
        for (j =0; j < m; j++)
         X[i][j] = 0;
      }
}
           
void
Matrix_mul (int m, double X[][m],double Y[][m],double Z[][m])
{
    int i,j, k;
   for (i = 0; i < m; i++)
    {
      for (j =0; j <m; j++)
        {
          Z[i][j] = 0;
          for (k =0; k <m; k++)
            {
              Z[i][j] += X[i][k] * Y[k][j];
            }
         }
     }

}



void
Matrix_init (int m , int n, double A[][n])
{
   int i,j;
   for (i=0; i < m ; i++)
    {
      for (j=0; j < n; j++)
        {
           A[m][n] = (double) (rand () %100);
        }
     }
//    printf ("in verify\n");
//    print_Mat (m, A);
}


void
print_Mat (int m, double X[][m])
{
   int i,j;
   for (i = 0 ; i <m; i++)
    {
       printf ("\n");
      for (j =0; j < m; j++)
        printf("%f\t", X[i][j]);
        
     }

}



double
Verify (int m, double X[][m], double Y[][m], double Z[][m])
{
  double Q[m][m];
  double start, end;
  start = MPI_Wtime ();
 
  int i,j,k;

   for (i=0; i<m; i++)
    {
     for (j = 0; j <m; j++)
       { Q[i][j] = 0;
         for (k=0; k<m; k++)
           Q[i][j] = Q[i][j] + X[i][k]*Y[k][j];
       }

    }

   end = MPI_Wtime ();
   int sum = 0;
   for (j = 0; j <m; j++)
     {
       for (k=0; k<m; k++)
          sum  =(int)(Q[j][k] - Z[j][k]);
     }
  //     print_Mat (m, Z);
  //     print_Mat (m, Q);
       printf ("error:%d", sum);
    return end - start;
}





MPI_Status status1;
MPI_Status status2;
MPI_Status status3;



int
main (int argc, char **argv)
{
   int row_receive, col_receive;
   int world_rank, world_size;
   int  source, destination;
   double start_time, end_time;
   int i,j;
   int row_i, column_i, cycle;
   int rank2;  
 
   MPI_Comm comm2;
      

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
   MPI_Comm_size (MPI_COMM_WORLD, &world_size);
   
   double root_p;                  /* sqrt of no of processors */
   root_p = sqrt ((double) world_size);
    
   if (NRA % (int) root_p != 0)
     {
      printf ("Please enter a processor count which a perfect square and a multiple ");
      MPI_Abort (MPI_COMM_WORLD, 1);
     }
    
   int sub_matrix = NRA / root_p;
   
   /* Need to create a grid of root p x root p processors */
   dims [0] = (int) root_p;
   dims [1] = (int) root_p;

   period [0] = 1;
   period [1] = 1;

   /* Now Matrix is made up of sub-matrices of size n/root p */

   double sub_A [sub_matrix][sub_matrix];   
   double sub_B [sub_matrix][sub_matrix];   
   double sub_C [sub_matrix][sub_matrix];   
   double sub_CT [sub_matrix][sub_matrix];   
   
   /* Now creating a cartesian topology */
   
   MPI_Cart_create (MPI_COMM_WORLD, 2, dims, period, 0, &comm2);     
   
   /* NOw getting a new rank */
   MPI_Comm_rank (comm2, &world_rank);
   
   /*Determine process co ordinate based on rank */
   MPI_Cart_coords (comm2, world_rank, 2, coordinates);

   

   Init_zero_Mat (sub_matrix, sub_C);

   if (world_rank == 0)
     {
        Matrix_init (NRA, NRA, A);
        Matrix_init (NRA, NRA, B);
  
        //print_Mat (NRA, A);
        //print_Mat (NRA, B);
     
        Init_zero_Mat (sub_matrix, sub_C); 

        /* Let us send each portion of A and B and start multiplying */         
         start_time = MPI_Wtime ();

         for (i = 0; i < root_p; i++)
           {
              for (j = 0; j < root_p; j++)
                 {
                    if ( i != 0 || j != 0)
                      {
                         send_coordinates [0] = i;
                         send_coordinates [1] = j;
                         row_i = -1;
                         int k;
                         for (k = i * sub_matrix; k < i * sub_matrix + sub_matrix; k++)
                           {
                             column_i = 0;
                             row_i ++;
                             int l;
                             for (l = j *sub_matrix; l < j * sub_matrix + sub_matrix; l++)
                              {
                                       
                                 sub_A[row_i][column_i] = A[k][l];
                                 sub_B[row_i][column_i] = B[k][l];
                                 column_i++;
                              } 
                            }
                            
                            /* Make the co ordinate reference to column and send it to processor pij */
                            send_coordinates [0] = i;
                            send_coordinates [1] = ((j - i) < 0) ? (j-i) + root_p : (j-i);
                            MPI_Cart_rank (comm2, send_coordinates, &rank2);
                            MPI_Send (sub_A, sub_matrix * sub_matrix, MPI_DOUBLE, rank2, 1, comm2);
                            send_coordinates [0] = ((i-j) < 0) ? (i-j) + root_p : i-j;
                            send_coordinates [1] = j;
                        
                            
                            MPI_Cart_rank (comm2, send_coordinates, &rank2);
                            MPI_Send (sub_B, sub_matrix * sub_matrix, MPI_DOUBLE, rank2, 2, comm2);
 

                      }
                   }
              } 

             
          /* NOws send  to process 0 */
          for (i =0 ; i<sub_matrix; i++)
           {
            for ( j = 0; j < sub_matrix; j++)
              {
                 sub_A[i][j] = A[i][j];
                 sub_B[i][j] = B[i][j];
               }
            }
     
           /* calculate c for matrix in process 0 */
           /* Todo: use in function */
           for (cycle = 0; cycle < sub_matrix; cycle++)
            {
                  
                Matrix_mul (sub_matrix, sub_A, sub_B, sub_C);      
            
            
               MPI_Cart_shift (comm2, 1, -1, &source, &destination);
               MPI_Sendrecv_replace (sub_A, sub_matrix * sub_matrix, MPI_DOUBLE,destination, 1, source, 1, comm2, &status1);
               
               MPI_Cart_shift (comm2, 0, -1, &source, &destination);
               MPI_Sendrecv_replace (sub_B, sub_matrix * sub_matrix, MPI_DOUBLE,destination, 2, source, 2, comm2, &status2);
            
     
             }
             
             }
/*end of master */


    else 
      {
         MPI_Recv (sub_A, sub_matrix * sub_matrix, MPI_DOUBLE, 0, 1, comm2, &status1);
         MPI_Recv (sub_B, sub_matrix * sub_matrix, MPI_DOUBLE, 0, 2, comm2, &status2);

         
       for (cycle = 0; cycle < sub_matrix; cycle ++)
         {
            Matrix_mul (sub_matrix, sub_A, sub_B, sub_C);
         
            MPI_Cart_shift (comm2, 1, -1, &source, &destination);
            MPI_Sendrecv_replace (sub_A, sub_matrix * sub_matrix, MPI_DOUBLE,destination, 1, source, 1, comm2, &status1);

            MPI_Cart_shift (comm2, 0, -1, &source, &destination);
            MPI_Sendrecv_replace (sub_B, sub_matrix * sub_matrix, MPI_DOUBLE,destination, 2, source, 2, comm2, &status2);
         }

        /* send final result to process 0 */
  
        MPI_Send (sub_C, sub_matrix * sub_matrix, MPI_DOUBLE, 0 , world_rank, comm2);
    
      }

    if (world_rank == 0)
      {
         Init_zero_Mat (NRA , C);
      
        int k;
     for (i =1; i < world_size; i++)
      {
         MPI_Recv (sub_CT, sub_matrix * sub_matrix, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,comm2, &status3);
         MPI_Cart_coords (comm2, status3.MPI_TAG, 2, send_coordinates);
          
         row_receive = send_coordinates [0];
         col_receive = send_coordinates [1];
         row_i = -1;
         column_i = 0;
         
         for ( j = row_receive * sub_matrix; j < row_receive * sub_matrix + sub_matrix; j++)
          {
            row_i ++;
            for ( k = col_receive * sub_matrix; k < col_receive * sub_matrix + sub_matrix; k++)
               {
                  C[j][k] = sub_CT[row_i][column_i];
                  column_i ++;
               }
             column_i = 0;
           }
        }
 
     /* On process 0 */
      for (i = 0; i < sub_matrix; i++)
        {
          for (j =0 ; j <sub_matrix; j++)
            {
                C[i][j] = sub_C[i][j];
            }
         }
      end_time = MPI_Wtime ();
      double serial = Verify (NRA, A,B,C);
      printf ("speedup :%f s",serial/( end_time - start_time));
      }

    MPI_Finalize ();
    return 0;
  }

   

           
