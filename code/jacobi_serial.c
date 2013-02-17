// Serial implementation of jacobi method
// compile: gcc jacobi_serial.c time.c -lpthread -o jacobi_serial
// execute: ./jacobi_serial size_of_fine num_iterations
// example: ./multi_parallel 200 500

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAXSIZE 200  /* maximum coarse matrix size */
#define MAXITER 1000   // Maximum number of iterations
#define EPSILON 0.005

//*********************global variables and function definitions***************


int fine_dim;
int coarse_dim_with_boundary;
int fine_dim_with_boundary;


float *grid_fine;
float *newMatrix_fine;


void printMatrix(float *matrix, int fine_dim);
float max(float a, float b);
float absolute(float a);
void jacobi(float *old, float *new, int sent_dim, int num_iter); 



int main (int argc, const char * argv[]) 

{
    
    int i, j, kk;
    
    int x, y;
    
    
    float maxdiff;
    float Finalmaxdiff = 0.0;
    int solution_iter;
    
    float time;
    FILE *fp;
    
    
    
    //get command line arguments
    fine_dim = (argc > 1)? atoi(argv[1]) : MAXSIZE;
    
    
    solution_iter = (argc > 2)? atoi(argv[2]) : MAXITER;
    
    if (fine_dim > MAXSIZE) fine_dim = MAXSIZE;
    
    if ((solution_iter > MAXITER)||(solution_iter <= 0)) solution_iter = MAXITER;
    
    
    //accomodate the boundary conditions in size
    fine_dim_with_boundary=fine_dim+2;
    
    //calculate the coarse grid with double spacing
    
    
    
    printf("\n\n******* Fine Grid Size: %d and Number of iterations: %d *******\n\n", fine_dim, solution_iter);
    
    
    
    //create the matrices
    grid_fine = (float* )malloc(fine_dim_with_boundary * fine_dim_with_boundary * sizeof(float));
    newMatrix_fine = (float* )malloc(fine_dim_with_boundary * fine_dim_with_boundary * sizeof(float));
    
    //Set inner values
    for (i=1; i<=fine_dim; i++) 
    {
        
        for (j=1; j<=fine_dim; j++) 
        {
            
            grid_fine[(i*fine_dim+j)-1]=0;
            
        }
        
    }
    
    
    //Set boundary conditions
    for (i=0; i<fine_dim_with_boundary; i++)
    {
        
        grid_fine[i]=1;// First row
        
        grid_fine[i*fine_dim_with_boundary]=1; // First column
        
        grid_fine[(i*fine_dim_with_boundary)+(fine_dim+1)]=1; // Last column
        
        grid_fine[(fine_dim_with_boundary*(fine_dim_with_boundary-1))+i]=1; // Last Row
        
        
        newMatrix_fine[i]=1;// First row
        
        newMatrix_fine[i*fine_dim_with_boundary]=1; // First column
        
        newMatrix_fine[(i*fine_dim_with_boundary)+(fine_dim+1)]=1; // Last column
        
        newMatrix_fine[(fine_dim_with_boundary*(fine_dim_with_boundary-1))+i]=1; // Last Row         
    }
    
  
    
    
    
    
    time = timer();
    
        
    jacobi(grid_fine, newMatrix_fine, fine_dim_with_boundary, solution_iter);
        
    
    Finalmaxdiff = 0.0;
    
    
    for (i=1; i<fine_dim_with_boundary; i++) 
        
    {
        
        for (j=1; j<fine_dim_with_boundary; j++) 
            
        {
            Finalmaxdiff = max(Finalmaxdiff, absolute(1 - grid_fine[(i*fine_dim_with_boundary)+j]));
                      
        }
    }
    
    printf("\nFinal maxdiff: %f\n\n",Finalmaxdiff);
    
    time=timer()-time;
    printf("Elapsed time: %f\n",time/1000000.0);
    fp=fopen("jacobi_serial_data.txt", "wb");
    
    if(fp==NULL) 
    {
        
        printf("Error: can't open file.\n");
        
        exit(0);
        
    }
    
    
    //save in file
    for (i=0; i<fine_dim_with_boundary; i++)
        
    {
        for(j=0; j<fine_dim_with_boundary; j++)
            
        {
            
            fprintf(fp, "%f ", grid_fine[i*fine_dim_with_boundary+j]); 
            
        }
        
        fputs("\n", fp);
        
    }
    
    printf("data saved in jacobi_serial_data.txt\n");
    fclose(fp);
    
    
    
    free(newMatrix_fine);
    free(grid_fine);
    
    
    return 0;
    
}//end Main




float max(float a, float b) 
{
    return a > b ? a : b;
}

float absolute(float a)
{
    return a >= 0 ? a : a*(-1);
}

void printMatrix(float *matrix, int s)
{
    int i,j;
    printf("\n\n");
    for (i=0; i<s; i++)
    {
        for(j=0; j<s; j++)
        {
            
            printf("%f ",matrix[(i*s)+j]);
            
        }
        printf("\n");
    }
}


void jacobi(float *old, float *new, int sent_dim, int num_iter) 
{
    
    int j;
    int i, k;
    
    for (k=0; k<(num_iter/2); k++) 
    {
        
        //compute the new inner points
        for (i=1; i<sent_dim-1; i++) 
            
        {
            for (j=1; j<sent_dim-1; j++) 
                
                
            {
                
                new[(i*sent_dim)+j] = (old[((i+1)*sent_dim)+j] + old[((i-1)*sent_dim)+j] + old[(i*sent_dim)+j-1] + old[(i*sent_dim)+j+1])*0.25;
                
            }
            
            
        }
        
        //loop unrolling
        for (i=1; i<sent_dim-1; i++) 
            
        {
            for (j=1; j<sent_dim-1; j++) 
                
            {
                
                old[(i*sent_dim)+j] = (new[((i+1)*sent_dim)+j] + new[((i-1)*sent_dim)+j] + new[(i*sent_dim)+j-1] + new[(i*sent_dim)+j+1])*0.25;
                
            }
            
        }
        
    }
    
}


