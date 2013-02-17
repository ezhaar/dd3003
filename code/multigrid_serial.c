// Serial implementation of multigrid method (2-grid)
// compile: gcc multigrid_serial.c time.c -lpthread -o multi_serial
// execute: ./multi_serial size_of_fine num_iterations
// example: ./multi_parallel 99 500

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAXSIZE 99  /* maximum coarse matrix size */
#define MAXITER 1000   // Maximum number of iterations
#define EPSILON 0.005

//*********************global variables and function definitions***************


int waiting=0, state=0;       /* condition variable for leaving */
int numWorkers;
int numArrived = 0;       /* number who have arrived */
int globalCounter = 0;

float epsilon = EPSILON;

int v_cycles = 8;
int smoothing_iter = 4;
int solution_iter;

int fine_dim;
int coarse_dim;
int coarse_dim_with_boundary;
int fine_dim_with_boundary;


float *grid_fine;
float *newMatrix_fine;
float *grid_coarse;
float *newMatrix_coarse;


void printMatrix(float *matrix, int fine_dim);
float max(float a, float b);
float absolute(float a);
void jacobi(float *old, float *new, int sent_dim, int num_iter); 
void restriction (void);
void interpolate(void);


int main (int argc, const char * argv[]) 

{
    
    int i, j, kk;
    
    int x, y;
    
    
    float maxdiff;
    float Finalmaxdiff = 0.0;
    
    float time;
    FILE *fp;
    
    
    
    //get command line arguments
    coarse_dim = (argc > 1)? atoi(argv[1]) : MAXSIZE;
    
        
    solution_iter = (argc > 2)? atoi(argv[2]) : MAXITER;
    
    if (coarse_dim > MAXSIZE) coarse_dim = MAXSIZE;
    
    if ((solution_iter > MAXITER)||(solution_iter <= 0)) solution_iter = MAXITER;
    
    
    //accomodate the boundary conditions in size
    coarse_dim_with_boundary = coarse_dim + 2;
    fine_dim = (coarse_dim*2)+1;
    fine_dim_with_boundary=fine_dim+2;
    
    //calculate the coarse grid with double spacing
    
    
    
    printf("\n\n******* Fine Grid Size: %d and Number of coarse iterations: %d *******\n\n", fine_dim, solution_iter);
    
    
    
    //create the matrices
    grid_fine = (float* )malloc(fine_dim_with_boundary * fine_dim_with_boundary * sizeof(float));
    newMatrix_fine = (float* )malloc(fine_dim_with_boundary * fine_dim_with_boundary * sizeof(float));
    grid_coarse = (float* )malloc(coarse_dim_with_boundary * coarse_dim_with_boundary * sizeof(float));
    newMatrix_coarse = (float* )malloc(coarse_dim_with_boundary * coarse_dim_with_boundary * sizeof(float));
    
    
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
    
    i = 0;
    j = 0;
    for (x=0; x<coarse_dim_with_boundary; x++)
    {
        
        for (y=0; y<coarse_dim_with_boundary; y++)
        {
            
            grid_coarse[x*coarse_dim_with_boundary+y] = grid_fine[i*fine_dim_with_boundary+j];
            newMatrix_coarse[x*coarse_dim_with_boundary+y] = grid_fine[i*fine_dim_with_boundary+j];
            j=j+2;
            
        }
        i=i+2;
        j=0;
    }
    
    
    
    
    time = timer();
    
    
    //******************************* STEP 1 SMOOTHING **********************************************
    
    // Step1:Smoothing via jacobi    
    
    //********************************************** //**********************************************    
    
    for(kk = 0; kk < v_cycles; kk++)
    {
        //printf("\nStep1 Smoothing on fine matrix: DONE\n");
        jacobi(grid_fine, newMatrix_fine, fine_dim_with_boundary, smoothing_iter);
        
        //printMatrix(grid_fine, fine_dim_with_boundary);
        
        
        
        //***************************** STEP 2 RESTRICTION **********************************************
        
        //step2: Restrict the fine grid to a coarser grid in which the points are twice as far apart
        //restriction operator
        //coarse[x][y] = fine[i][j]*0.5 + (fine[i-1][j] + fine[i][j-1] + fine[i][j+1] + fine[i+1][j])* 0.125
        
        
        //********************************************** //**********************************************
        
        restriction();
        
        //printf("\nStep2 coarse grid restriction: DONE\n");
        // printMatrix(grid_coarse, coarse_dim_with_boundary);
        
        
        
        //******************************** STEP 3 SOLUTION **********************************************
        
        //step3: compute the solution to desired accuracy
        
        //********************************************** //**********************************************
        
        

        jacobi(grid_coarse, newMatrix_coarse, coarse_dim_with_boundary, solution_iter);
        
        
        //printf("\n\n\nStep3 %d iterations on coarse: DONE\n", solution_iter);
        //printMatrix(grid_coarse, coarse_dim_with_boundary);
        
        //*************************** STEP 4 INTERPOLATION **********************************************
        
        //step4: Interpolate the coarse grid back to fine grid
        
        //********************************************** //**********************************************
        
        interpolate();        
        //printf("\n\n\nStep4 matrix Interpolation back to fine grid: DONE\n"); 
        //printMatrix(grid_fine, fine_dim_with_boundary);
        
        //****************************** STEP 5: SMOOTHING **********************************************
        
        //step5: update the fine grid for a few iterations
        
        //********************************************** //**********************************************
        
        jacobi(grid_fine, newMatrix_fine, fine_dim_with_boundary, smoothing_iter);
        
        //printf("\n\n\nStep5 Final Smoothing: DONE\n"); 
        // printMatrix(grid_fine, fine_dim_with_boundary);
    }
    
    Finalmaxdiff = 0.0;
    
    
    for (i=1; i<fine_dim_with_boundary; i++) 
        
    {
        
        for (j=1; j<fine_dim_with_boundary; j++) 
            
        {
            Finalmaxdiff = max(Finalmaxdiff, absolute(1 - grid_fine[(i*fine_dim_with_boundary)+j]));
            //Finalmaxdiff = max(Finalmaxdiff, absolute(newMatrix_fine[(i*fine_dim_with_boundary)+j] - grid_fine[(i*fine_dim_with_boundary)+j]));            
        }
    }
    
    printf("\nFinal maxdiff: %f after %d V-cycles\n\n",Finalmaxdiff, v_cycles);
    
    time=timer()-time;
    printf("Elapsed time: %f\n",time/1000000.0);
    fp=fopen("multigrid_serial_data.txt", "wb");
    
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
    
    printf("data saved in multigrid_serial_data.txt\n");
    fclose(fp);
    
    
    
    free(newMatrix_fine);
    free(grid_fine);
    free(grid_coarse);
    free(newMatrix_coarse);
    
    return 0;
    
}//end Main


void restriction(void)
{
    int i, j;
    int x,y;
    i = 2;
    j = 2;
    for (x=1; x<coarse_dim_with_boundary-1; x++)
    {
        
        for (y=1; y<coarse_dim_with_boundary-1; y++)
        {
            
            grid_coarse[x*coarse_dim_with_boundary+y] = (grid_fine[i*fine_dim_with_boundary+j]*0.5) + (grid_fine[((i-1)*fine_dim_with_boundary)+j] + grid_fine[((i+1)*fine_dim_with_boundary)+j] + grid_fine[(i*fine_dim_with_boundary)+j-1] + grid_fine[(i*fine_dim_with_boundary)+j+1])*0.125;
            j=j+2;
            
        }
        i=i+2;
        j=2;
    }
    
    
}


void interpolate(void)
{
    
    int i,j,x,y;
    i = 2;
    j = 2;
    for (x=1; x<coarse_dim_with_boundary; x++)
    {
        
        for (y=1; y<coarse_dim_with_boundary; y++)
        {
            
            grid_fine[i*fine_dim_with_boundary+j] = grid_coarse[x*coarse_dim_with_boundary+y];
            
            grid_fine[((i-1)*fine_dim_with_boundary)+j] = ((grid_fine[((i-2)*fine_dim_with_boundary)+j]) + (grid_fine[(i*fine_dim_with_boundary)+j]))*0.5;//above
            
            grid_fine[(i*fine_dim_with_boundary)+(j-1)] = ((grid_fine[(i*fine_dim_with_boundary)+(j-2)]) + (grid_fine[(i*fine_dim_with_boundary)+j]))*0.5;//left8
            
            grid_fine[((i-1)*fine_dim_with_boundary)+(j-1)] = ((grid_fine[((i-1)*fine_dim_with_boundary)+(j-2)]) + (grid_fine[((i-1)*fine_dim_with_boundary)+j]))*0.5;
            
            j=j+2;
            
        }
        i=i+2;
        j=2;
    }
    
}

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


