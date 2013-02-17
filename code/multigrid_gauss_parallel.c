// parallel implementation of multigrid method using pthreads
// compile: gcc multigrid_parallel.c time.c -lpthread -o multi_parallel
// execute: ./multi_parallel size_of_fine num_threads num_iterations
// example: ./multi_parallel 5 3 50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#define MAXSIZE 99  /* maximum coarse matrix size */
#define MAXITER 1000   // Maximum number of iterations
#define MAXWORKERS 8   /* maximum number of workers */
#define EPSILON 0.005

//*********************global variables and function definitions***************

pthread_mutex_t lock;
pthread_cond_t mysignal;
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
float *grid_coarse;
float *maxDiffArrary;



void printMatrix(float *matrix, int fine_dim);
float max(float a, float b);
float absolute(float a);
void jacobi(float *old, int sent_dim, int num_iter); 
void restriction (void);
void interpolate(void);
void *jacobiP(void *);

//*********************barrier function****************

void Barrier(){  
    int mystate; 
    pthread_mutex_lock (&lock);
    mystate=state;
    waiting++;
    if (waiting==numWorkers){
		waiting=0; state=1-mystate;
        pthread_cond_broadcast(&mysignal);}
    while (mystate==state)
		pthread_cond_wait(&mysignal,&lock);
    pthread_mutex_unlock (&lock);
}
//***************************************************************************

int main (int argc, const char * argv[]) 

{
    
    int i, j, kk;
    
    int x, y;
    
    
    float maxdiff;
    float Finalmaxdiff = 0.0;
    
    float time;
    FILE *fp;
    
    pthread_attr_t attr;
    pthread_t workerid[MAXWORKERS];
    
    //get command line arguments
    coarse_dim = (argc > 1)? atoi(argv[1]) : MAXSIZE;
    
    numWorkers = (argc > 2)? atoi(argv[2]) : MAXWORKERS;
    
    solution_iter = (argc > 3)? atoi(argv[3]) : MAXITER;
    
    if (coarse_dim > MAXSIZE) coarse_dim = MAXSIZE;
    
    if ((solution_iter > MAXITER)||(solution_iter <= 0)) solution_iter = MAXITER;
    
    if (numWorkers > MAXWORKERS) numWorkers = MAXWORKERS;
    
    
    
    //accomodate the boundary conditions in size
    coarse_dim_with_boundary = coarse_dim + 2;
    fine_dim = (coarse_dim*2)+1;
    fine_dim_with_boundary=fine_dim+2;
    
    //calculate the coarse grid with double spacing
    
    
    
    printf("\n\n******* Fine Grid Size: %d and Number of coarse iterations: %d and Number of Threads: %d *******\n\n", fine_dim, solution_iter, numWorkers);
    
    
    
    //set global thread attributes 
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    pthread_cond_init(&mysignal,NULL);
	pthread_mutex_init(&lock,NULL);    
    
    //create the matrices
    grid_fine = (float* )malloc(fine_dim_with_boundary * fine_dim_with_boundary * sizeof(float));

    grid_coarse = (float* )malloc(coarse_dim_with_boundary * coarse_dim_with_boundary * sizeof(float));

    
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
        
        
    }
    
    i = 0;
    j = 0;
    for (x=0; x<coarse_dim_with_boundary; x++)
    {
        
        for (y=0; y<coarse_dim_with_boundary; y++)
        {
            
            grid_coarse[x*coarse_dim_with_boundary+y] = grid_fine[i*fine_dim_with_boundary+j];
            j=j+2;
            
        }
        i=i+2;
        j=0;
    }
    
    
    
    
    time = timer();
    
    
        
    for(kk = 0; kk < v_cycles; kk++)
    {
        
        //******************************* STEP 1 SMOOTHING **********************************************
        
        // Step1:Smoothing via jacobi    
        
        //********************************************** //**********************************************    

        //printf("\nStep1 Smoothing on fine matrix: DONE\n");
        jacobi(grid_fine, fine_dim_with_boundary, smoothing_iter);
        
        //printMatrix(grid_fine, fine_dim_with_boundary);
        
        
        
        //***************************** STEP 2 RESTRICTION **********************************************
        
        //step2: Restrict the fine grid to a coarser grid in which the points are twice as far apart
        //coarse[x][y] = fine[i][j]*0.5 + (fine[i-1][j] + fine[i][j-1] + fine[i][j+1] + fine[i+1][j])* 0.125
        
        
        //********************************************** //**********************************************
        
        restriction();
        
       // printf("\nStep2 coarse grid restriction: DONE\n");
        // printMatrix(grid_coarse, coarse_dim_with_boundary);
        
        
        
        //******************************** STEP 3 SOLUTION **********************************************
        
        //step3: compute the solution to desired accuracy
        
        //********************************************** //**********************************************
        
        
        //create the workers, then exit 
        for (i = 0; i < numWorkers; i++)
            pthread_create(&workerid[i], &attr, jacobiP, (void *) i);
        
        
        
        for (i=0; i<numWorkers; i++)
            pthread_join(workerid[i],NULL);
        
        
       // printf("\n\n\nStep3 %d iterations on coarse: DONE\n", solution_iter);
        
        //*************************** STEP 4 INTERPOLATION **********************************************
        
        //step4: Interpolate the coarse grid back to fine grid
        
        //********************************************** //**********************************************
        
        interpolate();        
        //printf("\n\n\nStep4 matrix Interpolation back to fine grid: DONE\n"); 
        //printMatrix(grid_fine, fine_dim_with_boundary);
        
        //****************************** STEP 5: SMOOTHING **********************************************
        
        //step5: update the fine grid for a few iterations
        
        //********************************************** //**********************************************
        
        jacobi(grid_fine, fine_dim_with_boundary, smoothing_iter);
        
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
    fp=fopen("multigrid_gauss_parallel_data.txt", "wb");
    
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
    
    printf("data saved in multigrid_gauss_parallel_data.txt\n");
    fclose(fp);
    
    
    
    free(grid_fine);
    free(grid_coarse);
    
    pthread_attr_destroy(&attr);
    
    pthread_exit(NULL);
    
    
    return 0;
    
}//end Main




void *jacobiP(void *arg)
{
    
    int height, myFirst, myLast, residue;
    int i,j;
    int x,y;
    int kk;
    
    float diff = 0.0;
    float myDiff = 0.0;
    float largest = 1.0;
    
    int thread_id = (int) arg;
    
    
    /************************ Work Distribution ****************************/
    
    residue = coarse_dim % numWorkers;
    height = (coarse_dim-residue) / numWorkers;
    myFirst = (thread_id * height) +1;
    myLast = (thread_id * height)+height;
    if(thread_id==numWorkers-1)
        myLast = myLast+residue;
    
    
    
    for (kk=0; kk<solution_iter; kk++) 
    {
        
        //compute the new inner points
        for (i=myFirst; i<=myLast; i++) 
            
        {
            for (j=1; j<coarse_dim_with_boundary-1; j++) 
                
            {
                
                grid_coarse[(i*coarse_dim_with_boundary)+j] = (grid_coarse[((i-1)*coarse_dim_with_boundary)+j] + grid_coarse[((i+1)*coarse_dim_with_boundary)+j] + grid_coarse[(i*coarse_dim_with_boundary)+j-1] + grid_coarse[(i*coarse_dim_with_boundary)+j+1])*0.25;
                
            }
            
        }
        
        
        Barrier();
                
    }
    
    
}


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


void jacobi(float *old, int sent_dim, int num_iter) 
{
    
    int j;
    int i, k;
    
    for (k=0; k<num_iter; k++) 
    {
        
        //compute the new inner points
        for (i=1; i<sent_dim-1; i++) 
            
        {
            for (j=1; j<sent_dim-1; j++) 
                
                
            {
                
                old[(i*sent_dim)+j] = (old[((i+1)*sent_dim)+j] + old[((i-1)*sent_dim)+j] + old[(i*sent_dim)+j-1] + old[(i*sent_dim)+j+1])*0.25;
                
            }
            
            
        }
        
    }
    
}


