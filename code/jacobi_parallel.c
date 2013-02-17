// parallel implementation of jacobi using pthreads
// compile: gcc multigrid_parallel.c time.c -lpthread -o jacobi_parallel
// execute: ./jabobi_parallel size_of_grid num_threads num_iterations
// example: ./multi_parallel 200 3 500

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#define MAXSIZE 200  /* maximum fine matrix size */
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


int solution_iter;

int fine_dim;

int fine_dim_with_boundary;


float *grid_fine;
float *newMatrix_fine;
float *maxDiffArrary;



void printMatrix(float *matrix, int fine_dim);
float max(float a, float b);
float absolute(float a);
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
    fine_dim = (argc > 1)? atoi(argv[1]) : MAXSIZE;
    
    numWorkers = (argc > 2)? atoi(argv[2]) : MAXWORKERS;
    
    solution_iter = (argc > 3)? atoi(argv[3]) : MAXITER;
    
    if (fine_dim > MAXSIZE) fine_dim = MAXSIZE;
    
    if ((solution_iter > MAXITER)||(solution_iter <= 0)) solution_iter = MAXITER;
    
    if (numWorkers > MAXWORKERS) numWorkers = MAXWORKERS;
    
    
    
    //accomodate the boundary conditions in size
   
    fine_dim_with_boundary=fine_dim+2;
    
    
    printf("\n\n******* Fine Grid Size: %d and Number of fine iterations: %d and Number of Threads: %d *******\n\n", fine_dim, solution_iter, numWorkers);
    
    
    
    //set global thread attributes 
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    pthread_cond_init(&mysignal,NULL);
	pthread_mutex_init(&lock,NULL);    
    
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

        
        
        //create the workers, then exit 
        for (i = 0; i < numWorkers; i++)
            pthread_create(&workerid[i], &attr, jacobiP, (void *) i);
        
        
        
        for (i=0; i<numWorkers; i++)
            pthread_join(workerid[i],NULL);
    
    
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
    fp=fopen("jacobi_parallel_data.txt", "wb");
    
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
    
    printf("data saved in jacobi_parallel_data.txt\n");
    fclose(fp);
    
    
    
    free(newMatrix_fine);
    free(grid_fine);
    
    
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
    
    residue = fine_dim % numWorkers;
    height = (fine_dim-residue) / numWorkers;
    myFirst = (thread_id * height) +1;
    myLast = (thread_id * height)+height;
    if(thread_id==numWorkers-1)
        myLast = myLast+residue;
    
    
    
    for (kk=0; kk<(solution_iter/2); kk++) 
    {
        
        //compute the new inner points
        for (i=myFirst; i<=myLast; i++) 
            
        {
            for (j=1; j<fine_dim_with_boundary-1; j++) 
                
            {
                
                newMatrix_fine[(i*fine_dim_with_boundary)+j] = (grid_fine[((i-1)*fine_dim_with_boundary)+j] + grid_fine[((i+1)*fine_dim_with_boundary)+j] + grid_fine[(i*fine_dim_with_boundary)+j-1] + grid_fine[(i*fine_dim_with_boundary)+j+1])*0.25;
                
            }
            
        }
        
        
        Barrier();
        //loop unrolling
        for (i=myFirst; i<=myLast; i++) 
            
        {
            for (j=1; j<fine_dim_with_boundary-1; j++) 
                
            {
                
                grid_fine[(i*fine_dim_with_boundary)+j] = (newMatrix_fine[((i-1)*fine_dim_with_boundary)+j] + newMatrix_fine[((i+1)*fine_dim_with_boundary)+j] + newMatrix_fine[(i*fine_dim_with_boundary)+j-1] + newMatrix_fine[(i*fine_dim_with_boundary)+j+1])*0.25;
                
            }
            
        }
        
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

