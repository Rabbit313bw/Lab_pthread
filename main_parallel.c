#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>

#define DT 0.05
#define eps 0.0001

#include <sys/time.h>


#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
}


typedef struct
{
    double x, y;
} vector;


int bodies, timeSteps, thread_count;
double *masses, GravConstant;
double total_time = 0.0;
vector *positions, *velocities, *accelerations, *F_q_k;
pthread_t *thread_handles;
pthread_barrier_t barrier;

vector addVectors(vector a, vector b)
{
    vector c = {a.x + b.x, a.y + b.y};

    return c;
}

vector scaleVector(double b, vector a)
{
    vector c = {b * a.x, b * a.y};

    return c;
}

vector subtractVectors(vector a, vector b)
{
    vector c = {a.x - b.x, a.y - b.y};

    return c;
}

double mod(vector a)
{
    return sqrt(a.x * a.x + a.y * a.y);
}

void initiateSystem(char *fileName)
{
    int i;
    FILE *fp = fopen(fileName, "r");

    fscanf(fp, "%lf%d%d", &GravConstant, &bodies, &timeSteps);

    masses = (double *)malloc(bodies * sizeof(double));
    positions = (vector *)malloc(bodies * sizeof(vector));
    velocities = (vector *)malloc(bodies * sizeof(vector));
    accelerations = (vector *)malloc(bodies * sizeof(vector));

    for (i = 0; i < bodies; i++)
    {
        fscanf(fp, "%lf", &masses[i]);
        fscanf(fp, "%lf%lf", &positions[i].x, &positions[i].y);
        fscanf(fp, "%lf%lf", &velocities[i].x, &velocities[i].y);
    }

    fclose(fp);
}


void computeVelocities()
{
    int i;

    for (i = 0; i < bodies; i++)
        velocities[i] = addVectors(velocities[i], scaleVector(DT, accelerations[i]));
}


void computePositions()
{
    int i;

    for (i = 0; i < bodies; i++)
        positions[i] = addVectors(positions[i], scaleVector(DT,velocities[i]));
}

void *computeAccelerations(void* rank)
{
    int local_n = bodies / thread_count;
    long my_rank = (long)rank;
    int local_start = my_rank * local_n;
    int local_finish = local_start + local_n;
    if (my_rank == thread_count - 1)
        local_finish = bodies;
    int j;
    for (int i = local_start; i < local_finish; i++)
    {
        accelerations[i].x = 0;
        accelerations[i].y = 0;
        for (j = 0; j < local_finish; j++)
        {
            if (i != j)
            {
                double dist = pow(mod(subtractVectors(positions[i], positions[j])), 3);
                if (dist >= eps)
                {
                    F_q_k[i * bodies + j] = scaleVector(GravConstant * masses[j] / dist, subtractVectors(positions[j], positions[i]));
                    F_q_k[j * bodies + i] = scaleVector(-1, F_q_k[i * bodies + j]);
                }
                else{
                    F_q_k[i * bodies + j] = scaleVector(GravConstant * masses[j] / eps, subtractVectors(positions[j], positions[i]));
                    F_q_k[j * bodies + i] = scaleVector(-1, F_q_k[i * bodies + j]);
                }
                accelerations[i] = addVectors(accelerations[i], F_q_k[i * bodies + j]);
            }
        }
    }
    pthread_barrier_wait(&barrier);
    if (my_rank != thread_count - 1){
        for (int i = local_start; i < local_finish; i++)
        {
            for (int j = local_finish; j < bodies; j++)
            {
                if (i != j)
                {
                    accelerations[i] = addVectors(accelerations[i], F_q_k[i * bodies + j]);
                }
            }
        }
    }
    return NULL;
}

// void resolveCollisions()
// {
//     int i, j;

//     for (i = 0; i < bodies - 1; i++)
//         for (j = i + 1; j < bodies; j++)
//         {
//             if (positions[i].x == positions[j].x && positions[i].y == positions[j].y)
//             {
//                 vector temp = velocities[i];
//                 velocities[i] = velocities[j];
//                 velocities[j] = temp;
//             }
//         }
// }


void simulate()
{
    for (long i = 0; i < thread_count; ++i)
    {
        pthread_create(&thread_handles[i], NULL, computeAccelerations, (void *)i);
    }  
    for (long i = 0; i < thread_count; ++i)
    {
        pthread_join(thread_handles[i], NULL);
    }
    computePositions();
    computeVelocities();
}



int main(int argC, char *argV[])
{
    int i, j;

    if (argC != 3)
    {
        printf("%d", argC);
        printf("Usage : %s <file name containing system configuration data>", argV[0]);
    }
    else
    {
        thread_count = strtol(argV[2], NULL, 10);
        initiateSystem(argV[1]);

        
        thread_handles = malloc(thread_count * sizeof(pthread_t));


        F_q_k = (vector *)malloc(bodies * bodies * sizeof(vector));




        pthread_barrier_init(&barrier, NULL, thread_count);
        initiateSystem(argV[1]);
        //FILE *fw = fopen("output_par", "w");
        //fputs("Body   :     x              y           vx              vy   ", fw);
        double begin;
        GET_TIME(begin);
        char str_num_bodies[10];
        char str_num_timesteps[100];
        for (i = 0; i < timeSteps; i++)
        {
            char cycl[100];
            //printf("\nCycle %d\n", i + 1);
            //fputs(cycl, fw);
            simulate();
            for (j = 0; j < bodies; j++)
            {
                char str[10000];
                if (i == timeSteps - 1)
                    printf("Body %d : %lf\t%lf\t%lf\t%lf\n", j + 1, positions[j].x, positions[j].y, velocities[j].x, velocities[j].y);
                //fputs(str, fw);
            }
        }
        double end;
        GET_TIME(end);
        total_time = end - begin;
        printf("\nTime is : %lf\n", total_time);
        pthread_barrier_destroy(&barrier);
        free(thread_handles);
        free(masses);
        free(positions);
        free(velocities);
        free(F_q_k);
    }
    return 0;
}