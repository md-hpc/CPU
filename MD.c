//
// Created by sorak on 11/30/2024.
//

#include "MD.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#define NUM_PARTICLES_UNIVERSE 10000
#define NUM_TIMESTEP 1
#define NUM_THREADS 16
#define BSIZE 1
#define BLOCKED 1

void velocityUpdateBlocked(Particle *particleList);
void velocityUpdate(Particle *particleList);
void printData(Particle *particleList);
void plot_particles(Particle *particleList,FILE *fp);

int main()
{
    double start, end;
    double cpu_time;

    

    Particle *particleList = malloc(sizeof(Particle) * NUM_PARTICLES_UNIVERSE);
    init_ParticleList(particleList);
    Particle current;

    //FILE *fpBefore = fopen("Before.dat", "w");
    //plot_particles(particleList, fpBefore);

    
    start = omp_get_wtime();
    for (int i = 0; i < NUM_TIMESTEP; i++)
    {

        if (BLOCKED)
        {
            velocityUpdateBlocked(particleList);
           
        }
        else
        {
            velocityUpdate(particleList);
        }

        for (int ref = 0; ref < NUM_PARTICLES_UNIVERSE; ref++)
        {
            current = particleList[ref];

            current.x = fmod(current.x + current.vX * DT, L);
            current.y = fmod(current.y + current.vY * DT, L);
            current.z = fmod(current.z + current.vZ * DT, L);

            current.x = current.x < 0 ? current.x + L : current.x;
            current.y = current.y < 0 ? current.y + L : current.y;
            current.z = current.z < 0 ? current.z + L : current.z;

            current.x = current.x > L ? current.x - L : current.x;
            current.y = current.y > L ? current.y - L : current.y;
            current.z = current.z > L ? current.z - L : current.z;
        }
    }
    end = omp_get_wtime();
    cpu_time = (double)(end - start);
    //FILE *fpAfter = fopen("After.dat", "w");
    //plot_particles(particleList, fpAfter);


    printf("Time Taken: %f seconds\n", cpu_time);
    return 0;
}
void init_ParticleList(Particle *particleList)
{
    srand(SEED);
    int min = -10;
    int max = 10;

    for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++)
    {
        particleList[i].x = ((float)rand() / (float)(RAND_MAX)) * UNIVERSE_SIZE;
        particleList[i].y = ((float)rand() / (float)(RAND_MAX)) * UNIVERSE_SIZE;
        particleList[i].z = ((float)rand() / (float)(RAND_MAX)) * UNIVERSE_SIZE;

        particleList[i].vX = ((float)(min + rand() % (max - min + 1)) / (float)max) * UNIVERSE_SIZE / 2;
        particleList[i].vY = ((float)(min + rand() % (max - min + 1)) / (float)max) * UNIVERSE_SIZE / 2;
        particleList[i].vZ = ((float)(min + rand() % (max - min + 1)) / (float)max) * UNIVERSE_SIZE / 2;

        particleList[i].particleId = i;
    }
}
void velocityUpdateBlocked(Particle *particleList)
{
    Particle vec;
    float r, f;
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int ii = 0; ii < NUM_PARTICLES_UNIVERSE; ii += BSIZE)
    {
        for (int jj = 0; jj < NUM_PARTICLES_UNIVERSE; jj += BSIZE)
        {
            for (int ref = ii; ref < ii + BSIZE; ref++)
            {
                for (int neighbor = jj; neighbor < jj + BSIZE; neighbor++)
                {
                    if (ref == neighbor)
                    {
                        continue;
                    }
                    else
                    {
                        modr(&vec, &particleList[ref], &particleList[neighbor]);
                        r = norm(&vec);
                        if (r > CUTOFF)
                        {
                            continue;
                        }
                        if (vec.x < 0)
                        {
                            continue;
                        }
                        f = lj(r);
                        scalar_mul(&vec, DT * f / r);
                        vec_add(&particleList[ref], &vec);
                        scalar_mul(&vec, -1.);
                        vec_add(&particleList[neighbor], &vec);
                    }
                }
            }
        }
    }
}
void velocityUpdate(Particle *particleList)
{
    Particle vec;
    float r, f;
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int ref = 0; ref < NUM_PARTICLES_UNIVERSE; ref++)
    {
        for (int neighbor = 0; neighbor < NUM_PARTICLES_UNIVERSE; neighbor++)
        {
            
            if (ref == neighbor)
            {
                continue;
            }
            else
            {

                modr(&vec, &particleList[ref], &particleList[neighbor]);
                r = norm(&vec);

                if (r > CUTOFF)
                {
                    continue;
                }

                if (vec.x < 0)
                {
                    continue;
                }

                f = lj(r);
                scalar_mul(&vec, DT * f / r);
                vec_add(&particleList[ref], &vec);
                scalar_mul(&vec, -1.);
                vec_add(&particleList[neighbor], &vec);
                
            }
        }
    }
}
void printData(Particle *particleList){
    printf("Printing Particle List!!\n --------------------\n");
    for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++)
    {
        printf("Particle:%d\n", particleList[i].particleId);
        printf("X: %f\n", particleList[i].x);
        printf("Y: %f\n", particleList[i].y);
        printf("Z: %f\n", particleList[i].z);
        printf("\n");
    }
    printf("Printing Particle Velocity!!\n --------------------\n");
    for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++)
    {
        printf("Particle:%d\n", particleList[i].particleId);
        printf("X: %f\n", particleList[i].vX);
        printf("Y: %f\n", particleList[i].vY);
        printf("Z: %f\n", particleList[i].vZ);
        printf("\n");
    }
}
void plot_particles(Particle *particleList, FILE *fp)
{

    for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++)
    {
        fprintf(fp, "%lf %lf %lf\n", particleList[i].x, particleList[i].y, particleList[i].z);
    }
    fclose(fp);
}

