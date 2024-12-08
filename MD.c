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
#define NUM_THREADS 16
void plot_particles(Particle *particleList,FILE *fp){

    for(int i = 0;i <NUM_PARTICLES_UNIVERSE;i++){
        fprintf(fp, "%lf %lf %lf\n",particleList[i].x,particleList[i].y,particleList[i].z);
    }
    fclose(fp);
}
int main(){
    double start,end;
    double cpu_time;

    start = omp_get_wtime();

    Particle *particleList = malloc(sizeof(Particle)*NUM_PARTICLES_UNIVERSE);
    Acceleration *accelerationCache = malloc(sizeof(Acceleration)*NUM_PARTICLES_UNIVERSE);
    init_ParticleList(particleList);

    FILE *fpBefore = fopen("Before.dat","w");
    plot_particles(particleList,fpBefore);
    printf("Printing Particle List!!\n --------------------\n");
    for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++) {
        printf("Particle:%d\n", particleList[i].particleId);
        printf("X: %f\n", particleList[i].x);
        printf("Y: %f\n", particleList[i].y);
        printf("Z: %f\n", particleList[i].z);
        printf("\n");
    }
    printf("Printing Particle Velocity!!\n --------------------\n");
    for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++) {
        printf("Particle:%d\n", particleList[i].particleId);
        printf("X: %f\n", particleList[i].vX);
        printf("Y: %f\n", particleList[i].vY);
        printf("Z: %f\n", particleList[i].vZ);
        printf("\n");
    }

    for(int i=0;i<NUM_TIMESTEP;i++) {

        init_AccelerationCache(accelerationCache);

        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int ref = 0; ref < NUM_PARTICLES_UNIVERSE; ref++) {
            for (int neighbor = 0; neighbor < NUM_PARTICLES_UNIVERSE; neighbor++) {
                if (particleList[ref].particleId == particleList[neighbor].particleId) {
                    continue;
                } else {
                    Acceleration curr = LJ_3D(&particleList[ref], &particleList[neighbor]);
                    accelerationCache[ref].aX += curr.aX*TIMESTEP;
                    accelerationCache[ref].aY += curr.aY*TIMESTEP;
                    accelerationCache[ref].aZ += curr.aZ*TIMESTEP;
                }
            }
        }
        /*printf("Printing Particle Acceleration!!\n --------------------\n");
        for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++) {
            printf("Particle:%d\n", accelerationCache[i].particleId);
            printf("X: %f\n", accelerationCache[i].aX);
            printf("Y: %f\n", accelerationCache[i].aY);
            printf("Z: %f\n", accelerationCache[i].aZ);
            printf("\n");
        }
        */
        for (int ref = 0; ref < NUM_PARTICLES_UNIVERSE; ref++) {
            particleList[ref].vX += accelerationCache[ref].aX * TIMESTEP;
            particleList[ref].vY += accelerationCache[ref].aY * TIMESTEP;
            particleList[ref].vZ += accelerationCache[ref].aZ * TIMESTEP;
            /*
            float newVX = particleList[ref].vX + accelerationCache[ref].aX * TIMESTEP;
            float newVY = particleList[ref].vY + accelerationCache[ref].aY * TIMESTEP;
            float newVZ = particleList[ref].vZ + accelerationCache[ref].aZ * TIMESTEP;


            particleList[ref].vX = velocityCheck(newVX);
            particleList[ref].vY = velocityCheck(newVY);
            particleList[ref].vZ = velocityCheck(newVZ);
             */
        }
        for (int ref = 0; ref < NUM_PARTICLES_UNIVERSE; ref++) {
            particleList[ref].x = fmod((particleList[ref].x + particleList[ref].vX * TIMESTEP), L);
            particleList[ref].y = fmod((particleList[ref].y + particleList[ref].vY * TIMESTEP), L);
            particleList[ref].z = fmod((particleList[ref].z + particleList[ref].vZ * TIMESTEP), L);


            /*
            float newX = fmod((particleList[ref].x + particleList[ref].vX * TIMESTEP), L);
            float newY = fmod((particleList[ref].y + particleList[ref].vY * TIMESTEP), L);
            float newZ = fmod((particleList[ref].z + particleList[ref].vZ * TIMESTEP), L);

            particleList[ref].x = positionCheck(newX);
            particleList[ref].y = positionCheck(newY);
            particleList[ref].z = positionCheck(newZ);
            */
        }
    }
    printf("Printing Particle List!!\n --------------------\n");
    for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++) {
        printf("Particle:%d\n", particleList[i].particleId);
        printf("X: %f\n", particleList[i].x);
        printf("Y: %f\n", particleList[i].y);
        printf("Z: %f\n", particleList[i].z);
        printf("\n");
    }
    printf("Printing Particle Velocity!!\n --------------------\n");
    for (int i = 0; i < NUM_PARTICLES_UNIVERSE; i++) {
        printf("Particle:%d\n", particleList[i].particleId);
        printf("X: %f\n", particleList[i].vX);
        printf("Y: %f\n", particleList[i].vY);
        printf("Z: %f\n", particleList[i].vZ);
        printf("\n");
    }
    FILE *fpAfter = fopen("After.dat","w");
    plot_particles(particleList,fpAfter);

    printf("LJ MIN: %f\n", LJ_MIN);

    end = omp_get_wtime();
    cpu_time =  (double)(end - start);

    printf("Time Taken: %f seconds\n",cpu_time);
    return 0;

}