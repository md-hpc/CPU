//
// Created by sorak on 11/30/2024.
//

#ifndef CPU_MD_H
#define CPU_MD_H

#define EPSILON 40
#define SIGMA 1
#define UNIVERSE_SIZE  10
#define L (CUTOFF * UNIVERSE_SIZE)
#define CUTOFF (SIGMA * 2.5)
#define TIMESTEP (1*pow(10,-7))
#define SEED 246
#define DT 0.1
#define LJ_MIN (-4*24*EPSILON/SIGMA*(powf(7./26.,7./6.)-2*powf(7./26.,13./6.)))


#include <math.h>
#include <time.h>
#include <stdlib.h>



typedef struct{
    float x,y,z;
    float vX,vY,vZ;
    int particleId;
    //Cell Coordinate Here might be useful

}Particle;

typedef struct{
    float aX,aY,aZ;
    int particleId;
}Acceleration;

typedef struct{
    float vX,vY,vZ;
    int particleId;
}Velocity;
void init_ParticleList(Particle *particleList);
void init_AccelerationCache(Acceleration *accelerationCache);

/*
float distance(float a, float b){
    float opts [] = {((b - L) - a),
                  (b - a),
                  ((b + L) - a)};
    float m = abs(opts[0]);
    int mi = 0;
    for(int i = 0;i<3;i++){
        if(abs(opts[i]) < m){
            m = abs(opts[i]);
            mi = i;
        }
    }
    return opts[mi];
}
void distance3D (Particle reference, Particle neighbor, float * dist){
    dist[0] = distance(reference.x,neighbor.x);
    dist[1] = distance(reference.y,neighbor.y);
    dist[2] = distance(reference.z,reference.y);
}

void normalize(float *distance3D){
    double sum =0,mean,std_dev=0;

    int i;

    for(i = 0;i<3;i++){
        sum+= distance3D[i];
    }
    mean = sum/3;

    for(i = 0; i < 3; i++){
        std_dev += pow(distance3D[i] - mean,2);
    }
    std_dev = sqrt(std_dev/3);

    for(i = 0 ;i < 3;i++){
        distance3D[i] = (distance3D[i]-mean / std_dev);
    }

}
*/
float dummy_LJ_1D(float  ref, float  neighbor){
    float r = ref - neighbor;

    float LJ;
    if(fabs(r) > CUTOFF){
        return 0;
    }
    else if(r == 0){
        return 0;
    }
    else{
         LJ = 100/r ;
    }
    return LJ;

}

Acceleration dummy_LJ_3D(Particle *ref,Particle *neighbor){
    float aX = dummy_LJ_1D(ref->x,neighbor->x);
    float aY = dummy_LJ_1D(ref->y,neighbor->y);
    float aZ = dummy_LJ_1D(ref->z,neighbor->z);

    Acceleration output = {aX,aY,aZ,ref->particleId};
    return output;
}
float positionCheck(float coordinate){
    if (coordinate > UNIVERSE_SIZE){
        return coordinate - UNIVERSE_SIZE;
    }
    else if(coordinate<0){
        return coordinate + UNIVERSE_SIZE;
    }
    else{
        return coordinate;
    }
}
#define MAX_VELOCITY 1000
float velocityCheck(float velocity){
    if(velocity<0){
        return fmax(-MAX_VELOCITY,velocity);
    }
    else if(velocity>0){
        return fmin(MAX_VELOCITY,velocity);
    }
    else{
        return 0;
    }
}
float calculate_Velocity(float a2);
static inline float subm(float a,float b){
    float c,d;
    c = b-a;
    if(c>0){
        d=b-L-a;
        if(c<-d){
            return c;
        }
        else{
            return d;
        }
    }
    else{
        d=b+L-a;
        if(-c<d){
            return c;
        }
        else{
            return d;
        }
    }
}

static inline void modr(Particle *c, Particle *a,Particle *b){
    c->x=subm(a->x,b->x);
    c->y=subm(a->y,b->y);
    c->z=subm(a->z,b->z);
}

static inline float norm(Particle *r){
    return sqrt(powf(r->x,2)+powf(r->y,2)+powf(r->z,2));
}

static inline float lj(float r){
    float f=4*EPSILON*(6*powf(SIGMA,6)/powf(r,7)-12*powf(SIGMA,12)/powf(r,13));

    if(f<LJ_MIN){
        return LJ_MIN;
    }
    else{
        return f;
    }

} 

static inline void scalar_mul(Particle *v,float c){
    v->vX*=c;
    v->vY*=c;
    v->vY*=c;
}

static inline void vec_add(Particle *a, Particle *b){
    a->x+=b->x;
    a->y+=b->y;
    a->z+=b->z;
}

#endif //CPU_MD_H
