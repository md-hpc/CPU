#ifndef I_COMMON
#define I_COMMON 1

using namespace std;

#define R_MAX (SIGMA*powf(26/7,1/6))
#define LJ(r) 4*EPSILON*(6*powf(SIGMA,6)/powf(r,7)-12*powf(SIGMA,12)/powf(r,13))

class vec {
public:
    float x;
    float y;
    float z;
    
    vec(float x, float y, float z);
    vec();

    vec& operator=(const vec &other);

#ifdef DEBUG
    char *str();
private:
    char strbuf[32];
#endif
};

class timer {
public:
    timer();

    void start();
    void stop();
    unsigned long get();

private:
    unsigned long time;
    unsigned long last;
    bool running;
};

unsigned long rdtsc();
float frand();

#endif 
