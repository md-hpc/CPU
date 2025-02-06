#include <cmath>
using namespace std;
 
extern float SIGMA;
extern float EPSILON;
extern float CUTOFF;
extern int UNIVERSE_SIZE;
extern int N_PARTICLE;
extern int N_TIMESTEP;
extern float DT;
extern int SEED;
extern char *path;
extern int RESOLUTION;

#define N_CELL (UNIVERSE_SIZE*UNIVERSE_SIZE*UNIVERSE_SIZE)
#define L (CUTOFF * UNIVERSE_SIZE)

#define LJ_MIN (-4*24*EPSILON/SIGMA*(powf(7./26.,7./6.)-2*powf(7./26.,13./6.)))

#define DIR_MODE (S_IRWXU|S_IRWXG|S_IROTH)
#define FILE_MODE (S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH)

int linear_idx(int, int, int);
void cubic_idx(int *, int);
float subm(float, float);
float lj(float);
float frand(void);
void thread(void (*)(long), long);
int parse_cli(int, char **);

class vec {
public:
    float x;
    float y;
    float z;
    
    vec(float x, float y, float z);
    vec();

    vec operator+(const vec &other);
    vec operator*(const float c);

    vec& operator*=(const float c);
    vec& operator+=(const float c);
    vec& operator+=(const vec &other);

    vec& operator=(const vec &other);
    
    float norm();

    void read(float *buf);
};

class mod_vec : public vec {
public:
    float x;
    float y;
    float z;

    mod_vec(float x, float y, float z);
    mod_vec &operator+=(const vec &other);
    mod_vec operator+(const vec &other);
    void add(const vec &other);
    int cell();
    mod_vec modr(mod_vec& other);

private:
    void normalize();
};

class particle {
public:
    particle();
    particle(mod_vec r);
    
    int cell();

    mod_vec r;
    vec v;
    int idx;
}; 
