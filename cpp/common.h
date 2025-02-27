#include <vector>
#include <cmath>
using namespace std;

typedef enum {
	ALGO_CELLS,
	ALGO_LISTS,
	ALGO_REF_LISTS,
	ALGO_NONE,

} mdalgo_t;

extern float SIGMA;
extern float EPSILON;
extern float CUTOFF;
extern float M;

extern int UNIVERSE_SIZE;
extern int N_PARTICLE;
extern int N_TIMESTEP;
extern float DT;
extern int SEED;
extern char *path;
extern int RESOLUTION;
extern int NEIGHBOR_REFRESH_RATE;
extern int BR;
extern int BN;
extern int THREADS;
extern char *LOG_PATH;
extern mdalgo_t ALGO;

extern float LJ_MIN;
extern float R;

extern int t;

#define N_CELL (UNIVERSE_SIZE*UNIVERSE_SIZE*UNIVERSE_SIZE)
#define L (R * UNIVERSE_SIZE)

// The potential well is centered at (26/7)^(1/6)*SIGMA
#define R_MAX (SIGMA*powf(26/7,1/6))
#define LJ(r) 4*EPSILON*(6*powf(SIGMA,6)/powf(r,7)-12*powf(SIGMA,12)/powf(r,13))



#define DIR_MODE (S_IRWXU|S_IRWXG|S_IROTH)
#define FILE_MODE (S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH)

class vec {
public:
    float x;
    float y;
    float z;
    
    vec(float x, float y, float z);
    vec();

    vec operator+(const vec &other);
    vec operator*(const float c);
    vec operator%(const vec &other);
    
    vec& operator*=(const float c);
    vec& operator+=(const float c);
    vec& operator+=(const vec &other);

    vec& operator=(const vec &other);

    void apbc();
    int cell();

    float norm();
    float normsq();

    void read(float *buf);

#ifdef DEBUG
    void sprint(char *buf);
    void print();
    char *str();
private:
    char strbuf[32];
#endif
};

class particle {
public:
    particle();
    particle(vec r);

	int interact(particle *pn);	
	void update_position();

    vec r;
    vec v;
    int id;
    int cell;
    static int counter;
#ifdef DEBUG
    char *str();

private:
    char dbstr[16]; 
    int old_cell;
#endif

}; 

template <typename T>
class vector {
public:
	vector() : i(0), n(10) {
		array = aligned_alloc(sizeof(T), sizeof(T) * size);
	}

	void push_back(T x) {
		if (i == size) {
			T *tmp = aligned_alloc(sizeof(T), 2 * sizeof(T) * size);
			memcpy(tmp,array, sizeof(T) * size);
			free(array);
			array = tmp;
		}

		array[i++] = x;
	}

	int size() {
		return n;
	}

private:
	int i;
	int n;
	T *array;
}


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

int linear_idx(int, int, int);
void cubic_idx(int *, int);
float subm(float, float);
float lj(float);
float frand(void);
void thread(void (*)(int,int), int);
int parse_cli(int, char **);
void init_particles(vector<particle> &);
void save(vector<particle> &, int);
void save(vector<vector<particle>> &, int);

unsigned long rdtsc();
