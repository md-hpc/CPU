using namespace std;

#define N_CELL (UNIVERSE_SIZE*UNIVERSE_SIZE*UNIVERSE_SIZE)
#define L (R * UNIVERSE_SIZE)

// The potential well is centered at (26/7)^(1/6)*SIGMA
#define R_MAX (SIGMA*powf(26/7,1/6))
#define LJ(r) 4*EPSILON*(6*powf(SIGMA,6)/powf(r,7)-12*powf(SIGMA,12)/powf(r,13))

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
