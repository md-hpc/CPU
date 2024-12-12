typedef struct {
    float x;
    float y;
    float z;
} vector_t;

typedef struct {
    vector_t r;
    vector_t v;
} particle_t;

typedef struct {
    long n;
    particle_t particles[128];
} cell_t;

typedef struct {
    long cidx;
    particle_t p;
} outbound_particle_t;

typedef struct {
    long n;
    outbound_particle_t particles[32];
} outbound_t;

static inline long linear_idx(long i, long j, long k);
static inline void cubic_idx(long *res, long idx);
static inline long position_to_cell(vector_t *r);
static inline float subm(float a, float b);
static inline void scalar_mul(vector_t *v, float c);
static inline void vec_add(vector_t *a, vector_t* b);
static inline void modr(vector_t *c, vector_t *a, vector_t *b);
static inline float norm(vector_t *r);
static inline float lj(float r);
static inline float frand();
