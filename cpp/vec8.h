class vec8 {
public:
    fvec x;
    fvec y;
    fvec z;
    
    vec8(fvec x, fvec y, fvec z);
    vec8();

    vec8 operator+(const vec8 &other);
    vec8 operator*(const fvec c);
    vec8 operator%(const vec8 &other);
    
    vec8& operator*=(const fvec c);
    vec8& operator+=(const fvec c);
    vec8& operator+=(const vec8 &other);

    vec8& operator=(const vec8 &other);
};
