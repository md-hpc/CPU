#include <stdio.h>
#include "common.h"

int main() {
    CUTOFF = 1.;
    UNIVERSE_SIZE = 3;

    vec a = vec(.5,.5,.5);
    vec i = vec(1,0,0);
    vec j = vec(0,1,0);
    vec k = vec(0,0,1);

    vec b = a + i*2;

    b.print();
    vec c = a % b;
    c.print();

}
