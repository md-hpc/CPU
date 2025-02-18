#include <cstdio>
#include <atomic>

using namespace std;

atomic<unsigned long> a{0};

int main() {
    a = 15;

    a.fetch_add(5);
    
    unsigned long b = a;
    printf("%d\n",b);
}
