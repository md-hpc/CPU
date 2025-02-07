#include <vector>
#include <stdio.h>

using namespace std;

int main() {
    vector<vector<int>> a;

    vector<int> *b;

    a.resize(16);


    b = &a[5];

    b->push_back(5);

    printf("%d\n",(*b)[0]);
}
