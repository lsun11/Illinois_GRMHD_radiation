#include <cassert>
#include <cmath>
#include <limits>

#include <utils.hh>

#define SIZE 10

int main(void) {
    double x[SIZE];
    double y[SIZE];
    double st;
    double sp;

    for(int i = 0; i < SIZE; ++i) {
        x[i] = i;
        y[i] = SIZE - i;
        st += x[i] * y[i];
        sp += x[i] * y[i] * (i % 2 == 0);
    }

    double dott = utils::dot<double, SIZE>::eval(x, 1, y, 1);
    double dotp = utils::dot<double, SIZE/2>::eval(x, 2, y, 2);

    assert(std::abs(st - dott) < std::numeric_limits<double>::epsilon());
    assert(std::abs(sp - dotp) < std::numeric_limits<double>::epsilon());
}
