#include <cassert>
#include <cmath>
#include <limits>

#include <utils.hh>

int main(void) {
    int A[3][3] = {{1, 2, 3},
                   {4, 5, 6},
                   {7, 8, 9}};
    int x[3]    = {1, 2, 3};
    int y[3];
    int z[3];

    utils::gemv<int, false, 3, 3>::eval(1, &A[0][0], 3, &x[0], 1, 0, &y[0], 1);
    utils::gemv<int, true,  3, 3>::eval(1, &A[0][0], 3, &x[0], 1, 0, &z[0], 1);

    assert(y[0] == 14);
    assert(y[1] == 32);
    assert(y[2] == 50);

    assert(z[0] == 30);
    assert(z[1] == 36);
    assert(z[2] == 42);

    int B[4][3] = {{ 1,  2,  3},
                   { 4,  5,  6},
                   { 7,  8,  9},
                   {10, 11, 12}};
    int u[4];

    utils::gemv<int, false, 4, 3>::simple(&B[0][0], &x[0], &u[0]);

    assert(u[0] == 14);
    assert(u[1] == 32);
    assert(u[2] == 50);
    assert(u[3] == 68);

    int C[3][4] = {{1,  2,  3, 4},
                   {5,  6,  7, 8},
                   {9, 10, 11, 12}};

    utils::gemv<int, true, 3, 4>::simple(&C[0][0], &x[0], &u[0]);

    assert(u[0] == 38);
    assert(u[1] == 44);
    assert(u[2] == 50);
    assert(u[3] == 56);

    double D[15][10];
    for(int i = 0; i < 15; ++i) {
        for(int j = 0; j < 10; ++j) {
            D[i][j] = i*j;
        }
    }

    double w[10];
    for(int j = 0; j < 10; ++j) {
        w[j] = j*j;
    }

    double Dw[15];
    utils::gemv<double, false, 15, 10>::simple(&D[0][0], &w[0], &Dw[0]);

    for(int i = 0; i < 15; ++i) {
        double check = 0;
        for(int j = 0; j < 10; ++j) {
            check += D[i][j]*w[j];
        }
        assert(std::abs(check - Dw[i]) < std::numeric_limits<double>
                ::epsilon());
    }
};
