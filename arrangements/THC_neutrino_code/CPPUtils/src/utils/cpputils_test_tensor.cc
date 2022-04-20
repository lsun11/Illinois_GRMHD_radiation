#include <cassert>
#include <cmath>
#include <iostream>

#ifndef CPPUTILS_DEBUG
#define CPPUTILS_DEBUG
#include <utils.hh>
#endif

#define DELTA(a,b) ((a)==(b) ? 1 : 0)
#define TEST_EPSILON 1e-14
//#define TEST_VERBOSE

int main(void) {
    // Test rank 2 symmetric tensor
    utils::tensor::symmetric2<CCTK_REAL, 4, 2> A;
    for(int i = 0; i < utils::tensor::symmetric2<CCTK_REAL, 4, 2>::ndof; ++i) {
        A[i] = i;
    }
    for(int a = 0; a < 4; ++a)
    for(int b = 0; b < 4; ++b) {
        assert(A(a,b) == A(b,a));
    }
#ifdef TEST_VERBOSE
    std::cout << std::endl << "symmetric2<CCTK_REAL, 4, 2>" << std::endl;
    for(int a = 0; a < 4; ++a) {
        for(int b = 0; b < 4; ++b) {
            int const i = static_cast<int>(&A(a,b) - &A[0]);
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << "symmetric2<CCTK_REAL, 3, 2>" << std::endl;
    utils::tensor::symmetric2<3, 2> B;
    for(int a = 0; a < 3; ++a) {
        for(int b = 0; b < 3; ++b) {
            int const i = static_cast<int>(&B(a,b) - &B[0]);
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
#endif

    // Test rank 3 symmetric tensor
    utils::tensor::symmetric2<CCTK_REAL, 4, 3> Gamma;
    for(int i = 0; i < utils::tensor::symmetric2<CCTK_REAL, 4, 3>::ndof; ++i) {
        Gamma[i] = i;
    }
    for(int a = 0; a < 4; ++a)
    for(int b = 0; b < 4; ++b)
    for(int c = 0; c < 4; ++c) {
        assert(Gamma(a,b,c) == Gamma(a,c,b));
    }
#ifdef TEST_VERBOSE
    std::cout << std::endl << "symmetric2<CCTK_REAL, 4, 3>" << std::endl;
    for(int a = 0; a < 4; ++a) {
        for(int b = 0; b < 4; ++b) {
            for(int c = 0; c < 4; ++c) {
                int const i = static_cast<int>(&Gamma(a,b,c) - &Gamma[0]);
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
#endif

    // Test spacetime metric and inverse
    CCTK_REAL const alp = 0.8;
    CCTK_REAL const betax = 0.1;
    CCTK_REAL const betay = 0.0;
    CCTK_REAL const betaz = -0.1;
    CCTK_REAL const gxx = 0.9;
    CCTK_REAL const gxy = 0.6;
    CCTK_REAL const gxz = 0.0;
    CCTK_REAL const gyy = 1.0;
    CCTK_REAL const gyz = 0.1;
    CCTK_REAL const gzz = 1.0;
    utils::tensor::metric<4> g;
    g.from_adm(alp, betax, betay, betaz, gxx, gxy, gxz, gyy, gyz, gzz);
    utils::tensor::inv_metric<4> u;
    u.from_adm(alp, betax, betay, betaz, gxx, gxy, gxz, gyy, gyz, gzz);

    utils::tensor::symmetric2<CCTK_REAL, 4, 2> id;
    for(int a = 0; a < 4; ++a) {
        for(int b = a; b < 4; ++b) {
            id(a, b) = 0.0;
            for(int c = 0; c < 4; ++c) {
                id(a, b) += g(a, c) * u(c, b);
            }
            assert(std::abs(id(a,b) - DELTA(a,b)) < TEST_EPSILON);
        }
    }
}
