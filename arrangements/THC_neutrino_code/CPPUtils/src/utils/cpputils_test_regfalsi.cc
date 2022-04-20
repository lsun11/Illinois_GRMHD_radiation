#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include <utils.hh>

class Function: public utils::unary_function<CCTK_REAL, CCTK_REAL> {
    public:
        CCTK_REAL operator()(CCTK_REAL x) const {
            return std::atan(10*x);
        }
};

int main(void) {
    Function f;
    CCTK_REAL xp, res;
    int iter;

    utils::regfalsi::status const err = utils::regfalsi::rootfinder(
            f, -1.0, M_PI, std::numeric_limits<double>::epsilon(),
            100, xp, res, iter);

    assert(err == utils::regfalsi::ok);
}
