#include "biort.hpp"

SoilParameters::SoilParameters() {
    ssa = { BADVAL };
    sw_thld = { BADVAL };
    sw_exp =  { BADVAL };
    q10 = { BADVAL };
    n_alpha = { BADVAL };
}

SoilParameters SoilParameters::copy() const {
    return SoilParameters(
        this->ssa,
        this->sw_thld,
        this->sw_exp,
        this->q10,
        this->n_alpha
    );
}