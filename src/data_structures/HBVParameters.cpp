#include "biort.hpp"

HBVParameters::HBVParameters() {
    k1 = BADVAL;
    k2 = BADVAL;
    maxbas = BADVAL;
    perc = BADVAL;
    sfcf = BADVAL;
    tt = BADVAL;
}

HBVParameters HBVParameters::copy() {
    return HBVParameters(
        this->k1,
        this->k2,
        this->maxbas,
        this->perc,
        this->sfcf,
        this->tt
    );
}