#include "biort.hpp"

SoilConstants::SoilConstants() {
    porosity = BADVAL;
    depth = BADVAL;
    ws_passive = BADVAL;
}

SoilConstants SoilConstants::copy() {
    return SoilConstants(
        this->porosity,
        this->depth,
        this->ws_passive
    );
}