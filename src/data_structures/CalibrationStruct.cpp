#include "biort.hpp"

CalibrationStruct::CalibrationStruct() {
    xsorption = BADVAL;
    rate = BADVAL;
    ssa = BADVAL;
}

CalibrationStruct CalibrationStruct::copy() {
    return CalibrationStruct(
        this->xsorption,
        this->rate,
        this->ssa
    );
}