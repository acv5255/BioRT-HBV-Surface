#include "biort.hpp"

ControlData::ControlData() {
    this->recycle = 0;
    this->use_activity = 1;
    this->transport_only = 0;
    this->variable_precipchem = 1;
    this->precipchem_numexp = 1;
}

ControlData ControlData::copy() {
    return ControlData(
        this->recycle,
        this->use_activity,
        this->transport_only,
        this->variable_precipchem,
        this->precipchem_numexp
    );
}