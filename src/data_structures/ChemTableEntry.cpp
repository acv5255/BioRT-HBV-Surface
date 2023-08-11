#include "biort.hpp"

ChemTableEntry::ChemTableEntry() {
    this->molar_mass = BADVAL;
    this->molar_vol = BADVAL;
    this->charge = BADVAL;
    this->size_fac = BADVAL;
    this->itype = -1;
    this->mtype = -1;
}

ChemTableEntry::ChemTableEntry(const string& name, const f64 molar_mass, const f64 molar_vol, const f64 charge, const f64 size_fac, const int itype, const int mtype) {
    this->name = name;
    this->molar_mass = molar_mass;
    this->molar_vol = molar_vol;
    this->charge = charge;
    this->size_fac = size_fac;
    this->itype = itype;
    this->mtype = mtype;
}

ChemTableEntry ChemTableEntry::copy() {
    return ChemTableEntry(
        this->name,
        this->molar_mass,
        this->molar_vol,
        this->charge,
        this->size_fac,
        this->itype,
        this->mtype
    );
}