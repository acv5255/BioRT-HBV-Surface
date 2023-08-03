#include "biort.hpp"

ChemicalState::ChemicalState() {
    this->tot_conc = { BADVAL };
    this->prim_conc = { BADVAL };
    this->sec_conc = { BADVAL };
    this->prim_conc = { BADVAL };
    this->tot_mol = { BADVAL };
    this->soil_parameters = SoilParameters();
}

ChemicalState ChemicalState::copy() {
    return ChemicalState(
        this->tot_conc,
        this->prim_conc,
        this->sec_conc,
        this->prim_actv,
        this->tot_mol,
        this->soil_parameters
    );
}