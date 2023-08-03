#include "biort.hpp"

ReactionNetwork::ReactionNetwork() {
    this->num_stc = 0;
    this->num_spc = 0;
    this->num_ssc = 0;
    this->num_sdc = 0;
    this->num_min = 0;
    this->num_ads = 0;
    this->num_cex = 0;
    this->num_mkr = 0;
    this->num_akr = 0;
    this->tmp = BADVAL;
    this->dep_mtx = { { 0.0 } };
    this->dep_kin = { { 0.0 } };
    this->conc_contrib = { { 0.0 } };
    this->keq = { 0.0 };
    this->keq_kin = {0.0};
    this->adh = BADVAL;
    this->bdh = BADVAL;
    this->bdt = BADVAL;
}