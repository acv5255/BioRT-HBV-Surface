#include "biort.hpp"

void CopyConstSubcatchProp(const Subcatchment& subcatch, Subcatchment& subcatch_numexp)
{
    subcatch_numexp.hbv_parameters.k1 = subcatch.hbv_parameters.k1;
    subcatch_numexp.hbv_parameters.k2 = subcatch.hbv_parameters.k2;
    subcatch_numexp.hbv_parameters.maxbas = subcatch.hbv_parameters.maxbas;
    subcatch_numexp.hbv_parameters.perc = subcatch.hbv_parameters.perc;
    subcatch_numexp.hbv_parameters.tt = subcatch.hbv_parameters.tt;
    subcatch_numexp.hbv_parameters.sfcf = subcatch.hbv_parameters.sfcf;
    
    subcatch_numexp.soil_surface = subcatch.soil_surface;
    subcatch_numexp.soil_sz = subcatch.soil_sz;
    subcatch_numexp.soil_dz = subcatch.soil_dz;
}

void CopyInitChemSubcatch(ReactionNetwork *rttbl, const Subcatchment& subcatch, Subcatchment& subcatch_numexp)
{
    for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
    {
        subcatch_numexp.chms[SNOW].soil_parameters.ssa[kspc] = subcatch.chms[SNOW].soil_parameters.ssa[kspc];
        subcatch_numexp.chms[SURFACE].soil_parameters.ssa[kspc] = subcatch.chms[SURFACE].soil_parameters.ssa[kspc];
        subcatch_numexp.chms[UZ].soil_parameters.ssa[kspc] = subcatch.chms[UZ].soil_parameters.ssa[kspc];
        subcatch_numexp.chms[LZ].soil_parameters.ssa[kspc] = subcatch.chms[LZ].soil_parameters.ssa[kspc];
        
        subcatch_numexp.chms[SURFACE].soil_parameters.q10[kspc] = subcatch.chms[SURFACE].soil_parameters.q10[kspc];
        subcatch_numexp.chms[UZ].soil_parameters.q10[kspc] = subcatch.chms[UZ].soil_parameters.q10[kspc];
        subcatch_numexp.chms[LZ].soil_parameters.q10[kspc] = subcatch.chms[LZ].soil_parameters.q10[kspc];
        
        subcatch_numexp.chms[SURFACE].soil_parameters.sw_thld[kspc] = subcatch.chms[SURFACE].soil_parameters.sw_thld[kspc];
        subcatch_numexp.chms[UZ].soil_parameters.sw_thld[kspc] = subcatch.chms[UZ].soil_parameters.sw_thld[kspc];
        subcatch_numexp.chms[LZ].soil_parameters.sw_thld[kspc] = subcatch.chms[LZ].soil_parameters.sw_thld[kspc];
        
        subcatch_numexp.chms[SURFACE].soil_parameters.sw_exp[kspc] = subcatch.chms[SURFACE].soil_parameters.sw_exp[kspc];
        subcatch_numexp.chms[UZ].soil_parameters.sw_exp[kspc] = subcatch.chms[UZ].soil_parameters.sw_exp[kspc];
        subcatch_numexp.chms[LZ].soil_parameters.sw_exp[kspc] = subcatch.chms[LZ].soil_parameters.sw_exp[kspc];
        
        subcatch_numexp.chms[SURFACE].soil_parameters.n_alpha[kspc] = subcatch.chms[SURFACE].soil_parameters.n_alpha[kspc];
        subcatch_numexp.chms[UZ].soil_parameters.n_alpha[kspc] = subcatch.chms[UZ].soil_parameters.n_alpha[kspc];
        subcatch_numexp.chms[LZ].soil_parameters.n_alpha[kspc] = subcatch.chms[LZ].soil_parameters.n_alpha[kspc];

        subcatch_numexp.chms[SNOW].prim_conc[kspc] = subcatch.chms[SNOW].prim_conc[kspc];
        subcatch_numexp.chms[SURFACE].prim_conc[kspc] = subcatch.chms[SURFACE].prim_conc[kspc];
        subcatch_numexp.chms[UZ].prim_conc[kspc] = subcatch.chms[UZ].prim_conc[kspc];
        subcatch_numexp.chms[LZ].prim_conc[kspc] = subcatch.chms[LZ].prim_conc[kspc];
        subcatch_numexp.chms[STREAM].prim_conc[kspc] = subcatch.chms[STREAM].prim_conc[kspc];

        subcatch_numexp.chms[SNOW].prim_actv[kspc] = subcatch.chms[SNOW].prim_actv[kspc];
        subcatch_numexp.chms[SURFACE].prim_actv[kspc] = subcatch.chms[SURFACE].prim_actv[kspc];
        subcatch_numexp.chms[UZ].prim_actv[kspc] = subcatch.chms[UZ].prim_actv[kspc];
        subcatch_numexp.chms[LZ].prim_actv[kspc] = subcatch.chms[LZ].prim_actv[kspc];
        subcatch_numexp.chms[STREAM].prim_actv[kspc] = subcatch.chms[STREAM].prim_actv[kspc];

        subcatch_numexp.chms[SNOW].tot_conc[kspc] = subcatch.chms[SNOW].tot_conc[kspc];
        subcatch_numexp.chms[SURFACE].tot_conc[kspc] = subcatch.chms[SURFACE].tot_conc[kspc];
        subcatch_numexp.chms[UZ].tot_conc[kspc] = subcatch.chms[UZ].tot_conc[kspc];
        subcatch_numexp.chms[LZ].tot_conc[kspc] = subcatch.chms[LZ].tot_conc[kspc];
        subcatch_numexp.chms[STREAM].tot_conc[kspc] = subcatch.chms[STREAM].tot_conc[kspc];

        subcatch_numexp.chms[SNOW].tot_mol[kspc] = subcatch.chms[SNOW].tot_mol[kspc];
        subcatch_numexp.chms[SURFACE].tot_mol[kspc] = subcatch.chms[SURFACE].tot_mol[kspc];
        subcatch_numexp.chms[UZ].tot_mol[kspc] = subcatch.chms[UZ].tot_mol[kspc];
        subcatch_numexp.chms[LZ].tot_mol[kspc] = subcatch.chms[LZ].tot_mol[kspc];
        subcatch_numexp.chms[STREAM].tot_mol[kspc] = subcatch.chms[STREAM].tot_mol[kspc];
    }

    for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)
    {
        subcatch_numexp.chms[SNOW].sec_conc[kspc] = subcatch.chms[SNOW].sec_conc[kspc];
        subcatch_numexp.chms[SURFACE].sec_conc[kspc] = subcatch.chms[SURFACE].sec_conc[kspc];
        subcatch_numexp.chms[UZ].sec_conc[kspc] = subcatch.chms[UZ].sec_conc[kspc];
        subcatch_numexp.chms[LZ].sec_conc[kspc] = subcatch.chms[LZ].sec_conc[kspc];
    }

}