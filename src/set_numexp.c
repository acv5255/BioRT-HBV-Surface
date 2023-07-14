#include "biort.h"

void CopyConstSubcatchProp(const Subcatchment* subcatch, Subcatchment* subcatch_numexp)
{
    subcatch_numexp->k1 = subcatch->k1;
    subcatch_numexp->k2 = subcatch->k2;
    subcatch_numexp->maxbas = subcatch->maxbas;
    subcatch_numexp->perc = subcatch->perc;
    subcatch_numexp->tt = subcatch->tt;
    subcatch_numexp->sfcf = subcatch->sfcf;
    
    subcatch_numexp->soil_surface = subcatch->soil_surface;
    subcatch_numexp->soil_sz = subcatch->soil_sz;
    subcatch_numexp->soil_dz = subcatch->soil_dz;
}

void CopyInitChemSubcatch(ReactionNetwork *rttbl, const Subcatchment subcatch[], Subcatchment subcatch_numexp[])
{
    for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
    {
        subcatch_numexp->chms[SNOW].ssa[kspc] = subcatch->chms[SNOW].ssa[kspc];
        subcatch_numexp->chms[SURFACE].ssa[kspc] = subcatch->chms[SURFACE].ssa[kspc];
        subcatch_numexp->chms[UZ].ssa[kspc] = subcatch->chms[UZ].ssa[kspc];
        subcatch_numexp->chms[LZ].ssa[kspc] = subcatch->chms[LZ].ssa[kspc];
        
        subcatch_numexp->chms[SURFACE].q10[kspc] = subcatch->chms[SURFACE].q10[kspc];
        subcatch_numexp->chms[UZ].q10[kspc] = subcatch->chms[UZ].q10[kspc];
        subcatch_numexp->chms[LZ].q10[kspc] = subcatch->chms[LZ].q10[kspc];
        
        subcatch_numexp->chms[SURFACE].sw_thld[kspc] = subcatch->chms[SURFACE].sw_thld[kspc];
        subcatch_numexp->chms[UZ].sw_thld[kspc] = subcatch->chms[UZ].sw_thld[kspc];
        subcatch_numexp->chms[LZ].sw_thld[kspc] = subcatch->chms[LZ].sw_thld[kspc];
        
        subcatch_numexp->chms[SURFACE].sw_exp[kspc] = subcatch->chms[SURFACE].sw_exp[kspc];
        subcatch_numexp->chms[UZ].sw_exp[kspc] = subcatch->chms[UZ].sw_exp[kspc];
        subcatch_numexp->chms[LZ].sw_exp[kspc] = subcatch->chms[LZ].sw_exp[kspc];
        
        subcatch_numexp->chms[SURFACE].n_alpha[kspc] = subcatch->chms[SURFACE].n_alpha[kspc];
        subcatch_numexp->chms[UZ].n_alpha[kspc] = subcatch->chms[UZ].n_alpha[kspc];
        subcatch_numexp->chms[LZ].n_alpha[kspc] = subcatch->chms[LZ].n_alpha[kspc];

        subcatch_numexp->chms[SNOW].prim_conc[kspc] = subcatch->chms[SNOW].prim_conc[kspc];
        subcatch_numexp->chms[SURFACE].prim_conc[kspc] = subcatch->chms[SURFACE].prim_conc[kspc];
        subcatch_numexp->chms[UZ].prim_conc[kspc] = subcatch->chms[UZ].prim_conc[kspc];
        subcatch_numexp->chms[LZ].prim_conc[kspc] = subcatch->chms[LZ].prim_conc[kspc];
        subcatch_numexp->chms[STREAM].prim_conc[kspc] = subcatch->chms[STREAM].prim_conc[kspc];

        subcatch_numexp->chms[SNOW].prim_actv[kspc] = subcatch->chms[SNOW].prim_actv[kspc];
        subcatch_numexp->chms[SURFACE].prim_actv[kspc] = subcatch->chms[SURFACE].prim_actv[kspc];
        subcatch_numexp->chms[UZ].prim_actv[kspc] = subcatch->chms[UZ].prim_actv[kspc];
        subcatch_numexp->chms[LZ].prim_actv[kspc] = subcatch->chms[LZ].prim_actv[kspc];
        subcatch_numexp->chms[STREAM].prim_actv[kspc] = subcatch->chms[STREAM].prim_actv[kspc];

        subcatch_numexp->chms[SNOW].tot_conc[kspc] = subcatch->chms[SNOW].tot_conc[kspc];
        subcatch_numexp->chms[SURFACE].tot_conc[kspc] = subcatch->chms[SURFACE].tot_conc[kspc];
        subcatch_numexp->chms[UZ].tot_conc[kspc] = subcatch->chms[UZ].tot_conc[kspc];
        subcatch_numexp->chms[LZ].tot_conc[kspc] = subcatch->chms[LZ].tot_conc[kspc];
        subcatch_numexp->chms[STREAM].tot_conc[kspc] = subcatch->chms[STREAM].tot_conc[kspc];

        subcatch_numexp->chms[SNOW].tot_mol[kspc] = subcatch->chms[SNOW].tot_mol[kspc];
        subcatch_numexp->chms[SURFACE].tot_mol[kspc] = subcatch->chms[SURFACE].tot_mol[kspc];
        subcatch_numexp->chms[UZ].tot_mol[kspc] = subcatch->chms[UZ].tot_mol[kspc];
        subcatch_numexp->chms[LZ].tot_mol[kspc] = subcatch->chms[LZ].tot_mol[kspc];
        subcatch_numexp->chms[STREAM].tot_mol[kspc] = subcatch->chms[STREAM].tot_mol[kspc];
    }

    for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)
    {
        subcatch_numexp->chms[SNOW].sec_conc[kspc] = subcatch->chms[SNOW].sec_conc[kspc];
        subcatch_numexp->chms[SURFACE].sec_conc[kspc] = subcatch->chms[SURFACE].sec_conc[kspc];
        subcatch_numexp->chms[UZ].sec_conc[kspc] = subcatch->chms[UZ].sec_conc[kspc];
        subcatch_numexp->chms[LZ].sec_conc[kspc] = subcatch->chms[LZ].sec_conc[kspc];
    }

}