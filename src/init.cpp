#include "biort.hpp"

// Initialize RT structures
void InitChem(const char dir[], const CalibrationStruct *calib, const ControlData ctrl, array<ChemTableEntry, MAXSPS>& chemtbl,
    array<KineticTableEntry, MAXSPS>& kintbl, ReactionNetwork *rttbl, Subcatchment* subcatch)
{
    char            file_name[MAXSTRING];
    FILE           *file_pointer;

    // LOOK UP DATABASE TO FIND REQUIRED PARAMETERS AND DEPENDENCIES FOR CHEMICAL SPECIES IN CDBS.TXT
    sprintf(file_name, "input/%s/cdbs.txt", dir);
    file_pointer = fopen(file_name, "r");

    Lookup(file_pointer, calib, chemtbl, kintbl, rttbl);
    fclose(file_pointer);

    for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
    {
        // Apply calibration
        // subcatch->chms[SNOW].ssa[kspc] *= (chemtbl[kspc].itype == MINERAL) ? calib->ssa : 1.0;   // 2021-05-07
        subcatch->chms[UZ].soil_parameters.ssa[kspc] *= (chemtbl[kspc].itype == MINERAL) ? calib->ssa : 1.0;
        subcatch->chms[LZ].soil_parameters.ssa[kspc] *= (chemtbl[kspc].itype == MINERAL) ? calib->ssa : 1.0;

        // Snow and soil moisture zone should have the same concentrations as the upper zone at the beginning
        subcatch->chms[SNOW].tot_conc[kspc] = subcatch->chms[UZ].tot_conc[kspc];
        subcatch->chms[SNOW].prim_conc[kspc] = subcatch->chms[SNOW].tot_conc[kspc];
        subcatch->chms[SNOW].soil_parameters.ssa[kspc] = subcatch->chms[UZ].soil_parameters.ssa[kspc];
        subcatch->chms[SNOW].tot_mol[kspc] =
            subcatch->chms[SNOW].tot_conc[kspc] * subcatch->ws[0][SNOW];
    }

    // Initialize upper and lower zone concentrations taking into account speciation
    InitChemState(subcatch->soil_surface.porosity, subcatch->ws[0][SURFACE], chemtbl, rttbl, ctrl,
        &subcatch->chms[SURFACE]);
    CheckChmsForNonFinite(&subcatch->chms[SURFACE], "init.c", 35);
    InitChemState(subcatch->soil_sz.porosity, subcatch->ws[0][UZ], chemtbl, rttbl, ctrl,
        &subcatch->chms[UZ]);
    CheckChmsForNonFinite(&subcatch->chms[UZ], "init.c", 38);
    InitChemState(subcatch->soil_dz.porosity, subcatch->ws[0][LZ], chemtbl, rttbl, ctrl,
        &subcatch->chms[LZ]);
    CheckChmsForNonFinite(&subcatch->chms[LZ], "init.c", 41);
}

void InitChemState(double smcmax, double vol, const array<ChemTableEntry, MAXSPS>& chemtbl, const ReactionNetwork *rttbl,
    const ControlData ctrl, ChemicalState *chms)
{
    for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
    {
        if (strcmp(chemtbl[kspc].name, "'H+'") == 0)
        {
            chms->prim_actv[kspc] = chms->tot_conc[kspc];
            chms->prim_conc[kspc] = chms->tot_conc[kspc];
        }
        else if (chemtbl[kspc].itype == MINERAL)
        {
            // Update the concentration of mineral using molar volume
            //concentration (mole of mineral/L of porous media) = Absolute mineral volume fraction (cm3 of mineral/cm3 of porous media)*
          //1000 (cm3/l)/molar volume (cm3 of mineral/mole of mineral)

            chms->tot_conc[kspc] *= 1000.0 / chemtbl[kspc].molar_vol ;     // Absolute mineral volume fraction
            chms->prim_actv[kspc] = 1.0;
            chms->prim_conc[kspc] = chms->tot_conc[kspc];
        }
        else if ((chemtbl[kspc].itype == CATION_ECHG) || (chemtbl[kspc].itype == ADSORPTION))
        {
            chms->prim_actv[kspc] = chms->tot_conc[kspc] * 0.5;
            chms->tot_conc[kspc] *= (1.0 - smcmax) * 2650.0;    // Change unit of CEC (eq g-1) into C(ion site)
                                                                // (eq L-1 porous space), assuming density of solid is
                                                                // always 2650 g L-1
            chms->prim_conc[kspc] = chms->tot_conc[kspc];
        }
        else
        {
            chms->prim_actv[kspc] = chms->tot_conc[kspc] * 0.5;
            chms->prim_conc[kspc] = chms->tot_conc[kspc] * 0.5;
        }
    }

    for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)
    {
        chms->sec_conc[kspc] = ZERO_CONC;
    }

    // Speciation
    SolveSpeciation(chemtbl, ctrl, rttbl, 1, chms);

    // Total moles should be calculated after speciation
    for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
    {
        chms->tot_mol[kspc] = chms->tot_conc[kspc] * vol;
    }
}
