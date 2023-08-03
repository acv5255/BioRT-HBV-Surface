#include "biort.hpp"

const int MAX_SPECIATION_ITERS = 100;
const double MAX_STEP = 1.0;

void GetSecondarySpeciesRange(const ReactionNetwork* rttbl, array<f64, MAXSPS>& tmpconc, double gamma[MAXSPS]) {
    for (int i = 0; i < rttbl->num_ssc; i++)
    {
        tmpconc[i + rttbl->num_stc] = 0.0;
        for (int j = 0; j < rttbl->num_sdc; j++)
        {
            tmpconc[i + rttbl->num_stc] += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
        }
        tmpconc[i + rttbl->num_stc] -= rttbl->keq[i] + gamma[i + rttbl->num_stc];
    }
    return;
}

double GetIonicStrength(const ReactionNetwork* rttbl, const ChemTableEntry chemtbl[], const array<f64, MAXSPS>& conc) {
    double imat = 0.0;
    for (int i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        imat += pow(10, conc[i]) * chemtbl[i].charge * chemtbl[i].charge;
    }

    return 0.5 * imat;
}

int SolveSpeciation(const ChemTableEntry chemtbl[], const ControlData ctrl, const ReactionNetwork *rttbl, int speciation_flg,
    ChemicalState *chms)
{
    double          residue[MAXSPS];
    array<f64, MAXSPS> tmpconc;
    double          tot_conc[MAXSPS];
    double          log10gamma[MAXSPS];
    double          maxerror;
    const double    TMPPRB = 1E-2;

    // If speciation flg = 1, pH is defined. Total concentration is calculated from the activity of H+. Dependency is
    // the same but the total concentration for H+ does not need to be solved.
    // If speciation flg = 0, all defined value is total concentration
    SetZero(residue);
    tmpconc.fill(0.0);
    SetZero(tot_conc);
    SetZero(log10gamma);

    Log10Arr(chms->prim_conc, tmpconc, rttbl->num_stc);

    ComputeDependence(tmpconc, rttbl->dep_mtx, rttbl->keq, rttbl->num_ssc, rttbl->num_sdc, rttbl->num_stc);

    const int jcb_dim = (speciation_flg == 1) ? rttbl->num_stc - 1: rttbl->num_stc;

    int cur_iter = 0;
    do
    {
        sunindextype    pivots[MAXSPS];     // Whether any rows had to be swapped (???)
        realtype        x[MAXSPS];          // Solution to the problem
        int             row, col;
        realtype** jcb = newDenseMat(jcb_dim, jcb_dim);

        if (ctrl.use_activity == true)
        {
            double ionic_strength = GetIonicStrength(rttbl, chemtbl, tmpconc);
            double root_ion_str = sqrt(ionic_strength);

            for (int i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
            {
                // Aqueous species in the unit of mol/L, however the solids are in the unit of mol/L porous media.
                // Activity of solid is 1, log10 of activity is 0. By assigning gamma[minerals] to negative of the
                // tmpconc[minerals], we ensured the log10 of activity of solids are 0. gamma stores log10gamma[i]
                log10gamma[i] = (chemtbl[i].itype == MINERAL) ? -tmpconc[i] :
                    (-rttbl->adh * root_ion_str * chemtbl[i].charge * chemtbl[i].charge) /
                    (1.0 + rttbl->bdh * chemtbl[i].size_fac * root_ion_str) + rttbl->bdt * ionic_strength;

                if (speciation_flg == 1 && strcmp(chemtbl[i].name, "'H+'") == 0)
                {
                    tmpconc[i] = log10(chms->prim_actv[i]) - log10gamma[i];
                }
            }
        }

        GetLogActivity(tmpconc, log10gamma,rttbl->dep_mtx, rttbl->keq,  rttbl->num_ssc, rttbl->num_sdc, rttbl->num_stc);

        for (int i = 0; i < rttbl->num_stc; i++)
        {
            // Update the total concentration of H+ for later stage RT at initialization
            tot_conc[i] = 0.0;
            for (int j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
            {
                tot_conc[i] += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
            }

            if (speciation_flg == 1 && strcmp(chemtbl[i].name, "'H+'") == 0)
            {
                chms->tot_conc[i] = tot_conc[i];
            }

            residue[i] = tot_conc[i] - chms->tot_conc[i];
        }


        col = 0;
        for (int k = 0; k < rttbl->num_stc; k++)
        {
            if (speciation_flg == 1 && strcmp(chemtbl[k].name, "'H+'") == 0)
            {
                continue;
            }

            tmpconc[k] += TMPPRB;
            for (int i = 0; i < rttbl->num_ssc; i++)
            {
                tmpconc[i + rttbl->num_stc] = 0.0;
                for (int j = 0; j < rttbl->num_sdc; j++)
                {
                    tmpconc[i + rttbl->num_stc] += (tmpconc[j] + log10gamma[j]) * rttbl->dep_mtx[i][j];
                }
                tmpconc[i + rttbl->num_stc] -= rttbl->keq[i] + log10gamma[i + rttbl->num_stc];
            }

            row = 0;
            for (int i = 0; i < rttbl->num_stc; i++)
            {
                double          tmpval;

                if (speciation_flg == 1 && strcmp(chemtbl[i].name, "'H+'") == 0)
                {
                    continue;
                }

                tmpval = 0.0;
                for (int j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
                {
                    tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
                }

                jcb[col][row] = (tmpval - chms->tot_conc[i] - residue[i]) / TMPPRB;

                row++;
            }

            tmpconc[k] -= TMPPRB;

            col++;
        }

        row = 0;
        for (int i = 0; i < rttbl->num_stc; i++)
        {
            if (speciation_flg == 1 && strcmp(chemtbl[i].name, "'H+'") == 0)
            {
                continue;
            }

            x[row++] = -residue[i];
        }

        if (denseGETRF(jcb, jcb_dim, jcb_dim, pivots) != 0)
        {
            biort_printf(VL_ERROR, "Speciation error.\n");
            exit(EXIT_FAILURE);
        }

        denseGETRS(jcb, jcb_dim, pivots, x);

        maxerror = 0.0;
        row = 0;
        for (int i = 0; i < rttbl->num_stc; i++)
        {
            if (speciation_flg == 1 && strcmp(chemtbl[i].name, "'H+'") == 0)
            {
                continue;
            }

            tmpconc[i] += x[row] > 0 ? std::min(x[row], MAX_STEP) : std::max(x[row], -MAX_STEP);
            maxerror = std::max(fabs(residue[i] / tot_conc[i]), maxerror);
            row += 1;
        }

        destroyMat(jcb);

        // printf("Speciation errors: \n");
        // for (int i = 0; i < rttbl->num_stc; i++) {
        //     printf("%g ", fabs(residue[i] / tot_conc[i]));
        // }
        // printf("\n");

        // printf("Speciation error: %g\n", maxerror);
        cur_iter += 1;
    } while ((maxerror > TOLERANCE) && (cur_iter < MAX_SPECIATION_ITERS));
    // printf("\n");

    if (cur_iter == MAX_SPECIATION_ITERS) {
        printf("Speciation exceeded maximum iterations and failed to converge...\n");
        exit(-1);
    }

    GetSecondarySpeciesRange(rttbl, tmpconc, log10gamma);    // Compute secondary species

    for (int i = 0; i < rttbl->num_stc; i++) {
        tot_conc[i] = 0.0;
        for (int j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
        {
            tot_conc[i] += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
        }
        residue[i] = tot_conc[i] - chms->tot_conc[i];
    }

    // Set species concentrations now
    for (int i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        if (i < rttbl->num_stc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                chms->prim_conc[i] = pow(10, tmpconc[i]);
                chms->prim_actv[i] = 1.0;
            }
            else
            {
                chms->prim_conc[i] = pow(10, tmpconc[i]);
                chms->prim_actv[i] = pow(10, (tmpconc[i] + log10gamma[i]));
            }
        }
        else
        {
            chms->sec_conc[i - rttbl->num_stc] = pow(10, tmpconc[i]);
        }
    }
    CheckChmsForNonFinite(chms, "speciation.c", 217);

    return 0;
}

void Speciation(const ChemTableEntry chemtbl[], const ControlData ctrl, const ReactionNetwork *rttbl,
    Subcatchment* subcatch)
{
    SolveSpeciation(chemtbl, ctrl, rttbl, 0, &subcatch->chms[UZ]);
    SolveSpeciation(chemtbl, ctrl, rttbl, 0, &subcatch->chms[LZ]);
}

void StreamSpeciation(int step, const ChemTableEntry chemtbl[], const ControlData ctrl,
    const ReactionNetwork *rttbl, Subcatchment* subcatch)
{
    static int      init_flag = 1;

    // In the first model step, initialize primary concentrations, activities in river. This cannot be done during
    // initialization step because stream concentrations are unknown
    if (init_flag == 1)
    {

        for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            if (chemtbl[kspc].itype == AQUEOUS)
            {
                subcatch->chms[STREAM].prim_actv[kspc] = subcatch->chms[STREAM].tot_conc[kspc];
                subcatch->chms[STREAM].prim_conc[kspc] = subcatch->chms[STREAM].tot_conc[kspc];
            }
            else
            {
                subcatch->chms[STREAM].tot_conc[kspc] = ZERO_CONC;
                subcatch->chms[STREAM].prim_conc[kspc] = ZERO_CONC;
                subcatch->chms[STREAM].prim_actv[kspc] = ZERO_CONC;
                subcatch->chms[STREAM].tot_mol[kspc] = 0.0;
            }
        }

        for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            subcatch->chms[STREAM].sec_conc[kspc] = ZERO_CONC;
        }
    }

    init_flag = 0;

    if (subcatch->q[step][Q0] + subcatch->q[step][Q1] + subcatch->q[step][Q2] <= 0.0)
    {
        for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            if (chemtbl[kspc].itype == AQUEOUS)
            {
                subcatch->chms[STREAM].prim_actv[kspc] = subcatch->chms[STREAM].tot_conc[kspc];
                subcatch->chms[STREAM].prim_conc[kspc] = subcatch->chms[STREAM].tot_conc[kspc];
            }
            else
            {
                subcatch->chms[STREAM].tot_conc[kspc] = ZERO_CONC;
                subcatch->chms[STREAM].prim_conc[kspc] = ZERO_CONC;
                subcatch->chms[STREAM].prim_actv[kspc] = ZERO_CONC;
                subcatch->chms[STREAM].tot_mol[kspc] = 0.0;
            }
        }

        for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            subcatch->chms[STREAM].sec_conc[kspc] = ZERO_CONC;
        }
    }
    else
    {
        SolveSpeciation(chemtbl, ctrl, rttbl, 0, &subcatch->chms[STREAM]);
    }
}
