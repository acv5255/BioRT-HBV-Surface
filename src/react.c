#include "biort.h"

const double SATN_MINIMUM = 1e-4;       // Minimum saturation to initiate kinetic reactions
const int MAX_ITERATIONS = 25;
const double MINIMUM_SUBSTEP = 1e-20;

void GetIAP(double iap[MAXSPS], const double activity[MAXSPS], const double dep_kin[MAXSPS][MAXSPS], int num_reactions, int num_stc) {
    for (int i = 0; i < num_reactions; i++) {
        iap[i] = 0.0;
        for (int j = 0; j < num_stc; j++) {
            iap[i] += log10(activity[j]) * dep_kin[i][j];
        }
        iap[i] = pow(10.0, iap[i]);
    }
    return;
}

double GetMonodTerm(const KineticTableEntry* entry, const ChemicalState* chms) {
    double monodterm = 1.0;

    for (int kmonod = 0; kmonod < entry->nmonod; kmonod++) {
        const int id = entry->monod_index[kmonod];
        const double k_i = entry->monod_para[kmonod];
        const double c_i = chms->prim_conc[id];
        monodterm *= c_i / (c_i + k_i);
    }

    return monodterm;
}

double GetInhibTerm(const KineticTableEntry* entry, const ChemicalState* chms) {
    double inhibterm = 1.0;

    for (int kinhib = 0; kinhib < entry->ninhib; kinhib++) {
        const int id = entry->inhib_index[kinhib];
        const double k_i = entry->inhib_para[kinhib];
        const double c_i = chms->prim_conc[id];
        inhibterm *= k_i / (c_i + k_i);
    }

    return inhibterm;
}

double GetDependenceTerm(const KineticTableEntry* entry, const ChemicalState* chms) {
    double dep_term = 1.0;
    for (int k = 0; k < entry->ndep; k++) {
        const int dep_ind = entry->dep_index[k];
        dep_term *= pow(chms->prim_actv[dep_ind], entry->dep_power[k]);
    }
    return dep_term;
}

void GetSecondarySpecies(double conc[MAXSPS], const double gamma[MAXSPS], const ReactionNetwork* rttbl) {
    for (int i = 0; i < rttbl->num_ssc; i++) {
        double tmpval = 0.0;
        for (int j = 0; j < rttbl->num_sdc; j++) {
            tmpval += (conc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
        }
        tmpval -= rttbl->keq[i] + gamma[i + rttbl->num_stc];
        if (!isfinite(tmpval)) {
            printf("Value in 'GetSecondarySpecies is infinite, exiting...\n"); 
            printf("Gamma for secondary species: %g\n", gamma[i + rttbl->num_stc]);
            printf("Values leading up to this:\n");
            for (int j = 0; j < rttbl->num_sdc; j++) {
                printf("(conc, gamma) = (%g, %g)\n", conc[j], gamma[j]);
            }
            exit(-1);
        }
        conc[i + rttbl->num_stc] = tmpval;
    }

    return;
}

void ReactZone(const int kzone, const double temp, const SoilConstants soil, const double tot_water, const ChemTableEntry chemtbl[],
        const KineticTableEntry kintbl[], const ReactionNetwork* rttbl, double stepsize, Subcatchment* subcatch) {
    
    const double Zw = soil.depth - tot_water / soil.porosity;     // UZ or LZ or SURFACE
    double substep;

    for (int kspc = 0; kspc < MAXSPS; kspc++)
    {
        subcatch->react_rate[kzone][kspc] = BADVAL;   // Set reaction rate to -999
    }

    double satn = tot_water / (soil.depth * soil.porosity);  // add porosity for saturation calculation
    satn = MIN(satn, 1.0);

    CheckChmsForNonFinite(&subcatch->chms[kzone], "react.c", 80);
    if (!CheckArrayForNan(subcatch->react_rate[kzone])) {
        printf("ChemicalState->react_rate contains nan in 'react.c' near line 82 with kzone = %d", kzone);
    }

    if (satn > SATN_MINIMUM)
    {
        substep = ReactControl(chemtbl, kintbl, rttbl, stepsize, soil.porosity, soil.depth, satn, temp, Zw,
            subcatch->react_rate[kzone], &subcatch->chms[kzone]);
        CheckChmsForNonFinite(&subcatch->chms[kzone], "react.c", 89);
        
        if (substep < 0.0) {
            printf("React failed to converge...\n");
        }
    }
}

void GetRates(double rate[MAXSPS], double rate_spe[MAXSPS], const double area[MAXSPS], const double ftemp[MAXSPS], const double fsw[MAXSPS], const double fzw[MAXSPS],
    const ReactionNetwork* rttbl, const KineticTableEntry kintbl[MAXSPS], const ChemicalState* chms) {
        
    for (int i = 0; i < rttbl->num_mkr; i++) {
        int min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

        if (kintbl[i].type == TST)
        {
            double iap[MAXSPS];
            double dependency[MAXSPS];
            GetIAP(iap, chms->prim_actv, rttbl->dep_kin, rttbl->num_mkr, rttbl->num_stc);
            double temp_keq = pow(10, rttbl->keq_kin[i]);
            dependency[i] = GetDependenceTerm(&kintbl[i], chms);

            // Calculate the predicted rate depending on the type of rate law
            //   rate_pre: rate per reaction (mol L-1 porous media s-1)
            //   area: m2 of min/ L of porous media
            //   rate: mol/m2 of min/s
            //   dependency: dimensionless
            rate[i] = area[min_pos] * pow(10, kintbl[i].rate) * dependency[i] * (1.0 - iap[i] / temp_keq) * ftemp[min_pos] * fsw[min_pos] * fzw[min_pos];
        }
        else if (kintbl[i].type == MONOD)
        {
            double monodterm = GetMonodTerm(&kintbl[i], chms);
            double inhibterm = GetInhibTerm(&kintbl[i], chms);

            // Based on CrunchTope
            rate[i] = area[min_pos] * pow(10, kintbl[i].rate) * monodterm * inhibterm * ftemp[min_pos] * fsw[min_pos] * fzw[min_pos];
        }
        for (int j = 0; j < rttbl->num_stc; j++)
        {
            rate_spe[j] += rate[i] * rttbl->dep_kin[i][j];
        }
    }

}

void Reaction(int kstep, double stepsize, const ChemTableEntry chemtbl[],
    const KineticTableEntry kintbl[], const ReactionNetwork *rttbl, Subcatchment* subcatch)
{
    const double temp = subcatch->tmp[kstep];
    if (subcatch->q[kstep][Q0] >= 0.1) {
        // printf("Starting surface reactions with a water value of %g mm\n", subcatch->q[kstep][Q0]);
        // printf("Surface zone chemical state before reactions: \n");
        // PrintChemicalState(&subcatch->chms[SURFACE]);
        ErrOnZeroRanged("react.c", "subcatch->chms[SURFACE].sec_conc", 150, subcatch->chms[SURFACE].sec_conc, rttbl->num_ssc);
        ReactZone(SURFACE, temp, subcatch->soil_surface, subcatch->q[kstep][Q0], chemtbl, kintbl, rttbl, stepsize, subcatch);
    }
    CheckChmsForNonFinite(&subcatch->chms[SURFACE], "react.c", 153);

    ReactZone(UZ, temp, subcatch->soil_sz, subcatch->ws[kstep][UZ], chemtbl, kintbl, rttbl, stepsize, subcatch);
    ReactZone(LZ, temp, subcatch->soil_dz, subcatch->ws[kstep][LZ], chemtbl, kintbl, rttbl, stepsize, subcatch);
}

int SolveReact(double stepsize, const ChemTableEntry chemtbl[], const KineticTableEntry kintbl[], const ReactionNetwork *rttbl,
    double satn, double temp, double porosity, double Zw, ChemicalState *chms)
{
    double          tmpconc[MAXSPS];
    double          tot_conc[MAXSPS];
    double          area[MAXSPS];
    double          ftemp[MAXSPS];
    double          fsw[MAXSPS];
    double          fzw[MAXSPS];
    double          gamma[MAXSPS];
    double          rate_pre[MAXSPS];
    double          rate_spe[MAXSPS];
    double          rate_spet[MAXSPS];
    const double    TMPPRB = 1.0E-4;
    const double    TMPPRB_INV = 1.0 / TMPPRB;
    const double    inv_sat = 1.0 / satn;//inv_sat(L of porous space/L of water)= 1/sat(L of water/L of porous space)
    SetZero(tmpconc); SetZero(tot_conc);
    CheckChmsForNonFinite(chms, "react.c", 174);
    ErrOnZeroRanged("react.c", "chms->sec_conc", 175, chms->sec_conc, rttbl->num_ssc);
    {
        int start = rttbl->num_stc - rttbl->num_min;
        int end = rttbl->num_stc;
        int offset = start;
        GetSurfaceAreaRange(area, chms->prim_conc, chms->ssa, chemtbl, rttbl->num_stc - rttbl->num_min, rttbl->num_stc, rttbl->num_stc - rttbl->num_min);
        GetTempFactorRange(ftemp, chms->q10, temp, start, end, offset);
        GetWTDepthFactorRange(fzw, Zw, chms->n_alpha, start, end, offset);
    }

    SoilMoistFactorRange(fsw, satn, chms->sw_thld, chms->sw_exp, rttbl->num_stc - rttbl->num_min, rttbl->num_stc, rttbl->num_stc - rttbl->num_min);
    SetZeroRange(rate_spe, 0, rttbl->num_stc);
    GetRates(rate_pre, rate_spe, area, ftemp, fsw, fzw, rttbl, kintbl, chms);

    for (int i = 0; i < rttbl->num_mkr; i++)
    {
        int min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

        if (rate_pre[i] < 0.0)
        {
            // Mineral cutoff when mineral is disappearing
            area[min_pos] = (chms->prim_conc[kintbl[i].position] < 1.0E-8) ? 0.0 : area[min_pos];
        }
    }

    for (int i = 0; i < rttbl->num_spc; i++)
    {
        // Aqueous species, saturation term for aqueous volume
        // rate(mol/m2 water/s)= rate(mol/m2 pm/s)*inv_sat(L of porous space/L of water)/porosity(L of porous space/L of pm)
        rate_spe[i] *= (chemtbl[i].itype == AQUEOUS) ? (inv_sat/porosity) : 1.0;
    }

    const double adh = rttbl->adh;
    const double bdh = rttbl->bdh;
    const double bdt = rttbl->bdt;

    double tot_cec = 0.0;
    double ionic_strength = 0.0;
    for (int i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        ErrorOnArrayNan(tmpconc, "tmpconc", "react.c", 215);
        double conc_val = (i < rttbl->num_stc) ? log10(chms->prim_conc[i]) : log10(chms->sec_conc[i - rttbl->num_stc]);
        if (!isfinite(conc_val)) {
            printf("Got non-finite value for log10(conc) for species '%s' with index %d\n", chemtbl[i].name, i);
        }
        tmpconc[i] = conc_val;

        tot_cec += (chemtbl[i].itype == CATION_ECHG) ? pow(10.0, tmpconc[i]) : 0.0;

        ionic_strength += 0.5 * pow(10, tmpconc[i]) * chemtbl[i].charge * chemtbl[i].charge;
        ErrorOnArrayNan(tmpconc, "tmpconc", "react.c", 221);
    }
    double root_ionic_strength = sqrt(ionic_strength);
    
    for (int i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        switch (chemtbl[i].itype)
        {
            case AQUEOUS:
                gamma[i] = (-adh * chemtbl[i].charge * chemtbl[i].charge * root_ionic_strength) /
                    (1.0 + bdh * chemtbl[i].size_fac * root_ionic_strength) + bdt * ionic_strength;
                break;
            case ADSORPTION:
                gamma[i] = log10(satn);
                break;
            case CATION_ECHG:
                gamma[i] = -log10(tot_cec);
                break;
            case MINERAL:
                gamma[i] = -tmpconc[i];
                break;
        }
    }
    ErrorOnArrayNan(tmpconc, "tmpconc", "react.c", 244);
    ErrorOnArrayNan(gamma, "gamma", "react.c", 245);
    CheckChmsForNonFinite(chms, "react.c", 246);

    int control = 0;
    double max_error = 0.0;
    do
    {
        const int matrix_dimension = rttbl->num_stc - rttbl->num_min;
        realtype** jcb = newDenseMat(matrix_dimension, matrix_dimension);
        double          residue[MAXSPS];
        double          residue_t[MAXSPS];
        sunindextype    p[MAXSPS];
        realtype        x[MAXSPS];

        ErrorOnArrayNanIter(tmpconc, "tmpconc", "react.c", 249, control);
        ErrorOnArrayNanIter(gamma, "gamma", "react.c", 250, control);
        GetSecondarySpecies(tmpconc, gamma, rttbl);
        ErrorOnArrayNanIter(tmpconc, "tmpconc", "react.c", 252, control);

        SetZeroRange(rate_spet, 0, rttbl->num_stc);
        ErrorOnArrayNanIter(tmpconc, "tmpconc", "react.c", 255, control);
        for (int i = 0; i < rttbl->num_mkr; i++)
        {
            double iap[MAXSPS];
            double dependency[MAXSPS];
            int min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

            if (kintbl[i].type == TST)
            {
                iap[i] = 0.0;
                for (int j = 0; j < rttbl->num_stc; j++)
                {
                    iap[i] += (chemtbl[j].itype != MINERAL) ? (tmpconc[j] + gamma[j]) * rttbl->dep_kin[i][j] : 0.0;
                }
                iap[i] = pow(10, iap[i]);

                double temp_keq = pow(10, rttbl->keq_kin[i]);

                dependency[i] = 0.0;
                for (int k = 0; k < kintbl[i].ndep; k++)
                {
                    dependency[i] += (tmpconc[kintbl[i].dep_index[k]] +
                        gamma[kintbl[i].dep_index[k]]) * kintbl[i].dep_power[k];
                }
                dependency[i] = pow(10.0, dependency[i]);

                // Calculate predicted rate depending on type of rate law
                // rate_pre: in mol / L water / s
                // area: m2/L water
                // rate: mol/m2/s
                // dependency: dimensionless
                rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * dependency[i] * (1.0 - (iap[i] / temp_keq)) * ftemp[min_pos] * fsw[min_pos] * fzw[min_pos];
            }
            else if (kintbl[i].type == MONOD)
            {
                const double monodterm = GetMonodTerm(&kintbl[i], chms);
                const double inhibterm = GetInhibTerm(&kintbl[i], chms);

                // Based on CrunchTope
                rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * monodterm * inhibterm * ftemp[min_pos] * fsw[min_pos] * fzw[min_pos];
            }

            for (int j = 0; j < rttbl->num_stc; j++)
            {
                rate_spet[j] += rate_pre[i] * rttbl->dep_kin[i][j];
            }
            // Adjust the unit of the calculated rate. Note that for mineral, the unit of rate and the unit of
            // concentration are mol/L porous media. For the aqueous species, the unit of the rate and the unit of the
            // concentration are mol/L pm and mol/L water respectively.
        }
        ErrorOnArrayNanIter(tmpconc, "tmpconc", "react.c", 301, control);
        for (int i = 0; i < rttbl->num_spc; i++)
        {
            // rate(mol/m2 water/s)= rate(mol/m2 pm/s)*inv_sat(L of porous space/L of water)/porosity(L of porous space/L of pm)
            rate_spet[i] *= (chemtbl[i].itype == AQUEOUS) ? (inv_sat/porosity) : 1.0;
        }

        // Calculate the residual for each aqueous primary species
        ErrorOnArrayNanIter(tmpconc, "tmpconc", "react.c", 309, control);
        for (int i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            double tmpval = 0.0;
            for (int j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
            {
                tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
            }
            tot_conc[i] = tmpval;
            residue[i] = tmpval - (chms->tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
        }
        ErrorOnArrayNan(tmpconc, "tmpconc", "react.c", 320);

        // Calculate the concentrations with a perturbation TMPPRB to calculate the numerical Jacobian
        for (int k = 0; k < rttbl->num_stc - rttbl->num_min; k++)
        {
            tmpconc[k] += TMPPRB;
            for (int i = 0; i < rttbl->num_ssc; i++)
            {
                double tmpval = 0.0;
                for (int j = 0; j < rttbl->num_sdc; j++)
                {
                    tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
                }
                tmpval -= rttbl->keq[i] + gamma[i + rttbl->num_stc];
                tmpconc[i + rttbl->num_stc] = tmpval;
            }
            for (int i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
            {
                double tmpval = 0.0;
                for (int j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
                {
                    tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
                }
                residue_t[i] = tmpval - (chms->tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
                jcb[k][i] = (residue_t[i] - residue[i]) * TMPPRB_INV;
            }
            tmpconc[k] -= TMPPRB;
        }
        
        for (int i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            x[i] = -residue[i];
        }

        int pivot_flg = denseGETRF(jcb, rttbl->num_stc - rttbl->num_min, rttbl->num_stc - rttbl->num_min, p);
        if (pivot_flg != 0)
        {
            destroyMat(jcb);
            return 1;
        }

        denseGETRS(jcb, rttbl->num_stc - rttbl->num_min, p, x);

        max_error = 0.0;
        for (int i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            // Take the step, limited to 0.3
            if (fabs(x[i]) < 0.3)
            {
                tmpconc[i] += x[i];
            }
            else
            {
                tmpconc[i] += (x[i] < 0) ? -0.3 : 0.3;
            }

            max_error = MAX(fabs(residue[i] / tot_conc[i]), max_error);
        }

        control++;
        if (control > MAX_ITERATIONS)   // Limit the model steps
        {
            // biort_printf(VL_NORMAL, "React failed to converge...\n");
            // printf("Failed to converge with stepsize of %g\n", stepsize);
            // printf("Jacobian matrix: \n");
            // PrintMatrix((const realtype**)jcb, matrix_dimension, matrix_dimension);
            destroyMat(jcb);
            // printf("\n\n");
            return 1;
            // exit(-1);
        }
        destroyMat(jcb);
    } while (max_error > TOLERANCE);

    ErrorOnArrayNan(tmpconc, "tmpconc", "react.c", 385);
    ErrorOnArrayNan(gamma, "gamma", "react.c", 386);
    ErrorOnArrayNan(tot_conc, "tot_conc", "react.c", 388);
    CheckChmsForNonFinite(chms, "react.c", 389);

    // Update secondary species here
    for (int i = 0; i < rttbl->num_ssc; i++)
    {
        double tmpval = 0.0;
        for (int j = 0; j < rttbl->num_sdc; j++)
        {
            tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
        }
        tmpval -= rttbl->keq[i] + gamma[i + rttbl->num_stc];
        tmpconc[i + rttbl->num_stc] = tmpval;
    }
    CheckChmsForNonFinite(chms, "react.c", 402);
    
    for (int i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
    {
        double tmpval = 0.0;
        for (int j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
        {
            tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
        }
        tot_conc[i] = tmpval;
    }
    CheckChmsForNonFinite(chms, "react.c", 413);

    for (int i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        if (i < rttbl->num_stc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                chms->tot_conc[i] += (rate_spe[i] + rate_spet[i]) * stepsize * 0.5;
                chms->prim_actv[i] = 1.0;
                chms->prim_conc[i] = chms->tot_conc[i];
            }
            else
            {
                chms->prim_conc[i] = pow(10, tmpconc[i]);
                chms->prim_actv[i] = pow(10, tmpconc[i] + gamma[i]);
                chms->tot_conc[i] = tot_conc[i];
            }
        }
        else
        {
            chms->sec_conc[i - rttbl->num_stc] = pow(10, tmpconc[i]);
        }
    }

    CheckChmsForNonFinite(chms, "react.c", 438);
    return 0;
}

double ReactControl(const ChemTableEntry chemtbl[], const KineticTableEntry kintbl[], const ReactionNetwork *rttbl,
    double stepsize, double porosity, double depth, double satn, double temp, double Zw, double react_rate[],
    ChemicalState *chms)
{
    double          substep;
    double          step_counter = 0.0;
    double          conc0[MAXSPS];

    // Copy initial mineral concentration to array
    for (int kspc = 0; kspc < rttbl->num_min; kspc++)
    {
        conc0[kspc] = chms->tot_conc[kspc + rttbl->num_stc - rttbl->num_min];
    }

    substep = stepsize;

    while (1.0 - step_counter / stepsize > 1.0E-10 && substep > MINIMUM_SUBSTEP)
    {
        CheckChmsForNonFinite(chms, "react.c", 451);
        int flag = SolveReact(substep, chemtbl, kintbl, rttbl, satn, temp, porosity, Zw, chms);
        CheckChmsForNonFinite(chms, "react.c", 453);

        if (flag == 0)
        {
            // Reaction passed with current step
            step_counter += substep;
        }
        else
        {
            substep *= 0.5;
        }
    }

    if (roundi(step_counter) != roundi(stepsize))   // Reactions fail
    {
        printf("Failed to converge with stepsize %g\n", substep);
        // return -substep;
        exit(-1);
    }
    else    // Reactions succeed
    {
        for (int kspc = 0; kspc < rttbl->num_min; kspc++)
        {
            // Calculate reaction rate (mole/m2-pm/day) = mol/L-pm/day * mm of depth   * (1m/1000mm) * (1000L/1m3)
            react_rate[kspc] =
                (chms->tot_conc[kspc + rttbl->num_stc - rttbl->num_min] - conc0[kspc]) * depth;
            CheckChmsForNonFinite(chms, "react.c", 477);
        }

        for (int kspc = 0; kspc <rttbl->num_spc; kspc++)
        {
            chms->tot_mol[kspc] = chms->tot_conc[kspc] * porosity * satn * depth; // tot_mol (moles-mm of water/L water ) =tot_conc (mol/L water) * satn * porosity * depth
            CheckChmsForNonFinite(chms, "react.c", 483);
        }
        if (roundi(substep) != roundi(stepsize))
        {
            return substep;
        }
        else
        {
            return 0.0;
        }
    }
}

double SoilTempFactor(double q10, double stc)
{
    return pow(q10, (stc - 20.0) / 10.0);
}

double SoilMoistFactor(double satn, double sw_thld, double sw_exp)
{
    double fsw;
    fsw     =   (satn <= sw_thld) ?
                    pow(satn/sw_thld, sw_exp) : pow((1.0 - satn) / (1.0 - sw_thld), sw_exp);
    return fsw;
}

double WTDepthFactor(double Zw, double n_alpha){
    double fzw;
    fzw     =   (n_alpha==0) ? 1 : exp(-fabs(n_alpha)*pow(Zw,n_alpha/fabs(n_alpha)));
    return    fzw;
}

void SoilMoistFactorRange(double dst[MAXSPS], double satn, const double sw_threshold[MAXSPS], const double sw_exponent[MAXSPS], int start, int end, int offset) {
    /* Calculate the soil moisture factor over an array */
    if (satn < 1.0) {
        for (int i = start; i < end; i++) {
            dst[i - offset] = SoilMoistFactor(satn, sw_threshold[i], sw_exponent[i]);
        }
    }
    return;
}

void GetSurfaceAreaRange(double area[MAXSPS], const double prim_conc[MAXSPS], const double ssa[MAXSPS], const ChemTableEntry chemtbl[], int start, int end, int offset) {
    for (int i = start; i < end; i++) {
        area[i - offset] = prim_conc[i] * ssa[i] * chemtbl[i].molar_mass;
    }
}

void GetTempFactorRange(double ftemp[MAXSPS], const double q10[MAXSPS], double temperature, int start, int end, int offset) {
    for (int i = start; i < end; i++) {
        ftemp[i - offset] = SoilTempFactor(q10[i], temperature);
    }
}

void GetWTDepthFactorRange(double fzw[MAXSPS], double Zw, const double n_alpha[MAXSPS], int start, int end, int offset) {
    for (int i = start; i < end; i++) {
        fzw[i - offset] = WTDepthFactor(Zw, n_alpha[i]);
    }
}