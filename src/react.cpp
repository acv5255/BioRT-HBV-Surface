#include "biort.hpp"

void GetIAP(double iap[MAXSPS], const array<f64, MAXSPS>& activity, const array<array<f64, MAXSPS>, MAXSPS>& dep_kin, int num_reactions, int num_stc) {
    for (int i = 0; i < num_reactions; i++) {
        iap[i] = 0.0;
        for (int j = 0; j < num_stc; j++) {
            iap[i] += log10(activity[j]) * dep_kin[i][j];
        }
        iap[i] = pow(10.0, iap[i]);
    }
    return;
}

double GetMonodTerm(const KineticTableEntry* entry, const ChemicalState& chms) {
    double monodterm = 1.0;

    for (int kmonod = 0; kmonod < entry->nmonod; kmonod++) {
        const int id = entry->monod_index[kmonod];
        const double k_i = entry->monod_para[kmonod];
        const double c_i = chms.prim_conc[id];
        monodterm *= c_i / (c_i + k_i);
    }

    return monodterm;
}

double GetInhibTerm(const KineticTableEntry* entry, const ChemicalState& chms) {
    double inhibterm = 1.0;

    for (int kinhib = 0; kinhib < entry->ninhib; kinhib++) {
        const int id = entry->inhib_index[kinhib];
        const double k_i = entry->inhib_para[kinhib];
        const double c_i = chms.prim_conc[id];
        inhibterm *= k_i / (c_i + k_i);
    }

    return inhibterm;
}

double GetDependenceTerm(const KineticTableEntry* entry, const ChemicalState& chms) {
    double dep_term = 1.0;
    for (int k = 0; k < entry->ndep; k++) {
        const int dep_ind = entry->dep_index[k];
        dep_term *= pow(chms.prim_actv[dep_ind], entry->dep_power[k]);
    }
    return dep_term;
}

void GetSecondarySpecies(array<f64, MAXSPS>& conc, const array<f64, MAXSPS>& gamma, const ReactionNetwork& rttbl) {
    for (int i = 0; i < rttbl.num_ssc; i++) {
        double tmpval = 0.0;
        for (int j = 0; j < rttbl.num_sdc; j++) {
            tmpval += (conc[j] + gamma[j]) * rttbl.dep_mtx[i][j];
        }
        tmpval -= rttbl.keq[i] + gamma[i + rttbl.num_stc];
        if (std::isinf(tmpval)) {
            printf("Value in 'GetSecondarySpecies is infinite, exiting...\n"); 
            printf("Gamma for secondary species: %g\n", gamma[i + rttbl.num_stc]);
            printf("Values leading up to this:\n");
            for (int j = 0; j < rttbl.num_sdc; j++) {
                printf("(conc, gamma) = (%g, %g)\n", conc[j], gamma[j]);
            }
            exit(-1);
        }
        conc[i + rttbl.num_stc] = tmpval;
    }

    return;
}

void ReactZone(const int kzone, const double temp, const SoilConstants soil, const double tot_water, const array<ChemTableEntry, MAXSPS>& chemtbl,
        const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl, double stepsize, Subcatchment& subcatch) {
    
    const double Zw = soil.depth - tot_water / soil.porosity;     // UZ or LZ or SURFACE
    double substep;

    for (int kspc = 0; kspc < MAXSPS; kspc++)
    {
        subcatch.react_rate[kzone][kspc] = BADVAL;   // Set reaction rate to -999
    }

    double satn = tot_water / (soil.depth * soil.porosity);  // add porosity for saturation calculation
    satn = std::min(satn, 1.0);

    CheckChmsForNonFinite(subcatch.chms[kzone], "react.c", 80);
    if (!CheckArrayForNan(subcatch.react_rate[kzone])) {
        printf("ChemicalState->react_rate contains nan in 'react.c' near line 82 with kzone = %d", kzone);
    }

    if (satn > SATN_MINIMUM)
    {
        substep = ReactControl(chemtbl, kintbl, rttbl, stepsize, soil.porosity, soil.depth, satn, temp, Zw,
            subcatch.react_rate[kzone], subcatch.chms[kzone]);
        CheckChmsForNonFinite(subcatch.chms[kzone], "react.c", 89);
        
        if (substep < 0.0) {
            printf("React failed to converge...\n");
        }
    }
}

void ReactSurfaceZone(const double temp, const SoilConstants soil, const double tot_water, const array<ChemTableEntry, MAXSPS>& chemtbl,
        const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl, double stepsize, Subcatchment& subcatch) {
    
    double substep;

    for (int kspc = 0; kspc < MAXSPS; kspc++)
    {
        subcatch.react_rate[SURFACE][kspc] = BADVAL;   // Set reaction rate to -999
    }

    CheckChmsForNonFinite(subcatch.chms[SURFACE], "react.c", 80);
    if (!CheckArrayForNan(subcatch.react_rate[SURFACE])) {
        printf("ChemicalState->react_rate contains nan in 'react.c' near line 82 with SURFACE = %d", SURFACE);
    }

    substep = ReactSurfaceControl(chemtbl, kintbl, rttbl, stepsize, soil.porosity, soil.depth, tot_water, temp,
        subcatch.react_rate[SURFACE], subcatch.chms[SURFACE]);
    CheckChmsForNonFinite(subcatch.chms[SURFACE], "react.c", 89);
    
    if (substep < 0.0) {
        printf("React failed to converge...\n");
    }
}

void GetRates(array<f64, MAXSPS>& rate, array<f64, MAXSPS>& rate_spe, const array<f64, MAXSPS>& area, const array<f64, MAXSPS>& ftemp, const array<f64, MAXSPS>& fsw, const array<f64, MAXSPS>& fzw,
    const ReactionNetwork& rttbl, const array<KineticTableEntry, MAXSPS>& kintbl, const ChemicalState& chms) {
        
    for (int i = 0; i < rttbl.num_mkr; i++) {
        int min_pos = kintbl[i].position - rttbl.num_stc + rttbl.num_min;

        if (kintbl[i].type == TST)
        {
            double iap[MAXSPS];
            double dependency[MAXSPS];
            GetIAP(iap, chms.prim_actv, rttbl.dep_kin, rttbl.num_mkr, rttbl.num_stc);
            double temp_keq = pow(10, rttbl.keq_kin[i]);
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
        for (int j = 0; j < rttbl.num_stc; j++)
        {
            rate_spe[j] += rate[i] * rttbl.dep_kin[i][j];
        }
    }

}

void Reaction(int kstep, double stepsize, const array<ChemTableEntry, MAXSPS>& chemtbl,
    const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl, Subcatchment& subcatch)
{
    const double temp = subcatch.tmp[kstep];
    if (subcatch.q[kstep][Q0] >= 0.1) {
        const double tot_water_surf = subcatch.ws[kstep][SURFACE] + subcatch.q[kstep][Q0];
        ReactSurfaceZone(temp, subcatch.soil_surface, tot_water_surf, chemtbl, kintbl, rttbl, stepsize, subcatch);
    }
    CheckChmsForNonFinite(subcatch.chms[SURFACE], "react.c", 153);

    ReactZone(UZ, temp, subcatch.soil_sz, subcatch.ws[kstep][UZ], chemtbl, kintbl, rttbl, stepsize, subcatch);
    ReactZone(LZ, temp, subcatch.soil_dz, subcatch.ws[kstep][LZ], chemtbl, kintbl, rttbl, stepsize, subcatch);
}

int SolveReact(double stepsize, const array<ChemTableEntry, MAXSPS>& chemtbl, const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl,
    double satn, double temp, double porosity, double Zw, ChemicalState&chms)
{
    array<f64, MAXSPS> tmpconc = { 0.0 };
    array<f64, MAXSPS> tot_conc = { 0.0 };
    array<f64, MAXSPS> area;
    array<f64, MAXSPS> ftemp;
    array<f64, MAXSPS> fsw;
    array<f64, MAXSPS> fzw;
    array<f64, MAXSPS> gamma;
    array<f64, MAXSPS>rate_pre;
    array<f64, MAXSPS> rate_spe;
    array<f64, MAXSPS> rate_spet;
    const double    TMPPRB = 1.0E-4;
    const double    TMPPRB_INV = 1.0 / TMPPRB;
    const double    inv_sat = 1.0 / satn;//inv_sat(L of porous space/L of water)= 1/sat(L of water/L of porous space)
    CheckChmsForNonFinite(chms, "react.c", 174);
    // ErrOnZeroRanged("react.c", "chms.sec_conc", 175, chms.sec_conc, rttbl.num_ssc);
    {
        int start = rttbl.num_stc - rttbl.num_min;
        int end = rttbl.num_stc;
        int offset = start;
        GetSurfaceAreaRange(area, chms.prim_conc, chms.soil_parameters.ssa, chemtbl, rttbl.num_stc - rttbl.num_min, rttbl.num_stc, rttbl.num_stc - rttbl.num_min);
        GetTempFactorRange(ftemp, chms.soil_parameters.q10, temp, start, end, offset);
        GetWTDepthFactorRange(fzw, Zw, chms.soil_parameters.n_alpha, start, end, offset);
    }

    SoilMoistFactorRange(fsw, satn, chms.soil_parameters.sw_thld, chms.soil_parameters.sw_exp, rttbl.num_stc - rttbl.num_min, rttbl.num_stc, rttbl.num_stc - rttbl.num_min);
    rate_spe.fill(0.0);
    GetRates(rate_pre, rate_spe, area, ftemp, fsw, fzw, rttbl, kintbl, chms);

    for (int i = 0; i < rttbl.num_mkr; i++)
    {
        int min_pos = kintbl[i].position - rttbl.num_stc + rttbl.num_min;

        if (rate_pre[i] < 0.0)
        {
            // Mineral cutoff when mineral is disappearing
            area[min_pos] = (chms.prim_conc[kintbl[i].position] < 1.0E-8) ? 0.0 : area[min_pos];
        }
    }

    for (int i = 0; i < rttbl.num_spc; i++)
    {
        // Aqueous species, saturation term for aqueous volume
        // rate(mol/m2 water/s)= rate(mol/m2 pm/s)*inv_sat(L of porous space/L of water)/porosity(L of porous space/L of pm)
        rate_spe[i] *= (chemtbl[i].itype == AQUEOUS) ? (inv_sat/porosity) : 1.0;
    }

    const double adh = rttbl.adh;
    const double bdh = rttbl.bdh;
    const double bdt = rttbl.bdt;

    double tot_cec = 0.0;
    double ionic_strength = 0.0;
    for (int i = 0; i < rttbl.num_stc + rttbl.num_ssc; i++)
    {
        double conc_val = (i < rttbl.num_stc) ? log10(chms.prim_conc[i]) : log10(chms.sec_conc[i - rttbl.num_stc]);
        if (std::isinf(conc_val)) {
            printf("Got non-finite value for log10(conc) for species '%s' with index %d\n", chemtbl[i].name.c_str(), i);
        }
        tmpconc[i] = conc_val;

        tot_cec += (chemtbl[i].itype == CATION_ECHG) ? pow(10.0, tmpconc[i]) : 0.0;

        ionic_strength += 0.5 * pow(10, tmpconc[i]) * chemtbl[i].charge * chemtbl[i].charge;
    }
    double root_ionic_strength = sqrt(ionic_strength);
    
    for (int i = 0; i < rttbl.num_stc + rttbl.num_ssc; i++)
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
    CheckChmsForNonFinite(chms, "react.c", 246);

    int control = 0;
    double max_error = 0.0;
    do
    {
        const int matrix_dimension = rttbl.num_stc - rttbl.num_min;
        realtype** jcb = newDenseMat(matrix_dimension, matrix_dimension);
        array<f64, MAXSPS> residue;
        array<f64, MAXSPS> residue_t;
        array<sunindextype, MAXSPS> p;
        array<realtype, MAXSPS> x;

        GetSecondarySpecies(tmpconc, gamma, rttbl);

        rate_spet.fill(0.0);
        for (int i = 0; i < rttbl.num_mkr; i++)
        {
            double iap[MAXSPS];
            double dependency[MAXSPS];
            int min_pos = kintbl[i].position - rttbl.num_stc + rttbl.num_min;

            if (kintbl[i].type == TST)
            {
                iap[i] = 0.0;
                for (int j = 0; j < rttbl.num_stc; j++)
                {
                    iap[i] += (chemtbl[j].itype != MINERAL) ? (tmpconc[j] + gamma[j]) * rttbl.dep_kin[i][j] : 0.0;
                }
                iap[i] = pow(10, iap[i]);

                double temp_keq = pow(10, rttbl.keq_kin[i]);

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

            for (int j = 0; j < rttbl.num_stc; j++)
            {
                rate_spet[j] += rate_pre[i] * rttbl.dep_kin[i][j];
            }
            // Adjust the unit of the calculated rate. Note that for mineral, the unit of rate and the unit of
            // concentration are mol/L porous media. For the aqueous species, the unit of the rate and the unit of the
            // concentration are mol/L pm and mol/L water respectively.
        }
        for (int i = 0; i < rttbl.num_spc; i++)
        {
            // rate(mol/m2 water/s)= rate(mol/m2 pm/s)*inv_sat(L of porous space/L of water)/porosity(L of porous space/L of pm)
            rate_spet[i] *= (chemtbl[i].itype == AQUEOUS) ? (inv_sat/porosity) : 1.0;
        }

        // Calculate the residual for each aqueous primary species
        for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
        {
            double tmpval = 0.0;
            for (int j = 0; j < rttbl.num_stc + rttbl.num_ssc; j++)
            {
                tmpval += rttbl.conc_contrib[i][j] * pow(10, tmpconc[j]);
            }
            tot_conc[i] = tmpval;
            residue[i] = tmpval - (chms.tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
        }

        // Calculate the concentrations with a perturbation TMPPRB to calculate the numerical Jacobian
        for (int k = 0; k < rttbl.num_stc - rttbl.num_min; k++)
        {
            tmpconc[k] += TMPPRB;
            for (int i = 0; i < rttbl.num_ssc; i++)
            {
                double tmpval = 0.0;
                for (int j = 0; j < rttbl.num_sdc; j++)
                {
                    tmpval += (tmpconc[j] + gamma[j]) * rttbl.dep_mtx[i][j];
                }
                tmpval -= rttbl.keq[i] + gamma[i + rttbl.num_stc];
                tmpconc[i + rttbl.num_stc] = tmpval;
            }
            for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
            {
                double tmpval = 0.0;
                for (int j = 0; j < rttbl.num_stc + rttbl.num_ssc; j++)
                {
                    tmpval += rttbl.conc_contrib[i][j] * pow(10, tmpconc[j]);
                }
                residue_t[i] = tmpval - (chms.tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
                jcb[k][i] = (residue_t[i] - residue[i]) * TMPPRB_INV;
            }
            tmpconc[k] -= TMPPRB;
        }
        
        for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
        {
            x[i] = -residue[i];
        }

        int pivot_flg = denseGETRF(jcb, rttbl.num_stc - rttbl.num_min, rttbl.num_stc - rttbl.num_min, p.data());
        if (pivot_flg != 0)
        {
            destroyMat(jcb);
            return 1;
        }

        denseGETRS(jcb, rttbl.num_stc - rttbl.num_min, p.data(), x.data());

        max_error = 0.0;
        for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
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

            max_error = std::max(fabs(residue[i] / tot_conc[i]), max_error);
        }

        control++;
        if (control > MAX_ITERATIONS)   // Limit the model steps
        {
            destroyMat(jcb);
            return 1;
        }
        destroyMat(jcb);
    } while (max_error > TOLERANCE);

    CheckChmsForNonFinite(chms, "react.c", 389);

    // Update secondary species here
    for (int i = 0; i < rttbl.num_ssc; i++)
    {
        double tmpval = 0.0;
        for (int j = 0; j < rttbl.num_sdc; j++)
        {
            tmpval += (tmpconc[j] + gamma[j]) * rttbl.dep_mtx[i][j];
        }
        tmpval -= rttbl.keq[i] + gamma[i + rttbl.num_stc];
        tmpconc[i + rttbl.num_stc] = tmpval;
    }
    CheckChmsForNonFinite(chms, "react.c", 402);
    
    for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
    {
        double tmpval = 0.0;
        for (int j = 0; j < rttbl.num_stc + rttbl.num_ssc; j++)
        {
            tmpval += rttbl.conc_contrib[i][j] * pow(10, tmpconc[j]);
        }
        tot_conc[i] = tmpval;
    }
    CheckChmsForNonFinite(chms, "react.c", 413);

    for (int i = 0; i < rttbl.num_stc + rttbl.num_ssc; i++)
    {
        if (i < rttbl.num_stc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                chms.tot_conc[i] += (rate_spe[i] + rate_spet[i]) * stepsize * 0.5;
                chms.prim_actv[i] = 1.0;
                chms.prim_conc[i] = chms.tot_conc[i];
            }
            else
            {
                chms.prim_conc[i] = pow(10, tmpconc[i]);
                chms.prim_actv[i] = pow(10, tmpconc[i] + gamma[i]);
                chms.tot_conc[i] = tot_conc[i];
            }
        }
        else
        {
            chms.sec_conc[i - rttbl.num_stc] = pow(10, tmpconc[i]);
        }
    }

    return 0;
}

int SolveSurfaceReact(double stepsize, const array<ChemTableEntry, MAXSPS>& chemtbl, const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl,
    double tot_water, double temp, double porosity, ChemicalState&chms)
{
    array<double, MAXSPS> tmpconc = { 0.0 };
    array<double, MAXSPS> tot_conc = { 0.0 };
    array<double, MAXSPS> area;
    array<double, MAXSPS> ftemp;
    array<double, MAXSPS> fsw;
    array<double, MAXSPS> fzw;
    array<f64, MAXSPS> gamma;
    array<double, MAXSPS> rate_pre;
    array<double, MAXSPS> rate_spe;
    array<double, MAXSPS> rate_spet;
    const double    TMPPRB = 1.0E-4;
    const double    TMPPRB_INV = 1.0 / TMPPRB;
    const double    inv_sat = 1.0 / tot_water;//inv_sat(L of porous space/L of water)= 1/sat(L of water/L of porous space)
    CheckChmsForNonFinite(chms, "react.c", 174);
    {
        int start = rttbl.num_stc - rttbl.num_min;
        int end = rttbl.num_stc;
        int offset = start;
        GetSurfaceAreaRange(area, chms.prim_conc, chms.soil_parameters.ssa, chemtbl, rttbl.num_stc - rttbl.num_min, rttbl.num_stc, rttbl.num_stc - rttbl.num_min);
        GetTempFactorRange(ftemp, chms.soil_parameters.q10, temp, start, end, offset);
        // GetWTDepthFactorRange(fzw, Zw, chms.n_alpha, start, end, offset);
    }

    // Set the factors for the surface zone
    for (int i = 0; i < rttbl.num_min; i++) {
        fzw[i] = 1.0;
    }

    for (int i = rttbl.num_stc - rttbl.num_min; i < rttbl.num_stc; i++) {
        int offset = rttbl.num_stc - rttbl.num_min;
        fsw[i - offset] = pow(tot_water, chms.soil_parameters.sw_exp[i]);
        // fsw[i - offset] = 1.0;
        // printf("fsw for species %d: %g\n", i - offset, fsw[i - offset]);
    }

    // SoilMoistFactorRange(fsw, tot_water, chms.sw_thld, chms.sw_exp, rttbl.num_stc - rttbl.num_min, rttbl.num_stc, rttbl.num_stc - rttbl.num_min);
    rate_spe.fill(0.0);
    GetRates(rate_pre, rate_spe, area, ftemp, fsw, fzw, rttbl, kintbl, chms);

    for (int i = 0; i < rttbl.num_mkr; i++)
    {
        int min_pos = kintbl[i].position - rttbl.num_stc + rttbl.num_min;

        if (rate_pre[i] < 0.0)
        {
            // Mineral cutoff when mineral is disappearing
            area[min_pos] = (chms.prim_conc[kintbl[i].position] < 1.0E-8) ? 0.0 : area[min_pos];
        }
    }

    for (int i = 0; i < rttbl.num_spc; i++)
    {
        // Aqueous species, saturation term for aqueous volume
        // rate(mol/m2 water/s)= rate(mol/m2 pm/s)*inv_sat(L of porous space/L of water)/porosity(L of porous space/L of pm)
        rate_spe[i] *= (chemtbl[i].itype == AQUEOUS) ? (inv_sat/porosity) : 1.0;
    }

    const double adh = rttbl.adh;
    const double bdh = rttbl.bdh;
    const double bdt = rttbl.bdt;

    double tot_cec = 0.0;
    double ionic_strength = 0.0;
    for (int i = 0; i < rttbl.num_stc + rttbl.num_ssc; i++)
    {
        double conc_val = (i < rttbl.num_stc) ? log10(chms.prim_conc[i]) : log10(chms.sec_conc[i - rttbl.num_stc]);
        if (std::isinf(conc_val)) {
            printf("Got non-finite value for log10(conc) for species '%s' with index %d\n", chemtbl[i].name.c_str(), i);
        }
        tmpconc[i] = conc_val;

        tot_cec += (chemtbl[i].itype == CATION_ECHG) ? pow(10.0, tmpconc[i]) : 0.0;

        ionic_strength += 0.5 * pow(10, tmpconc[i]) * chemtbl[i].charge * chemtbl[i].charge;
    }
    double root_ionic_strength = sqrt(ionic_strength);
    
    for (int i = 0; i < rttbl.num_stc + rttbl.num_ssc; i++)
    {
        switch (chemtbl[i].itype)
        {
            case AQUEOUS:
                gamma[i] = (-adh * chemtbl[i].charge * chemtbl[i].charge * root_ionic_strength) /
                    (1.0 + bdh * chemtbl[i].size_fac * root_ionic_strength) + bdt * ionic_strength;
                break;
            case ADSORPTION:
                gamma[i] = log10(tot_water);
                break;
            case CATION_ECHG:
                gamma[i] = -log10(tot_cec);
                break;
            case MINERAL:
                gamma[i] = -tmpconc[i];
                break;
        }
    }
    CheckChmsForNonFinite(chms, "react.c", 246);

    int control = 0;
    double max_error = 0.0;
    do
    {
        const int matrix_dimension = rttbl.num_stc - rttbl.num_min;
        realtype** jcb = newDenseMat(matrix_dimension, matrix_dimension);
        array<f64, MAXSPS> residue;
        array<f64, MAXSPS> residue_t;
        array<sunindextype, MAXSPS> p;
        array<realtype, MAXSPS> x;

        GetSecondarySpecies(tmpconc, gamma, rttbl);

        rate_spet.fill(0.0);
        for (int i = 0; i < rttbl.num_mkr; i++)
        {
            double iap[MAXSPS];
            double dependency[MAXSPS];
            int min_pos = kintbl[i].position - rttbl.num_stc + rttbl.num_min;

            if (kintbl[i].type == TST)
            {
                iap[i] = 0.0;
                for (int j = 0; j < rttbl.num_stc; j++)
                {
                    iap[i] += (chemtbl[j].itype != MINERAL) ? (tmpconc[j] + gamma[j]) * rttbl.dep_kin[i][j] : 0.0;
                }
                iap[i] = pow(10, iap[i]);

                double temp_keq = pow(10, rttbl.keq_kin[i]);

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

            for (int j = 0; j < rttbl.num_stc; j++)
            {
                rate_spet[j] += rate_pre[i] * rttbl.dep_kin[i][j];
            }
            // Adjust the unit of the calculated rate. Note that for mineral, the unit of rate and the unit of
            // concentration are mol/L porous media. For the aqueous species, the unit of the rate and the unit of the
            // concentration are mol/L pm and mol/L water respectively.
        }
        for (int i = 0; i < rttbl.num_spc; i++)
        {
            // rate(mol/m2 water/s)= rate(mol/m2 pm/s)*inv_sat(L of porous space/L of water)/porosity(L of porous space/L of pm)
            rate_spet[i] *= (chemtbl[i].itype == AQUEOUS) ? (inv_sat/porosity) : 1.0;
        }

        // Calculate the residual for each aqueous primary species
        for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
        {
            double tmpval = 0.0;
            for (int j = 0; j < rttbl.num_stc + rttbl.num_ssc; j++)
            {
                tmpval += rttbl.conc_contrib[i][j] * pow(10, tmpconc[j]);
            }
            tot_conc[i] = tmpval;
            residue[i] = tmpval - (chms.tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
        }

        // Calculate the concentrations with a perturbation TMPPRB to calculate the numerical Jacobian
        for (int k = 0; k < rttbl.num_stc - rttbl.num_min; k++)
        {
            tmpconc[k] += TMPPRB;
            for (int i = 0; i < rttbl.num_ssc; i++)
            {
                double tmpval = 0.0;
                for (int j = 0; j < rttbl.num_sdc; j++)
                {
                    tmpval += (tmpconc[j] + gamma[j]) * rttbl.dep_mtx[i][j];
                }
                tmpval -= rttbl.keq[i] + gamma[i + rttbl.num_stc];
                tmpconc[i + rttbl.num_stc] = tmpval;
            }
            for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
            {
                double tmpval = 0.0;
                for (int j = 0; j < rttbl.num_stc + rttbl.num_ssc; j++)
                {
                    tmpval += rttbl.conc_contrib[i][j] * pow(10, tmpconc[j]);
                }
                residue_t[i] = tmpval - (chms.tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
                jcb[k][i] = (residue_t[i] - residue[i]) * TMPPRB_INV;
            }
            tmpconc[k] -= TMPPRB;
        }
        
        for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
        {
            x[i] = -residue[i];
        }

        int pivot_flg = denseGETRF(jcb, rttbl.num_stc - rttbl.num_min, rttbl.num_stc - rttbl.num_min, p.data());
        if (pivot_flg != 0)
        {
            destroyMat(jcb);
            return 1;
        }

        denseGETRS(jcb, rttbl.num_stc - rttbl.num_min, p.data(), x.data());

        max_error = 0.0;
        for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
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

            max_error = std::max(fabs(residue[i] / tot_conc[i]), max_error);
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

    CheckChmsForNonFinite(chms, "react.c", 389);

    // Update secondary species here
    for (int i = 0; i < rttbl.num_ssc; i++)
    {
        double tmpval = 0.0;
        for (int j = 0; j < rttbl.num_sdc; j++)
        {
            tmpval += (tmpconc[j] + gamma[j]) * rttbl.dep_mtx[i][j];
        }
        tmpval -= rttbl.keq[i] + gamma[i + rttbl.num_stc];
        tmpconc[i + rttbl.num_stc] = tmpval;
    }
    CheckChmsForNonFinite(chms, "react.c", 402);
    
    for (int i = 0; i < rttbl.num_stc - rttbl.num_min; i++)
    {
        double tmpval = 0.0;
        for (int j = 0; j < rttbl.num_stc + rttbl.num_ssc; j++)
        {
            tmpval += rttbl.conc_contrib[i][j] * pow(10, tmpconc[j]);
        }
        tot_conc[i] = tmpval;
    }
    CheckChmsForNonFinite(chms, "react.c", 413);

    for (int i = 0; i < rttbl.num_stc + rttbl.num_ssc; i++)
    {
        if (i < rttbl.num_stc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                chms.tot_conc[i] += (rate_spe[i] + rate_spet[i]) * stepsize * 0.5;
                chms.prim_actv[i] = 1.0;
                chms.prim_conc[i] = chms.tot_conc[i];
            }
            else
            {
                chms.prim_conc[i] = pow(10, tmpconc[i]);
                chms.prim_actv[i] = pow(10, tmpconc[i] + gamma[i]);
                chms.tot_conc[i] = tot_conc[i];
            }
        }
        else
        {
            chms.sec_conc[i - rttbl.num_stc] = pow(10, tmpconc[i]);
        }
    }

    CheckChmsForNonFinite(chms, "react.c", 438);
    return 0;
}

double ReactControl(const array<ChemTableEntry, MAXSPS>& chemtbl, const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl,
    double stepsize, double porosity, double depth, double satn, double temp, double Zw, array<f64, MAXSPS>& react_rate,
    ChemicalState&chms)
{
    double          substep;
    double          step_counter = 0.0;
    double          conc0[MAXSPS];

    // Copy initial mineral concentration to array
    for (int kspc = 0; kspc < rttbl.num_min; kspc++)
    {
        conc0[kspc] = chms.tot_conc[kspc + rttbl.num_stc - rttbl.num_min];
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
        for (int kspc = 0; kspc < rttbl.num_min; kspc++)
        {
            // Calculate reaction rate (mole/m2-pm/day) = mol/L-pm/day * mm of depth   * (1m/1000mm) * (1000L/1m3)
            react_rate[kspc] =
                (chms.tot_conc[kspc + rttbl.num_stc - rttbl.num_min] - conc0[kspc]) * depth;
            CheckChmsForNonFinite(chms, "react.c", 477);
        }

        for (int kspc = 0; kspc <rttbl.num_spc; kspc++)
        {
            chms.tot_mol[kspc] = chms.tot_conc[kspc] * porosity * satn * depth; // tot_mol (moles-mm of water/L water ) =tot_conc (mol/L water) * satn * porosity * depth
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

double ReactSurfaceControl(const array<ChemTableEntry, MAXSPS>& chemtbl, const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl,
    double stepsize, double porosity, double depth, double tot_water, double temp, array<f64, MAXSPS>& react_rate,
    ChemicalState&chms)
{
    double          substep;
    double          step_counter = 0.0;
    double          conc0[MAXSPS];

    // Copy initial mineral concentration to array
    for (int kspc = 0; kspc < rttbl.num_min; kspc++)
    {
        conc0[kspc] = chms.tot_conc[kspc + rttbl.num_stc - rttbl.num_min];
    }

    substep = stepsize;

    while (1.0 - step_counter / stepsize > 1.0E-10 && substep > MINIMUM_SUBSTEP)
    {
        CheckChmsForNonFinite(chms, "react.c", 892);
        int flag = SolveSurfaceReact(substep, chemtbl, kintbl, rttbl, tot_water, temp, porosity, chms);
        CheckChmsForNonFinite(chms, "react.c", 894);

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
        for (int kspc = 0; kspc < rttbl.num_min; kspc++)
        {
            // Calculate reaction rate (mole/m2-pm/day) = mol/L-pm/day * mm of depth   * (1m/1000mm) * (1000L/1m3)
            react_rate[kspc] =
                (chms.tot_conc[kspc + rttbl.num_stc - rttbl.num_min] - conc0[kspc]) * depth;
            CheckChmsForNonFinite(chms, "react.c", 477);
        }

        for (int kspc = 0; kspc <rttbl.num_spc; kspc++)
        {
            chms.tot_mol[kspc] = chms.tot_conc[kspc] * tot_water; // tot_mol (moles-mm of water/L water ) =tot_conc (mol/L water) * satn * porosity * depth
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

void SoilMoistFactorRange(array<f64, MAXSPS>& dst, double satn, const array<f64, MAXSPS>& sw_threshold, const array<f64, MAXSPS>& sw_exponent, int start, int end, int offset) {
    /* Calculate the soil moisture factor over an array */
    if (satn < 1.0) {
        for (int i = start; i < end; i++) {
            dst[i - offset] = SoilMoistFactor(satn, sw_threshold[i], sw_exponent[i]);
        }
    }
    return;
}

void GetSurfaceAreaRange(array<f64, MAXSPS>& area, const array<f64, MAXSPS>& prim_conc, const array<f64, MAXSPS>& ssa, const array<ChemTableEntry, MAXSPS>& chemtbl, int start, int end, int offset) {
    for (int i = start; i < end; i++) {
        area[i - offset] = prim_conc[i] * ssa[i] * chemtbl[i].molar_mass;
    }
}

void GetTempFactorRange(array<f64, MAXSPS>& ftemp, const array<f64, MAXSPS>& q10, double temperature, int start, int end, int offset) {
    for (int i = start; i < end; i++) {
        ftemp[i - offset] = SoilTempFactor(q10[i], temperature);
    }
}

void GetWTDepthFactorRange(array<f64, MAXSPS>& fzw, double Zw, const array<f64, MAXSPS>& n_alpha, int start, int end, int offset) {
    for (int i = start; i < end; i++) {
        fzw[i - offset] = WTDepthFactor(Zw, n_alpha[i]);
    }
}