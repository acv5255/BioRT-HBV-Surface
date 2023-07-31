#include "biort.h"

// Calculates mass from transport
// To increase stability, transport is calculated using backward-in-time scheme, in which
//   d(CV) / dt = C_in * q_in - C * q_out + f
// i.e.,
//   C * V - C_0 * V_0 = C_in * Q_in - C * Q_out + F.
// Thus,
//   C = (C_in * Qin + F + C0 * V0) / (V + Q_out).
void Transport(int step, const ChemTableEntry *chemtbl, ReactionNetwork *rttbl, const ControlData ctrl, Subcatchment* subcatch)
{
    /* Transport for all species in each of the zones */
    TransportSurfaceZone(rttbl, ctrl, subcatch, step);
    
    // printf("Concentration after 'TransportSurfaceZone': \n");
    // PrintChemicalState(&subcatch->chms[SURFACE]);
    // printf("\n\n");

    // printf("UZ Concentration after 'TransportSurfaceZone': \n");
    // PrintChemicalState(&subcatch->chms[UZ]);
    // printf("\n\n");

    SolveSpeciation(chemtbl, ctrl, rttbl, 0, &subcatch->chms[UZ]);
    
    TransportShallowZone(rttbl, chemtbl, subcatch, step);

    SolveSpeciation(chemtbl, ctrl, rttbl, 0, &subcatch->chms[LZ]);
    
    TransportDeepZone(rttbl, chemtbl, subcatch, step);
}

void UpdatePrimConc(const ReactionNetwork *rttbl, const ControlData ctrl, Subcatchment* subcatch)
{
    /* Set concentrations */

    for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
    {
        subcatch->chms[SNOW].prim_conc[kspc] = subcatch->chms[SNOW].tot_conc[kspc];
        subcatch->chms[SURFACE].prim_conc[kspc] = subcatch->chms[SURFACE].tot_conc[kspc];
        if (ctrl.transport_only == true)
        {
            subcatch->chms[UZ].prim_conc[kspc] = subcatch->chms[UZ].tot_conc[kspc];
            subcatch->chms[LZ].prim_conc[kspc] = subcatch->chms[LZ].tot_conc[kspc];
        }
        subcatch->chms[STREAM].prim_conc[kspc] = subcatch->chms[STREAM].tot_conc[kspc];
    }
}

void TransportSurfaceZone(const ReactionNetwork* rttbl, const ControlData ctrl, Subcatchment* subcatch, const int step) {
    const double q_snow = subcatch->q[step][Psnow];
    const double q_snowmelt = subcatch->q[step][snowmelt];
    const double q_rain = subcatch->q[step][Prain];
    const double q_q0 = subcatch->q[step][Q0];
    const double q_q1 = subcatch->q[step][Q1];
    const double q_perc = subcatch->q[step][PERC];

    const double ws_surf = subcatch->ws[step][SURFACE];
    const double ws_snow = subcatch->ws[step][SNOW];
    const double ws_sz = subcatch->ws[step][UZ];       // Water storage in the shallow zone

    const double ws_snow_tot = ws_snow + q_snowmelt;
    const double ws_surf_tot = q_q0 + ws_surf;
    const double ws_sz_tot = ws_sz + q_q1 + q_perc;

    for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
    {

        double c_prcp_i;        // Concentration in precipitation of species (i)
        if (ctrl.variable_precipchem == false) {
            c_prcp_i = subcatch->prcp_conc[kspc];
        }
        else if (ctrl.variable_precipchem == true) {
            c_prcp_i = subcatch->prcp_conc_time[step][kspc];
        }
        // Precipitation
        double dm_rain = c_prcp_i * q_rain;             // Input of mass of chemical species i from rain - note: 1 mm water = 1 L water
        double dm_snow = c_prcp_i * q_snow;             // Change in mass of chemical species i in snow reservoir
        subcatch->chms[SNOW].tot_mol[kspc] += dm_snow;

        double c_snow = ZERO_CONC;
        double dm_snowmelt = 0.0;
        // Temporary concentration in snow
        if (ws_snow_tot == 0.0)
        {
            subcatch->chms[SNOW].tot_mol[kspc] = 0.0;
        }
        else
        {
            c_snow = subcatch->chms[SNOW].tot_mol[kspc] / ws_snow_tot;
            dm_snowmelt = c_snow * q_snowmelt;
 
            // Infiltration
            subcatch->chms[SNOW].tot_mol[kspc] -= dm_snowmelt;  // Mass flux (mol/d)
        }

        const double m_surf_0 = subcatch->chms[SURFACE].tot_mol[kspc];
        double c_in = (dm_rain + dm_snowmelt) / (q_rain + q_snowmelt);
        if (!isfinite(c_in)) c_in = ZERO_CONC;
        const double v_0 = ws_surf + q_q0;
        const double q_in = (q_rain + q_snowmelt);
        const double q_out = q_q0;
        const double delta_q = q_in - q_out;
        const double ws_tot = MAX(v_0 + delta_q, 0.001);
        const double delta_m = m_surf_0 - (c_in * ws_tot + (m_surf_0 - c_in * ws_tot) * exp(-(q_in * c_in) / ws_tot));

        // printf("m_surf_0: %g\n", m_surf_0);
        // printf("c_in: %g\n", c_in);
        // printf("q_in: %g\n", q_in);
        // printf("v_0: %g\n", v_0);
        // printf("q_out: %g\n", q_out);
        // printf("ws_tot: %g\n", ws_tot);
        // printf("delta_q: %g\n", delta_q);
        // printf("delta_m: %g\n", delta_m);
        // printf("\n");

        // subcatch->chms[SURFACE].tot_mol[kspc] += dm_rain + dm_snowmelt;
        // subcatch->chms[SURFACE].tot_conc[kspc] = subcatch->chms[SURFACE].tot_mol[kspc] / ws_surf_tot;

        // double c_surf = subcatch->chms[SURFACE].tot_conc[kspc];
        // double dm_surf = c_surf * q_q0;     // Mass mass transfer to the stream
        // double dm_inf = c_surf * (q_rain + q_snowmelt);
        const double m_surf_to_stream = delta_m * (q_q0) / (q_q0 + q_rain + q_snowmelt);        // Mass transfer to stream
        const double m_inf = delta_m * (q_rain + q_snowmelt) / (q_q0 + q_rain + q_snowmelt);    // Mass transfer of infiltration

        // Q0
        // subcatch->chms[SURFACE].tot_mol[kspc] -= dm_surf;
        // subcatch->chms[STREAM].tot_mol[kspc] = dm_surf;
        subcatch->chms[SURFACE].tot_mol[kspc] -= delta_m;
        subcatch->chms[STREAM].tot_mol[kspc] += m_surf_to_stream;
        subcatch->chms[UZ].tot_mol[kspc] += m_inf;

        // Infiltration
        // subcatch->chms[SURFACE].tot_mol[kspc] -= dm_inf;
        // subcatch->chms[UZ].tot_mol[kspc] += dm_inf;
        /*
            Note from Andrew - since HBV does not have an infiltration term, I just use all precipitation + snowmelt
            entering from the surface zone to the shallow zone. That way, the entire water volume can react both at the
            surface and then move into the shallow zone and because I'm not sure how to remove the dependence on the surface zone
         */
        subcatch->chms[SURFACE].tot_conc[kspc] = subcatch->chms[SURFACE].tot_mol[kspc] / ws_surf_tot;
        subcatch->chms[UZ].tot_conc[kspc] = subcatch->chms[UZ].tot_mol[kspc] / ws_sz_tot;
        
        // Just in case
        // subcatch->chms[SURFACE].tot_conc[kspc] = MIN(subcatch->chms[SURFACE].tot_conc[kspc], ZERO_CONC);

        // Newer version:
    }

    return;
}

void TransportShallowZone(const ReactionNetwork* rttbl, const ChemTableEntry chemtbl[], Subcatchment* subcatch, const int step) {
    const double q_q1 = subcatch->q[step][Q1];
    const double q_q2 = subcatch->q[step][Q2];
    const double q_perc = subcatch->q[step][PERC];
    const double ws_dz = subcatch->ws[step][LZ];
    const double ws_dz_tot = ws_dz + q_q2;

    for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
    {
        // Temporary concentration in UZ
        
        // Q1 -- transport from the upper zone to the stream
        subcatch->chms[UZ].tot_mol[kspc] -= subcatch->chms[UZ].tot_conc[kspc] * q_q1;
        subcatch->chms[STREAM].tot_mol[kspc] += subcatch->chms[UZ].tot_conc[kspc] * q_q1;
        
        // Percolation
        subcatch->chms[UZ].tot_mol[kspc] -= subcatch->chms[UZ].tot_conc[kspc] * q_perc;
        subcatch->chms[LZ].tot_mol[kspc] += subcatch->chms[UZ].tot_conc[kspc] * q_perc;
        
        if (chemtbl[kspc].mtype == MIXED_MA)
        {
            for (int kssc = 0; kssc < rttbl->num_ssc; kssc++)
            {
                const double nu_i = rttbl->conc_contrib[kspc][rttbl->num_stc+kssc]; // Stoichiometric coefficient in each species
                const double c_i = subcatch->chms[UZ].sec_conc[kssc];
                if ((nu_i != 0) && 
                    (chemtbl[rttbl->num_stc+kssc].itype != AQUEOUS))
                    {
                        subcatch->chms[UZ].tot_mol[kspc] += c_i * q_q1 * nu_i;
                        subcatch->chms[STREAM].tot_mol[kspc] -= c_i * q_q1 * nu_i;
                        subcatch->chms[UZ].tot_mol[kspc] += c_i * q_perc * nu_i;
                        subcatch->chms[LZ].tot_mol[kspc] -= c_i * q_perc * nu_i;
                        //biort_printf(VL_NORMAL, "Working\n");
                    }
            }   
        }

        // Temporary concentration in LZ
        subcatch->chms[LZ].tot_conc[kspc] = subcatch->chms[LZ].tot_mol[kspc] / ws_dz_tot;   
    }

    return;
}

void TransportDeepZone(const ReactionNetwork* rttbl, const ChemTableEntry chemtbl[], Subcatchment* subcatch, const int step) {
    const double q_q0 = subcatch->q[step][Q0];
    const double q_q1 = subcatch->q[step][Q1];
    const double q_q2 = subcatch->q[step][Q2];
    const double ws_surf = subcatch->ws[step][SURFACE];
    const double ws_snow = subcatch->ws[step][SNOW];
    const double ws_sz = subcatch->ws[step][UZ];       // Water storage in the shallow zone
    const double ws_dz = subcatch->ws[step][LZ];       // Water storage in the deep zone

    for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
    {
        {
            // Q2
            subcatch->chms[LZ].tot_mol[kspc] -= subcatch->chms[LZ].tot_conc[kspc] * q_q2;
            subcatch->chms[STREAM].tot_mol[kspc] += subcatch->chms[LZ].tot_conc[kspc] * q_q2;
            
            if (chemtbl[kspc].mtype == MIXED_MA)
            {
                for (int kssc = 0; kssc < rttbl->num_ssc; kssc++)
                {
                    double nu_i = rttbl->conc_contrib[kspc][rttbl->num_stc+kssc];   // Contribution of each species to this value
                    double c_i = subcatch->chms[LZ].sec_conc[kssc];     // Concentration of species 'kspc'
                    double dm_i = nu_i * c_i * q_q2;                    // Mass transfer between zones
                    if ((nu_i != 0) && 
                        (chemtbl[rttbl->num_stc+kssc].itype != AQUEOUS))
                        {
                            subcatch->chms[LZ].tot_mol[kspc] += dm_i;
                            subcatch->chms[STREAM].tot_mol[kspc] -= dm_i;
                        }
                    
                }   
            }
        }
        

        // UPDATE CONCENTRATIONS
        const double q_stream = q_q0 + q_q1 + q_q2;
        subcatch->chms[SNOW].tot_conc[kspc] =
            (ws_snow == 0) ? ZERO_CONC : subcatch->chms[SNOW].tot_mol[kspc] / ws_snow;
        subcatch->chms[SURFACE].tot_conc[kspc] = subcatch->chms[SURFACE].tot_mol[kspc] / (ws_surf + q_q0);
        subcatch->chms[UZ].tot_conc[kspc] = subcatch->chms[UZ].tot_mol[kspc] / ws_sz;
        subcatch->chms[LZ].tot_conc[kspc] = subcatch->chms[LZ].tot_mol[kspc] / ws_dz;

        subcatch->chms[STREAM].tot_conc[kspc] = (q_stream > 0.0) ?
            subcatch->chms[STREAM].tot_mol[kspc] / q_stream : ZERO_CONC;
    }
}