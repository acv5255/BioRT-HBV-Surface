#include "biort.hpp"

const double SURFACE_MIN = 1e-3;

double WellMixedConcentration(double c_in, double c_0, double q, double v_0, double t) {
    /*
        Calculate the concentration exiting a well-mixed bucket
     */
    return c_in + (c_0 - c_in) * exp(-t * (q / (q + v_0)));
}

double TransportMassChange(const double c_in, const double c_0, const double q, const double v_0, const double t) {
    /*
        Calculate mass change after transport
        q: flux in and out of zone (must be equal)
        c_in: concentration entering the bucket
        c_out: concentration exiting the bucket
        v_0: passive volume of bucket
        t: time step to use
     */ 
    return WellMixedConcentration(c_in, c_0, q, v_0, t) * (q + v_0);
}

// Calculates mass from transport
// To increase stability, transport is calculated using backward-in-time scheme, in which
//   d(CV) / dt = C_in * q_in - C * q_out + f
// i.e.,
//   C * V - C_0 * V_0 = C_in * Q_in - C * Q_out + F.
// Thus,
//   C = (C_in * Qin + F + C0 * V0) / (V + Q_out).
void Transport(int step, const array<ChemTableEntry, MAXSPS>& chemtbl, ReactionNetwork *rttbl, const ControlData ctrl, Subcatchment* subcatch)
{
    /* Transport for all species in each of the zones */
    
    TransportSurfaceZone(rttbl, ctrl, subcatch, step);
    
    SolveSpeciation(chemtbl, ctrl, rttbl, 0, &subcatch->chms[SURFACE]);
    SolveSpeciation(chemtbl, ctrl, rttbl, 0, &subcatch->chms[UZ]);

    CheckChmsForNonFinite(&subcatch->chms[UZ], "react.c", 54);
    
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
    const double t = 1.0;     // Time step for transportation
    const double q_snow = subcatch->q[step][Psnow];
    const double q_snowmelt = subcatch->q[step][snowmelt];
    const double q_rain = subcatch->q[step][Prain];
    const double q_q0 = subcatch->q[step][Q0];
    const double q_q1 = subcatch->q[step][Q1];
    const double q_perc = subcatch->q[step][PERC];

    const double ws_surf = std::max(subcatch->ws[step][SURFACE], 0.001);
    const double v_0 = ws_surf;
    const double ws_snow = subcatch->ws[step][SNOW];
    const double ws_sz = subcatch->ws[step][UZ];       // Water storage in the shallow zone

    const double ws_snow_tot = ws_snow + q_snowmelt;
    const double ws_surf_tot = q_q0 + ws_surf;
    const double ws_sz_tot = ws_sz + q_q1 + q_perc;

    const double q = (q_rain + q_snowmelt);
    const bool do_surface_transport = q > SURFACE_MIN;   // Only do surface transport of there is a threshold
    const double q_inf = (q_rain + q_snowmelt) - q_q0;      // The residual water enters the shallow zone
    const double x_0 = do_surface_transport ? q_q0 / q : 0.0;    // Fraction of water input going to the stream
    const double x_inf = do_surface_transport ? q_inf / q : 1.0; // Fraction of water input entering as infiltration
        
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

        // Mass transfer section
        const double m_0 = subcatch->chms[SURFACE].tot_mol[kspc];  // Initial number of moles
        double c_in = (dm_rain + dm_snowmelt) / (q_rain + q_snowmelt);  // Weighted average of snowmelt and rain concentrations
        if (std::isinf(c_in)) c_in = ZERO_CONC;
        
        const double c_0 = m_0 / (v_0 + q);

        // const double dm = TransportMassChange(q, c_in, c_0, v_0, t);
        const double c_f = WellMixedConcentration(c_in, c_0, q, v_0, t);
        const double m = c_f * (q + v_0);
        const double dm = m - m_0;
        const double m_in = c_in * q * t;
        const double m_out = m_in - dm;
                
        // Diagnostic values //
        const double m_f_surf = m_0 + dm;
        const double c_final = m_f_surf / (v_0 + q);

        const double new_mass_stream = subcatch->chms[STREAM].tot_mol[kspc] + x_0 * fabs(m_out);
        const double new_mass_sz = subcatch->chms[UZ].tot_mol[kspc] + x_inf * fabs(m_out);

        bool should_increase = (c_in > c_0);
        bool did_increase = (c_final > c_0);
        bool expected_behavior = (should_increase == did_increase);



        if (!expected_behavior) {
            // This means that the concentration should have increased but really decreased, or the opposite
            printf("Failed to get expected behavior in transport in surface zone...\n");
            exit(-1);
        }

        // Transport moles between reservoirs
        subcatch->chms[SURFACE].tot_mol[kspc] = m_f_surf;
        subcatch->chms[STREAM].tot_mol[kspc] = new_mass_stream;
        subcatch->chms[UZ].tot_mol[kspc] = new_mass_sz;
        // Infiltration
        /*
            Note from Andrew - since HBV does not have an infiltration term, I just use all precipitation + snowmelt
            entering from the surface zone to the shallow zone. That way, the entire water volume can react both at the
            surface and then move into the shallow zone and because I'm not sure how to remove the dependence on the surface zone
         */
        subcatch->chms[SURFACE].tot_conc[kspc] = subcatch->chms[SURFACE].tot_mol[kspc] / ws_surf_tot;
        subcatch->chms[UZ].tot_conc[kspc] = subcatch->chms[UZ].tot_mol[kspc] / ws_sz_tot;
    }

    return;
}

void TransportShallowZone(const ReactionNetwork* rttbl, const array<ChemTableEntry, MAXSPS>& chemtbl, Subcatchment* subcatch, const int step) {
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
                const double nu_i = rttbl->conc_contrib[kspc][rttbl->num_stc+kssc];     // Stoichiometric coefficient in each species
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

    CheckChmsForNonFinite(&subcatch->chms[UZ], "transpt.c", 220);

    return;
}

void TransportDeepZone(const ReactionNetwork* rttbl, const array<ChemTableEntry, MAXSPS>& chemtbl, Subcatchment* subcatch, const int step) {
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
        
        CheckChmsForNonFinite(&subcatch->chms[LZ], "transpt.c", 259);

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

    CheckChmsForNonFinite(&subcatch->chms[LZ], "transpt.c", 274);

}