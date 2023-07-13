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

    const double ws_snow = subcatch->ws[step][SNOW];
    const double ws_sz = subcatch->ws[step][UZ];       // Water storage in the shallow zone

    const double ws_snow_tot = ws_snow + q_snowmelt;
    const double ws_sz_tot = ws_sz + q_q1 + q_perc;

    double conc_temp;

    for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
    {
        double c_prcp_i;
        if (ctrl.variable_precipchem == false) {
            c_prcp_i = subcatch->prcp_conc[kspc];
        }
        else if (ctrl.variable_precipchem == true) {
            c_prcp_i = subcatch->prcp_conc_time[step][kspc];
        }
        // Precipitation
        subcatch->chms[SNOW].tot_mol[kspc] += c_prcp_i * q_snow;

        // Temporary concentration in snow
        if (ws_snow_tot == 0.0)
        {
            // conc_temp = 0;
            subcatch->chms[SNOW].tot_mol[kspc] = 0;
        }
        else
        {
            conc_temp = subcatch->chms[SNOW].tot_mol[kspc] / ws_snow_tot;

        // Recharge
            subcatch->chms[SNOW].tot_mol[kspc] -= conc_temp * q_snowmelt;
        }

        if ((q_snowmelt + q_rain) == 0)
        {
            conc_temp = ZERO_CONC;
        } else // variable or constant precipitation chemistry mode
        {
            subcatch->chms[UZ].tot_mol[kspc] += ((conc_temp * q_snowmelt) + (c_prcp_i * q_rain));
            conc_temp = ((conc_temp * q_snowmelt) + (c_prcp_i * q_rain)) / (q_snowmelt + q_rain);
        }

        subcatch->chms[SURFACE].tot_conc[kspc] = conc_temp;
        
        // Q0
        subcatch->chms[UZ].tot_mol[kspc] -= conc_temp * q_q0;
        subcatch->chms[STREAM].tot_mol[kspc] = conc_temp * q_q0;
                
        subcatch->chms[UZ].tot_conc[kspc] =  subcatch->chms[UZ].tot_mol[kspc] / ws_sz_tot;
    }
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
                if ((rttbl->conc_contrib[kspc][rttbl->num_stc+kssc] != 0) && 
                    (chemtbl[rttbl->num_stc+kssc].itype != AQUEOUS))
                    {
                        subcatch->chms[UZ].tot_mol[kspc] += subcatch->chms[UZ].sec_conc[kssc] * q_q1 * rttbl->conc_contrib[kspc][rttbl->num_stc+kssc];
                        subcatch->chms[STREAM].tot_mol[kspc] -= subcatch->chms[UZ].sec_conc[kssc] * q_q1 * rttbl->conc_contrib[kspc][rttbl->num_stc+kssc];
                        subcatch->chms[UZ].tot_mol[kspc] += subcatch->chms[UZ].sec_conc[kssc] * q_perc * rttbl->conc_contrib[kspc][rttbl->num_stc+kssc];
                        subcatch->chms[LZ].tot_mol[kspc] -= subcatch->chms[UZ].sec_conc[kssc] * q_perc * rttbl->conc_contrib[kspc][rttbl->num_stc+kssc];
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
                    if ((rttbl->conc_contrib[kspc][rttbl->num_stc+kssc] != 0) && 
                        (chemtbl[rttbl->num_stc+kssc].itype != AQUEOUS))
                        {
                            subcatch->chms[LZ].tot_mol[kspc] += subcatch->chms[LZ].sec_conc[kssc] * q_q2 * rttbl->conc_contrib[kspc][rttbl->num_stc+kssc];
                            subcatch->chms[STREAM].tot_mol[kspc] -= subcatch->chms[LZ].sec_conc[kssc] * q_q2 * rttbl->conc_contrib[kspc][rttbl->num_stc+kssc];
                        }
                    
                }   
            }
        }
        

        // UPDATE CONCENTRATIONS
        const double q_stream = q_q0 + q_q1 + q_q2;
        subcatch->chms[SNOW].tot_conc[kspc] =
            (ws_snow == 0) ? ZERO_CONC : subcatch->chms[SNOW].tot_mol[kspc] / ws_snow;
        //subcatch->chms[SURFACE].tot_conc[kspc] =
        //    (ws_snow == 0) ? ZERO_CONC : subcatch->chms[SURFACE].tot_mol[kspc] / subcatch->q[step][SURFACE];
        subcatch->chms[UZ].tot_conc[kspc] =
            subcatch->chms[UZ].tot_mol[kspc] / ws_sz;
        subcatch->chms[LZ].tot_conc[kspc] =
            subcatch->chms[LZ].tot_mol[kspc] / ws_dz;

        subcatch->chms[STREAM].tot_conc[kspc] = (q_stream > 0.0) ?
            subcatch->chms[STREAM].tot_mol[kspc] / q_stream : ZERO_CONC;
        
        
        //subcatch->chms[STREAM].tot_conc[kspc] =
        //    (q_stream > 0.0) ?
        //    (subcatch->chms[SURFACE].tot_conc[kspc] * q_q0 +  subcatch->chms[UZ].tot_conc[kspc] * q_q1 + subcatch->chms[LZ].tot_conc[kspc] * q_q2) / (q_q0 + q_q1 + q_q2) : ZERO_CONC;
    }
}