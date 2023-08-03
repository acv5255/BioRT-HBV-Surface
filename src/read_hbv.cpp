#include "biort.hpp"

const double SURFACE_STORAGE_MIN = 1e-5; // 1 mm water, minimum water allowed in the surface zone

void ReadHbvResults(const char dir[], int *nsteps, int *steps[], Subcatchment* subcatch, int mode)
{
    FILE           *file_pointer;
    char            file_name[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             line_number = 0;
    int             len_numexp=1;
    int             numexp_file_flag=0;

    if (mode==0)
    {
        sprintf(file_name, "input/%s/Results.txt", dir);
        file_pointer = fopen(file_name, "r");
    } else if (mode == 1)
    {
        sprintf(file_name, "input/%s/Numexp_Results.txt", dir);
        if (access(file_name , F_OK ) != -1)
        {
            file_pointer = fopen(file_name, "r");
            numexp_file_flag = 1;
            biort_printf(VL_NORMAL, "\nHydrology in %s will be used for numerical experiment.\n", file_name);
        } else
        {
            sprintf(file_name, "input/%s/Numexp_precipchem.txt", dir);
            file_pointer = fopen(file_name, "r");
            len_numexp = CountLines(file_pointer, cmdstr, 0) - 1;
            fclose(file_pointer);
            sprintf(file_name, "input/%s/Results.txt", dir);
            file_pointer = fopen(file_name, "r");
            numexp_file_flag = 0;
            biort_printf(VL_NORMAL, "\nHydrology in %s will be used for numerical experiment.\n", file_name);
        }
    }


    *nsteps = CountLines(file_pointer, cmdstr, 0) - 1;

    rewind(file_pointer);

    *steps = (int *)malloc(*nsteps * sizeof(int));

    if ((mode == 1) & (numexp_file_flag == 0))
    {
      if (len_numexp % *nsteps != 0){
          biort_printf(VL_ERROR,"\nNumber of time steps in \"Numexp_precipchem.txt\" should be a multiple of time steps in \"Results.txt\" file, if \"Numexp_precipchem.txt\" file is not provided.\n");
          exit(EXIT_FAILURE);
    }
        len_numexp /= *nsteps;
    }


    subcatch->ws = (double (*)[NWS])malloc(len_numexp * *nsteps * sizeof(double[NWS]));
    subcatch->q = (double (*)[NQ])malloc(len_numexp * *nsteps * sizeof(double[NQ]));
    subcatch->tmp = (double *)malloc(len_numexp * *nsteps * sizeof(double));

    NextLine(file_pointer, cmdstr, &line_number);     // Skip header line

    for (int kstep = 0; kstep < *nsteps; kstep++)
    {
        double          snow, sm;

        fscanf(file_pointer, "%d", &((*steps)[kstep]));    // Read model steps

        fscanf(file_pointer, "%*f %*f");  // Skip "Qsim" and "Qobs"
        fscanf(file_pointer, "%lf %lf", &subcatch->q[kstep][PRECIP], &subcatch->tmp[kstep]);    // Read precip. and
                                                                                                // air temperature
        fscanf(file_pointer, "%*f %*f");   // Skip "AET" and PET"
        fscanf(file_pointer, "%lf %*f %lf", &snow, &sm);     // Read snow and soil moisture
        subcatch->ws[kstep][SNOW] = snow;  // 2021-05-14
        subcatch->ws[kstep][SM] = sm;
        fscanf(file_pointer, "%lf", &subcatch->q[kstep][RECHG]);  // Read recharge
        fscanf(file_pointer, "%*f %*f");   // Skip upper and lower zone storages
        fscanf(file_pointer, "%lf %lf %lf", &subcatch->q[kstep][Q0],
            &subcatch->q[kstep][Q1], &subcatch->q[kstep][Q2]);  // Read Q0, Q1, and Q2
        fscanf(file_pointer, "%*f %*f");   // Skip "Qsim_rain" and "Qsim_snow"

        //Incoming precipitation might become snow or enter the soil zone directly.
        //  Psnow = PRECIP * SFCF (if temp < TT)
        //  Prain = PRECIP (if temp >=TT)
        //  Psnow + snow0 = snowmelt + snow
        //
        // Solving the above equations,
        //
        subcatch->q[kstep][Prain] = (subcatch->tmp[kstep] < subcatch->tt) ?
            0.0 : subcatch->q[kstep][PRECIP];
        
        subcatch->q[kstep][Psnow] = (subcatch->tmp[kstep] < subcatch->tt) ?
            subcatch->q[kstep][PRECIP] * subcatch->sfcf : 0.0;
        
        subcatch->q[kstep][snowmelt] = (kstep == 0) ?
            0.0 : MAX(0, subcatch->q[kstep][Psnow] + subcatch->ws[kstep - 1][SNOW] - subcatch->ws[kstep][SNOW]);
        //biort_printf(VL_NORMAL, "  snowmelt %.2f m3 m-3.\n", subcatch->q[kstep][snowmelt]);
        // In the upper zone, HBV first calculates percolation to the lower zone.
        //   percolation = MIN(perc0, SUZ0 + recharge).                         (1)
        // If there is water left in the upper zone, calculate Q0 and Q1:
        //   Q0 = MAX(SUZ0 + recharge - percoloation - UZL, 0) * K0,            (2)
        //   Q1 = (SUZ0 + recharge - percoloation) * K1,                        (3)
        //   SUZ = SUZ0 + recharge - percolation - Q0 - Q1,                     (4)
        // i.e.,
        //   SUZ = Q1 / K1 - Q0 - Q1.                                           (5)
        // In the lower zone, percoloation is first added to the lower zone storage. Then Q2 is calculated:
        //   Q2 = (SLZ0 + percolation) * K2,                                    (6)
        //   SLZ = SLZ0 + percololation - Q2,                                   (7)
        // i.e.,
        //   SLZ = Q2 / k2 - Q2.                                                (8)
        //
        // HBV Light outputs SUZ and SLZ with low precision (one decimal digit), and does not output percolation
        // rate. Therefore, SUZ, SLZ, and percolation rate are calculated using Equations (1), (4), and (8).

        subcatch->ws[kstep][UZ] = subcatch->q[kstep][Q1] / subcatch->k1 -
            (subcatch->q[kstep][Q0] + subcatch->q[kstep][Q1]);
        subcatch->ws[kstep][LZ] = subcatch->q[kstep][Q2] / subcatch->k2 -
            subcatch->q[kstep][Q2];
        subcatch->q[kstep][PERC] = (kstep == 0) ?
            0.0 : MIN(subcatch->perc, subcatch->ws[kstep - 1][UZ] + subcatch->q[kstep][RECHG]);

    }
    
    //Change PERC value for 1st time step using water storage for last timestep
    
    subcatch->q[0][PERC] = MIN(subcatch->perc, subcatch->ws[*nsteps - 1][UZ] + subcatch->q[0][RECHG]);

    // Add 1. residual moisture to LZ & UZ and 2. SM to UZ
    for (int kstep = 0; kstep < *nsteps; kstep++)
    {
        subcatch->ws[kstep][SURFACE]+= MIN(subcatch->soil_surface.ws_passive, STORAGE_MIN);
        // if (subcatch->q[kstep][Q0] > 0.0) subcatch->ws[kstep][SURFACE] = subcatch->q[kstep][Q0];
        subcatch->ws[kstep][UZ] += subcatch->soil_sz.ws_passive;
        subcatch->ws[kstep][LZ] += subcatch->soil_dz.ws_passive;
        subcatch->ws[kstep][UZ] += subcatch->ws[kstep][SM];
        
    }

    for (int kstep = 0; kstep < *nsteps; kstep++)
    {
        if (subcatch->ws[kstep][UZ] > (subcatch->soil_sz.depth * subcatch->soil_sz.porosity))
        {
            biort_printf(VL_ERROR, "\nWater storage in UZ exceeds maximum water storage capacity at line %d.\n", kstep);
            exit(EXIT_FAILURE);
        }
        if (subcatch->ws[kstep][LZ] > (subcatch->soil_dz.depth * subcatch->soil_dz.porosity))
        {
            biort_printf(VL_ERROR, "\nWater storage in LZ exceeds maximum water storage capacity at line %d.\n", kstep);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file_pointer);

    if ((mode == 1) & (numexp_file_flag == 0))
    {
        if (len_numexp > 1)
        {
            for (int i = 1; i < len_numexp; i++)
            {
                for (int kstep = 0; kstep < *nsteps; kstep++)
                {
                    subcatch->ws[(i * *nsteps)+kstep][SNOW] = subcatch->ws[kstep][SNOW];
                    subcatch->ws[(i * *nsteps)+kstep][SURFACE] = subcatch->ws[kstep][SURFACE];
                    subcatch->ws[(i * *nsteps)+kstep][UZ] = subcatch->ws[kstep][UZ];
                    subcatch->ws[(i * *nsteps)+kstep][LZ] = subcatch->ws[kstep][LZ];
                    subcatch->q[(i * *nsteps)+kstep][PRECIP] = subcatch->q[kstep][PRECIP];
                    subcatch->q[(i * *nsteps)+kstep][RECHG] = subcatch->q[kstep][RECHG];
                    subcatch->q[(i * *nsteps)+kstep][PERC] = subcatch->q[kstep][PERC];
                    subcatch->q[(i * *nsteps)+kstep][Q0] = subcatch->q[kstep][Q0];
                    subcatch->q[(i * *nsteps)+kstep][Q1] = subcatch->q[kstep][Q1];
                    subcatch->q[(i * *nsteps)+kstep][Q2] = subcatch->q[kstep][Q2];
                    subcatch->tmp[(i * *nsteps)+kstep] = subcatch->tmp[kstep];
                }
            }
        }
    }
    *nsteps *= len_numexp;
}

void ReadHbvParam(const char dir[], Subcatchment* subcatch)
{
    FILE           *file_pointer;
    char            file_name[MAXSTRING];
    char            cmdstr[MAXSTRING];
    char            tag[MAXSTRING];
    int             line_number = 0;
    double          value;

    sprintf(file_name, "input/%s/Parameter.xml", dir);
    file_pointer = fopen(file_name, "r");

    biort_printf(VL_NORMAL, "\nHBV MODEL PARAMETERS\n");

    FindLine(file_pointer, "<SubCatchmentParameters>", &line_number, cmdstr);

    while (!feof(file_pointer))
    {
        NextLine(file_pointer, cmdstr, &line_number);
        ParseLine(cmdstr, tag, &value);
        if (strcmp(tag, "PERC") == 0)
        {
            subcatch->perc = value;
            biort_printf(VL_NORMAL, "  Percolation rate is %.2lf mm day-1\n", subcatch->perc);
        }
        else if (strcmp(tag, "K1") == 0)
        {
            subcatch->k1 = value;
            biort_printf(VL_NORMAL, "  K1 is %.2lf day-1\n", subcatch->k1);
        }
        else if (strcmp(tag, "K2") == 0)
        {
            subcatch->k2 = value;
            biort_printf(VL_NORMAL, "  K2 is %.2lf day-1\n", subcatch->k2);
        }
        else if (strcmp(tag, "MAXBAS") == 0)
        {
            subcatch->maxbas = value;
            biort_printf(VL_NORMAL, "  Routing parameter is %.2lf\n", subcatch->maxbas);
            break;
        }
    }
    FindLine(file_pointer, "<SubCatchmentVegetationZoneParameters>", &line_number, cmdstr);

    while (!feof(file_pointer))
    {
        NextLine(file_pointer, cmdstr, &line_number);
        ParseLine(cmdstr, tag, &value);
        if (strcmp(tag,"TT") == 0)
        {
            subcatch->tt = value;
            biort_printf(VL_NORMAL, "  TT is %.2lf\n", subcatch->tt);
        }
        else if (strcmp(tag,"SFCF") == 0)
        {
            subcatch->sfcf = value;
            biort_printf(VL_NORMAL, "  SFCF is %.2lf\n", subcatch->sfcf);
            break;
        }
    }

    fclose(file_pointer);
}

void ParseLine(const char cmdstr[], char tag[], double *value)
{
    char           *temp;
    int             tag_index1, tag_index2, tag_index3;
    int             k;

    temp = strchr((char*)cmdstr, '<');
    tag_index1 = (int)(temp - cmdstr);

    temp = strchr((char*)cmdstr, '>');
    tag_index2 = (int)(temp - cmdstr);

    temp = strchr(temp, '<');
    tag_index3 = (int)(temp - cmdstr);

    for (k = tag_index1 + 1; k < tag_index2; k++)
    {
        tag[k - tag_index1 - 1] = cmdstr[k];
    }
    tag[tag_index2 - tag_index1 - 1] = '\0';

    for (k = tag_index2 + 1; k < tag_index3; k++)
    {
        temp[k - tag_index2 - 1] = cmdstr[k];
    }
    temp[tag_index3 - tag_index2 - 1] = '\0';
    sscanf(temp, "%lf", value);
}
