#include "biort.h"

int             verbose_mode;

int main(int argc, char *argv[])
{
    char                dir[MAXSTRING];         // name of input directory
    char                file_name[MAXSTRING];          // name of output file
    char                timestr[MAXSTRING];     // time stamp
    int                 nsteps;                 // number of simulation steps
    int                 *steps;                  // simulation steps
    int                 nsteps_numexp;                 // number of simulation steps
    int                 *steps_numexp;                  // simulation steps
    time_t              rawtime;
    struct tm           *timestamp;
    ReactionNetwork     rttbl;
    ChemTableEntry      chemtbl[MAXSPS];
    KineticTableEntry   kintbl[MAXSPS];
    CalibrationStruct   calib;
    ControlData         ctrl;

    // Read command line arguments
    ParseCmdLineParam(argc, argv, dir);

    // Allocate
    Subcatchment subcatch;
    Subcatchment subcatch_numexp;

    calib.rate = 0.0;
    calib.xsorption = 0.0;
    calib.ssa = 1.0;

    // Read subcatchment soil parameters
    ReadSoil(dir, &subcatch);

    // Read HBV simulation steps, water states and fluxes
    ReadHbvParam(dir, &subcatch);
    ReadHbvResults(dir, &nsteps, &steps, &subcatch, 0);

    // Read chemistry control file
    ReadChem(dir, &ctrl, &rttbl, chemtbl, kintbl);

    // Read chemistry initial conditions
    ReadCini(dir, chemtbl, &rttbl, &subcatch);

    // printf("Surface zone chemistry after 'ReadCini': \n");
    // PrintChemicalState(&subcatch.chms[SURFACE]);

    // Read time-series precipitation chemistry if defined in chem.txt  2021-05-20
    if (ctrl.variable_precipchem == 1) {
        ReadPrecipChem(dir, &nsteps, &steps, &subcatch, rttbl.num_stc, chemtbl, 0);
    }

    if (ctrl.precipchem_numexp == 1){
        ctrl.recycle = ctrl.recycle - 1;
        CopyConstSubcatchProp(&subcatch, &subcatch_numexp);
        ReadHbvResults(dir, &nsteps_numexp,  &steps_numexp, &subcatch_numexp, 1);
        ReadPrecipChem(dir, &nsteps_numexp, &steps_numexp, &subcatch_numexp, rttbl.num_stc, chemtbl, 1);
    }

    // Initialize RT structures
    InitChem(dir, &calib, ctrl, chemtbl, kintbl, &rttbl, &subcatch);

    // Check on secondary concentrations in the surface zone
    ErrOnZeroRanged("main.c", "subcatch.chms[SURFACE].sec_conc", 63, subcatch.chms[SURFACE].sec_conc, 2);

    // Create output directory when necessary
    mkdir("output");

    // {
    //     printf("Surface total concentration at before starting: \n");
    //     PrintArray(subcatch.chms[SURFACE].tot_conc);
    // }

    //
    time(&rawtime);
    timestamp = localtime(&rawtime);
    strftime(timestr, 11, "%y%m%d%H%M", timestamp);

    // printf("Surface zone chemistry before starting loops: \n");
    // PrintChemicalState(&subcatch.chms[SURFACE]);

    // Open output file
    {
        FILE* file_pointer;
        char dir_substr[SUBSTRING_SIZE];
        char time_substr[SUBSTRING_SIZE];
        strncpy(dir_substr, dir, SUBSTRING_SIZE);
        strncpy(time_substr, timestr, SUBSTRING_SIZE);

        // sprintf(file_name, "output/%s_results_%s.txt", dir, timestr);
        sprintf(file_name, "output/%s_results_%s.txt", dir_substr, time_substr);
        file_pointer = fopen(file_name, "w");
        PrintHeader(file_pointer, ctrl.transport_only, &rttbl, chemtbl);
    
        biort_printf(VL_NORMAL, "\nHBV-BioRT %s simulation started.\n", dir);

        // Loop through forcing cycles
        for (int kcycle = 0; kcycle <= ctrl.recycle; kcycle++)
        {
            // Loop through model steps to calculate reactive transport
            for (int kstep = 0; kstep < nsteps; kstep++)
            {
                // Transport and routing
                Transport(kstep, chemtbl, &rttbl, ctrl, &subcatch); // 2021-05-20
                // printf("Finished Transport for step %d\n", kstep);
                // Transport changes total concentrations. Primary concentrations needs to be updated using total
                // concentrations
                
                UpdatePrimConc(&rttbl, ctrl, &subcatch);
                
                StreamSpeciation(kstep, chemtbl, ctrl, &rttbl, &subcatch);

                if (ctrl.transport_only == KIN_REACTION)
                {
                    // In reaction mode, simulate reaction for soil, and speciation for stream
                    Reaction(kstep, 86400.0, chemtbl, kintbl, &rttbl, &subcatch);
                }
                else
                {
                    Speciation(chemtbl, ctrl, &rttbl, &subcatch);
                }
                PrintDailyResults(file_pointer, ctrl.transport_only, steps[kstep], &rttbl, &subcatch);
                ResetReactionRates(&subcatch);
            }
        }

        fclose(file_pointer);
    }

    biort_printf(VL_NORMAL, "\nHBV-BioRT %s simulation succeeded.\n", dir);


    if (ctrl.precipchem_numexp == 1){
        CopyInitChemSubcatch(&rttbl, &subcatch, &subcatch_numexp);
        {
            FILE* file_pointer;
            char dir_substr[SUBSTRING_SIZE];
            char time_substr[SUBSTRING_SIZE];
            strncpy(dir_substr, dir, SUBSTRING_SIZE);
            strncpy(time_substr, timestr, SUBSTRING_SIZE);

            // sprintf(file_name, "output/%s_results_numexp_%s.txt", dir, timestr);
            sprintf(file_name, "output/%s_results_numexp_%s.txt", dir_substr, time_substr);
            file_pointer = fopen(file_name, "w");
            PrintHeader(file_pointer, ctrl.transport_only, &rttbl, chemtbl);

            biort_printf(VL_NORMAL, "\nHBV-BioRT %s numerical experiment started.\n", dir);

            for (int kstep = 0; kstep < nsteps_numexp; kstep++){


                Transport(kstep, chemtbl, &rttbl, ctrl, &subcatch_numexp); // 2021-05-20

                // Transport changes total concentrations. Primary concentrations needs to be updated using total
                // concentrations
                UpdatePrimConc(&rttbl, ctrl, &subcatch_numexp);
                
                StreamSpeciation(kstep, chemtbl, ctrl, &rttbl, &subcatch_numexp);


                if (ctrl.transport_only == KIN_REACTION)
                {
                    // In reaction mode, simulate reaction for soil, and speciation for stream
                    Reaction(kstep, 86400.0, chemtbl, kintbl, &rttbl, &subcatch_numexp);
                }
                else
                {
                    Speciation(chemtbl, ctrl, &rttbl, &subcatch_numexp);
                }
                
                PrintDailyResults(file_pointer, ctrl.transport_only, steps_numexp[kstep], &rttbl, &subcatch_numexp);
                ResetReactionRates(&subcatch_numexp);

            }
            biort_printf(VL_NORMAL, "\nHBV-BioRT %s numerical experiment succeeded.\n", dir);

            fclose(file_pointer);
        }
    }

    FreeStruct(&steps, &subcatch);
    FreeStruct(&steps_numexp, &subcatch_numexp);
}
