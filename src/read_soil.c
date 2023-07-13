#include "biort.h"

void ReadSoil(const char dir[], Subcatchment* sub)
{
    char            cmdstr[MAXSTRING];
    char            file_name[MAXSTRING];
    int             line_number = 0;
    FILE           *file_pointer;

    // READ SOIL.TXT FILE
    sprintf(file_name, "input/%s/soil.txt", dir);
    file_pointer = fopen(file_name, "r");

    biort_printf(VL_NORMAL, "SOIL PARAMTERS\n");

    //NextLine(fp, cmdstr, &lno);   // 2021-05-14
    //ReadParam(cmdstr, "POROSITY_SURFACE", 'd', fn, lno, &sub[0].porosity_surface);
    //biort_printf(VL_NORMAL, "  Surface porosity is %.2f m3 m-3.\n", sub[0].porosity_surface);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "POROSITY_UZ", 'd', file_name, line_number, &sub->soil_sz.porosity);
    biort_printf(VL_NORMAL, "  Upper zone porosity is %.2f m3 m-3.\n", sub->soil_sz.porosity);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "POROSITY_LZ", 'd', file_name, line_number, &sub->soil_dz.porosity);
    biort_printf(VL_NORMAL, "  Lower zone porosity is %.2f m3 m-3.\n", sub->soil_dz.porosity);

    //NextLine(fp, cmdstr, &lno);  // 2021-05-14
    //ReadParam(cmdstr, "WS_PASSIVE_SURFACE", 'd', fn, lno, &sub[0].res_surface);
    //sub[0].res_surface = MAX(sub[0].res_surface, STORAGE_MIN);
    //biort_printf(VL_NORMAL, "  Surface passive water storage is %.2f mm.\n", sub[0].res_surface);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "WS_PASSIVE_UZ", 'd', file_name, line_number, &sub->soil_sz.ws_passive);
    sub->soil_sz.ws_passive = MAX(sub->soil_sz.ws_passive, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Upper zone passive water storage is %.2f mm.\n", sub->soil_sz.ws_passive);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "WS_PASSIVE_LZ", 'd', file_name, line_number, &sub->soil_dz.ws_passive);
    sub->soil_dz.ws_passive = MAX(sub->soil_dz.ws_passive, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Lower zone passive water storage is %.2f mm.\n", sub->soil_dz.ws_passive);

    //NextLine(fp, cmdstr, &lno);  // 2021-05-14
    //ReadParam(cmdstr, "DEPTH_SURFACE", 'd', fn, lno, &sub[0].d_surface);
    //biort_printf(VL_NORMAL, "  Depth of surface zone is %.2f mm.\n", sub[0].d_surface);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "DEPTH_UZ", 'd', file_name, line_number, &sub->soil_sz.depth);
    biort_printf(VL_NORMAL, "  Depth of upper zone is %.2f mm.\n", sub->soil_sz.depth);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "DEPTH_LZ", 'd', file_name, line_number, &sub->soil_dz.depth);
    biort_printf(VL_NORMAL, "  Depth of lower zone is %.2f mm.\n", sub->soil_dz.depth);

    fclose(file_pointer);
}
