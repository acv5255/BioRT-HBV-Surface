#include "biort.h"

void ReadSoil(const char dir[], subcatch_struct* sub)
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
    ReadParam(cmdstr, "POROSITY_UZ", 'd', file_name, line_number, &sub->porosity_uz);
    biort_printf(VL_NORMAL, "  Upper zone porosity is %.2f m3 m-3.\n", sub->porosity_uz);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "POROSITY_LZ", 'd', file_name, line_number, &sub->porosity_lz);
    biort_printf(VL_NORMAL, "  Lower zone porosity is %.2f m3 m-3.\n", sub->porosity_lz);

    //NextLine(fp, cmdstr, &lno);  // 2021-05-14
    //ReadParam(cmdstr, "WS_PASSIVE_SURFACE", 'd', fn, lno, &sub[0].res_surface);
    //sub[0].res_surface = MAX(sub[0].res_surface, STORAGE_MIN);
    //biort_printf(VL_NORMAL, "  Surface passive water storage is %.2f mm.\n", sub[0].res_surface);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "WS_PASSIVE_UZ", 'd', file_name, line_number, &sub->res_uz);
    sub->res_uz = MAX(sub->res_uz, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Upper zone passive water storage is %.2f mm.\n", sub->res_uz);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "WS_PASSIVE_LZ", 'd', file_name, line_number, &sub->res_lz);
    sub->res_lz = MAX(sub->res_lz, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Lower zone passive water storage is %.2f mm.\n", sub->res_lz);

    //NextLine(fp, cmdstr, &lno);  // 2021-05-14
    //ReadParam(cmdstr, "DEPTH_SURFACE", 'd', fn, lno, &sub[0].d_surface);
    //biort_printf(VL_NORMAL, "  Depth of surface zone is %.2f mm.\n", sub[0].d_surface);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "DEPTH_UZ", 'd', file_name, line_number, &sub->d_uz);
    biort_printf(VL_NORMAL, "  Depth of upper zone is %.2f mm.\n", sub->d_uz);

    NextLine(file_pointer, cmdstr, &line_number);
    ReadParam(cmdstr, "DEPTH_LZ", 'd', file_name, line_number, &sub->d_lz);
    biort_printf(VL_NORMAL, "  Depth of lower zone is %.2f mm.\n", sub->d_lz);

    fclose(file_pointer);
}
