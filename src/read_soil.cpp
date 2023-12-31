#include "biort.hpp"

void ReadSoil(const string& dir, Subcatchment& subcatch)
{
    char            cmdstr[MAXSTRING];
    char            file_name[MAXSTRING];
    int             line_number = 0;
    FILE           *file_pointer;

    // READ SOIL.TXT FILE
    sprintf(file_name, "input/%s/soil.txt", dir.c_str());
    file_pointer = fopen(file_name, "r");

    biort_printf(VL_NORMAL, "SOIL PARAMTERS\n");

    //===== Porosity =====//
    NextLine(file_pointer, cmdstr, &line_number);   // 2021-05-14
    subcatch.soil_surface.porosity = ReadParamToDouble(cmdstr, "POROSITY_SURFACE", file_name, line_number);
    biort_printf(VL_NORMAL, "  Surface zone porosity is %.2f m3 m-3.\n", subcatch.soil_surface.porosity);

    NextLine(file_pointer, cmdstr, &line_number);
    subcatch.soil_sz.porosity = ReadParamToDouble(cmdstr, "POROSITY_UZ", file_name, line_number);
    biort_printf(VL_NORMAL, "  Upper zone porosity is %.2f m3 m-3.\n", subcatch.soil_sz.porosity);

    NextLine(file_pointer, cmdstr, &line_number);
    subcatch.soil_dz.porosity = ReadParamToDouble(cmdstr, "POROSITY_LZ", file_name, line_number);
    biort_printf(VL_NORMAL, "  Lower zone porosity is %.2f m3 m-3.\n", subcatch.soil_dz.porosity);

    //===== Passive water storage =====//
    NextLine(file_pointer, cmdstr, &line_number);
    subcatch.soil_surface.ws_passive = ReadParamToDouble(cmdstr, "WS_PASSIVE_SURFACE", file_name, line_number);
    subcatch.soil_surface.ws_passive = std::max(subcatch.soil_surface.ws_passive, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Surface zone passive water storage is %.2f mm.\n", subcatch.soil_surface.ws_passive);

    NextLine(file_pointer, cmdstr, &line_number);
    subcatch.soil_sz.ws_passive = ReadParamToDouble(cmdstr, "WS_PASSIVE_UZ", file_name, line_number);
    subcatch.soil_sz.ws_passive = std::max(subcatch.soil_sz.ws_passive, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Upper zone passive water storage is %.2f mm.\n", subcatch.soil_sz.ws_passive);

    NextLine(file_pointer, cmdstr, &line_number);
    subcatch.soil_dz.ws_passive = ReadParamToDouble(cmdstr, "WS_PASSIVE_LZ", file_name, line_number);
    subcatch.soil_dz.ws_passive = std::max(subcatch.soil_dz.ws_passive, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Lower zone passive water storage is %.2f mm.\n", subcatch.soil_dz.ws_passive);

    //===== Depth =====//
    NextLine(file_pointer, cmdstr, &line_number);
    subcatch.soil_surface.depth = ReadParamToDouble(cmdstr, "DEPTH_SURFACE", file_name, line_number);
    biort_printf(VL_NORMAL, "  Depth of surface zone is %.2f mm.\n", subcatch.soil_surface.depth);

    NextLine(file_pointer, cmdstr, &line_number);
    subcatch.soil_sz.depth = ReadParamToDouble(cmdstr, "DEPTH_UZ", file_name, line_number);
    biort_printf(VL_NORMAL, "  Depth of upper zone is %.2f mm.\n", subcatch.soil_sz.depth);

    NextLine(file_pointer, cmdstr, &line_number);
    subcatch.soil_dz.depth = ReadParamToDouble(cmdstr, "DEPTH_LZ", file_name, line_number);
    biort_printf(VL_NORMAL, "  Depth of lower zone is %.2f mm.\n", subcatch.soil_dz.depth);

    fclose(file_pointer);
}
