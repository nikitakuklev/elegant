#include "../track.h"

void sdds_strength_output(char *outputfile, LINE_LIST *beamline, char *input);
void convert_to_patpet(char *output, LINE_LIST *beamline, long flip_k, double angle_tolerance,
                    char *header_file, char *ender_file);
void convert_to_parmela(char *output, LINE_LIST *beamline, long flip_k, double angle_tolerance, char *header_file, char *ender_file,
    double quad_aperture, double sext_aperture, double pc);
void convert_to_patricia(char *output, LINE_LIST *beamline, long flip_k, double angle_tolerance,
                    char *header_file, char *ender_file);
void convert_to_transport(char *output, LINE_LIST *beamline, long flip_k, double angle_tolerance, char *header_file, char *ender_file,
    double quad_aperture, double sext_aperture, double pc);

void do_output_transport(FILE *fp, char *s);
void section(char *s, int n);
long nearly_equal(double x1, double x2, double frac);
char *quoted_label(char *s);
void do_output(FILE *fp, char *s);
