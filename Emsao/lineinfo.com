# File Emsao/lineinfo.com
# April 23, 2008

common/lineinfo/ velline, errline, eqwline, htline, widline, eqeline

double	velline[MAXREF], errline[MAXREF], eqwline[MAXREF]
double	htline[MAXREF], widline[MAXREF], eqeline[MAXREF]
