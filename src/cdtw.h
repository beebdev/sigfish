#include <stdlib.h>
#include <stdint.h>
#include "shared.h"

typedef struct Path
{
  int k;
  int *px;
  int *py;
} Path;


float std_dtw(float *x, float *y, int n, int m, float *cost, int squared);
int path(float *cost, int n, int m, int startx, int starty, Path *p);
void subsequence(float *x, float *y, int n, int m, float *cost);
int subsequence_path(float *cost, int n, int m, int starty, Path *p);

void _hw_sdtw(SIG_DTYPE *scaled_x, SIG_DTYPE *scaled_y, int n, int m, COST_DTYPE *scaled_cost);
void _sw_sdtw(SIG_DTYPE *scaled_x, SIG_DTYPE *scaled_y, int n, int m, COST_DTYPE *scaled_cost);