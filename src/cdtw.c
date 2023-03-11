/*
    This code is written by Davide Albanese <davide.albanese@gmail.com>.
    (C) mlpy Developers.

    This program is free software: you can redistribute it and/or modify
    it underthe terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "cdtw.h"

float
min3(float a, float b, float c)
{
  float min;

  min = a;
  if (b < min)
    min = b;
  if (c < min)
    min = c;
  return min;
}

// Paliwal adjustment window used for restricting the warping function
// r: window length
int
paliwal_window(int i, int j, int n, int m, int r)
{
  float s, f;

  s = ((float) m) / n;
  f = fabs(i - (((float) j) / s));

  if (f <= r)
    return 1;
  else
    return 0;
}

// euclidean distance
float e_dist(float x, float y)
{
  return fabs(x - y);
}

// squared euclidean distance
float se_dist(float x, float y)
{
  return pow(x - y, 2);
}


// Returns the (unnormalized) minimum-distance warp path
// between time series x and y and the cost matrix C
float
std_dtw(float *x, float *y, int64_t n, int64_t m, float *cost, int64_t squared)
{
  int i, j;
  float (*dist)(float, float);

  if (squared == 0)
    dist = &e_dist;
  else
    dist = &se_dist;

  cost[0] = (*dist)(x[0], y[0]);

  for (i=1; i<n; i++)
    cost[i*m] = (*dist)(x[i], y[0]) + cost[(i-1)*m];

  for (j=1; j<m; j++)
    cost[j] = (*dist)(x[0], y[j]) + cost[(j-1)];

  for (i=1; i<n; i++)
    for (j=1; j<m; j++)
      cost[i*m+j] = (*dist)(x[i], y[j]) +
	min3(cost[(i-1)*m+j], cost[(i-1)*m+(j-1)], cost[i*m+(j-1)]);

  return cost[n*m-1];
}

// Compute the warp path starting at cost[startx, starty]
// If startx = -1 -> startx = n-1; if starty = -1 -> starty = m-1
int
path(float *cost, int64_t n, int64_t m, int64_t startx, int64_t starty, Path *p)
{
  int i, j, k, z1, z2;
  int *px;
  int *py;
  float min_cost;

  if ((startx >= n) || (starty >= m))
    return 0;

  if (startx < 0)
    startx = n - 1;

  if (starty < 0)
    starty = m - 1;

  i = startx;
  j = starty;
  k = 1;

  // allocate path for the worst case
  px = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));
  py = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));

  px[0] = i;
  py[0] = j;

  while ((i > 0) || (j > 0))
    {
      if (i == 0)
	j--;
      else if (j == 0)
	i--;
      else
	{
	  min_cost = min3(cost[(i-1)*m+j],
			  cost[(i-1)*m+(j-1)],
			  cost[i*m+(j-1)]);

	  if (cost[(i-1)*m+(j-1)] == min_cost)
	    {
	      i--;
	      j--;
	    }
	  else if (cost[i*m+(j-1)] == min_cost)
	    j--;
	  else
	    i--;
	}

      px[k] = i;
      py[k] = j;
      k++;
    }

  p->px = (int *) malloc (k * sizeof(int));
  p->py = (int *) malloc (k * sizeof(int));
  for (z1=0, z2=k-1; z1<k; z1++, z2--)
    {
      p->px[z1] = px[z2];
      p->py[z1] = py[z2];
    }
  p->k = k;

  free(px);
  free(py);

  return 1;
}


//
void
subsequence(float *x, float *y, int64_t n, int64_t m, float *cost)
{
  int i, j;

  cost[0] = fabs(x[0]-y[0]);

  for (i=1; i<n; i++)
    cost[i*m] = fabs(x[i]-y[0]) + cost[(i-1)*m];

  for (j=1; j<m; j++)
    cost[j] = fabs(x[0]-y[j]); // subsequence variation: D(0,j) := c(x0, yj)

  for (i=1; i<n; i++)
    for (j=1; j<m; j++)
      cost[i*m+j] = fabs(x[i]-y[j]) +
	min3(cost[(i-1)*m+j], cost[(i-1)*m+(j-1)], cost[i*m+(j-1)]);

}


int
subsequence_path(float *cost, int64_t n, int64_t m, int64_t starty, Path *p)
{
  int i, z1, z2;
  int a_star;
  int *tmpx, *tmpy;

  // find path
  if (!path(cost, n, m, -1, starty, p))
    return 0;

  // find a_star
  a_star = 0;
  for (i=1; i<p->k; i++)
    if (p->px[i] == 0)
      a_star++;
    else
      break;

  // rebuild path
  tmpx = p->px;
  tmpy = p->py;
  p->px = (int *) malloc ((p->k-a_star) * sizeof(int));
  p->py = (int *) malloc ((p->k-a_star) * sizeof(int));
  for (z1=0, z2=a_star; z2<p->k; z1++, z2++)
    {
      p->px[z1] = tmpx[z2];
      p->py[z1] = tmpy[z2];
    }
  p->k = p->k-a_star;

  free(tmpx);
  free(tmpy);

  return 1;
}



// ======================================================================

COST_DTYPE hw_pe(SIG_DTYPE x, SIG_DTYPE y, COST_DTYPE n, COST_DTYPE m, COST_DTYPE nw) {
    COST_DTYPE min = n;
    if (m < min) min = m;
    if (nw < min) min = nw;

    COST_DTYPE cost;
    if (x > y) cost = (COST_DTYPE) (x - y);
    else cost = (COST_DTYPE) (y - x);

    return cost + min;
}

void hw_subsequence(float *x, float *y, int n, int m, float *cost)
{
    SIG_DTYPE *scaled_x = (SIG_DTYPE *) malloc(n * sizeof(SIG_DTYPE));
    SIG_DTYPE *scaled_y = (SIG_DTYPE *) malloc(m * sizeof(SIG_DTYPE));
    COST_DTYPE *scaled_cost = (COST_DTYPE *) malloc(n * m * sizeof(COST_DTYPE));

    int i, j;
    for (i = 0; i < n; i++) {
        scaled_x[i] = (SIG_DTYPE) (x[i] * SCALING);
    }

    for (j = 0; j < m; j++) {
        scaled_y[j] = (SIG_DTYPE) (y[j] * SCALING);
    }

    // ===
    _sw_sdtw(scaled_x, scaled_y, n, m, scaled_cost);

    // ===

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            cost[i * m + j] = ((float) scaled_cost[i * m + j]) / SCALING;
        }
    }
}

void _hw_sdtw(SIG_DTYPE *scaled_x, SIG_DTYPE *scaled_y, int64_t n, int64_t m, COST_DTYPE *scaled_cost) {
    int i, j;

    // First row
    for (j = 0; j < m; j++) {
        scaled_cost[j] = hw_pe(scaled_x[0], scaled_y[j], 0, 0, 0);
        // printf("%d ", scaled_cost[j]);
    }

    // First column
    for (i = 1; i < n; i++) {
        scaled_cost[i*m] = hw_pe(scaled_x[i], scaled_y[0], scaled_cost[(i - 1)*m], COST_DTYPE_MAX, COST_DTYPE_MAX);
    }

    // Rest of the matrix
    for (i = 1; i < n; i++) {
        for (j = 1; j < m; j++) {
            scaled_cost[(i * m) + j] = hw_pe(scaled_x[i], scaled_y[j], scaled_cost[(i - 1) * m + j], scaled_cost[(i * m) + j - 1], scaled_cost[(i - 1) * m  + j - 1]);
        }
    }
}

// ===

COST_DTYPE _cost(SIG_DTYPE x, SIG_DTYPE y) {
    COST_DTYPE cost;
    if (x > y) cost = (COST_DTYPE) (x - y);
    else cost = (COST_DTYPE) (y - x);
    return cost;
}

COST_DTYPE _min(COST_DTYPE a, COST_DTYPE b, COST_DTYPE c) {
    COST_DTYPE min = a;
    if (b < min) min = b;
    if (c < min) min = c;
    return min;
}

void _sw_sdtw(SIG_DTYPE *scaled_x, SIG_DTYPE *scaled_y, int64_t n, int64_t m, COST_DTYPE *scaled_cost) {
    int i, j;
    scaled_cost[0] = _cost(scaled_x[0], scaled_y[0]);
    for (i = 1; i < n; i++) {
        scaled_cost[i*m] = _cost(scaled_x[i], scaled_y[0]) + scaled_cost[(i - 1)*m];
    }

    for (j = 1; j < m; j++) {
        scaled_cost[j] = _cost(scaled_x[0], scaled_y[j]);
    }

    for (i = 1; i < n; i++) {
        for (j = 1; j < m; j++) {
            scaled_cost[(i * m) + j] = _cost(scaled_x[i], scaled_y[j])
                + _min(scaled_cost[(i - 1) * m + j], scaled_cost[(i * m) + j - 1], scaled_cost[(i - 1) * m  + j - 1]);
        }
    }
}

int
path_DS(COST_DTYPE *cost, int64_t n, int64_t m, int64_t startx, int64_t starty, Path *p)
{
  int i, j, k, z1, z2;
  int *px;
  int *py;
  COST_DTYPE min_cost;

  if ((startx >= n) || (starty >= m))
    return 0;

  if (startx < 0)
    startx = n - 1;

  if (starty < 0)
    starty = m - 1;

  i = startx;
  j = starty;
  k = 1;

  // allocate path for the worst case
  px = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));
  py = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));

  px[0] = i;
  py[0] = j;

  while ((i > 0) || (j > 0))
    {
      if (i == 0)
	j--;
      else if (j == 0)
	i--;
      else
	{
    // min3
    min_cost = cost[(i-1)*m+j];
    if (cost[(i-1)*m+(j-1)] < min_cost) min_cost = cost[(i-1)*m+(j-1)];
    if (cost[i*m+(j-1)] < min_cost) min_cost = cost[i*m+(j-1)];

	  if (cost[(i-1)*m+(j-1)] == min_cost)
	    {
	      i--;
	      j--;
	    }
	  else if (cost[i*m+(j-1)] == min_cost)
	    j--;
	  else
	    i--;
	}

      px[k] = i;
      py[k] = j;
      k++;
    }

  p->px = (int *) malloc (k * sizeof(int));
  p->py = (int *) malloc (k * sizeof(int));
  for (z1=0, z2=k-1; z1<k; z1++, z2--)
    {
      p->px[z1] = px[z2];
      p->py[z1] = py[z2];
    }
  p->k = k;

  free(px);
  free(py);

  return 1;
}

int
subsequence_path_DS(COST_DTYPE *cost, int64_t n, int64_t m, int64_t starty, Path *p)
{
  int i, z1, z2;
  int a_star;
  int *tmpx, *tmpy;

  // find path
  if (!path_DS(cost, n, m, -1, starty, p))
    return 0;

  // find a_star
  a_star = 0;
  for (i=1; i<p->k; i++)
    if (p->px[i] == 0)
      a_star++;
    else
      break;

  // rebuild path
  tmpx = p->px;
  tmpy = p->py;
  p->px = (int *) malloc ((p->k-a_star) * sizeof(int));
  p->py = (int *) malloc ((p->k-a_star) * sizeof(int));
  for (z1=0, z2=a_star; z2<p->k; z1++, z2++)
    {
      p->px[z1] = tmpx[z2];
      p->py[z1] = tmpy[z2];
    }
  p->k = p->k-a_star;

  free(tmpx);
  free(tmpy);

  return 1;
}