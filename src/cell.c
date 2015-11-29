/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "trconfig.h"
#include "transsys.h"


int add_neighbor(CELL *cell1, CELL *cell2, double weight)
{
  CELL **ca;
  double *w;

  if (cell1 == cell2)
  {
    fprintf(stderr, "add_neighbor: cells are identical\n");
    return (-1);
  }
  ca = (CELL **) realloc(cell1->neighbor, (cell1->num_neighbors + 1) * sizeof(CELL *));
  if (ca == NULL)
    return (-1);
  else
    cell1->neighbor = ca;
  w = realloc(cell1->contact_weight, (cell1->num_neighbors + 1) * sizeof(double));
  if (w == NULL)
    return (-1);
  else
    cell1->contact_weight = w;
  cell1->num_neighbors++;
  ca = (CELL **) realloc(cell2->neighbor, (cell2->num_neighbors + 1) * sizeof(CELL *));
  if (ca == NULL)
    return (-1);
  else
    cell2->neighbor = ca;
  w = realloc(cell2->contact_weight, (cell2->num_neighbors + 1) * sizeof(double));
  if (w == NULL)
    return (-1);
  else
    cell2->contact_weight = w;
  cell2->num_neighbors++;
  cell1->neighbor[cell1->num_neighbors - 1] = cell2;
  cell1->contact_weight[cell1->num_neighbors - 1] = weight;
  cell2->neighbor[cell2->num_neighbors - 1] = cell1;
  cell2->contact_weight[cell2->num_neighbors - 1] = weight;
  return (0);
}


int set_contact_weight(CELL *cell1, CELL *cell2, double weight)
{
  int n1, n2;

  for (n1 = 0; n1 < cell1->num_neighbors; n1++)
  {
    if (cell1->neighbor[n1] == cell2)
      break;
  }
  for (n2 = 0; n2 < cell2->num_neighbors; n2++)
  {
    if (cell2->neighbor[n2] == cell1)
      break;
  }
  if ((n1 < cell1->num_neighbors) && (n2 < cell2->num_neighbors))
  {
    cell1->contact_weight[n1] = weight;
    cell2->contact_weight[n2] = weight;
    return (0);
  }
  else if ((n1 == cell1->num_neighbors) && (n2 == cell2->num_neighbors))
    return (add_neighbor(cell1, cell2, weight));
  else
  {
    fprintf(stderr, "set_contact_weight: inconsistent contacts\n");
    return (-1);
  }
}


int cell_expression(CELL *cell)
{
  int return_value;

  if (!cell->existing || !cell->alive)
    return (0);
  return_value = process_expression(&(cell->transsys_instance));
  return (return_value);
}


int diffuse(int num_cells, CELL *cell)
{
  int cn, i, n;
  double d, d1, diff;
  const TRANSSYS_INSTANCE *ti;

  /* fprintf(stderr, "diffuse: starting\n"); */
  for (cn = 0; cn < num_cells; cn++)
  {
    if (cell[cn].existing)
    {
      for (i = 0; i < cell[cn].transsys_instance.transsys->num_factors; i++)
      {
	cell[cn].transsys_instance.new_concentration[i] = cell[cn].transsys_instance.factor_concentration[i];
      }
    }
  }
  /* fprintf(stderr, "diffuse: initialised new_concentration\n"); */
  for (cn = 0; cn < num_cells; cn++)
  {
    /* fprintf(stderr, "diffuse: processing cell %d\n", cn); */
    if (cell[cn].existing)
    {
      /* fprintf(stderr, "diffuse: processing existing cell %d\n", cn); */
      for (i = 0; i < cell[cn].transsys_instance.transsys->num_factors; i++)
      {
	d = cell[cn].transsys_instance.factor_concentration[i] / (cell[cn].num_neighbors + 1);
	for (n = 0; n < cell[cn].num_neighbors; n++)
	{
          /* fprintf(stderr, "diffuse: processing cell %d, neighbour %d of %d\n", cn, n, cell[cn].num_neighbors); */
	  ti = &(cell[cn].transsys_instance);
          /* fprintf(stderr, "diffuse:   got transsys instance\n"); */
          if (ti->transsys->factor_list[i].diffusibility_expression == NULL)
          {
            diff = 0.0;
            fprintf(stderr, "diffuse: transsys %s, factor %s has no diffusibility expression, assuming %f\n", ti->transsys->name, ti->transsys->factor_list[i].name, diff);
          }
          else
          {
            diff = evaluate_expression(ti->transsys->factor_list[i].diffusibility_expression, &ti);
          }
          /* fprintf(stderr, "diffuse:   evaluated diffusibility, %f\n", diff); */
	  d1 = d * diff * cell[cn].contact_weight[n];
	  ti->new_concentration[i] -= d1;
          /* fprintf(stderr, "diffuse:   updated cell's new_concentration\n"); */
	  cell[cn].neighbor[n]->transsys_instance.new_concentration[i] += d1;
          /* fprintf(stderr, "diffuse:   updated neighbour's new_concentration\n"); */
	}
      }
    }
  }
  /* fprintf(stderr, "diffuse: updated new_concentration\n"); */
  for (cn = 0; cn < num_cells; cn++)
  {
    if (cell[cn].existing)
    {
      for (i = 0; i < cell[cn].transsys_instance.transsys->num_factors; i++)
      {
	cell[cn].transsys_instance.factor_concentration[i] = cell[cn].transsys_instance.new_concentration[i];
      }
    }
  }
  /* fprintf(stderr, "diffuse: finishing\n"); */
  return (0);
}


double total_concentration(int num_cells, const CELL *cell, int factor_no)
{
  int cn;
  double s;

  if (factor_no >= cell->transsys_instance.transsys->num_factors)
  {
    fprintf(stderr, "total_concentration: factor number %d out of range (%d)\n",
            factor_no, cell->transsys_instance.transsys->num_factors);
    return (0.0);
  }
  s = 0.0;
  for (cn = 0; cn < num_cells; cn++)
  {
    if (cell[cn].existing)
      s += cell[cn].transsys_instance.factor_concentration[factor_no];
  }
  return (s);
}


double max_concentration(int num_cells, const CELL *cell, int factor_no)
{
  int cn;
  double s;

  if (factor_no >= cell->transsys_instance.transsys->num_factors)
  {
    fprintf(stderr, "max_concentration: factor number %d out of range (%d)\n",
            factor_no, cell->transsys_instance.transsys->num_factors);
    return (0.0);
  }
  s = 0.0;
  for (cn = 0; cn < num_cells; cn++)
  {
    if (cell[cn].existing && (cell[cn].transsys_instance.factor_concentration[factor_no] > s))
      s = cell[cn].transsys_instance.factor_concentration[factor_no];
  }
  return (s);
}
