/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.2  2002/01/25 03:35:03  kim
 * Added gnuplot link functionality to transexpr, transscatter, improved
 *     PostScript stuff
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "transsys.h"


static int postscript_header(FILE *f, const char *title, int eps, const POSTSCRIPT_BOX *bounding_box)
{
  if (eps)
    fprintf(f, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  else
    fprintf(f, "%%!PS-Adobe-3.0\n");
  fprintf(f, "%%%%Creator: %s\n", prgname);
  if (title)
    fprintf(f, "%%%%Title: %s\n", title);
  fprintf(f, "%%%%BoundingBox: %f %f %f %f\n",
	  bounding_box->x, bounding_box->y, bounding_box->x + bounding_box->w, bounding_box->y + bounding_box->h);
  fprintf(f, "%%%%EndComments\n");
  return (0);
}


static int postscript_page_start(FILE *f, const POSTSCRIPT_STYLE *style)
{
  fprintf(f, "%%%%Page: %d %d\n", style->page_num, style->page_num);
  fprintf(f, "save\n");
  fprintf(f, "/%s findfont %f scalefont setfont\n", style->fontname, style->fontheight);
}


static int postscript_page_end(FILE *f, const POSTSCRIPT_STYLE *style)
{
  fprintf(f, "restore\n");
  fprintf(f, "showpage\n");
}


static int procset_transsys_arcplot(FILE *f)
{
  fprintf(f, "%%%%BeginProcSet: transsys_arcplot\n");
  fprintf(f, "/B\n");
  fprintf(f, "{\n");
  fprintf(f, "  newpath moveto 0.0 rlineto 0.0 exch rlineto 0.0 rlineto closepath stroke\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "/S\n");
  fprintf(f, "{\n");
  fprintf(f, "  moveto show\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "/SC\n");
  fprintf(f, "{\n");
  fprintf(f, "  moveto show\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "/AH\n");
  fprintf(f, "{\n");
  fprintf(f, "  newpath moveto rlineto 0.0 exch rlineto closepath stroke\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "/AHF\n");
  fprintf(f, "{\n");
  fprintf(f, "  newpath moveto rlineto 0.0 exch rlineto closepath fill\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "/AC\n");
  fprintf(f, "{\n");
  fprintf(f, "  newpath moveto curveto stroke\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "%%%%EndProcSet\n");
  return (0);
}


static int procset_ltranssys_primitives(FILE *f)
{
  fprintf(f, "%%%%BeginProcSet: ltranssys_primitives\n");
  fprintf(f, "/fillstroke { gsave gsave fill grestore 0 setgray stroke grestore } bind def\n");
  fprintf(f, "/MOVE { 0 exch translate } bind def\n");
  fprintf(f, "/PUSH { save } bind def\n");
  fprintf(f, "/POP { restore } bind def\n");
  fprintf(f, "/SPHERE\n");
  fprintf(f, "{\n");
  fprintf(f, "  newpath 0 exch 0 exch 0 360 arc fillstroke\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "/CYLINDER\n");
  fprintf(f, "{\n");
  fprintf(f, "  newpath dup -0.5 mul exch 3 1 roll exch dup -0.5 mul exch 4 1 roll\n");
  fprintf(f, "  moveto dup 3 1 roll 0 rlineto 0 exch rlineto -1.0 mul 0 rlineto closepath\n");
  fprintf(f, "  fillstroke\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "/BOX\n");
  fprintf(f, "{\n");
  fprintf(f, "  newpath dup -0.5 mul exch 3 1 roll exch dup -0.5 mul exch 4 1 roll\n");
  fprintf(f, "  moveto dup 3 1 roll 0 rlineto 0 exch rlineto -1.0 mul 0 rlineto closepath\n");
  fprintf(f, "  fillstroke pop\n");
  fprintf(f, "} bind def\n");
  fprintf(f, "/TURN { rotate } bind def\n");
  fprintf(f, "/ROLL { pop } bind def\n");
  fprintf(f, "/BANK { pop } bind def\n");
  fprintf(f, "/COLOR { setrgbcolor } bind def\n");
  fprintf(f, "%%%%EndProcSet\n");
  return (0);
}


static int postscript_trailer(FILE *f, const POSTSCRIPT_STYLE *style)
{
  fprintf(f, "%%%%EOF\n");
  return (0);
}


static int postscript_graphics_primitive(FILE *f, const GRAPHICS_PRIMITIVE *gp, const TRANSSYS_INSTANCE *transsys_instance, double scale)
{
  double arg[GRAPHICS_PRIMITIVE_ARGMAX];
  int i;

  for (i = 0; i < gp->num_arguments; i++)
    arg[i] = evaluate_expression(gp->argument[i], &transsys_instance);
  switch (gp->type)
  {
  case GRAPHICS_PUSH:
    fprintf(f, "PUSH\n");
    break;
  case GRAPHICS_POP:
    fprintf(f, "POP\n");
    break;
  case GRAPHICS_MOVE:
    fprintf(f, "%e MOVE\n", arg[0] * scale);
    break;
  case GRAPHICS_TURN:
    fprintf(f, "%e TURN\n", arg[0]);
    break;
  case GRAPHICS_ROLL:
    fprintf(f, "%e ROLL\n", arg[0]);
    break;
  case GRAPHICS_BANK:
    fprintf(f, "%e BANK\n", arg[0]);
    break;
  case GRAPHICS_SPHERE:
    fprintf(f, "%e SPHERE\n", arg[0] * scale);
    break;
  case GRAPHICS_CYLINDER:
    fprintf(f, "%e %e CYLINDER\n", arg[0] * scale, arg[1] * scale);
    break;
  case GRAPHICS_BOX:
    fprintf(f, "%e %e %e BOX\n", arg[0] * scale, arg[1] * scale, arg[2] * scale);
    break;
  case GRAPHICS_COLOR:
    fprintf(f, "%e %e %e COLOR\n", arg[0], arg[1], arg[2]);
    break;
  default:
    fprintf(stderr, "postscript_graphics_primitive: primitive type %d not implemented\n", (int) gp->type);
    return (-1);
  }
  return (0);
}


int postscript_symbol_string(FILE *f, const SYMBOL_INSTANCE *symbol_string, const POSTSCRIPT_STYLE *style, double scale)
{
  const SYMBOL_INSTANCE *symbol;
  const GRAPHICS_PRIMITIVE *gp;
  int return_value = 0;

  fprintf(f, "%%%%Page: %d %d\n", style->page_num, style->page_num);
  fprintf(f, "%% symbol_strlen: %lu\n", (unsigned long) symbol_strlen(symbol_string));
  fprintf(f, "save\n");
  fprintf(f, "%e setlinewidth\n", style->linewidth);
  fprintf(f, "%e %e translate\n", style->lnd_x0, style->lnd_y0);
  for (symbol = symbol_string; symbol; symbol = symbol->next)
  {
    for (gp = symbol->lsys->symbol_list[symbol->symbol_index].graphics_primitive_list; gp; gp = gp->next)
    {
      return_value = postscript_graphics_primitive(f, gp, &(symbol->transsys_instance), scale);
      if (return_value)
	break;
    }
  }
  fprintf(f, "restore\n");
  fprintf(f, "showpage\n\n");
  return (return_value);
}


int postscript_lsys_prolog(FILE *f, const LSYS *lsys, const POSTSCRIPT_STYLE *style)
{
  postscript_header(f, lsys->name, style->eps, &(style->pagebox));
  procset_ltranssys_primitives(f);
  fprintf(f, "%%%%EndProlog\n");
  return (0);
}


static size_t max_factorname_length(const TRANSSYS *transsys)
{
  const FACTOR_ELEMENT *fe;
  size_t m, l;

  m = 0;
  for (fe = transsys->factor_list; fe; fe = fe->next)
  {
    l = strlen(fe->name);
    if (l > m)
      m = l;
  }
  return (m);
}


static size_t max_genename_length(const TRANSSYS *transsys)
{
  const GENE_ELEMENT *ge;
  size_t m, l;

  m = 0;
  for (ge = transsys->gene_list; ge; ge = ge->next)
  {
    l = strlen(ge->name);
    if (l > m)
      m = l;
  }
  return (m);
}


/*
static int draw_arrow(FILE *f, const SQL_FIELD *field1, const SQL_FIELD *field2, const PS_STYLE *style)
{
  double y1, y2, xb, xt, xt_min;

  y1 = field1->y + 0.5 * field1->h;
  y2 = field2->y + 0.5 * field2->h;
  xb = style->box_x - style->arr_h;
  xt_min = style->box_x - 3.0 * style->arr_h;
  if (xt_min < style->x0)
    xt_min = style->x0;
  xt = xt_min - (xt_min - style->x0) * (fabs(y2 - y1) / style->height);
  fprintf(f, "%f %f %f %f %f AHV\n", style->arr_w * 2, -style->arr_h, -style->arr_w, style->box_x, y2);
  fprintf(f, "%f %f %f %f %f %f %f %f AC\n", xt, y1, xt, y2, xb, y2, style->box_x, y1);
  return (0);
}
*/

static int draw_activation_arrows(FILE *f, const TRANSSYS *transsys, const PROMOTER_ELEMENT *a, double xp, double yp, const POSTSCRIPT_STYLE *style)
{
  const GENE_ELEMENT *ge;
  int i;
  double y1, xt, xt_min;

  for (ge = transsys->gene_list; ge; ge = ge->next)
  {
    for (i = 0; i < a->num_binding_factors; i++)
    {
      if (a->factor_index[i] == ge->product_index)
      {
	fprintf(f, "%% activation arrow from %s\n", ge->name);
	y1 = ge->box.y + 1.0 * style->fontheight;
	xt_min = ge->box.x - 3.0 * style->arrow_length;
	if (xt_min < style->pagebox.x)
	  xt_min = style->pagebox.x;
	xt = xt_min - (xt_min - style->pagebox.x) * (fabs(yp - y1) / style->pagebox.h);
	/* fprintf(f, "%% xt_min = %f, xt = %f, y1 = %f, xp = %f, yp = %f\n", xt_min, xt, y1, xp, yp); */
	fprintf(f, "%f %f %f %f %f %f %f %f AC\n", xt, y1, xt, yp, xp, yp, ge->box.x, y1);
      }
    }
  }
  return (0);
}


static int draw_repression_arrows(FILE *f, const TRANSSYS *transsys, const PROMOTER_ELEMENT *a, double xp, double yp, const POSTSCRIPT_STYLE *style)
{
  const GENE_ELEMENT *ge;
  int i;
  double y1, xt, xt_min;

  for (ge = transsys->gene_list; ge; ge = ge->next)
  {
    fprintf(f, "%% checking gene %s (product index %d)\n", ge->name, ge->product_index);
    for (i = 0; i < a->num_binding_factors; i++)
    {
      fprintf(f, "%% ... %d\n", a->factor_index[i]);
      if (a->factor_index[i] == ge->product_index)
      {
	fprintf(f, "%% repression arrow from %s\n", ge->name);
	y1 = ge->box.y + 1.0 * style->fontheight;
	xt_min = ge->box.x + ge->box.w + 3.0 * style->arrow_length;
	if (xt_min > style->pagebox.x + style->pagebox.w)
	  xt_min = style->pagebox.x + style->pagebox.w;
	xt = xt_min + (style->pagebox.x + style->pagebox.w - xt_min) * (fabs(yp - y1) / style->pagebox.h);
	/* fprintf(f, "%% xt_min = %f, xt = %f, y1 = %f, xp = %f, yp = %f\n", xt_min, xt, y1, xp, yp); */
	fprintf(f, "%f %f %f %f %f %f %f %f AC\n", xt, y1, xt, yp, xp, yp, ge->box.x + ge->box.w, y1);
      }
    }
  }
  return (0);
}


int postscript_gene(FILE *f, const TRANSSYS *transsys, const GENE_ELEMENT *gene, const POSTSCRIPT_STYLE *style)
{
  double y;
  const PROMOTER_ELEMENT *a;

  fprintf(f, "%f setlinewidth\n", style->box_linewidth);
  fprintf(f, "%f %f %f %f %f B\n", -gene->box.w, gene->box.h, gene->box.w, gene->box.x, gene->box.y);
  y = gene->box.y + gene->box.h - style->fontheight * 1.2;
  /* fprintf(f, "(%s) %f %f SC\n", gene->name, x, y); */
  fprintf(f, "(%s) %f %f SC\n", gene->name, gene->box.x + 0.7 * style->fontheight, y);
  y = gene->box.y + 0.8 * style->fontheight;
  /* fprintf(f, "(%s) %f %f SC\n", transsys->factor_list[gene->product_index].name, x, y); */
  fprintf(f, "(%s) %f %f SC\n", transsys->factor_list[gene->product_index].name, gene->box.x + 0.7 * style->fontheight, y);
  y = gene->box.y + gene->box.h - style->fontheight * 0.5;
  fprintf(f, "%f setlinewidth\n", style->arrow_linewidth);
  for (a = gene->promoter_list; a; a = a->next)
  {
    switch (a->type)
    {
    case ACT_CONSTITUTIVE:
      break;
    case ACT_ACTIVATE:
      fprintf(f, "%f %f %f %f %f AHF\n", style->arrow_width, -style->arrow_length, -style->arrow_width * 0.5, gene->box.x, y);
      fprintf(f, "%% activation arrows to %s\n", gene->name);
      draw_activation_arrows(f, transsys, a, gene->box.x - style->arrow_length, y, style);
      break;
    case ACT_REPRESS:
      fprintf(f, "%f %f %f %f %f AH\n", style->arrow_width, style->arrow_length, -style->arrow_width * 0.5, gene->box.x + gene->box.w, y);
      fprintf(f, "%% repression arrows to %s\n", gene->name);
      draw_repression_arrows(f, transsys, a, gene->box.x + gene->box.w + style->arrow_length, y, style);
      break;
    default:
      fprintf(stderr, "postscript_gene: cannot draw activation type %d\n", (int) a->type);
      break;
    }
    y -= 2.0 * style->arrow_width;
  }
  return (0);
}


int set_genebox_coordinates(const TRANSSYS *transsys, const POSTSCRIPT_STYLE *style)
{
  GENE_ELEMENT *ge;
  const PROMOTER_ELEMENT *a;
  double x0, y0, h;
  double boxwidth;
  size_t cw, cwg;
  int n;

  cw = max_factorname_length(transsys);
  cwg = max_genename_length(transsys);
  if (cw < cwg)
    cw = cwg;
  boxwidth = (cw + 2) * style->fontheight * 0.7;
  x0 = style->pagebox.x + (style->pagebox.w - boxwidth) * 0.5;
  y0 = style->pagebox.y + style->pagebox.h - 2.0 * style->fontheight;
  for (ge = transsys->gene_list; ge; ge = ge->next)
  {
    ge->box.x = x0;
    ge->box.w = boxwidth;
    n = 0;
    for (a = ge->promoter_list; a; a = a->next)
      n++;
    h = n * 2 * style->arrow_width;
    if (h > style->fontheight)
      ge->box.h = 2.0 * style->fontheight + h;
    else
      ge->box.h = 3.0 * style->fontheight;
    ge->box.y = y0 - ge->box.h;
    y0 = ge->box.y - style->fontheight;
  }
  return (0);
}


int postscript_transsys_init(FILE *f, const char *title, const POSTSCRIPT_STYLE *style)
{
  postscript_header(f, title, style->eps, (&style->pagebox));
  procset_transsys_arcplot(f);
  fprintf(f, "%%%%EndProlog\n");
  fprintf(f, "%%%%BeginSetup\n");
  fprintf(f, "%%%%PaperSize: A4\n");
  fprintf(f, "/%s findfont %f scalefont setfont\n", style->fontname, style->fontheight);
  fprintf(f, "%%%%EndSetup\n");
  return (0);
}


int postscript_transsys_finish(FILE *f, const POSTSCRIPT_STYLE *style)
{
  postscript_trailer(f, style);
}


int postscript_transsys(FILE *f, const TRANSSYS *transsys, const POSTSCRIPT_STYLE *style)
{
  GENE_ELEMENT *ge;

  postscript_page_start(f, style);
  set_genebox_coordinates(transsys, style);
  fprintf(f, "%f (%s) stringwidth pop 0.5 mul sub %f moveto\n", style->pagebox.x + 0.5 * style->pagebox.w, transsys->name, style->pagebox.y + style->pagebox.h - style->fontheight);
  fprintf(f, "(%s) show\n", transsys->name);
  for (ge = transsys->gene_list; ge; ge = ge->next)
  {
    postscript_gene(f, transsys, ge, style);
  }
  postscript_page_end(f, style);
  return (0);
}


int set_default_postscript_style(int eps, const TRANSSYS *transsys, POSTSCRIPT_STYLE *style)
{
  GENE_ELEMENT *ge;

  if (eps)
  {
    if (transsys == NULL)
    {
      fprintf(stderr, "set_default_postscript_style: transsys needed for default eps style\n");
      return (-1);
    }
    style->eps = 1;
    style->pagebox.x = 0.0;
    style->pagebox.y = 0.0;
    style->pagebox.w = 600.0;
    style->pagebox.h = 850.0;
    style->box_linewidth = 0.5;
    style->fontheight = 8.0;
    style->arrow_width = 3.0;
    style->arrow_length = 9.0;
    style->arrow_linewidth = 0.3;
    strncpy(style->fontname, "Courier-Bold", IDENTIFIER_MAX);
    set_genebox_coordinates(transsys, style);
    for (ge = transsys->gene_list; ge; ge = ge->next)
    {
      if (ge->next == NULL)
	style->pagebox.h -= ge->box.y;
    }
  }
  else
  {
    style->eps = 0;
    style->pagebox.x = 50.0;
    style->pagebox.y = 50.0;
    style->pagebox.w = 500.0;
    style->pagebox.h = 750.0;
    style->page_num = 0;
    style->lnd_x0 = style->pagebox.x + 0.5 * style->pagebox.w;
    style->lnd_y0 = style->pagebox.y + 0.2 * style->pagebox.h;
    style->linewidth = 1.0;
    style->box_linewidth = 0.5;
    style->fontheight = 8.0;
    style->arrow_width = 3.0;
    style->arrow_length = 9.0;
    style->arrow_linewidth = 0.3;
    strncpy(style->fontname, "Courier-Bold", IDENTIFIER_MAX);
  }
  return (0);
}

