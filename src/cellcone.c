#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <GL/glut.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "transsys.h"


#define NO_FACTOR -1


typedef struct
{
  double x, y, z;
} ANGLES;

typedef struct
{
  int show_extrusion;
  int solid_flag;
  int lighting;
  int extrusion_box_flag;
  int xmin, xmax, ymin, ymax, zmin, zmax;
  int red_factor, green_factor, blue_factor;
  double red_max, green_max, blue_max;
  int current_axis, current_factor, current_minmax;
  double background;
  double boxsize;
  double extrusion_boxheight;
  ANGLES angles, spin_angles;
} DISPLAY_PARAMETERS;


typedef struct
{
  int num_cells;
  int xwidth, ywidth, zwidth;
  TRANSSYS *transsys;
  CELL ***cell;
  int extrusion_length;
  int extrusion_index, extrusion_planes;
  double ****extrusion;
} MERISTEM_CONE;



static MERISTEM_CONE cone;
static DISPLAY_PARAMETERS display_parameters;
static int window_width, window_height;
static double boxsize = 0.8;
static const int xoff[6] = {-1,  0,  0,  1,  0,  0};
static const int yoff[6] = { 0, -1,  0,  0,  1,  0};
static const int zoff[6] = { 0,  0, -1,  0,  0,  1};
static int time_running = 0, trace_flag = 0;
static unsigned long time_step = 0;
static double rendering_time = 0.0;
static long num_renderings = 0;

void display_eps(void);

/* cannibalized from Mark Kilgard's rendereps.c */

typedef struct _Feedback3Dcolor {
  GLfloat x;
  GLfloat y;
  GLfloat z;
  GLfloat red;
  GLfloat green;
  GLfloat blue;
  GLfloat alpha;
} Feedback3Dcolor;


/* Write contents of one vertex to stdout. */
void
print3DcolorVertex(GLint size, GLint * count,
  GLfloat * buffer)
{
  int i;

  printf("  ");
  for (i = 0; i < 7; i++) {
    printf("%4.2f ", buffer[size - (*count)]);
    *count = *count - 1;
  }
  printf("\n");
}

void
printBuffer(GLint size, GLfloat * buffer)
{
  GLint count;
  int token, nvertices;

  count = size;
  while (count) {
    token = buffer[size - count];
    count--;
    switch (token) {
    case GL_PASS_THROUGH_TOKEN:
      printf("GL_PASS_THROUGH_TOKEN\n");
      printf("  %4.2f\n", buffer[size - count]);
      count--;
      break;
    case GL_POINT_TOKEN:
      printf("GL_POINT_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_LINE_TOKEN:
      printf("GL_LINE_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_LINE_RESET_TOKEN:
      printf("GL_LINE_RESET_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_POLYGON_TOKEN:
      printf("GL_POLYGON_TOKEN\n");
      nvertices = buffer[size - count];
      count--;
      for (; nvertices > 0; nvertices--) {
        print3DcolorVertex(size, &count, buffer);
      }
    }
  }
}

GLfloat pointSize;

static char *gouraudtriangleEPS[] =
{
  "/bd{bind def}bind def /triangle { aload pop   setrgbcolor  aload pop 5 3",
  "roll 4 2 roll 3 2 roll exch moveto lineto lineto closepath fill } bd",
  "/computediff1 { 2 copy sub abs threshold ge {pop pop pop true} { exch 2",
  "index sub abs threshold ge { pop pop true} { sub abs threshold ge } ifelse",
  "} ifelse } bd /computediff3 { 3 copy 0 get 3 1 roll 0 get 3 1 roll 0 get",
  "computediff1 {true} { 3 copy 1 get 3 1 roll 1 get 3 1 roll 1 get",
  "computediff1 {true} { 3 copy 2 get 3 1 roll  2 get 3 1 roll 2 get",
  "computediff1 } ifelse } ifelse } bd /middlecolor { aload pop 4 -1 roll",
  "aload pop 4 -1 roll add 2 div 5 1 roll 3 -1 roll add 2 div 3 1 roll add 2",
  "div 3 1 roll exch 3 array astore } bd /gouraudtriangle { computediff3 { 4",
  "-1 roll aload 7 1 roll 6 -1 roll pop 3 -1 roll pop add 2 div 3 1 roll add",
  "2 div exch 3 -1 roll aload 7 1 roll exch pop 4 -1 roll pop add 2 div 3 1",
  "roll add 2 div exch 3 -1 roll aload 7 1 roll pop 3 -1 roll pop add 2 div 3",
  "1 roll add 2 div exch 7 3 roll 10 -3 roll dup 3 index middlecolor 4 1 roll",
  "2 copy middlecolor 4 1 roll 3 copy pop middlecolor 4 1 roll 13 -1 roll",
  "aload pop 17 index 6 index 15 index 19 index 6 index 17 index 6 array",
  "astore 10 index 10 index 14 index gouraudtriangle 17 index 5 index 17",
  "index 19 index 5 index 19 index 6 array astore 10 index 9 index 13 index",
  "gouraudtriangle 13 index 16 index 5 index 15 index 18 index 5 index 6",
  "array astore 12 index 12 index 9 index gouraudtriangle 17 index 16 index",
  "15 index 19 index 18 index 17 index 6 array astore 10 index 12 index 14",
  "index gouraudtriangle 18 {pop} repeat } { aload pop 5 3 roll aload pop 7 3",
  "roll aload pop 9 3 roll 4 index 6 index 4 index add add 3 div 10 1 roll 7",
  "index 5 index 3 index add add 3 div 10 1 roll 6 index 4 index 2 index add",
  "add 3 div 10 1 roll 9 {pop} repeat 3 array astore triangle } ifelse } bd",
  NULL
};

GLfloat *
spewPrimitiveEPS(FILE * file, GLfloat * loc)
{
  int token;
  int nvertices, i;
  GLfloat red, green, blue;
  int smooth;
  GLfloat dx, dy, dr, dg, db, absR, absG, absB, colormax;
  int steps;
  Feedback3Dcolor *vertex;
  GLfloat xstep, ystep, rstep, gstep, bstep;
  GLfloat xnext = 0.0, ynext = 0.0, rnext = 0.0, gnext = 0.0, bnext = 0.0, distance = 0.0;

  token = *loc;
  loc++;
  switch (token) {
  case GL_LINE_RESET_TOKEN:
  case GL_LINE_TOKEN:
    vertex = (Feedback3Dcolor *) loc;

    dr = vertex[1].red - vertex[0].red;
    dg = vertex[1].green - vertex[0].green;
    db = vertex[1].blue - vertex[0].blue;

    if (dr != 0 || dg != 0 || db != 0) {
      /* Smooth shaded line. */
      dx = vertex[1].x - vertex[0].x;
      dy = vertex[1].y - vertex[0].y;
      distance = sqrt(dx * dx + dy * dy);

      absR = fabs(dr);
      absG = fabs(dg);
      absB = fabs(db);

#define Max(a,b) (((a)>(b))?(a):(b))

#define EPS_SMOOTH_LINE_FACTOR 0.06  /* Lower for better smooth 

                                        lines. */

      colormax = Max(absR, Max(absG, absB));
      steps = Max(1.0, colormax * distance * EPS_SMOOTH_LINE_FACTOR);

      xstep = dx / steps;
      ystep = dy / steps;

      rstep = dr / steps;
      gstep = dg / steps;
      bstep = db / steps;

      xnext = vertex[0].x;
      ynext = vertex[0].y;
      rnext = vertex[0].red;
      gnext = vertex[0].green;
      bnext = vertex[0].blue;

      /* Back up half a step; we want the end points to be
         exactly the their endpoint colors. */
      xnext -= xstep / 2.0;
      ynext -= ystep / 2.0;
      rnext -= rstep / 2.0;
      gnext -= gstep / 2.0;
      bnext -= bstep / 2.0;
    } else {
      /* Single color line. */
      steps = 0;
    }

    fprintf(file, "%g %g %g setrgbcolor\n",
      vertex[0].red, vertex[0].green, vertex[0].blue);
    fprintf(file, "%g %g moveto\n", vertex[0].x, vertex[0].y);

    for (i = 0; i < steps; i++) {
      xnext += xstep;
      ynext += ystep;
      rnext += rstep;
      gnext += gstep;
      bnext += bstep;
      fprintf(file, "%g %g lineto stroke\n", xnext, ynext);
      fprintf(file, "%g %g %g setrgbcolor\n", rnext, gnext, bnext);
      fprintf(file, "%g %g moveto\n", xnext, ynext);
    }
    fprintf(file, "%g %g lineto stroke\n", vertex[1].x, vertex[1].y);

    loc += 14;          /* Each vertex element in the feedback
                           buffer is 7 GLfloats. */

    break;
  case GL_POLYGON_TOKEN:
    nvertices = *loc;
    loc++;

    vertex = (Feedback3Dcolor *) loc;

    if (nvertices > 0) {
      red = vertex[0].red;
      green = vertex[0].green;
      blue = vertex[0].blue;
      smooth = 0;
      for (i = 1; i < nvertices; i++) {
        if (red != vertex[i].red || green != vertex[i].green || blue != vertex[i].blue) {
          smooth = 1;
          break;
        }
      }
      if (smooth) {
        /* Smooth shaded polygon; varying colors at vetices. */
        int triOffset;

        /* Break polygon into "nvertices-2" triangle fans. */
        for (i = 0; i < nvertices - 2; i++) {
          triOffset = i * 7;
          fprintf(file, "[%g %g %g %g %g %g]",
            vertex[0].x, vertex[i + 1].x, vertex[i + 2].x,
            vertex[0].y, vertex[i + 1].y, vertex[i + 2].y);
          fprintf(file, " [%g %g %g] [%g %g %g] [%g %g %g] gouraudtriangle\n",
            vertex[0].red, vertex[0].green, vertex[0].blue,
            vertex[i + 1].red, vertex[i + 1].green, vertex[i + 1].blue,
            vertex[i + 2].red, vertex[i + 2].green, vertex[i + 2].blue);
        }
      } else {
        /* Flat shaded polygon; all vertex colors the same. */
        fprintf(file, "newpath\n");
        fprintf(file, "%g %g %g setrgbcolor\n", red, green, blue);

        /* Draw a filled triangle. */
        fprintf(file, "%g %g moveto\n", vertex[0].x, vertex[0].y);
        for (i = 1; i < nvertices; i++) {
          fprintf(file, "%g %g lineto\n", vertex[i].x, vertex[i].y);
        }
        fprintf(file, "closepath fill\n\n");
      }
    }
    loc += nvertices * 7;  /* Each vertex element in the
                              feedback buffer is 7 GLfloats. */
    break;
  case GL_POINT_TOKEN:
    vertex = (Feedback3Dcolor *) loc;
    fprintf(file, "%g %g %g setrgbcolor\n", vertex[0].red, vertex[0].green, vertex[0].blue);
    fprintf(file, "%g %g %g 0 360 arc fill\n\n", vertex[0].x, vertex[0].y, pointSize / 2.0);
    loc += 7;           /* Each vertex element in the feedback
                           buffer is 7 GLfloats. */
    break;
  default:
    /* XXX Left as an excersie to the reader. */
    printf("Incomplete implementation.  Unexpected token (%d).\n", token);
    exit(1);
  }
  return loc;
}

void
spewUnsortedFeedback(FILE * file, GLint size, GLfloat * buffer)
{
  GLfloat *loc, *end;

  loc = buffer;
  end = buffer + size;
  while (loc < end) {
    loc = spewPrimitiveEPS(file, loc);
  }
}

typedef struct _DepthIndex {
  GLfloat *ptr;
  GLfloat depth;
} DepthIndex;

static int
compare(const void *a, const void *b)
{
  DepthIndex *p1 = (DepthIndex *) a;
  DepthIndex *p2 = (DepthIndex *) b;
  GLfloat diff = p2->depth - p1->depth;

  if (diff > 0.0) {
    return 1;
  } else if (diff < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

void
spewSortedFeedback(FILE * file, GLint size, GLfloat * buffer)
{
  int token;
  GLfloat *loc, *end;
  Feedback3Dcolor *vertex;
  GLfloat depthSum;
  int nprimitives, item;
  DepthIndex *prims;
  int nvertices, i;

  end = buffer + size;

  /* Count how many primitives there are. */
  nprimitives = 0;
  loc = buffer;
  while (loc < end) {
    token = *loc;
    loc++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      loc += 14;
      nprimitives++;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = *loc;
      loc++;
      loc += (7 * nvertices);
      nprimitives++;
      break;
    case GL_POINT_TOKEN:
      loc += 7;
      nprimitives++;
      break;
    default:
      /* XXX Left as an excersie to the reader. */
      printf("Incomplete implementation.  Unexpected token (%d).\n",
        token);
      exit(1);
    }
  }

  /* Allocate an array of pointers that will point back at
     primitives in the feedback buffer.  There will be one
     entry per primitive.  This array is also where we keep the
     primitive's average depth.  There is one entry per
     primitive  in the feedback buffer. */
  prims = (DepthIndex *) malloc(sizeof(DepthIndex) * nprimitives);

  item = 0;
  loc = buffer;
  while (loc < end) {
    prims[item].ptr = loc;  /* Save this primitive's location. */
    token = *loc;
    loc++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      vertex = (Feedback3Dcolor *) loc;
      depthSum = vertex[0].z + vertex[1].z;
      prims[item].depth = depthSum / 2.0;
      loc += 14;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = *loc;
      loc++;
      vertex = (Feedback3Dcolor *) loc;
      depthSum = vertex[0].z;
      for (i = 1; i < nvertices; i++) {
        depthSum += vertex[i].z;
      }
      prims[item].depth = depthSum / nvertices;
      loc += (7 * nvertices);
      break;
    case GL_POINT_TOKEN:
      vertex = (Feedback3Dcolor *) loc;
      prims[item].depth = vertex[0].z;
      loc += 7;
      break;
    default:
      /* XXX Left as an excersie to the reader. */
      assert(1);
    }
    item++;
  }
  assert(item == nprimitives);

  /* Sort the primitives back to front. */
  qsort(prims, nprimitives, sizeof(DepthIndex), compare);

  /* Understand that sorting by a primitives average depth
     doesn't allow us to disambiguate some cases like self
     intersecting polygons.  Handling these cases would require
     breaking up the primitives.  That's too involved for this
     example.  Sorting by depth is good enough for lots of
     applications. */

  /* Emit the Encapsulated PostScript for the primitives in
     back to front order. */
  for (item = 0; item < nprimitives; item++) {
    (void) spewPrimitiveEPS(file, prims[item].ptr);
  }

  free(prims);
}

#define EPS_GOURAUD_THRESHOLD 0.1  /* Lower for better (slower) 

                                      smooth shading. */

void
spewWireFrameEPS(FILE * file, int doSort, GLint size, GLfloat * buffer, char *creator)
{
  GLfloat clearColor[4], viewport[4];
  GLfloat lineWidth;
  int i;

  /* Read back a bunch of OpenGL state to help make the EPS
     consistent with the OpenGL clear color, line width, point
     size, and viewport. */
  glGetFloatv(GL_VIEWPORT, viewport);
  glGetFloatv(GL_COLOR_CLEAR_VALUE, clearColor);
  glGetFloatv(GL_LINE_WIDTH, &lineWidth);
  glGetFloatv(GL_POINT_SIZE, &pointSize);

  /* Emit EPS header. */
  fputs("%!PS-Adobe-2.0 EPSF-2.0\n", file);
  /* Notice %% for a single % in the fprintf calls. */
  fprintf(file, "%%%%Creator: %s (using OpenGL feedback)\n", creator);
  fprintf(file, "%%%%BoundingBox: %g %g %g %g\n",
    viewport[0], viewport[1], viewport[2], viewport[3]);
  fputs("%%EndComments\n", file);
  fputs("\n", file);
  fputs("gsave\n", file);
  fputs("\n", file);

  /* Output Frederic Delhoume's "gouraudtriangle" PostScript
     fragment. */
  fputs("% the gouraudtriangle PostScript fragement below is free\n", file);
  fputs("% written by Frederic Delhoume (delhoume@ilog.fr)\n", file);
  fprintf(file, "/threshold %g def\n", EPS_GOURAUD_THRESHOLD);
  for (i = 0; gouraudtriangleEPS[i]; i++) {
    fprintf(file, "%s\n", gouraudtriangleEPS[i]);
  }

  fprintf(file, "\n%g setlinewidth\n", lineWidth);

  /* Clear the background like OpenGL had it. */
  fprintf(file, "%g %g %g setrgbcolor\n",
    clearColor[0], clearColor[1], clearColor[2]);
  fprintf(file, "%g %g %g %g rectfill\n\n",
    viewport[0], viewport[1], viewport[2], viewport[3]);

  if (doSort) {
    spewSortedFeedback(file, size, buffer);
  } else {
    spewUnsortedFeedback(file, size, buffer);
  }

  /* Emit EPS trailer. */
  fputs("grestore\n\n", file);
  fputs("%Add `showpage' to the end of this file to be able to print to a printer.\n",
    file);

  fclose(file);
}

void
outputEPS(size_t size, int doSort, char *filename)
{
  GLfloat *feedbackBuffer;
  GLint returned;
  FILE *file;

  feedbackBuffer = calloc(size, sizeof(GLfloat));
  glFeedbackBuffer(size, GL_3D_COLOR, feedbackBuffer);
  (void) glRenderMode(GL_FEEDBACK);
  display_eps();
  returned = glRenderMode(GL_RENDER);
  fprintf(stderr, "outputEPS: allocated %lu, returned: %d\n", (unsigned long) size, returned);
  if (filename) {
    file = fopen(filename, "w");
    if (file) {
      spewWireFrameEPS(file, doSort, returned, feedbackBuffer, "rendereps");
    } else {
      printf("Could not open %s\n", filename);
    }
  } else {
    /* Helps debugging to be able to see the decode feedback
       buffer as text. */
    printBuffer(returned, feedbackBuffer);
  }
  free(feedbackBuffer);
}

/* end cannibalized section */


/* cloned from ltransgl.c... GL stuff should really go to some lib... */

static void render_wireframe_box(double x, double y, double z)
{
  fprintf(stderr, "render_wireframe_box: currently not implemented\n");
}


static void render_solid_box(double x, double y, double z)
{
  glBegin(GL_QUADS);
  glNormal3f(0.0, -1.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(x, 0.0, 0.0);
  glVertex3f(x, 0.0, z);
  glVertex3f(0.0, 0.0, z);
  glNormal3f(0.0, 1.0, 0.0);
  glVertex3f(0.0, y, 0.0);
  glVertex3f(0.0, y, z);
  glVertex3f(x, y, z);
  glVertex3f(x, y, 0.0);
  glNormal3f(0.0, 0.0, -1.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, y, 0.0);
  glVertex3f(x, y, 0.0);
  glVertex3f(x, 0.0, 0.0);
  glNormal3f(0.0, 0.0, 1.0);
  glVertex3f(0.0, 0.0, z);
  glVertex3f(x, 0.0, z);
  glVertex3f(x, y, z);
  glVertex3f(0.0, y, z);
  glNormal3f(-1.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, z);
  glVertex3f(0.0, y, z);
  glVertex3f(0.0, y, 0.0);
  glNormal3f(1.0, 0.0, 0.0);
  glVertex3f(x, 0.0, 0.0);
  glVertex3f(x, y, 0.0);
  glVertex3f(x, y, z);
  glVertex3f(x, 0.0, z);
  glEnd();
/*   glBegin(GL_QUAD_STRIP); */
/*   glVertex3f(0.0, 0.0, z); */
/*   glVertex3f(0.0, y, z); */
/*   glVertex3f(x, 0.0, z); */
/*   glVertex3f(x, y, z); */
/*   glVertex3f(x, 0.0, 0.0); */
/*   glVertex3f(x, y, 0.0); */
/*   glVertex3f(0.0, 0.0, 0.0); */
/*   glVertex3f(0.0, y, 0.0); */
/*   glVertex3f(0.0, 0.0, z); */
/*   glVertex3f(0.0, y, z); */
/*   glEnd(); */
}


size_t cone_polygon_complexity(const MERISTEM_CONE *cone, const DISPLAY_PARAMETERS *dp)
{
  int x, y, z;
  size_t c;

  c = 0;
  for (x = 0; x < cone->xwidth; x++)
  {
    for (y = 0; y < cone->ywidth; y++)
    {
      for (z = 0; z < cone->zwidth; z++)
      {
	if (cone->cell[x][y][z].existing)
	  c += 276;
      }
    }
  }
  if (dp->show_extrusion)
  {
    for (x = 0; x < cone->xwidth; x++)
    {
      for (y = 0; y < cone->extrusion_length; y++)
      {
	for (z = 0; z < cone->zwidth; z++)
	{
	  if (cone->extrusion[x][y][z][0] != -1.0)
	  {
	    if (dp->extrusion_box_flag)
	      c += 276;
	    else
	      c += 46;
	  }
	}
      }
    }
  }
  return (c);
}


void xz_boxf(GLfloat xw, GLfloat zw, int polygon_flag)
{
  static GLuint xz_boxlist = 0;
  static GLfloat xz_xwidth = -1.0, xz_zwidth = -1.0;
  GLfloat x, z;

  if (xz_boxlist == 0)
  {
    xz_boxlist = glGenLists(1);
    if (xz_boxlist == 0)
    {
      fprintf(stderr, "xz_boxf: failed to get display list\n");
      return;
    }
  }
  if ((xw != xz_xwidth) || (zw != xz_zwidth))
  {
    x = xw * 0.5;
    z = zw * 0.5;
    glNewList(xz_boxlist, GL_COMPILE);
    glNormal3f(0.0, -1.0, 0.0);
    glVertex3f(-x, 0.0, -z);
    glVertex3f(x, 0.0, -z);
    glVertex3f(x, 0.0, z);
    glVertex3f(-x, 0.0, z);
    glEnd();
    glEndList();
    xz_xwidth = xw;
    xz_zwidth = zw;
  }
  if (polygon_flag)
    glBegin(GL_POLYGON);
  else
    glBegin(GL_LINE_LOOP);
  glCallList(xz_boxlist);
}


void free_double4(double ****d)
{
    free(d[0][0][0]);
    free(d[0][0]);
    free(d[0]);
    free(d);
}


double **** alloc_double4(int d1, int d2, int d3, int d4, double v)
{
  double ****d;
  int i1, i2, i3, i4;

  d = (double ****) malloc(d1 * sizeof(double ***));
  if (d == NULL)
    return (NULL);
  d[0] = (double ***) malloc(d1 * d2 * sizeof(double **));
  if (d[0] == NULL)
  {
    free(d);
    return (NULL);
  }
  d[0][0] = (double **) malloc(d1 * d2 * d3 * sizeof(double *));
  if (d[0][0] == NULL)
  {
    free(d[0]);
    free(d);
    return (NULL);
  }
  d[0][0][0] = (double *) malloc(d1 * d2 * d3 * d4 * sizeof(double));
  if (d[0][0][0] == NULL)
  {
    free(d[0][0]);
    free(d[0]);
    free(d);
    return (NULL);
  }
  for (i1 = 0; i1 < d1; i1++)
  {
    d[i1] = d[0] + i1 * d2;
    d[i1][0] = d[0][0] + i1 * d2 * d3;
    d[i1][0][0] = d[0][0][0] + i1 * d2 * d3 * d4;
    for (i2 = 0; i2 < d2; i2++)
    {
      d[i1][i2] = d[i1][0] + i2 * d3;
      d[i1][i2][0] = d[i1][0][0] + i2 * d3 * d4;
      for (i3 = 0; i3 < d3; i3++)
      {
	d[i1][i2][i3] = d[i1][i2][0] + i3 * d4;
	for (i4 = 0; i4 < d4; i4++)
	  d[i1][i2][i3][i4] = v;
      }
    }
  }
  return (d);
}


void init_cone(MERISTEM_CONE *cone)
{
  cone->cell = NULL;
  cone->transsys = NULL;
  cone->xwidth = 0;
  cone->ywidth = 0;
  cone->zwidth = 0;
  cone->num_cells = 0;
}


void free_cone(MERISTEM_CONE *cone)
{
  if (cone->cell)
  {
    if (cone->cell[0])
    {
      if (cone->cell[0][0])
	free_cells(cone->num_cells, cone->cell[0][0]);
    }
    free(cone->cell);
  }
  init_cone(cone);
}


int create_cone(int xwidth, int ywidth, int zwidth, TRANSSYS *transsys, int extrusion_length, MERISTEM_CONE *cone)
{
  int x, y, z, xn, yn, zn, n;
  double q, d, xx, yy, zz;

  cone->extrusion = alloc_double4(xwidth, extrusion_length, zwidth, transsys->num_factors, -1.0);
  if (cone->extrusion == NULL)
    return (-1);
  cone->extrusion_length = extrusion_length;
  cone->extrusion_index = 0;
  cone->extrusion_planes = 0;
  cone->cell = (CELL ***) malloc(xwidth * sizeof(CELL **));
  if (cone->cell == NULL)
  {
    free_double4(cone->extrusion);
    return (-1);
  }
  cone->cell[0] = (CELL **) malloc(xwidth * ywidth * sizeof(CELL *));
  if (cone->cell[0] == NULL)
  {
    free_double4(cone->extrusion);
    free(cone->cell);
    return (-1);
  }
  for (x = 1; x < xwidth; x++)
    cone->cell[x] = cone->cell[0] + x * ywidth;
  cone->cell[0][0] = new_cells(xwidth * ywidth * zwidth, transsys);
  if (cone->cell[0][0] == NULL)
  {
    free_double4(cone->extrusion);
    free(cone->cell[0]);
    free(cone->cell);
    return (-1);
  }
  for (x = 1; x < xwidth; x++)
    cone->cell[x][0] = cone->cell[0][0] + x * ywidth * zwidth;
  for (x = 0; x < xwidth; x++)
  {
    for (y = 1; y < ywidth; y++)
      cone->cell[x][y] = cone->cell[x][0] + y * zwidth;
  }
  cone->transsys = transsys;
  cone->num_cells = xwidth * ywidth * zwidth;
  cone->xwidth = xwidth;
  cone->ywidth = ywidth;
  cone->zwidth = zwidth;
  for (y = 0; y < ywidth; y++)
  {
    yy = ((double) y) / (ywidth - 1);
    q = sqrt(1.0 - yy * yy) * 0.5;
    for (x = 0; x < xwidth; x++)
    {
      xx = ((double) x) / (xwidth - 1) - 0.5;
      for (z = 0; z < zwidth; z++)
      {
	zz = ((double) z) / (zwidth - 1) - 0.5;
        d = sqrt(xx * xx + zz * zz);
	if (d <= q)
	{
	  cone->cell[x][y][z].existing = 1;
	  cone->cell[x][y][z].alive = 1;
	}
      }
    }
  }
  for (y = 0; y < ywidth; y++)
  {
    for (x = 0; x < xwidth; x++)
    {
      for (z = 0; z < zwidth; z++)
      {
	for (n = 0; n < 6; n++)
	{
	  xn = x + xoff[n];
	  yn = y + yoff[n];
	  zn = z + zoff[n];
	  if ((0 <= xn) &&  (xn < xwidth) && (0 <= yn) && (yn < ywidth) && (0 <= zn) && (zn < zwidth)
		  && cone->cell[xn][yn][zn].existing)
	  {
	    set_contact_weight(&(cone->cell[x][y][z]), &(cone->cell[xn][yn][zn]), 1.0 / 6.0);
	  }
        }
      }
    }
  }
  return (0);
}


void printf_cone(const MERISTEM_CONE *cone)
{
  int x, y, z;

  for (z = 0; z < cone->zwidth; z++)
  {
    for (y = 0; y < cone->ywidth; y++)
    {
      for (x = 0; x < cone->xwidth; x++)
      {
        if (cone->cell[x][y][z].existing)
	  printf("*");
	else
	  printf(" ");
      }
      printf("\n");
    }
    printf("\n");
  }
}


int extend_extrusion(MERISTEM_CONE *cone)
{
  int x, z, i;

  if (cone->extrusion_planes < cone->extrusion_length)
    cone->extrusion_planes++;
  cone->extrusion_index = (cone->extrusion_index + 1) % cone->extrusion_planes;
  for (x = 0; x < cone->xwidth; x++)
  {
    for (z = 0; z < cone->zwidth; z++)
    {
      if (cone->cell[x][0][z].alive)
      {
	for (i = 0; i < cone->transsys->num_factors; i++)
	  cone->extrusion[x][cone->extrusion_index][z][i] = cone->cell[x][0][z].transsys_instance.factor_concentration[i];
      }
    }
  }
  return (0);
}


void display_model(void)
{
  int i, x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;
  double x0, y0, z0, r, g, b, rmax, gmax, bmax;
  GLfloat v[4];

  x0 = (cone.xwidth - 1) * 0.5;
  y0 = (cone.ywidth - 1) * 0.5;
  z0 = (cone.zwidth - 1) * 0.5;
  if (display_parameters.red_factor != NO_FACTOR)
    rmax = max_concentration(cone.num_cells, cone.cell[0][0], display_parameters.red_factor);
  if (display_parameters.green_factor != NO_FACTOR)
    gmax = max_concentration(cone.num_cells, cone.cell[0][0], display_parameters.green_factor);
  if (display_parameters.blue_factor != NO_FACTOR)
    bmax = max_concentration(cone.num_cells, cone.cell[0][0], display_parameters.blue_factor);
  if (display_parameters.red_max < rmax)
    display_parameters.red_max = rmax;
  if (display_parameters.green_max < gmax)
    display_parameters.green_max = gmax;
  if (display_parameters.blue_max < bmax)
    display_parameters.blue_max = bmax;
  xmin = display_parameters.xmin > 0 ? display_parameters.xmin : 0;
  xmax = display_parameters.xmax < cone.xwidth ? display_parameters.xmax : cone.xwidth;
  ymin = display_parameters.ymin > 0 ? display_parameters.ymin : 0;
  ymax = display_parameters.ymax < cone.ywidth ? display_parameters.ymax : cone.ywidth;
  zmin = display_parameters.zmin > 0 ? display_parameters.zmin : 0;
  zmax = display_parameters.zmax < cone.zwidth ? display_parameters.zmax : cone.zwidth;
  glRotatef(display_parameters.angles.x, 1.0, 0.0, 0.0);
  glRotatef(display_parameters.angles.z, 0.0, 0.0, 1.0);
  glRotatef(display_parameters.angles.y, 0.0, 1.0, 0.0);
  glTranslatef(-x0, -y0, -z0);
  for (x = xmin; x < xmax; x++)
  {
    for (y = ymin; y < ymax; y++)
    {
      for (z = zmin; z < zmax; z++)
      {
	if (cone.cell[x][y][z].existing)
	{
	  r = display_parameters.background;
	  g = display_parameters.background;
	  b = display_parameters.background;
	  if ((display_parameters.red_factor != NO_FACTOR) && (display_parameters.red_max > 0.0))
	    r += cone.cell[x][y][z].transsys_instance.factor_concentration[display_parameters.red_factor] / display_parameters.red_max * (1.0 - display_parameters.background);
	  if ((display_parameters.green_factor != NO_FACTOR) && (display_parameters.green_max > 0.0))
	    g += cone.cell[x][y][z].transsys_instance.factor_concentration[display_parameters.green_factor] / display_parameters.green_max * (1.0 - display_parameters.background);
	  if ((display_parameters.blue_factor != NO_FACTOR) && (display_parameters.blue_max > 0.0))
	    b += cone.cell[x][y][z].transsys_instance.factor_concentration[display_parameters.blue_factor] / display_parameters.blue_max * (1.0 - display_parameters.background);
	  if (display_parameters.lighting)
	  {
	    v[0] = r;
	    v[1] = g;
	    v[2] = b;
	    v[3] = 1.0;
	    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, v);
	  }
	  else
	    glColor3f(r, g, b);
	  glPushMatrix();
	  glTranslatef(x, y, z);
	  if (display_parameters.solid_flag)
	    render_solid_box(display_parameters.boxsize, display_parameters.boxsize, display_parameters.boxsize);
	  else
	    render_wireframe_box(display_parameters.boxsize, display_parameters.boxsize, display_parameters.boxsize);
	  glPopMatrix();
	}
      }
    }
  }
  if (display_parameters.show_extrusion)
  {
    glTranslatef(0.0, -0.5, 0.0);
    for (i = cone.extrusion_planes; i > 0; i--)
    {
      glTranslatef(0.0, -display_parameters.extrusion_boxheight, 0.0);
      y = (cone.extrusion_index + i) % cone.extrusion_planes;
      for (x = xmin; x < xmax; x++)
      {
	for (z = zmin; z < zmax; z++)
	{
	  if (cone.extrusion[x][y][z][0] != -1.0)
	  {
	    r = display_parameters.background;
	    g = display_parameters.background;
	    b = display_parameters.background;
	    if ((display_parameters.red_factor != NO_FACTOR) && (display_parameters.red_max > 0.0))
	      r += cone.extrusion[x][y][z][display_parameters.red_factor] / display_parameters.red_max * (1.0 - display_parameters.background);
	    if ((display_parameters.green_factor != NO_FACTOR) && (display_parameters.green_max > 0.0))
	      g += cone.extrusion[x][y][z][display_parameters.green_factor] / display_parameters.green_max * (1.0 - display_parameters.background);
	    if ((display_parameters.blue_factor != NO_FACTOR) && (display_parameters.blue_max > 0.0))
	      b += cone.extrusion[x][y][z][display_parameters.blue_factor] / display_parameters.blue_max * (1.0 - display_parameters.background);
	    if (display_parameters.lighting)
	    {
	      v[0] = r;
	      v[1] = g;
	      v[2] = b;
	      v[3] = 1.0;
	      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, v);
	    }
	    else
	      glColor3f(r, g, b);
	    glPushMatrix();
	    glTranslatef(x, 0.0, z);
	    if (display_parameters.extrusion_box_flag)
	    {
	      if (display_parameters.solid_flag)
		render_solid_box(display_parameters.boxsize, display_parameters.boxsize * display_parameters.extrusion_boxheight, display_parameters.boxsize);
	      else
		render_wireframe_box(display_parameters.boxsize, display_parameters.boxsize * display_parameters.extrusion_boxheight, display_parameters.boxsize);
	    }
	    else
	      xz_boxf(display_parameters.boxsize, display_parameters.boxsize, display_parameters.solid_flag);
	    glPopMatrix();
	  }
	}
      }
    }
  }
}


void display_init(void)
{
  double maxwidth;
  GLfloat v[4];

  maxwidth = cone.xwidth;
  if (display_parameters.show_extrusion)
  {
    if (maxwidth < cone.ywidth + display_parameters.extrusion_boxheight * cone.extrusion_length)
      maxwidth = cone.ywidth + display_parameters.extrusion_boxheight * cone.extrusion_length;
  }
  else
  {
    if (maxwidth < cone.ywidth)
      maxwidth = cone.ywidth;
  }
  if (maxwidth < cone.zwidth)
    maxwidth = cone.zwidth;
  if (display_parameters.solid_flag)
    glEnable(GL_DEPTH_TEST);
  else
    glDisable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, (GLfloat) window_width / (GLfloat) window_height, 1.0, 5000.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (display_parameters.lighting)
  {
    glEnable(GL_LIGHTING);
    if (display_parameters.extrusion_box_flag)
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    else
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    v[0] = -0.5;
    v[1] = 1.0;
    v[2] = 1.0;
    v[3] = 0.0;
    glLightfv(GL_LIGHT0, GL_POSITION, v);
    glEnable(GL_LIGHT0);
  }
  else
    glDisable(GL_LIGHTING);
  if (display_parameters.show_extrusion)
    gluLookAt(0.0, -0.5 * display_parameters.extrusion_boxheight * cone.extrusion_length, -maxwidth * 1.8, 0.0, -0.5 * display_parameters.extrusion_boxheight * cone.extrusion_length, 1.0, 0.0, 1.0, 0.0);
  else
    gluLookAt(0.0, 0.0, -maxwidth * 1.8, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0);
}


void display_eps(void)
{
  display_init();
  display_model();
}


void display(void)
{
  double t;
  time_t starttime;

  starttime = time(NULL);
  glViewport(0, 0, window_width, window_height);
  display_init();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  display_model();
  glFlush();
  glutSwapBuffers();
  t = difftime(time(NULL), starttime);
  rendering_time += t;
  num_renderings++;
}


void myinit (void) 
{
  glShadeModel(GL_FLAT);
}


void myReshape(int w, int h)
{
  window_width = w;
  window_height = h;
}


void tick(int t)
{
  int i;

  if (time_running)
  {
    for (i = 0; i < cone.num_cells; i++)
      process_expression(&(cone.cell[0][0][i].transsys_instance));
    diffuse(cone.num_cells, cone.cell[0][0]);
    time_step++;
    extend_extrusion(&cone);
    if (trace_flag)
    {
      trace_flag = 0;
      time_running = 0;
    }
  }
  if (time_running)
    glutTimerFunc(10, tick, t);
  display_parameters.angles.x += display_parameters.spin_angles.x;
  display_parameters.angles.y += display_parameters.spin_angles.y;
  display_parameters.angles.z += display_parameters.spin_angles.z;
  glutPostRedisplay();
}


void toggle_running(void)
{
  time_running = !time_running;
  if (time_running)
  {
    fprintf(stderr, "continuing at time step %lu\n", time_step);
    glutTimerFunc(10, tick, 1);
  }
  else
    fprintf(stderr, "pausing at time step %lu\n", time_step);
}


void toggle_solid(void)
{
  display_parameters.solid_flag = !display_parameters.solid_flag;
}


void run_one_step(void)
{
  time_running = 1;
  trace_flag = 1;
  fprintf(stderr, "running time step %lu\n", time_step);
}


int zero_display_rgb(DISPLAY_PARAMETERS *d)
{
  d->red_max = 0.0;
  d->green_max = 0.0;
  d->blue_max = 0.0;
  return (0);
}


int init_display_parameters(DISPLAY_PARAMETERS *d, MERISTEM_CONE *cone)
{
  d->show_extrusion = 1;
  d->solid_flag = 1;
  d->lighting = 1;
  d->extrusion_box_flag = 1;
  d->xmin = 0;
  d->xmax = cone->xwidth;
  d->ymin = 0;
  d->ymax = cone->ywidth;
  d->zmin = 0;
  d->zmax = cone->zwidth;
  if (cone->transsys->num_factors > 0)
    d->red_factor = 0;
  else
    d->red_factor = NO_FACTOR;
  if (cone->transsys->num_factors > 1)
    d->green_factor = 1;
  else
    d->green_factor = NO_FACTOR;
  if (cone->transsys->num_factors > 2)
    d->blue_factor = 2;
  else
    d->blue_factor = NO_FACTOR;
  d->background = 0.05;
  d->boxsize = boxsize;
  d->extrusion_boxheight = 0.1;
  d->current_axis = 0;
  d->current_factor = 0;
  d->current_minmax = 0;
  d->angles.x = 0.0;
  d->angles.y = 0.0;
  d->angles.z = 0.0;
  d->spin_angles.x = 0.0;
  d->spin_angles.y = 0.0;
  d->spin_angles.z = 0.0;
  zero_display_rgb(d);
  return (0);
}


void reset_cone(void)
{
  int x, y, z, i;

  for (x = 0; x < cone.xwidth; x++)
  {
    for (y = 0; y < cone.ywidth; y++)
    {
      for (z = 0; z < cone.zwidth; z++)
      {
	for (i = 0; i < cone.transsys->num_factors; i++)
	  cone.cell[x][y][z].transsys_instance.factor_concentration[i] = 0.0;
      }
    }
  }
  for (x = 0; x < cone.xwidth; x++)
  {
    for (y = 0; y < cone.extrusion_length; y++)
    {
      for (z = 0; z < cone.zwidth; z++)
      {
	for (i = 0; i < cone.transsys->num_factors; i++)
	  cone.extrusion[x][y][z][i] = -1.0;
      }
    }
  }
  cone.extrusion_index = 0;
  cone.extrusion_planes = 0;
  time_step = 0;
}


void reset_display(void)
{
  init_display_parameters(&display_parameters, &cone);
}


void reset_rgb(void)
{
  zero_display_rgb(&display_parameters);
}


void cutplane_up(void)
{
  switch (display_parameters.current_axis)
  {
  case 0:
    if (display_parameters.current_minmax)
      display_parameters.xmax++;
    else
    {
      if (display_parameters.xmin < (display_parameters.xmax - 1))
	display_parameters.xmin++;
    }
    break;
  case 1:
    if (display_parameters.current_minmax)
      display_parameters.ymax++;
    else
    {
      if (display_parameters.ymin < (display_parameters.ymax - 1))
	display_parameters.ymin++;
    }
    break;
  case 2:
    if (display_parameters.current_minmax)
      display_parameters.zmax++;
    else
    {
      if (display_parameters.zmin < (display_parameters.zmax - 1))
	display_parameters.zmin++;
    }
    break;
  default:
    fprintf(stderr, "bad current axis %d\n", display_parameters.current_axis);
    break;
  }
}


void cutplane_down(void)
{
  switch (display_parameters.current_axis)
  {
  case 0:
    if (display_parameters.current_minmax)
    {
      if (display_parameters.xmin < (display_parameters.xmax - 1))
	display_parameters.xmax--;
    }
    else
      display_parameters.xmin--;
    break;
  case 1:
    if (display_parameters.current_minmax)
    {
	if (display_parameters.zmin < (display_parameters.zmax - 1))
	  display_parameters.ymax--;
    }
    else
      display_parameters.ymin--;
    break;
  case 2:
    if (display_parameters.current_minmax)
    {
      if (display_parameters.zmin < (display_parameters.zmax - 1))
	display_parameters.zmax--;
    }
    else
      display_parameters.zmin--;
    break;
  default:
    fprintf(stderr, "bad current axis %d\n", display_parameters.current_axis);
    break;
  }
}


void toggle_extrusion(void)
{
  display_parameters.show_extrusion = !display_parameters.show_extrusion;
}


void toggle_extrusion_box_flag(void)
{
  display_parameters.extrusion_box_flag = !display_parameters.extrusion_box_flag;
}


void toggle_lighting(void)
{
  display_parameters.lighting = !display_parameters.lighting;
}

void toggle_minmax(void)
{
  display_parameters.current_minmax = !display_parameters.current_minmax;
}


void set_red(void)
{
  display_parameters.red_factor = display_parameters.current_factor;
  fprintf(stderr, "red is factor \"%s\"\n", cone.transsys->factor_list[display_parameters.current_factor].name);
}


void set_green(void)
{
  display_parameters.green_factor = display_parameters.current_factor;
  fprintf(stderr, "green is factor \"%s\"\n", cone.transsys->factor_list[display_parameters.current_factor].name);
}


void set_blue(void)
{
  display_parameters.blue_factor = display_parameters.current_factor;
  fprintf(stderr, "blue is factor \"%s\"\n", cone.transsys->factor_list[display_parameters.current_factor].name);
}


void red_off(void)
{
  display_parameters.red_factor = NO_FACTOR;
  fprintf(stderr, "red off\n");
}


void green_off(void)
{
  display_parameters.green_factor = NO_FACTOR;
  fprintf(stderr, "green off\n");
}


void blue_off(void)
{
  display_parameters.blue_factor = NO_FACTOR;
  fprintf(stderr, "blue off\n");
}


void spin_backward(void)
{
  /* fprintf(stderr, "spinning backward\n"); */
  switch (display_parameters.current_axis)
  {
  case 0:
    display_parameters.angles.x += 355.0;
    if (display_parameters.angles.x >= 360.0)
      display_parameters.angles.x -= 360.0;
    break;
  case 1:
    display_parameters.angles.y += 355.0;
    if (display_parameters.angles.y >= 360.0)
      display_parameters.angles.y -= 360.0;
    break;
  case 2:
    display_parameters.angles.z += 355.0;
    if (display_parameters.angles.z >= 360.0)
      display_parameters.angles.z -= 360.0;
    break;
  default:
    fprintf(stderr, "bad current axis %d\n", display_parameters.current_axis);
    break;
  }
}


void spin_forward(void)
{
  /* fprintf(stderr, "spinning forward\n"); */
  switch (display_parameters.current_axis)
  {
  case 0:
    display_parameters.angles.x += 5.0;
    if (display_parameters.angles.x >= 360.0)
      display_parameters.angles.x -= 360.0;
    break;
  case 1:
    display_parameters.angles.y += 5.0;
    if (display_parameters.angles.y >= 360.0)
      display_parameters.angles.y -= 360.0;
    break;
  case 2:
    display_parameters.angles.z += 5.0;
    if (display_parameters.angles.z >= 360.0)
      display_parameters.angles.z -= 360.0;
    break;
  default:
    fprintf(stderr, "bad current axis %d\n", display_parameters.current_axis);
    break;
  }
}


void advance_factor(void)
{
  display_parameters.current_factor++;
  display_parameters.current_factor %= cone.transsys->num_factors;
  fprintf(stderr, "current factor is \"%s\"\n", cone.transsys->factor_list[display_parameters.current_factor].name);
}


void advance_axis(void)
{
  display_parameters.current_axis++;
  display_parameters.current_axis%= 3;
  fprintf(stderr, "current axis is: %c\n", 'x' + display_parameters.current_axis);
}


void toggle_spin(void)
{
  switch (display_parameters.current_axis)
  {
  case 0:
    if (display_parameters.spin_angles.x)
      display_parameters.spin_angles.x = 0.0;
    else
      display_parameters.spin_angles.x = 1.0;
    break;
  case 1:
    if (display_parameters.spin_angles.y)
      display_parameters.spin_angles.y = 0.0;
    else
      display_parameters.spin_angles.y = 1.0;
    break;
  case 2:
    if (display_parameters.spin_angles.z)
      display_parameters.spin_angles.z = 0.0;
    else
      display_parameters.spin_angles.z = 1.0;
    break;
  default:
    fprintf(stderr, "bad current axis %d\n", display_parameters.current_axis);
    break;
  }
}


void write_postscript(void)
{
  static int postscript_fileno = 0;
  char fname[256];

  sprintf(fname, "cellcone%d.eps", postscript_fileno++);
  outputEPS(cone_polygon_complexity(&cone, &display_parameters) + 4096, 1, fname);
}


void print_info(void)
{
  int i;

  printf("time step: %lu\n", time_step);
  printf("%f seconds for %ld renderings (%f renderings / second)\n", rendering_time, num_renderings, num_renderings / rendering_time);
  for (i = 0; i < cone.transsys->num_factors; i++)
  {
    printf("total concentration of \"%s\": %f\n",
            cone.transsys->factor_list[i].name,
	    total_concentration(cone.num_cells, cone.cell[0][0], i));
  }
  printf("current axis is: %c\n", 'x' + display_parameters.current_axis);
  printf("current factor is \"%s\"\n", cone.transsys->factor_list[display_parameters.current_factor].name);
  if (display_parameters.red_factor != NO_FACTOR)
    printf("colors: red is factor \"%s\"\n", cone.transsys->factor_list[display_parameters.red_factor].name);
  if (display_parameters.green_factor != NO_FACTOR)
    printf("colors: green is factor \"%s\"\n", cone.transsys->factor_list[display_parameters.green_factor].name);
  if (display_parameters.blue_factor != NO_FACTOR)
    printf("colors: blue is factor \"%s\"\n", cone.transsys->factor_list[display_parameters.blue_factor].name);
  printf("angles: x: %f, y: %f, z: %f\n", display_parameters.angles.x,
	  display_parameters.angles.y, display_parameters.angles.z);
  printf("visibility: x: %d ... %d, y: %d ... %d, z: %d ... %d\n",
	  display_parameters.xmin, display_parameters.xmax,
	  display_parameters.ymin, display_parameters.ymax,
	  display_parameters.zmin, display_parameters.zmax);
}


int init_factor_concentrations(const MERISTEM_CONE *cone, double max_concentration)
{
  int x, y, z, i;

  for (x = 0; x < cone->xwidth; x++)
  {
    for (y = 0; y < cone->ywidth; y++)
    {
      for (z = 0; z < cone->zwidth; z++)
      {
	for (i = 0; i < cone->transsys->num_factors; i++)
	  cone->cell[x][y][z].transsys_instance.factor_concentration[i] = urandom_double() * max_concentration;
      }
    }
  }
  return (0);
}


static void process_keystroke(unsigned char key, int x, int y)
{
  if (key == 27)
    exit(EXIT_SUCCESS);
  else if (key == 'r')
    set_red();
  else if (key == 'g')
    set_green();
  else if (key == 'b')
    set_blue();
  else if (key == 'R')
    red_off();
  else if (key == 'G')
    green_off();
  else if (key == 'B')
    green_off();
  else if (key == 'm')
    toggle_minmax();
  else if (key == 'a')
    advance_axis();
  else if (key == 'f')
    advance_factor();
  else if (key == 's')
    toggle_spin();
  else if (key == 't')
    run_one_step();
  else if (key == 'p')
    toggle_running();
  else if (key == 'w')
    toggle_solid();
  else if (key == 'l')
    toggle_lighting();
  else if (key == 'z')
    reset_rgb();
  else if (key == '0')
    reset_cone();
  else if (key == 'e')
    toggle_extrusion();
  else if (key == 'E')
    toggle_extrusion_box_flag();
  else if (key == '\r')
    reset_display();
  else if (key == 'P')
    write_postscript();
  else if (key == ' ')
    print_info();
  else
    fprintf(stderr, "process_keystroke: unknown key \'%c\' (%d)\n", key, key);
}


static void process_special_keystroke(int key, int x, int y)
{
  if (key == GLUT_KEY_LEFT)
    spin_backward();
  else if (key == GLUT_KEY_RIGHT)
    spin_forward();
  else if (key == GLUT_KEY_UP)
    cutplane_up();
  else if (key == GLUT_KEY_DOWN)
    cutplane_down();
}


int main(int argc, char** argv)
{
  int xwidth = 10, ywidth = 6, zwidth = 10, extrusion_length = 30;
  double max_concentration = 0.0;
  int yyreturn;
  char *outfile_name = NULL;
  FILE *outfile;
  int oc;
  extern int optind;
  extern char *optarg;

  while ((oc = getopt(argc, argv, "m:l:s:x:y:z:h")) != -1)
  {
    switch(oc)
    {
    case 'l':
      extrusion_length = strtol(optarg, NULL, 10);
      break;
    case 's':
      boxsize = strtod(optarg, NULL);
      break;
    case 'x':
      xwidth = strtol(optarg, NULL, 10);
      break;
    case 'y':
      ywidth = strtol(optarg, NULL, 10);
      break;
    case 'z':
      zwidth = strtol(optarg, NULL, 10);
      break;
    case 'm':
      max_concentration = strtod(optarg, NULL);
      break;
    case 'h':
      printf("-l <num>: specify extrusion length\n");
      printf("-m <num>: specify maximum of randomly generated initial concenctration\n");
      printf("-s <num>: specify relative box size\n");
      printf("-x <num>: specify x width\n");
      printf("-y <num>: specify y width\n");
      printf("-z <num>: specify z width\n");
      printf("-h: print this help and exit\n");
      exit(EXIT_SUCCESS);
    default:
      fprintf(stderr, "unknown option -%c\n", oc);
      exit(EXIT_FAILURE);
    }
  }
  if (optind < argc)
    yyin_name = argv[optind++];
  if (optind < argc)
    outfile_name = argv[optind++];
  if (yyin_name)
  {
    if ((yyin = fopen(yyin_name, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open \"%s\" for input -- exit\n", yyin_name);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    yyin = stdin;
    yyin_name = "stdin";
  }
  if (outfile_name)
  {
    if ((outfile = fopen(outfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open \"%s\" for output -- exit\n", outfile_name);
      if (yyin != stdin)
        fclose(yyin);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    outfile = stdout;
    outfile_name = "stdout";
  }
  yyreturn = yyparse();
  if (yyreturn)
  {
    fprintf(stderr, "error parsing transsys from \"%s\" -- exit\n", yyin_name);
    exit(EXIT_FAILURE);
  }
  arrange_transsys_arrays(parsed_transsys);
  init_cone(&cone);
  create_cone(xwidth, ywidth, zwidth, parsed_transsys, extrusion_length, &cone);
  init_factor_concentrations(&cone, max_concentration);
  /* printf_cone(&cone); */
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(640, 480);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  if (glutCreateWindow(argv[0]) == GL_FALSE)
    exit(1);
  myinit ();
  reset_display();
  glutKeyboardFunc(process_keystroke);
  glutSpecialFunc(process_special_keystroke);
  glutReshapeFunc(myReshape);
  glutDisplayFunc(display);
  glutMainLoop();
  free_cone(&cone);
  free_transsys_list(parsed_transsys);
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  return (EXIT_SUCCESS);
}

