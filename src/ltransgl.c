/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.7  2005/05/16 12:02:10  jtk
 * in transition from distance matrices to contact graphs
 *
 * Revision 1.6  2005/04/12 17:49:46  jtk
 * minor changes to make string and state reporting more convenient
 *
 * Revision 1.5  2005/04/05 10:12:39  jtk
 * made diffusion consistent (no oscillation due to overshooting), small fixes
 *
 * Revision 1.4  2005/04/04 21:30:54  jtk
 * minor changes / comments
 *
 * Revision 1.3  2005/04/04 09:39:54  jtk
 * added lsys capabilities to transexpr, various small changes
 *
 * Revision 1.2  2005/03/29 17:33:02  jtk
 * introduced arrayed lsys string, with symbol distance matrix.
 *
 * Revision 1.1.1.1  2005/03/08 17:12:02  jtk
 * new cvs after loss at INB
 *
 * Revision 1.2  2003/02/26 17:50:04  kim
 * fixed bug of returning success when parse errors occurred
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>

#include "trconfig.h"
#include "trtypes.h"
#include "transsys.h"

/* FIXME: should grow dynamically */
#define LSTR_BLOCK_STRING_MAX 5000


typedef struct
{
  double x, y, z;
} VECTOR3D;


typedef struct
{
  double w, x, y, z;
} QUATERNION;

typedef struct
{
  const LSYS *lsys;
  int num_strings;
  int string_index;
  LSYS_STRING *lstr[LSTR_BLOCK_STRING_MAX];
} LSTR_BLOCK;


typedef struct tag_turtle_state
{
  struct tag_turtle_state *next;
  VECTOR3D position;
  QUATERNION orientation;
  double red, green, blue;
} TURTLE_STATE;


typedef struct
{
  int polygon_mode;
  int animate_mode;
  int lighting_mode;
  int width, height;
  VECTOR3D view_position;
  QUATERNION view_orientation;
  VECTOR3D model_position;
  QUATERNION model_orientation;
  const char *ppm_filename;
} RENDER_PARAMETERS;


typedef struct
{
  int button, state, modifiers, x, y;
} MOUSE_BLOCK;


char *prgname = "ltransgl";
int verbose = 0;

static MOUSE_BLOCK mouse_block = { 0, GLUT_UP, 0, 0, 0 };

static RENDER_PARAMETERS render_parameters =
{
  GL_FILL,
  0, 1,
  0, 0,
  { 0.0, 0.0, 10.0 },
  { 1.0, 0.0, 0.0, 0.0 },
  { 0.0, 0.0, 0.0 },
  { 1.0, 0.0, 0.0, 0.0 },
  NULL
};
static LSTR_BLOCK lstr_block;
static GLUquadricObj *glu_quadric_object = NULL;


QUATERNION rotation_quaternion(double angle, const VECTOR3D *axis)
{
  QUATERNION q;
  double s = sin(angle * 0.5);

  q.w = cos(angle * 0.5);
  q.x = axis->x * s;
  q.y = axis->y * s;
  q.z = axis->z * s;
  return (q);
}


QUATERNION quaternion_conjugate(QUATERNION q)
{
  QUATERNION qc;

  qc.w = q.w;
  qc.x = -q.x;
  qc.y = -q.y;
  qc.z = -q.z;
  return (qc);
}


QUATERNION quaternion_add(QUATERNION q1, QUATERNION q2)
{
  QUATERNION q;

  q.w = q1.w + q2.w;
  q.x = q1.x + q2.x;
  q.y = q1.y + q2.y;
  q.z = q1.z + q2.z;
  return (q);
}


QUATERNION quaternion_mult(QUATERNION q1, QUATERNION q2)
{
  QUATERNION q;

  q.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
  q.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
  q.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
  q.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
  return (q);
}


QUATERNION quaternion_mult_real(double x, QUATERNION q1)
{
  QUATERNION q;

  q.w = x * q1.w;
  q.x = x * q1.x;
  q.y = x * q1.y;
  q.z = x * q1.z;
  return (q);
}


double quaternion_abs(QUATERNION q)
{
  return (sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z));
}


QUATERNION quaternion_normalized(QUATERNION q)
{
  return (quaternion_mult_real(1.0 / quaternion_abs(q), q));
}


int quaternion_to_glmatrixd(GLdouble *m, const QUATERNION *q)
{
  m[0] = q->w * q->w + q->x * q->x - q->y * q->y - q->z * q->z;
  m[1] = 2.0 * q->x * q->y - 2.0 * q->w * q->z;
  m[2] = 2.0 * q->x * q->z + 2.0 * q->w * q->y;
  m[3] = 0.0;
  m[4] = 2.0 * q->x * q->y + 2.0 * q->w * q->z;
  m[5] = q->w * q->w - q->x * q->x + q->y * q->y - q->z * q->z;
  m[6] = 2.0 * q->y * q->z - 2.0 * q->w * q->x;
  m[7] = 0.0;
  m[8] = 2.0 * q->x * q->z - 2.0 * q->w * q->y;
  m[9] = 2.0 * q->y * q->z + 2.0 * q->w * q->x;
  m[10] = q->w * q->w - q->x * q->x - q->y * q->y + q->z * q->z;
  m[11] = 0.0;
  m[12] = 0.0;
  m[13] = 0.0;
  m[14] = 0.0;
  m[15] = q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z;
  /* fprintf(stderr, "rotation matrix check: %f %f %f\n", m[0] * m[0] + m[1] * m[1] + m[2] * m[2], m[4] * m[4] + m[5] * m[5] + m[6] * m[6], m[8] * m[8] + m[9] * m[9] + m[10] * m[10]); */
  return (0);
}


VECTOR3D vector3d_rotation(VECTOR3D v, const QUATERNION *q)
{
  QUATERNION qv, qr;
  VECTOR3D vr;

  qv.w = 0.0;
  qv.x = v.x;
  qv.y = v.y;
  qv.z = v.z;
  qr = quaternion_mult(quaternion_mult(*q, qv), quaternion_conjugate(*q));
  vr.x = qr.x;
  vr.y = qr.y;
  vr.z = qr.z;
  return (vr);
}


static void init_lstr_block(LSTR_BLOCK *lstr_block, const LSYS *lsys)
{
  int i;

  lstr_block->lsys = lsys;
  lstr_block->num_strings = 0;
  lstr_block->string_index = 0;
  for (i = 0; i < LSTR_BLOCK_STRING_MAX; i++)
    lstr_block->lstr[i] = NULL;
}


static void free_lstr_block(LSTR_BLOCK *lstr_block)
{
  int i;

  for (i = 0; i < lstr_block->num_strings; i++)
    free_lsys_string(lstr_block->lstr[i]);
}


static int compute_derived_strings(LSTR_BLOCK *lstr_block)
{
  if (lstr_block->string_index >= LSTR_BLOCK_STRING_MAX)
  {
    fprintf(stderr, "compute_derived_strings: clipping %d to %d\n", lstr_block->string_index, LSTR_BLOCK_STRING_MAX - 1);
    lstr_block->string_index = LSTR_BLOCK_STRING_MAX - 1;
  }
  if (lstr_block->string_index < 0)
    lstr_block->string_index = 0;
  if (lstr_block->string_index < lstr_block->num_strings)
    return (0);
  if (lstr_block->num_strings == 0)
  {
    lstr_block->lstr[0] = axiom_string(lstr_block->lsys);
    if (lstr_block->lstr[0] == NULL)
    {
      fprintf(stderr, "compute_derived_strings: axiom_string() returned NULL\n");
      return (-1);
    }
    lstr_block->num_strings++;
  }
  for ( ; lstr_block->num_strings <= lstr_block->string_index; lstr_block->num_strings++)
  {
    /* fprintf(stderr, "compute_derived_strings: string #%d\n", lstr_block->num_strings); */
    lsys_string_expression(lstr_block->lstr[lstr_block->num_strings - 1]);
    lsys_string_diffusion(lstr_block->lstr[lstr_block->num_strings - 1]);
    lstr_block->lstr[lstr_block->num_strings] = derived_string(lstr_block->lstr[lstr_block->num_strings - 1]);
    if (lstr_block->lstr[lstr_block->num_strings] == NULL)
    {
      fprintf(stderr, "compute_derived_strings: derived_string() returned NULL\n");
      return (-1);
    }
  }
  return (0);
}


static void print_render_parameters(FILE *f, const RENDER_PARAMETERS *rp)
{
  fprintf(f, "polygon_mode: %d\n", (int) rp->polygon_mode);
  fprintf(f, "view position:   (%f, %f, %f)\n", rp->view_position.x, rp->view_position.y, rp->view_position.z);
  fprintf(f, "model origin:    (%f, %f, %f)\n", rp->model_position.x, rp->model_position.y, rp->model_position.z);
  fprintf(f, "model orientation:  (%f %f, %f, %f)\n", rp->model_orientation.w, rp->model_orientation.x, rp->model_orientation.y, rp->model_orientation.z);
}


static void print_lstr_block(FILE *f, const LSTR_BLOCK *lstr_block)
{
  int i, lastlen_pos = 0;
  size_t len, lastlen = 0;

  fprintf(f, "lsys \"%s\", %d strings derived\n", lstr_block->lsys->name, lstr_block->num_strings);
  for (i = 0; i < lstr_block->num_strings; i++)
  {
    len = lsys_string_length(lstr_block->lstr[i]);
    if (len != lastlen)
    {
      fprintf(f, "#%4d - %4d: %lu symbols\n", lastlen_pos, i, (unsigned long) lsys_string_length(lstr_block->lstr[i]));
      lastlen = len;
      lastlen_pos = i + 1;
    }
  }
  if (lstr_block->string_index < lstr_block->num_strings)
  {
    fprintf(f, "current string is #%d:\n", lstr_block->string_index);
    fprint_lsys_string(f, lstr_block->lstr[lstr_block->string_index], "\n");
    if (verbose)
    {
      fprint_lsys_string_contact_graph(f, lstr_block->lstr[lstr_block->string_index]);
    }
  }
}


static void print_lstrings(FILE *f)
{
  print_lstr_block(f, &lstr_block);
}


static void print_globals(FILE *f)
{
  fprintf(f, "current string is #%d out of %d, length = %lu\n", lstr_block.string_index, lstr_block.num_strings, (unsigned long) lstr_block.lstr[lstr_block.string_index]->num_symbols);
  print_render_parameters(f, &render_parameters);
}


static TURTLE_STATE *pop_turtle_state(TURTLE_STATE *s)
{
  TURTLE_STATE *sn;

  sn = s->next;
  free(s);
  return (sn);
}


static TURTLE_STATE *push_turtle_state(TURTLE_STATE *s)
{
  TURTLE_STATE *sn;

  sn = (TURTLE_STATE *) malloc(sizeof(TURTLE_STATE));
  if (sn == NULL)
    return (NULL);
  if (s)
    *sn = *s;
  else
  {
    sn->position.x = 0.0;
    sn->position.y = 0.0;
    sn->position.z = 0.0;
    sn->orientation.w = 1.0;
    sn->orientation.x = 0.0;
    sn->orientation.y = 0.0;
    sn->orientation.z = 0.0;
    sn->red = 0.7;
    sn->green = 0.7;
    sn->blue = 0.7;
  }
  sn->next = s;
  return (sn);
}


static int turtle_move(TURTLE_STATE *turtle_state, double d)
{
  VECTOR3D v = { 0.0, 0.0, 1.0 };
  QUATERNION q;

  if (turtle_state == NULL)
    return (-1);
  q = quaternion_conjugate(turtle_state->orientation);
  v = vector3d_rotation(v, &q);
  /* fprintf(stderr, "turtle_move: rotated vector is (%f, %f, %f)\n", v.x, v.y, v.z); */
  turtle_state->position.x += d * v.x;
  turtle_state->position.y += d * v.y;
  turtle_state->position.z += d * v.z;
  return (0);
}


static int turtle_turn(TURTLE_STATE *ts, double angle)
{
  VECTOR3D axis = { 0.0, 1.0, 0.0 };

  if (ts == NULL)
    return (-1);
  ts->orientation = quaternion_normalized(quaternion_mult(rotation_quaternion(-angle * M_PI / 180.0, &axis), ts->orientation));
  return (0);
}


static int turtle_roll(TURTLE_STATE *ts, double angle)
{
  VECTOR3D axis = { 0.0, 0.0, 1.0 };

  if (ts == NULL)
    return (-1);
  ts->orientation = quaternion_normalized(quaternion_mult(rotation_quaternion(-angle * M_PI / 180.0, &axis), ts->orientation));
  return (0);
}


static int turtle_bank(TURTLE_STATE *ts, double angle)
{
  VECTOR3D axis = { 1.0, 0.0, 0.0 };

  if (ts == NULL)
    return (-1);
  ts->orientation = quaternion_normalized(quaternion_mult(rotation_quaternion(-angle * M_PI / 180.0, &axis), ts->orientation));
  return (0);
}


static void ltgl_box(double x, double y, double z)
{
  glTranslatef(-x * 0.5, -y * 0.5, -z * 0.5);
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


static int ltgl_setcolor(const TURTLE_STATE *ts, const RENDER_PARAMETERS *r)
{
  GLfloat v[4];

  if (ts == NULL)
    return (-1);
  if (r->lighting_mode)
  {
    v[0] = ts->red;
    v[1] = ts->green;
    v[2] = ts->blue;
    v[3] = 1.0;
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, v);
  }
  else
    glColor3f(ts->red, ts->green, ts->blue);
  return (0);
}


static void ltgl_rotate_q(const QUATERNION *q)
{
  GLdouble m[16];

  quaternion_to_glmatrixd(m, q);
  glMultMatrixd(m);
}


static int ltgl_graphics_primitive(const GRAPHICS_PRIMITIVE *gp, const TRANSSYS_INSTANCE *transsys_instance, TURTLE_STATE **ts, const RENDER_PARAMETERS *rp)
{
  double arg[GRAPHICS_PRIMITIVE_ARGMAX];
  int i;

  for (i = 0; i < gp->num_arguments; i++)
    arg[i] = evaluate_expression(gp->argument[i], &transsys_instance);
  switch (gp->type)
  {
  case GRAPHICS_PUSH:
    *ts = push_turtle_state(*ts);
    break;
  case GRAPHICS_POP:
    *ts = pop_turtle_state(*ts);
    ltgl_setcolor(*ts, rp);
    break;
  case GRAPHICS_MOVE:
    turtle_move(*ts, arg[0]);
    break;
  case GRAPHICS_TURN:
    turtle_turn(*ts, arg[0]);
    break;
  case GRAPHICS_ROLL:
    turtle_roll(*ts, arg[0]);
    break;
  case GRAPHICS_BANK:
    turtle_bank(*ts, arg[0]);
    break;
  case GRAPHICS_SPHERE:
    glPushMatrix();
    glTranslatef((*ts)->position.x, (*ts)->position.y, (*ts)->position.z);
    /* fprintf(stderr, "ltgl_graphics_primitive: sphere at (%f, %f, %f)\n", (*ts)->position.x, (*ts)->position.y, (*ts)->position.z); */
    ltgl_rotate_q(&((*ts)->orientation));
    gluSphere(glu_quadric_object, arg[0], 10, 7);
    glPopMatrix();
    break;
  case GRAPHICS_CYLINDER:
    glPushMatrix();
    glTranslatef((*ts)->position.x, (*ts)->position.y, (*ts)->position.z);
    ltgl_rotate_q(&((*ts)->orientation));
    glTranslatef(0.0, 0.0, -0.5 * arg[1]);
    gluQuadricOrientation(glu_quadric_object, GLU_INSIDE);
    gluDisk(glu_quadric_object, 0.0, arg[0], 10, 1);
    gluQuadricOrientation(glu_quadric_object, GLU_OUTSIDE);
    if (arg[1] > 0.0)
    {
      gluCylinder(glu_quadric_object, arg[0], arg[0], arg[1], 10, 1);
    }
    glTranslatef(0.0, 0.0, arg[1]);
    gluDisk(glu_quadric_object, 0.0, arg[0], 10, 1);
    glPopMatrix();
    break;
  case GRAPHICS_BOX:
    glPushMatrix();
    glTranslatef((*ts)->position.x, (*ts)->position.y, (*ts)->position.z);
    ltgl_rotate_q(&((*ts)->orientation));
    ltgl_box(arg[0], arg[1], arg[2]);
    glPopMatrix();
    break;
  case GRAPHICS_COLOR:
    (*ts)->red = arg[0];
    (*ts)->green = arg[1];
    (*ts)->blue = arg[2];
    ltgl_setcolor(*ts, rp);
    break;
  default:
    fprintf(stderr, "postscript_graphics_primitive: primitive type %d not implemented\n", (int) gp->type);
    return (-1);
  }
  return (0);
}


static void draw_model(LSTR_BLOCK *lstr_block, const RENDER_PARAMETERS *rp)
{
  const SYMBOL_INSTANCE *symbol;
  const GRAPHICS_PRIMITIVE *gp;
  TURTLE_STATE *turtle_state_list;

  compute_derived_strings(lstr_block);
  turtle_state_list = push_turtle_state(NULL);
  for (symbol = lstr_block->lstr[lstr_block->string_index]->symbol; symbol; symbol = symbol->next) /* migrate to arrayed string */
  {
    for (gp = symbol->lsys_string->lsys->symbol_list[symbol->symbol_index].graphics_primitive_list; gp; gp = gp->next) /* migrate */
    {
      ltgl_graphics_primitive(gp, &(symbol->transsys_instance), &turtle_state_list, rp);
      if (turtle_state_list == NULL)
	break;
    }
    if (turtle_state_list == NULL)
      break;
  }
  while (turtle_state_list)
    turtle_state_list = pop_turtle_state(turtle_state_list);
}


static void ltgl_reshape(int width, int height)
{
  /* fprintf(stderr, "ltgl_reshape: start\n"); */
  /* Darstellung auf gesamten Clientbereich des Fensters zulassen */
  render_parameters.width = width;
  render_parameters.height = height;
  glViewport(0, 0, (GLint)width, (GLint)height);

  /* Projektionsmatix initialisieren auf 60 Grad horizontales */
  /* Sichtfeld, Verhaeltnis Breite:Hoehe = 1:1, Clipping fuer z<1 */
  /* und z>200 */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  /* angle, aspect, near Clip, far Clip */
  gluPerspective(60.0, ((double) width) / height, 1.0, 200.0);

  /* Modelview Matrix wieder zur aktuellen Matrix machen */
  glMatrixMode(GL_MODELVIEW);
}


static void write_ppm()
{
  FILE *f;
  GLubyte *pixels;
  size_t arrsize, x, y, base;

  arrsize = render_parameters.width * render_parameters.height * 3;
  pixels = (GLubyte *) malloc(arrsize * sizeof(GLubyte));
  if (pixels == NULL)
  {
    fprintf(stderr, "write_ppm: malloc failed\n");
    return;
  }
  glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, render_parameters.width, render_parameters.height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  if (render_parameters.ppm_filename)
    f = fopen(render_parameters.ppm_filename, "w");
  else
    f = fopen("ltransgl.ppm", "w");
  if (f == NULL)
  {
    fprintf(stderr, "write_ppm: fopen failed\n");
    free(pixels);
    return;
  }
  fprintf(f, "P3\n");
  fprintf(f, "%d %d 255\n", render_parameters.width, render_parameters.height);
  for (y = 0; y < render_parameters.height; y++)
  {
    for (x = 0; x < render_parameters.width; x++)
    {
      base = ((render_parameters.height - y - 1) * render_parameters.width + x) * 3,
      fprintf(f, "%d %d %d\n", pixels[base], pixels[base + 1], pixels[base + 2]);
    }
  }
  fclose(f);
  free(pixels);
}


static void ltgl_animate_func(int dt)
{
  /* fprintf(stderr, "ltgl_animate_func(%d)\n", dt); */
  if (!render_parameters.animate_mode)
    return;
  if (((dt > 0) && (lstr_block.string_index == lstr_block.num_strings - 1)) || ((dt < 0) && (lstr_block.string_index == 0)))
    return;
  lstr_block.string_index += dt;
  glutTimerFunc(40, ltgl_animate_func, dt);
  glutPostRedisplay();
}


static void ltgl_key(unsigned char key, int x, int y)
{
  if (key==27)
  {
    free_lstr_block(&lstr_block);
    free_lsys_list(parsed_lsys);
    free_transsys_list(parsed_transsys);
    exit(0);
  }
  if (key==' ')
  {
    if (render_parameters.polygon_mode == GL_FILL)
      render_parameters.polygon_mode = GL_LINE;
    else
      render_parameters.polygon_mode = GL_FILL;
    glPolygonMode(GL_FRONT_AND_BACK, render_parameters.polygon_mode);
  }
  else if (key == 'n')
    lstr_block.string_index++;
  else if (key == 'N')
    lstr_block.string_index += 50;
  else if (key == '<')
  {
    render_parameters.animate_mode = 1;
    ltgl_animate_func(-1);
  }
  else if (key == '>')
  {
    render_parameters.animate_mode = 1;
    ltgl_animate_func(1);
  }
  else if (key == '.')
    render_parameters.animate_mode = 0;
  else if (key == 'l')
    render_parameters.lighting_mode = !render_parameters.lighting_mode;
  else if (key == 'p')
    lstr_block.string_index--;
  else if (key == 'X')
    render_parameters.view_position.x += 0.2;
  else if (key == 'x')
    render_parameters.view_position.x -= 0.2;
  else if (key == 'Y')
    render_parameters.view_position.y += 0.2;
  else if (key == 'y')
    render_parameters.view_position.y -= 0.2;
  else if (key == 'Z')
    render_parameters.view_position.z += 0.2;
  else if (key == 'z')
    render_parameters.view_position.z -= 0.2;
  else if (key == '|')
  {
    render_parameters.view_orientation.w = 1.0;
    render_parameters.view_orientation.x = 0.0;
    render_parameters.view_orientation.y = 0.0;
    render_parameters.view_orientation.z = 0.0;
    render_parameters.model_orientation.w = 1.0;
    render_parameters.model_orientation.x = 0.0;
    render_parameters.model_orientation.y = 0.0;
    render_parameters.model_orientation.z = 0.0;
  }
  else if (key == 'i')
    write_ppm();
  else if (key == '?')
    print_globals(stdout);
  else if (key == 's')
    print_lstrings(stdout);
  else if (key == 'v')
    verbose = !verbose;
  else
  {
    if (isprint(key))
      fprintf(stderr, "ltgl_key: unrecognized key '%c'\n", key);
    else if (iscntrl(key))
      fprintf(stderr, "ltgl_key: unrecognized key '^%c'\n", key + 64);
    else
      fprintf(stderr, "ltgl_key: unrecognized key #%d\n", key);
  }
  glutPostRedisplay ();
}


static int rotate_model(RENDER_PARAMETERS *rp, int axis_no, double angle)
{
  VECTOR3D axis = { 0.0, 0.0, 0.0 };

  switch(axis_no)
  {
  case 1:
    axis.x = 1.0;
    break;
  case 2:
    axis.y = 1.0;
    break;
  case 3:
    axis.z = 1.0;
    break;
  }
  rp->model_orientation = quaternion_normalized(quaternion_mult(rotation_quaternion(angle / 180.0 * M_PI, &axis), render_parameters.model_orientation));
  return (0);
}


static void ltgl_special_key(int key, int x, int y)
{
  switch (key)
  {
  case GLUT_KEY_LEFT:
    rotate_model(&render_parameters, 2, -5.0);
    break;
  case GLUT_KEY_RIGHT:
    rotate_model(&render_parameters, 2, 5.0);
    break;
  case GLUT_KEY_UP:
    rotate_model(&render_parameters, 1, -5.0);
    break;
  case GLUT_KEY_DOWN:
    rotate_model(&render_parameters, 1, 5.0);
    break;
  case GLUT_KEY_PAGE_UP:
    rotate_model(&render_parameters, 3, -5.0);
    break;
  case GLUT_KEY_PAGE_DOWN:
    rotate_model(&render_parameters, 3, 5.0);
    break;
  case GLUT_KEY_HOME:
    lstr_block.string_index = 0;
    break;
  case GLUT_KEY_END:
    lstr_block.string_index = lstr_block.num_strings - 1;
    break;
  default:
    fprintf(stderr, "ltgl_special_key: unknown special key %d\n", key);
    return;
  }
  glutPostRedisplay();
}


static void ltgl_mouse(int button, int state, int x, int y)
{
  mouse_block.button = button;
  mouse_block.state = state;
  mouse_block.modifiers = glutGetModifiers();
  mouse_block.x = x;
  mouse_block.y = y;
  /* fprintf(stderr, "mouse: button = %d, state = %d, modifiers = %d, x = %d, y = %d\n", mouse_block.button, mouse_block.state, mouse_block.modifiers, mouse_block.x, mouse_block.y); */
}


static void ltgl_motion(int x, int y)
{
  int dx, dy;
  VECTOR3D axis = { 0.0, 0.0, 0.0 };

  if (mouse_block.state == GLUT_DOWN)
  {
    dx = x - mouse_block.x;
    dy = y - mouse_block.y;
    mouse_block.x = x;
    mouse_block.y = y;
    if (mouse_block.button == GLUT_LEFT_BUTTON)
    {
      if ((mouse_block.modifiers & GLUT_ACTIVE_CTRL))
	render_parameters.view_position.z -= dy * 0.1;
      else if ((mouse_block.modifiers & GLUT_ACTIVE_SHIFT))
	render_parameters.view_position.y += dy * 0.1;
      else
	render_parameters.view_position.x -= dx * 0.1;
    }
    if (mouse_block.button == GLUT_RIGHT_BUTTON)
    {
      if ((mouse_block.modifiers & GLUT_ACTIVE_CTRL))
      {
	axis.z = 1.0;
	render_parameters.model_orientation = quaternion_normalized(quaternion_mult(rotation_quaternion(-dx / 180.0 * M_PI, &axis), render_parameters.model_orientation));
      }
      else if ((mouse_block.modifiers & GLUT_ACTIVE_SHIFT))
      {
	axis.y = 1.0;
	render_parameters.model_orientation = quaternion_normalized(quaternion_mult(rotation_quaternion(dx / 180.0 * M_PI, &axis), render_parameters.model_orientation));
      }
      else
      {
	axis.x = 1.0;
	render_parameters.model_orientation = quaternion_normalized(quaternion_mult(rotation_quaternion(-dy / 180.0 * M_PI, &axis), render_parameters.model_orientation));
      }
    }
    glutPostRedisplay();
  }
  /* fprintf(stderr, "motion: %d %d\n", x, y); */
}


static void ltgl_display(void)
{
  GLfloat v[4];

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPolygonMode(GL_FRONT_AND_BACK, render_parameters.polygon_mode);
  /* initialize modelview matrix */
  glLoadIdentity ();
  if (render_parameters.lighting_mode)
  {
    glEnable(GL_LIGHTING);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    v[0] = 0.2;
    v[1] = 0.2;
    v[2] = 0.2;
    v[3] = 1.0;
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, v);
    v[0] = -0.5;
    v[1] = 1.0;
    v[2] = 1.0;
    v[3] = 0.0;
    glLightfv(GL_LIGHT0, GL_POSITION, v);
    glEnable(GL_LIGHT0);
  }
  else
    glDisable(GL_LIGHTING);
  /* View transformation */
  glTranslatef(-render_parameters.view_position.x, -render_parameters.view_position.y, -render_parameters.view_position.z);
  ltgl_rotate_q(&(render_parameters.view_orientation));
  /* Model transformation */
  glTranslatef(render_parameters.model_position.x, render_parameters.model_position.y, render_parameters.model_position.z);
  ltgl_rotate_q(&(render_parameters.model_orientation));
  /* draw the plant */
  draw_model(&lstr_block, &render_parameters);
  glutSwapBuffers();
}


static void glu_quadric_error(GLenum horror)
{
  fprintf(stderr, "GLU quadric error #%d: %s\n", (int) horror, gluErrorString(horror));
}


static void ltgl_init(void)
{
  glu_quadric_object = gluNewQuadric();
  if (glu_quadric_object == NULL)
  {
    fprintf(stderr, "ltgl_init: gluNewQuadric() failed\n");
    exit(EXIT_FAILURE);
  }
  gluQuadricCallback(glu_quadric_object, GLU_ERROR, glu_quadric_error);
  gluQuadricNormals(glu_quadric_object, GLU_FLAT);
  /* Ausgabefenster definieren */
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(640, 480);
  /* Renderkontext mit Z-Buffer, Doublebuffer fuer RGB-Modus anfordern. */
  glutInitDisplayMode(GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE);

  if (glutCreateWindow(prgname) == GL_FALSE)
    exit(1);
  
  /* Callback Funktionen vereinbaren */
  glutReshapeFunc(ltgl_reshape);
  glutKeyboardFunc(ltgl_key);
  glutSpecialFunc(ltgl_special_key);
  glutDisplayFunc(ltgl_display);
  glutMouseFunc(ltgl_mouse);
  glutMotionFunc(ltgl_motion);
  /* enable z-buffer */
  glEnable(GL_DEPTH_TEST);
  /* Vorder- u. Rueckseite der Polygone nur als Randlinien darstellen */
  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
}


int get_triple(double triple[3], char *arg)
{
  char *s, *s1;

  s = arg;
  triple[0] = strtod(s, &s1);
  if ((s == s1) || (*s1 != ','))
    return (-1);
  s = s1 + 1;
  triple[1] = strtod(s, &s1);
  if ((s == s1) || (*s1 != ','))
    return (-1);
  s = s1 + 1;
  triple[2] = strtod(s, &s1);
  if (s == s1)
    return (-1);
  return (0);
}


int main(int argc, char **argv)
{
  int oc;
  GLint mdepth;
  extern char *optarg;
  extern int optind;
  char *outfile_name = NULL;
  FILE *outfile;
  int yyreturn;
  double triple[3];
  int num_initial_derivations = 0;
  unsigned int rndseed = 1;

  glutInit(&argc, argv);
  /* printf("%s\n", gluGetString(GLU_VERSION)); */
  while ((oc = getopt(argc, argv, "d:i:p:r:vh")) != -1)
  {
    switch(oc)
    {
    case 'v':
      verbose = 1;
      fprintf(stderr, "verbose mode activated\n");
      break;
    case 'd':
      num_initial_derivations = strtol(optarg, NULL, 10);
      break;
    case 'i':
      render_parameters.ppm_filename = optarg;
      break;
    case 's':
      rndseed = strtoul(optarg, NULL, 10);
      break;
    case 'r':
      if (get_triple(triple, optarg) == 0)
      {
	rotate_model(&render_parameters, 1, triple[0]);
	rotate_model(&render_parameters, 2, triple[1]);
	rotate_model(&render_parameters, 3, triple[2]);
      }
      else
	fprintf(stderr, "bad angle triple \"%s\"\n", optarg);
      break;
    case 'p':
      if (get_triple(triple, optarg) == 0)
      {
	render_parameters.view_position.x = triple[0];
	render_parameters.view_position.y = triple[1];
	render_parameters.view_position.z = triple[2];
      }
      else
	fprintf(stderr, "bad coordinate triple \"%s\"\n", optarg);
      break;
    case 'h':
      printf("ltransps -- make PostScript ltranssys derivation series\n");
      printf("usage: %s [options] <infile> <outfile>\n", argv[0]);
      printf("options:\n");
      printf("-d <num>: specify initial number of derivations\n");
      printf("-i <filename>: specify image file name (PPM format)\n");
      printf("-p <num>,<num>,<num>: specify initial viewer coordinates\n");
      printf("-r <num>,<num>,<num>: specify initial model rotation angles\n");
      printf("-h: print this help and exit\n");
      exit(EXIT_SUCCESS);
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
    fprintf(stderr, "Failed to parse from \"%s\" -- exit\n", yyin_name);
    if (yyin != stdin)
      fclose(yyin);
    if (outfile != stdout)
      fclose(outfile);
    exit(EXIT_FAILURE);
  }
  if (parsed_lsys == NULL)
  {
    fprintf(stderr, "No lsys found in \"%s\" -- exit\n", yyin_name);
    if (yyin != stdin)
      fclose(yyin);
    if (outfile != stdout)
      fclose(outfile);
    free_transsys_list(parsed_transsys);
    exit(EXIT_FAILURE);
  }    
  ulong_srandom(rndseed);
  init_lstr_block(&lstr_block, parsed_lsys);
  lstr_block.string_index = num_initial_derivations;
  ltgl_init();
  glGetIntegerv(GL_MAX_MODELVIEW_STACK_DEPTH, &mdepth);
  printf("max matrix depth is %d\n", (int) mdepth);
  glutMainLoop();
  fprintf(stderr, "glutMainLoop() finished\n");
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  return (EXIT_SUCCESS);
}

