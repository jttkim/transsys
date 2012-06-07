#include <Python.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include <transsys.h>

/*
 * The API version must be changed manually each time the API is
 * changed.
 */
static char clib_api_version[] = "350";

/*
 * Struct to contain pointers to Python classes that directly correspond
 * to C level classes.
 */
typedef struct
{
  int initialised;
  PyObject *ExpressionNode;
  PyObject *ExpressionNodeValue;
  PyObject *ExpressionNodeIdentifier;
  PyObject *ExpressionNodeBinary;
  PyObject *ExpressionNodeMult;
  PyObject *ExpressionNodeDiv;
  PyObject *ExpressionNodeAdd;
  PyObject *ExpressionNodeSubtract;
  PyObject *ExpressionNodeLower;
  PyObject *ExpressionNodeLowerEqual;
  PyObject *ExpressionNodeGreater;
  PyObject *ExpressionNodeGreaterEqual;
  PyObject *ExpressionNodeEqual;
  PyObject *ExpressionNodeUnequal;
  PyObject *ExpressionNodeNot;
  PyObject *ExpressionNodeAnd;
  PyObject *ExpressionNodeOr;
  PyObject *ExpressionNodeFunction;
  PyObject *ExpressionNodeUniformRandom;
  PyObject *ExpressionNodeGaussianRandom;
  PyObject *ExpressionNodePower;
  PyObject *ExpressionNodeLogarithm;
  PyObject *ExpressionNodeAtan;
  PyObject *PromoterElement;
  PyObject *PromoterElementLink;
  PyObject *PromoterElementConstitutive;
  PyObject *PromoterElementActivate;
  PyObject *PromoterElementRepress;
  PyObject *Factor;
  PyObject *Gene;
  PyObject *TranssysProgram;
  PyObject *GraphicsPrimitive;
  PyObject *GraphicsPrimitiveMove;
  PyObject *GraphicsPrimitivePush;
  PyObject *GraphicsPrimitivePop;
  PyObject *GraphicsPrimitiveTurn;
  PyObject *GraphicsPrimitiveRoll;
  PyObject *GraphicsPrimitiveBank;
  PyObject *GraphicsPrimitiveSphere;
  PyObject *GraphicsPrimitiveCylinder;
  PyObject *GraphicsPrimitiveBox;
  PyObject *GraphicsPrimitiveColor;
  PyObject *Symbol;
  PyObject *Assignment;
  PyObject *LhsSymbol;
  PyObject *ProductionElement;
  PyObject *Rule;
  PyObject *LsysProgram;
  PyObject *TranssysInstance;
  PyObject *SymbolInstance;
  PyObject *LsysSymbolString;
} PYTHON_CLASSES;


typedef struct
{
  PyObject *pythonClass;
  EXPR_NODE_TYPE node_type;
} EXPRESSION_NODETYPE_MAPENTRY;


typedef struct
{
  PyObject *pythonClass;
  GRAPHICS_PRIMITIVE_TYPE graphics_primitive_type;
} GRAPHICS_PRIMITIVETYPE_MAPENTRY;


static EXPRESSION_NODETYPE_MAPENTRY expressionNodetypeMap[] = {
  {NULL, NT_NONE}, /* value */
  {NULL, NT_NONE}, /* identifier */
  {NULL, NT_NONE}, /* add */
  {NULL, NT_NONE}, /* subtract */
  {NULL, NT_NONE}, /* mult */
  {NULL, NT_NONE}, /* div */
  {NULL, NT_NONE}, /* lower */
  {NULL, NT_NONE}, /* lowerEqual*/
  {NULL, NT_NONE}, /* greater */
  {NULL, NT_NONE}, /* greaterEqual */
  {NULL, NT_NONE}, /* equal */
  {NULL, NT_NONE}, /* unequal */
  {NULL, NT_NONE}, /* not */
  {NULL, NT_NONE}, /* and */
  {NULL, NT_NONE}, /* or */
  {NULL, NT_NONE}, /* random */
  {NULL, NT_NONE}, /* gauss */
  {NULL, NT_NONE}, /* power */
  {NULL, NT_NONE}, /* log */
  {NULL, NT_NONE}, /* atan */
  {NULL, NT_NONE}  /* end flagged by NULL */
};


static GRAPHICS_PRIMITIVETYPE_MAPENTRY graphicsPrimitivetypeMap[] = {
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE},
  {NULL, GRAPHICS_NONE}
};


static PYTHON_CLASSES pythonClasses = {
  0,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL
};


typedef enum
{
  CLIB_MSG_TRACE,
  CLIB_MSG_DEVMESSAGE,
  CLIB_MSG_WARNING,
  CLIB_MSG_ERROR,
  CLIB_MSG_FATAL
} CLIB_MSG_IMPORTANCE;


static CLIB_MSG_IMPORTANCE message_importance_threshold = CLIB_MSG_WARNING;


/* #define REFCOUNTDEBUG */

/***** reference count monitoring *********************************************/

#ifdef REFCOUNTDEBUG

#define REFCOUNTDEBUG_STRINGLENGTH 200


typedef struct tag_refcountdebug_record
{
  struct tag_refcountdebug_record *next;
  int refcount;
  PyObject *python_object;
  char name[REFCOUNTDEBUG_STRINGLENGTH];
  int line;
} REFCOUNTDEBUG_RECORD;


typedef struct
{
  REFCOUNTDEBUG_RECORD *record_list;
} REFCOUNTDEBUG_STATE;


REFCOUNTDEBUG_STATE _refcountdebug_state = {NULL};


static REFCOUNTDEBUG_RECORD *_refcountDebug_findRecord(PyObject *o)
{
  REFCOUNTDEBUG_RECORD *r;

  for (r = _refcountdebug_state.record_list; r; r = r->next)
  {
    if (r->python_object == o)
    {
      return (r);
    }
  }
  return (NULL);
}


static void _refcountDebug_setString(char *str, const char *source)
{
  strncpy(str, source, REFCOUNTDEBUG_STRINGLENGTH - 1);
  str[REFCOUNTDEBUG_STRINGLENGTH - 1] = '\0';
}


static void _refcountDebug_setInfo(REFCOUNTDEBUG_RECORD *r, const char *name, int line)
{
  _refcountDebug_setString(r->name, name);
  r->line = line;
}


static REFCOUNTDEBUG_RECORD *_refcountDebug_newRecord(PyObject *o)
{
  REFCOUNTDEBUG_RECORD *r;

  r = (REFCOUNTDEBUG_RECORD *) malloc(sizeof(REFCOUNTDEBUG_RECORD));
  if (r == NULL)
  {
    return (NULL);
  }
  r->next = _refcountdebug_state.record_list;
  r->python_object = o;
  r->refcount = 0;
  _refcountdebug_state.record_list = r;
  return (r);
}


static void _refcountDebug_findOrNewAndIncrement(PyObject *o, const char *name, int line)
{
  REFCOUNTDEBUG_RECORD *r;

  if (o == NULL)
  {
    return;
  }
  r = _refcountDebug_findRecord(o);
  if (r == NULL)
  {
    r = _refcountDebug_newRecord(o);
  }
  if (r != NULL)
  {
    if (r->refcount == 0)
    {
      _refcountDebug_setInfo(r, name, line);
    }
    r->refcount++;
    /* fprintf(stderr, "_refcountDebug_findOrNewAndIncrement: record for %p: refcount = %d\n", (void *) o, r->refcount); */
  }
  else
  {
    fprintf(stderr, "_refcountDebug_findOrNewAndIncrement: newRecord failed\n");
  }
}


static void _refcountDebug_Py_INCREF(PyObject *o, const char *name, int line)
{
  Py_XINCREF(o);
  _refcountDebug_findOrNewAndIncrement(o, name, line);
}


static void _refcountDebug_checkAndDecrement(PyObject *o, const char *funcname, const char *name, int line)
{
  REFCOUNTDEBUG_RECORD *r = _refcountDebug_findRecord(o);

  if (r == NULL)
  {
    fprintf(stderr, "refcountDebug: %s on unknown object at %p\n", funcname, (void *) o);
    fprintf(stderr, "  name: \"%s\"\n", name);
    fprintf(stderr, "  line: %d\n", line);
  }
  else
  {
    if (r->refcount > 0)
    {
      r->refcount--;
      /* fprintf(stderr, "_refcountDebug_Py_DECREF: decremented refcount of %p to %d\n", (void *) o, r->refcount); */
    }
    else
    {
      fprintf(stderr, "refcountDebug: Py_DECREF beyond 0\n");
      fprintf(stderr, "  name: \"%s\"\n", name);
      fprintf(stderr, "  line: %d\n", line);
      fprintf(stderr, "  record name: \"%s\"\n", r->name);
      fprintf(stderr, "  record line: %d\n", r->line);
    }
  }
}


static void _refcountDebug_Py_DECREF(PyObject *o, const char *name, int line)
{
  if (o != NULL)
  {
    _refcountDebug_checkAndDecrement(o, "Py_DECREF", name, line);
  }
  Py_XDECREF(o);
}


static PyObject *_refcountDebug_PyObject_GetAttrString(PyObject *o, char *attr, const char *name, int line)
{
  PyObject *python_object = PyObject_GetAttrString(o, attr);
  char aname[REFCOUNTDEBUG_STRINGLENGTH];

  if (strlen(name) + strlen(attr) + 1 < REFCOUNTDEBUG_STRINGLENGTH)
  {
    sprintf(aname, "%s.%s", name, attr);
  }
  else
  {
    sprintf(aname, "attribute");
  }
  _refcountDebug_findOrNewAndIncrement(python_object, aname, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyTuple_New(int size, const char *name, int line)
{
  PyObject *python_object = PyTuple_New(size);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyList_New(int size, const char *name, int line)
{
  PyObject *python_object = PyList_New(size);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyDict_New(const char *name, int line)
{
  PyObject *python_object = PyDict_New();

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyInt_FromLong(long v, const char *name, int line)
{
  PyObject *python_object = PyInt_FromLong(v);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyFloat_FromDouble(double v, const char *name, int line)
{
  PyObject *python_object = PyFloat_FromDouble(v);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyInstance_New(PyObject *pyClass, PyObject *args, PyObject *kw, const char *name, int line)
{
  PyObject *python_object = PyInstance_New(pyClass, args, kw);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyObject_Call(PyObject *callable_object, PyObject *args, PyObject *kw, const char *name, int line)
{
  PyObject *python_object = PyObject_Call(callable_object, args, kw);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static int _refcountDebug_PyTuple_SetItem(PyObject *tuple, int i, PyObject *item, const char *name, int line)
{
  _refcountDebug_checkAndDecrement(item, "PyTuple_SetItem", name, line);
  return (PyTuple_SetItem(tuple, i, item));
}


static int _refcountDebug_PyList_SetItem(PyObject *list, int i, PyObject *item, const char *name, int line)
{
  _refcountDebug_checkAndDecrement(item, "PyList_SetItem", name, line);
  return (PyList_SetItem(list, i, item));
}


static int _refcountDebug_numDisplayableRecords(int showall)
{
  REFCOUNTDEBUG_RECORD *r;
  int n = 0;

  for (r = _refcountdebug_state.record_list; r; r = r->next)
  {
    if ((r->refcount > 0) || showall)
    {
      n++;
    }
  }
  return (n);
}


static void _refcountDebug_report(const char *title, int showall, const char *filename, int line)
{
  REFCOUNTDEBUG_RECORD *r;

  fprintf(stderr, "+----- refcountDebug: %-45s -----+\n", title);
  fprintf(stderr, "| file: %-64s |\n", filename);
  fprintf(stderr, "| line: %-6d                                                           |\n", line);
  fprintf(stderr, "+------------------------------------------------------------------------+\n");
  if (_refcountDebug_numDisplayableRecords(showall) == 0)
  {
  fprintf(stderr, "| no records                                                             |\n");
  }
  else
  {
    fprintf(stderr, "| %-41s  %10s  %6s   %6s |\n", "reference name", "PyObject", "line", "count");
    fprintf(stderr, "+------------------------------------------------------------------------+\n");
    for (r = _refcountdebug_state.record_list; r; r = r->next)
    {
      if ((r->refcount > 0) || showall)
      {
	fprintf(stderr, "| %-41s  %10p  %6d   %6d |\n", r->name, (void *) r->python_object, r->line, r->refcount);
      }
    }
  }
  fprintf(stderr, "+------------------------------------------------------------------------+\n\n");
}


static void refcountDebug_init(void)
{
  REFCOUNTDEBUG_RECORD *r = _refcountdebug_state.record_list, *r1;

  while (r)
  {
    r1 = r;
    r = r->next;
    free(r1);
  }
  _refcountdebug_state.record_list = NULL;
}


/* should also intercept the functions that steal references */

#define PyObject_GetAttrString(o, s) _refcountDebug_PyObject_GetAttrString(o, s, #o, __LINE__)
#define PyTuple_New(size) _refcountDebug_PyTuple_New(size, "new tuple", __LINE__)
#define PyList_New(size) _refcountDebug_PyList_New(size, "new list", __LINE__)
#define PyDict_New() _refcountDebug_PyDict_New("new dict", __LINE__)
#define PyInt_FromLong(v) _refcountDebug_PyInt_FromLong(v, #v, __LINE__)
#define PyFloat_FromDouble(v) _refcountDebug_PyFloat_FromDouble(v, #v, __LINE__)
#define PyInstance_New(pyclass, args, kw) _refcountDebug_PyInstance_New(pyclass, args, kw, #pyclass, __LINE__)
#define PyObject_Call(callable_object, args, kw) _refcountDebug_PyObject_Call(callable_object, args, kw, #callable_object, __LINE__)
#define PyTuple_SetItem(tuple, i, item) _refcountDebug_PyTuple_SetItem(tuple, i, item, #item, __LINE__)
#define PyList_SetItem(list, i, item) _refcountDebug_PyList_SetItem(list, i, item, #item, __LINE__)
#undef Py_INCREF
#undef Py_DECREF
#undef Py_XINCREF
#undef Py_XDECREF
#define Py_INCREF(o) _refcountDebug_Py_INCREF(o, #o, __LINE__)
#define Py_DECREF(o) _refcountDebug_Py_DECREF(o, #o, __LINE__)
#define Py_XINCREF(o) _refcountDebug_Py_INCREF(o, #o, __LINE__)
#define Py_XDECREF(o) _refcountDebug_Py_DECREF(o, #o, __LINE__)

#define refcountDebug_report(title, showall) _refcountDebug_report(title, showall, __FILE__, __LINE__)

#else /* REFCOUNTDEBUG */

#define refcountDebug_report(title, showall)
#define refcountDebug_init()

#endif /* REFCOUNTDEBUG */

/******************************************************************************/



static void clib_message(CLIB_MSG_IMPORTANCE importance, const char *format, ...)
{
  va_list arglist;
  int imp = (int) importance;

  if (imp >= ((int) message_importance_threshold))
  {
    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
  }
}


/*
 * Get a Python class specified by module and by classname.
 * Returns a pointer to the Python class, NULL in case of failure.
 */
static PyObject *getPythonClass(char *modulename, char *classname)
{
  PyObject *sys_modules = PyImport_GetModuleDict();
  PyObject *module = NULL;
  PyObject *pythonClass = NULL;

  if (sys_modules == NULL)
  {
    return (NULL);
  }
  Py_INCREF(sys_modules);
  module = PyDict_GetItemString(sys_modules, modulename);
  Py_DECREF(sys_modules);
  if (module == NULL)
  {
    return (NULL);
  }
  Py_INCREF(module);
  pythonClass = PyObject_GetAttrString(module, classname);
  Py_DECREF(module);
  return (pythonClass);
}


static void freePythonClasses(void)
{
  Py_XDECREF(pythonClasses.ExpressionNode);
  Py_XDECREF(pythonClasses.ExpressionNodeValue);
  Py_XDECREF(pythonClasses.ExpressionNodeIdentifier);
  Py_XDECREF(pythonClasses.ExpressionNodeBinary);
  Py_XDECREF(pythonClasses.ExpressionNodeMult);
  Py_XDECREF(pythonClasses.ExpressionNodeDiv);
  Py_XDECREF(pythonClasses.ExpressionNodeAdd);
  Py_XDECREF(pythonClasses.ExpressionNodeSubtract);
  Py_XDECREF(pythonClasses.ExpressionNodeLower);
  Py_XDECREF(pythonClasses.ExpressionNodeLowerEqual);
  Py_XDECREF(pythonClasses.ExpressionNodeGreater);
  Py_XDECREF(pythonClasses.ExpressionNodeGreaterEqual);
  Py_XDECREF(pythonClasses.ExpressionNodeEqual);
  Py_XDECREF(pythonClasses.ExpressionNodeUnequal);
  Py_XDECREF(pythonClasses.ExpressionNodeNot);
  Py_XDECREF(pythonClasses.ExpressionNodeAnd);
  Py_XDECREF(pythonClasses.ExpressionNodeOr);
  Py_XDECREF(pythonClasses.ExpressionNodeFunction);
  Py_XDECREF(pythonClasses.ExpressionNodeUniformRandom);
  Py_XDECREF(pythonClasses.ExpressionNodeGaussianRandom);
  Py_XDECREF(pythonClasses.ExpressionNodePower);
  Py_XDECREF(pythonClasses.ExpressionNodeLogarithm);
  Py_XDECREF(pythonClasses.ExpressionNodeAtan);
  Py_XDECREF(pythonClasses.PromoterElement);
  Py_XDECREF(pythonClasses.PromoterElementLink);
  Py_XDECREF(pythonClasses.PromoterElementConstitutive);
  Py_XDECREF(pythonClasses.PromoterElementActivate);
  Py_XDECREF(pythonClasses.PromoterElementRepress);
  Py_XDECREF(pythonClasses.Factor);
  Py_XDECREF(pythonClasses.Gene);
  Py_XDECREF(pythonClasses.TranssysProgram);
  Py_XDECREF(pythonClasses.GraphicsPrimitive);
  Py_XDECREF(pythonClasses.GraphicsPrimitiveMove);
  Py_XDECREF(pythonClasses.GraphicsPrimitivePush);
  Py_XDECREF(pythonClasses.GraphicsPrimitivePop);
  Py_XDECREF(pythonClasses.GraphicsPrimitiveTurn);
  Py_XDECREF(pythonClasses.GraphicsPrimitiveRoll);
  Py_XDECREF(pythonClasses.GraphicsPrimitiveBank);
  Py_XDECREF(pythonClasses.GraphicsPrimitiveSphere);
  Py_XDECREF(pythonClasses.GraphicsPrimitiveCylinder);
  Py_XDECREF(pythonClasses.GraphicsPrimitiveBox);
  Py_XDECREF(pythonClasses.GraphicsPrimitiveColor);
  Py_XDECREF(pythonClasses.Symbol);
  Py_XDECREF(pythonClasses.Assignment);
  Py_XDECREF(pythonClasses.LhsSymbol);
  Py_XDECREF(pythonClasses.ProductionElement);
  Py_XDECREF(pythonClasses.Rule);
  Py_XDECREF(pythonClasses.LsysProgram);
  Py_XDECREF(pythonClasses.TranssysInstance);
  Py_XDECREF(pythonClasses.SymbolInstance);
  Py_XDECREF(pythonClasses.LsysSymbolString);
  pythonClasses.ExpressionNode = NULL;
  pythonClasses.ExpressionNodeValue = NULL;
  pythonClasses.ExpressionNodeIdentifier = NULL;
  pythonClasses.ExpressionNodeBinary = NULL;
  pythonClasses.ExpressionNodeMult = NULL;
  pythonClasses.ExpressionNodeDiv = NULL;
  pythonClasses.ExpressionNodeAdd = NULL;
  pythonClasses.ExpressionNodeSubtract = NULL;
  pythonClasses.ExpressionNodeLower = NULL;
  pythonClasses.ExpressionNodeLowerEqual = NULL;
  pythonClasses.ExpressionNodeGreater = NULL;
  pythonClasses.ExpressionNodeGreaterEqual = NULL;
  pythonClasses.ExpressionNodeEqual = NULL;
  pythonClasses.ExpressionNodeUnequal = NULL;
  pythonClasses.ExpressionNodeNot = NULL;
  pythonClasses.ExpressionNodeAnd = NULL;
  pythonClasses.ExpressionNodeOr = NULL;
  pythonClasses.ExpressionNodeFunction = NULL;
  pythonClasses.ExpressionNodeUniformRandom = NULL;
  pythonClasses.ExpressionNodeGaussianRandom = NULL;
  pythonClasses.ExpressionNodePower = NULL;
  pythonClasses.ExpressionNodeLogarithm = NULL;
  pythonClasses.ExpressionNodeAtan = NULL;
  pythonClasses.PromoterElement = NULL;
  pythonClasses.PromoterElementLink = NULL;
  pythonClasses.PromoterElementConstitutive = NULL;
  pythonClasses.PromoterElementActivate = NULL;
  pythonClasses.PromoterElementRepress = NULL;
  pythonClasses.Factor = NULL;
  pythonClasses.Gene = NULL;
  pythonClasses.TranssysProgram = NULL;
  pythonClasses.GraphicsPrimitive = NULL;
  pythonClasses.GraphicsPrimitiveMove = NULL;
  pythonClasses.GraphicsPrimitivePush = NULL;
  pythonClasses.GraphicsPrimitivePop = NULL;
  pythonClasses.GraphicsPrimitiveTurn = NULL;
  pythonClasses.GraphicsPrimitiveRoll = NULL;
  pythonClasses.GraphicsPrimitiveBank = NULL;
  pythonClasses.GraphicsPrimitiveSphere = NULL;
  pythonClasses.GraphicsPrimitiveCylinder = NULL;
  pythonClasses.GraphicsPrimitiveBox = NULL;
  pythonClasses.GraphicsPrimitiveColor = NULL;
  pythonClasses.Symbol = NULL;
  pythonClasses.Assignment = NULL;
  pythonClasses.LhsSymbol = NULL;
  pythonClasses.ProductionElement = NULL;
  pythonClasses.Rule = NULL;
  pythonClasses.LsysProgram = NULL;
  pythonClasses.TranssysInstance = NULL;
  pythonClasses.SymbolInstance = NULL;
  pythonClasses.LsysSymbolString = NULL;
  pythonClasses.initialised = 0;
}


static void initExpressionNodetypeMap(const PYTHON_CLASSES *pc, EXPRESSION_NODETYPE_MAPENTRY map[])
{
  int i = 0;

  map[i].pythonClass = pc->ExpressionNodeValue;
  map[i++].node_type = NT_VALUE;
  map[i].pythonClass = pc->ExpressionNodeIdentifier;
  map[i++].node_type = NT_IDENTIFIER;
  map[i].pythonClass = pc->ExpressionNodeAdd;
  map[i++].node_type = NT_ADD;
  map[i].pythonClass = pc->ExpressionNodeSubtract;
  map[i++].node_type = NT_SUBTRACT;
  map[i].pythonClass = pc->ExpressionNodeMult;
  map[i++].node_type = NT_MULT;
  map[i].pythonClass = pc->ExpressionNodeDiv;
  map[i++].node_type = NT_DIV;
  map[i].pythonClass = pc->ExpressionNodeLower;
  map[i++].node_type = NT_LOWER;
  map[i].pythonClass = pc->ExpressionNodeLowerEqual;
  map[i++].node_type = NT_LOWER_EQUAL;
  map[i].pythonClass = pc->ExpressionNodeGreater;
  map[i++].node_type = NT_GREATER;
  map[i].pythonClass = pc->ExpressionNodeGreaterEqual;
  map[i++].node_type = NT_GREATER_EQUAL;
  map[i].pythonClass = pc->ExpressionNodeEqual;
  map[i++].node_type = NT_EQUAL;
  map[i].pythonClass = pc->ExpressionNodeUnequal;
  map[i++].node_type = NT_UNEQUAL;
  map[i].pythonClass = pc->ExpressionNodeNot;
  map[i++].node_type = NT_NOT;
  map[i].pythonClass = pc->ExpressionNodeAnd;
  map[i++].node_type = NT_LOGICAL_AND;
  map[i].pythonClass = pc->ExpressionNodeOr;
  map[i++].node_type = NT_LOGICAL_OR;
  map[i].pythonClass = pc->ExpressionNodeUniformRandom;
  map[i++].node_type = NT_RANDOM;
  map[i].pythonClass = pc->ExpressionNodeGaussianRandom;
  map[i++].node_type = NT_GAUSS;
  map[i].pythonClass = pc->ExpressionNodePower;
  map[i++].node_type = NT_POW;
  map[i].pythonClass = pc->ExpressionNodeLogarithm;
  map[i++].node_type = NT_LOG;
  map[i].pythonClass = pc->ExpressionNodeAtan;
  map[i++].node_type = NT_ATAN;
  map[i].pythonClass = NULL;
  map[i++].node_type = NT_NONE;
}


static void initGraphicsPrimitivetypeMap(const PYTHON_CLASSES *pc, GRAPHICS_PRIMITIVETYPE_MAPENTRY map[])
{
  int i = 0;

  map[i].pythonClass = pc->GraphicsPrimitiveMove;
  map[i++].graphics_primitive_type = GRAPHICS_MOVE;
  map[i].pythonClass = pc->GraphicsPrimitivePush;
  map[i++].graphics_primitive_type = GRAPHICS_PUSH;
  map[i].pythonClass = pc->GraphicsPrimitivePop;
  map[i++].graphics_primitive_type = GRAPHICS_POP;
  map[i].pythonClass = pc->GraphicsPrimitiveTurn;
  map[i++].graphics_primitive_type = GRAPHICS_TURN;
  map[i].pythonClass = pc->GraphicsPrimitiveRoll;
  map[i++].graphics_primitive_type = GRAPHICS_ROLL;
  map[i].pythonClass = pc->GraphicsPrimitiveBank;
  map[i++].graphics_primitive_type = GRAPHICS_BANK;
  map[i].pythonClass = pc->GraphicsPrimitiveSphere;
  map[i++].graphics_primitive_type = GRAPHICS_SPHERE;
  map[i].pythonClass = pc->GraphicsPrimitiveCylinder;
  map[i++].graphics_primitive_type = GRAPHICS_CYLINDER;
  map[i].pythonClass = pc->GraphicsPrimitiveBox;
  map[i++].graphics_primitive_type = GRAPHICS_BOX;
  map[i].pythonClass = pc->GraphicsPrimitiveColor;
  map[i++].graphics_primitive_type = GRAPHICS_COLOR;
  map[i].pythonClass = NULL;
  map[i++].graphics_primitive_type = GRAPHICS_NONE;
}


/*
 * Initialise the pythonClass struct
 */
static int initPythonClasses(void)
{
  if (pythonClasses.initialised)
  {
    return (0);
  }
  pythonClasses.ExpressionNode = getPythonClass("transsys", "ExpressionNode");
  if (pythonClasses.ExpressionNode == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeValue = getPythonClass("transsys", "ExpressionNodeValue");
  if (pythonClasses.ExpressionNodeValue == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeIdentifier = getPythonClass("transsys", "ExpressionNodeIdentifier");
  if (pythonClasses.ExpressionNodeIdentifier == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeBinary = getPythonClass("transsys", "ExpressionNodeBinary");
  if (pythonClasses.ExpressionNodeBinary == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeMult = getPythonClass("transsys", "ExpressionNodeMult");
  if (pythonClasses.ExpressionNodeMult == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeDiv = getPythonClass("transsys", "ExpressionNodeDiv");
  if (pythonClasses.ExpressionNodeDiv == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeAdd = getPythonClass("transsys", "ExpressionNodeAdd");
  if (pythonClasses.ExpressionNodeAdd == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeSubtract = getPythonClass("transsys", "ExpressionNodeSubtract");
  if (pythonClasses.ExpressionNodeSubtract == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeLower = getPythonClass("transsys", "ExpressionNodeLower");
  if (pythonClasses.ExpressionNodeLower == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeLowerEqual = getPythonClass("transsys", "ExpressionNodeLowerEqual");
  if (pythonClasses.ExpressionNodeLowerEqual == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeGreater = getPythonClass("transsys", "ExpressionNodeGreater");
  if (pythonClasses.ExpressionNodeGreater == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeGreaterEqual = getPythonClass("transsys", "ExpressionNodeGreaterEqual");
  if (pythonClasses.ExpressionNodeGreaterEqual == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeEqual = getPythonClass("transsys", "ExpressionNodeEqual");
  if (pythonClasses.ExpressionNodeEqual == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeUnequal = getPythonClass("transsys", "ExpressionNodeUnequal");
  if (pythonClasses.ExpressionNodeUnequal == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeNot = getPythonClass("transsys", "ExpressionNodeNot");
  if (pythonClasses.ExpressionNodeNot == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeAnd = getPythonClass("transsys", "ExpressionNodeAnd");
  if (pythonClasses.ExpressionNodeAnd == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeOr = getPythonClass("transsys", "ExpressionNodeOr");
  if (pythonClasses.ExpressionNodeOr == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeFunction = getPythonClass("transsys", "ExpressionNodeFunction");
  if (pythonClasses.ExpressionNodeFunction == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeUniformRandom = getPythonClass("transsys", "ExpressionNodeUniformRandom");
  if (pythonClasses.ExpressionNodeUniformRandom == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeGaussianRandom = getPythonClass("transsys", "ExpressionNodeGaussianRandom");
  if (pythonClasses.ExpressionNodeGaussianRandom == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodePower = getPythonClass("transsys", "ExpressionNodePower");
  if (pythonClasses.ExpressionNodePower == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeLogarithm = getPythonClass("transsys", "ExpressionNodeLogarithm");
  if (pythonClasses.ExpressionNodeLogarithm == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeAtan = getPythonClass("transsys", "ExpressionNodeAtan");
  if (pythonClasses.ExpressionNodeAtan == NULL)
  {
    return (-1);
  }
  pythonClasses.PromoterElement = getPythonClass("transsys", "PromoterElement");
  if (pythonClasses.PromoterElement == NULL)
  {
    return (-1);
  }
  pythonClasses.PromoterElementLink = getPythonClass("transsys", "PromoterElementLink");
  if (pythonClasses.PromoterElementLink == NULL)
  {
    return (-1);
  }
  pythonClasses.PromoterElementConstitutive = getPythonClass("transsys", "PromoterElementConstitutive");
  if (pythonClasses.PromoterElementConstitutive == NULL)
  {
    return (-1);
  }
  pythonClasses.PromoterElementActivate = getPythonClass("transsys", "PromoterElementActivate");
  if (pythonClasses.PromoterElementActivate == NULL)
  {
    return (-1);
  }
  pythonClasses.PromoterElementRepress = getPythonClass("transsys", "PromoterElementRepress");
  if (pythonClasses.PromoterElementRepress == NULL)
  {
    return (-1);
  }
  pythonClasses.Factor = getPythonClass("transsys", "Factor");
  if (pythonClasses.Factor == NULL)
  {
    return (-1);
  }
  pythonClasses.Gene = getPythonClass("transsys", "Gene");
  if (pythonClasses.Gene == NULL)
  {
    return (-1);
  }
  pythonClasses.TranssysProgram = getPythonClass("transsys", "TranssysProgram");
  if (pythonClasses.TranssysProgram == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitive = getPythonClass("transsys", "GraphicsPrimitive");
  if (pythonClasses.GraphicsPrimitive == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveMove = getPythonClass("transsys", "GraphicsPrimitiveMove");
  if (pythonClasses.GraphicsPrimitiveMove == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitivePush = getPythonClass("transsys", "GraphicsPrimitivePush");
  if (pythonClasses.GraphicsPrimitivePush == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitivePop = getPythonClass("transsys", "GraphicsPrimitivePop");
  if (pythonClasses.GraphicsPrimitivePop == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveTurn = getPythonClass("transsys", "GraphicsPrimitiveTurn");
  if (pythonClasses.GraphicsPrimitiveTurn == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveRoll = getPythonClass("transsys", "GraphicsPrimitiveRoll");
  if (pythonClasses.GraphicsPrimitiveRoll == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveBank = getPythonClass("transsys", "GraphicsPrimitiveBank");
  if (pythonClasses.GraphicsPrimitiveBank == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveSphere = getPythonClass("transsys", "GraphicsPrimitiveSphere");
  if (pythonClasses.GraphicsPrimitiveSphere == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveCylinder = getPythonClass("transsys", "GraphicsPrimitiveCylinder");
  if (pythonClasses.GraphicsPrimitiveCylinder == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveBox = getPythonClass("transsys", "GraphicsPrimitiveBox");
  if (pythonClasses.GraphicsPrimitiveBox == NULL)
  {
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveColor = getPythonClass("transsys", "GraphicsPrimitiveColor");
  if (pythonClasses.GraphicsPrimitiveColor == NULL)
  {
    return (-1);
  }
  pythonClasses.Symbol = getPythonClass("transsys", "Symbol");
  if (pythonClasses.Symbol == NULL)
  {
    return (-1);
  }
  pythonClasses.Assignment = getPythonClass("transsys", "Assignment");
  if (pythonClasses.Assignment == NULL)
  {
    return (-1);
  }
  pythonClasses.LhsSymbol = getPythonClass("transsys", "LhsSymbol");
  if (pythonClasses.LhsSymbol == NULL)
  {
    return (-1);
  }
  pythonClasses.ProductionElement = getPythonClass("transsys", "ProductionElement");
  if (pythonClasses.ProductionElement == NULL)
  {
    return (-1);
  }
  pythonClasses.Rule = getPythonClass("transsys", "Rule");
  if (pythonClasses.Rule == NULL)
  {
    return (-1);
  }
  pythonClasses.LsysProgram = getPythonClass("transsys", "LsysProgram");
  if (pythonClasses.LsysProgram == NULL)
  {
    return (-1);
  }
  pythonClasses.TranssysInstance = getPythonClass("transsys", "TranssysInstance");
  if (pythonClasses.TranssysInstance == NULL)
  {
    return (-1);
  }
  pythonClasses.SymbolInstance = getPythonClass("transsys", "SymbolInstance");
  if (pythonClasses.SymbolInstance == NULL)
  {
    return (-1);
  }
  pythonClasses.LsysSymbolString = getPythonClass("transsys", "LsysSymbolString");
  if (pythonClasses.LsysSymbolString == NULL)
  {
    return (-1);
  }
  initExpressionNodetypeMap(&pythonClasses, expressionNodetypeMap);
  initGraphicsPrimitivetypeMap(&pythonClasses, graphicsPrimitivetypeMap);
  pythonClasses.initialised = 1;
  return (0);
}


/*
 * Convenience function to get an attribute of string type as
 * a char *.
 * The char * returned is that belonging to the python string
 * object. Therefore, it can be expected to be valid only as long
 * as the calling function holds a reference to o, and it does
 * not alter o's attribute attr.
 * In practice, the returned char * should be used immediately e.g.
 * for creating C level transsys struct instances or for name
 * lookups etc., but it should not passed around or returned further, or
 * used after more extensive use of python library functions.
 */
static const char *getStringAttribute(PyObject *o, char *attr)
{
  PyObject *p = PyObject_GetAttrString(o, attr);
  const char *s = NULL;

  if (p == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "getStringAttribute: PyObject_GetAttrString failed for \"%s\"\n", attr);
    return (NULL);
  }
  if (!PyString_Check(p))
  {
    clib_message(CLIB_MSG_ERROR, "attribute \"%s\" is not a string\n", attr);
    PyErr_SetString(PyExc_TypeError, "getStringAttribute: attribute found but is not a string instance");
  }
  s = PyString_AsString(p);
  Py_DECREF(p);
  return (s);
}


/*
 * Checks whether object o is an instance of pyClass. If the check
 * fails, an exception is set and 0 is returned. The exception is
 * either set by PyObject_IsInstance or, if the check ran ok but
 * returns false, a TypeError exception with message msg is set.
 *
 * Callers can thus return with an error code (NULL etc.), and
 * without worrying about setting an exception.
 */
static int checkClass(PyObject *o, PyObject *pyClass, const char *msg)
{
  int chk = PyObject_IsInstance(o, pyClass);

  if (chk == -1)
  {
    return (0);
  }
  if (!chk)
  {
    PyErr_SetString(PyExc_TypeError, msg);
    return (0);
  }
  return (1);
}


/*
 * Checks whether o is None or an instance of pyClass. See checkClass
 * for further details.
 */
static int checkClassOrNone(PyObject *o, PyObject *pyClass, const char *msg)
{
  if (o == Py_None)
  {
    return (1);
  }
  return (checkClass(o, pyClass, msg));
}



static int initTranssysInstance(PyObject *python_ti, const TRANSSYS_INSTANCE *ti)
{
  long i;
  PyObject *python_fc = PyObject_GetAttrString(python_ti, "factor_concentration");

  if (python_fc == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "initTranssysInstance: PyObject_GetAttrString failed for \"factor_concentration\"\n");
    return (-1);
  }
  if (ti->transsys->num_factors != PyList_Size(python_fc))
  {
    PyErr_SetString(PyExc_TypeError, "initTranssysInstance: ti and python_ti num_factors mismatch");
  }
  for (i = 0; i < ti->transsys->num_factors; i++)
  {
    PyObject *python_c = PyFloat_FromDouble(ti->factor_concentration[i]);

    if (python_c == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "initTranssysInstance: PyFloat_FromDouble failed\n");
      Py_DECREF(python_fc);
      return (-1);
    }
    if (PyList_SetItem(python_fc, i, python_c) != 0)
    {
      clib_message(CLIB_MSG_TRACE, "initTranssysInstance: PyList_SetItem failed for index %d\n", i);
      Py_DECREF(python_fc);
      return (-1);
    }
  }
  Py_DECREF(python_fc);
  return (0);
}


static PyObject *newTranssysInstance(PyObject *python_transsys, const TRANSSYS_INSTANCE *ti)
{
  PyObject *a;
  PyObject *kw;
  PyObject *python_ti = NULL;

  a = PyTuple_New(1);
  if (a == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newTranssysInstance: PyTuple_New failed\n");
    return (NULL);
  }
  /* keep a reference to python_transsys (PyTuple_SetItem steals) */
  Py_INCREF(python_transsys);
  PyTuple_SetItem(a, 0, python_transsys);
  kw = PyDict_New();
  if (kw == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newTranssysInstance: PyDict_New failed\n");
    Py_DECREF(a);
    return (NULL);
  }
  python_ti = PyObject_Call(pythonClasses.TranssysInstance, a, kw);
  Py_DECREF(a);
  Py_DECREF(kw);
  if (python_ti == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newTranssysInstance: PyInstance_New failed\n");
    return (NULL);
  }
  if (ti != NULL)
  {
    if (initTranssysInstance(python_ti, ti) != 0)
    {
      clib_message(CLIB_MSG_TRACE, "newTranssysInstance: initTranssysInstance failed\n");
      Py_DECREF(python_ti);
      return (NULL);
    }
  }
  clib_message(CLIB_MSG_TRACE, "newTranssysInstance: returning %p\n", ti);
  return (python_ti);
}


/*
 * Construct a new Python SymbolInstance object.
 */
static PyObject *newSymbolInstance(PyObject *python_symbol, PyObject *python_ti, PyObject *python_rule)
{
  PyObject *a;
  PyObject *kw;
  PyObject *python_si = NULL;

  if (python_ti == NULL)
  {
    Py_INCREF(Py_None);
    python_ti = Py_None;
  }
  if (python_rule == NULL)
  {
    Py_INCREF(Py_None);
    python_rule = Py_None;
  }
  a = PyTuple_New(3);
  if (a == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newSymbolInstance: PyTuple_New failed\n");
    return (NULL);
  }
  /* keep references to python_symbol, python_ti, python_rule (PyTuple_SetItem steals) */
  Py_INCREF(python_symbol);
  PyTuple_SetItem(a, 0, python_symbol);
  Py_INCREF(python_ti);
  PyTuple_SetItem(a, 1, python_ti);
  Py_INCREF(python_rule);
  PyTuple_SetItem(a, 2, python_rule);
  kw = PyDict_New();
  if (kw == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newSymbolInstance: PyDict_New failed\n");
    Py_DECREF(a);
    return (NULL);
  }
  python_si = PyObject_Call(pythonClasses.SymbolInstance, a, kw);
  if (python_si == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newSymbolInstance: PyObject_Call failed\n");
  }
  Py_DECREF(a);
  Py_DECREF(kw);
  return (python_si);
}


static PyObject *newLsysSymbolString(PyObject *python_lsys)
{
  PyObject *a;
  PyObject *kw;
  PyObject *python_lstr = NULL;

  a = PyTuple_New(1);
  if (a == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newLsysSymbolString: PyTuple_New failed\n");
    return (NULL);
  }
  /* keep a reference to python_lsys (PyTuple_SetItem steals) */
  Py_INCREF(python_lsys);
  PyTuple_SetItem(a, 0, python_lsys);
  kw = PyDict_New();
  if (kw == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newLsysSymbolString: PyDict_New failed\n");
    Py_DECREF(a);
    return (NULL);
  }
  python_lstr = PyObject_Call(pythonClasses.LsysSymbolString, a, kw);
  if (python_lstr == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "newLsysSymbolString: PyObject_Call failed\n");
  }
  Py_DECREF(a);
  Py_DECREF(kw);
  return (python_lstr);
}


static EXPRESSION_NODE *extract_expression(PyObject *python_node);


static EXPR_NODE_TYPE identify_node_type(PyObject *python_node, EXPRESSION_NODETYPE_MAPENTRY map[])
{
  int i, return_value;

  for (i = 0; map[i].pythonClass != NULL; i++)
  {
    return_value = PyObject_IsInstance(python_node, map[i].pythonClass);
    if (return_value == -1)
    {
      clib_message(CLIB_MSG_TRACE, "identify_node_type: PyObject_IsInstance failed\n");
      return (NT_NONE);
    }
    if (return_value)
    {
      return (map[i].node_type);
    }
  }
  PyErr_SetString(PyExc_TypeError, "identify_node_type: unknown node type or unsuitable python object");
  return (NT_NONE);
}


static GRAPHICS_PRIMITIVE_TYPE identify_graphics_primitive_type(PyObject *python_gp, GRAPHICS_PRIMITIVETYPE_MAPENTRY map[])
{
  int i, return_value;

  for (i = 0; map[i].pythonClass != NULL; i++)
  {
    return_value = PyObject_IsInstance(python_gp, map[i].pythonClass);
    if (return_value == -1)
    {
      clib_message(CLIB_MSG_TRACE, "identify_graphics_primitive_type: PyObject_IsInstance failed\n");
      return (GRAPHICS_NONE);
    }
    if (return_value)
    {
      return (map[i].graphics_primitive_type);
    }
  }
  PyErr_SetString(PyExc_TypeError, "identify_node_type: unknown graphics primitive type or unsuitable python object");
  return (GRAPHICS_NONE);
}


static EXPRESSION_NODE *extract_expression_value(PyObject *python_node)
{
  PyObject *v_obj;
  double v;
  EXPRESSION_NODE *expression;

  clib_message(CLIB_MSG_TRACE, "extract_expression_value: start\n");
  if (!checkClass(python_node, pythonClasses.ExpressionNodeValue, "extract_expression_value: python_node is not an instance of ExpressionNodeValue"))
  {
    return (NULL);
  }
  v_obj = PyObject_GetAttrString(python_node, "value");
  if (v_obj == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_value: PyObject_GetAttrString failed for \"value\"\n");
    return (NULL);
  }
  v = PyFloat_AsDouble(v_obj);
  if (PyErr_Occurred() != NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_value: exception in PyFloat_AsDouble\n");
    Py_DECREF(v_obj);
    return (NULL);
  }
  Py_DECREF(v_obj);
  expression = new_expression_node(NT_VALUE, v);
  if (expression == NULL)
  {
    clib_message(CLIB_MSG_ERROR, "extract_expression_value: new_expression_node failed\n");
  }
  return (expression);
}


static EXPRESSION_NODE *extract_expression_identifier(PyObject *python_node)
{
/*
    self.factor = factor
    self.transsys_label = transsys_label
*/
  const char *factor_name, *transsys_label;
  PyObject *python_factor, *python_transsys_label;
  EXPRESSION_NODE *node;

  clib_message(CLIB_MSG_TRACE, "extract_expression_identifier: start\n");
  if (!checkClass(python_node, pythonClasses.ExpressionNodeIdentifier, "extract_expression_identifier: not an instance of ExpressionNodeIdentifier"))
  {
    return (NULL);
  }
  python_factor = PyObject_GetAttrString(python_node, "factor");
  if (python_factor == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_identifier: PyObject_GetAttrString failed for \"factor\"\n");
    return (NULL);
  }
  if (!checkClass(python_factor, pythonClasses.Factor, "extract_expression_identifier: factor attribute is not an instance of Factor"))
  {
    Py_DECREF(python_factor);
    return (NULL);
  }
  python_transsys_label = PyObject_GetAttrString(python_node, "transsys_label");
  if (python_transsys_label == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_identifier: PyObject_GetAttrString failed for \"transsys_label\"\n");
    Py_DECREF(python_factor);
    return (NULL);
  }
  if (python_transsys_label == Py_None)
  {
    transsys_label = NULL;
  }
  else
  {
    transsys_label = getStringAttribute(python_node, "transsys_label");
  }
  Py_DECREF(python_transsys_label);
  factor_name = getStringAttribute(python_factor, "name");
  node = new_expression_node(NT_RAW_IDENTIFIER, transsys_label, factor_name);
  Py_DECREF(python_factor);
  if (node == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_expression_identifier: new_expression_node failed");
  }
  return (node);
}


static EXPRESSION_NODE *extract_binary_expression(EXPR_NODE_TYPE node_type, PyObject *operand1, PyObject *operand2)
{
  EXPRESSION_NODE *arg1, *arg2, *expression;

  arg1 = extract_expression(operand1);
  if (arg1 == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_binary_expression: extract_expression failed for arg1\n");
    return (NULL);
  }
  arg2 = extract_expression(operand2);
  if (arg2 == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_binary_expression: extract_expression failed for arg2\n");
    free_expression_tree(arg1);
    return (NULL);
  }
  expression = new_expression_node(node_type, arg1, arg2);
  if (expression == NULL)
  {
    clib_message(CLIB_MSG_ERROR, "extract_binary_expression: new_expression_node failed\n");
  }
  return (expression);
}


static EXPRESSION_NODE *extract_expression_binary(PyObject *python_node, EXPR_NODE_TYPE node_type)
{
  PyObject *operand1, *operand2;
  EXPRESSION_NODE *expression;

  if (!checkClass(python_node, pythonClasses.ExpressionNodeBinary, "extract_expression_binary: python_node is not an instance of ExpressionNodeBinary"))
  {
    return (NULL);
  }
  operand1 = PyObject_GetAttrString(python_node, "operand1");
  if (operand1 == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_binary: PyObject_GetAttrString failed for \"operand1\"\n");
    return (NULL);
  }
  operand2 = PyObject_GetAttrString(python_node, "operand2");
  if (operand2 == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_binary: PyObject_GetAttrString failed for \"operand2\"\n");
    Py_DECREF(operand1);
    return (NULL);
  }
  expression = extract_binary_expression(node_type, operand1, operand2);
  Py_DECREF(operand1);
  Py_DECREF(operand2);
  return (expression);
}


/*
 * FIXME: this is good enough for unary and binary functions, but will
 * have to be implemented more properly if functions with higher arities
 * should be introduced.
 */
static EXPRESSION_NODE *extract_expression_function(PyObject *python_node, EXPR_NODE_TYPE node_type)
{
  PyObject *operand_list, *operand1, *operand2;
  EXPRESSION_NODE *arg, *expression = NULL;
  int list_size, arity;

  switch (node_type)
  {
  case NT_ATAN :
    arity = 1;
    break;
  case NT_RANDOM :
  case NT_GAUSS :
  case NT_POW :
  case NT_LOG :
    arity = 2;
    break;
  default :
    PyErr_SetString(PyExc_TypeError, "extract_expression_function: could not determine arity");
    return (NULL);
  }
  if (!checkClass(python_node, pythonClasses.ExpressionNodeFunction, "extract_expression_function: python_node is not an instance of ExpressionNodeFunction"))
  {
    return (NULL);
  }
  operand_list = PyObject_GetAttrString(python_node, "operand");
  if (operand_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_function: PyObject_GetAttrString failed for \"operand\"\n");
    return (NULL);
  }
  list_size = PyList_Size(operand_list);
  if (list_size != arity)
  {
    Py_DECREF(operand_list);
    PyErr_SetString(PyExc_TypeError, "extract_expression_function: incorrect arity");
    return (NULL);
  }
  operand1 = PyList_GetItem(operand_list, 0);
  if (operand1 == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_function: PyList_GetItem failed for index 0\n");
    Py_DECREF(operand_list);
    return (NULL);
  }
  Py_INCREF(operand1);
  if (arity == 1)
  {
    arg = extract_expression(operand1);
    if (arg == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_expression_function: extract_expression failed for unary arg\n");
      return (NULL);
    }
    expression = new_expression_node(node_type, arg);
    if (expression == NULL)
    {
      clib_message(CLIB_MSG_ERROR, "extract_expression_function: new_expression_node failed\n");
    }
  }
  else if (arity == 2)
  {
    operand2 = PyList_GetItem(operand_list, 1);
    if (operand2 == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_expression_function: PyList_GetItem failed for index 1\n");
      Py_DECREF(operand1);
      Py_DECREF(operand_list);
      return (NULL);
    }
    Py_INCREF(operand2);
    expression = extract_binary_expression(node_type, operand1, operand2);
    if (expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_expression_function: extract_binary_expression failed\n");
    }
    Py_DECREF(operand2);
  }
  else
  {
    PyErr_SetString(PyExc_StandardError, "extract_expression_function: cannot handle arity larger than 2");
  }
  Py_DECREF(operand_list);
  Py_DECREF(operand1);
  return (expression);
}


static EXPRESSION_NODE *extract_expression_not(PyObject *python_node)
{
  EXPRESSION_NODE *arg, *expression;
  PyObject *operand;

  if (!checkClass(python_node, pythonClasses.ExpressionNodeNot, "extract_expression_not: python_node is not an instance of ExpressionNodeNot"))
  {
    return (NULL);
  }
  operand = PyObject_GetAttrString(python_node, "operand");
  if (operand == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_not: PyObject_GetAttrString failed for \"operand\"\n");
    return (NULL);
  }
  arg = extract_expression(operand);
  Py_DECREF(operand);
  if (arg == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_not: extract_expression failed for arg\n");
    return (NULL);
  }
  expression = new_expression_node(NT_NOT, arg);
  if (expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression_not: new_expression_node failed\n");
  }
  return (expression);
}


static EXPRESSION_NODE *extract_expression(PyObject *python_node)
{
  EXPR_NODE_TYPE node_type;
  EXPRESSION_NODE *expression;

  clib_message(CLIB_MSG_TRACE, "extract_expression: starting\n");
  if (!checkClass(python_node, pythonClasses.ExpressionNode, "extract_expression: python_node is not an instance of ExpressionNode"))
  {
    return (NULL);
  }
  node_type = identify_node_type(python_node, expressionNodetypeMap);
  clib_message(CLIB_MSG_TRACE, "extract_expression: node type identified as %d\n", (int) node_type);
  if (node_type == NT_NONE)
  {
    clib_message(CLIB_MSG_TRACE, "extract_expression: identify_node_type failed\n");
    return (NULL);
  }
  switch (node_type)
  {
  case NT_VALUE :
    clib_message(CLIB_MSG_TRACE, "extract_expression: extracting value\n");
    expression = extract_expression_value(python_node);
    if (expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_expression: extract_expression_value failed\n");
    }
    return (expression);
  case NT_IDENTIFIER :
    clib_message(CLIB_MSG_TRACE, "extract_expression: extracting identifier\n");
    expression = extract_expression_identifier(python_node);
    if (expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_expression: extract_expression_identifier failed\n");
    }
    return (expression);
  case NT_MULT :
  case NT_DIV :
  case NT_ADD :
  case NT_SUBTRACT :
  case NT_LOWER :
  case NT_LOWER_EQUAL :
  case NT_GREATER :
  case NT_GREATER_EQUAL :
  case NT_EQUAL :
  case NT_UNEQUAL :
  case NT_LOGICAL_AND :
  case NT_LOGICAL_OR :
    clib_message(CLIB_MSG_TRACE, "extract_expression: extracting binary (%d)\n", (int) node_type);
    expression = extract_expression_binary(python_node, node_type);
    if (expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_expression: extract_expression_binary failed\n");
    }
    return (expression);
  case NT_RANDOM :
  case NT_GAUSS :
  case NT_POW :
  case NT_LOG :
  case NT_ATAN :
    clib_message(CLIB_MSG_TRACE, "extract_expression: extracting function (type %d)\n", (int) node_type);
    expression = extract_expression_function(python_node, node_type);
    if (expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_expression: extract_expression_function failed\n");
    }
    return (expression);
  case NT_NOT :
    clib_message(CLIB_MSG_TRACE, "extract_expression: extracting !\n");
    expression = extract_expression_not(python_node);
    if (expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_expression: extract_expression_not failed\n");
    }
    return (expression);
  default :
    clib_message(CLIB_MSG_ERROR, "extract_expression: unrecognised type of expression node\n");
    PyErr_SetString(PyExc_TypeError, "extract_expression: unrecognised type of expression node");
    /* set some exception indicating that no matching type was found */
    return (NULL);
  }
  PyErr_SetString(PyExc_SystemError, "extract_expression: fell through switch (half-implemented expression type?)");
  return (NULL);
}


static EXPRESSION_NODE *getExpressionString(PyObject *o, char *attr)
{
  PyObject *python_expression = PyObject_GetAttrString(o, attr);
  EXPRESSION_NODE *expression;

  if (python_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "getExpressionString: PyObject_GetAttrString failed for \"%s\"\n", attr);
    return (NULL);
  }
  expression = extract_expression(python_expression);
  Py_DECREF(python_expression);
  if (expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "getExpressionString: extract_expression failed for \"%s\"\n", attr);
  }
  return (expression);
}


static FACTOR_ELEMENT *extract_factor(PyObject *python_factor)
{
  PyObject *python_expression;
  EXPRESSION_NODE *decay_expression = NULL, *diffusibility_expression = NULL, *synthesis_expression = NULL;
  FACTOR_ELEMENT *factor = NULL;
  const char *name;

  clib_message(CLIB_MSG_TRACE, "extract_factor: getting decay\n");
  if (!checkClass(python_factor, pythonClasses.Factor, "extract_factor: python_factor is not an instance of Factor"))
  {
    return (NULL);
  }
  python_expression = PyObject_GetAttrString(python_factor, "decay_expression");
  if (python_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_factor: PyObject_GetAttrString failed for \"decay_expression\"\n");
    return (NULL);
  }
  if (python_expression == Py_None)
  {
    decay_expression = NULL;
  }
  else
  {
    decay_expression = extract_expression(python_expression);
    if (decay_expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_factor: extract_expression failed for \"decay_expression\"\n");
      Py_DECREF(python_expression);
      return (NULL);
    }
    clib_message(CLIB_MSG_TRACE, "extract_factor: extracted decay expression\n");
  }
  Py_DECREF(python_expression);
  clib_message(CLIB_MSG_TRACE, "extract_factor: got decay expression: %p\n", (void *) decay_expression);
  python_expression = PyObject_GetAttrString(python_factor, "diffusibility_expression");
  if (python_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_factor: PyObject_GetAttrString failed for \"diffusibility_expression\"\n");
    return (NULL);
  }
  if (python_expression == Py_None)
  {
    diffusibility_expression = NULL;
  }
  else
  {
    diffusibility_expression = extract_expression(python_expression);
    if (diffusibility_expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_factor: extract_expression failed for \"diffusibility_expression\"\n");
      free_expression_tree(decay_expression);
      Py_DECREF(python_expression);
      return (NULL);
    }
    clib_message(CLIB_MSG_TRACE, "extract_factor: extracted diffusibility expression\n");
  }
  Py_DECREF(python_expression);
  clib_message(CLIB_MSG_TRACE, "extract_factor: got diffusibility expression: %p\n", (void *) diffusibility_expression);


  python_expression = PyObject_GetAttrString(python_factor, "synthesis_expression");
  if (python_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_factor: PyObject_GetAttrString failed for \"synthesis_expression\"\n");
    return (NULL);
  }
  if (python_expression == Py_None)
  {
    synthesis_expression = NULL;
  }
  else
  {
    synthesis_expression = extract_expression(python_expression);
    if (synthesis_expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_factor: extract_expression failed for \"synthesis_expression\"\n");
      free_expression_tree(decay_expression);
      free_expression_tree(diffusibility_expression);
      Py_DECREF(python_expression);
      return (NULL);
    }
    clib_message(CLIB_MSG_TRACE, "extract_factor: extracted synthesis expression\n");
  }
  Py_DECREF(python_expression);
  clib_message(CLIB_MSG_TRACE, "extract_factor: got synthesis expression: %p\n", (void *) synthesis_expression);
  name = getStringAttribute(python_factor, "name");
  if (name == NULL)
  {
    clib_message(CLIB_MSG_ERROR, "extract_factor: getStringAttribute failed for \"name\"\n");
    free_expression_tree(decay_expression);
    free_expression_tree(diffusibility_expression);
    free_expression_tree(synthesis_expression);
    return (NULL);
  }
  factor = new_factor_element(name, decay_expression, diffusibility_expression, synthesis_expression);
  if (factor == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_factor: new_factor_element failed");
  }
  clib_message(CLIB_MSG_TRACE, "returning factor\n");
  return (factor);
}


static int extract_factor_elements(PyObject *python_tp, TRANSSYS *tp)
{
  int factor_index, num_factors;
  PyObject *factor_list, *python_factor;
  FACTOR_ELEMENT *factor_element;

  clib_message(CLIB_MSG_TRACE, "extract_factor_elements: starting\n");
  factor_list = PyObject_GetAttrString(python_tp, "factor_list");
  if (factor_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_factor_elements: PyObject_GetAttrString failed for \"factor_list\"\n");
    return (-1);
  }
  if (!PyList_Check(factor_list))
  {
    PyErr_SetString(PyExc_TypeError, "extract_factor_elements: factor_list is not a list");
    Py_DECREF(factor_list);
    return (-1);
  }
  num_factors = PyList_Size(factor_list);
  clib_message(CLIB_MSG_TRACE, "extract_factor_elements: %d factors\n", num_factors);
  for (factor_index = 0; factor_index < num_factors; factor_index++)
  {
    clib_message(CLIB_MSG_TRACE, "extracting factor #%d\n", factor_index);
    python_factor = PyList_GetItem(factor_list, factor_index);
    if (python_factor == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_factor_elements: PyList_GetItem failed for index %d\n", factor_index);
      Py_DECREF(factor_list);
      return (-1);
    }
    Py_INCREF(python_factor);
    factor_element = extract_factor(python_factor);
    if (factor_element)
    {
      clib_message(CLIB_MSG_TRACE, "got factor element \"%s\"\n", factor_element->name);
      add_factor_definition(tp, factor_element);
    }
    Py_DECREF(python_factor);
    if (factor_element == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_factor_elements: extract_factor failed\n");
      Py_DECREF(factor_list);
      return (-1);
    }
  }
  Py_DECREF(factor_list);
  return (0);
}


static PROMOTERELEMENT_TYPE identify_promoter_element(PyObject *python_pe)
{
  int return_value;

  return_value = PyObject_IsInstance(python_pe, pythonClasses.PromoterElementConstitutive);
  if (return_value == -1)
  {
    clib_message(CLIB_MSG_TRACE, "identify_promoter_element: PyObject_IsInstance failed for PromoterElementConstitutive\n");
    return (PROMOTERELEMENT_NONE);
  }
  if (return_value)
  {
    return (PROMOTERELEMENT_CONSTITUTIVE);
  }
  return_value = PyObject_IsInstance(python_pe, pythonClasses.PromoterElementActivate);
  if (return_value == -1)
  {
    clib_message(CLIB_MSG_TRACE, "identify_promoter_element: PyObject_IsInstance failed for PromoterElementActivate\n");
    return (PROMOTERELEMENT_NONE);
  }
  if (return_value)
  {
    return (PROMOTERELEMENT_ACTIVATE);
  }
  return_value = PyObject_IsInstance(python_pe, pythonClasses.PromoterElementRepress);
  if (return_value == -1)
  {
    clib_message(CLIB_MSG_TRACE, "identify_promoter_element: PyObject_IsInstance failed for PromoterElementRepress\n");
    return (PROMOTERELEMENT_NONE);
  }
  if (return_value)
  {
    return (PROMOTERELEMENT_REPRESS);
  }
  PyErr_SetString(PyExc_TypeError, "identify_promoter_element: unrecognised promoter element type");
  return (PROMOTERELEMENT_NONE);
}


static PROMOTER_ELEMENT *extract_promoterelement_constitutive(PyObject *python_pe)
{
  PROMOTER_ELEMENT *pe;
  EXPRESSION_NODE *expression;

  if (!checkClass(python_pe, pythonClasses.PromoterElementConstitutive, "extract_promoterelement_constitutive: python_pe is not an instance of PromoterElementConstitutive"))
  {
    return (NULL);
  }
  expression = getExpressionString(python_pe, "expression");
  if (expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_promoterelement_constitutive: getExpressionString failed\n");
    return (NULL);
  }
  pe = new_promoter_element(PROMOTERELEMENT_CONSTITUTIVE, 0, NULL, expression, NULL);
  if (pe == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_promoterelement_constitutive: new_promoter_element failed\n");
  }
  return (pe);
}


static int extract_factor_index(PyObject *python_factor, const TRANSSYS *tp)
{
  int factor_index;
  const char *name = NULL;

  if (!checkClass(python_factor, pythonClasses.Factor, "extract_factor_index: python_factor is not an instance of Factor"))
  {
    clib_message(CLIB_MSG_TRACE, "extract_factor_index: python_factor is not an instance of Factor\n");
    return (NO_INDEX);
  }
  name = getStringAttribute(python_factor, "name");
  if (name == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_factor_index: getStringAttribute failed\n");
    return (NO_INDEX);
  }
  factor_index = find_factor_index(tp, name);
  if (factor_index == NO_INDEX)
  {
    clib_message(CLIB_MSG_TRACE, "extract_factor_index: find_factor_index failed\n");
    PyErr_SetString(PyExc_StandardError, "extract_factor_index: find_factor_index failed");
  }
  return (factor_index);
}


static INTEGER_ARRAY *extract_factor_index_array(PyObject *factor_list, const TRANSSYS *tp)
{
  INTEGER_ARRAY *factor_index_array = NULL;
  int i, num_factors;

  if (!PyList_Check(factor_list))
  {
    PyErr_SetString(PyExc_TypeError, "extract_factor_index_array: factor_list is not a list");
    return (NULL);
  }
  num_factors = PyList_Size(factor_list);
  for (i = 0; i < num_factors; i++)
  {
    PyObject *factor = PyList_GetItem(factor_list, i);
    int factor_index;

    if (factor == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "find_factor_index_array: PyList_GetItem failed for index %d\n", i);
      free_integer_array(factor_index_array);
      return (NULL);
    }
    Py_INCREF(factor);
    factor_index = extract_factor_index(factor, tp);
    if (factor_index == NO_INDEX)
    {
      clib_message(CLIB_MSG_TRACE, "find_factor_index_array: extract_factor_index failed\n");
      free_integer_array(factor_index_array);
      clib_message(CLIB_MSG_TRACE, "find_factor_index_array: free_integer_array done\n");
      Py_DECREF(factor);
      return (NULL);
    }
    factor_index_array = extend_integer_array(factor_index_array, factor_index);
    if (factor_index_array == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "find_factor_index_array: extend_integer_array failed\n");
      free_integer_array(factor_index_array);
      Py_DECREF(factor);
      PyErr_SetString(PyExc_StandardError, "extend_integer_array failed");
      return (NULL);
    }
    Py_DECREF(factor);
  }
  return (factor_index_array);
}


static PROMOTER_ELEMENT *extract_promoterelement_link(PyObject *python_pe, PROMOTERELEMENT_TYPE pe_type, const TRANSSYS *tp)
{
  PyObject *factor_list;
  PROMOTER_ELEMENT *pe = NULL;
  INTEGER_ARRAY *factor_index_array;
  EXPRESSION_NODE *expression1, *expression2;

  if (!checkClass(python_pe, pythonClasses.PromoterElementLink, "extract_promoterelement_link: python_pe is not an instance of PromoterElementLink"))
  {
    return (NULL);
  }
  factor_list = PyObject_GetAttrString(python_pe, "factor_list");
  if (factor_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_promoterelement_link: PyObject_GetAttrString failed for \"factor_list\"\n");
    return (NULL);
  }
  factor_index_array = extract_factor_index_array(factor_list, tp);
  Py_DECREF(factor_list);
  if (factor_index_array == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_promoterelement_link: extract_factor_index_array failed\n");
    return (NULL);
  }
  expression1 = getExpressionString(python_pe, "expression1");
  if (expression1 == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_promoterelement_link: getExpressionString failed for \"expression1\"\n");
    return (NULL);
  }
  expression2 = getExpressionString(python_pe, "expression2");
  if (expression2 == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_promoterelement_link: getExpressionString failed for \"expression2\"\n");
    return (NULL);
  }
  pe = create_promoter(pe_type, factor_index_array, expression1, expression2);
  if (pe == NULL)
  {
    free_expression_tree(expression1);
    free_expression_tree(expression2);
    PyErr_SetString(PyExc_StandardError, "extract_promoterelement_link: create_promoter failed");
  }
  return (pe);
}


static PROMOTER_ELEMENT *extract_promoter_element(PyObject *python_pe, const TRANSSYS *tp)
{
  PROMOTER_ELEMENT *pe;
  PROMOTERELEMENT_TYPE pe_type = identify_promoter_element(python_pe);

  switch (pe_type)
  {
  case PROMOTERELEMENT_CONSTITUTIVE :
    pe = extract_promoterelement_constitutive(python_pe);
    if (pe == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_promoter_element: extract_promoterelement_constitutive failed\n");
    }
    return (pe);
  case PROMOTERELEMENT_ACTIVATE :
  case PROMOTERELEMENT_REPRESS :
    pe = extract_promoterelement_link(python_pe, pe_type, tp);
    if (pe == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_promoter_element: extract_promoterelement_link failed\n");
    }
    return (pe);
  default :
    clib_message(CLIB_MSG_TRACE, "extract_promoter_element: unidentified promoter element\n");
    /* identify_promoter_element has set the exception */
    return (NULL);
  }
}


static PROMOTER_ELEMENT *extract_promoter(PyObject *python_promoter, const TRANSSYS *tp)
{
  int num_promoter_elements, i;
  PROMOTER_ELEMENT *promoter = NULL, *pe;
  PyObject *python_pe;

  if (!PyList_Check(python_promoter))
  {
    PyErr_SetString(PyExc_TypeError, "extract_promoter: python_promoter is not a list");
    return (NULL);
  }
  num_promoter_elements = PyList_Size(python_promoter);
  for (i = 0; i < num_promoter_elements; i++)
  {
    python_pe = PyList_GetItem(python_promoter, i);
    if (python_pe == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_promoter: PyList_GetItem failed for index %d\n", i);
      return (NULL);
    }
    Py_INCREF(python_pe);
    pe = extract_promoter_element(python_pe, tp);
    Py_DECREF(python_pe);
    promoter = extend_promoter_list(promoter, pe);
    clib_message(CLIB_MSG_TRACE, "extract_promoter: list extended by element %d\n", i);
  }
  return (promoter);
}


static GENE_ELEMENT *extract_gene(PyObject *python_gene, const TRANSSYS *tp)
{
  GENE_ELEMENT *gene = NULL;
  PyObject *python_promoter, *python_product;
  PROMOTER_ELEMENT *promoter;
  int product_index;
  const char *name;

  clib_message(CLIB_MSG_TRACE, "extract_gene: starting\n");
  if (!checkClass(python_gene, pythonClasses.Gene, "extract_gene: python_gene is not an instance of Gene"))
  {
    return (NULL);
  }
  name = getStringAttribute(python_gene, "name");
  if (name == NULL)
  {
    clib_message(CLIB_MSG_ERROR, "extract_gene: getStringAttribute failed for \"name\"\n");
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "gene: \"%s\"\n", name);
  python_promoter = PyObject_GetAttrString(python_gene, "promoter");
  if (python_promoter == NULL)
  {
    clib_message(CLIB_MSG_ERROR, "extract_gene: PyObject_GetAttrString failed for \"promoter\"\n");
    return (NULL);
  }
  promoter = extract_promoter(python_promoter, tp);
  Py_DECREF(python_promoter);
  python_product = PyObject_GetAttrString(python_gene, "product");
  if (python_product == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_gene: PyObject_GetAttrString failed for \"product\"\n");
    free_promoter_list(promoter);
    return (NULL);
  }
  product_index = extract_factor_index(python_product, tp);
  Py_DECREF(python_product);
  if (product_index == NO_INDEX)
  {
    clib_message(CLIB_MSG_TRACE, "extract_gene: extract_factor_index failed\n");
    return (NULL);
  }
  gene = create_gene(name, promoter, product_index);
  if (gene == NULL)
  {
    PyErr_SetString(PyExc_StandardError, "extract_gene: create_gene failed");
    return (NULL);
  }
  return (gene);
}


static int extract_gene_elements(PyObject *python_tp, TRANSSYS *tp)
{
  int gene_index, num_genes;
  PyObject *gene_list, *python_gene;
  GENE_ELEMENT *gene_element;

  clib_message(CLIB_MSG_TRACE, "extract_gene_elements: starting\n");
  gene_list = PyObject_GetAttrString(python_tp, "gene_list");
  if (gene_list == NULL)
  {
    clib_message(CLIB_MSG_ERROR, "extract_gene_elements: PyObject_GetAttrString failed for \"gene_list\"\n");
    return (-1);
  }
  if (!PyList_Check(gene_list))
  {
    PyErr_SetString(PyExc_TypeError, "extract_gene_elements: gene_list attribute is not a list");
    return (-1);
  }
  num_genes = PyList_Size(gene_list);
  clib_message(CLIB_MSG_TRACE, "extract_gene_elements: %d genes\n", num_genes);
  for (gene_index = 0; gene_index < num_genes; gene_index++)
  {
    clib_message(CLIB_MSG_TRACE, "extracting gene #%d\n", gene_index);
    python_gene = PyList_GetItem(gene_list, gene_index);
    if (python_gene == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_gene_elements: PyList_GetItem failed for index %d\n", gene_index);
      Py_DECREF(gene_list);
      return (-1);
    }
    Py_INCREF(python_gene);
    gene_element = extract_gene(python_gene, tp);
    Py_DECREF(python_gene);
    if (gene_element == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_gene_elements: extract_gene failed\n");
      Py_DECREF(gene_list);
      return (-1);
    }
    clib_message(CLIB_MSG_TRACE, "got gene element \"%s\"\n", gene_element->name);
    add_gene_definition(tp, gene_element);
  }
  Py_DECREF(gene_list);
  clib_message(CLIB_MSG_TRACE, "extract_gene_elements: returning, arrayed = %d\n", tp->arrayed);
  return (0);
}


static TRANSSYS *extract_transsys(PyObject *python_tp)
{
  PyObject *tp_name;
  TRANSSYS *tp;

  /*
   * Notice that initPythonClasses must have been successfully called
   * before calling extract_transsys.
   */
  if (!checkClass(python_tp, pythonClasses.TranssysProgram, "extract_transsys: not an instance of TranssysProgram"))
  {
    return (NULL);
  }
  tp_name = PyObject_GetAttrString(python_tp, "name");
  if (tp_name == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_transsys: PyObject_GetAttrString failed for \"name\"\n");
    return (NULL);
  }
  if (!PyString_Check(tp_name))
  {
    PyErr_SetString(PyExc_TypeError, "extract_transsys: name is not a string");
    return (NULL);
  }
  tp = new_transsys(PyString_AsString(tp_name));
  Py_DECREF(tp_name);
  if (extract_factor_elements(python_tp, tp) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "extract_transsys: extract_factor_elements failed\n");
    free_transsys_list(tp);
    return (NULL);
  }
  if (extract_gene_elements(python_tp, tp) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "extract_transsys: extracting gene elements failed, freeing transsys\n");
    free_transsys_list(tp);
    return (NULL);
  }
  if (resolve_transsys(tp) != 0)
  {
    PyErr_SetString(PyExc_SystemError, "extract_transsys: resolve_transsys failed");
    free_transsys_list(tp);
    return (NULL);
  }
  return (tp);
}


static GRAPHICS_PRIMITIVE *extract_primitive_cylinder(PyObject *python_gp)
{
  EXPRESSION_NODE *diameter_expression, *length_expression;
  GRAPHICS_PRIMITIVE *gp;

  if (!checkClass(python_gp, pythonClasses.GraphicsPrimitiveCylinder, "extract_primitive_cylinder: python_gp is not an instance of GraphicsPrimitiveCylinder"))
  {
    return (NULL);
  }
  diameter_expression = getExpressionString(python_gp, "diameterExpression");
  if (diameter_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_primitive_cylinder: getExpressionString failed for \"diameterExpression\"\n");
    return (NULL);
  }
  length_expression = getExpressionString(python_gp, "lengthExpression");
  if (length_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_primitive_cylinder: getExpressionString failed for \"lengthExpression\"\n");
    free_expression_tree(diameter_expression);
    return (NULL);
  }
  gp = new_graphics_primitive(GRAPHICS_CYLINDER, diameter_expression, length_expression);
  if (gp == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_primitive_cylinder: new_graphics_primitive failed");
  }
  return (gp);
}


static GRAPHICS_PRIMITIVE *extract_primitive_box(PyObject *python_gp)
{
  EXPRESSION_NODE *x_expression, *y_expression, *z_expression;
  GRAPHICS_PRIMITIVE *gp;

  if (!checkClass(python_gp, pythonClasses.GraphicsPrimitiveBox, "extract_primitive_box: python_gp is not an instance of GraphicsPrimitiveBox"))
  {
    return (NULL);
  }
  x_expression = getExpressionString(python_gp, "xExpression");
  if (x_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_primitive_box: getExpressionString failed for \"xExpressoion\"\n");
    return (NULL);
  }
  y_expression = getExpressionString(python_gp, "yExpression");
  if (y_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_primitive_box: getExpressionString failed for \"yExpressoion\"\n");
    free_expression_tree(x_expression);
    return (NULL);
  }
  z_expression = getExpressionString(python_gp, "zExpression");
  if (z_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_primitive_box: getExpressionString failed for \"zExpressoion\"\n");
    free_expression_tree(x_expression);
    free_expression_tree(y_expression);
    return (NULL);
  }
  gp = new_graphics_primitive(GRAPHICS_BOX, x_expression, y_expression, z_expression);
  if (gp == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_primitive_box: new_graphics_primitive failed");
  }
  return (gp);
}


static GRAPHICS_PRIMITIVE *extract_primitive_color(PyObject *python_gp)
{
  EXPRESSION_NODE *red_expression, *green_expression, *blue_expression;
  GRAPHICS_PRIMITIVE *gp;

  if (!checkClass(python_gp, pythonClasses.GraphicsPrimitiveColor, "extract_primitive_color: python_gp is not an instance of GraphicsPrimitiveColor"))
  {
    return (NULL);
  }
  red_expression = getExpressionString(python_gp, "redExpression");
  if (red_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_primitive_color: getExpressionString failed for \"redExpression\"\n");
    return (NULL);
  }
  green_expression = getExpressionString(python_gp, "greenExpression");
  if (green_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_primitive_color: getExpressionString failed for \"greenExpression\"\n");
    free_expression_tree(red_expression);
    return (NULL);
  }
  blue_expression = getExpressionString(python_gp, "blueExpression");
  if (blue_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_primitive_color: getExpressionString failed for \"blueExpression\"\n");
    free_expression_tree(red_expression);
    free_expression_tree(green_expression);
    return (NULL);
  }
  gp = new_graphics_primitive(GRAPHICS_COLOR, red_expression, green_expression, blue_expression);
  if (gp == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_primitive_color: new_graphics_primitive failed");
  }
  return (gp);
}


static GRAPHICS_PRIMITIVE *extract_graphics_primitive(PyObject *python_gp)
{
  GRAPHICS_PRIMITIVE *gp;
  GRAPHICS_PRIMITIVE_TYPE gp_type;
  EXPRESSION_NODE *expression;

  if (!checkClass(python_gp, pythonClasses.GraphicsPrimitive, "extract_graphics_primitive: python_gp is not an instance of GraphicsPrimitive"))
  {
    return (NULL);
  }
  gp_type = identify_graphics_primitive_type(python_gp, graphicsPrimitivetypeMap);
  switch (gp_type)
  {
  case GRAPHICS_NONE :
    clib_message(CLIB_MSG_TRACE, "extract_graphics_primitive: identify_graphics_primitive_type failed\n");
    return (NULL);
  case GRAPHICS_PUSH :
  case GRAPHICS_POP :
    gp = new_graphics_primitive(gp_type);
    if (gp == NULL)
    {
      PyErr_SetString(PyExc_MemoryError, "extract_graphics_primitive: new_graphics_primitive failed for push / pop");
    }
    return (gp);
  case GRAPHICS_MOVE :
  case GRAPHICS_TURN :
  case GRAPHICS_ROLL :
  case GRAPHICS_BANK :
  case GRAPHICS_SPHERE :
    expression = getExpressionString(python_gp, "expression");
    if (expression == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_graphics_primitive: getExpressionString failed for \"expression\"\n");
      return (NULL);
    }
    gp = new_graphics_primitive(gp_type, expression);
    if (gp == NULL)
    {
      PyErr_SetString(PyExc_MemoryError, "extract_graphics_primitive: new_graphics_primitive failed");
      return (NULL);
    }
    return (gp);
  case GRAPHICS_CYLINDER :
    gp = extract_primitive_cylinder(python_gp);
    if (gp == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_graphics_primitive: extract_primitive_cylinder failed\n");
    }
    return (gp);
  case GRAPHICS_BOX :
    gp = extract_primitive_box(python_gp);
    if (gp == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_graphics_primitive: extract_primitive_box failed\n");
    }
    return (gp);
  case GRAPHICS_COLOR :
    gp =extract_primitive_color(python_gp);
    if (gp == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_graphics_primitive: extract_primitive_color failed\n");
    }
    return (gp);
  }
  PyErr_SetString(PyExc_SystemError, "extract_graphics_primitive: fell through switch (half-implemented primitive?)");
  return (NULL);
}


static int add_extracted_graphics_primitives(PyObject *python_symbol, SYMBOL_ELEMENT *symbol)
{
  PyObject *python_gplist, *python_gp;
  GRAPHICS_PRIMITIVE *gp_list = NULL, *gp;
  int num_primitives, i;

  python_gplist = PyObject_GetAttrString(python_symbol, "graphics");
  if (python_gplist == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "add_graphics_primitives: PyObject_GetAttrString failed for \"graphics\"\n");
    return (-1);
  }
  if (!PyList_Check(python_gplist))
  {
    Py_DECREF(python_gplist);
    PyErr_SetString(PyExc_TypeError, "add_graphics_primitives: graphics attribute is not a list");
    return (-1);
  }
  num_primitives = PyList_Size(python_gplist);
  for (i = 0; i < num_primitives; i++)
  {
    python_gp = PyList_GetItem(python_gplist, i);
    Py_INCREF(python_gp);
    gp = extract_graphics_primitive(python_gp);
    Py_DECREF(python_gp);
    if (gp == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "add_graphics_primitives: extract_graphics_primitive failed\n");
      Py_DECREF(python_gplist);
      free_graphics_primitive_list(gp_list);
      return (-1);
    }
    gp_list = add_graphics_primitive(gp_list, gp);
  }
  Py_DECREF(python_gplist);
  add_graphics_to_symbol(symbol, gp_list);
  return (0);
}


/*
 * Extract a SYMBOL from an instance of transsys.Symbol, and as a side
 * effect, extract the associated transsys program and append that
 * to transsys_list.
 * If the transsys program is already in transsys_list, it is not
 * extracted again; the pointer from the list is used. Identification
 * of the transsys program is done by its name.
 * Notice that it is therefore a severe bug to create different transsys
 * program instances with the same name on the python level.
 */
static SYMBOL_ELEMENT *extract_symbol(PyObject *python_symbol, TRANSSYS **transsys_list)
{
/*
    self.name = name
    self.transsys = transsys
    self.graphics = graphics
*/
  SYMBOL_ELEMENT *symbol;
  const TRANSSYS *tp;
  PyObject *python_transsys;
  const char *transsys_name, *symbol_name;

  if (!checkClass(python_symbol, pythonClasses.Symbol, "extract_symbol: python_symbol is not an instance of Symbol"))
  {
    return (NULL);
  }
  symbol_name = getStringAttribute(python_symbol, "name");
  if (symbol_name == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_symbol: getStringAttribute failed for \"name\"\n");
    return (NULL);
  }
  python_transsys = PyObject_GetAttrString(python_symbol, "transsys");
  if (python_transsys == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_symbol: PyObject_GetAttrString failed for \"transsys\"\n");
    return (NULL);
  }
  if (python_transsys == Py_None)
  {
    tp = NULL;
    clib_message(CLIB_MSG_TRACE, "extract_symbol: symbol \"%s\" has no transsys\n", symbol_name);
  }
  else
  {
    if (!checkClass(python_transsys, pythonClasses.TranssysProgram, "extract_symbol: transsys attribute is not a TranssysProgram instance"))
    {
      Py_DECREF(python_transsys);
      return (NULL);
    }
    transsys_name = getStringAttribute(python_transsys, "name");
    if (transsys_name == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_symbol: getStringAttribute failed for \"name\" (of transsys)\n");
      Py_DECREF(python_transsys);
      return (NULL);
    }
    clib_message(CLIB_MSG_TRACE, "extract_symbol: symbol \"%s\" has transsys \"%s\"\n", symbol_name, transsys_name);
    tp = find_transsys(*transsys_list, transsys_name);
    if (tp == NULL)
    {
      tp = extract_transsys(python_transsys);
      if (tp == NULL)
      {
	clib_message(CLIB_MSG_TRACE, "extract_symbol: extract_transsys failed\n");
	Py_DECREF(python_transsys);
	return (NULL);
      }
      *transsys_list = add_transsys(*transsys_list, (TRANSSYS *) tp);
      clib_message(CLIB_MSG_TRACE, "extract_symbol: add_transsys returned %p\n", (void *) *transsys_list);
      if (*transsys_list == NULL)
      {
	PyErr_SetString(PyExc_MemoryError, "extract_symbol: add_transsys failed");
	Py_DECREF(python_transsys);
	return (NULL);
      }
      clib_message(CLIB_MSG_TRACE, "extract_symbol: extracted and added transsys \"%s\"\n", transsys_name);
    }
    else
    {
      clib_message(CLIB_MSG_TRACE, "extract_symbol: transsys \"%s\" already in list\n", transsys_name);
    }
  }
  Py_DECREF(python_transsys);
  symbol = new_symbol_element(symbol_name, tp);
  if (add_extracted_graphics_primitives(python_symbol, symbol) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "extract_symbol: add_extracted_graphics_primitives failed\n");
    free_symbol_element_list(symbol);
    return (NULL);
  }
  return (symbol);
}


static int extract_symbol_elements(PyObject *python_lsys, LSYS *lsys, TRANSSYS **transsys_list)
{
  PyObject *symbol_list, *python_symbol;
  SYMBOL_ELEMENT *symbol;
  int num_symbols, symbol_index;

  symbol_list = PyObject_GetAttrString(python_lsys, "symbols");
  if (symbol_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_symbol_elements: PyObject_GetAttrString failed for \"symbols\"\n");
    return (-1);
  }
  if (!PyList_Check(symbol_list))
  {
    PyErr_SetString(PyExc_TypeError, "extract_symbol_elements: symbols attribute is not a list");
    return (-1);
  }
  num_symbols = PyList_Size(symbol_list);
  for (symbol_index = 0; symbol_index < num_symbols; symbol_index++)
  {
    python_symbol = PyList_GetItem(symbol_list, symbol_index);
    if (python_symbol == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_symbol_elements: PyList_GetItem failed for index %d\n", symbol_index);
      Py_DECREF(symbol_list);
      return (-1);
    }
    Py_INCREF(python_symbol);
    symbol = extract_symbol(python_symbol, transsys_list);
    Py_DECREF(python_symbol);
    if (symbol == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_symbol_elements: extract_symbol failed\n");
      Py_DECREF(symbol_list);
      return (-1);
    }
    add_symbol_definition(lsys, symbol);
  }
  Py_DECREF(symbol_list);
  return (0);
}


static int extract_diffusionrange(PyObject *python_lsys, LSYS *lsys)
{
  PyObject *python_diffusionrange;
  long diffusion_range;

  python_diffusionrange = PyObject_GetAttrString(python_lsys, "diffusionrange");
  if (python_diffusionrange == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_diffusionrange: PyObject_GetAttrString failed for \"diffusionrange\"\n");
    return (-1);
  }
  diffusion_range = PyInt_AsLong(python_diffusionrange);
  if (PyErr_Occurred() != NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_diffusionrange: exception in PyInt_AsLong\n");
    Py_DECREF(python_diffusionrange);
    return (-1);
  }
  Py_DECREF(python_diffusionrange);
  lsys->diffusion_range = diffusion_range;
  return (0);
}


RAW_ASSIGNMENT *extract_raw_assignment(PyObject *python_assignment)
{
/*
    self.transsys = transsys
    self.factor = factor
    self.expression = expression
*/
  PyObject *python_factor, *python_expression;
  RAW_ASSIGNMENT *ra;
  EXPRESSION_NODE *expression;
  const char *factor_name;

  if (!checkClass(python_assignment, pythonClasses.Assignment, "extract_raw_assignment: python_assignment is not an instance of Assignment"))
  {
    return (NULL);
  }
  python_factor = PyObject_GetAttrString(python_assignment, "factor");
  if (python_factor == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_raw_assignment: PyObject_GetAttrString failed for \"factor\"\n");
    return (NULL);
  }
  if (!checkClass(python_factor, pythonClasses.Factor, "extract_raw_assignment: factor attribute of assignment is not an instance of Factor"))
  {
    Py_DECREF(python_factor);
    return (NULL);
  }
  python_expression = PyObject_GetAttrString(python_assignment, "expression");
  if (python_expression == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_raw_assignment: PyObject_GetAttrString failed for \"expression\"\n");
    Py_DECREF(python_factor);
    return (NULL);
  }
  if (!checkClass(python_expression, pythonClasses.ExpressionNode, "extract_raw_assignment: expression attribute of assignment is not an instance of Expression"))
  {
    Py_DECREF(python_factor);
    Py_DECREF(python_expression);
    return (NULL);
  }
  factor_name = getStringAttribute(python_factor, "name");
  if (factor_name == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_raw_assignment: PyObject_GetAttrString failed for \"name\" (of factor)\n");
    Py_DECREF(python_factor);
    Py_DECREF(python_expression);
    return (NULL);
  }
  expression = extract_expression(python_expression);
  clib_message(CLIB_MSG_TRACE, "extract_raw_assignment: calling create_raw_assignment(factor_name = \"%s\")\n", factor_name);
  ra = create_raw_assignment(factor_name, NULL, expression);
  Py_DECREF(python_factor);
  Py_DECREF(python_expression);
  if (ra == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_raw_assignment: create_raw_assignment failed");
  }
  clib_message(CLIB_MSG_TRACE, "extract_raw_assignment: returning %p\n", (void *) ra);
  return (ra);
}


static PRODUCTION_ELEMENT *extract_production_element(PyObject *python_pe, LSYS *lsys, LHS_SYMBOL *lhs_symbol_list)
{
  PyObject *python_symbol, *python_template_label, *python_assignment_list, *python_assignment;
  PRODUCTION_ELEMENT *pe = NULL;
  int num_assignments, i;
  RAW_ASSIGNMENT *raw_assignment_list = NULL, *ra;

/*
    self.symbol = symbol
    self.template_label = template_label
    self.assignments = assignments
*/

  if (!checkClass(python_pe, pythonClasses.ProductionElement, "extract_production_element: not an instance of ProductionElement"))
  {
    return (NULL);
  }
  python_symbol = PyObject_GetAttrString(python_pe, "symbol");
  if (python_symbol == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_production_element: PyObject_GetAttrString failed for \"symbol\"\n");
    return (NULL);
  }
  python_template_label = PyObject_GetAttrString(python_pe, "template_label");
  if (python_template_label == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_production_element: PyObject_GetAttrString failed for \"template_label\"\n");
    Py_DECREF(python_symbol);
    return (NULL);
  }
  if (PyString_Check(python_template_label))
  {
    clib_message(CLIB_MSG_TRACE, "extract_production_element: raw assignment with template label \"%s\"\n", PyString_AsString(python_template_label));
    ra = create_raw_assignment(NULL, PyString_AsString(python_template_label), NULL);
    raw_assignment_list = add_raw_assignment(raw_assignment_list, ra);
  }
  else if (python_template_label != Py_None)
  {
    Py_DECREF(python_symbol);
    Py_DECREF(python_template_label);
    PyErr_SetString(PyExc_TypeError, "extract_production_element: template_label is neither string nor None");
    return (NULL);
  }
  Py_DECREF(python_template_label);
  python_assignment_list = PyObject_GetAttrString(python_pe, "assignments");
  if (python_assignment_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_production_element: PyObject_GetAttrString failed for \"assignments\"\n");
    Py_DECREF(python_symbol);
    return (NULL);
  }
  if (!PyList_Check(python_assignment_list))
  {
    Py_DECREF(python_symbol);
    Py_DECREF(python_assignment_list);
    PyErr_SetString(PyExc_TypeError, "extract_production_element: assignments attribute is not a list");
    return (NULL);
  }
  num_assignments = PyList_Size(python_assignment_list);
  for (i = 0; i < num_assignments; i++)
  {
    python_assignment = PyList_GetItem(python_assignment_list, i);
    Py_INCREF(python_assignment);
    ra = extract_raw_assignment(python_assignment);
    Py_DECREF(python_assignment);
    if (ra == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_production_element: extract_raw_assignment failed\n");
      Py_DECREF(python_symbol);
      Py_DECREF(python_assignment_list);
      return (NULL);
    }
    raw_assignment_list = add_raw_assignment(raw_assignment_list, ra);
    clib_message(CLIB_MSG_TRACE, "add_raw_assignment done, i = %d\n", i);
  }
  pe = create_production_element(getStringAttribute(python_symbol, "name"), raw_assignment_list, lsys, lhs_symbol_list);
  Py_DECREF(python_symbol);
  Py_DECREF(python_assignment_list);
  if (pe == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_production_element: create_production_element failed");
  }
  return (pe);
}


static SYMBOL_PRODUCTION *extract_symbol_production(PyObject *python_sp, LSYS *lsys, LHS_SYMBOL *lhs_symbol_list)
{
  PyObject *python_pe;
  SYMBOL_PRODUCTION *sp = NULL;
  PRODUCTION_ELEMENT *pe_list = NULL, *pe;
  int num_production_elements, i;

  if (!PyList_Check(python_sp))
  {
    PyErr_SetString(PyExc_TypeError, "extract_symbol_production: python_sp is not a list");
    return (NULL);
  }
  num_production_elements = PyList_Size(python_sp);
  for (i = 0; i < num_production_elements; i++)
  {
    python_pe = PyList_GetItem(python_sp, i);
    if (python_pe == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_symbol_production: PyList_GetItem failed for index %d\n", i);
      free_production_element_list(pe_list);
      return (NULL);
    }
    Py_INCREF(python_pe);
    pe = extract_production_element(python_pe, lsys, lhs_symbol_list);
    Py_DECREF(python_pe);
    if (pe == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_symbol_production: extract_production_element failed\n");
      free_production_element_list(pe_list);
      return (NULL);
    }
    pe_list = add_production_element(pe_list, pe);
    clib_message(CLIB_MSG_TRACE, "extract_symbol_production: add_production_element done\n");
  }
  sp = new_symbol_production(pe_list);
  if (sp == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_symbol_production: new_symbol_production failed");
  }
  return (sp);
}


static int extract_axiom(PyObject *python_lsys, LSYS *lsys)
{
  PyObject *python_axiom;
  SYMBOL_PRODUCTION *axiom;

  python_axiom = PyObject_GetAttrString(python_lsys, "axiom");
  if (python_axiom == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_axiom: PyObject_GetAttrString failed for \"axiom\"\n");
    return (-1);
  }
  axiom = extract_symbol_production(python_axiom, lsys, NULL);
  if (axiom == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_axiom: extract_symbol_production failed\n");
    Py_DECREF(python_axiom);
    return (-1);
  }
  Py_DECREF(python_axiom);
  add_axiom_definition(lsys, axiom);
  clib_message(CLIB_MSG_TRACE, "extract_axiom: added axiom\n");
  return (0);
}


static LHS_SYMBOL *extract_lhs_symbol(PyObject *python_lhs_symbol, LSYS *lsys)
{
/*
    self.symbol = symbol
    self.transsys_label = transsys_label
*/
  PyObject *python_symbol, *python_symbol_transsys, *python_transsys_label;
  LHS_SYMBOL *lhs_symbol;
  const char *symbol_name, *transsys_label = NULL;

  if (!checkClass(python_lhs_symbol, pythonClasses.LhsSymbol, "extract_lhs_symbol: python_lhs_symbol is not an instance of LhsSymbol"))
  {
    return (NULL);
  }
  python_symbol = PyObject_GetAttrString(python_lhs_symbol, "symbol");
  if (python_symbol == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lhs_symbol: PyObject_GetAttrString failed for \"symbol\"\n");
    return (NULL);
  }
  if (!checkClass(python_symbol, pythonClasses.Symbol, "extract_lhs_symbol: symbol is not an instance of Symbol"))
  {
    Py_DECREF(python_symbol);
    return (NULL);
  }
  symbol_name = getStringAttribute(python_symbol, "name");
  if (symbol_name == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lhs_symbol: getStringAttribute failed for \"name\" (of symbol)\n");
    Py_DECREF(python_symbol);
    return (NULL);
  }
  python_symbol_transsys = PyObject_GetAttrString(python_symbol, "transsys");
  if (python_symbol_transsys == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lhs_symbol: PyObject_GetAttrString failed for \"transsys\" (of symbol)\n");
    Py_DECREF(python_symbol);
    return (NULL);
  }
  if (!checkClassOrNone(python_symbol_transsys, pythonClasses.TranssysProgram, "extract_lhs_symbol: transsys attribute of symbol is neither None nor an instance of TranssysProgram"))
  {
    Py_DECREF(python_symbol);
    Py_DECREF(python_symbol_transsys);
    return (NULL);
  }
  python_transsys_label = PyObject_GetAttrString(python_lhs_symbol, "transsys_label");
  if (python_transsys_label == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lhs_symbol: PyObject_GetAttrString failed for \"transsys_label\"\n");
    Py_DECREF(python_symbol);
    Py_DECREF(python_symbol_transsys);
    return (NULL);
  }
  if (PyString_Check(python_transsys_label))
  {
    if (python_symbol_transsys == Py_None)
    {
      PyErr_SetString(PyExc_RuntimeError, "extract_lhs_symbol: transsys label given but symbol has no transsys");
      Py_DECREF(python_symbol);
      Py_DECREF(python_symbol_transsys);
      Py_DECREF(python_transsys_label);
      return (NULL);
    }
    transsys_label = PyString_AsString(python_transsys_label);
  }
  else if (python_transsys_label == Py_None)
  {
    if (python_symbol_transsys != Py_None)
    {
      PyErr_SetString(PyExc_RuntimeError, "extract_lhs_symbol: no transsys label given but symbol does have transsys");
      Py_DECREF(python_symbol);
      Py_DECREF(python_symbol_transsys);
      Py_DECREF(python_transsys_label);
      return (NULL);
    }
    transsys_label = NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "extract_lhs_symbol: transsys_label attribute is neither string nor None");
    Py_DECREF(python_symbol);
    Py_DECREF(python_symbol_transsys);
    Py_DECREF(python_transsys_label);
    return (NULL);
  }
  lhs_symbol = create_lhs_symbol(symbol_name, transsys_label, lsys);
  Py_DECREF(python_symbol);
  Py_DECREF(python_symbol_transsys);
  Py_DECREF(python_transsys_label);
  if (lhs_symbol == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_lhs_symbol: create_lhs_symbol failed");
  }
  return (lhs_symbol);
}


static LHS_DESCRIPTOR *extract_lhs_descriptor(PyObject *python_lhs, LSYS *lsys)
{
  PyObject *python_lhs_symbol;
  int num_lhs_symbols, i;
  LHS_SYMBOL *lhs_symbol_list = NULL, *lhs_symbol;
  LHS_DESCRIPTOR *lhs_descriptor;

  if (!PyList_Check(python_lhs))
  {
    PyErr_SetString(PyExc_TypeError, "extract_lhs_descriptor: python_lhs is not a list");
    return (NULL);
  }
  num_lhs_symbols = PyList_Size(python_lhs);
  if (num_lhs_symbols < 1)
  {
    PyErr_SetString(PyExc_RuntimeError, "extract_lhs_descriptor: lhs has 0 lhs symbols");
    return (NULL);
  }
  for (i = 0; i < num_lhs_symbols; i++)
  {
    python_lhs_symbol = PyList_GetItem(python_lhs, i);
    if (python_lhs_symbol == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_lhs_descriptor: PyList_GetItem failed for index %d\n", i);
      return (NULL);
    }
    Py_INCREF(python_lhs_symbol);
    lhs_symbol = extract_lhs_symbol(python_lhs_symbol, lsys);
    Py_DECREF(python_lhs_symbol);
    if (lhs_symbol == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_lhs_descriptor: extract_lhs_symbol failed\n");
      free_lhs_symbol_list(lhs_symbol_list);
      return (NULL);
    }
    lhs_symbol_list = add_lhs_symbol(lhs_symbol_list, lhs_symbol);
    clib_message(CLIB_MSG_TRACE, "extract_lhs_descriptor: added lhs_symbol = %p (%d / %d)\n", (void *) lhs_symbol, i, num_lhs_symbols);
  }
  lhs_descriptor = create_lhs_descriptor(lsys, lhs_symbol_list);
  if (lhs_descriptor == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_lhs_descriptor: create_lhs_descriptor failed");
  }
  return (lhs_descriptor);
}


static RULE_ELEMENT *extract_rule(PyObject *python_rule, LSYS *lsys)
{
/*
    self.name = name
    self.lhs = lhs
    self.condition = condition
    self.rhs = rhs
*/
  PyObject *python_lhs, *python_condition, *python_rhs;
  const char *name;
  RULE_ELEMENT *rule;
  LHS_DESCRIPTOR *lhs_descriptor;

  if (!checkClass(python_rule, pythonClasses.Rule, "extract_rule: python_rule is not an instance of Rule"))
  {
    return (NULL);
  }
  name = getStringAttribute(python_rule, "name");
  if (name == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_rule: getStringAttribute failed for \"name\"\n");
    return (NULL);
  }
  rule = new_rule_element(name, NULL, NULL, NULL);
  if (rule == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "extract_rule: new_rule_element failed");
  }
  python_lhs = PyObject_GetAttrString(python_rule, "lhs");
  if (python_lhs == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_rule: PyObject_GetAttrString failed for \"lhs\"\n");
    free_rule_element_list(rule);
    return (NULL);
  }
  lhs_descriptor = extract_lhs_descriptor(python_lhs, lsys);
  Py_DECREF(python_lhs);
  if (lhs_descriptor == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_rule: extract_lhs_descriptor failed\n");
    free_rule_element_list(rule);
    return (NULL);
  }
  rule->lhs = lhs_descriptor;
  python_condition = PyObject_GetAttrString(python_rule, "condition");
  if (python_condition == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_rule: PyObject_GetAttrString failed for \"condition\"\n");
    free_rule_element_list(rule);
    return (NULL);
  }
  if (python_condition != Py_None)
  {
    rule->condition = extract_expression(python_condition);
    if (rule->condition == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_rule: extract expression failed (for condition)\n");
      Py_DECREF(python_condition);
      free_rule_element_list(rule);
      return (NULL);
    }
  }
  Py_DECREF(python_condition);
  python_rhs = PyObject_GetAttrString(python_rule, "rhs");
  if (python_rhs == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_rule: PyObject_GetAttrString failed for \"rhs\"\n");
    free_rule_element_list(rule);
    return (NULL);
  }
  rule->rhs = extract_symbol_production(python_rhs, lsys, rule->lhs->symbol_list);
  Py_DECREF(python_rhs);
  if (rule->rhs == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_rule: extract_symbol_production failed\n");
    free_rule_element_list(rule);
    return (NULL);
  }
  if (arrange_symbol_production_arrays(rule->rhs) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "extract_rule: arrange_symbol_production_arrays failed\n");
    free_rule_element_list(rule);
    return (NULL);
  }
  return (rule);
}


static int extract_rule_elements(PyObject *python_lsys, LSYS *lsys)
{
  PyObject *python_rule_list, *python_rule;
  int num_rules, i;
  RULE_ELEMENT *rule;

  python_rule_list = PyObject_GetAttrString(python_lsys, "rules");
  if (python_rule_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_rule_elements: PyObject_GetAttrString failed for \"rules\"\n");
    return (-1);
  }
  if (!PyList_Check(python_rule_list))
  {
    PyErr_SetString(PyExc_TypeError, "extract_rule_elements: rules attribute is not a list");
    Py_DECREF(python_rule_list);
    return (-1);
  }
  num_rules = PyList_Size(python_rule_list);
  for (i = 0; i < num_rules; i++)
  {
    python_rule = PyList_GetItem(python_rule_list, i);
    Py_INCREF(python_rule);
    rule = extract_rule(python_rule, lsys);
    Py_DECREF(python_rule);
    if (rule == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_rule_elements: extract_rule failed\n");
      Py_DECREF(python_rule_list);
      return (-1);
    }
    clib_message(CLIB_MSG_TRACE, "extract_rule_elements: adding rule definition\n");
    if (add_rule_definition(lsys, rule) != 0)
    {
      clib_message(CLIB_MSG_TRACE, "extract_rule_elements: add_rule_definition failed\n");
      Py_DECREF(python_rule_list);
      return (-1);
    }
    clib_message(CLIB_MSG_TRACE, "extract_rule_elements: added rule %d / %d\n", i, num_rules);
  }
  Py_DECREF(python_rule_list);
  return (0);
}


static LSYS *extract_lsys(PyObject *python_lsys)
{
  const char *lsys_name;
  LSYS *lsys;
  TRANSSYS *transsys_list = NULL;

  clib_message(CLIB_MSG_TRACE, "extract_lsys: starting\n");
  /*
   * Notice that initPythonClasses must have been successfully called
   * before calling extract_transsys.
   */
  if (!checkClass(python_lsys, pythonClasses.LsysProgram, "extract_lsys: not an instance of LsysProgram"))
  {
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "extract_lsys: type check ok\n");
  lsys_name = getStringAttribute(python_lsys, "name");
  if (lsys_name == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lsys: getStringAttribute failed for \"name\"\n");
    return (NULL);
  }
  lsys = new_lsys(lsys_name);
  clib_message(CLIB_MSG_TRACE, "extract_lsys: new_lsys done\n");
  if (lsys == NULL)
  {
    PyErr_SetString(PyExc_TypeError, "extract_lsys: new_lsys failed");
    return (NULL);
  }
  if (extract_symbol_elements(python_lsys, lsys, &transsys_list) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lsys: extract_symbol_elements failed\n");
    free_lsys_with_transsys(lsys);
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "extract_lsys: extract_symbol_elements done\n");
  if (extract_diffusionrange(python_lsys, lsys) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lsys: extract_diffusionrange failed\n");
    free_lsys_with_transsys(lsys);
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "extract_lsys: extract_diffusionrange done\n");
  if (extract_axiom(python_lsys, lsys) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lsys: extract_axiom failed\n");
    free_lsys_with_transsys(lsys);
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "extract_lsys: extract_axiom done\n");
  if (extract_rule_elements(python_lsys, lsys) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "extract_lsys: extract_rule_elements failed\n");
    free_lsys_with_transsys(lsys);
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "extract_lsys: extract_rule_elements done\n");
  if (resolve_lsys(lsys) != 0)
  {
    PyErr_SetString(PyExc_SystemError, "extract_lsys: resolve_lsys failed");
    free_lsys_with_transsys(lsys);
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "extract_lsys: resolve_lsys succeeded\n");
  return (lsys);
}


/*
 * Notice that the python_ti and the ti must actually correspond
 * to the same transsys program on the python and the C level,
 * respectively.
 */
static int extract_initial_factor_concentrations(PyObject *python_ti, TRANSSYS_INSTANCE *ti)
{
  PyObject *python_fc_start;
  long num_factors, i;

  if (!checkClass(python_ti, pythonClasses.TranssysInstance, "extract_initial_factor_concentrations: python_ti is not an instance of TranssysInstance"))
  {
    return (-1);
  }
  python_fc_start = PyObject_GetAttrString(python_ti, "factor_concentration");
  if (python_fc_start == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "extract_initial_factor_concentrations: PyObject_GetAttrString failed for \"factor_concentration\"\n");
    return (-1);
  }
  num_factors = PyList_Size(python_fc_start);
  if (num_factors != ti->transsys->num_factors)
  {
    clib_message(CLIB_MSG_ERROR, "extract_initial_factor_concentrations: python_ti and ti num_factors mismatch");
    Py_DECREF(python_fc_start);
    return (-1);
  }
  for (i = 0; i < num_factors; i++)
  {
    PyObject *python_c = PyList_GetItem(python_fc_start, i);

    if (python_c == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_initial_factor_concentration: PyList_GetItem failed for index %d\n", i);
      Py_DECREF(python_fc_start);
      return (-1);
    }
    Py_INCREF(python_c);
    ti->factor_concentration[i] = PyFloat_AsDouble(python_c);
    if (PyErr_Occurred() != NULL)
    {
      clib_message(CLIB_MSG_TRACE, "extract_initial_factor_concentration: exception in PyFloat_AsDouble\n");
      Py_DECREF(python_c);
      return (-1);
    }
    Py_DECREF(python_c);
  }
  Py_DECREF(python_fc_start);
  return (0);
}


/*
 * Construct a time series of transsys instances
 */
static PyObject *transsysInstanceTimeSeries(PyObject *python_ti_start, long num_timesteps, long sampling_period)
{
  TRANSSYS *tp;
  TRANSSYS_INSTANCE *ti;
  PyObject *python_tp, *python_ti_list, *python_ts;
  long i;

  if (!checkClass(python_ti_start, pythonClasses.TranssysInstance, "transsysInstanceTimeSeries: python_ti_start is not an instance of TranssysInstance"))
  {
    return (NULL);
  }
  python_tp = PyObject_GetAttrString(python_ti_start, "transsys_program");
  if (python_tp == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: PyObject_GetAttrString failed for \"transsys_program\"\n");
    return (NULL);
  }
  tp = extract_transsys(python_tp);
  if (tp == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: extract_transsys failed\n");
    Py_DECREF(python_tp);
    return (NULL);
  }
  ti = new_transsys_instance(tp);
  if (ti == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: new_transsys_instance failed\n");
    free_transsys_list(tp);
    Py_DECREF(python_tp);
    return (NULL);
  }
  if (extract_initial_factor_concentrations(python_ti_start, ti) != 0)
  {
    clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: extract_initial_factor_concentrations failed\n");
    free_transsys_instance(ti);
    free_transsys_list(tp);
    Py_DECREF(python_tp);
    return (NULL);
  }
  python_ti_list = PyList_New(0);
  if (python_ti_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: PyList_New failed\n");
    free_transsys_instance(ti);
    free_transsys_list(tp);
    Py_DECREF(python_tp);
    return (NULL);
  }
  for (i = 0; i < num_timesteps; i++)
  {
    if (i % sampling_period == 0)
    {
      PyObject *python_ti;
      clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: time step %d\n", i);
      python_ti = newTranssysInstance(python_tp, ti);
      if (python_ti == NULL)
      {
	clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: newTranssysInstance failed\n");
	free_transsys_instance(ti);
	free_transsys_list(tp);
	Py_DECREF(python_tp);
	Py_DECREF(python_ti_list);
	return (NULL);
      }
      python_ts = PyInt_FromLong(i);
      if (PyObject_SetAttrString(python_ti, "timestep", python_ts) == -1)
      {
	clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: PyObject_SetAttrString failed for \"timestep\"\n");
	free_transsys_instance(ti);
	free_transsys_list(tp);
	Py_DECREF(python_tp);
	Py_DECREF(python_ti);
	Py_DECREF(python_ti_list);
	return (NULL);
      }
      Py_DECREF(python_ts);
      if (PyList_Append(python_ti_list, python_ti) == -1)
      {
	clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries: PyList_Append failed for timestep %d\n", i);
	free_transsys_instance(ti);
	free_transsys_list(tp);
	Py_DECREF(python_tp);
	Py_DECREF(python_ti);
	Py_DECREF(python_ti_list);
	return (NULL);
      }
      Py_DECREF(python_ti);
    }
    process_expression(ti);
  }
  free_transsys_instance(ti);
  free_transsys_list(tp);
  Py_DECREF(python_tp);
  return (python_ti_list);
}


static PyObject *clib_timeseries(PyObject *self, PyObject *args)
{
  PyObject *python_ti_start, *python_int, *python_ti_list;
  long num_timesteps, sampling_period;

  refcountDebug_init();
  if (initPythonClasses() != 0)
  {
    clib_message(CLIB_MSG_TRACE, "clib_timeseries: initPythonClasses failed\n");
    /* Error should be set by getPythonClass */
    refcountDebug_report("clib_timeseries returning exception", 0);
    return (NULL);
  }
  if (PyTuple_Size(args) != 3)
  {
    PyErr_SetString(PyExc_TypeError, "exactly 3 arguments required: self (transsys instance), num_timesteps, sampling_period");
    refcountDebug_report("clib_timeseries returning exception", 0);
    return (NULL);
  }
  python_int = PyTuple_GetItem(args, 1);
  if (python_int == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "clib_timeseries: PyTuple_GetItem failed for index 1 (num_timesteps)\n");
    refcountDebug_report("clib_timeseries returning exception", 0);
    return (NULL);
  }
  Py_INCREF(python_int);
  if (!PyInt_Check(python_int))
  {
    PyErr_SetString(PyExc_TypeError, "number of timesteps (arg 2) must be an int");
    Py_DECREF(python_int);
    refcountDebug_report("clib_timeseries returning exception", 0);
    return (NULL);
  }
  num_timesteps = PyInt_AsLong(python_int);
  Py_DECREF(python_int);
  python_int = PyTuple_GetItem(args, 2);
  if (python_int == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "clib_timeseries: PyTuple_GetItem failed for index 2 (sampling_period)\n");
    refcountDebug_report("clib_timeseries returning exception", 0);
    return (NULL);
  }
  Py_INCREF(python_int);
  if (!PyInt_Check(python_int))
  {
    PyErr_SetString(PyExc_TypeError, "sampling period (arg 3) must be an int");
    Py_DECREF(python_int);
    refcountDebug_report("clib_timeseries returning exception", 0);
    return (NULL);
  }
  sampling_period = PyInt_AsLong(python_int);
  Py_DECREF(python_int);
  python_ti_start = PyTuple_GetItem(args, 0);
  if (python_ti_start == NULL)
  {
    clib_message(CLIB_MSG_ERROR, "PyTuple_GetItem failed for index 0 (transsys instance)\n");
    refcountDebug_report("clib_timeseries returning exception", 0);
    return (NULL);
  }
  Py_INCREF(python_ti_start);
  clib_message(CLIB_MSG_TRACE, "num_timesteps = %ld, sampling_period = %ld\n", num_timesteps, sampling_period);
  python_ti_list = transsysInstanceTimeSeries(python_ti_start, num_timesteps, sampling_period);
  if (python_ti_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "transsysInstanceTimeSeries failed\n");
    return (NULL);
  }
  Py_DECREF(python_ti_start);
  if (PyErr_Occurred() != NULL)
  {
    clib_message(CLIB_MSG_ERROR, "clib_timeseries: exception raised at regular return\n");
    Py_DECREF(python_ti_list);
    refcountDebug_report("clib_timeseries returning unhandled exception", 0);
    return (NULL);
  }
  refcountDebug_report("clib_timeseries returning successfully", 0);
  return (python_ti_list);
}


static PyObject *make_SymbolInstance(PyObject *python_symbol, const TRANSSYS_INSTANCE *ti, PyObject *python_rule)
{
  PyObject *python_transsys, *python_ti, *python_si;
  const char *transsys_name;

  python_transsys = PyObject_GetAttrString(python_symbol, "transsys");
  if (python_transsys == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: PyObject_GetAttrString failed for \"transsys\"\n");
    return (NULL);
  }
  if (python_transsys == Py_None)
  {
    if (ti->transsys != NULL)
    {
      clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: incompatible: python_transsys is None, transsys is \"%s\"\n", ti->transsys->name);
      PyErr_SetString(PyExc_RuntimeError, "make_SymbolInstance: python_transsys is None but symbol has transsys");
      Py_DECREF(python_transsys);
      return (NULL);
    }
    Py_INCREF(Py_None);
    python_ti = Py_None;
  }
  else
  {
    transsys_name = getStringAttribute(python_transsys, "name");
    if (transsys_name == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: getStringAttribute failed on transsys for \"name\"\n");
      Py_DECREF(python_transsys);
      return (NULL);
    }
    if (strcmp(transsys_name, ti->transsys->name))
    {
      clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: incompatible: python_transsys \"%s\", transsys \"%s\"\n", transsys_name, ti->transsys->name);
      PyErr_SetString(PyExc_RuntimeError, "make_SymbolInstance: python_transsys and transsys incompatible");
      Py_DECREF(python_transsys);
      return (NULL);
    }
    python_ti = newTranssysInstance(python_transsys, ti);
    if (python_ti == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: newTranssysInstance failed\n");
      Py_DECREF(python_transsys);
      return (NULL);
    }
  }
  Py_DECREF(python_transsys);
  clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: calling newSymbolInstance\n");
  python_si = newSymbolInstance(python_symbol, python_ti, python_rule);
  Py_DECREF(python_ti);
  clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: got python_si\n");
  if (python_si == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: newSymbolInstance failed\n");
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "make_SymbolInstance: done\n");
  return (python_si);
}


static PyObject *make_LsysSymbolString(PyObject *python_lsys, LSYS_STRING *lstr)
{
  int i;
  PyObject *python_lstr;
  PyObject *python_symbol_list, *python_symbol;
  PyObject *python_rule_list, *python_rule;
  PyObject *python_si_list, *python_si;
  const char *lsys_name;

  lsys_name = getStringAttribute(python_lsys, "name");
  if (lsys_name == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "make_LsysSmbolString: getStringAttribute failed for \"name\"\n");
    return (NULL);
  }
  if (strcmp(lsys_name, lstr->lsys->name))
  {
    PyErr_SetString(PyExc_RuntimeError, "make_LsysSymbolString: python_lsys and lstr are incompatible");
    clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: incompatible: python lsys \"%s\", lstr->lsys \"%s\"\n", lsys_name, lstr->lsys->name);
    return (NULL);
  }
  python_symbol_list = PyObject_GetAttrString(python_lsys, "symbols");
  if (python_symbol_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: PyObject_GetAttrString failed for \"symbols\"\n");
    return (NULL);
  }
  python_rule_list = PyObject_GetAttrString(python_lsys, "rules");
  if (python_rule_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: PyObject_GetAttrString failed for \"rules\"\n");
    Py_DECREF(python_symbol_list);
    return (NULL);
  }
  python_lstr = newLsysSymbolString(python_lsys);
  if (python_lstr == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: newLsysSymbolString failed\n");
    Py_DECREF(python_symbol_list);
    Py_DECREF(python_rule_list);
    return (NULL);
  }
  python_si_list = PyObject_GetAttrString(python_lstr, "symbol_list");
  if (python_si_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: PyObject_GetAttriString failed on python_lstr for \"symbol_list\"\n");
    Py_DECREF(python_symbol_list);
    Py_DECREF(python_rule_list);
    return (NULL);
  }
  for (i = 0; i < lstr->num_symbols; i++)
  {
    python_symbol = PyList_GetItem(python_symbol_list, lstr->symbol[i].symbol_index);
    if (python_symbol == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: PyList_GetItem failed on symbol list for index %d\n", lstr->symbol[i].symbol_index);
      Py_DECREF(python_symbol_list);
      Py_DECREF(python_rule_list);
      Py_DECREF(python_si_list);
      Py_DECREF(python_lstr);
      return (NULL);
    }
    Py_INCREF(python_symbol);
    if (lstr->symbol[i].rule_index == NO_INDEX)
    {
      Py_INCREF(Py_None);
      python_rule = Py_None;
    }
    else
    {
      clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: rule index is %d\n", lstr->symbol[i].rule_index);
      python_rule = PyList_GetItem(python_rule_list, lstr->symbol[i].rule_index);
      if (python_rule == NULL)
      {
	clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: PyList_GetItem failed on rule list for index %d\n", lstr->symbol[i].rule_index);
	Py_DECREF(python_symbol_list);
	Py_DECREF(python_symbol);
	Py_DECREF(python_rule_list);
	Py_DECREF(python_si_list);
	Py_DECREF(python_lstr);
	return (NULL);
      }
      Py_INCREF(python_rule);
      clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: got rule with index %d\n", lstr->symbol[i].rule_index);
    }
    python_si = make_SymbolInstance(python_symbol, &(lstr->symbol[i].transsys_instance), python_rule);
    clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: make_SymbolInstance done\n");
    Py_DECREF(python_symbol);
    Py_DECREF(python_rule);
    if (python_si == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: make_SymbolInstance failed\n");
      Py_DECREF(python_symbol_list);
      Py_DECREF(python_rule_list);
      Py_DECREF(python_si_list);
      Py_DECREF(python_lstr);
      return (NULL);
    }
    if (PyList_Append(python_si_list, python_si) == -1)
    {
      clib_message(CLIB_MSG_TRACE, "make_LsysSymbolString: PyList_Append failed for symbol #%d\n", i);
      Py_DECREF(python_symbol_list);
      Py_DECREF(python_rule_list);
      Py_DECREF(python_lstr);
      Py_DECREF(python_si_list);
      Py_DECREF(python_si);
      return (NULL);
    }
    Py_DECREF(python_si);
  }
  Py_DECREF(python_si_list);
  Py_DECREF(python_symbol_list);
  Py_DECREF(python_rule_list);
  return (python_lstr);
}


static PyObject *lsysStringSeries(PyObject *python_lsys, int num_timesteps, int sampling_period)
{
  PyObject *python_lstring_list, *python_lstr, *python_timestep;
  LSYS *lsys;
  LSYS_STRING *lstr, *lstr_next;
  int t;

  python_lstring_list = PyList_New(0);
  clib_message(CLIB_MSG_TRACE, "lsysStringSeries: python_lstring_list = %p\n", (void *) python_lstring_list);
  if (python_lstring_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "lsysStringSeries: PyList_New failed\n");
    return (NULL);
  }
  lsys = extract_lsys(python_lsys);
  if (lsys == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "lsysStringSeries: extract_lsys failed\n");
    Py_DECREF(python_lstring_list);
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "lsysStringSeries: extract_lsys succeeded\n");
  lstr = axiom_string(lsys);
  if (lstr == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "lsysStringSeries: axiom_string failed\n");
    Py_DECREF(python_lstring_list);
    free_lsys_with_transsys(lsys);
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "lsysStringSeries: axiom_string succeeded\n");
  for (t = 0; t < num_timesteps; t++)
  {
    clib_message(CLIB_MSG_TRACE, "lsysStringSeries: time step %d, %lu symbols\n", t, (unsigned long) lstr->num_symbols);
    lsys_string_expression(lstr);
    lsys_string_diffusion(lstr);
    lstr_next = derived_string(lstr);
    if (lstr_next == NULL)
    {
      clib_message(CLIB_MSG_TRACE, "lsysStringSeries: derived_string failed\n");
      Py_DECREF(python_lstring_list);
      free_lsys_string(lstr);
      free_lsys_with_transsys(lsys);
      return (NULL);
    }
    if ((t % sampling_period) == 0)
    {
      python_lstr = make_LsysSymbolString(python_lsys, lstr);
      if (python_lstr == NULL)
      {
	clib_message(CLIB_MSG_TRACE, "lsysStringSeries: make_LsysSymbolString failed\n");
	Py_DECREF(python_lstring_list);
	free_lsys_string(lstr);
	free_lsys_with_transsys(lsys);
	return (NULL);
      }
      python_timestep = PyInt_FromLong(t);
      if (python_timestep == NULL)
      {
	clib_message(CLIB_MSG_TRACE, "lsysStringSeries: PyInt_FromLong failed\n");
	Py_DECREF(python_lstring_list);
	Py_DECREF(python_lstr);
	free_lsys_string(lstr);
	free_lsys_with_transsys(lsys);
	return (NULL);
      }
      if (PyObject_SetAttrString(python_lstr, "timestep", python_timestep) == -1)
      {
	clib_message(CLIB_MSG_TRACE, "lsysStringSeries: PyObject_SetAttrString failed \"timestep\"\n");
	Py_DECREF(python_lstring_list);
	Py_DECREF(python_timestep);
	Py_DECREF(python_lstr);
	free_lsys_string(lstr);
	free_lsys_with_transsys(lsys);
	return (NULL);
      }
      Py_DECREF(python_timestep);
      if (PyList_Append(python_lstring_list, python_lstr) != 0)
      {
	clib_message(CLIB_MSG_TRACE, "lsysStringSeries: PyList_Append failed\n");
	Py_DECREF(python_lstring_list);
	Py_DECREF(python_lstr);
	free_lsys_string(lstr);
	free_lsys_with_transsys(lsys);
	return (NULL);
      }
      Py_DECREF(python_lstr);
    }
    free_lsys_string(lstr);
    lstr = lstr_next;
  }
  free_lsys_string(lstr);
  free_lsys_with_transsys(lsys);
  return (python_lstring_list);
}


static PyObject *clib_stringseries(PyObject *self, PyObject *args)
{
  PyObject *python_lsys, *python_int, *python_lstring_list;
  int num_timesteps, sampling_period;

  refcountDebug_report("clib_stringseries starting", 0);
  /* refcountDebug_init(); */
  if (initPythonClasses() != 0)
  {
    /* Error should be set by getPythonClass */
    refcountDebug_report("clib_stringseries returning exception", 0);
    return (NULL);
  }
  if (PyTuple_Size(args) != 3)
  {
    PyErr_SetString(PyExc_TypeError, "exactly 3 arguments required: lsys, num_timesteps, sampling_period");
    return (NULL);
  }
  python_lsys = PyTuple_GetItem(args, 0);
  if (python_lsys == NULL)
  {
    refcountDebug_report("clib_stringseries returning exception", 0);
    return (NULL);
  }
  Py_INCREF(python_lsys);
  python_int = PyTuple_GetItem(args, 1);
  if (python_int == NULL)
  {
    Py_DECREF(python_lsys);
    return (NULL);
  }
  Py_INCREF(python_int);
  if (!PyInt_Check(python_int))
  {
    PyErr_SetString(PyExc_TypeError, "argument 2 (num_timesteps) must be an int");
    Py_DECREF(python_lsys);
    Py_DECREF(python_int);
    return (NULL);
  }
  num_timesteps = PyInt_AsLong(python_int);
  Py_DECREF(python_int);
  python_int = PyTuple_GetItem(args, 2);
  if (python_int == NULL)
  {
    Py_DECREF(python_lsys);
    return (NULL);
  }
  Py_INCREF(python_int);
  if (!PyInt_Check(python_int))
  {
    PyErr_SetString(PyExc_TypeError, "argument 3 (sampling_period) must be an int");
    Py_DECREF(python_lsys);
    Py_DECREF(python_int);
    return (NULL);
  }
  sampling_period = PyInt_AsLong(python_int);
  Py_DECREF(python_int);
  python_lstring_list = lsysStringSeries(python_lsys, num_timesteps, sampling_period);
  Py_DECREF(python_lsys);
  if (python_lstring_list == NULL)
  {
    clib_message(CLIB_MSG_TRACE, "clib_stringseries: lsysStringSeries failed\n");
    return (NULL);
  }
  freePythonClasses();
  if (PyErr_Occurred() != NULL)
  {
    clib_message(CLIB_MSG_ERROR, "clib_stringseries: exception raised at regular return\n");
    Py_DECREF(python_lstring_list);
    refcountDebug_report("clib_stringseries returning unhandled exception", 0);
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "clib_stringseries: python_lstring_list = %p\n", (void *) python_lstring_list);
  refcountDebug_report("clib_stringseries returning successfully", 0);
  return (python_lstring_list);
}


static PyObject *clib_setverbose(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "i", &message_importance_threshold))
  {
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "transsys.clib: verbosity level set to %d\n", message_importance_threshold);
  Py_INCREF(Py_None);
  return (Py_None);
}


static PyObject *clib_srandom(PyObject *self, PyObject *args)
{
  unsigned long rndseed;
  if (!PyArg_ParseTuple(args, "k", &rndseed))
  {
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "clib_srandom: setting random seed %lu\n", rndseed);
  ulong_srandom(rndseed);
  Py_INCREF(Py_None);
  return (Py_None);
}


static PyObject *clib_dummy(PyObject *self, PyObject *args)
{
  Py_INCREF(Py_None);
  return (Py_None);
}


static PyMethodDef clib_methods[] = {
  {"timeseries", clib_timeseries, METH_VARARGS, "compute time series from a transsys instance"},
  {"stringseries", clib_stringseries, METH_VARARGS, "compute derivation series from a lsys program"},
  {"srandom", clib_srandom, METH_VARARGS, "set the random seed for clib transsys computations"},
  {"dummy", clib_dummy, METH_VARARGS, "dummy test function for clib development"},
  {"setverbose", clib_setverbose, METH_VARARGS, "set verbosity level for transsys.clib module"},
  {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initclib(void)
{
  PyObject *clib_module;
  clib_module = Py_InitModule("transsys.clib", clib_methods);
  /* FIXME: should not ignore return value */
  PyModule_AddStringConstant(clib_module, "clib_api_version", clib_api_version);
}

/* don't forget to change the clib_api_version */
