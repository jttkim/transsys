#include <Python.h>

#include <stdlib.h>
#include <stdio.h>

#include <transsys.h>

/*
 * FIXME: Determine which functions are to set exceptions and which
 * ones are to just pass them upwards.
 */

/*
 * Struct to contain pointers to Python classes that directly correspond
 * to C level classes.
 *
 * Note: The ExpressionNodeValue is used to detect whether the struct has
 * successfully been initialised.
 */
typedef struct
{
  PyObject *ExpressionNodeValue;
  PyObject *ExpressionNodeIdentifier;
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
  PyObject *ExpressionNodeUniformRandom;
  PyObject *ExpressionNodeGaussianRandom;
  PyObject *PromoterElementConstitutive;
  PyObject *PromoterElementActivate;
  PyObject *PromoterElementRepress;
  PyObject *Factor;
  PyObject *Gene;
  PyObject *TranssysProgram;
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
} PYTHON_CLASSES;


typedef struct
{
  PyObject *pythonClass;
  EXPR_NODE_TYPE node_type;
} EXPRESSION_NODETYPE_MAPENTRY;


static EXPRESSION_NODETYPE_MAPENTRY expressionNodetypeMap[] = {
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE},
  {NULL, NT_NONE}
};

static PYTHON_CLASSES pythonClasses = {
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
  CLIB_ERROR_MESSAGE,
  CLIB_ERROR_DEVMESSAGE,
  CLIB_ERROR_WARNING,
  CLIB_ERROR_ERROR,
  CLIB_ERROR_FATAL
} CLIB_ERROR_IMPORTANCE;


static CLIB_ERROR_IMPORTANCE error_importance_threshold = CLIB_ERROR_DEVMESSAGE;


static void clib_message(CLIB_ERROR_IMPORTANCE importance, const char *format, ...)
{
  va_list arglist;
  int imp = (int) importance;

  if (imp >= ((int) error_importance_threshold))
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
    clib_message(CLIB_ERROR_ERROR, "PyDict_GetItemString failed\n");
    return (NULL);
  }
  Py_INCREF(module);
  pythonClass = PyObject_GetAttrString(module, classname);
  Py_DECREF(module);
  return (pythonClass);
}


static void freePythonClasses(void)
{
  Py_DECREF(pythonClasses.ExpressionNodeValue);
  Py_DECREF(pythonClasses.ExpressionNodeIdentifier);
  Py_DECREF(pythonClasses.ExpressionNodeMult);
  Py_DECREF(pythonClasses.ExpressionNodeDiv);
  Py_DECREF(pythonClasses.ExpressionNodeAdd);
  Py_DECREF(pythonClasses.ExpressionNodeSubtract);
  Py_DECREF(pythonClasses.ExpressionNodeLower);
  Py_DECREF(pythonClasses.ExpressionNodeLowerEqual);
  Py_DECREF(pythonClasses.ExpressionNodeGreater);
  Py_DECREF(pythonClasses.ExpressionNodeGreaterEqual);
  Py_DECREF(pythonClasses.ExpressionNodeEqual);
  Py_DECREF(pythonClasses.ExpressionNodeUnequal);
  Py_DECREF(pythonClasses.ExpressionNodeNot);
  Py_DECREF(pythonClasses.ExpressionNodeAnd);
  Py_DECREF(pythonClasses.ExpressionNodeOr);
  Py_DECREF(pythonClasses.ExpressionNodeUniformRandom);
  Py_DECREF(pythonClasses.ExpressionNodeGaussianRandom);
  Py_DECREF(pythonClasses.PromoterElementConstitutive);
  Py_DECREF(pythonClasses.PromoterElementActivate);
  Py_DECREF(pythonClasses.PromoterElementRepress);
  Py_DECREF(pythonClasses.Factor);
  Py_DECREF(pythonClasses.Gene);
  Py_DECREF(pythonClasses.TranssysProgram);
  Py_DECREF(pythonClasses.GraphicsPrimitiveMove);
  Py_DECREF(pythonClasses.GraphicsPrimitivePush);
  Py_DECREF(pythonClasses.GraphicsPrimitivePop);
  Py_DECREF(pythonClasses.GraphicsPrimitiveTurn);
  Py_DECREF(pythonClasses.GraphicsPrimitiveRoll);
  Py_DECREF(pythonClasses.GraphicsPrimitiveBank);
  Py_DECREF(pythonClasses.GraphicsPrimitiveSphere);
  Py_DECREF(pythonClasses.GraphicsPrimitiveCylinder);
  Py_DECREF(pythonClasses.GraphicsPrimitiveBox);
  Py_DECREF(pythonClasses.GraphicsPrimitiveColor);
  Py_DECREF(pythonClasses.Symbol);
  Py_DECREF(pythonClasses.Assignment);
  Py_DECREF(pythonClasses.LhsSymbol);
  Py_DECREF(pythonClasses.ProductionElement);
  Py_DECREF(pythonClasses.Rule);
  Py_DECREF(pythonClasses.LsysProgram);
  Py_DECREF(pythonClasses.TranssysInstance);
  pythonClasses.ExpressionNodeValue = NULL;
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
  map[i].pythonClass = NULL;
  map[i++].node_type = NT_NONE;
}


/*
 * Initialise the pythonClass struct
 */
static int initPythonClasses(void)
{
  if (pythonClasses.ExpressionNodeValue != NULL)
  {
    return (0);
  }
  pythonClasses.ExpressionNodeValue = getPythonClass("transsys", "ExpressionNodeValue");
  if (pythonClasses.ExpressionNodeValue == NULL)
  {
    return (-1);
  }
  pythonClasses.ExpressionNodeIdentifier = getPythonClass("transsys", "ExpressionNodeIdentifier");
  if (pythonClasses.ExpressionNodeIdentifier == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeMult = getPythonClass("transsys", "ExpressionNodeMult");
  if (pythonClasses.ExpressionNodeMult == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeDiv = getPythonClass("transsys", "ExpressionNodeDiv");
  if (pythonClasses.ExpressionNodeDiv == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeAdd = getPythonClass("transsys", "ExpressionNodeAdd");
  if (pythonClasses.ExpressionNodeAdd == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeSubtract = getPythonClass("transsys", "ExpressionNodeSubtract");
  if (pythonClasses.ExpressionNodeSubtract == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeLower = getPythonClass("transsys", "ExpressionNodeLower");
  if (pythonClasses.ExpressionNodeLower == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeLowerEqual = getPythonClass("transsys", "ExpressionNodeLowerEqual");
  if (pythonClasses.ExpressionNodeLowerEqual == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeGreater = getPythonClass("transsys", "ExpressionNodeGreater");
  if (pythonClasses.ExpressionNodeGreater == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeGreaterEqual = getPythonClass("transsys", "ExpressionNodeGreaterEqual");
  if (pythonClasses.ExpressionNodeGreaterEqual == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeEqual = getPythonClass("transsys", "ExpressionNodeEqual");
  if (pythonClasses.ExpressionNodeEqual == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeUnequal = getPythonClass("transsys", "ExpressionNodeUnequal");
  if (pythonClasses.ExpressionNodeUnequal == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeNot = getPythonClass("transsys", "ExpressionNodeNot");
  if (pythonClasses.ExpressionNodeNot == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeAnd = getPythonClass("transsys", "ExpressionNodeAnd");
  if (pythonClasses.ExpressionNodeAnd == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeOr = getPythonClass("transsys", "ExpressionNodeOr");
  if (pythonClasses.ExpressionNodeOr == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeUniformRandom = getPythonClass("transsys", "ExpressionNodeUniformRandom");
  if (pythonClasses.ExpressionNodeUniformRandom == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ExpressionNodeGaussianRandom = getPythonClass("transsys", "ExpressionNodeGaussianRandom");
  if (pythonClasses.ExpressionNodeGaussianRandom == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.PromoterElementConstitutive = getPythonClass("transsys", "PromoterElementConstitutive");
  if (pythonClasses.PromoterElementConstitutive == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.PromoterElementActivate = getPythonClass("transsys", "PromoterElementActivate");
  if (pythonClasses.PromoterElementActivate == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.PromoterElementRepress = getPythonClass("transsys", "PromoterElementRepress");
  if (pythonClasses.PromoterElementRepress == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.Factor = getPythonClass("transsys", "Factor");
  if (pythonClasses.Factor == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.Gene = getPythonClass("transsys", "Gene");
  if (pythonClasses.Gene == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.TranssysProgram = getPythonClass("transsys", "TranssysProgram");
  if (pythonClasses.TranssysProgram == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveMove = getPythonClass("transsys", "GraphicsPrimitiveMove");
  if (pythonClasses.GraphicsPrimitiveMove == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitivePush = getPythonClass("transsys", "GraphicsPrimitivePush");
  if (pythonClasses.GraphicsPrimitivePush == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitivePop = getPythonClass("transsys", "GraphicsPrimitivePop");
  if (pythonClasses.GraphicsPrimitivePop == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveTurn = getPythonClass("transsys", "GraphicsPrimitiveTurn");
  if (pythonClasses.GraphicsPrimitiveTurn == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveRoll = getPythonClass("transsys", "GraphicsPrimitiveRoll");
  if (pythonClasses.GraphicsPrimitiveRoll == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveBank = getPythonClass("transsys", "GraphicsPrimitiveBank");
  if (pythonClasses.GraphicsPrimitiveBank == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveSphere = getPythonClass("transsys", "GraphicsPrimitiveSphere");
  if (pythonClasses.GraphicsPrimitiveSphere == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveCylinder = getPythonClass("transsys", "GraphicsPrimitiveCylinder");
  if (pythonClasses.GraphicsPrimitiveCylinder == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveBox = getPythonClass("transsys", "GraphicsPrimitiveBox");
  if (pythonClasses.GraphicsPrimitiveBox == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.GraphicsPrimitiveColor = getPythonClass("transsys", "GraphicsPrimitiveColor");
  if (pythonClasses.GraphicsPrimitiveColor == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.Symbol = getPythonClass("transsys", "Symbol");
  if (pythonClasses.Symbol == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.Assignment = getPythonClass("transsys", "Assignment");
  if (pythonClasses.Assignment == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.LhsSymbol = getPythonClass("transsys", "LhsSymbol");
  if (pythonClasses.LhsSymbol == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.ProductionElement = getPythonClass("transsys", "ProductionElement");
  if (pythonClasses.ProductionElement == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.Rule = getPythonClass("transsys", "Rule");
  if (pythonClasses.Rule == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.LsysProgram = getPythonClass("transsys", "LsysProgram");
  if (pythonClasses.LsysProgram == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  pythonClasses.TranssysInstance = getPythonClass("transsys", "TranssysInstance");
  if (pythonClasses.TranssysInstance == NULL)
  {
    pythonClasses.ExpressionNodeValue = NULL;
    return (-1);
  }
  Py_INCREF(pythonClasses.ExpressionNodeValue);
  Py_INCREF(pythonClasses.ExpressionNodeIdentifier);
  Py_INCREF(pythonClasses.ExpressionNodeMult);
  Py_INCREF(pythonClasses.ExpressionNodeDiv);
  Py_INCREF(pythonClasses.ExpressionNodeAdd);
  Py_INCREF(pythonClasses.ExpressionNodeSubtract);
  Py_INCREF(pythonClasses.ExpressionNodeLower);
  Py_INCREF(pythonClasses.ExpressionNodeLowerEqual);
  Py_INCREF(pythonClasses.ExpressionNodeGreater);
  Py_INCREF(pythonClasses.ExpressionNodeGreaterEqual);
  Py_INCREF(pythonClasses.ExpressionNodeEqual);
  Py_INCREF(pythonClasses.ExpressionNodeUnequal);
  Py_INCREF(pythonClasses.ExpressionNodeNot);
  Py_INCREF(pythonClasses.ExpressionNodeAnd);
  Py_INCREF(pythonClasses.ExpressionNodeOr);
  Py_INCREF(pythonClasses.ExpressionNodeUniformRandom);
  Py_INCREF(pythonClasses.ExpressionNodeGaussianRandom);
  Py_INCREF(pythonClasses.PromoterElementConstitutive);
  Py_INCREF(pythonClasses.PromoterElementActivate);
  Py_INCREF(pythonClasses.PromoterElementRepress);
  Py_INCREF(pythonClasses.Factor);
  Py_INCREF(pythonClasses.Gene);
  Py_INCREF(pythonClasses.TranssysProgram);
  Py_INCREF(pythonClasses.GraphicsPrimitiveMove);
  Py_INCREF(pythonClasses.GraphicsPrimitivePush);
  Py_INCREF(pythonClasses.GraphicsPrimitivePop);
  Py_INCREF(pythonClasses.GraphicsPrimitiveTurn);
  Py_INCREF(pythonClasses.GraphicsPrimitiveRoll);
  Py_INCREF(pythonClasses.GraphicsPrimitiveBank);
  Py_INCREF(pythonClasses.GraphicsPrimitiveSphere);
  Py_INCREF(pythonClasses.GraphicsPrimitiveCylinder);
  Py_INCREF(pythonClasses.GraphicsPrimitiveBox);
  Py_INCREF(pythonClasses.GraphicsPrimitiveColor);
  Py_INCREF(pythonClasses.Symbol);
  Py_INCREF(pythonClasses.Assignment);
  Py_INCREF(pythonClasses.LhsSymbol);
  Py_INCREF(pythonClasses.ProductionElement);
  Py_INCREF(pythonClasses.Rule);
  Py_INCREF(pythonClasses.LsysProgram);
  Py_INCREF(pythonClasses.TranssysInstance);
  initExpressionNodetypeMap(&pythonClasses, expressionNodetypeMap);
  return (0);
}


static PyObject *newTranssysInstance(PyObject *transsysProgram)
{
  PyObject *a;
  PyObject *kw;
  PyObject *transsysInstance = NULL;

  a = PyTuple_New(1);
  if (a == NULL)
  {
    return (NULL);
  }
  /* keep a reference to transsysProgram (PyTuple_SetItem steals) */
  Py_INCREF(transsysProgram);
  PyTuple_SetItem(a, 0, transsysProgram);
  kw = PyDict_New();
  if (kw == NULL)
  {
    Py_DECREF(a);
    return (NULL);
  }
  transsysInstance = PyInstance_New(pythonClasses.TranssysInstance, a, kw);
  Py_DECREF(a);
  Py_DECREF(kw);
  return (transsysInstance);
}


static const char *getStringAttribute(PyObject *o, char *attr)
{
  PyObject *p = PyObject_GetAttrString(o, attr);
  const char *s;

  if (p == NULL)
  {
    return (NULL);
  }
  s = PyString_AsString(p);
  Py_DECREF(p);
  return (s);
}


static EXPRESSION_NODE *extract_expression(PyObject *python_node);


static EXPR_NODE_TYPE identify_node_type(PyObject *python_node, EXPRESSION_NODETYPE_MAPENTRY map[])
{
  int i, return_value;

  for (i = 0; expressionNodetypeMap[i].pythonClass != NULL; i++)
  {
    return_value = PyObject_IsInstance(python_node, map[i].pythonClass);
    if (return_value == -1)
    {
      return (NT_NONE);
    }
    if (return_value)
    {
      return (map[i].node_type);
    }
  }
  return (NT_NONE);
}


static EXPRESSION_NODE *extract_expression_value(PyObject *python_node)
{
  PyObject *v_obj = PyObject_GetAttrString(python_node, "value");
  double v;

  clib_message(CLIB_ERROR_MESSAGE, "extract_expression_value: start\n");
  if (v_obj == NULL)
  {
    return (NULL);
  }
  v = PyFloat_AsDouble(v_obj);
  Py_DECREF(v_obj);
  return (new_expression_node(NT_VALUE, v));
}


static EXPRESSION_NODE *extract_expression_identifier(PyObject *python_node)
{
  int return_value;
  const char *s;
  PyObject *id_obj = PyObject_GetAttrString(python_node, "factor");
  EXPRESSION_NODE *node;

  clib_message(CLIB_ERROR_MESSAGE, "extract_expression_identifier: start\n");
  if (id_obj == NULL)
  {
    return (NULL);
  }
  return_value = PyObject_IsInstance(id_obj, pythonClasses.Factor);
  if (return_value == -1)
  {
    Py_DECREF(id_obj);
    return (NULL);
  }
  if (return_value)
  {
    PyObject *id_str = PyObject_GetAttrString(id_obj, "name");

    Py_DECREF(id_obj);
    if (id_str == NULL)
    {
      return (NULL);
    }
    id_obj = id_str;
  }
  if (PyString_Check(id_obj))
  {
    s = PyString_AsString(id_obj);
    node = new_expression_node(NT_RAW_IDENTIFIER, NULL, s);
    Py_DECREF(id_obj);
    return (node);
  }
  Py_DECREF(id_obj);
  PyErr_SetString(PyExc_TypeError, "bad factor type (neither Factor nor string)");
  return (NULL);
}


static EXPRESSION_NODE *extract_expression_binary(EXPR_NODE_TYPE node_type, PyObject *operand1, PyObject *operand2)
{
  EXPRESSION_NODE *arg1, *arg2;

  if ((operand1 == NULL) || (operand2 == NULL))
  {
    return (NULL);
  }
  arg1 = extract_expression(operand1);
  if (arg1 == NULL)
  {
    return (NULL);
  }
  arg2 = extract_expression(operand2);
  if (arg2 == NULL)
  {
    free_expression_tree(arg1);
    return (NULL);
  }
  return(new_expression_node(node_type, arg1, arg2));
}


static EXPRESSION_NODE *extract_expression_binary_operator(PyObject *python_node, EXPR_NODE_TYPE node_type)
{
  PyObject *operand1, *operand2;
  EXPRESSION_NODE *expression;

  operand1 = PyObject_GetAttrString(python_node, "operand1");
  if (operand1 == NULL)
  {
    return (NULL);
  }
  Py_INCREF(operand1);
  operand2 = PyObject_GetAttrString(python_node, "operand2");
  if (operand2 == NULL)
  {
    Py_DECREF(operand1);
    return (NULL);
  }
  Py_INCREF(operand2);
  expression = extract_expression_binary(node_type, operand1, operand2);
  Py_DECREF(operand1);
  Py_DECREF(operand2);
  return (expression);
}


static EXPRESSION_NODE *extract_expression_binary_function(PyObject *python_node, EXPR_NODE_TYPE node_type)
{
  PyObject *operand_list, *operand1, *operand2;
  EXPRESSION_NODE *expression;
  int list_size;

  operand_list = PyObject_GetAttrString(python_node, "operand");
  if (operand_list == NULL)
  {
    return (NULL);
  }
  list_size = PyList_Size(operand_list);
  if (list_size != 2)
  {
    Py_DECREF(operand_list);
    PyErr_SetString(PyExc_TypeError, "operand list does not have 2 elements");
    return (NULL);
  }
  operand1 = PyList_GetItem(operand_list, 0);
  if (operand1 == NULL)
  {
    Py_DECREF(operand_list);
    return (NULL);
  }
  Py_INCREF(operand1);
  operand2 = PyList_GetItem(operand_list, 1);
  if (operand2 == NULL)
  {
    Py_DECREF(operand1);
    Py_DECREF(operand_list);
    return (NULL);
  }
  Py_INCREF(operand2);
  Py_DECREF(operand_list);
  expression = extract_expression_binary(node_type, operand1, operand2);
  Py_DECREF(operand1);
  Py_DECREF(operand2);
  return (expression);
}


static EXPRESSION_NODE *extract_expression_not(PyObject *python_node)
{
  EXPRESSION_NODE *arg;
  PyObject *operand;
  

  operand = PyObject_GetAttrString(python_node, "operand");
  if (operand == NULL)
  {
    return (NULL);
  }
  arg = extract_expression(operand);
  Py_DECREF(operand);
  if (arg == NULL)
  {
    return (NULL);
  }
  return (new_expression_node(NT_NOT, arg));
}


static EXPRESSION_NODE *extract_expression(PyObject *python_node)
{
  EXPR_NODE_TYPE node_type;

  clib_message(CLIB_ERROR_MESSAGE, "extract_expression: starting\n");
  node_type = identify_node_type(python_node, expressionNodetypeMap);
  clib_message(CLIB_ERROR_MESSAGE, "extract_expression: node type identified as %d\n", (int) node_type);
  if (node_type == NT_NONE)
  {
    return (NULL);
  }
  switch (node_type)
  {
  case NT_VALUE :
    clib_message(CLIB_ERROR_MESSAGE, "extract_expression: extracting value\n");
    return (extract_expression_value(python_node));
  case NT_IDENTIFIER :
    clib_message(CLIB_ERROR_MESSAGE, "extract_expression: extracting identifier\n");
    return (extract_expression_identifier(python_node));
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
    clib_message(CLIB_ERROR_MESSAGE, "extract_expression: extracting binary (%d)\n", (int) node_type);
    return (extract_expression_binary_operator(python_node, node_type));
  case NT_RANDOM :
  case NT_GAUSS :
    clib_message(CLIB_ERROR_MESSAGE, "extract_expression: extracting binary (%d)\n", (int) node_type);
    return (extract_expression_binary_function(python_node, node_type));
  case NT_NOT :
    clib_message(CLIB_ERROR_MESSAGE, "extract_expression: extracting !\n");
    return (extract_expression_not(python_node));
  default :
    PyErr_SetString(PyExc_TypeError, "unrecognised type of expression node");
    /* set some exception indicating that no matching type was found */
    return (NULL);
  }
  clib_message(CLIB_ERROR_MESSAGE, "clib internal node extraction error\n");
  PyErr_SetString(PyExc_RuntimeError, "clib internal node extraction error");
  return (NULL);
}


static FACTOR_ELEMENT *extract_factor(PyObject *python_factor)
{
  PyObject *python_name, *python_decay, *python_diffusibility;
  EXPRESSION_NODE *decay_expression = NULL, *diffusibility_expression = NULL;
  FACTOR_ELEMENT *factor = NULL;
  const char *name;

  clib_message(CLIB_ERROR_MESSAGE, "extract_factor: getting decay\n");
  python_decay = PyObject_GetAttrString(python_factor, "decay_expression");
  clib_message(CLIB_ERROR_MESSAGE, "extract_factor: got decay_expression attribute\n");
  if (python_decay == NULL)
  {
    clib_message(CLIB_ERROR_ERROR, "failed to extract decay expression\n");
    return (NULL);
  }
  decay_expression = extract_expression(python_decay);
  clib_message(CLIB_ERROR_MESSAGE, "extract_factor: extracted decay expression: %p\n", (void *) decay_expression);
  Py_DECREF(python_decay);
  if (decay_expression == NULL)
  {
    clib_message(CLIB_ERROR_ERROR, "extract_factor: decay expression extraction failed\n");
    return (NULL);
  }
  clib_message(CLIB_ERROR_MESSAGE, "extract_factor: got decay\n");
  python_diffusibility = PyObject_GetAttrString(python_factor, "diffusibility_expression");
  if (python_diffusibility == NULL)
  {
    clib_message(CLIB_ERROR_ERROR, "failed to extract diffusibility expression\n");
    free_expression_tree(decay_expression);
    return (NULL);
  }
  diffusibility_expression = extract_expression(python_diffusibility);
  Py_DECREF(python_diffusibility);
  if (diffusibility_expression == NULL)
  {
    free_expression_tree(decay_expression);
    return (NULL);
  }
  clib_message(CLIB_ERROR_MESSAGE, "returning factor\n");
  python_name = PyObject_GetAttrString(python_factor, "name");
  if (python_name == NULL)
  {
    return (NULL);
  }
  name = PyString_AsString(python_name);
  factor = new_factor_element(name, decay_expression, diffusibility_expression);
  Py_DECREF(python_name);
  return (factor);
}


static int extract_factor_elements(PyObject *python_tp, TRANSSYS *tp)
{
  int factor_index, num_factors;
  PyObject *factor_list, *python_factor;
  FACTOR_ELEMENT *factor_element;

  clib_message(CLIB_ERROR_MESSAGE, "extracting factor elements\n");
  factor_list = PyObject_GetAttrString(python_tp, "factor_list");
  if (factor_list == NULL)
  {
    clib_message(CLIB_ERROR_MESSAGE, "no factor list found\n");
    return (-1);
  }
  num_factors = PyList_Size(factor_list);
  clib_message(CLIB_ERROR_MESSAGE, "extract_factor_elements: %d factors\n", num_factors);
  for (factor_index = 0; factor_index < num_factors; factor_index++)
  {
    clib_message(CLIB_ERROR_MESSAGE, "extracting factor #%d\n", factor_index);
    python_factor = PyList_GetItem(factor_list, factor_index);
    if (python_factor == NULL)
    {
      Py_DECREF(factor_list);
      return (-1);
    }
    Py_INCREF(python_factor);
    factor_element = extract_factor(python_factor);
    if (factor_element)
    {
      clib_message(CLIB_ERROR_MESSAGE, "got factor element \"%s\"\n", factor_element->name);
      add_factor_definition(tp, factor_element);
    }
    Py_DECREF(python_factor);
    if (factor_element == NULL)
    {
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
    return (PROMOTERELEMENT_NONE);
  }
  if (return_value)
  {
    return (PROMOTERELEMENT_CONSTITUTIVE);
  }
  return_value = PyObject_IsInstance(python_pe, pythonClasses.PromoterElementActivate);
  if (return_value == -1)
  {
    return (PROMOTERELEMENT_NONE);
  }
  if (return_value)
  {
    return (PROMOTERELEMENT_ACTIVATE);
  }
  return_value = PyObject_IsInstance(python_pe, pythonClasses.PromoterElementRepress);
  if (return_value == -1)
  {
    return (PROMOTERELEMENT_NONE);
  }
  if (return_value)
  {
    return (PROMOTERELEMENT_REPRESS);
  }
  PyErr_SetString(PyExc_TypeError, "unrecognised promoter element type");
  return (PROMOTERELEMENT_NONE);
}


static PROMOTER_ELEMENT *extract_promoterelement_constitutive(PyObject *python_pe)
{
  PyObject *python_expression;
  PROMOTER_ELEMENT *pe;
  EXPRESSION_NODE *expression;

  python_expression = PyObject_GetAttrString(python_pe, "expression");
  if (python_expression == NULL)
  {
    return (NULL);
  }
  expression = extract_expression(python_expression);
  Py_DECREF(python_expression);
  if (expression == NULL)
  {
    return (NULL);
  }
  pe = new_promoter_element(PROMOTERELEMENT_CONSTITUTIVE, 0, NULL, expression, NULL);
  return (pe);
}


static int extract_factor_index(PyObject *python_product, const TRANSSYS *tp)
{
  int return_value, factor_index;
  PyObject *python_name = NULL;
  const char *name = NULL;

  return_value = PyObject_IsInstance(python_product, pythonClasses.Factor);
  if (return_value == -1)
  {
    return (NO_INDEX);
  }
  if (return_value)
  {
    python_name = PyObject_GetAttrString(python_product, "name");
    if (python_name == NULL)
    {
      return (NO_INDEX);
    }
    name = PyString_AsString(python_name);
  }
  else
  {
    name = PyString_AsString(python_product);
  }
  if (name == NULL)
  {
    if (python_name != NULL)
    {
      Py_DECREF(python_name);
    }
    return (NO_INDEX);
  }
  factor_index = find_factor_index(tp, name);
  if (python_name != NULL)
  {
    Py_DECREF(python_name);
  }
  if (factor_index == NO_INDEX)
  {
    PyErr_SetString(PyExc_StandardError, "find_factor_index failed");
  }
  return (factor_index);
}


static INTEGER_ARRAY *extract_factor_index_array(PyObject *factor_list, const TRANSSYS *tp)
{
  INTEGER_ARRAY *factor_index_array = NULL;
  int i, num_factors;

  num_factors = PyList_Size(factor_list);
  for (i = 0; i < num_factors; i++)
  {
    PyObject *factor = PyList_GetItem(factor_list, i);
    int factor_index;

    if (factor == NULL)
    {
      free_integer_array(factor_index_array);
      return (NULL);
    }
    Py_INCREF(factor);
    factor_index = extract_factor_index(factor, tp);
    if (factor_index == NO_INDEX)
    {
      free_integer_array(factor_index_array);
      Py_DECREF(factor);
      return (NULL);
    }
    factor_index_array = extend_integer_array(factor_index_array, factor_index);
    if (factor_index_array == NULL)
    {
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
  PyObject *factor_list, *python_exp1, *python_exp2;
  PROMOTER_ELEMENT *pe = NULL;
  INTEGER_ARRAY *factor_index_array;
  EXPRESSION_NODE *expression1, *expression2;

  factor_list = PyObject_GetAttrString(python_pe, "factor_list");
  if (factor_list == NULL)
  {
    return (NULL);
  }
  factor_index_array = extract_factor_index_array(factor_list, tp);
  Py_DECREF(factor_list);
  if (factor_index_array == NULL)
  {
    return (NULL);
  }
  python_exp1 = PyObject_GetAttrString(python_pe, "expression1");
  if (python_exp1 == NULL)
  {
    free_integer_array(factor_index_array);
    return (NULL);
  }
  expression1 = extract_expression(python_exp1);
  Py_DECREF(python_exp1);
  python_exp2 = PyObject_GetAttrString(python_pe, "expression2");
  if (python_exp2 == NULL)
  {
    free_integer_array(factor_index_array);
    free_expression_tree(expression1);
    return (NULL);
  }
  expression2 = extract_expression(python_exp2);
  Py_DECREF(python_exp2);
  pe = create_promoter(pe_type, factor_index_array, expression1, expression2);
  if (pe == NULL)
  {
    free_expression_tree(expression1);
    free_expression_tree(expression2);
    PyErr_SetString(PyExc_StandardError, "create_promoter failed");
  }
  return (pe);
}


static PROMOTER_ELEMENT *extract_promoter_element(PyObject *python_pe, const TRANSSYS *tp)
{
  PROMOTERELEMENT_TYPE pe_type = identify_promoter_element(python_pe);

  switch (pe_type)
  {
  case PROMOTERELEMENT_CONSTITUTIVE :
    return (extract_promoterelement_constitutive(python_pe));
  case PROMOTERELEMENT_ACTIVATE :
  case PROMOTERELEMENT_REPRESS :
    return (extract_promoterelement_link(python_pe, pe_type, tp));
  default :
    /* identify_promoter_element has set the exception */
    return (NULL);
  }
}


static PROMOTER_ELEMENT *extract_promoter(PyObject *python_promoter, const TRANSSYS *tp)
{
  int num_promoter_elements, i;
  PROMOTER_ELEMENT *promoter = NULL, *pe;
  PyObject *python_pe;

  num_promoter_elements = PyList_Size(python_promoter);
  for (i = 0; i < num_promoter_elements; i++)
  {
    python_pe = PyList_GetItem(python_promoter, i);
    if (python_pe == NULL)
    {
      return (NULL);
    }
    Py_INCREF(python_pe);
    pe = extract_promoter_element(python_pe, tp);
    Py_DECREF(python_pe);
    promoter = extend_promoter_list(promoter, pe);
    clib_message(CLIB_ERROR_MESSAGE, "extract_promoter: list extended by element %d\n", i);
  }
  return (promoter);
}


static GENE_ELEMENT *extract_gene(PyObject *python_gene, const TRANSSYS *tp)
{
  GENE_ELEMENT *gene = NULL;
  PyObject *python_promoter, *python_product;
  PROMOTER_ELEMENT *promoter;
  int product_index;
  const char *name = getStringAttribute(python_gene, "name");

  clib_message(CLIB_ERROR_MESSAGE, "extract_gene: starting\n");
  if (name == NULL)
  {
    clib_message(CLIB_ERROR_ERROR, "extract_gene: no gene name found\n");
    return (NULL);
  }
  clib_message(CLIB_ERROR_MESSAGE, "gene: \"%s\"\n", name);
  python_promoter = PyObject_GetAttrString(python_gene, "promoter");
  if (python_promoter == NULL)
  {
    clib_message(CLIB_ERROR_ERROR, "extract_gene: no promoter found\n");
    return (NULL);
  }
  promoter = extract_promoter(python_promoter, tp);
  Py_DECREF(python_promoter);
  python_product = PyObject_GetAttrString(python_gene, "product");
  if (python_product == NULL)
  {
    free_promoter_list(promoter);
    return (NULL);
  }
  product_index = extract_factor_index(python_product, tp);
  Py_DECREF(python_product);
  if (product_index == NO_INDEX)
  {
    return (NULL);
  }
  gene = create_gene(name, promoter, product_index);
  if (gene == NULL)
  {
    PyErr_SetString(PyExc_StandardError, "create_gene failed");
    return (NULL);
  }
  return (gene);
}


static int extract_gene_elements(PyObject *python_tp, TRANSSYS *tp)
{
  int gene_index, num_genes;
  PyObject *gene_list, *python_gene;
  GENE_ELEMENT *gene_element;

  clib_message(CLIB_ERROR_MESSAGE, "extract_gene_elements: starting\n");
  gene_list = PyObject_GetAttrString(python_tp, "gene_list");
  if (gene_list == NULL)
  {
    clib_message(CLIB_ERROR_ERROR, "no gene list found\n");
    return (-1);
  }
  num_genes = PyList_Size(gene_list);
  clib_message(CLIB_ERROR_MESSAGE, "extract_gene_elements: %d genes\n", num_genes);
  for (gene_index = 0; gene_index < num_genes; gene_index++)
  {
    clib_message(CLIB_ERROR_MESSAGE, "extracting gene #%d\n", gene_index);
    python_gene = PyList_GetItem(gene_list, gene_index);
    if (python_gene == NULL)
    {
      Py_DECREF(gene_list);
      return (-1);
    }
    Py_INCREF(python_gene);
    gene_element = extract_gene(python_gene, tp);
    Py_DECREF(python_gene);
    if (gene_element == NULL)
    {
      Py_DECREF(gene_list);
      return (-1);
    }
    clib_message(CLIB_ERROR_MESSAGE, "got gene element \"%s\"\n", gene_element->name);
    add_gene_definition(tp, gene_element);
  }
  Py_DECREF(gene_list);
  clib_message(CLIB_ERROR_MESSAGE, "extract_gene_elements: returning, arrayed = %d\n", tp->arrayed);
  return (0);
}


static TRANSSYS *extract_transsys(PyObject *python_tp)
{
  PyObject *tp_name;
  TRANSSYS *tp;
  int isTranssysProgram = PyObject_IsInstance(python_tp, pythonClasses.TranssysProgram);

  if (isTranssysProgram != 1)
  {
    if (isTranssysProgram == 0)
    {
      PyErr_SetString(PyExc_TypeError, "not an instance of TranssysProgram");
    }
    return (NULL);
  }
  tp_name = PyObject_GetAttrString(python_tp, "name");
  if (tp_name == NULL)
  {
    return (NULL);
  }
  tp = new_transsys(PyString_AsString(tp_name));
  Py_DECREF(tp_name);
  if (extract_factor_elements(python_tp, tp) != 0)
  {
    free_transsys_list(tp);
    return (NULL);
  }
  if (extract_gene_elements(python_tp, tp) != 0)
  {
    clib_message(CLIB_ERROR_ERROR, "extract_transsys: extracting gene elements failed, freeing transsys\n");
    free_transsys_list(tp);
    clib_message(CLIB_ERROR_ERROR, "extract_transsys: transsys freed, returning NULL\n");
    return (NULL);
  }
  resolve_transsys(tp);
  return (tp);
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

  python_fc_start = PyObject_GetAttrString(python_ti, "factor_concentration");
  if (python_fc_start == NULL)
  {
    return (-1);
  }
  num_factors = PyList_Size(python_fc_start);
  if (num_factors != ti->transsys->num_factors)
  {
    clib_message(CLIB_ERROR_ERROR, "extract_initial_factor_concentrations: python_ti and ti do not match");
    Py_DECREF(python_fc_start);
    return (-1);
  }
  for (i = 0; i < num_factors; i++)
  {
    PyObject *python_c = PyList_GetItem(python_fc_start, i);

    if (python_c == NULL)
    {
      Py_DECREF(python_fc_start);
      return (-1);
    }
    Py_INCREF(python_c);
    ti->factor_concentration[i] = PyFloat_AsDouble(python_c);
    Py_DECREF(python_c);
  }
  Py_DECREF(python_fc_start);
  return (0);
}


int make_factor_concentration(TRANSSYS_INSTANCE *ti, PyObject *python_ti)
{
  long i;
  PyObject *python_fc = PyObject_GetAttrString(python_ti, "factor_concentration");

  if (python_fc == NULL)
  {
    return (-1);
  }
  if (ti->transsys->num_factors != PyList_Size(python_fc))
  {
    PyErr_SetString(PyExc_TypeError, "make_factor_concentration: ti and python_ti do not match");
  }
  for (i = 0; i < ti->transsys->num_factors; i++)
  {
    PyObject *python_c = PyFloat_FromDouble(ti->factor_concentration[i]);
    if (python_c == NULL)
    {
      Py_DECREF(python_fc);
      return (-1);
    }
    if (PyList_SetItem(python_fc, i, python_c) != 0)
    {
      Py_DECREF(python_fc);
      return (-1);
    }
  }
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
  int isTranssysInstance = PyObject_IsInstance(python_ti_start, pythonClasses.TranssysInstance);

  if (isTranssysInstance != 1)
  {
    if (isTranssysInstance == 0)
    {
      PyErr_SetString(PyExc_TypeError, "not an instance of TranssysInstance");
    }
    return (NULL);
  }
  python_tp = PyObject_GetAttrString(python_ti_start, "transsys_program");
  if (python_tp == NULL)
  {
    return (NULL);
  }
  tp = extract_transsys(python_tp);
  if (tp == NULL)
  {
    Py_DECREF(python_tp);
    return (NULL);
  }
  ti = new_transsys_instance(tp);
  if (ti == NULL)
  {
    free_transsys_list(tp);
    Py_DECREF(python_tp);
    return (NULL);
  }
  if (extract_initial_factor_concentrations(python_ti_start, ti) != 0)
  {
    free_transsys_instance(ti);
    free_transsys_list(tp);
    Py_DECREF(python_tp);
    return (NULL);
  }
  python_ti_list = PyList_New(0);
  if (python_ti_list == NULL)
  {
    free_transsys_instance(ti);
    free_transsys_list(tp);
    Py_DECREF(python_tp);
    return (NULL);
  }
  for (i = 0; i < num_timesteps; i++)
  {
    if (i % sampling_period == 0)
    {
      PyObject *python_ti = newTranssysInstance(python_tp);
      if (python_ti == NULL)
      {
	free_transsys_instance(ti);
	free_transsys_list(tp);
	Py_DECREF(python_tp);
	Py_DECREF(python_ti_list);
	return (NULL);
      }
      python_ts = PyInt_FromLong(i);
      if (PyObject_SetAttrString(python_ti, "timestep", python_ts) == -1)
      {
	free_transsys_instance(ti);
	free_transsys_list(tp);
	Py_DECREF(python_tp);
	Py_DECREF(python_ti_list);
	return (NULL);
      }
      if (make_factor_concentration(ti, python_ti) != 0)
      {
	free_transsys_instance(ti);
	free_transsys_list(tp);
	Py_DECREF(python_tp);
	Py_DECREF(python_ti_list);
	return (NULL);
      }
      if (PyList_Append(python_ti_list, python_ti) != 0)
      {
	free_transsys_instance(ti);
	free_transsys_list(tp);
	Py_DECREF(python_tp);
	Py_DECREF(python_ti_list);
	return (NULL);
      }
    }
    process_expression(ti);
  }
  Py_DECREF(python_tp);
  return (python_ti_list);
}


static PyObject *clib_timeseries(PyObject *self, PyObject *args)
{
  PyObject *python_ti_start, *python_int, *python_ti_list;
  long num_timesteps, sampling_period;


  /* clib_message(CLIB_ERROR_MESSAGE, "dummy: starting\n"); */
  if (initPythonClasses() != 0)
  {
    /* Error should be set by getPythonClass */
    return (NULL);
  }
  if (PyTuple_Size(args) != 3)
  {
    PyErr_SetString(PyExc_TypeError, "exactly 3 arguments required");
    return (NULL);
  }
  python_int = PyTuple_GetItem(args, 1);
  if (python_int == NULL)
  {
    return (NULL);
  }
  Py_INCREF(python_int);
  if (!PyInt_Check(python_int))
  {
    Py_DECREF(python_int);
    PyErr_SetString(PyExc_TypeError, "number of timesteps (arg 2) must be an int");
    return (NULL);
  }
  num_timesteps = PyInt_AsLong(python_int);
  Py_DECREF(python_int);
  python_int = PyTuple_GetItem(args, 2);
  if (python_int == NULL)
  {
    return (NULL);
  }
  Py_INCREF(python_int);
  if (!PyInt_Check(python_int))
  {
    Py_DECREF(python_int);
    PyErr_SetString(PyExc_TypeError, "sampling period (arg 3) must be an int");
    return (NULL);
  }
  sampling_period = PyInt_AsLong(python_int);
  Py_DECREF(python_int);
  python_ti_start = PyTuple_GetItem(args, 0);
  if (python_ti_start == NULL)
  {
    clib_message(CLIB_ERROR_ERROR, "PyTuple_GetItem failed\n");
    return (NULL);
  }
  Py_INCREF(python_ti_start);
  clib_message(CLIB_ERROR_MESSAGE, "got tp\n");
  clib_message(CLIB_ERROR_MESSAGE, "num_timesteps = %ld, sampling_period = %ld\n", num_timesteps, sampling_period);
  python_ti_list = transsysInstanceTimeSeries(python_ti_start, num_timesteps, sampling_period);
  Py_DECREF(python_ti_start);
  /* clib_message(CLIB_ERROR_MESSAGE, "dummy: returning\n"); */
  return (python_ti_list);
}


static PyObject *clib_dummy(PyObject *self, PyObject *args)
{
  Py_INCREF(Py_None);
  return (Py_None);
}


static PyObject *clib_setverbose(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "i", &error_importance_threshold))
  {
    return (NULL);
  }
  clib_message(CLIB_ERROR_MESSAGE, "transsys.clib: verbosity level set to %d\n", error_importance_threshold);
  Py_INCREF(Py_None);
  return (Py_None);
}


static PyMethodDef clib_methods[] = {
  {"timeseries", clib_timeseries, METH_VARARGS, "compute time series from a transsys instance"},
  {"dummy", clib_dummy, METH_VARARGS, "dummy test function for clib development"},
  {"setverbose", clib_setverbose, METH_VARARGS, "set verbosity level for transsys.clib module"},
  {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initclib(void)
{
  (void) Py_InitModule("transsys.clib", clib_methods);
}
