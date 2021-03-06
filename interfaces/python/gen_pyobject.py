#! /usr/bin/python

""" This madness of code is used to generate the C code of the python interface
to aubio. Don't try this at home.

The list of typedefs and functions is obtained from the command line 'cpp
aubio.h'. This list is then used to parse all the functions about this object.

I hear the ones asking "why not use swig, or cython, or something like that?"

The requirements for this extension are the following:

    - aubio vectors can be viewed as numpy arrays, and vice versa
    - aubio 'object' should be python classes, not just a bunch of functions

I haven't met any python interface generator that can meet both these
requirements. If you know of one, please let me know, it will spare me
maintaining this bizarre file.
"""

# TODO
# do function: for now, only the following pattern is supported:
# void aubio_<foo>_do (aubio_foo_t * o, 
#       [input1_t * input, [output1_t * output, ..., output3_t * output]]);
# There is no way of knowing that output1 is actually input2. In the future,
# const could be used for the inputs in the C prototypes.

def split_type(arg):
    """ arg = 'foo *name' 
        return ['foo*', 'name'] """
    l = arg.split()
    # ['foo', '*name'] -> ['foo*', 'name']
    if l[-1].startswith('*'):
        return [l[0]+'*', l[1][1:]]
    # ['foo', '*', 'name'] -> ['foo*', 'name']
    if len(l) == 3:
        return [l[0]+l[1], l[2]]
    else:
        return l

def get_params(proto):
    """ get the list of parameters from a function prototype
    example: proto = "int main (int argc, char ** argv)"
    returns: ['int argc', 'char ** argv']
    """
    import re
    paramregex = re.compile('[\(, ](\w+ \*?\*? ?\w+)[, \)]')
    return paramregex.findall(proto)

def get_params_types_names(proto):
    """ get the list of parameters from a function prototype
    example: proto = "int main (int argc, char ** argv)"
    returns: [['int', 'argc'], ['char **','argv']]
    """
    return map(split_type, get_params(proto)) 

def get_return_type(proto):
    import re
    paramregex = re.compile('(\w+ ?\*?).*')
    outputs = paramregex.findall(proto)
    assert len(outputs) == 1
    return outputs[0].replace(' ', '')

def get_name(proto):
    name = proto.split()[1].split('(')[0]
    return name.replace('*','')

# the important bits: the size of the output for each objects. this data should
# move into the C library at some point.
defaultsizes = {
    'resampler':    'input->length * self->ratio',
    'specdesc':     '1',
    'onset':        '1',
    'pitchyin':     '1',
    'pitchyinfft':  '1',
    'pitchschmitt': '1',
    'pitchmcomb':   '1',
    'pitchfcomb':   '1',
    'pitch':        '1',
    'tss':          'self->hop_size',
    'mfcc':         'self->n_coeffs',
    'beattracking': 'self->hop_size',
    'tempo':        '1',
    'peakpicker':   '1',
}

# default value for variables
aubioinitvalue = {
    'uint_t': 0,
    'smpl_t': 0,
    'lsmp_t': 0.,
    'char_t*': 'NULL',
    }

aubiodefvalue = {
    # we have some clean up to do
    'buf_size': 'Py_default_vector_length', 
    # and here too
    'hop_size': 'Py_default_vector_length / 2', 
    # these should be alright
    'samplerate': 'Py_aubio_default_samplerate', 
    # now for the non obvious ones
    'n_filters': '40', 
    'n_coeffs': '13', 
    'nelems': '10',
    'flow': '0.', 
    'fhig': '1.', 
    'ilow': '0.', 
    'ihig': '1.', 
    'thrs': '0.5',
    'ratio': '0.5',
    'method': '"default"',
    }

# aubio to python
aubio2pytypes = {
    'uint_t': 'I',
    'smpl_t': 'f',
    'lsmp_t': 'd',
    'fvec_t*': 'O',
    'cvec_t*': 'O',
    'char_t*': 's',
}

# python to aubio
aubiovecfrompyobj = {
    'fvec_t*': 'PyAubio_ArrayToCFvec',
    'cvec_t*': 'PyAubio_ArrayToCCvec',
}

# aubio to python
aubiovectopyobj = {
    'fvec_t*': 'PyAubio_CFvecToArray',
    'cvec_t*': 'PyAubio_CCvecToPyCvec',
    'smpl_t': 'PyFloat_FromDouble',
}

def gen_new_init(newfunc, name):
    newparams = get_params_types_names(newfunc)
    # self->param1, self->param2, self->param3
    if len(newparams):
        selfparams = ', self->'+', self->'.join([p[1] for p in newparams])
    else:
        selfparams = '' 
    # "param1", "param2", "param3"
    paramnames = ", ".join(["\""+p[1]+"\"" for p in newparams])
    pyparams = "".join(map(lambda p: aubio2pytypes[p[0]], newparams))
    paramrefs = ", ".join(["&" + p[1] for p in newparams])
    s = """\
// WARNING: this file is generated, DO NOT EDIT

// WARNING: if you haven't read the first line yet, please do so
#include "aubiowraphell.h"

typedef struct
{
  PyObject_HEAD
  aubio_%(name)s_t * o;
""" % locals()
    for ptype, pname in newparams:
        s += """\
  %(ptype)s %(pname)s;
""" % locals()
    s += """\
} Py_%(name)s;

static char Py_%(name)s_doc[] = "%(name)s object";

static PyObject *
Py_%(name)s_new (PyTypeObject * pytype, PyObject * args, PyObject * kwds)
{
  Py_%(name)s *self;
""" % locals()
    for ptype, pname in newparams:
        initval = aubioinitvalue[ptype]
        s += """\
  %(ptype)s %(pname)s = %(initval)s;
""" % locals()
    # now the actual PyArg_Parse
    if len(paramnames):
        s += """\
  static char *kwlist[] = { %(paramnames)s, NULL };

  if (!PyArg_ParseTupleAndKeywords (args, kwds, "|%(pyparams)s", kwlist,
          %(paramrefs)s)) {
    return NULL;
  }
""" % locals()
    s += """\

  self = (Py_%(name)s *) pytype->tp_alloc (pytype, 0);

  if (self == NULL) {
    return NULL;
  }
""" % locals()
    for ptype, pname in newparams:
        defval = aubiodefvalue[pname]
        if ptype == 'char_t*':
            s += """\

  self->%(pname)s = %(defval)s;
  if (%(pname)s != NULL) {
    self->%(pname)s = %(pname)s;
  }
""" % locals()
        elif ptype == 'uint_t':
            s += """\

  self->%(pname)s = %(defval)s;
  if (%(pname)s > 0) {
    self->%(pname)s = %(pname)s;
  } else if (%(pname)s < 0) {
    PyErr_SetString (PyExc_ValueError,
        "can not use negative value for %(pname)s");
    return NULL;
  }
""" % locals()
        elif ptype == 'smpl_t':
            s += """\

  self->%(pname)s = %(defval)s;
  if (%(pname)s != %(defval)s) {
    self->%(pname)s = %(pname)s;
  }
""" % locals()
        else:
            print "ERROR, unknown type of parameter %s %s" % (ptype, pname)
    s += """\

  return (PyObject *) self;
}

AUBIO_INIT(%(name)s %(selfparams)s)

AUBIO_DEL(%(name)s)

""" % locals()
    return s

def gen_do(dofunc, name):
    funcname = dofunc.split()[1].split('(')[0]
    doparams = get_params_types_names(dofunc) 
    # make sure the first parameter is the object
    assert doparams[0][0] == "aubio_"+name+"_t*", \
        "method is not in 'aubio_<name>_t"
    # and remove it
    doparams = doparams[1:]
    # guess the input/output params, assuming we have less than 3
    assert len(doparams) > 0, \
        "no parameters for function do in object %s" % name
    #assert (len(doparams) <= 2), \
    #    "more than 3 parameters for do in object %s" % name

    # build strings for inputs, assuming there is only one input 
    inputparams = [doparams[0]]
    # build the parsing string for PyArg_ParseTuple
    pytypes = "".join([aubio2pytypes[p[0]] for p in doparams[0:1]])
    inputdefs = "\n  ".join(["PyObject * " + p[-1] + "_obj;" for p in inputparams])
    inputvecs = "\n  ".join(map(lambda p: \
                p[0] + p[-1] + ";", inputparams))
    parseinput = ""
    for p in inputparams:
        inputvec = p[-1]
        inputdef = p[-1] + "_obj"
        converter = aubiovecfrompyobj[p[0]]
        parseinput += """%(inputvec)s = %(converter)s (%(inputdef)s);

  if (%(inputvec)s == NULL) {
    return NULL;
  }""" % locals()
    # build the string for the input objects references
    inputrefs = ", ".join(["&" + p[-1] + "_obj" for p in inputparams])
    # end of inputs strings

    # build strings for outputs
    outputparams = doparams[1:]
    if len(outputparams) >= 1:
        #assert len(outputparams) == 1, \
        #    "too many output parameters"
        outputvecs = "\n  ".join([p[0] + p[-1] + ";" for p in outputparams])
        outputcreate = "\n  ".join(["""\
  %(name)s = new_%(autype)s (%(length)s);""" % \
    {'name': p[-1], 'pytype': p[0], 'autype': p[0][:-3],
        'length': defaultsizes[name]} \
        for p in outputparams]) 
        if len(outputparams) > 1:
            returnval = "PyObject *outputs = PyList_New(0);\n"
            for p in outputparams:
                returnval += "  PyList_Append( outputs, (PyObject *)" + aubiovectopyobj[p[0]] + " (" + p[-1] + ")" +");\n"
            returnval += "  return outputs;"
        else:
            returnval = "return (PyObject *)" + aubiovectopyobj[p[0]] + " (" + p[-1] + ")"
    else:
        # no output
        outputvecs = ""
        outputcreate = ""
        #returnval = "Py_None";
        returnval = "return (PyObject *)" + aubiovectopyobj[p[0]] + " (" + p[-1] + ")"
    # end of output strings

    # build the parameters for the  _do() call
    doparams_string = "self->o, " + ", ".join([p[-1] for p in doparams])

    # put it all together
    s = """\
static PyObject * 
Py_%(name)s_do(Py_%(name)s * self, PyObject * args)
{
  %(inputdefs)s
  %(inputvecs)s
  %(outputvecs)s

  if (!PyArg_ParseTuple (args, "%(pytypes)s", %(inputrefs)s)) {
    return NULL;
  }

  %(parseinput)s
  
  %(outputcreate)s

  /* compute _do function */
  %(funcname)s (%(doparams_string)s);

  %(returnval)s;
}
""" % locals()
    return s

def gen_members(new_method, name):
    newparams = get_params_types_names(new_method)
    s = """
AUBIO_MEMBERS_START(%(name)s)""" % locals()
    for param in newparams:
        if param[0] == 'char_t*':
            s += """
  {"%(pname)s", T_STRING, offsetof (Py_%(name)s, %(pname)s), READONLY, ""},""" \
        % { 'pname': param[1], 'ptype': param[0], 'name': name}
        elif param[0] == 'uint_t':
            s += """
  {"%(pname)s", T_INT, offsetof (Py_%(name)s, %(pname)s), READONLY, ""},""" \
        % { 'pname': param[1], 'ptype': param[0], 'name': name}
        elif param[0] == 'smpl_t':
            s += """
  {"%(pname)s", T_FLOAT, offsetof (Py_%(name)s, %(pname)s), READONLY, ""},""" \
        % { 'pname': param[1], 'ptype': param[0], 'name': name}
        else:
            print "-- ERROR, unknown member type ", param
    s += """
AUBIO_MEMBERS_STOP(%(name)s)

""" % locals()
    return s


def gen_methods(get_methods, set_methods, name):
    s = ""
    method_defs = ""
    for method in set_methods:
        method_name = get_name(method)
        params = get_params_types_names(method)
        out_type = get_return_type(method)
        assert params[0][0] == "aubio_"+name+"_t*", \
            "get method is not in 'aubio_<name>_t"
        print method
        print params[1:]
        setter_args = "self->o, " +",".join([p[1] for p in params[1:]])
        parse_args = ""
        for p in params[1:]:
            parse_args += p[0] + " " + p[1] + ";\n"
        argmap = "".join([aubio2pytypes[p[0]] for p in params[1:]])
        arglist = ", ".join(["&"+p[1] for p in params[1:]])
        parse_args += """
  if (!PyArg_ParseTuple (args, "%(argmap)s", %(arglist)s)) {
    return NULL;
  } """ % locals()
        s += """
static PyObject *
Py%(funcname)s (Py_%(objname)s *self, PyObject *args)
{
  uint_t err = 0;

  %(parse_args)s

  err = %(funcname)s (%(setter_args)s);

  if (err > 0) {
    PyErr_SetString (PyExc_ValueError,
        "error running %(funcname)s");
    return NULL;
  }
  return Py_None;
}
""" % {'funcname': method_name, 'objname': name, 
        'out_type': out_type, 'setter_args': setter_args, 'parse_args': parse_args }
        shortname = method_name.split(name+'_')[-1]
        method_defs += """\
  {"%(shortname)s", (PyCFunction) Py%(method_name)s,
    METH_VARARGS, ""},
""" % locals()

    for method in get_methods:
        method_name = get_name(method)
        params = get_params_types_names(method)
        out_type = get_return_type(method)
        assert params[0][0] == "aubio_"+name+"_t*", \
            "get method is not in 'aubio_<name>_t %s" % params[0][0]
        assert len(params) == 1, \
            "get method has more than one parameter %s" % params
        getter_args = "self->o" 
        returnval = "(PyObject *)" + aubiovectopyobj[out_type] + " (tmp)"
        shortname = method_name.split(name+'_')[-1]
        method_defs += """\
  {"%(shortname)s", (PyCFunction) Py%(method_name)s,
    METH_NOARGS, ""},
""" % locals()
        s += """
static PyObject *
Py%(funcname)s (Py_%(objname)s *self, PyObject *unused)
{
  %(out_type)s tmp = %(funcname)s (%(getter_args)s);
  return %(returnval)s;
}
""" % {'funcname': method_name, 'objname': name, 
        'out_type': out_type, 'getter_args': getter_args, 'returnval': returnval }

    s += """
static PyMethodDef Py_%(name)s_methods[] = {
""" % locals() 
    s += method_defs 
    s += """\
  {NULL} /* sentinel */
};
""" % locals() 
    return s

def gen_finish(name):
    s = """\

AUBIO_TYPEOBJECT(%(name)s, "aubio.%(name)s")
""" % locals()
    return s
