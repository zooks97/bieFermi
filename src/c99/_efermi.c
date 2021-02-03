#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include "efermi.c"

// Docstrings
static char module_docstring[] =
    "Interface to C99 implementation of a bisection method for the Fermi energy.";
static char efermi_docstring[] =
    "Find the Fermi energy using bisection.";

// Available functions
static PyMethodDef ModuleMethods[] = {
    {"efermi", efermi_efermi, METH_VARARGS, efermi_docstring},
    {NULL, NULL, 0, NULL} // Sentinel
};

// Module definition
static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "efermi", // name
    module_docstring, // doc
    -1, // size
    ModuleMethods, // methods
    NULL,
    NULL,
    NULL,
    NULL
};

// Module initialization
PyMODINIT_FUNC PyInit_efermi(void)
{
    // Load `numpy`
    import_array();

    return PyModule_Create(&module_def)
}

static PyObject *efermi_efermi(PyObject *self, PyObject *args)
{
    PyObject *bands_obj, *weights_obj;
    int nelec, stype;
    double swidth;

    // Parse the input tuple
    if (!PyArg_ParseTuple(args, "OOidi", &bands_obj, &weights_obj, &nelec, &swidth, &stype))
        return NULL;

    // Interpret the input objects as numpy arrays
    PyObject *bands_array = PyArray_FROM_OTF(bands_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *weights_array = PyArray_FROM_OTF(weights_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    // If that didn't work, throw an exception
    if (bands_array == NULL || weights_array == NULL)
    {
        Py_XDECREF(bands_array);
        Py_XDECREF(weights_array);
        return NULL;
    }

    // Get dimensions
    size_t nkpt = (size_t)PyArray_DIM(bands_array, 0);
    size_t nbnd = (size_t)PyArray_DIM(bands_array, 1);

    // Get pointers to the data as C-types
    double *bands = (double *)PyArray_DATA(bands_array);
    double *weights = (double *)PyArray_DATA(weights_array);

    double ef = efermi(nkpt, nbnd, bands, weights, nelec, swidth, stype);

    // Clean up
    Py_DECREF(bands_array);
    Py_DECREF(weights_array);

    // Build the output tuple
    PyObject *ret = Py_BuildValue("d", ef);
    return ret;
}