#include <Python.h>
#include <numpy/arrayobject.h>
#include "efermi.c"

// Docstrings
static char module_docstring[] =
    "Interface to C99 implementation of a bisection method for the Fermi energy.";
static char efermi_docstring[] =
    "Find the Fermi energy using bisection.";

// Available functions
static PyMethodDef module_methods[] = {
    {"efermi", efermi_efermi, METH_VARARGS, efermi_docstring},
    {NULL, NULL, 0, NULL}}

// Initialize the module
PyMODINIT_FUNC init_efermi(void)
{
    PyObject m * = Py_InitModule3("_efermi", module_methods, module_docstring);
    if (m == NULL)
        return;

    // Load `numpy`
    import_array();
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
    int nkpt = (int)PyArray_DIM(bands_array, 0);
    int nbnd = (int)PyArray_DIM(bands_array, 1);

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