#include "Python.h"
static PyMethodDef SpamMethods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};
PyMODINIT_FUNC
initext_falcon(void)
{
    PyObject *m;

    m = Py_InitModule("ext_falcon", SpamMethods);
    if (m == NULL)
        return;
}
