/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//<LicenseText>
//
// CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
// Copyright (C) 2002-2005, California Institute of Technology.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//</LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

#include <Python.h>

#include "exceptions.h"
#include "mesher.h"

#include "global_defs.h"
#include "parallel_related.h"


extern void initial_mesh_solver_setup(struct All_variables *);



char pyCitcom_set_3dsphere_defaults__doc__[] = "";
char pyCitcom_set_3dsphere_defaults__name__[] = "set_3dsphere_defaults";

PyObject * pyCitcom_set_3dsphere_defaults(PyObject *self, PyObject *args)
{
    PyObject *obj;
    struct All_variables* E;

    if (!PyArg_ParseTuple(args, "O:set_3dsphere_defaults", &obj))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));

    (E->solver.set_3dsphere_defaults)(E);

    Py_INCREF(Py_None);
    return Py_None;
}



char pyCitcom_full_sphere_launch__doc__[] = "";
char pyCitcom_full_sphere_launch__name__[] = "full_sphere_launch";

PyObject * pyCitcom_full_sphere_launch(PyObject *self, PyObject *args)
{
    PyObject *obj;
    struct All_variables* E;

    if (!PyArg_ParseTuple(args, "O:full_sphere_launch", &obj))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));

    initial_mesh_solver_setup(E);

    Py_INCREF(Py_None);
    return Py_None;
}



char pyCitcom_regional_sphere_launch__doc__[] = "";
char pyCitcom_regional_sphere_launch__name__[] = "regional_sphere_launch";

PyObject * pyCitcom_regional_sphere_launch(PyObject *self, PyObject *args)
{
    PyObject *obj;
    struct All_variables* E;

    if (!PyArg_ParseTuple(args, "O:regional_sphere_launch", &obj))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));

    initial_mesh_solver_setup(E);

    Py_INCREF(Py_None);
    return Py_None;
}



/* $Id: mesher.c 14352 2009-03-16 22:50:06Z tan2 $ */

/* End of file */
