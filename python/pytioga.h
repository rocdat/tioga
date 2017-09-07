#ifndef PYTIOGA_H
#define PYTIOGA_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <boost/python.hpp>
#include "boost/python/extract.hpp"
#include <numpy/ndarrayobject.h>
#include "tioga.h"

class PyTioga {

  tioga *tg;

  Py_buffer b_coord;
  Py_buffer b_wallnode;
  Py_buffer b_obcnode;
  Py_buffer b_iblank;
  Py_buffer b_tetConn;
  Py_buffer b_pyraConn;
  Py_buffer b_prismConn;
  Py_buffer b_hexaConn;
  Py_buffer b_qvars;

  int initialized;

  int mpi_size, mpi_rank;
  int nv2[4],ncell[4];
  int *connptr[4];
  int *ndc4, *ndc5, *ndc6, *ndc8;


  MPI_Comm comm;

 public:
  PyTioga();
  PyTioga(boost::python::object pycomm);
  ~PyTioga();

  void setup(MPI_Comm comm);
  void register_data(boost::python::dict data);
  void preprocess_grids();
  void perform_connectivity();
  void register_solution(int block_id, boost::python::object qo);
  void data_update(int nq);
  void write_data(int nq);

};

#endif
