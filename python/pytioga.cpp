// ##################################################################
//
// Python interface for TIOGA through BOOST
//
// ##################################################################
#include "pytioga.h"
#include "mpi.h"
#include "tioga.h"
#include <mpi4py/mpi4py.h>

// ==================================================================
// constructor
// ==================================================================

PyTioga::PyTioga(){
  this->setup(MPI_COMM_WORLD);
}

PyTioga::PyTioga(boost::python::object py_comm){
  PyObject* py_obj = py_comm.ptr();
  if(import_mpi4py()<0) return;
  MPI_Comm *comm = PyMPIComm_Get(py_obj);
  if(comm == NULL){
    printf("Error with the python mpi communicator\n");
    throw 28;
  }
  this->initialized = 0;
  this->setup(*comm);
}

void PyTioga::setup(MPI_Comm comm){
  int dummy = 0;
  int size, rank;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  this->comm     = comm;
  this->mpi_size = size;
  this->mpi_rank = rank;

  tg = new tioga[1];

  tg->setCommunicator(comm,rank,size);

}

// ==================================================================
// destructor
// ==================================================================

PyTioga::~PyTioga(){

  //
  // Release hold of the python data buffers
  //
  if(initialized){
    PyBuffer_Release(&b_hexaConn);
    PyBuffer_Release(&b_wallnode);
    PyBuffer_Release(&b_obcnode);
    PyBuffer_Release(&b_coord);
    PyBuffer_Release(&b_iblank);
  }
  //
  delete [] tg;
  
}

// ==================================================================
// check dictionary
// ==================================================================

using namespace boost::python;

bool check_python_data(dict data){
  if(not data.has_key("tetConn")){           printf("missing tetconn key\n");          return 1; }
  if(not data.has_key("pyraConn")){          printf("missing pyraconn key\n");         return 1; }
  if(not data.has_key("prismConn")){         printf("missing prismconn key\n");        return 1; }
  if(not data.has_key("hexaConn")){          printf("missing hexaconn key\n");         return 1; }
  if(not data.has_key("wallnode")){          printf("missing wallnode key\n");         return 1; }
  if(not data.has_key("obcnode")){           printf("missing obcnode key\n");          return 1; }
  if(not data.has_key("bodyTag")){           printf("missing btag key\n");             return 1; }
  if(not data.has_key("grid-coordinates")){  printf("missing grid-coordinates key\n"); return 1; }
  if(not data.has_key("iblanking")){         printf("missing iblanking key\n");        return 1; }
  // do more checks here, like perhaps sizing checks?
  return 0;
}

void PyTioga::register_data(dict data){

  if( check_python_data(data) ){
    printf("Error in Dictionary Inputs\n");
    MPI_Abort(this->comm,789);
  }

  int btag, bid, iblk, nv, nwbc, nobc;
  double* xyz;
  int *wbcnode, *obcnode, *iblank;
  int *tmpbtag;
  int ntypes, tmpcount;
  int dbg, nq;

  // list tmplist; // things are stored in the dictionary as lists.

  // list test = extract<list>(data["hexaConn"]);
  // list test = (list)data["hexaConn"];
  // test[0];

  // cast as list, get first item, and cast as object, then get pointer to data
  PyObject_GetBuffer(((object)((list)data["tetConn"])[0]).ptr(),          &b_tetConn,     PyBUF_SIMPLE);
  PyObject_GetBuffer(((object)((list)data["pyraConn"])[0]).ptr(),         &b_pyraConn,    PyBUF_SIMPLE);
  PyObject_GetBuffer(((object)((list)data["prismConn"])[0]).ptr(),        &b_prismConn,   PyBUF_SIMPLE);
  PyObject_GetBuffer(((object)((list)data["hexaConn"])[0]).ptr(),         &b_hexaConn,    PyBUF_SIMPLE);
  PyObject_GetBuffer(((object)((list)data["wallnode"])[0]).ptr(),         &b_wallnode,    PyBUF_SIMPLE);
  PyObject_GetBuffer(((object)((list)data["obcnode"])[0]).ptr(),          &b_obcnode,     PyBUF_SIMPLE);
  PyObject_GetBuffer(((object)((list)data["grid-coordinates"])[0]).ptr(), &b_coord,       PyBUF_SIMPLE);
  PyObject_GetBuffer(((object)((list)data["iblanking"])[0]).ptr(),        &b_iblank,      PyBUF_SIMPLE);
  initialized = 1;
  xyz         = reinterpret_cast<double*>(b_coord.buf);
  wbcnode     = reinterpret_cast<int*   >(b_wallnode.buf);
  obcnode     = reinterpret_cast<int*   >(b_obcnode.buf);
  iblank      = reinterpret_cast<int*   >(b_iblank.buf);
  ndc4        = reinterpret_cast<int*   >(b_tetConn.buf);
  ndc5        = reinterpret_cast<int*   >(b_pyraConn.buf);
  ndc6        = reinterpret_cast<int*   >(b_prismConn.buf);
  ndc8        = reinterpret_cast<int*   >(b_hexaConn.buf);

  // temporary buffer usage for bodyTag, which is an array of ints
  PyObject_GetBuffer(((object)((list)data["bodyTag"])[0]).ptr(),          &b_qvars,      PyBUF_SIMPLE);
  tmpbtag = reinterpret_cast<int*   >(b_qvars.buf);
  btag    = tmpbtag[0];
  PyBuffer_Release(&b_qvars);    

  bid        = 0;
  nwbc       = extract<int>(((list)data["wallnode"])[0].attr("size"));
  nobc       = extract<int>(((list)data["obcnode"])[0].attr("size"));
  nv         = extract<int>(((list)data["grid-coordinates"])[0].attr("size"))/3;

  // nq         = extract<int>(((list)data["q-variables"])[0].attr("size"));
  // printf("nq, nv, nq/nv: %d %d %f\n", nq, nv, 1.0*nq/nv);

  // Now we'll look through all the connectivity data and see which
  // one are / aren't used. Keep track of the number of types of
  // connectivity.
  ntypes = 0;

  // 4-vertex Cells
  tmpcount   = extract<int>(((list)data["tetConn"])[0].attr("size"))/4;
  if(tmpcount > 0){
    ncell[ntypes]   = tmpcount;
    nv2[ntypes]     = 4;
    connptr[ntypes] = ndc4;
    ntypes++;
  }
  // 5-vertex Cells
  tmpcount   = extract<int>(((list)data["pyraConn"])[0].attr("size"))/5;
  if(tmpcount > 0){
    ncell[ntypes]   = tmpcount;
    nv2[ntypes]     = 5;
    connptr[ntypes] = ndc5;
    ntypes++;
  }
  // 6-vertex Cells
  tmpcount   = extract<int>(((list)data["prismConn"])[0].attr("size"))/6;
  if(tmpcount > 0){
    ncell[ntypes]   = tmpcount;
    nv2[ntypes]     = 6;
    connptr[ntypes] = ndc6;
    ntypes++;
  }
  // 8-vertex Cells
  tmpcount   = extract<int>(((list)data["hexaConn"])[0].attr("size"))/8;
  if(tmpcount > 0){
    ncell[ntypes]   = tmpcount;
    nv2[ntypes]     = 8;
    connptr[ntypes] = ndc8;
    ntypes++;
  }

  if(ntypes == 0){
    printf("Error: using tioga without providing data\n");
  } else {
    tg->registerGridData(btag,nv,xyz,iblank,nwbc,nobc,wbcnode,obcnode,ntypes,nv2,ncell,connptr);
  } 

}

void PyTioga::preprocess_grids(){
  tg->profile();
}

void PyTioga::perform_connectivity(){
  tg->performConnectivity();
}

void PyTioga::register_solution(int block_id, object qo){

  double* q;

  PyObject_GetBuffer(qo.ptr(), &b_qvars, PyBUF_SIMPLE);

  q = reinterpret_cast<double*>(b_qvars.buf);

  tg->registerSolution(block_id,q);

  PyBuffer_Release(&b_qvars);
}

void PyTioga::data_update(int nq){
  tg->dataUpdate(nq,0);   // only "row" update implemented
}

void PyTioga::write_data(int nvar){
  tg->writeData(nvar, 0); // only "row" update implemented
}
