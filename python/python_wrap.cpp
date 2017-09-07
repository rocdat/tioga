
#include "pytioga.h"

BOOST_PYTHON_MODULE(libpytioga){

  using namespace boost::python;

  import_array();

  class_<PyTioga>("PyTioga", init<object>())
    .def(init<>())
    .def("register_data", &PyTioga::register_data)
    .def("preprocess_grids", &PyTioga::preprocess_grids)
    .def("perform_connectivity", &PyTioga::perform_connectivity)
    .def("register_solution", &PyTioga::register_solution)
    .def("write_data", &PyTioga::write_data)
    .def("data_update", &PyTioga::data_update);
  

}
