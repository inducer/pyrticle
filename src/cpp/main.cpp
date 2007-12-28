#include <stdexcept>
#include <vector>
#include <list>
#include <climits>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <hedge/base.hpp>
#include <boost/foreach.hpp> 




using namespace boost::assign;
namespace python = boost::python;
namespace ublas = boost::numeric::ublas;




namespace {
  template<unsigned dimensions_pos, unsigned dimensions_velocity>
  class particle_cloud : boost::noncopyable
  {
    public:
      // member data ----------------------------------------------------------
      python::dict                      m_vis_info;

      const double                      m_vacuum_c;




      // setup ----------------------------------------------------------------
      particle_cloud(
          unsigned dimensions_mesh,
          unsigned vertices_sizehint,
          unsigned elements_sizehint, 
          unsigned discretizations_sizehint,
          double vacuum_c)
        : m_mesh_info(dimensions_mesh, 
            vertices_sizehint,
            elements_sizehint, 
            discretizations_sizehint),
        m_vacuum_c(vacuum_c)
      {
      }





      // operation ------------------------------------------------------------








  };
}




BOOST_PYTHON_MODULE(_internal)
{
  using python::arg;

}
