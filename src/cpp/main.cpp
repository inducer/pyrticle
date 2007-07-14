#include <vector>
#include <boost/python.hpp>
#include <hedge/base.hpp>




namespace {
  typedef boost::numeric::ublas::vector<double> vector;

  class point_cloud
  {
    void update_containing_elements()
    {

    }

    protected:
      std::vector<int> m_containing_elements;
      vector m_positions;
      vector m_velocities;

      struct element_info
      {
        hedge::affine_map       m_map;
      };

      std::vector<element_info> m_elements;
  };

}
BOOST_PYTHON_MODULE(_internal)
{
}
