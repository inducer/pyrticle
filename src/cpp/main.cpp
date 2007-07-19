#include <stdexcept>
#include <vector>
#include <list>
#include <climits>
#include <numeric>
#include <algorithm>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <hedge/base.hpp>
#include <boost/assign/list_of.hpp> 
#include <boost/foreach.hpp> 
#include <boost/foreach.hpp> 
#include <boost/ptr_container/ptr_vector.hpp>




using namespace boost::assign;
namespace python = boost::python;
namespace ublas = boost::numeric::ublas;




#define COPY_PY_LIST(CPP_LIST, PY_LIST, EXTRACT_TYPE) \
        for (unsigned i = 0; i < unsigned(len(PY_LIST)); i++) \
          (CPP_LIST).push_back( \
              python::extract<EXTRACT_TYPE>((PY_LIST)[i]));

#define DEF_SIMPLE_METHOD(NAME) \
  def(#NAME, &cl::NAME)
#define DEF_SIMPLE_FUNCTION(NAME) \
  def(#NAME, &NAME)




namespace {
  typedef ublas::compressed_matrix<
    double, ublas::column_major, 0, 
    ublas::unbounded_array<int> >
      csr_matrix;




  struct monomial_basis_function
  {
    std::vector<unsigned> m_exponents;

    monomial_basis_function(python::list exponents)
    { 
      COPY_PY_LIST(m_exponents, exponents, unsigned);
    }

    monomial_basis_function(unsigned i, unsigned j)
    { m_exponents = list_of(i)(j); }

    monomial_basis_function(unsigned i, unsigned j, unsigned k)
    { m_exponents = list_of(i)(j)(k); }

    const double operator()(const hedge::vector &v) const
    {
      double result = 1;
      unsigned i = 0;
      BOOST_FOREACH(unsigned exp, m_exponents)
        result *= pow(v[i++], exp);

      return result;
    }
  };




  bool is_in_unit_simplex(const hedge::vector &unit_coords)
  {
    BOOST_FOREACH(hedge::vector::value_type ri, unit_coords)
      if (ri < -1)
        return false;

    return std::accumulate(unit_coords.begin(), unit_coords.end(), 
        (double) 0) <= -(signed(unit_coords.size())-2);
  }




  class zero_vector { };




  class interpolator
  {
    private:
      unsigned m_el_start, m_el_end;
      hedge::vector m_interpolation_coefficients;

    public:
      interpolator(unsigned el_start, unsigned el_end, 
          const hedge::vector &intp_coeff)
        : m_el_start(el_start), m_el_end(el_end),
        m_interpolation_coefficients(intp_coeff)
      { }

      const double operator()(const hedge::vector &data) const
      {
        return inner_prod(m_interpolation_coefficients,
            subrange(data, m_el_start, m_el_end));
      }

      const double operator()(zero_vector) const
      { return 0; }
  };




  class mesh_info : boost::noncopyable
  {
    public:
      // data structures ------------------------------------------------------
      typedef unsigned element_number;
      typedef unsigned vertex_number;
      typedef std::vector<element_number> el_id_vector;
      typedef std::vector<vertex_number> vtx_id_vector;

      static const element_number INVALID_ELEMENT = UINT_MAX;
      static const vertex_number INVALID_VERTEX = UINT_MAX;

      struct element_info
      {
        element_number                 m_id;
        hedge::affine_map              m_inverse_map;
        unsigned                       m_ldis_index;

        unsigned                       m_start, m_end;

        vtx_id_vector                  m_vertices;

        // the indices for the following two lists match up:
        // say at matching indices you find normal n
        // and element index i, then n is the normal
        // of the face leading to element i.

        std::vector<hedge::vector>      m_normals;
        std::vector<element_number>     m_neighbors;
      };

      struct local_discretization
      {
        std::vector<monomial_basis_function> m_basis;
        hedge::matrix m_l_mon_vandermonde_t;
        hedge::matrix m_u_mon_vandermonde_t;
        csr_matrix m_p_mon_vandermonde_t;
      };

      // setup ----------------------------------------------------------------
      mesh_info(
          unsigned dimensions,
          unsigned vertices_sizehint, 
          unsigned elements_sizehint, 
          unsigned discretizations_sizehint)
        : m_dimensions(dimensions)
      {
        m_vertices.reserve(vertices_sizehint);
        m_element_info.reserve(elements_sizehint);
        m_local_discretizations.reserve(discretizations_sizehint);
      }

      void add_local_discretization(python::list basis,
          const hedge::matrix &l_vdmt, 
          const hedge::matrix &u_vdmt, 
          const csr_matrix &p_vdmt)
      {
        local_discretization ldis;

        COPY_PY_LIST(ldis.m_basis, basis, const monomial_basis_function &);

        ldis.m_l_mon_vandermonde_t = l_vdmt;
        ldis.m_u_mon_vandermonde_t = u_vdmt;
        ldis.m_p_mon_vandermonde_t = p_vdmt;
        m_local_discretizations.push_back(ldis);
      }

      void add_element(const hedge::affine_map &inverse_map, 
          unsigned ldis_index, unsigned start, unsigned end,
          python::object vertices,
          python::object normals, 
          python::object neighbors
          )
      {
        element_info ei;
        ei.m_id = m_element_info.size();
        ei.m_inverse_map = inverse_map;
        ei.m_ldis_index = ldis_index;

        ei.m_start = start;
        ei.m_end = end;

        COPY_PY_LIST(ei.m_vertices, vertices, vertex_number);
        COPY_PY_LIST(ei.m_normals, normals, const hedge::vector &);
        COPY_PY_LIST(ei.m_neighbors, neighbors, element_number);

        m_element_info.push_back(ei);
      }

      void add_vertex(
          vertex_number vn, 
          const hedge::vector &pos,
          python::list adjacent_elements)
      {
        if (vn != m_vertex_adj_elements.size())
          throw std::runtime_error("vertices must be added in increasing order");

        m_vertices.push_back(pos);

        std::auto_ptr<el_id_vector> adj_els_cpp(new el_id_vector);
        COPY_PY_LIST(*adj_els_cpp, adjacent_elements, element_number);

        m_vertex_adj_elements.push_back(adj_els_cpp);
      }




      // operations -----------------------------------------------------------
      const bool is_in_element(element_number en, const hedge::vector &pt) const
      {
        const element_info &el(m_element_info[en]);
        hedge::vector uc = el.m_inverse_map(pt);
        return is_in_unit_simplex(uc);
      }

      const element_number find_containing_element(const hedge::vector &pt) const
      {
        BOOST_FOREACH(const element_info &el, m_element_info)
          if (is_in_unit_simplex(el.m_inverse_map(pt)))
            return el.m_id;
        return INVALID_ELEMENT;
      }

      const interpolator make_interpolator(
          const hedge::vector &pt, 
          element_number in_element) const
      {
        const element_info &el_inf = m_element_info[in_element];
        const local_discretization &ldis = 
          m_local_discretizations[el_inf.m_ldis_index];

        unsigned basis_length = el_inf.m_end-el_inf.m_start;
        hedge::vector mon_basis_values_at_pt(basis_length);

        for (unsigned i = 0; i < basis_length; i++)
          mon_basis_values_at_pt[i] = ldis.m_basis[i](pt);

        hedge::vector permuted_basis_values = prod(
                  ldis.m_p_mon_vandermonde_t,
                  mon_basis_values_at_pt);

        hedge::vector coeff = 
          solve(
              ldis.m_u_mon_vandermonde_t,
              solve(
                ldis.m_l_mon_vandermonde_t,
                permuted_basis_values,
                ublas::lower_tag()),
              ublas::upper_tag());

        return interpolator(el_inf.m_start, el_inf.m_end, coeff);
      }

      // data members ---------------------------------------------------------
      unsigned m_dimensions;

      hedge::vector m_nodes;
      std::vector<local_discretization> m_local_discretizations;
      std::vector<element_info> m_element_info;
      std::vector<hedge::vector> m_vertices;
      boost::ptr_vector<el_id_vector> m_vertex_adj_elements;
  };




  class particle_cloud : boost::noncopyable
  {
    public:
      // member data ----------------------------------------------------------
      mesh_info                         m_mesh_info;

      mesh_info::el_id_vector           m_containing_elements;
      hedge::vector                     m_positions;
      hedge::vector                     m_velocities;
      hedge::vector                     m_charges;
      hedge::vector                     m_masses;

      std::vector<unsigned>             m_deadlist;
      python::dict                      m_vis_info;

      // setup ----------------------------------------------------------------
      particle_cloud(
          unsigned dimensions,
          unsigned vertices_sizehint,
          unsigned elements_sizehint, 
          unsigned discretizations_sizehint)
        : m_mesh_info(dimensions, 
            vertices_sizehint,
            elements_sizehint, 
            discretizations_sizehint)
      {
      }

      // operation ------------------------------------------------------------
      void update_containing_elements()
      {
        for (unsigned i = 0; i < m_containing_elements.size(); i++)
        {
          unsigned pstart = i*m_mesh_info.m_dimensions;
          unsigned pend = (i+1)*m_mesh_info.m_dimensions;

          mesh_info::element_number prev = m_containing_elements[i];
          if (prev == mesh_info::INVALID_ELEMENT)
            continue;

          mesh_info::element_info &prev_el = m_mesh_info.m_element_info[prev];
          hedge::vector pt = subrange(m_positions, pstart, pend);

          // check if we're still in the same element -------------------------
          if (is_in_unit_simplex(prev_el.m_inverse_map(pt)))
            continue;

          // we're not: lookup via normal -------------------------------------
          {
            int closest_normal_idx = -1;
            double max_ip = 0;
            unsigned normal_idx = 0;

            BOOST_FOREACH(hedge::vector &n, prev_el.m_normals)
            {
              double ip = inner_prod(n, subrange(m_velocities, pstart, pend));
              if (ip > max_ip)
              {
                closest_normal_idx = normal_idx;
                max_ip = ip;
              }
              ++normal_idx;
            }

            if (closest_normal_idx == -1)
              throw std::runtime_error("no best normal found--weird");

            mesh_info::element_number possible_idx =
              prev_el.m_neighbors[closest_normal_idx];

            if (possible_idx != mesh_info::INVALID_ELEMENT)
            {
              mesh_info::element_info &possible = 
                m_mesh_info.m_element_info[possible_idx];

              if (is_in_unit_simplex(possible.m_inverse_map(pt)))
              {
                m_containing_elements[i] = possible.m_id;
                continue;
              }
            }
          }

          // look up via closest vertex ---------------------------------------
          {
            mesh_info::vertex_number closest_vertex = 
              mesh_info::INVALID_VERTEX;

            {
              double min_dist = std::numeric_limits<double>::infinity();

              BOOST_FOREACH(mesh_info::vertex_number vi, prev_el.m_vertices)
              {
                double dist = norm_2(m_mesh_info.m_vertices[vi] - pt);
                if (dist < min_dist)
                {
                  closest_vertex = vi;
                  min_dist = dist;
                }
              }
            }

            // found closest vertex, go through adjacent elements
            bool found = false;

            BOOST_FOREACH(mesh_info::element_number possible_idx, 
                m_mesh_info.m_vertex_adj_elements[closest_vertex])
            {
              mesh_info::element_info &possible = 
                m_mesh_info.m_element_info[possible_idx];

              if (is_in_unit_simplex(possible.m_inverse_map(pt)))
              {
                m_containing_elements[i] = possible.m_id;
                found = true;
                break;
              }
            }

            if (found)
              continue;
          }

          // last resort: global search ---------------------------------------
          m_containing_elements[i] = m_mesh_info.find_containing_element(pt);

          // no element found? kill the particle ------------------------------
          m_containing_elements[i] = mesh_info::INVALID_ELEMENT;
          m_deadlist.push_back(i);

          std::cout << "KILL" << i << std::endl;
          continue;
        }
      }

      // why all these template arguments? In 2D and 1D,
      // instead of passing a hedge::vector, you may simply
      // pass a zero_vector, and interpolation will know to
      // not even compute anything, but just return zero.
      template <class EX, class EY, class EZ, 
               class HX, class HY, class HZ>
      hedge::vector accelerations(
          const EX &ex, const EY &ey, const EZ &ez,
          const HX &hx, const HY &hy, const HZ &hz,
          bool update_vis_info
          )
      {
        hedge::vector result(m_positions.size());
        return result;
      }
  };



  // Python wrap helpers ------------------------------------------------------
  template <class Vec>
  void vector_extend(Vec &dest, const Vec &src)
  {
    std::copy(src.begin(), src.end(), back_inserter(dest));
  }
}




BOOST_PYTHON_MODULE(_internal)
{
  {
    typedef monomial_basis_function cl;
    python::class_<cl>("MonomialBasisFunction", python::init<unsigned, unsigned>())
      .def(python::init<unsigned, unsigned, unsigned>())
      .def(python::init<python::list>())
      .def("__call__", &cl::operator())
      ;
  }

  {
    typedef zero_vector cl;
    python::class_<cl>("ZeroVector");
  }

  {
    typedef mesh_info cl;
    python::class_<cl, boost::noncopyable>("MeshInfo", 
        python::init<unsigned, unsigned, unsigned, unsigned>())
      .def_readonly("INVALID_ELEMENT", &cl::INVALID_ELEMENT)
      .def_readwrite("dimensions", &cl::m_dimensions)
      .def_readwrite("nodes", &cl::m_nodes)

      .DEF_SIMPLE_METHOD(add_local_discretization)
      .DEF_SIMPLE_METHOD(add_element)
      .DEF_SIMPLE_METHOD(add_vertex)

      .DEF_SIMPLE_METHOD(is_in_element)
      .DEF_SIMPLE_METHOD(find_containing_element)
      ;
  }
  {
    typedef particle_cloud cl;
    python::class_<cl, boost::noncopyable>("ParticleCloud", 
        python::init<unsigned, unsigned, unsigned, unsigned>())
      .def_readonly("mesh_info", &cl::m_mesh_info)

      .def_readonly("containing_elements", &cl::m_containing_elements)

      .def_readwrite("positions", &cl::m_positions)
      .def_readwrite("velocities", &cl::m_velocities)
      .def_readwrite("charges", &cl::m_charges)
      .def_readwrite("masses", &cl::m_masses)

      .def_readonly("deadlist", &cl::m_deadlist)

      .def_readonly("vis_info", &cl::m_deadlist)

      .DEF_SIMPLE_METHOD(update_containing_elements)
      .def("accelerations", &cl::accelerations<
          hedge::vector, hedge::vector, hedge::vector,
          hedge::vector, hedge::vector, hedge::vector>)
      ;
  }

  {
    typedef std::vector<unsigned> cl;
    python::class_<cl, boost::noncopyable>("UnsignedVector")
      .def(python::vector_indexing_suite<cl>())
      .DEF_SIMPLE_METHOD(clear)
      .def("append", &cl::push_back)
      ;
  }
}
