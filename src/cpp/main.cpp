#include <stdexcept>
#include <vector>
#include <list>
#include <climits>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/acosh.hpp>
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
#define DEF_RO_MEMBER(NAME) \
  def_readonly(#NAME, &cl::m_##NAME)
#define DEF_RW_MEMBER(NAME) \
  def_readwrite(#NAME, &cl::m_##NAME)




namespace {
  typedef ublas::compressed_matrix<
    double, ublas::column_major, 0, 
    ublas::unbounded_array<int> >
      csr_matrix;
  typedef ublas::zero_vector<
    hedge::vector::value_type>
    zero_vector;
  typedef unsigned particle_number;





  template <class T>
  PyObject *manage_new_object(T *obj)
  {
    typename python::manage_new_object::apply<T *>::type 
      result_converter;
    return result_converter(obj);
  }





  template <class T>
  inline const T square(T x)
  {
    return x*x;
  }





  template <class VecType>
  inline hedge::vector::value_type entry_or_zero(const VecType &v, int i)
  {
    if (i >= v.size())
      return 0;
    else
      return v[i];
  }




  template <class VecType1, class VecType2>
  inline
  const hedge::vector cross(
      const VecType1 &a, 
      const VecType2 &b)
  {
    hedge::vector result(3);
    result[0] = entry_or_zero(a,1)*entry_or_zero(b,2) - entry_or_zero(a,2)*entry_or_zero(b,1);
    result[1] = entry_or_zero(a,2)*entry_or_zero(b,0) - entry_or_zero(a,0)*entry_or_zero(b,2);
    result[2] = entry_or_zero(a,0)*entry_or_zero(b,1) - entry_or_zero(a,1)*entry_or_zero(b,0);
    return result;
  }




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




  class shape_function
  {
    public:
      shape_function(
          double radius,
          unsigned dimensions,
          double alpha=2)
        : m_alpha(alpha), m_l(radius), 
        m_l_squared(square(radius))
      {
        using boost::math::tgamma;
        using boost::math::beta;

        double n = dimensions;

        // see doc/notes.tm
        double sphere_area = 2*pow(M_PI, n/2) / tgamma(n/2);
        m_normalizer = 2 / (sphere_area *
          pow(m_l, n+alpha)*beta(n/2, alpha+1)
          );
      }

      const double operator()(const hedge::vector &r) const
      {
        double r_squared = inner_prod(r, r);
        if (r_squared > m_l_squared)
          return 0;
        else
          return m_normalizer * pow(m_l-r_squared/m_l, m_alpha);
      }

    private:
      double m_normalizer;
      double m_alpha;
      double m_l, m_l_squared;
  };




  /** The ReconstructionTarget protocol:
   *
   * template <class Scaler>
   * class reconstruction_target
   * {
   *   void begin_particle(particle_number pn);
   *   void add_shape_at_point(unsigned i, double shape_factor)
   * };
   *
   * Note: this is a stateful protocol.
   */

  class rho_reconstruction_target
  {
    private:
      hedge::vector m_target_vector;
      const hedge::vector &m_charges;
      double m_scale_factor;

    public:
      rho_reconstruction_target(unsigned points, const hedge::vector &charges)
        : m_target_vector(points), m_charges(charges)
      { 
        m_target_vector.clear();
      }

      void begin_particle(particle_number pn)
      {
        m_scale_factor = m_charges[pn];
      }

      void add_shape_at_point(unsigned i, double shape_factor)
      {
        m_target_vector[i] += shape_factor * m_scale_factor;
      }

      const hedge::vector &result() const
      {
        return m_target_vector;
      }
  };




  /** Reconstruction Target for the current density.
   */
  template<unsigned dimensions_velocity>
  class j_reconstruction_target
  {
    private:
      hedge::vector m_target_vector;
      const hedge::vector &m_charges;
      const hedge::vector &m_velocities;
      double m_scale_factors[dimensions_velocity];

    public:
      j_reconstruction_target(unsigned points, const hedge::vector &charges,
          const hedge::vector &velocities)
        : m_target_vector(3*points), m_charges(charges), m_velocities(velocities)
      { 
        m_target_vector.clear();
      }

      void begin_particle(particle_number pn)
      {
        const double charge = m_charges[pn];
        for (unsigned axis = 0; axis < dimensions_velocity; axis++)
          m_scale_factors[axis] = charge * m_velocities[pn*dimensions_velocity+axis];
      }

      void add_shape_at_point(unsigned i, double shape_factor)
      {
        const unsigned base = i*dimensions_velocity;
        for (unsigned axis = 0; axis < dimensions_velocity; axis++)
          m_target_vector[base+axis] += shape_factor * m_scale_factors[axis];
      }

      const hedge::vector &result() const
      {
        return m_target_vector;
      }
  };




  template <class T1, class T2>
  class chained_reconstruction_target
  {
    private:
      T1 &m_target1;
      T2 &m_target2;

    public:
      chained_reconstruction_target(T1 &target1, T2 &target2)
        : m_target1(target1), m_target2(target2)
      { }

      void begin_particle(particle_number pn)
      {
        m_target1.begin_particle(pn);
        m_target2.begin_particle(pn);
      }

      void add_shape_at_point(unsigned i, double shape_factor)
      {
        m_target1.add_shape_at_point(i, shape_factor);
        m_target2.add_shape_at_point(i, shape_factor);
      }
  };



  template <class T1, class T2>
  inline
  chained_reconstruction_target<T1, T2> 
  make_chained_reconstruction_target(T1 target1, T2 target2)
  {
    return chained_reconstruction_target<T1, T2>(target1, target2);
  }




  const bool is_in_unit_simplex(const hedge::vector &unit_coords)
  {
    const double eps = 1e-10;

    BOOST_FOREACH(hedge::vector::value_type ri, unit_coords)
      if (ri < -1-eps)
        return false;

    return std::accumulate(unit_coords.begin(), unit_coords.end(), 
        (double) 0) <= -(signed(unit_coords.size())-2)+eps;
  }




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

      // data members ---------------------------------------------------------
      const unsigned m_dimensions;

      std::vector<local_discretization> m_local_discretizations;
      std::vector<element_info> m_element_info;
      std::vector<hedge::vector> m_vertices, m_nodes;
      boost::ptr_vector<el_id_vector> m_vertex_adj_elements;




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




      void add_nodes(unsigned sizehint, python::object iterable)
      {
        m_nodes.reserve(sizehint);
        python::stl_input_iterator<const hedge::vector &> 
          first(iterable), last;
        std::copy(first, last, std::back_inserter(m_nodes));
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

        hedge::vector unit_pt = el_inf.m_inverse_map(pt);

        for (unsigned i = 0; i < basis_length; i++)
          mon_basis_values_at_pt[i] = ldis.m_basis[i](unit_pt);

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
  };




  class event_counter
  {
    private:
      unsigned          m_count;

    public:
      event_counter()
        : m_count(0)
        { }

      unsigned get()
      { return m_count; }

      unsigned pop()
      { 
        unsigned result = m_count; 
        m_count = 0;
        return result;
      }

      void tick()
      { ++m_count; }
  };




  template<unsigned dimensions_pos, unsigned dimensions_velocity>
  class particle_cloud : boost::noncopyable
  {
    public:
      // member data ----------------------------------------------------------
      mesh_info                         m_mesh_info;

      static const unsigned             m_dimensions_pos = dimensions_pos;
      static const unsigned             m_dimensions_velocity = dimensions_velocity;

      mesh_info::el_id_vector           m_containing_elements;
      hedge::vector                     m_positions;
      hedge::vector                     m_momenta;
      hedge::vector                     m_charges;
      hedge::vector                     m_masses;

      std::vector<particle_number>      m_deadlist;
      python::dict                      m_vis_info;

      mutable event_counter             m_same_searches,
                                        m_normal_searches,
                                        m_vertex_searches,
                                        m_global_searches,
                                        m_vertex_shape_adds,
                                        m_neighbor_shape_adds,
                                        m_periodic_hits;

      const double                      m_epsilon, m_mu, m_c;

      python::object                    m_boundary_hit_cb;




      // setup ----------------------------------------------------------------
      particle_cloud(
          unsigned dimensions_mesh,
          unsigned vertices_sizehint,
          unsigned elements_sizehint, 
          unsigned discretizations_sizehint,
          double epsilon, 
          double mu,
          python::object boundary_hit_cb)
        : m_mesh_info(dimensions_mesh, 
            vertices_sizehint,
            elements_sizehint, 
            discretizations_sizehint),
        m_epsilon(epsilon), m_mu(mu), m_c(1/sqrt(mu*epsilon)),
        m_boundary_hit_cb(boundary_hit_cb)
      {
      }




      // operation ------------------------------------------------------------
      const hedge::vector velocities() const
      {
        const unsigned vdim = m_dimensions_velocity;

        hedge::vector result(m_momenta.size());

        for (particle_number pn = 0; pn < m_containing_elements.size(); pn++)
        {
          unsigned vpstart = vdim*pn;
          unsigned vpend = vdim*(pn+1);

          mesh_info::element_number in_el = m_containing_elements[pn];
          if (in_el != mesh_info::INVALID_ELEMENT)
          {
            const double m = m_masses[pn];
            double p = norm_2(subrange(m_momenta, vpstart, vpend));
            double v = m_c*p/sqrt(m*m*m_c*m_c + p*p);
            subrange(result, vpstart, vpend) = v/p*subrange(m_momenta, vpstart, vpend);
          }
        }
        return result;
      }


      

      mesh_info::element_number find_new_containing_element(particle_number i,
          mesh_info::element_number prev) const
      {
        const unsigned xdim = m_dimensions_pos;

        const unsigned x_pstart = i*xdim;
        const unsigned x_pend = (i+1)*xdim;
        
        const hedge::vector pt = subrange(m_positions, x_pstart, x_pend);

        if (prev != mesh_info::INVALID_ELEMENT)
        {
          const mesh_info::element_info &prev_el = m_mesh_info.m_element_info[prev];

          // check if we're still in the same element -------------------------
          if (is_in_unit_simplex(prev_el.m_inverse_map(pt)))
          {
            m_same_searches.tick();
            return prev;
          }

          // we're not: lookup via normal -------------------------------------
          {
            int closest_normal_idx = -1;
            double max_ip = 0;
            unsigned normal_idx = 0;

            BOOST_FOREACH(const hedge::vector &n, prev_el.m_normals)
            {
              double ip = inner_prod(n, subrange(m_momenta, x_pstart, x_pend));
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
              const mesh_info::element_info &possible = 
                m_mesh_info.m_element_info[possible_idx];

              if (is_in_unit_simplex(possible.m_inverse_map(pt)))
              {
                m_normal_searches.tick();
                return possible.m_id;
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
            BOOST_FOREACH(mesh_info::element_number possible_idx, 
                m_mesh_info.m_vertex_adj_elements[closest_vertex])
            {
              const mesh_info::element_info &possible = 
                m_mesh_info.m_element_info[possible_idx];

              if (is_in_unit_simplex(possible.m_inverse_map(pt)))
              {
                m_vertex_searches.tick();
                return possible.m_id;
              }
            }
          }

        }

        // last resort: global search ---------------------------------------
        {
          m_global_searches.tick();

          mesh_info::element_number new_el = 
            m_mesh_info.find_containing_element(pt);
          if (new_el != mesh_info::INVALID_ELEMENT)
            return new_el;
        }

        return mesh_info::INVALID_ELEMENT;
      }




      void update_containing_elements()
      {
        for (particle_number i = 0; i < m_containing_elements.size(); i++)
        {
          mesh_info::element_number prev = m_containing_elements[i];
          if (prev == mesh_info::INVALID_ELEMENT)
            continue;

          mesh_info::element_number new_el = 
            find_new_containing_element(i, prev);

          // no element found? delegate to callback ---------------------------
          if (new_el == mesh_info::INVALID_ELEMENT)
            m_boundary_hit_cb(i);
          else
            m_containing_elements[i] = new_el;
        }
      }




      // why all these template arguments? In 2D and 1D,
      // instead of passing a hedge::vector, you may simply
      // pass a zero_vector, and interpolation will know to
      // not even compute anything, but just return zero.
      template <class EX, class EY, class EZ, 
               class HX, class HY, class HZ>
      hedge::vector forces(
          const EX &ex, const EY &ey, const EZ &ez,
          const HX &hx, const HY &hy, const HZ &hz,
          const hedge::vector &velocities,
          bool update_vis_info
          )
      {
        hedge::vector result(m_momenta.size());
        std::auto_ptr<hedge::vector> 
          vis_e, vis_h, vis_el_force, vis_lorentz_force;

        if (update_vis_info)
        {
          vis_e = std::auto_ptr<hedge::vector>(
              new hedge::vector(3*m_containing_elements.size()));
          vis_h = std::auto_ptr<hedge::vector>(
              new hedge::vector(3*m_containing_elements.size()));
          vis_el_force = std::auto_ptr<hedge::vector>(
              new hedge::vector(3*m_containing_elements.size()));
          vis_lorentz_force = std::auto_ptr<hedge::vector>(
              new hedge::vector(3*m_containing_elements.size()));
        }

        for (particle_number i = 0; i < m_containing_elements.size(); i++)
        {
          unsigned x_pstart = m_dimensions_pos*i;
          unsigned x_pend = m_dimensions_pos*(i+1);
          unsigned v_pstart = m_dimensions_velocity*i;
          unsigned v_pend = m_dimensions_velocity*(i+1);

          mesh_info::element_number in_el = m_containing_elements[i];
          if (in_el == mesh_info::INVALID_ELEMENT)
          {
            subrange(result, v_pstart, v_pend) = zero_vector(m_dimensions_pos);
            continue;
          }

          interpolator interp = m_mesh_info.make_interpolator(
              subrange(m_positions, x_pstart, x_pend), in_el);

          hedge::vector e(3);
          e[0] = interp(ex);
          e[1] = interp(ey);
          e[2] = interp(ez);

          hedge::vector h(3);
          h[0] = interp(hx);
          h[1] = interp(hy);
          h[2] = interp(hz);

          const double charge = m_charges[i];

          hedge::vector el_force(3);
          el_force[0] = charge*e[0];
          el_force[1] = charge*e[1];
          el_force[2] = charge*e[2];

          const hedge::vector v = subrange(velocities, v_pstart, v_pend);
          hedge::vector lorentz_force = cross(v, charge*m_mu*h);

          // truncate forces to m_dimensions_velocity entries
          subrange(result, v_pstart, v_pend) = subrange(
              el_force + lorentz_force, 0, m_dimensions_velocity);

          if (update_vis_info)
          {
            subrange(*vis_e, 3*i, 3*(i+1)) = e;
            subrange(*vis_h, 3*i, 3*(i+1)) = h;
            subrange(*vis_el_force, 3*i, 3*(i+1)) = el_force;
            subrange(*vis_lorentz_force, 3*i, 3*(i+1)) = lorentz_force;
          }
        }

        if (update_vis_info)
        {
          m_vis_info["pt_e"] = python::object(*vis_e);
          m_vis_info["pt_h"] = python::object(*vis_h);
          m_vis_info["el_force"] = python::object(*vis_el_force);
          m_vis_info["lorentz_force"] = python::object(*vis_lorentz_force);
        }

        return result;
      }




      template <class Target, class ShapeFunction>
      void add_shape_on_element(
          Target &tgt, 
          const hedge::vector &center,
          mesh_info::element_number en,
          const ShapeFunction &sf) const
      {
        const mesh_info::element_info &el = 
          m_mesh_info.m_element_info[en];

        for (unsigned i = el.m_start; i < el.m_end; i++)
          tgt.add_shape_at_point(i, sf(m_mesh_info.m_nodes[i]-center));
      }




      template <class Target, class ShapeFunction>
      void add_shape_by_neighbors(
          Target &target,
          const ShapeFunction &sf,
          particle_number pi) const
      {
        const unsigned dim = m_dimensions_pos;
        const hedge::vector pos = subrange(
            m_positions, pi*dim, (pi+1)*dim);
        const mesh_info::element_info &el(
            m_mesh_info.m_element_info[m_containing_elements[pi]]);

        add_shape_on_element(
            target, pos, m_containing_elements[pi], sf);
        BOOST_FOREACH(mesh_info::element_number en, el.m_neighbors)
          if (en != mesh_info::INVALID_ELEMENT)
            add_shape_on_element(target, pos, en, sf);
      }




      template <class Target, class ShapeFunction>
      void add_shape(
          Target &target, 
          const ShapeFunction &sf,
          double radius, particle_number pi) const
      {
        const unsigned dim = m_dimensions_pos;
        const hedge::vector pos = subrange(
            m_positions, pi*dim, (pi+1)*dim);
        const mesh_info::element_info &el(
            m_mesh_info.m_element_info[m_containing_elements[pi]]);

        // We're deciding between RULE A and RULE B below.
        // The decision is made by looking at the barycentric coordinates of
        // the center of the shape function. If that center has all barycentric
        // coordinates 1/2 or smaller, then use the faster RULE A.
        //
        // For this purpose, observe that the unit coordinates are the first
        // barycentric coordinates, and the remaining one is found by calculating
        // 1-sum(lambda_i)
        //
        // Note also that this is not purely a speed tradeoff. RULE B works fairly
        // well alone, but falls flat on its face for near the center of the hypotenuse
        // of a right triangle.

        // FIXME this assumes m_dimension_pos == m_dimension_mesh
        hedge::vector unit_pt = el.m_inverse_map(pos);
        
        bool can_use_rule_a = true;

        double uc_sum = 0;
        BOOST_FOREACH(double uc, unit_pt)
        {
          if (uc > 0)
          {
            can_use_rule_a = false;
            break;
          }
          uc_sum += uc;
        }

        if (can_use_rule_a && (1-0.5*(uc_sum+unit_pt.size()) < 0.5))
        {
          // RULE A: we're far enough away from vertices, 
          // weight onto neighboring elements
          m_neighbor_shape_adds.tick();

          add_shape_by_neighbors(target, sf, pi);
        }
        else
        {
          // RULE B: we're close to a vertex, weight onto all elements
          // adjoining that vertex

          // find closest vertex
          mesh_info::vertex_number closest_vertex = 
            mesh_info::INVALID_VERTEX;
          double min_dist = std::numeric_limits<double>::infinity();

          BOOST_FOREACH(mesh_info::vertex_number vi, el.m_vertices)
          {
            double dist = norm_2(m_mesh_info.m_vertices[vi] - pos);
            if (dist < min_dist)
            {
              closest_vertex = vi;
              min_dist = dist;
            }
          }
          // found a close vertex, go through adjacent elements
          m_vertex_shape_adds.tick();

          BOOST_FOREACH(mesh_info::element_number en, 
              m_mesh_info.m_vertex_adj_elements[closest_vertex])
            add_shape_on_element(target, pos, en, sf);
        }
      }




      template<class Target>
      void reconstruct_densities_on_target(Target &tgt, double radius) const
      {
        const shape_function sf(radius, m_mesh_info.m_dimensions);

        particle_number pn = 0;

        BOOST_FOREACH(mesh_info::element_number en,
            m_containing_elements)
        {
          tgt.begin_particle(pn);
          if (en != mesh_info::INVALID_ELEMENT)
            add_shape(tgt, sf, radius, pn);
          ++pn;
        }
      }




      void reconstruct_densities(
          hedge::vector &rho, 
          python::object py_j,
          double radius,
          const hedge::vector &velocities) const
      {
        rho_reconstruction_target rho_tgt(m_mesh_info.m_nodes.size(), m_charges);
        j_reconstruction_target<m_dimensions_velocity> j_tgt(m_mesh_info.m_nodes.size(), 
            m_charges, velocities);

        chained_reconstruction_target
          <rho_reconstruction_target, j_reconstruction_target<m_dimensions_velocity> >
          tgt(rho_tgt, j_tgt);
        reconstruct_densities_on_target(tgt, radius);

        rho = rho_tgt.result();

        int i = 0;
        BOOST_FOREACH(hedge::vector &j_i, 
            std::make_pair(
              python::stl_input_iterator<hedge::vector &>(py_j),
              python::stl_input_iterator<hedge::vector &>()))
        {
          if (i >= m_dimensions_velocity)
            throw std::runtime_error("j field must have v_dim dimensions");
          j_i = subslice(j_tgt.result(), 
              i++, m_dimensions_velocity, m_mesh_info.m_nodes.size());
        }
        if (i < m_dimensions_velocity)
            throw std::runtime_error("j field must have v_dim dimensions");
      }




      void reconstruct_rho(hedge::vector &rho, double radius) const
      {
        rho_reconstruction_target rho_tgt(m_mesh_info.m_nodes.size(), m_charges);

        reconstruct_densities_on_target(rho_tgt, radius);

        rho = rho_tgt.result();
      }
  };




  // Python wrap helpers ------------------------------------------------------
  template <class Cloud>
  void expose_particle_cloud(const std::string &name)
  {
    using python::arg;

    typedef Cloud cl;

    python::class_<cl, boost::noncopyable> 
      cloud_wrap(name.c_str(), 
          python::init<
          unsigned, 
          unsigned, unsigned, unsigned, 
          double, double, python::object>());
    cloud_wrap
      .DEF_RO_MEMBER(mesh_info)

      .DEF_RO_MEMBER(dimensions_pos)
      .DEF_RO_MEMBER(dimensions_velocity)

      .DEF_RO_MEMBER(containing_elements)
      .DEF_RW_MEMBER(positions)
      .DEF_RW_MEMBER(momenta)
      .DEF_RW_MEMBER(charges)
      .DEF_RW_MEMBER(masses)

      .DEF_RO_MEMBER(deadlist)

      .DEF_RO_MEMBER(vis_info)

      .DEF_RO_MEMBER(same_searches)
      .DEF_RO_MEMBER(normal_searches)
      .DEF_RO_MEMBER(vertex_searches)
      .DEF_RO_MEMBER(global_searches)
      .DEF_RO_MEMBER(vertex_shape_adds)
      .DEF_RO_MEMBER(neighbor_shape_adds)
      .DEF_RO_MEMBER(periodic_hits)

      .DEF_RO_MEMBER(epsilon)
      .DEF_RO_MEMBER(mu)
      .DEF_RO_MEMBER(c)

      .DEF_SIMPLE_METHOD(velocities)
      .DEF_SIMPLE_METHOD(find_new_containing_element)
      .DEF_SIMPLE_METHOD(update_containing_elements)
      ;

    if (Cloud::m_dimensions_velocity == 3)
    {
      cloud_wrap
        .def("forces", &cl::template forces< // full-field case
            hedge::vector, hedge::vector, hedge::vector,
            hedge::vector, hedge::vector, hedge::vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("hx"), arg("hy"), arg("hz"),
             arg("velocities"), arg("verbose_vis")))
        ;
    }
    else if (Cloud::m_dimensions_velocity == 2)
    {
      cloud_wrap
        .def("forces", &cl::template forces< // TM case
            zero_vector, zero_vector, hedge::vector,
            hedge::vector, hedge::vector, zero_vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("hx"), arg("hy"), arg("hz"),
             arg("velocities"), arg("verbose_vis")))
        .def("forces", &cl::template forces< // TE case
            hedge::vector, hedge::vector, zero_vector,
            zero_vector, zero_vector, hedge::vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("hx"), arg("hy"), arg("hz"),
             arg("velocities"), arg("verbose_vis")))
        ;
    }

    cloud_wrap
      .DEF_SIMPLE_METHOD(reconstruct_densities)
      .DEF_SIMPLE_METHOD(reconstruct_rho)
      ;
  }
}




BOOST_PYTHON_MODULE(_internal)
{
  using python::arg;

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
      .DEF_RW_MEMBER(nodes)
      .DEF_RO_MEMBER(dimensions)

      .DEF_SIMPLE_METHOD(add_local_discretization)
      .DEF_SIMPLE_METHOD(add_element)
      .DEF_SIMPLE_METHOD(add_vertex)
      .DEF_SIMPLE_METHOD(add_nodes)

      .DEF_SIMPLE_METHOD(is_in_element)
      .DEF_SIMPLE_METHOD(find_containing_element)
      ;
  }

  expose_particle_cloud<
    particle_cloud<3,3> >("ParticleCloud33");
  expose_particle_cloud<
    particle_cloud<2,2> >("ParticleCloud22");

  {
    typedef event_counter cl;
    python::class_<cl>("EventCounter")
      .DEF_SIMPLE_METHOD(get)
      .DEF_SIMPLE_METHOD(pop)
      .DEF_SIMPLE_METHOD(tick)
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

  python::def("asinh", (double (*)(double)) boost::math::asinh);
  python::def("acosh", (double (*)(double)) boost::math::acosh);
}
