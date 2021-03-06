// Pyrticle - Particle in Cell in Python
// Little bits of helpfulness
// Copyright (C) 2007 Andreas Kloeckner
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.





#ifndef _BADFJAH_PYRTICLE_TOOLS_HPP_INCLUDED
#define _BADFJAH_PYRTICLE_TOOLS_HPP_INCLUDED




#include <utility>
#include <functional>
#include <pyublas/numpy.hpp>
#include <pyublas/elementwise_op.hpp>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/python/slice.hpp>



#define FOR_ALL_SLICE_INDICES_PREP(SLICE, LEN) \
  Py_ssize_t \
    fsi__start, \
    fsi__stop, \
    fsi__step, \
    fsi__length; \
  if (PySlice_GetIndicesEx( \
        reinterpret_cast<PySliceObject *>((SLICE).ptr()), (LEN), \
      &fsi__start, &fsi__stop, &fsi__step, &fsi__length)) \
    throw boost::python::error_already_set(); \



#define FOR_ALL_SLICE_INDICES_LOOP \
  for (Py_ssize_t fsi__cnt = 0; fsi__cnt < fsi__length;  ++fsi__cnt)



#define FOR_ALL_SLICE_INDICES(SLICE, LEN) \
  FOR_ALL_SLICE_INDICES_PREP(SLICE, LEN) \
  FOR_ALL_SLICE_INDICES_LOOP

#define FOR_ALL_SLICE_INDICES_INNER(ITYPE, IVAR) \
  ITYPE IVAR = fsi__start + fsi__cnt*fsi__step;



namespace pyrticle 
{
  // common types -------------------------------------------------------------
  typedef unsigned particle_number;
  static const particle_number INVALID_PARTICLE = UINT_MAX;




  // vector / matrix types ----------------------------------------------------
  typedef pyublas::numpy_vector<int> py_int_vector;
  typedef pyublas::numpy_vector<double> py_vector;
  typedef pyublas::numpy_matrix<double> py_matrix;
  typedef pyublas::numpy_matrix<double,
          boost::numeric::ublas::column_major> py_fortran_matrix;
  typedef boost::numeric::ublas::vector<double> dyn_vector;
  typedef boost::numeric::ublas::matrix<double> dyn_matrix;
  typedef boost::numeric::ublas::matrix<double,
          boost::numeric::ublas::column_major> dyn_fortran_matrix;
  static const unsigned bounded_max_dims = 3;
  typedef boost::numeric::ublas::bounded_vector<double, bounded_max_dims> bounded_vector;
  typedef boost::numeric::ublas::bounded_vector<npy_int32, bounded_max_dims> bounded_int_vector;

  typedef boost::numeric::ublas::compressed_matrix<
    double, boost::numeric::ublas::column_major, 0, 
    boost::numeric::ublas::unbounded_array<int> >
      csr_matrix;
  typedef boost::numeric::ublas::zero_vector<double> zero_vector;




  // box utilties -------------------------------------------------------------
  template <class VecT>
  struct box
  {
    typedef VecT vector_type;
    vector_type m_lower, m_upper;

    box(vector_type const &lower, vector_type const &upper)
      : m_lower(lower), m_upper(upper)
    { }

    bool operator==(box const &b2) const
    {
      return 
        std::equal(m_lower.begin(), m_lower.end(), b2.m_lower.begin())
        &&
        std::equal(m_upper.begin(), m_upper.end(), b2.m_upper.begin());
    }

    bool operator!=(box const &b2) const
    { return !operator==(b2); }

    bool is_empty() const
    {
      for (unsigned i = 0; i < m_lower.size(); ++i)
        if (m_lower[i] >= m_upper[i])
          return true;
      return false;
    }

    template <class VecType>
    bool contains(VecType const &pt) const
    { return contains(pt, 1e-10); }

    template <class VecType>
    bool contains(VecType const &pt, double threshold=1e-10) const
    {
      for (unsigned i = 0; i < m_lower.size(); ++i)
        if (pt[i] < m_lower[i] - threshold 
            || pt[i] >= m_upper[i]+threshold)
          return false;
      return true;
    }

    template <class VecType2>
    box intersect(box<VecType2> const &b2) const
    {
      const unsigned dims = m_lower.size();

      const vector_type d_vector_lower(dims);
      const vector_type d_vector_upper(dims);
      box<vector_type> result(d_vector_lower, d_vector_upper);

      for (unsigned i = 0; i < dims; ++i)
      {
        result.m_lower[i] = std::max(m_lower[i], b2.m_lower[i]);
        result.m_upper[i] = std::min(m_upper[i], b2.m_upper[i]);
      }

      return result;
    }

    template <class VecType2>
    box enlarged(VecType2 const &vec) const
    { return box(m_lower-vec, m_upper+vec); }

    box enlarged(typename vector_type::value_type scalar) const
    { 
      return enlarged(
          boost::numeric::ublas::scalar_vector<typename vector_type::value_type>(
            m_lower.size(), scalar));
    }
  };



  typedef box<bounded_vector> bounded_box;
  typedef box<bounded_int_vector> bounded_int_box;




  // utilities ----------------------------------------------------------------
  template <class Map>
  typename Map::mapped_type const &map_get(
      Map const &map, 
      typename Map::key_type const &key)
  {
    typename Map::const_iterator it = map.find(key);
    if (it == map.end())
      throw std::logic_error("item not found in map");
    else
      return it->second;
  }




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




  template <class T>
  inline const T square(T x)
  {
    return x*x;
  }





  template <class VecType>
  inline typename VecType::value_type entry_or_zero(const VecType &v, unsigned i)
  {
    if (i >= v.size())
      return 0;
    else
      return v[i];
  }




  template <class T>
  inline T entry_or_zero(const T *v, int i)
  {
    return v[i];
  }




  template <class VecType1, class VecType2>
  inline
  const VecType1 cross(
      const VecType1 &a, 
      const VecType2 &b)
  {
    VecType1 result(3);
    result[0] = entry_or_zero(a,1)*entry_or_zero(b,2) - entry_or_zero(a,2)*entry_or_zero(b,1);
    result[1] = entry_or_zero(a,2)*entry_or_zero(b,0) - entry_or_zero(a,0)*entry_or_zero(b,2);
    result[2] = entry_or_zero(a,0)*entry_or_zero(b,1) - entry_or_zero(a,1)*entry_or_zero(b,0);
    return result;
  }




  template <class T> 
  struct identity : std::unary_function<T,T> {
    const T operator()(const T& x) const
    {
      return x;
    }
  };




  template <class InputIterator, class Function>
  inline double average(InputIterator first, InputIterator last, 
      Function fn = identity<double>())
  {
    double result = 0;
    unsigned count = 0;

    BOOST_FOREACH(double value, std::make_pair(first,last))
    {
      result += fn(value);
      ++count;
    }

    if (count == 0)
      throw std::runtime_error("attempted to take empty average");

    return result/count;
  }





  template <class InputIterator>
  inline double std_dev(InputIterator first, InputIterator last)
  {
    double square_sum = 0;
    double sum = 0;
    unsigned count = 0;

    BOOST_FOREACH(double value, std::make_pair(first,last))
    {
      sum += value;
      square_sum += square(value);
      ++count;
    }

    if (count == 0)
      throw std::runtime_error("attempted to take empty average");

    double mean = sum/count;
    return sqrt(square_sum/count - square(mean));
  }




  template <class VectorType>
  bool isnan_any(VectorType const &vec)
  {
    BOOST_FOREACH(typename VectorType::value_type value, vec)
      if (isnan(value))
        return true;
    return false;
  }




  template <class T>
  class stats_gatherer
  {
    // http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    private:
      unsigned m_count;
      T m_mean, m_m2;
      T m_min, m_max;

    public:
      stats_gatherer()
        : m_count(0), m_mean(0), m_m2(0)
      { }

      void add(T x)
      {
        if (m_count == 0 || x < m_min)
          m_min = x;

        if (m_count == 0 || x > m_max)
          m_max = x;

        ++m_count;
        T delta = x - m_mean;
        m_mean += delta/T(m_count);
        m_m2 += delta*(x - m_mean);
      }

      void reset()
      {
        m_count = 0;
        m_mean = 0;
        m_m2 = 0;
      }

      unsigned count() const
      { return m_count; }

      T minimum() const
      { return m_min; }

      T maximum() const
      { return m_max; }

      T mean() const
      {
        if (m_count == 0)
          throw std::runtime_error("attempted to take empty mean");

        return m_mean;
      }

      T variance_sample() const
      {
        if (m_count < 2)
          throw std::runtime_error("sample too small for sample variance");

        return m_m2/(m_count-1);
      }

      T standard_deviation_sample() const
      {
        return sqrt(variance_sample());
      }

      T variance() const
      {
        if (m_count == 0)
          throw std::runtime_error("attempted to take empty variance");

        return m_m2/m_count;
      }

      T standard_deviation() const
      {
        return sqrt(variance());
      }

  };




  class visualization_listener
  {
    public:
      virtual ~visualization_listener()
      { }

      virtual void store_mesh_vis_vector(
          const char *name, const py_vector &vec) const = 0;
      virtual void store_particle_vis_vector(
          const char *name,
          const py_vector &vec) const = 0;
  };




  class number_shift_listener
  {
    public:
      virtual ~number_shift_listener()
      { }

      virtual void note_change_size(unsigned new_size) const 
      { }
      virtual void note_move(unsigned orig, unsigned dest, unsigned size) const 
      { }
      virtual void note_reset(unsigned start, unsigned size) const 
      { }
  };



  class boundary_hit_listener
  {
    public:
      virtual ~boundary_hit_listener()
      { }

      virtual void note_boundary_hit(particle_number pn) const 
      { }
  };




  class warning_listener
  {
    private:
      static warning_listener   *m_singleton;

    public:
      warning_listener()
      {
        if (m_singleton)
          throw std::runtime_error("warning listener singleton already exists");
        m_singleton = this;
      }

      virtual ~warning_listener()
      { 
        m_singleton = 0;
      }

      static void warn(
          std::string const &message,
          std::string const &filename,
          unsigned lineno
          ) 
      {
        if (m_singleton)
          m_singleton->note_warning(message, filename, lineno);
        else
          throw std::runtime_error("warning raised, but no listener registered");
      }

      virtual void note_warning(
          std::string const &message,
          std::string const &filename,
          unsigned lineno
          ) const = 0;
  };




#define WARN(MESSAGE) \
  warning_listener::warn(MESSAGE, __FILE__, __LINE__);




  // shape functions ----------------------------------------------------------
  class polynomial_shape_function
  {
    public:
      polynomial_shape_function(double radius=1, unsigned dimensions=1, double alpha=2);

      template <class VecType>
      const double operator()(const VecType &r) const
      {
        const double r_squared = pyublas::square_sum(r);
        if (r_squared > m_radius_squared)
          return 0;
        else
        {
          double radius_term = m_radius-r_squared/m_radius;

          if (m_alpha_is_2)
            return m_normalizer*radius_term*radius_term;
          else
            return m_normalizer * pow(radius_term, m_alpha);
        }
      }

      const double normalizer() const
      { return m_normalizer; }

      const double radius() const
      { return m_radius; }

      const double exponent() const
      { return m_alpha; }

      static const std::string name()
      {
        return "polynomial";
      }

    private:
      double m_normalizer;
      double m_alpha;
      double m_radius, m_radius_squared;

      bool m_alpha_is_2;
  };




  class c_infinity_shape_function
  {
    public:
      c_infinity_shape_function()
      { }

      c_infinity_shape_function(double radius, unsigned dimensions, 
          double integral_for_rad1);

      template <class VecType>
      const double operator()(const VecType &r) const
      {
        const double r_squared = pyublas::square_sum(r);
        if (r_squared > m_radius_squared)
          return 0;
        else
        {
          double sr_squared_m_1 = r_squared/m_radius_squared-1;
          return m_normalizer*exp(-1/(sr_squared_m_1*sr_squared_m_1));
        }
      }

      const double radius() const
      { return m_radius; }

      static const std::string name()
      {
        return "c_infinity";
      }

    private:
      double m_normalizer;
      double m_radius, m_radius_squared;
  };




}




#define PIC_THIS static_cast<PICAlgorithm *>(this)
#define CONST_PIC_THIS static_cast<const PICAlgorithm *>(this)




#endif
