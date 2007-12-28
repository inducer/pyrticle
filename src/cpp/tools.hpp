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




#include <hedge/base.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>




namespace pyrticle 
{
  // common types -------------------------------------------------------------
  typedef unsigned particle_number;

#define PIC_THIS static_cast<PICAlgorithm *>(this)
#define CONST_PIC_THIS static_cast<const PICAlgorithm *>(this)




  // common ublas types -------------------------------------------------------
  namespace ublas = boost::numeric::ublas;

  typedef ublas::compressed_matrix<
    double, ublas::column_major, 0, 
    ublas::unbounded_array<int> >
      csr_matrix;
  typedef ublas::zero_vector<
    hedge::vector::value_type>
    zero_vector;




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
}




#endif
