// Pyrticle - Particle in Cell in Python
// Python wrapper for PIC algorithm
// Copyright (C) 2007 Andreas Kloeckner
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or // (at your option) any later version.  // 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.





#include "wrap_pic.hpp"
#include "dep_shape.hpp"
#include "dep_normshape.hpp"
#include "dep_advective.hpp"
#include "dep_grid.hpp"
#include "dep_grid_find.hpp"




using namespace pyrticle;
using namespace boost::python;




namespace
{
  typedef polynomial_shape_function used_shape_function;
  // typedef c_infinity_shape_function used_shape_function;




  template <class Depositor>
  void expose_deposition_functions()
  {
    def("deposit_densities", deposit_densities<Depositor>);
    def("deposit_j", deposit_j<Depositor>);
    def("deposit_rho", deposit_rho<Depositor>);
  }




  template <class GridDep>
  py_vector get_extra_points(const GridDep &dep)
  { return dep.m_extra_points; }

  template <class GridDep>
  void set_extra_points(GridDep &dep, py_vector v)
  { dep.m_extra_points = v; }




  template <class ParticleState, class Brick, class GridDepBaseStateWrapper>
  void expose_grid_depositor(const std::string &brick_type, 
      GridDepBaseStateWrapper &gdbs_wrap)
  {
    typedef grid_depositor<ParticleState, used_shape_function, Brick> cl;
    class_<cl> wrp(
        (brick_type+"GridDepositor"+get_state_class_suffix<ParticleState>()).c_str(), 
        init<const mesh_data &>());

    wrp
      .DEF_RW_MEMBER(shape_function)

      .DEF_RW_MEMBER(bricks)
      .DEF_RW_MEMBER(elements_on_grid)

      .DEF_RW_MEMBER(first_extra_point)
      // PyUblas member-in-base-class issue: wrap by hand
      .add_property("extra_points", get_extra_points<cl>, set_extra_points<cl>)
      .DEF_RW_MEMBER(extra_point_brick_starts)

      .DEF_RW_MEMBER(average_groups)
      .DEF_RW_MEMBER(average_group_starts)

      .DEF_SIMPLE_METHOD(grid_node_count)
      .DEF_SIMPLE_METHOD(grid_node_count_with_extra)
      .DEF_SIMPLE_METHOD(find_points_in_element)

      .DEF_SIMPLE_METHOD(remap_grid_to_mesh)
      .DEF_SIMPLE_METHOD(remap_residual)

      .DEF_SIMPLE_METHOD(deposit_grid_densities)
      .DEF_SIMPLE_METHOD(deposit_grid_j)
      .DEF_SIMPLE_METHOD(deposit_grid_rho)
      ;

    wrp.attr("DepositorState") = gdbs_wrap;
  }




  template <class ParticleState, class GridDepBaseStateWrapper>
  void expose_depositors_for_pstate(GridDepBaseStateWrapper &gdbs_wrap)
  {
    {
      typedef shape_function_depositor<ParticleState, used_shape_function> cl;
      class_<cl> wrp(
        ("InterpolatingDepositor"+get_state_class_suffix<ParticleState>()).c_str(), 
        init<const mesh_data &>());

      wrp
        .DEF_RW_MEMBER(shape_function)
        ;

      scope cls_scope = wrp;
      {
        typedef typename cl::depositor_state cl;
        class_<cl> wrp("DepositorState")
          ;
      }
    }

    {
      typedef normalized_shape_function_depositor<
        ParticleState, used_shape_function> cl;
      class_<cl> wrp(
        ("NormalizingInterpolatingDepositor"+get_state_class_suffix<ParticleState>()).c_str(), 
        init<const mesh_data &, const py_matrix &>());

      wrp
        .DEF_RW_MEMBER(shape_function)
        ;

      scope cls_scope = wrp;
      {
        typedef typename cl::depositor_state cl;
        class_<cl>("DepositorState")
          .DEF_RW_MEMBER(stats)
          ;
      }
    }


    typedef advective_depositor<
      ParticleState, used_shape_function> adv_dep;
    {
      typedef adv_dep cl;
      class_<cl> wrp(
        ("AdvectiveDepositor"+get_state_class_suffix<ParticleState>()).c_str(), 
        init<const mesh_data &, unsigned, unsigned, 
        const py_matrix &,
        const py_matrix &,
        const py_matrix &,
        const py_matrix &,
        boost::shared_ptr<typename adv_dep::face_group_type>,
        boost::shared_ptr<typename adv_dep::face_group_type>,
        double, double, double
        >(args(
            "mesh_data", "faces_per_element", "dofs_per_element",
            "mass_matrix", "inverse_mass_matrix", "filter_matrix",
            "face_mass_matrix", "int_face_group", "bdry_face_group",
            "activation_threshold", "kill_threshold", "upwind_alpha"
            )));

      wrp
        .DEF_SIMPLE_METHOD(add_local_diff_matrix)
        .DEF_SIMPLE_METHOD(add_advective_particle)
        .DEF_SIMPLE_METHOD(get_debug_quantity_on_mesh)

        .DEF_SIMPLE_METHOD(get_advective_particle_rhs)
        .def("apply_advective_particle_rhs",
            &cl::apply_advective_particle_rhs,
            return_value_policy<manage_new_object>())
        
        .DEF_SIMPLE_METHOD(perform_depositor_upkeep)
        .DEF_SIMPLE_METHOD(kill_advected_particle)
        .DEF_SIMPLE_METHOD(note_move)
        .DEF_SIMPLE_METHOD(note_change_size)
        ;

      scope cls_scope = wrp;
      {
        typedef typename cl::depositor_state cl;
        class_<cl>("DepositorState")
          .DEF_RW_MEMBER(rho_dof_shift_listener)

          .DEF_RO_MEMBER(active_elements)
          .DEF_SIMPLE_METHOD(count_advective_particles)
          .DEF_RO_MEMBER(element_activation_counter)
          .DEF_RO_MEMBER(element_kill_counter)

          .DEF_SIMPLE_METHOD(resize_rho)
          .DEF_SIMPLE_METHOD(clear)
          ;
      }
    }

    EXPOSE_FOR_ALL_TARGET_RECONSTRUCTORS(
        expose_deposition_functions, 
        used_shape_function,
        ());

    expose_grid_depositor<ParticleState, brick>("Regular", gdbs_wrap);
    expose_grid_depositor<ParticleState, jiggly_brick>("Jiggly", gdbs_wrap);

    {
      typedef grid_find_depositor<ParticleState, used_shape_function, brick> cl;

      class_<cl> wrp(
        ("GridFindDepositor"+get_state_class_suffix<ParticleState>()).c_str(), 
        init<const mesh_data &>());

      wrp
        .DEF_RW_MEMBER(shape_function)
        .DEF_RW_MEMBER(bricks)
        .DEF_RW_MEMBER(node_number_list_starts)
        .DEF_RW_MEMBER(node_number_lists)

        .DEF_SIMPLE_METHOD(grid_node_count)

        .DEF_SIMPLE_METHOD(deposit_densities)
        .DEF_SIMPLE_METHOD(deposit_j)
        .DEF_SIMPLE_METHOD(deposit_rho)
        ;

      wrp.attr("DepositorState") = gdbs_wrap;
    }
  }
}




void expose_deposition()
{
  class_<grid_depositor_base_state> gdbs_wrap("GridDepositorBaseState");
  {
    typedef grid_depositor_base_state cl;
    gdbs_wrap
      .DEF_SIMPLE_METHOD(note_move)
      .DEF_SIMPLE_METHOD(note_change_size)
      ;
  }

  EXPOSE_FOR_ALL_STATE_TYPES(expose_depositors_for_pstate, (gdbs_wrap));

  {
    typedef element_on_grid cl;
    python::class_<cl>("ElementOnGrid")
      .DEF_RW_MEMBER(element_number)
      .DEF_RW_MEMBER(grid_nodes)
      .DEF_BYVAL_RW_MEMBER(weight_factors)
      .DEF_BYVAL_RW_MEMBER(interpolation_matrix)
      .DEF_BYVAL_RW_MEMBER(inverse_interpolation_matrix)
      ;
  }

  expose_std_vector<element_on_grid>("ElementOnGrid");

  python::def("get_shape_function_name", &used_shape_function::name);

  {
    typedef normalized_deposition_stats cl;
    class_<cl>("NormalizedDepositionStats")
      .DEF_RW_MEMBER(normalization_stats)
      .DEF_RW_MEMBER(centroid_distance_stats)
      .DEF_RW_MEMBER(el_per_particle_stats)
      ;
  }
}
