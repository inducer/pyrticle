"""Hyperbolic cleaning operator."""

from __future__ import division 

__copyright__ = "Copyright (C) 2008 Andreas Kloeckner"

__license__ = """
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see U{http://www.gnu.org/licenses/}.
"""




import pylinear.array as num




class CleaningMaxwellOperator(object):
    pass




class ECleaningMaxwellOperator(CleaningMaxwellOperator):
    def __init__(self, maxwell_op, chi=1, phi_decay=0):
        self.discr = maxwell_op.discr
        self.maxwell_op = maxwell_op
        self.chi = chi
        self.phi_decay = phi_decay

        assert chi > 0
        assert phi_decay >= 0

        from hedge.flux import make_normal, FluxVectorPlaceholder

        normal = make_normal(self.discr.dimensions)

        e_components = maxwell_op.count_subset(maxwell_op.get_eh_subset()[0:3])
        h_components = maxwell_op.count_subset(maxwell_op.get_eh_subset()[3:6])

        w = FluxVectorPlaceholder(
                maxwell_op.count_subset(maxwell_op.get_eh_subset())+1)
        e, h = self.split_eh(w)
        phi = w[e_components+h_components]

        from pytools.arithmetic_container import join_fields
        from hedge.tools import dot

        # see hedge/doc/maxima/eclean.mac for derivation
        strong_flux = 0.5*maxwell_op.c*chi*join_fields(
                # flux e
                normal*(phi.int-phi.ext - dot(normal, e.int-e.ext)),
                # flux h
                len(h)*[0],
                # flux phi
                dot(e.int-e.ext, normal)-(phi.int-phi.ext)
                )
        strong_flux += join_fields(maxwell_op.flux, 0)

        self.strong_flux_op = self.discr.get_flux_operator(strong_flux)

        self.pec_normals = self.discr.boundary_normals(self.maxwell_op.pec_tag)

    @property
    def count_subset(self): return self.maxwell_op.count_subset
    @property
    def epsilon(self): return self.maxwell_op.epsilon
    @property
    def mu(self): return self.maxwell_op.mu
    @property
    def c(self): return self.maxwell_op.c
    @property
    def e_cross(self): return self.maxwell_op.e_cross
    @property
    def h_cross(self): return self.maxwell_op.h_cross

    def get_eh_subset(self):
        return self.maxwell_op.get_eh_subset()

    def rhs(self, t, w, rho):
        from hedge.tools import dot
        from hedge.discretization import pair_with_boundary, cache_diff_results
        from pytools.arithmetic_container import join_fields
        
        pec_tag = self.maxwell_op.pec_tag
        nabla = self.maxwell_op.nabla

        c = self.maxwell_op.c
        chi = self.chi

        e, h, phi = self.split_ehphi(w)

        pec_e = self.discr.boundarize_volume_field(e, self.maxwell_op.pec_tag)
        pec_h = self.discr.boundarize_volume_field(h, self.maxwell_op.pec_tag)
        pec_phi = self.discr.boundarize_volume_field(phi, self.maxwell_op.pec_tag)
        pec_n = self.pec_normals

        from pytools.arithmetic_container import work_with_arithmetic_containers
        ac_multiply = work_with_arithmetic_containers(num.multiply)

        # see hedge/doc/maxima/eclean.mac for derivation
        pec_bc = join_fields(
                -pec_e
                +2*ac_multiply(pec_n, dot(pec_n, pec_e, num.multiply)),
                pec_h,
                -pec_phi)

        e_cache = cache_diff_results(e)
        h_cache = cache_diff_results(h)
        phi_cache = cache_diff_results(phi)

        # in conservation form: u_t + A u_x = 0
        max_local_op = join_fields(self.maxwell_op.local_op(e_cache, h_cache), 0)
        local_operator = max_local_op + join_fields(
                c*chi*(nabla*phi_cache),
                0*h,
                c*chi*(dot(nabla, e_cache) - rho/self.maxwell_op.epsilon) 
                - self.phi_decay * phi
                )

        return -local_operator + self.maxwell_op.m_inv*(
                    self.strong_flux_op * w
                    +self.strong_flux_op * pair_with_boundary(w, pec_bc, pec_tag)
                    )

    def assemble_fields(self, e=None, h=None, phi=None):
        if phi is None:
            phi = self.discr.volume_zeros()

        from pytools.arithmetic_container import join_fields
        return join_fields(
                self.maxwell_op.assemble_fields(e, h),
                phi)

    def split_eh(self, w):
        return self.maxwell_op.split_eh(w)

    def split_ehphi(self, w):
        e, h = self.split_eh(w)
                
        eh_components = self.maxwell_op.count_subset(self.maxwell_op.get_eh_subset())
        phi = w[eh_components]
        return e, h, phi

    def max_eigenvalue(self):
        return self.chi*self.maxwell_op.max_eigenvalue()
