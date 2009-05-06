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




import numpy
from pytools import memoize_method




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

        self.component_count = maxwell_op.count_subset(maxwell_op.get_eh_subset())+1
        w = FluxVectorPlaceholder(self.component_count)
        e, h, phi = self.split_ehphi(w)

        from hedge.tools import join_fields

        # see hedge/doc/maxima/eclean.mac for derivation
        strong_flux = 0.5*maxwell_op.c*chi*join_fields(
                # flux e
                normal*(phi.int-phi.ext - numpy.dot(normal, e.int-e.ext)),
                # flux h
                len(h)*[0],
                # flux phi
                numpy.dot(e.int-e.ext, normal)-(phi.int-phi.ext)
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

    @memoize_method
    def op_template(self):
        from hedge.tools import join_fields
        from hedge.optemplate import make_vector_field, pair_with_boundary

        w = make_vector_field("w", self.component_count)
        e, h, phi = self.split_ehphi(w)
        pec_bc = make_vector_field("pec_bc", self.component_count)

        nabla = self.discr.nabla
        m_inv = self.discr.inverse_mass_operator

        # in conservation form: u_t + A u_x = 0
        max_local_op = join_fields(self.maxwell_op.local_op(e, h), 0)

        c = self.maxwell_op.c
        chi = self.chi

        local_operator = max_local_op + join_fields(
                c*chi*(nabla*phi),
                0*h,
                c*chi*numpy.dot(nabla, e)
                )

        pec_tag = self.maxwell_op.pec_tag

        return self.discr.compile(
                -local_operator + m_inv*(
                    self.strong_flux_op * w
                    +self.strong_flux_op * pair_with_boundary(w, pec_bc, pec_tag)
                    )
                )

    def rhs(self, t, w, rho):
        from hedge.tools import join_fields

        e, h, phi = self.split_ehphi(w)

        pec_tag = self.maxwell_op.pec_tag
        pec_e = self.discr.boundarize_volume_field(e, self.maxwell_op.pec_tag)
        pec_h = self.discr.boundarize_volume_field(h, self.maxwell_op.pec_tag)
        pec_phi = self.discr.boundarize_volume_field(phi, self.maxwell_op.pec_tag)
        pec_n = self.pec_normals

        from hedge.tools import log_shape

        from hedge.tools import ptwise_dot
        # see hedge/doc/maxima/eclean.mac for derivation
        pec_bc = join_fields(
                -pec_e
                +2*pec_n * ptwise_dot(1, 1, pec_n, pec_e),
                pec_h,
                -pec_phi)

        c = self.maxwell_op.c
        chi = self.chi
        return (
                self.op_template()(w=w, pec_bc=pec_bc) 
                +
                self.assemble_fields(
                    phi=c*chi*rho/self.maxwell_op.epsilon - self.phi_decay * phi)
                )

    def assemble_fields(self, e=None, h=None, phi=None):
        if phi is None:
            phi = self.discr.volume_zeros()

        from hedge.tools import join_fields
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




class PhiFilter:
    def __init__(self, maxwell_op, filter):
        self.maxwell_op = maxwell_op
        self.filter = filter

    def __call__(self, em_fields):
        e, h, phi = self.maxwell_op.split_ehphi(em_fields)
        return self.maxwell_op.assemble_fields(e, h, self.filter(phi))
