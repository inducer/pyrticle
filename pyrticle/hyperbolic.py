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
from hedge.models import HyperbolicOperator




class CleaningMaxwellOperator(HyperbolicOperator):
    pass




class ECleaningMaxwellOperator(CleaningMaxwellOperator):
    def __init__(self, maxwell_op, chi=1, phi_decay=0):
        self.maxwell_op = maxwell_op
        self.chi = chi
        self.phi_decay = phi_decay

        assert chi > 0
        assert phi_decay >= 0

        from hedge.tools import count_subset
        self.component_count = count_subset(maxwell_op.get_eh_subset())+1

    @property
    def dimensions(self): return self.maxwell_op.dimensions
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

    def flux(self):
        from hedge.flux import make_normal, FluxVectorPlaceholder

        normal = make_normal(self.maxwell_op.dimensions)

        from hedge.tools import join_fields

        w = FluxVectorPlaceholder(self.component_count)
        e, h, phi = self.split_ehphi(w)

        # see hedge/doc/maxima/eclean.mac for derivation
        strong_flux = 0.5*self.c*self.chi*join_fields(
                # flux e
                normal*(phi.int-phi.ext - numpy.dot(normal, e.int-e.ext)),
                # flux h
                len(h)*[0],
                # flux phi
                numpy.dot(e.int-e.ext, normal)-(phi.int-phi.ext)
                )
        return strong_flux + join_fields(self.maxwell_op.flux(1), 0)

    @memoize_method
    def op_template(self):
        from hedge.tools import join_fields
        from hedge.optemplate import Field, make_vector_field, BoundaryPair, \
                BoundarizeOperator, make_normal, get_flux_operator, \
                make_nabla, InverseMassOperator

        w = make_vector_field("w", self.component_count)
        e, h, phi = self.split_ehphi(w)
        rho = Field("rho")


        # local part ----------------------------------------------------------
        nabla = make_nabla(self.maxwell_op.dimensions)

        # in conservation form: u_t + A u_x = 0
        # we're describing the A u_x part, the sign gets reversed
        # below.
        max_local_op = join_fields(self.maxwell_op.local_derivatives(w), 0)

        c = self.maxwell_op.c
        chi = self.chi

        hyp_local_operator = max_local_op + join_fields(
                c*chi*(nabla*phi),
                0*h,
                c*chi*numpy.dot(nabla, e) 
                - c*chi*rho/self.maxwell_op.epsilon

                # sign gets reversed below, so this is actually 
                # the decay it advertises to be.
                + self.phi_decay*phi
                )

        # BCs -----------------------------------------------------------------
        pec_tag = self.maxwell_op.pec_tag

        pec_e = BoundarizeOperator(pec_tag)(e)
        pec_h = BoundarizeOperator(pec_tag)(h)
        pec_phi = BoundarizeOperator(pec_tag)(phi)
        pec_n = make_normal(pec_tag, self.maxwell_op.dimensions)

        from hedge.tools import log_shape

        from hedge.tools import ptwise_dot
        bc = "invent"
        print "HYP CLEAN BC", bc
        if bc == "char":
            # see hedge/doc/maxima/eclean.mac for derivation
            pec_bc = join_fields(
                    -pec_e
                    + 3/2 * pec_n * numpy.dot(pec_n, pec_e)
                    + 1/2 * pec_phi * pec_n ,

                    pec_h,

                    1/2*(pec_phi+numpy.dot(pec_n, pec_e))
                    )
        if bc == "invent":
            # see hedge/doc/maxima/eclean.mac for derivation
            pec_bc = join_fields(
                    -pec_e
                    + 2 * pec_n * numpy.dot(pec_n, pec_e),
                    pec_h, 
                    pec_phi)
        elif bc == "munz":
            # Munz et al
            pec_bc = join_fields(
                    -pec_e, 
                    pec_h, 
                    pec_phi-numpy.dot(pec_n, pec_e))
        elif bc == "prev":
            # previous condition
            pec_bc = join_fields(
                    -pec_e
                    +2*pec_n * numpy.dot(pec_n, pec_e),

                    pec_h,

                    -pec_phi)

        # assemble operator ---------------------------------------------------
        flux_op = get_flux_operator(self.flux())
        return -hyp_local_operator + InverseMassOperator() * (
                    flux_op * w
                    + flux_op * BoundaryPair(w, pec_bc, pec_tag)
                    )

    def bind(self, discr):
        op = discr.compile(self.op_template())

        def rhs(t, w, rho):
            return op(w=w, rho=rho)

        return rhs

    def assemble_fields(self, e=None, h=None, phi=None, discr=None):
        if discr is None:
            def zero(): return 0
        else:
            def zero(): return discr.volume_zeros()

        if phi is None:
            phi = zero()

        from hedge.tools import join_fields
        return join_fields(
                self.maxwell_op.assemble_fields(e, h),
                phi)

    def split_eh(self, w):
        return self.maxwell_op.split_eh(w)

    def split_ehphi(self, w):
        e, h = self.split_eh(w)

        from hedge.tools import count_subset
        eh_components = count_subset(self.maxwell_op.get_eh_subset())
        phi = w[eh_components]
        return e, h, phi

    def max_eigenvalue(self, t, fields=None, discr=None):
        return self.chi*self.maxwell_op.max_eigenvalue(t, fields, discr)




class PhiFilter:
    def __init__(self, maxwell_op, filter):
        self.maxwell_op = maxwell_op
        self.filter = filter

    def __call__(self, em_fields):
        e, h, phi = self.maxwell_op.split_ehphi(em_fields)
        return self.maxwell_op.assemble_fields(e, h, self.filter(phi))
