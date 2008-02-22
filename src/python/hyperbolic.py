"""Hyperbolic cleaning operator."""

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




class CleaningMaxwellOperator(object):
    def __init__(self, maxwell_op, chi=1):
        self.discr = maxwell_op.discr
        self.maxwell_op = maxwell_op
        self.chi = chi

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

        addl_flux = 0.5*join_fields(
                # flux e
                (normal*(
                    dot(normal, e.int-e.ext) 
                    - maxwell_op.c**2*(phi.int - phi.ext))),
                # flux h
                [0] * h_components,
                # flux phi
                chi**2*(dot(normal, e.int-e.ext) - (phi.int - phi.ext)),
                )
        self.addl_flux = self.discr.get_flux_operator(addl_flux)

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
        from pytools.arithmetic_container import join_fields, ArithmeticList
        
        pec_tag = self.maxwell_op.pec_tag
        nabla = self.maxwell_op.nabla

        e, h, phi = self.split_ehphi(w)

        pec_bc = join_fields(
                -self.discr.boundarize_volume_field(e, pec_tag),
                self.discr.boundarize_volume_field(h, pec_tag),
                self.discr.boundarize_volume_field(phi, pec_tag)
                )

        rhs_e, rhs_h = self.split_eh(self.maxwell_op.rhs(t, w))

        phi_decay = 0
        return (join_fields(
                rhs_e - self.maxwell_op.c**2*(nabla*cache_diff_results(phi)),
                rhs_h,
                self.chi**2*(rho - dot(nabla, cache_diff_results(e))) - phi_decay*phi
                )
                +
                self.maxwell_op.m_inv*(
                    self.addl_flux * w
                    +self.addl_flux * pair_with_boundary(w, pec_bc, pec_tag)
                    )
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

