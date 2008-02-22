<TeXmacs|1.0.6.12>

<style|<tuple|article|maxima>>

<\body>
  <section|Hyperbolic Cleaning>

  The PDE solved in hyperbolic cleaning is as follows:

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>E-<frac|1|\<varepsilon\>>\<nabla\>\<times\>H+c<rsup|2>\<nabla\>\<Phi\>>|<cell|=>|<cell|-<frac|j|\<varepsilon\><rsub|0>>>>|<row|<cell|\<partial\><rsub|t>H+<frac|1|\<mu\>>\<nabla\>\<times\>E>|<cell|=>|<cell|0>>|<row|<cell|\<partial\><rsub|t>\<Phi\>+\<chi\><rsup|2>\<nabla\>\<cdot\>E>|<cell|=>|<cell|\<chi\><rsup|2><frac|\<rho\>|\<varepsilon\><rsub|0>>>>>>
  </eqnarray*>

  <with|prog-language|maxima|prog-session|default|<\session>
    <\input|<with|color|red|(<with|math-font-family|rm|%i>1)
    <with|color|black|>>>
      kill(all);
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o0>)
      <with|color|black|>><with|math-font-family|bf|done>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>1)
    <with|color|black|>>>
      load("eigen");
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o1>)
      <with|color|black|>><with|mode|text|/usr/share/maxima/5.13.0/share/matrix/eigen.mac>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>2)
    <with|color|black|>>>
      load("itensor");
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o2>)
      <with|color|black|>><with|mode|text|/usr/share/maxima/5.13.0/share/tensor/itensor.lisp>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>3)
    <with|color|black|>>>
      assume(c\<gtr\>0);assume(mu\<gtr\>0);assume(epsilon\<gtr\>0);
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o3>)
      <with|color|black|>><left|[>c\<gtr\>0<right|]>>

      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o4>)
      <with|color|black|>><left|[>\<mu\>\<gtr\>0<right|]>>

      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o5>)
      <with|color|black|>><left|[>\<varepsilon\>\<gtr\>0<right|]>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>6)
    <with|color|black|>>>
      /* ------------------------------------------- */

      /* Matrix creation helpers */

      /* ------------------------------------------- */

      \;

      /* A matrix resulting from a cross product */

      cpmat(coord):=genmatrix(

      \ \ lambda ([i,j], levi_civita([coord,i,j])),

      \ \ 3,3)$

      \;

      /* A constant matrix of size n x m */

      constmatrix(n,m,c):=genmatrix(lambda ([i,j], c),n,m)$

      \;

      vstack(a,b):=append(a,b)$

      hstack(a,b):=transpose(append(transpose(a),transpose(b)))$

      blockmat(a11,a12,a21,a22):=vstack(hstack(a11,a12),hstack(a21,a22))$
    </input>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>11)
    <with|color|black|>>>
      /* ------------------------------------------- */

      /* Begin treatment of Maxwell system */

      /* ------------------------------------------- */

      maxwellmat(i):=blockmat(

      \ \ zeromatrix(3,3),

      \ \ -epsinv*cpmat(i), /* epsinv = 1/epsilon */

      \ \ muinv*cpmat(i), /* muinv = 1/mu */

      \ \ zeromatrix(3,3))$

      \;

      n:[nx,ny,nz]$

      \;

      Amaxsimp:sum(n[i]*maxwellmat(i),i,1,3)$

      \;

      Amax:ev(Amaxsimp,epsinv=1/epsilon,muinv=1/mu)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o14>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|-<frac|<with|math-font-family|rm|nz>|\<varepsilon\>>>|<cell|<frac|<with|math-font-family|rm|ny>|\<varepsilon\>>>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|<frac|<with|math-font-family|rm|nz>|\<varepsilon\>>>|<cell|0>|<cell|-<frac|<with|math-font-family|rm|nx>|\<varepsilon\>>>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|-<frac|<with|math-font-family|rm|ny>|\<varepsilon\>>>|<cell|<frac|<with|math-font-family|rm|nx>|\<varepsilon\>>>|<cell|0>>|<row|<cell|0>|<cell|<frac|<with|math-font-family|rm|nz>|\<mu\>>>|<cell|-<frac|<with|math-font-family|rm|ny>|\<mu\>>>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|-<frac|<with|math-font-family|rm|nz>|\<mu\>>>|<cell|0>|<cell|<frac|<with|math-font-family|rm|nx>|\<mu\>>>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|<frac|<with|math-font-family|rm|ny>|\<mu\>>>|<cell|-<frac|<with|math-font-family|rm|nx>|\<mu\>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>15)
    <with|color|black|>>>
      evAmax:subst(

      \ \ [epsinv=1/epsilon, muinv=1/mu],

      \ \ ratsubst(1,n.n, eigenvectors(Amaxsimp))

      )$
    </input>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>16)
    <with|color|black|>>>
      Vmax:transpose(apply(matrix, makelist(evAmax[i],i,2,7)))$
    </input>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>17)
    <with|color|black|>>>
      invVmax:ratsubst(1,n.n,ratsimp(invert(Vmax)))$
    </input>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>18)
    <with|color|black|>>>
      ratsubst(1,n.n,Vmax.invVmax)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o18>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|1>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|1>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|1>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|1>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|1>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|1>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>19)
    <with|color|black|>>>
      ratsubst(1,n.n,invVmax.Amax.Vmax)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o19>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|-<frac|1|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|-<frac|1|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|<frac|1|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|<frac|1|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>20)
    <with|color|black|>>>
      /* ------------------------------------------- */

      /* Begin treatment of hyp. cleaning system */

      /* ------------------------------------------- */

      assume(chi\<gtr\>0);
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o20>)
      <with|color|black|>><left|[>\<chi\>\<gtr\>0<right|]>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>21)
    <with|color|black|>>>
      Acleansimp:blockmat(

      \ \ Amaxsimp,

      \ \ vstack(epsinv*muinv*covect(n),constmatrix(3,1,0)),

      \ \ hstack(chi^2*n,constmatrix(1,3,0)),

      \ \ zeromatrix(1,1)

      )$

      Aclean:subst([epsinv=1/epsilon, muinv=1/mu], Acleansimp)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o22>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|-<frac|<with|math-font-family|rm|nz>|\<varepsilon\>>>|<cell|<frac|<with|math-font-family|rm|ny>|\<varepsilon\>>>|<cell|<frac|<with|math-font-family|rm|nx>|\<varepsilon\>*\<mu\>>>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|<frac|<with|math-font-family|rm|nz>|\<varepsilon\>>>|<cell|0>|<cell|-<frac|<with|math-font-family|rm|nx>|\<varepsilon\>>>|<cell|<frac|<with|math-font-family|rm|ny>|\<varepsilon\>*\<mu\>>>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|-<frac|<with|math-font-family|rm|ny>|\<varepsilon\>>>|<cell|<frac|<with|math-font-family|rm|nx>|\<varepsilon\>>>|<cell|0>|<cell|<frac|<with|math-font-family|rm|nz>|\<varepsilon\>*\<mu\>>>>|<row|<cell|0>|<cell|<frac|<with|math-font-family|rm|nz>|\<mu\>>>|<cell|-<frac|<with|math-font-family|rm|ny>|\<mu\>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|-<frac|<with|math-font-family|rm|nz>|\<mu\>>>|<cell|0>|<cell|<frac|<with|math-font-family|rm|nx>|\<mu\>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|<frac|<with|math-font-family|rm|ny>|\<mu\>>>|<cell|-<frac|<with|math-font-family|rm|nx>|\<mu\>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|\<chi\><rsup|2>*<with|math-font-family|rm|nx>>|<cell|\<chi\><rsup|2>*<with|math-font-family|rm|ny>>|<cell|\<chi\><rsup|2>*<with|math-font-family|rm|nz>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>23)
    <with|color|black|>>>
      evAclean:subst(

      \ \ [epsinv=1/epsilon, muinv=1/mu],

      \ \ ratsubst(1,n.n,eigenvectors(Acleansimp))

      )$
    </input>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>24)
    <with|color|black|>>>
      Vclean:transpose(apply(matrix, makelist(evAclean[i],i,2,8)))$
    </input>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>25)
    <with|color|black|>>>
      invVclean:ratsubst(1,n.n,ratsimp(invert(Vclean)))$
    </input>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>26)
    <with|color|black|>>>
      ratsubst(1,n.n,Vclean.invVclean)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o26>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|1>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|1>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|1>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|1>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|1>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|1>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|1>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>27)
    <with|color|black|>>>
      Dclean:ratsubst(1,n.n,invVclean.Aclean.Vclean)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o27>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|-<frac|1|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|-<frac|1|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|<frac|1|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|<frac|1|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|-<frac|\<chi\>|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|<frac|\<chi\>|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>28)
    <with|color|black|>>>
      cleanwm:[E1m,E2m,E3m,B1m,B2m,B3m,phi];
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o28>)
      <with|color|black|>><left|[><with|math-font-family|rm|E1m>,<with|math-font-family|rm|E2m>,<with|math-font-family|rm|E3m>,<with|math-font-family|rm|B1m>,<with|math-font-family|rm|B2m>,<with|math-font-family|rm|B3m>,\<varphi\><right|]>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>29)
    <with|color|black|>>>
      cleansminw:invVclean.cleanwm$
    </input>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>30)
    <with|color|black|>>>
      /* ------------------------------------------- */

      /* radiation boundary condition for cleaning system */

      /* ------------------------------------------- */

      cleanradbdryspinw:makelist(

      \ \ if Dclean[i,i] \<gtr\>= 0 then cleansminw[i,1] else 0,\ 

      i, 1, 7)$

      \;

      ratsimp(

      ratsubst(nz^2,1-nx^2-ny^2,

      ratsubst(1,n.n,

      Vclean.cleanradbdryspinw)))
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o38>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|<frac|<with|math-font-family|rm|nx>*\<varphi\>-<with|math-font-family|rm|B2m>*\<chi\>*\<mu\>*<with|math-font-family|rm|nz>+<with|math-font-family|rm|B3m>*\<chi\>*\<mu\>*<with|math-font-family|rm|ny>+\<chi\>*<with|math-font-family|rm|E1m>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>|2*\<chi\>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>>|<row|<cell|<frac|<with|math-font-family|rm|ny>*\<varphi\>+<with|math-font-family|rm|B1m>*\<chi\>*\<mu\>*<with|math-font-family|rm|nz>-<with|math-font-family|rm|B3m>*\<chi\>*\<mu\>*<with|math-font-family|rm|nx>+\<chi\>*<with|math-font-family|rm|E2m>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>|2*\<chi\>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>>|<row|<cell|<frac|<with|math-font-family|rm|nz>*\<varphi\>-<with|math-font-family|rm|B1m>*\<chi\>*\<mu\>*<with|math-font-family|rm|ny>+<with|math-font-family|rm|B2m>*\<chi\>*\<mu\>*<with|math-font-family|rm|nx>+\<chi\>*<with|math-font-family|rm|E3m>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>|2*\<chi\>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>>|<row|<cell|<frac|<sqrt|\<mu\>>*<left|(><with|math-font-family|rm|B3m>*<with|math-font-family|rm|nx>*<with|math-font-family|rm|nz>+<with|math-font-family|rm|B2m>*<with|math-font-family|rm|nx>*<with|math-font-family|rm|ny>+<with|math-font-family|rm|B1m>*<with|math-font-family|rm|nx><rsup|2>+<with|math-font-family|rm|B1m><right|)>+<sqrt|\<varepsilon\>>*<left|(><with|math-font-family|rm|E2m>*<with|math-font-family|rm|nz>-<with|math-font-family|rm|E3m>*<with|math-font-family|rm|ny><right|)>|2*<sqrt|\<mu\>>>>>|<row|<cell|-<frac|<sqrt|\<mu\>>*<left|(><with|math-font-family|rm|B2m>*<with|math-font-family|rm|nz><rsup|2>-<with|math-font-family|rm|B3m>*<with|math-font-family|rm|ny>*<with|math-font-family|rm|nz>-<with|math-font-family|rm|B1m>*<with|math-font-family|rm|nx>*<with|math-font-family|rm|ny>+<with|math-font-family|rm|B2m>*<with|math-font-family|rm|nx><rsup|2>-2*<with|math-font-family|rm|B2m><right|)>+<sqrt|\<varepsilon\>>*<left|(><with|math-font-family|rm|E1m>*<with|math-font-family|rm|nz>-<with|math-font-family|rm|E3m>*<with|math-font-family|rm|nx><right|)>|2*<sqrt|\<mu\>>>>>|<row|<cell|<frac|<sqrt|\<mu\>>*<left|(><with|math-font-family|rm|B3m>*<with|math-font-family|rm|nz><rsup|2>+<left|(><with|math-font-family|rm|B2m>*<with|math-font-family|rm|ny>+<with|math-font-family|rm|B1m>*<with|math-font-family|rm|nx><right|)>*<with|math-font-family|rm|nz>+<with|math-font-family|rm|B3m><right|)>+<sqrt|\<varepsilon\>>*<left|(><with|math-font-family|rm|E1m>*<with|math-font-family|rm|ny>-<with|math-font-family|rm|E2m>*<with|math-font-family|rm|nx><right|)>|2*<sqrt|\<mu\>>>>>|<row|<cell|<frac|\<varphi\>+<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>*<left|(>\<chi\>*<with|math-font-family|rm|E3m>*<with|math-font-family|rm|nz>+\<chi\>*<with|math-font-family|rm|E2m>*<with|math-font-family|rm|ny>+\<chi\>*<with|math-font-family|rm|E1m>*<with|math-font-family|rm|nx><right|)>|2>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>35)
    <with|color|black|>>>
      /* ------------------------------------------- */

      /* radiation boundary only for "chi" waves */

      /* ------------------------------------------- */

      cleanchiradbdryspinw:makelist(

      \ \ if Dclean[i,i] \<gtr\>= 0 and not
      (ratsimp(diff(Dclean[i,i],chi))=0) then\ 

      \ \ \ \ cleansminw[i,1]\ 

      \ \ else\ 

      \ \ \ \ 0,\ 

      i, 1, 7)$

      \;

      ratsimp(

      ratsubst(nz^2,1-nx^2-ny^2,

      ratsubst(1,n.n,

      Vclean.cleanchiradbdryspinw)))
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o55>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|<frac|<with|math-font-family|rm|nx>*\<varphi\>+<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>*<left|(>\<chi\>*<with|math-font-family|rm|E3m>*<with|math-font-family|rm|nx>*<with|math-font-family|rm|nz>+\<chi\>*<with|math-font-family|rm|E2m>*<with|math-font-family|rm|nx>*<with|math-font-family|rm|ny>+\<chi\>*<with|math-font-family|rm|E1m>*<with|math-font-family|rm|nx><rsup|2><right|)>|2*\<chi\>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>>|<row|<cell|-<frac|<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>*<left|(>\<chi\>*<with|math-font-family|rm|E2m>*<with|math-font-family|rm|nz><rsup|2>-\<chi\>*<with|math-font-family|rm|E3m>*<with|math-font-family|rm|ny>*<with|math-font-family|rm|nz>-\<chi\>*<with|math-font-family|rm|E1m>*<with|math-font-family|rm|nx>*<with|math-font-family|rm|ny>+\<chi\>*<with|math-font-family|rm|E2m>*<with|math-font-family|rm|nx><rsup|2>-\<chi\>*<with|math-font-family|rm|E2m><right|)>-<with|math-font-family|rm|ny>*\<varphi\>|2*\<chi\>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>>|<row|<cell|<frac|<with|math-font-family|rm|nz>*\<varphi\>+<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>*<left|(>\<chi\>*<with|math-font-family|rm|E3m>*<with|math-font-family|rm|nz><rsup|2>+<left|(>\<chi\>*<with|math-font-family|rm|E2m>*<with|math-font-family|rm|ny>+\<chi\>*<with|math-font-family|rm|E1m>*<with|math-font-family|rm|nx><right|)>*<with|math-font-family|rm|nz><right|)>|2*\<chi\>*<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>>>>|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|<frac|\<varphi\>+<sqrt|\<varepsilon\>>*<sqrt|\<mu\>>*<left|(>\<chi\>*<with|math-font-family|rm|E3m>*<with|math-font-family|rm|nz>+\<chi\>*<with|math-font-family|rm|E2m>*<with|math-font-family|rm|ny>+\<chi\>*<with|math-font-family|rm|E1m>*<with|math-font-family|rm|nx><right|)>|2>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>50)
    <with|color|black|>>>
      \;
    </input>
  </session>>
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Hyperbolic
      Cleaning> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>