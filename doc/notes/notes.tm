<TeXmacs|1.0.6.14>

<style|<tuple|article|axiom|maxima>>

<\body>
  <doc-data|<doc-title|Notes on PIC>>

  <section|Shape function Integrals>

  <with|prog-language|maxima|prog-session|default|<\session>
    <\input|<with|color|red|(<with|math-font-family|rm|%i>7)
    <with|color|black|>>>
      omega(n) := 2*%pi^(n/2)/gamma(n/2)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o1>)
      <with|color|black|>>\<omega\><left|(>n<right|)>:=<frac|2*\<pi\><rsup|<frac|n|2>>|\<Gamma\><left|(><frac|n|2><right|)>>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>2)
    <with|color|black|>>>
      omega(10)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o11>)
      <with|color|black|>><frac|\<pi\><rsup|5>|12>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>12)
    <with|color|black|>>>
      intval:integrate((l-r^2/l)^alpha*spherearea*r^(n-1),r,0,l)
    </input>

    <\input|<with|color|red|><with|mode|math|math-display|true|<with|mode|text|Is
    >l<with|mode|text| positive or negative?>><with|color|black|>>
      positive
    </input>

    <\input|<with|color|red|><with|mode|math|math-display|true|<with|mode|text|Is
    >n<with|mode|text| positive, negative, or zero?>><with|color|black|>>
      positive
    </input>

    <\input|<with|color|red|><with|mode|math|math-display|true|<with|mode|text|Is
    >\<alpha\>+1<with|mode|text| positive, negative, or
    zero?>><with|color|black|>>
      positive
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o12>)
      <with|color|black|>><frac|l<rsup|n+\<alpha\>>*\<beta\><left|(><frac|n|2>,\<alpha\>+1<right|)>*<with|math-font-family|rm|spherearea>|2>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>13)
    <with|color|black|>>>
      ev(intval, n=2, alpha=2)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o14>)
      <with|color|black|>><frac|l<rsup|4>*<with|math-font-family|rm|spherearea>|6>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>15)
    <with|color|black|>>>
      \;
    </input>
  </session>>

  <section|Initial Conditions by Lorentz Transform>

  Before we start, a brief remark on notation: If we let <math|\<b-v\>> be a
  vector and <math|\<b-beta\>> another vector, then we define by

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-v\><rsub|\<\|\|\>\<b-beta\>>>|<cell|\<assign\>>|<cell|<frac|1|\<b-beta\>\<cdot\>\<b-beta\>>\<b-beta\>(\<b-beta\>\<cdot\>\<b-v\>),>>|<row|<cell|\<b-v\><rsub|\<perp\>\<b-beta\>>>|<cell|\<assign\>>|<cell|\<b-v\>-\<b-v\><rsub|\<\|\|\>\<b-beta\>>>>>>
  </eqnarray*>

  the components of <math|\<b-v\>> parallel and orthogonal to
  <math|\<b-beta\>>. Note that

  <\enumerate-alpha>
    <item>the length of <math|\<b-beta\>> plays no role in this. If
    <math|\<b-beta\>=\<b-0\>>, we define <math|\<b-v\><rsub|\<\|\|\>\<b-beta\>>\<assign\>\<b-0\>>.

    <item><math|\<b-v\><rsub|\<\|\|\>\<b-beta\>>> and
    <math|\<b-v\><rsub|\<perp\>\<b-beta\>>> remain vectors.

    <item><math|\<b-v\><rsub|\<\|\|\>-\<b-beta\>>=\<b-v\><rsub|\<\|\|\>\<b-beta\>>>,
    i.e. the sign of <math|\<b-beta\>> also plays no role. Think of the
    operation as a projection onto <math|span{\<b-beta\>}>.
  </enumerate-alpha>

  The situation is a charge density <math|\<rho\>> moving with velocity
  <math|\<b-v\><rsub|0>>. We fix that <math|K> is the rest frame of that
  charge carrier, thus <math|\<b-v\>=-\<b-v\><rsub|0>>, and
  <math|K<rprime|'>> is the laboratory frame.

  According to Jackson (11.19), the coordinate Lorentz Transform from a frame
  <math|K> to a frame <math|K<rprime|'>> moving with velocity <math|\<b-v\>>
  relative to <math|K> is:

  <\eqnarray*>
    <tformat|<table|<row|<cell|x<rsub|0><rprime|'>>|<cell|=>|<cell|\<gamma\>(x<rsub|0>-\<b-beta\>\<cdot\>\<b-x\>),>>|<row|<cell|\<b-x\><rprime|'>>|<cell|=>|<cell|\<b-x\><rsub|\<perp\>\<b-beta\>>+\<gamma\>\<b-x\><rsub|\<\|\|\>\<b-beta\>>-\<gamma\>\<b-beta\>x<rsub|0>,>>>>
  </eqnarray*>

  where <math|x<rsub|0>=c*t>. Inverting the coordinate formulas, we also
  obtain

  <\eqnarray*>
    <tformat|<table|<row|<cell|x<rsub|0>>|<cell|=>|<cell|\<gamma\>(x<rsub|0><rprime|'>+\<b-beta\>\<cdot\>\<b-x\><rprime|'>),>>|<row|<cell|\<b-x\>>|<cell|=>|<cell|\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>-\<gamma\>\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>+\<gamma\>\<b-beta\>x<rsub|0><rprime|'>.>>>>
  </eqnarray*>

  We only care about time <math|t<rprime|'>=0> in <math|K>, and we assume
  that the situation in <math|K> is in steady-state, i.e.
  <math|f(\<b-x\>,x<rsub|0>)=f(\<b-x\>,0)>. Therefore, wishing to transport a
  function from <math|K> to <math|K<rprime|'>>, we find

  <\eqnarray*>
    <tformat|<table|<row|<cell|f<rprime|'>(\<b-x\><rprime|'>,0)>|<cell|=>|<cell|f(\<b-x\>,x<rsub|0>)>>|<row|<cell|>|<cell|=>|<cell|f(\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>-\<gamma\>\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>,x<rsub|0>)>>|<row|<cell|>|<cell|=>|<cell|f(\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>-\<gamma\>\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>,0).>>>>
  </eqnarray*>

  On the other hand,

  <\eqnarray*>
    <tformat|<table|<row|<cell|f(\<b-x\>,0)>|<cell|=>|<cell|f(\<b-x\>,\<b-beta\>\<cdot\>\<b-x\>)>>|<row|<cell|>|<cell|=>|<cell|f<rprime|'>(\<b-x\><rprime|'>,[\<b-beta\>\<cdot\>\<b-x\>]<rsub|x<rsub|0>\<rightarrow\>x<rsub|0><rprime|'>>)>>|<row|<cell|>|<cell|=>|<cell|f<rprime|'>(\<b-x\><rsub|\<perp\>\<b-beta\>>+\<gamma\>\<b-x\><rsub|\<\|\|\>\<b-beta\>>-\<gamma\>\<b-beta\>(\<b-beta\>\<cdot\>\<b-x\>),0)>>|<row|<cell|>|<cell|=>|<cell|f<rprime|'>(\<b-x\><rsub|\<perp\>\<b-beta\>>+\<gamma\>(1-\<beta\><rsup|2>)\<b-x\><rsub|\<\|\|\>\<b-beta\>>,0)>>|<row|<cell|>|<cell|=>|<cell|f<rprime|'>(\<b-x\><rsub|\<perp\>\<b-beta\>>+<frac|\<b-x\><rsub|\<\|\|\>\<b-beta\>>|\<gamma\>>,0).>>>>
  </eqnarray*>

  We may therefore neglect time dependency in the following.

  Next, let us worry about how charge densities are affected by the Lorentz
  transform. <math|c\<rho\>> is part of the 4-vector
  <math|J<rsup|\<alpha\>>=(c\<rho\>,\<b-J\>)>, therefore transforms in
  timelike fashion, and we obtain

  <\equation*>
    \<rho\><rprime|'>(\<b-x\><rprime|'>)=\<gamma\>\<rho\>(\<b-x\>).
  </equation*>

  This makes sense simply because the same amount of charge is distributed
  over a larger amount of space. We can find
  <math|\<rho\><rprime|'>(\<b-x\><rprime|'>)> in <math|K<rprime|'>> by shape
  function reconstruction. Assume we have \ a shape function <math|S>
  satisfying

  <\equation*>
    <big|int><rsub|\<bbb-R\><rsup|d>>S<rprime|'>(\<b-x\><rprime|'>)\<mathd\>\<b-x\><rprime|'>=1.
  </equation*>

  Then

  <\equation*>
    \<rho\><rprime|'>(\<b-x\><rprime|'>)=<big|sum><rsub|i>q<rsub|i>S<left|(><frac|\<b-x\><rprime|'>-\<b-x\><rprime|'><rsub|i>|r<rsub|i>><right|)>
  </equation*>

  satisfies the defining relationship

  <\equation*>
    Q=<big|int><rsub|\<bbb-R\><rsup|d>>\<rho\><rprime|'>(\<b-x\><rprime|'>)*\<mathd\>\<b-x\><rprime|'>=<big|sum><rsub|i>q<rsub|i><big|int><rsub|\<bbb-R\><rsup|d>>S<left|(><frac|\<b-x\><rprime|'>-\<b-x\><rprime|'><rsub|i>|r<rsub|i>><right|)>=<big|sum><rsub|i>q<rsub|i>.
  </equation*>

  It is important to avoid confusion here--since the positions of the
  particles are known in terms of lab frame coordinates, we necessarily
  reconstruct <math|\<rho\><rprime|'>>, not <math|\<rho\>>. <math|\<rho\>>
  can then be found as

  <\equation*>
    \<rho\>(\<b-x\>)=<frac|\<rho\><rprime|'>(\<b-x\><rprime|'>)|\<gamma\>>.
  </equation*>

  Once we obtain <math|\<rho\>> in the rest frame <math|K>, where we may
  assume the requisite radial symmetry of the particle, the initial
  <math|\<b-E\>> can be found as <math|\<nabla\>\<Phi\>> from
  <math|\<Delta\>\<Phi\>(\<b-x\>)=\<rho\>(\<b-x\>)/\<varepsilon\>>. Figure
  <reference|fig:dilation-lab-flat> illustrates the situation.

  <big-figure|<postscript|pic-dilation-rest-round.fig|11cm|||||>|<label|fig:dilation-lab-flat>Illustration
  of the effects of Lorentz transformation on initial conditions.>

  <with|color|red|An important question arises here:> The scheme used in the
  code right now does not do what Figure <reference|fig:dilation-lab-flat>
  says. Instead, we deposit particle shapes as round blobs in the lab frame,
  and so they appear strangely elongated in the rest frame, as shown in
  Figure <reference|fig:dilation>. This does not seem to be the correct thing
  to do--if the ``electrons'' have proper rotational symmetry anywhere, then
  it should be in their rest frame. The treatment of Figure
  <reference|fig:dilation-lab-flat> seems more appropriate, where charge
  blobs shortened by <math|\<gamma\>> are deposited in the lab frame,
  resulting in the ``correct'' shape in the rest frame. However, option
  <reference|fig:dilation-lab-flat> obviously makes charge deposition in the
  lab frame harder, because it necessitates more mesh resolution in the
  <math|z> direction by a factor of <math|\<gamma\>>. Can this be helped
  somehow? Can we get away with option <reference|fig:dilation> and not worry
  about it? Am I correct in asserting that Picture
  <reference|fig:dilation-lab-flat> is more physically correct?

  <\big-figure>
    \ <postscript|pic-dilation-2.fig|11cm|||||>
  </big-figure|<label|fig:dilation>Present charge deposition scheme. Produces
  elongated particles in the rest frame.>

  \;

  Note that we never want to explicitly deal with a function defined on a
  mesh in the rest frame, since we do not have (and do not want to construct)
  such a mesh. Instead, we treat functions living in the rest frame on the
  (existing) lab frame mesh, and simply adapt any rest frame derivative
  operators to act on lab frame functions as if the function was ``properly''
  in the rest frame. To this end, we define

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|\<b-E\>|~>(\<b-x\><rprime|'>)>|<cell|\<assign\>>|<cell|\<b-E\>(\<b-x\>)=\<b-E\><left|(>\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>-\<gamma\>\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'><right|)>,>>|<row|<cell|<wide|\<Phi\>|~>(\<b-x\><rprime|'>)>|<cell|\<assign\>>|<cell|\<Phi\>(\<b-x\>)=\<Phi\>(\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>-\<gamma\>\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>),>>|<row|<cell|<wide|\<rho\>|~>(\<b-x\><rprime|'>)>|<cell|\<assign\>>|<cell|\<rho\>(\<b-x\>)=\<rho\>(\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>-\<gamma\>\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>),>>>>
  </eqnarray*>

  i.e. the tilde'd (semi-rest frame) quantities have their
  <math|\<b-x\>>-dependencies switchted to the <math|K<rprime|'>> frame, but
  still carry the values they do in <math|K>, i.e. no Lorentz transform has
  taken place.

  Therefore

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\><rsub|\<b-x\>>\<Phi\>(\<b-x\>)>|<cell|=>|<cell|\<nabla\><rsub|\<b-x\>>[<wide|\<Phi\>|~>(\<b-x\><rprime|'>(\<b-x\>))]>>|<row|<cell|>|<cell|=>|<cell|(\<nabla\><wide|\<Phi\>|~>)(\<b-x\><rprime|'>)\<cdot\>(\<nabla\>\<b-x\><rprime|'>)(\<b-x\>)>>|<row|<cell|>|<cell|=>|<cell|(\<nabla\><wide|\<Phi\>|~>)(\<b-x\><rprime|'>)\<cdot\><left|[>Id<rsub|\<perp\>\<b-beta\>>+<frac|Id<rsub|\<\|\|\>\<b-beta\>>|\<gamma\>><right|]>>>|<row|<cell|>|<cell|=>|<cell|<left|(>\<nabla\><rsub|\<b-x\><rprime|'><rsub|\<perp\>\<b-beta\>>>+<frac|\<nabla\><rsub|\<b-x\><rprime|'><rsub|\<\|\|\>\<b-beta\>>>|\<gamma\>><right|)><wide|\<Phi\>|~>(\<b-x\><rprime|'>).>>>>
  </eqnarray*>

  We obtain

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Delta\><rsub|\<b-x\>>\<Phi\>(\<b-x\>)>|<cell|=>|<cell|<frac|\<rho\>(\<b-x\>)|\<varepsilon\>>>>|<row|<cell|\<\|\|\>>|<cell|>|<cell|<space|1em>\<\|\|\>>>|<row|<cell|<left|(>\<Delta\><rsub|\<b-x\><rprime|'><rsub|\<perp\>\<b-beta\>>>+<frac|\<Delta\><rsub|\<b-x\><rprime|'><rsub|\<\|\|\>\<b-beta\>>>|\<gamma\><rsup|2>><right|)><wide|\<Phi\>|~>(\<b-x\><rprime|'>)>|<cell|=>|<cell|<frac|<wide|\<rho\>|~>(\<b-x\><rprime|'>)|\<varepsilon\>>,>>>>
  </eqnarray*>

  and we can compute <math|<wide|\<Phi\>|~>(\<b-x\><rprime|'>)>. Next, we
  have

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-E\>(\<b-x\>)>|<cell|=>|<cell|\<nabla\><rsub|\<b-x\>>\<Phi\>(\<b-x\>)>>|<row|<cell|\<\|\|\>>|<cell|>|<cell|\<\|\|\>>>|<row|<cell|<wide|\<b-E\>|~>(\<b-x\><rprime|'>)>|<cell|=>|<cell|<left|(>\<nabla\><rsub|\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>>+<frac|\<nabla\><rsub|\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>>|\<gamma\>><right|)><wide|\<Phi\>|~>(\<b-x\><rprime|'>).>>>>
  </eqnarray*>

  According to Jackson (11.149), the Lorentz Transform for transforming
  <math|\<b-E\>> and <math|\<b-B\>> in Gaussian units, reads:

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-E\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\>(\<b-E\>*(\<b-x\>,t)+\<b-beta\>\<times\>\<b-B\>(\<b-x\>,t))-<frac|\<gamma\><rsup|2>|\<gamma\>+1>\<b-beta\>(\<b-beta\>\<cdot\>\<b-E\>(\<b-x\>,t)),>>|<row|<cell|\<b-B\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\>(\<b-B\>(\<b-x\>,t)-\<b-beta\>\<times\>\<b-E\>(\<b-x\>,t))-<frac|\<gamma\><rsup|2>|\<gamma\>+1>\<b-beta\>(\<b-beta\>\<cdot\>\<b-B\>(x,t)).>>>>
  </eqnarray*>

  Using the identity

  <\equation*>
    \<gamma\>-<frac|\<gamma\><rsup|2>\<beta\><rsup|2>|\<gamma\>+1>=<frac|\<gamma\><rsup|2>+\<gamma\>-\<gamma\><rsup|2>\<beta\><rsup|2>|\<gamma\>+1>=<frac|\<gamma\><rsup|2>(1-\<beta\><rsup|2>)+\<gamma\>|\<gamma\>+1>=<frac|1+\<gamma\>|\<gamma\>+1>=1,
  </equation*>

  we find the (much simpler) expressions

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-E\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\>(\<b-E\>*<rsub|\<perp\>\<b-beta\>>(\<b-x\>,t)+\<b-beta\>\<times\>\<b-B\>(\<b-x\>,t))+\<b-E\><rsub|\<\|\|\>\<b-beta\>>(\<b-x\>,t),>>|<row|<cell|\<b-B\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\>(\<b-B\><rsub|\<perp\>\<b-beta\>>(\<b-x\>,t)-\<b-beta\>\<times\>\<b-E\>(\<b-x\>,t))+\<b-B\><rsub|\<\|\|\>\<b-beta\>>(x,t).>>>>
  </eqnarray*>

  Converted to SI, we get

  <\eqnarray*>
    <tformat|<table|<row|<cell|<sqrt|\<varepsilon\><rsub|0>>\<b-E\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\>(<sqrt|\<varepsilon\><rsub|0>>\<b-E\><rsub|\<perp\>\<b-beta\>>*(\<b-x\>,t)+\<b-beta\>\<times\><frac|\<b-B\>|<sqrt|\<mu\><rsub|0>>>(\<b-x\>,t))+<sqrt|\<varepsilon\><rsub|0>>\<b-E\><rsub|\<\|\|\>\<b-beta\>>(\<b-x\>,t),>>|<row|<cell|<frac|\<b-B\><rprime|'>|<sqrt|\<mu\><rsub|0>>>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\><left|(><frac|\<b-B\><rsub|\<perp\>\<b-beta\>>|<sqrt|\<mu\><rsub|0>>>(\<b-x\>,t)-\<b-beta\>\<times\><sqrt|\<varepsilon\><rsub|0>>\<b-E\>(\<b-x\>,t)<right|)>+<frac|\<b-B\><rsub|\<\|\|\>\<b-beta\>>|<sqrt|\<mu\><rsub|0>>>(x,t)<right|)>>>>>
  </eqnarray*>

  or

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-E\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\>(\<b-E\><rsub|\<perp\>\<b-beta\>>*(\<b-x\>,t)+c\<b-beta\>\<times\>\<b-B\>(\<b-x\>,t))-\<b-E\><rsub|\<\|\|\>\<b-beta\>>(\<b-x\>,t),>>|<row|<cell|\<b-B\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\><left|(>\<b-B\><rsub|\<perp\>\<b-beta\>>(\<b-x\>,t)-<frac|\<b-beta\>|c>\<times\>\<b-E\>(\<b-x\>,t)<right|)>-\<b-B\><rsub|\<\|\|\>\<b-beta\>>(x,t).>>>>
  </eqnarray*>

  Next, observe that the situation in the rest frame is purely electrostatic,
  i.e. there is no time dependency and no magnetic field. This simplifies the
  field Lorentz transform to

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-E\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\>\<b-E\><rsub|\<perp\>\<b-beta\>>(\<b-x\>)+\<b-E\><rsub|\<\|\|\>\<b-beta\>>(\<b-x\>),>>|<row|<cell|\<b-B\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\><frac|\<b-beta\>|c>\<times\>\<b-E\>(\<b-x\>),>>>>
  </eqnarray*>

  or, using <math|<wide|\<b-E\>|~>>, to the final form

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-E\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\><wide|\<b-E\>|~><rsub|\<perp\>\<b-beta\>>(\<b-x\><rprime|'>)+<wide|\<b-E\>|~><rsub|\<\|\|\>\<b-beta\>>(\<b-x\><rprime|'>),>>|<row|<cell|\<b-B\><rprime|'>(\<b-x\><rprime|'>,t<rprime|'>)>|<cell|=>|<cell|\<gamma\><frac|\<b-beta\>|c>\<times\><wide|\<b-E\>|~>(\<b-x\><rprime|'>).>>>>
  </eqnarray*>

  It remains to verify that this modified initial condition actually
  satisfies <math|div<rsub|\<b-x\><rprime|'>>\<b-E\><rprime|'>=\<rho\><rprime|'>/\<varepsilon\>>.

  <\eqnarray*>
    <tformat|<table|<row|<cell|div<rsub|\<b-x\><rprime|'>>\<b-E\><rprime|'>(\<b-x\><rprime|'>)>|<cell|=>|<cell|div<rsub|\<b-x\><rprime|'>><left|[>\<gamma\><wide|\<b-E\>|~><rsub|\<perp\>\<b-beta\>>(\<b-x\><rprime|'>)+<wide|\<b-E\>|~><rsub|\<\|\|\>\<b-beta\>>(\<b-x\><rprime|'>)<right|]>>>|<row|<cell|>|<cell|=>|<cell|div<rsub|\<b-x\><rprime|'>><left|[>\<gamma\>\<nabla\><rsub|\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>><wide|\<Phi\>|~>(\<b-x\><rprime|'>)+<frac|\<nabla\><rsub|\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>>|\<gamma\>><wide|\<Phi\>|~>(\<b-x\><rprime|'>)<right|]>>>|<row|<cell|>|<cell|=>|<cell|\<gamma\>\<nabla\><rsub|\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>>\<cdot\>(\<nabla\><rsub|\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>><wide|\<Phi\>|~>(\<b-x\><rprime|'>))+<frac|\<nabla\><rsub|\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>>\<cdot\>(\<nabla\><rsub|\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>><wide|\<Phi\>|~>(\<b-x\><rprime|'>))|\<gamma\>>>>|<row|<cell|>|<cell|=>|<cell|\<gamma\><left|[>\<nabla\><rsub|\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>>\<cdot\>(\<nabla\><rsub|\<b-x\><rsub|\<perp\>\<b-beta\>><rprime|'>><wide|\<Phi\>|~>(\<b-x\><rprime|'>))+<frac|\<nabla\><rsub|\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>>\<cdot\>(\<nabla\><rsub|\<b-x\><rsub|\<\|\|\>\<b-beta\>><rprime|'>><wide|\<Phi\>|~>(\<b-x\><rprime|'>))|\<gamma\><rsup|2>><right|]>>>|<row|<cell|>|<cell|=>|<cell|\<gamma\><left|[><frac|<wide|\<rho\>|~>(\<b-x\><rprime|'>)|\<varepsilon\>><right|]>=<frac|\<rho\><rprime|'>(\<b-x\><rprime|'>)|\<varepsilon\>>.>>>>
  </eqnarray*>

  As a final part of this derivation, we derive an expression for the current
  density seen in the lab frame <math|K<rprime|'>>. We note that in the
  moving frame we have <math|\<b-J\>=(c\<rho\>,\<b-0\>)>, and find

  <\eqnarray*>
    <tformat|<table|<row|<cell|c\<rho\><rprime|'>>|<cell|=>|<cell|\<gamma\>(c\<rho\>),>>|<row|<cell|\<b-J\><rprime|'>>|<cell|=>|<cell|-\<gamma\>\<b-beta\>c\<rho\>=\<gamma\><frac|\<b-v\><rsub|0>|c>c\<rho\>=\<gamma\>\<b-v\><rsub|0>\<rho\>.>>>>
  </eqnarray*> If the charge distribution <math|\<rho\>> only depends on the
  radial coordinate <math|r>, we may compute <math|\<b-E\>> outside the
  charge cloud exactly by using Gauss' law. We know <math|div
  \<b-E\><rprime|'>=\<rho\><rprime|'>/\<varepsilon\>=\<gamma\><wide|\<rho\>|~>/\<varepsilon\>>,
  and find

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|int><rsub|{r\<leqslant\>R}>div\<b-E\><rprime|'>*\<mathd\>V>|<cell|=>|<cell|<big|int><rsub|{r=R}>\<b-E\><rprime|'>\<cdot\>\<b-n\>\<mathd\>A>>|<row|<cell|<big|int><rsub|{r\<leqslant\>R}>\<gamma\><frac|<wide|\<rho\>|~>|\<varepsilon\>>*\<mathd\>V>|<cell|=>|<cell|\<b-E\><rprime|'><rsub|r>(R)\<cdot\>2\<pi\>R>>|<row|<cell|<frac|\<gamma\>|\<varepsilon\>><big|int><rsub|{r\<leqslant\>R}><wide|\<rho\>|~>*\<mathd\>V>|<cell|=>|<cell|\<b-E\><rprime|'><rsub|r>(R)\<cdot\>2\<pi\>R>>|<row|<cell|<frac|\<gamma\>|2\<pi\>R\<varepsilon\>>Q>|<cell|=>|<cell|\<b-E\><rprime|'><rsub|r>(R).>>>>
  </eqnarray*>

  <section|Hyperbolic Cleaning>

  Maxwell's in conservation form reads:

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>E-<frac|1|\<varepsilon\>>\<nabla\>\<times\>H>|<cell|=>|<cell|-<frac|j|\<varepsilon\><rsub|0>>>>|<row|<cell|\<partial\><rsub|t>H+<frac|1|\<mu\>>\<nabla\>\<times\>E>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  \;

  <section|Hyperbolic cleaning of <math|div E>>

  Maxwell's with hyperbolic cleaning is:

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>E-<frac|1|\<varepsilon\>>\<nabla\>\<times\>H+c<rsup|2>\<nabla\>\<Phi\>>|<cell|=>|<cell|-<frac|j|\<varepsilon\><rsub|0>>>>|<row|<cell|\<partial\><rsub|t>H+<frac|1|\<mu\>>\<nabla\>\<times\>E>|<cell|=>|<cell|0>>|<row|<cell|\<partial\><rsub|t>\<Phi\>+\<chi\><rsup|2>\<nabla\>\<cdot\>E>|<cell|=>|<cell|\<chi\><rsup|2><frac|\<rho\>|\<varepsilon\><rsub|0>>>>>>
  </eqnarray*>

  <section|RMS Emittance>

  <math|<with|mode|text|<assign|mean|<macro|arg|\<langle\><arg|arg>\<rangle\>>>>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<varepsilon\><rsub|y><rsup|2>>|<cell|=>|<cell|<mean|(y-<mean|y>)<rsup|2>><mean|(y<rprime|'>-<mean|y<rprime|'>>)<rsup|2>>-<mean|(y-<mean|y>)(y<rprime|'>-<mean|y<rprime|'>>)><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<mean|y<rsup|2>-2y<mean|y>+<mean|y><rsup|2>><mean|(y<rprime|'>-<mean|y<rprime|'>>)<rsup|2>>-<mean|(y-<mean|y>)(y<rprime|'>-<mean|y<rprime|'>>)><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|[<mean|y<rsup|2>>-<mean|y><rsup|2>][<mean|y<rprime|'><rsup|2>>-<mean|y<rprime|'>><rsup|2>]-<mean|y*y<rprime|'>-<mean|y>y<rprime|'>-y<mean|y<rprime|'>>+<mean|y><mean|y<rprime|'>>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|[<mean|y<rsup|2>><mean|y<rprime|'><rsup|2>>-<mean|y><rsup|2><mean|y<rprime|'><rsup|2>>-<mean|y<rsup|2>><mean|y<rprime|'>><rsup|2>+<mean|y><rsup|2><mean|y<rprime|'>><rsup|2>]-[<mean|y*y<rprime|'>>-<mean|y><mean|y<rprime|'>>]<rsup|2>>>|<row|<cell|>|<cell|=>|<cell|[<mean|y<rsup|2>><mean|y<rprime|'><rsup|2>>-<mean|y><rsup|2><mean|y<rprime|'><rsup|2>>-<mean|y<rsup|2>><mean|y<rprime|'>><rsup|2>+<mean|y><rsup|2><mean|y<rprime|'>><rsup|2>]-[<mean|y*y<rprime|'>><rsup|2>-2<mean|y*y<rprime|'>><mean|y><mean|y<rprime|'>>+<mean|y><rsup|2><mean|y<rprime|'>><rsup|2>]>>|<row|<cell|>|<cell|=>|<cell|[<mean|y<rsup|2>><mean|y<rprime|'><rsup|2>>-<mean|y><rsup|2><mean|y<rprime|'><rsup|2>>-<mean|y<rsup|2>><mean|y<rprime|'>><rsup|2>]-[<mean|y*y<rprime|'>><rsup|2>-2<mean|y*y<rprime|'>><mean|y><mean|y<rprime|'>>]>>>>
  </eqnarray*>
</body>

<\initial>
  <\collection>
    <associate|info-flag|detailed>
    <associate|page-type|letter>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|1|?>>
    <associate|auto-4|<tuple|2|?>>
    <associate|auto-5|<tuple|3|?>>
    <associate|auto-6|<tuple|4|?>>
    <associate|auto-7|<tuple|5|?>>
    <associate|fig:dilation|<tuple|2|?>>
    <associate|fig:dilation-lab-flat|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<label|fig:dilation-lab-flat>Illustration of the effects
      of Lorentz transformation on initial conditions.|<pageref|auto-3>>

      <tuple|normal|<label|fig:dilation>Present charge deposition scheme.
      Produces elongated particles in the rest frame.|<pageref|auto-4>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Shape
      function Integrals> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Initial
      Conditions by Lorentz Transform> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>