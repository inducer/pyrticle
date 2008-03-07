<TeXmacs|1.0.6.12>

<style|<tuple|article|maxima>>

<\body>
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
    <associate|page-orientation|portrait>
    <associate|page-type|letter>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|7>>
    <associate|auto-3|<tuple|3|14>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Hyperbolic
      Cleaning> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Hyperbolic
      cleaning of <with|mode|<quote|math>|div E>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Hyperbolic
      Cleaning of <with|mode|<quote|math>|div E> and
      <with|mode|<quote|math>|div H>> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>