<TeXmacs|1.0.6>

<style|<tuple|generic|maxima>>

<\body>
  Lorentz factor:

  <\equation*>
    \<beta\>=<frac|v|c>,<space|1em>\<gamma\>=<frac|1|<sqrt|1-\<beta\><rsup|2>>>
  </equation*>

  Momentum:

  <\equation*>
    p<rsub|x>=\<gamma\>*m*v<rsub|x><with|mode|text|<space|1em>><with|mode|text|or
    just><space|1em>\<gamma\>\<beta\>
  </equation*>

  \;

  <\equation*>
    x<rprime|'>=<frac|\<mathd\>x|\<mathd\>z>=<frac|v<rsub|x>|v<rsub|z>>
  </equation*>

  Emittance:

  <\equation*>
    \<varepsilon\>=<frac|A|\<pi\>>
  </equation*>

  where <with|mode|math|A> is the area in phase space
  <with|mode|math|(x,x<rprime|'>)>.

  <section|Solution of (1.72) in Chao>

  <with|prog-language|maxima|prog-session|default|<\session>
    <\output>
      \;

      Maxima 5.10.0 http://maxima.sourceforge.net

      Using Lisp GNU Common Lisp (GCL) GCL 2.6.7 (aka GCL)

      Distributed under the GNU Public License. See the file COPYING.

      Dedicated to the memory of William Schelter.

      This is a development version of Maxima. The function bug_report()

      provides bug reporting information.
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>1)
    <with|color|black|>>>
      Kx:0
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o1>)
      <with|color|black|>>0>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>2)
    <with|color|black|>>>
      xi:0
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o2>)
      <with|color|black|>>0>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>3)
    <with|color|black|>>>
      odesoln:ode2('diff(a,s,2)+Kx*a-epsilon[x]^2/a^3=xi/(a+b),a,s)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o17>)
      <with|color|black|>><left|[>-<frac|<sqrt|2*<with|math-font-family|rm|%k1>*a<rsup|2>-1>|2*<with|math-font-family|rm|%k1>*\<varepsilon\><rsub|x>>=s+<with|math-font-family|rm|%k2>,<frac|<sqrt|2*<with|math-font-family|rm|%k1>*a<rsup|2>-1>|2*<with|math-font-family|rm|%k1>*\<varepsilon\><rsub|x>>=s+<with|math-font-family|rm|%k2><right|]>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>21)
    <with|color|black|>>>
      \;
    </input>
  </session>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|<sqrt|2*A*a<rsup|2>-1>|2*A\<varepsilon\><rsub|x>>>|<cell|=>|<cell|s>>|<row|<cell|2*A*a<rsup|2>-1>|<cell|=>|<cell|(2s*A\<varepsilon\><rsub|x>)<rsup|2>>>|<row|<cell|*a<rsup|2>>|<cell|=>|<cell|<frac|(2s*A\<varepsilon\><rsub|x>)<rsup|2>+1|2A>>>>>
  </eqnarray*>
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
  </collection>
</references>