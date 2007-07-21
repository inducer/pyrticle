<TeXmacs|1.0.6>

<style|<tuple|generic|axiom|maxima>>

<\body>
  <section|Notes on PIC>

  \;

  <with|prog-language|axiom|prog-session|default|<\session>
    <\input|<with|color|red|<with|mode|math|\<rightarrow\>> >>
      integrate((l-r^2/l)^3*r^2,r=0..l)
    </input>

    <\output>
      <with|mode|math|math-display|true|<frac|16l<rsup|6>|315><leqno>(7)>

      <axiomtype|Union(f1: OrderedCompletion Expression Integer,...) >
    </output>

    <\input|<with|color|red|<with|mode|math|\<rightarrow\>> >>
      \;
    </input>
  </session>>

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

    <\input|<with|color|red|(<with|math-font-family|rm|%i>7)
    <with|color|black|>>>
      omega(n) := 2*%pi^(n/2)/gamma(n/2)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o7>)
      <with|color|black|>>\<omega\><left|(>n<right|)>:=<frac|2*\<pi\><rsup|<frac|n|2>>|\<Gamma\><left|(><frac|n|2><right|)>>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>8)
    <with|color|black|>>>
      omega(3)
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o13>)
      <with|color|black|>>4*\<pi\>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>14)
    <with|color|black|>>>
      integrate((l-r^2/l)^alpha*spherearea*r^(n-1),r,0,l)
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
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o16>)
      <with|color|black|>><frac|l<rsup|n+\<alpha\>>*\<beta\><left|(><frac|n|2>,\<alpha\>+1<right|)>*<with|math-font-family|rm|spherearea>|2>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>17)
    <with|color|black|>>>
      \;
    </input>
  </session>>

  \;
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Notes
      on PIC> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>