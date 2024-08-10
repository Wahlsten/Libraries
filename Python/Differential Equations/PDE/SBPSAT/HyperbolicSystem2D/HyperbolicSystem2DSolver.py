from GridClass import GridClass
from ProblemClass import ProblemClass
from SBPClass import SBPClass
from SchemeClass import SchemeClass

g = GridClass([0, 0, 0], [1, 1, 1], [10, 10, 1000])
#g.plot_grid()

p = ProblemClass('Maxwell', g)
p.createData(g)

sbp     = SBPClass(g)
sbp.acc = 2

s = SchemeClass(p)
s.ContinuousBoundaryOperators()
s.numericalBoundaryOperators(g, p)
s.createPenalties()
s.CreateMatrices(g, sbp, p)
s.Scheme(g, sbp, p)
s.Error(sbp, g, p, 1)
s.plotSolution(g, p)
