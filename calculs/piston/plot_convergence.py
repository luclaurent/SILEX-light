import numpy
import pylab

h = numpy.array([5.0,
                 3.0
                 ])

ComputationalTime = numpy.array([???,
                                 ???
                                 ])

nbnodes= numpy.array([???,
                ???
                ])

nbelem = numpy.array([???,
                      ???
                      ])

error = numpy.array([?????,
                     ?????
                     ])

MaxDisp=numpy.array([????,
                     ????
	             ])

VMmaxi=numpy.array([???,
                    ???
                    ])

VMmaxiSmooth=numpy.array([???,
                          ???
                          ])

pylab.figure(1)
pylab.plot(nbnodes,numpy.log10(MaxDisp))
pylab.scatter(nbnodes,numpy.log10(MaxDisp))
pylab.xlabel('Nb nodes')
pylab.ylabel('Disp. Maxi')
pylab.grid('on')

pylab.figure(2)
pylab.plot(nbnodes,VMmaxi,label='VM')
pylab.plot(nbnodes,VMmaxiSmooth,label='VM Smooth')
pylab.xlabel('Nb nodes')
pylab.ylabel('V.M. Maxi')
pylab.legend()
pylab.scatter(nbnodes,VMmaxi)
pylab.scatter(nbnodes,VMmaxiSmooth)
pylab.grid('on')

pylab.figure(3)
pylab.plot(numpy.log10(1.0/nbnodes),numpy.log10(error))
pylab.scatter(numpy.log10(1.0/nbnodes),numpy.log10(error))
pylab.xlabel('Average elemental length (log10)')
pylab.ylabel('Global error (log10)')
pylab.title('Global error convergence')
pylab.grid('on')

pylab.figure(4)
pylab.plot(numpy.log10(nbnodes),numpy.log10(ComputationalTime))
pylab.scatter(numpy.log10(nbnodes),numpy.log10(ComputationalTime))
pylab.xlabel('Nb nodes (log10)')
pylab.ylabel('Computational time (log10)')
pylab.grid('on')

pylab.show()

