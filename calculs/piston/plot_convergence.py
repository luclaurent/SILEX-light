import scipy
import pylab

h = scipy.array([5.0,
                 3.0
                 ])

ComputationalTime = scipy.array([???,
                                 ???
                                 ])

nbnodes= scipy.array([???,
                ???
                ])

nbelem = scipy.array([???,
                      ???
                      ])

error = scipy.array([?????,
                     ?????
                     ])

MaxDisp=scipy.array([????,
                     ????
	             ])

VMmaxi=scipy.array([???,
                    ???
                    ])

VMmaxiSmooth=scipy.array([???,
                          ???
                          ])

pylab.figure(1)
pylab.plot(nbnodes,scipy.log10(MaxDisp))
pylab.scatter(nbnodes,scipy.log10(MaxDisp))
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
pylab.plot(scipy.log10(1.0/nbnodes),scipy.log10(error))
pylab.scatter(scipy.log10(1.0/nbnodes),scipy.log10(error))
pylab.xlabel('Average elemental length (log10)')
pylab.ylabel('Global error (log10)')
pylab.title('Global error convergence')
pylab.grid('on')

pylab.figure(4)
pylab.plot(scipy.log10(nbnodes),scipy.log10(ComputationalTime))
pylab.scatter(scipy.log10(nbnodes),scipy.log10(ComputationalTime))
pylab.xlabel('Nb nodes (log10)')
pylab.ylabel('Computational time (log10)')
pylab.grid('on')

pylab.show()

