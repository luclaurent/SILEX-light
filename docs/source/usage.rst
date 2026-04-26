Usage
=====

Examples available in this repository:

- piston case: calculs/piston_readme.md
- landing gear fork case: calculs/fork_readme.md

Minimal example
---------------

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tri3_python as tri3

   nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
   elements = np.array([[1, 2, 3]])
   material = np.array([210e9, 0.3, 0.01])

   ik, jk, vk = tri3.stiffnessmatrix(nodes, elements, material)

