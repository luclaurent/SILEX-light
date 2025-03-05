
# Landing gear fork
*French version available [here](fork_readme_fr.md)*

The 3D model can be open in your browser using this [link](https://3dviewer.net#model=https://raw.githubusercontent.com/luclaurent/SILEX-light/main/calculs/landing_gear_fork/fork.step)
  
## Introduction:

*   The landing gear fork is a solid which links the airplane leg to the wheels and to the landing gear compass. The following pictures show it on 2 different legs.  
      
    <img src="../misc/jambe_de_train.png" alt="jambe_de_train" width="200"/>

*   A detailed sketch (A3 format) is given here: [A3 SKETCH](../misc/Plan_fourche_A3.pdf)   
*   A A4 format sketch is also given here: [3D A4 SKETCH](../misc/fourch.jpg)  
*   The CAD of the solid is given with a STEP format (file [`fork.step](landing_gear_fork/fork.step)). 
*   We ask to study the solid with the following conditions:

    *   the bottom plane surface is fixed  
    *   the load from the landing gear compasses is as follow, using the (X,Y,Z) coordinates of the CAD:  
         *   [1000 N, -2000 N , 0] on one "ear",  
         *   [-5000 N, -2000 N , 0] on the other "ear".  
    *   the used material is aluminum (AU-4G), the Young's modulus is 75000 MPa, the Poisson's coefficient is 0.27 and the elastic limit is 240 MPa.  
          
## Indications for the computation:
    
* The file [`fork.step](landing_gear_fork/fork.step) is the CAD model  
* Copy the file [`piston-tet4.geo`](piston/piston-tet4.geo) into `fork-tet4.geo`, then change with your file editor (`gedit`, `notepad++` etc...) the name of the `step` file and delete the other lines        
*   Open `fork-tet4.geo` with `gmsh`:
    *   Define the volume:  `geometry/physical group/add/volume` and select the yellow ball and valid by pressing the `e` key
    *   Define the useful surfaces; for each one, do: `geometry/physical group/add/Surface` and s elect the surface(s), and valid by pressing the `e` key  
    *   Define the size of the elements: `Mesh/Define/Characteristic length /Surface`  
*   Perform all the necessary tasks o have a mesh ready for the Python Main file.  
*   Copy the file [`Main-Piston.py`](piston/Main-Piston.py) (or [`Main-Piston.ipynb`](piston/Main-Piston.ipynb) into `Main-Fork.py` (or `Main-Fork.ipynb`), modify this file in order to perform the computation.  
*   Perform a convergence study  
*   Analyze the results.  


## Report  
*   Give briefly the formulation of the 4-node tetrahedral element.  
*   Explain the programs, specially the parts that you have completed.  
*  Analyze the 2 computations `piston` and `fork`  
*   Comment the curves coming from the convergence study
*   Format of the report is free: paper or pdf, by hand or with the help of a word-processing.  
 


## Results "Landing gear fork": displacements

<img src="../misc/resultats-fourche-disp.png" alt="resultats-fourche-disp" width="200"/>

## Results "Landing gear fork": smooth Von Mises stress

<img src="../misc/resultats-fourche-VMlissee.png" alt="resultats-fourche-VMlissee" width="200"/>

