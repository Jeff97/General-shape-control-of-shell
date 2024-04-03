Here, we provide some supplementary files for the manuscript. 

- Movie 1 shows the simulated growth process of different samples (2D, 3D and complex 3D cases). To download the movie, one can click the item, then click "View raw". 

- Each folder contains .off files of the parametric plane and target shape. One can open .off files using software [MeshLab](https://www.meshlab.net/). 



# Simulation(needs expertise)

Each folder also contains necessary files to run the simulation using ABAQUS. To run the simulation, please make sure 

- Visual Studio and Parallel Studio are linked to ABAQUS, such that it can be used for subroutine development. [Here](https://www.researchgate.net/publication/349991987_Linking_ABAQUS_20192020_and_Intel_oneAPI_Base_Toolkit_FORTRAN_Compiler) is a linking guide from the internet. 

- The path in UMAT subroutine (see line 135 to 200 in the .for file) has been modified to your own one, since it is set as absolute path. For example

  Before modification:

  ```fortran
          open(301,FILE='C:\Users\12872\Desktop\'//
       &  'Bunny\part1\b1.csv',status="old")
          read(301,*) b1Imp
          close(301)
  ```

  After modification:

  ```fortran
      open(301,FILE='C:\$YourPath$\'//
   &  'Bunny\part1\b1.CSV',status="old")
      read(301,*) b1Imp
      close(301)
  ```

- Submit the job through ABAQUS COMMAND window, for example

  ```fortran
  C:
  cd C:\$YourPath$\Abaqus_Files\Alex_Shocked
  
  abaqus job=Alex_P1 user=Growth-Alex.for cpus=6
  ```

# Simulated results

### Case 1: from 2D to 2D plane

![2Dto2D](https://github.com/Jeff97/General-shape-control-of-shell/blob/main/2Dto2D.jpg "2Dto2D")

### Case 2: from 2D to 3D surface

![2Dto3D](https://github.com/Jeff97/General-shape-control-of-shell/blob/main/2Dto3D.jpg "2Dto3D")

### Case 3: from 3D to 3D surface

![3Dto3D](https://github.com/Jeff97/General-shape-control-of-shell/blob/main/3Dto3D.jpg "3Dto3D")

### Case 4: the parametrization and simulated growth deformation of the complex surfaces

![(a) facial expression from a calm face to a shocked face; (b) growth deformation from a Beetle car to a Taxi car](https://github.com/Jeff97/General-shape-control-of-shell/blob/main/ComplexFace.jpg "(a) facial expression from a calm face to a shocked face; (b) growth deformation from a Beetle car to a Taxi car")

### Case 5: continuous face change

![the continuous change from a robot’s face to a woman’s face, and finally to a men’s face](https://github.com/Jeff97/General-shape-control-of-shell/blob/main/FaceChanging.jpg "the continuous change from a robot’s face to a woman’s face, and finally to a men’s face")
