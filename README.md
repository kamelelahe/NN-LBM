<h1>NN-LBM</h1>
<h2>library for simulating fluid flow in porous media</h2>
<p>This repository contains steps numerical simulation for polymer gel injection in porous media and steps to build ANN to reduce the compuational cost to a fraction of a second</p>

<h3>Methodology</h3>
<p>Since polymer gel in very complex fluid, the research divided into several steps as below to make sure about the reliability of the result.</p>

![alt text](https://github.com/kamelelahe/NN-LBM/blob/master/palabos-master/figs/NN-Methodology.JPG)
<p align=center>figure 1: Research methodology<p>
<h3>1st Validation: Model validation of non-Newtonian fluid flow in simple geometry (Poiseuille)</h3>
<p>The first step is to check the model's validity with an analytical solution. Therefore, we simulated the Poiseuille flow of non-Newtonian fluid in this study and checked the result with an analytical solution. As shown in below figure, simulation results have excellent match with analytical solution.</p>
<p align="center">
<img src="https://github.com/kamelelahe/NN-LBM/blob/master/palabos-master/figs/NN-poiseuille.JPG">
</p>
<p align=center> figure 2: Comparing the simulation result and analytical solution for normalized velocity profile of Poiseuille flow of non-Newtonian fluid. The comparison shows an excellent match</p>

<h3>2nd Validation: Verify the code for the flow of Newtonian fluid in porous media</h3>
<p>It is required to check whether the developed code gives us valid results in porous media. This part is done with the help of Darcy's law. We conducted various simulations with changing Δp  to investigate this relationship. For the geometry of the porous media, we used a square cross-section of the two-dimensional Berea micro-model proposed by Boek (Boek & Venturoli, 2010) presented in Figure 3. The boundary conditions are constant pressure (Zou & He, 1995) at the inlet and outlet (right and left side) and the bounce-back scheme applied in the rock and fluid domain interface.  As Figure 4 indicates, the results show an excellent relationship with the slope of about 0.214, which agrees with the code's estimated permeability.</p>

<p align="center">
 <img src="https://github.com/kamelelahe/NN-LBM/blob/master/palabos-master/figs/Berea2D.jpg" width:"50px" height:"100px">
</p>
<p align=center> Figure 3 Berea pattern used in the simulation. The grains are represented in black, and the pore spaces are in white. The size is 1418 μm by 1418 μm, and the etch depth of the physical unit is 24.54 μm<p>

<p align="center">
<img src="https://github.com/kamelelahe/NN-LBM/blob/master/palabos-master/figs/DarcyLaw.jpg">
</p>
<p align=center>Figure 4 Validation of flow of Newtonian fluid in 2D Berea sample. The slope of this plot indicates permeability.</p>
<h3>Simulation of polymer injection in porous media</h3>
<p>We are now ready to simulate polymer gel injection in porous media. In /palabos-master/source codes/ folder, the developed code for simulating polymer gel inejction in 2D and 3D cases are available.</p>
<p> as it is shown in below video, this code can be used for simulating diiferent fluid type just with varying input parameters. For Newtonian case, we can set power-law index and TD_factor to 1; which means fluid is not either shear depandant or time dependant. For simple Carreau fluid, TD_factor must be set to 1, which means gelation is not occured and we have no changes in pore structre. Finallym for polymer gel injectoion, as shown in the last row of this video, both power-law index and TD_factor are not 1, whic h indicates both shear-dependancy and time-dependancy. Clearly, only for the last case we have changes in pore structure </p>

https://user-images.githubusercontent.com/79846810/190846409-dc48ece4-007d-4e92-bb5e-3c16c081f012.mp4

<p>To take a closer look, Figure 6 summarizes the simulation procedure stages. This figure also indicated that polymer gel had changed the structure of porous media desirably.</P>
<p align="center">
<img src="https://github.com/kamelelahe/NN-LBM/blob/master/palabos-master/figs/changesGeom.JPG">
</p>

The simulation is also done for 3D cases in which 3D rock samples are extracted from https://www.digitalrocksportal.org/ . Here, we also graphically depicted some of the simulation results for 3D cases.
https://user-images.githubusercontent.com/79846810/190847592-c6ca8ef6-53f1-45af-b5ec-6366a6bd1f1f.mp4

<p align="center">
<img src="https://github.com/kamelelahe/NN-LBM/blob/master/palabos-master/figs/Berea3D.JPG">
</p>

<p align="center">
<img src="https://github.com/kamelelahe/NN-LBM/blob/master/palabos-master/figs/bentheimer3DFields.JPG">
</p>


https://user-images.githubusercontent.com/79846810/190848915-c6e7feba-82da-4e57-a1d1-a2f967acfea6.mp4



