================================================================================
NN-LBM 
================================================================================
This repository contains steps numerical simulation for polymer gel injection in porous media and steps to build ANN to reduce the compuational cost to a fraction of a second

################################################################################
Methodology
################################################################################
Since polymer gel in very complex fluid, the research divided into several steps as below to make sure about the reliability of the result.
--------------------------------------------------------------------------------
.. figure::  /palabos-master/figs/NN-Methodology.JPG
    :align: right
    :alt: alternate text
    :figclass: align-right
    
    figure 1: Research methodology
################################################################################
1st Validation: Model validation of non-Newtonian fluid flow in simple geometry (Poiseuille)
################################################################################
The first step is to check the model's validity with an analytical solution. Therefore, we simulated the Poiseuille flow of non-Newtonian fluid in this study and checked the result with an analytical solution. As shown in below figure, simulation results have excellent match with analytical solution.
--------------------------------------------------------------------------------
.. figure::  /palabos-master/figs/NN-poiseuille.JPG
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    figure 2: Comparing the simulation result and analytical solution for normalized velocity profile of Poiseuille flow of non-Newtonian fluid of table 1. The comparison shows an excellent match
################################################################################
2nd Validation: Verify the code for the flow of Newtonian fluid in porous media
################################################################################
It is required to check whether the developed code gives us valid results in porous media. This part is done with the help of Darcy's law. We conducted various simulations with changing Δp  to investigate this relationship. For the geometry of the porous media, we used a square cross-section of the two-dimensional Berea micro-model proposed by Boek (Boek & Venturoli, 2010) presented in Figure 3. The boundary conditions are constant pressure (Zou & He, 1995) at the inlet and outlet (right and left side) and the bounce-back scheme applied in the rock and fluid domain interface.  As Figure 4 indicates, the results show an excellent relationship with the slope of about 0.214, which agrees with the code's estimated permeability.
--------------------------------------------------------------------------------
.. figure::  /palabos-master/figs/Berea 2D sample.JPG
    :align: center
    :alt: alternate text
    :figclass: align-center
--------------------------------------------------------------------------------
.. figure::  /palabos-master/figs/DarcyLaw.jpg
    :align: center
    :alt: alternate text
    :figclass: align-center
