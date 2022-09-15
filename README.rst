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
################################################################################
1st Validation: Poiseuille
################################################################################
The first step is to check the model's validity with an analytical solution. Therefore, we simulated the Poiseuille flow of non-Newtonian fluid in this study and checked the result with an analytical solution. As shown in below figure, simulation results have excellent match with analytical solution.
--------------------------------------------------------------------------------
.. figure::  /palabos-master/figs/NN-poiseuille.JPG
    :align: center
    :alt: alternate text
    :figclass: align-right
----------------------------------------------------------------------------
