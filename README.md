# NeuralNetworkiPP

This repository contains supporting information for the main text: "Multi-Scale Constitutive Model of Polypropylene Fibers with Molecular Simulations"

Included are the neural network architecture and weights files for each of the significant material parameters to the Arruda-Boyce Plasticity constitutive model. The parameters are discussed in the main text. The constitutive model was implemented using the PolyUMod and MCalibration commercial software from Veryst Engineering. 

The python file in this repository contains functions which allow for utilization of the neural network model to predict the material parameters for isotactic polypropylene given inputs of Mn (number average molecular weight in g/mol), Mw (weight average molecular weight in g/mol), crystallinity (number between 0 and 1), and Hermans orientation factor (number between 0 and 1). 
