
# coding: utf-8

# In[7]:


import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from sklearn.metrics import mean_squared_error
from tensorflow.keras.optimizers import SGD
from tensorflow.keras.callbacks import EarlyStopping
import tensorflow.keras.backend as K
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import model_from_json
from scipy import optimize
import csv
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams


# In[8]:


def importNNmodel(iPPNNarchitectures=["model_architecture_Ee.json", "model_architecture_muA.json", "model_architecture_lambdaL.json", "model_architecture_Adk.json", "model_architecture_s0.json", "model_architecture_sss.json", "model_architecture_h.json"], iPPNNweights=["model_weights_Ee.h5", "model_weights_muA.h5", "model_weights_lambdaL.h5", "model_weights_Adk.h5", "model_weights_s0.h5", "model_weights_sss.h5", "model_weights_h.h5"]):
    """ Imports the seven neural network models developed in this work for isotactic polypropylene
    Each neural network outputs one material parameter for the Arruda-Boyce Plasticity constitutive equation"""
    assert len(iPPNNarchitectures) == len(iPPNNweights)
    numNN = len(iPPNNarchitectures)
    NNmodels = []
    for i in range(numNN):
        with open(iPPNNarchitectures[i], 'r') as f:
            model = model_from_json(f.read())
        model.load_weights(iPPNNweights[i])
        NNmodels.append(model)
    return NNmodels


# In[22]:


def useNNmodel(Mns, Mws, crystallinities, HermansOrientationFactors, sample_names=[], savefile=""):
    """This function takes in as input four 1 dimensional numpy arrays or lists: Mns, Mws, crystallinities, and HermansOrientationFactors
    These four arrays are used as inputs to each neural network, 
    Mns is an array of number-average molecular weight values in the units of g/mol
    Mws is an array of weight-average molecular weight values in the units of g/mol
    crystallinities is an array of crystallinity values in the range of 0 to 1
    HermansOrientationFactors is an array of Hermans orientation factors in the range of 0 to 1
    The length of each of the four arrays need to be the same
    The input variable sample_names is optional and if provided should be a list of strings with the length being the same as the length of the other arrays
    
    This function uses the scaling of the input and output values used in this work and gives the material parameters for the Arruda-Boyce Plasticity Constitutive Model 
    as implemented in PolyUMod and MCalibration software from Veryst Engineering
    
    The results are displayed to the screen in this function can also be saved with a csv format if the optional input variable savefile is supplied as a string with the filename
    """
    [model_Ee, model_muA, model_lambdaL, model_Adk, model_s0, model_sss, model_h] = importNNmodel()
    
    assert len(Mns) == len(Mws)
    assert len(Mns) == len(crystallinities)
    assert len(Mns) == len(HermansOrientationFactors)
    
    nPoints = len(Mns)
    input_data = np.zeros((nPoints, 4))
    
    input_data[:,0] = Mns
    input_data[:,1] = Mws
    input_data[:,2] = crystallinities
    input_data[:,3] = HermansOrientationFactors
    print("input_data", input_data)
    # scale the input values 
    
    input_dataS = np.copy(input_data)
    input_dataS[:,0] = input_data[:,0]/100000 
    input_dataS[:,1] = input_data[:,1]/1000000 
    input_dataS[:,2] = input_data[:,2]
    input_dataS[:,3] = input_data[:,3]
    
    
    predS_Ee = model_Ee.predict(input_dataS)
    predS_muA = model_muA.predict(input_dataS)
    predS_lambdaL = model_lambdaL.predict(input_dataS)
    predS_Adk = model_Adk.predict(input_dataS)
    predS_s0 = model_s0.predict(input_dataS)
    predS_sss = model_sss.predict(input_dataS)
    predS_h = model_h.predict(input_dataS)
    
    prediction_Ee = predS_Ee[:,0]*1000
    prediction_muA = predS_muA[:,0]*10 + 0.1 
    prediction_lambdaL = (predS_lambdaL[:,0]*10 - 4)+ 6
    prediction_Adk= predS_Adk[:,0]*1000 + 50
    prediction_s0 = predS_s0[:,0]*100 + 0.1 
    prediction_sss = predS_sss[:,0]*100 + 0.1 
    prediction_h = predS_h[:,0]*100 + 0.1
    
    #print to screen 
    
    if len(sample_names) == len(Mns):
        print ("sample_name input, Mn, Mw, Crystallinity, Hermans Orientation Factor, output, Ee, muA, lambdaL, Adk, s0, sss, h")
        for si,sn in enumerate(sample_names): 
            print(sn + " input", ",",Mns[si], ",", Mws[si], ",", crystallinities[si], ",", HermansOrientationFactors[si], "," , "predicted material parameters", ",", prediction_Ee[si], ",", prediction_muA[si], ",", prediction_lambdaL[si],",", prediction_Adk[si],",", prediction_s0[si],",", prediction_sss[si], ",", prediction_h[si])
            
    
    else:
        for si in range(len(Mns)):
            print("Sample " + str(si) + " input", ",", Mns[si], ",", Mws[si], ",", crystallinities[si], ",", HermansOrientationFactors[si], ",", "predicted material parameters", prediction_Ee[si], ",",prediction_muA[si], ",",prediction_lambdaL[si], ",",prediction_Adk[si], ",", prediction_s0[si], ",", prediction_sss[si], ",", prediction_h[si] )
    
    if savefile != "":
        fout = open(savefile, 'w')
        if len(sample_names) == len(Mns):
            for si,sn in enumerate(sample_names):
                fout.write(sn + " input, " + str(Mns[si]) + ", " + str(Mws[si]) + ", " + str(crystallinities[si]) + ", " + str(HermansOrientationFactors[si]) + ", predicted material parameters, " + str(prediction_Ee[si]) + ", " + str(prediction_muA[si]) + ", " + str(prediction_lambdaL[si]) + ", " + str(prediction_Adk[si]) + ", " + str(prediction_s0[si]) + ", " + str(prediction_sss[si]) + ", " + str(prediction_h[si]) + ",\n" )
        
        else:
            for si in range(len(Mns)):
                 fout.write("Sample " + str(si) + " input, " + str(Mns[si]) + ", " + str(Mws[si]) + ", " + str(crystallinities[si]) + ", " + str(HermansOrientationFactors[si]) + ", predicted material parameters, " + str(prediction_Ee[si]) + ", " + str(prediction_muA[si]) + ", " + str(prediction_lambdaL[si]) + ", " + str(prediction_Adk[si]) + ", " + str(prediction_s0[si]) + ", " + str(prediction_sss[si]) + ", " + str(prediction_h[si]) + ",\n" )

                
        fout.close()
        

    
    
    
    
    
    
    

