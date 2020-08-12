#from VarroaPy.VarroaPy.RunVarroaPop import VarroaPop

#Oh, ok, so these treatments when on for a long time (between 6/20-10/22/2014, but we only sample them on 4 days for the purposes of calibration )
#We initialize the model on 6/20 with the starting conditions, then we want to sample the model at 4 days.
#There are also varying #'s of replicates for each treatment condition, which we want to mimic with our model runs. 
#So, we set up the model to run for set number of times. We store all of the data in the appropriate places, 
#then we calculate summary statistics on everything which is returned out of here. This gets compared to the 
#the field data in the abc.new command( a method for ABCSMC), just prior to running the abc.run command (another method). 

#What happens is that for every list of parameter for the eight paramters that we're changing between simulations, a series of 96 simulations
#is carried out based on the treatments and numbers of replicats, each of which is sampled at 4 four dates. 

#So, I think that all we need to do for the beehave model is set up the model to run in series using a loop to go through all of the iterations of the model 
#and 

#Need to figure out what the output dictionary is supposed to look like 
#Need to make sure the function outputs what it's supposed to. 


from itertools import product
import numpy as np
import os
import pandas as pd
import datetime
import pyNetLogo
from random import sample


#Note, os.path.abspath finds the "normalized" absolute path of afile

#Start date for sims = CCA3, 06/20/2014
#Several of the commands below with os.path.foo are specifying particular filenames and directories; __file__ refers to the this script and is used to refer to it by the pyabc_run_BeeHave script
#DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','data')) 
#INITIAL_DF = os.path.join(DATA_DIR, 'bee_counts/initial_conditions.csv')
DATA_DIR = 'C:\\Users\\DDawso01\\repos\\bee_neonic_abc\\data'
INITIAL_DF = os.path.join(DATA_DIR, 'bee_counts/initial_conditions.csv')
#Initial_DF='C:\\Users\\DDawso01\\repos\\bee_neonic_abc\\data'
START_DATE = 171 #'06/20/2014' #These are the julian dates
END_DATE = 295  #'10/22/2014'


DATES = ['4', '5', '6', '7']
DATES_STR  = [197, 220, 253, 288]
DATES_STR_HIGH  = [197, 220, 253, 294]

TREATMENTS = ['0', '10', '20', '40', '80', '160'] 
REPS = [24, 12, 12, 12, 12, 12] #So, for each "replicate" of the model, there should be 96

#netlogo_home = "C:\\Program Files\\NetLogo 5.3.1" #Windows
#netlogo_version = "5"
#netlogo=pyNetLogo.NetLogoLink(gui=True, netlogo_home=netlogo_home, netlogo_version=netlogo_version)



def simulate(param_priors): #Note, need to add in com
#    netlogo_home = "C:\\Program Files\\NetLogo 5.3.1" #Windows
#    netlogo_version = "5"
#    netlogo=pyNetLogo.NetLogoLink(gui=True, netlogo_home=netlogo_home, netlogo_version=netlogo_version)
    #netlogo.load_model("/home/ddawso01/Documents/BeeProject/Beehave_PEEM_Nectar_DED033120.nlogo") #Linux Box
    #netlogo.load_model("C:\Users\DDawso01\DesktopRunFiles\DesktopRunFiles_070220\Beehave_PEEM_Nectar_Sensitivity_DED070220.nlogo")
  #  netlogo.load_model("Beehave_PEEM_Nectar_Sensitivity_DED070220.nlogo")
    #Set up problem set 
    RepExp=2 #This is the number of reps per experiment
    RepCon=2 #This is the number of reps per control
    NonZeroConcs=5
    RepTotal=RepExp*NonZeroConcs + RepCon  
    Experiment=np.concatenate((np.repeat(0, RepCon), np.repeat(10, RepExp),
                       np.repeat(20, RepExp),np.repeat(40,RepExp),
                       np.repeat(80, RepExp),np.repeat(160, RepExp )))

    pars=param_priors.copy()
    startday=np.repeat(171, RepTotal)
    runlength=np.repeat(124, RepTotal)
    #randomseed = sample(range(1,10000), k=RepTotal)
    RunRep = list(range(1, RepTotal + 1)) #Have to put list around range to produce the range of numbers; don't need list inside a loop becuase its an iterator
    ExpRep=np.concatenate((range(1,RepCon +1), range(1,RepExp + 1), range(1,RepExp + 1), range(1,RepExp + 1),range(1,RepExp + 1),range(1,RepExp + 1)))
    AAP1=np.repeat(pars["AdultAcutePar1"], RepTotal)
    AAP2=np.repeat(pars["AdultAcutePar2"], RepTotal)
    LAP1=np.repeat(pars["LarvaAcutePar1"], RepTotal)
    LAP2=np.repeat(pars["LarvaAcutePar2"], RepTotal)
    
    
    param_values = np.column_stack((Experiment,startday, runlength, AAP1, AAP2, LAP1, LAP2,ExpRep,RunRep))
    
    
    #Set up and Run the model 
    Runlist=[]
    for i in range(0,len(param_values)):          
        netlogo.command('setup')
        netlogo.command('set NectarPesticideInput' + " " + str(param_values[i][0]))
        netlogo.command('set StartDay' + " " + str(param_values[i][1]))
        netlogo.command('set X_Days' + " " + str(param_values[i][2]))
        netlogo.command('set AdultAcutePar1' + " " + str(param_values[i][3]))
        netlogo.command('set AdultAcutePar2' + " " + str(param_values[i][4]))
        netlogo.command('set LarvaAcutePar1' + " " + str(param_values[i][5]))
        netlogo.command('set LarvaAcutePar2' + " " + str(param_values[i][6]))
                 
        #netlogo.command('random-seed' + " " + str(param_values[i][3]))        
        X_Days=netlogo.report("X_Days")
        TE=netlogo.report("sum [number] of EggCohorts")
        TL=netlogo.report("sum [number] of LarvaeCohorts")
        TP=netlogo.report("sum [number] of PupaeCohorts")
        TIHB=netlogo.report("sum [number] of IHBeeCohorts") 
        TF=netlogo.report("sum [number] of ForagerSquadrons")
        AFF=netlogo.report("AFF")
        HoneyEnergyStore=netlogo.report("HoneyEnergyStore")
        InitialCond=[X_Days, TE,TL, TP, TIHB, TF, AFF, HoneyEnergyStore]

    
    #Run for the specified number of ticks and return the values of the variables specified at each timestep
    #Note that this outputs a dataframe, which is transferred into a list per output with parallel processing
        counts = netlogo.repeat_report(['TotalEggs', 'TotalIHBees', 
                                    'TotalForagers', 'AFF', 'HoneyEnergyStore', 
                                    'HoneyStorePesticideConc'], X_Days)
        counts.iloc[0,0:4] = InitialCond[0:4] #fills in initial numbers
        if param_values[i][0] == 160: #Note, need to make set this to the correct 
            dates = DATES_STR_HIGH
        else:
            dates = DATES_STR
        
        responses = counts.loc[dates,] #Ok, in Jeff's code, the output had not already been summed per category; we just need to pull out the right columns here

#Ok, now have to get it down to only Eggs and Total Adults 
        responses["TotalAdults"] = sum(responses["TotalIHBees"], responses["TotalForagers"])
        responses["Dates"] = responses.index #Adding a Dates column
        responses["Experiment"]=param_values[i][0]
        responses["Rep"]=param_values[i][4]
        Runlist.append(responses) #End of Loop


    #End Netlogo workspace
    #netlogo.kill_workspace()
    #Process the data into the format        
    Runtab = pd.concat(Runlist)
#Makes the columsn we're interested in numeric; not sure why they became "O" class
    Runtab["TotalAdults"]=pd.to_numeric(Runtab["TotalAdults"])
    Runtab["TotalEggs"]=pd.to_numeric(Runtab["TotalEggs"])
#Make Datelabels 
    Runtab["DateLabels"] = Runtab["Dates"].apply(datelabels)
        #Aggregates by groups
    MeanTabAdults=Runtab.groupby(['Experiment', 'DateLabels'])['TotalAdults'].mean().reset_index()
    SDTabAdults=Runtab.groupby(['Experiment', 'DateLabels'])['TotalAdults'].std().reset_index()
    MeanTabEggs=Runtab.groupby(['Experiment', 'DateLabels'])["TotalEggs"].mean().reset_index()
    SDTabEggs=Runtab.groupby(['Experiment', 'DateLabels'])["TotalEggs"].std().reset_index()
        
        #Adds in the labels we're using for the dictionary
    MeanTabAdults["Labels"]= ['_'.join(x) for x in product(TREATMENTS, DATES)]
    MeanTabAdults["Labels1"] = MeanTabAdults["Labels"] + "_Adults_mean"
    SDTabAdults["Labels"]= ['_'.join(x) for x in product(TREATMENTS, DATES)]
    SDTabAdults["Labels1"] = SDTabAdults["Labels"] + "_Adults_sd"
    MeanTabEggs["Labels"]= ['_'.join(x) for x in product(TREATMENTS, DATES)]
    MeanTabEggs["Labels1"] = MeanTabEggs["Labels"] + "_Eggs_mean"
    SDTabEggs["Labels"]= ['_'.join(x) for x in product(TREATMENTS, DATES)]
    SDTabEggs["Labels1"] = SDTabEggs["Labels"] + "_Eggs_sd"

#Renames columns so they can all be concatenated together
    MeanTabAdults.rename(columns={'TotalAdults': 'Values'}, inplace=True)
    SDTabAdults.rename(columns={'TotalAdults': 'Values'}, inplace=True)
    MeanTabEggs.rename(columns={'TotalEggs': 'Values'}, inplace=True)
    SDTabEggs.rename(columns={'TotalEggs': 'Values'}, inplace=True)

#SOrts data so it will align with the field data
    MeanTabAdults.sort_values(by=["DateLabels", "Experiment"])
    SDTabAdults.sort_values(by=["DateLabels", "Experiment"])
    MeanTabEggs.sort_values(by=["DateLabels", "Experiment"])
    SDTabEggs.sort_values(by=["DateLabels", "Experiment"])

#Concatenates it together
    DataTab=pd.concat([MeanTabAdults, SDTabAdults, MeanTabEggs, SDTabEggs])

#Outputs it as a dictionary 
    zipobj=zip(DataTab["Labels1"], DataTab["Values"])
    output_dict=dict(zipobj)

    return output_dict

def datelabels(series):
    if series == 197:
        return 1
    elif series == 220:
        return 2
    elif series == 253:
        return 3
    else:
        return 4