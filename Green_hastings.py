import numpy as np
from scipy.io import loadmat
from sklearn.metrics import pairwise_distances
from scipy.stats import spearmanr
import pickle
import argparse

exceptions = {}
exceptions['except1'] = Exception("File is not a matlab File, nor a .txt file ")
exceptions['except2'] = Exception("If conectome data is inside matlab file, you should specify the variable name inside the .mat file (see python Green-hastings.py -h)")
exceptions['except3'] = Exception("Matrix not found in matlab file")
exceptions['except4'] = Exception("Data is not a NxN Matrix")


class GreenbergHastings(object):
    
    def __init__(self, data_loc, total, r1, r2, temp_min, temp_max, temp_step, dynamics, matlab_matrix_name, output_dir,output_prefix):
        self.data_loc = data_loc
        self.total = total
        self.r1 = r1
        self.r2 = r2
        self.temp_range = np.round(np.arange(temp_min,temp_max + temp_step,temp_step),4) #avoid strange behaviour with floating points by rounding 
        self.dynamics = dynamics
        self.mmn = matlab_matrix_name
        self.data = []
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self._exceptions = exceptions
        
    def load_data(self):
        if not(self.data_loc[-4:] == '.txt' or self.data_loc[-4:] == '.mat'):
            raise self._exceptions['except1']
        if self.data_loc[-4:] == '.mat' and self.mmn == '':
            raise self._exceptions['except2']
            
        if self.data_loc[-4:] == '.txt':
            data = np.loadtxt(self.data_loc)
            if data.shape[0] != data.shape[1]:
                raise self._exceptions['except4']
            
        if self.data_loc[-4:] == '.mat':
            workspace = loadmat(self.data_loc)
            if not self.mmn in workspace:
                raise self._exceptions['except3']
            data = workspace[self.mmn] 
        
        self.data = data

    def compute_sim(self):
        M2 = self.data
        if self.dynamics == 'Maritan':
            temp = np.sum(M2,axis=1)
            temp[temp==0] =1
            M2 = M2 / temp[:,None]
            
        Ss = {}
        N = self.data.shape[0]
        
        for t in range(len(self.temp_range)):
            T = self.temp_range[t]
            S = np.zeros((self.total,N))
            randi = np.random.randint(0,N,10) #10 regions are activate randomly
            S[0,randi] = 1

            for it in range(1,self.total):
                S[it,:] = S[it-1,:]
                mask1 = S[it - 1,:] ==1
                S[it,mask1] = 2 #Change Active states to refractory period
                R = S[it-1,:] == 2 
                mask2 = (np.random.rand(1,N) < self.r2)[0]
                #Change nodes in Refractory Period to Quiescent State with prob r2
                S[it,R & mask2] = 0;
                Q = (S[it-1,:] == 0)
                mask3 = (np.random.rand(1,N) < self.r1)[0]
                #Spontaneous Activation of Quiescent nodes with prob r1
                S[it,Q & mask3] = 1
                E = S[it-1,:] == 1
                M1 = M2.copy()
                
                # keep non zero values only at rows with quiescent 
                # nodes and columns with excited nodes (at t-1)
                M1[E,:] = 0
                M1[R,:] = 0
                M1[:,Q] = 0
                M1[:,R] = 0
                # sum contributions from excited neighbors
                #(to determine if is going to fire)
                W = M1.sum(axis=1)
                # if sum is larger than threshold(temperature), then fire.
                S[it, (W > T)] = 1;

            Ss[T] = S
            
        self.sim = Ss
        
    def sim_to_files(self):
        for key in self.sim.keys():
            np.savetxt("{}/{}-T{}.txt".format(self.output_dir, self.output_prefix,key),self.sim[key])
            
    def sim_to_pickle(self, out_fname = 'simulation_dict.pkl'):
        with open(out_fname,"wb") as out_file:
            pickle.dump(self.sim,out_file)
        

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-location", required=True,
                      help=".txt or .mat file containing the Conectome Data (NxN) Matrix")
    parser.add_argument("--total", default=1000,
                      help="[Default: 1000] Length of the simulation")
    
    parser.add_argument("--r1", default=0.001,
                      help="[Default: 0.001] Spontaneous activation probability (it set a very small rate of background activity)")
    parser.add_argument("--r2", default=0.2,
                      help="[Default: 0.2] Refractory to Quiescent probability (it determine the lengh of time in which is active)")
    
    parser.add_argument("--temp-min", default=0.01,
                      help="[Default: 0.01] Minimum value for temperature range (activation threshold)")
    
    parser.add_argument("--temp-max", default=0.1,
                      help="[Default: 0.1] Maximum value for temperature range (activation threshold)")
    
    parser.add_argument("--temp-step", default=0.01,
                      help="[Default: 0.01] Increment for temperature range (activation threshold)")
    
    parser.add_argument("--dynamics", default="Haimovici",
                      help="[Default: Haimovici] Set the dynamics between 'Haimovici' and 'Maritan'")
    
    parser.add_argument("--output-dir", default="./",
                      help="[Default: './'] Output directory to export the simulations")
    parser.add_argument("--output-prefix", default="Sim",
                      help="[Default: 'Sim'] prefix for the exported files")
    
    parser.add_argument("--matlab-matrix-name", default="",
                      help="If --data-location is a .mat file, this argument is to parse the matlab cariable to python format")
    
    
    args = parser.parse_args()

    GH = GreenbergHastings(
        data_loc=args.data_location, 
        total=args.total, 
        r1=args.r1, 
        r2=args.r2, 
        temp_min=args.temp_min, 
        temp_max=args.temp_max, 
        temp_step=args.temp_step, 
        dynamics=args.dynamics, 
        matlab_matrix_name=args.matlab_matrix_name,
        output_dir = args.output_dir,
        output_prefix = args.output_prefix,
    )
    
    GH.load_data()
    
    GH.compute_sim()
    GH.sim_to_files()
    
        
        
if __name__ == "__main__":
    main()