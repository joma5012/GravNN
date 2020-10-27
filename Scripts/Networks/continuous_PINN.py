
import os, sys
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.io
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import pickle

from GravNN.Visualization.Grid import Grid
from GravNN.Visualization.VisualizationBase import VisualizationBase
from GravNN.Visualization.MapVisualization import MapVisualization
from GravNN.GravityModels.SphericalHarmonics import SphericalHarmonics
from GravNN.CelestialBodies.Planets import Earth
from GravNN.Trajectories.DHGridDist import DHGridDist
from GravNN.Trajectories.RandomDist import RandomDist
from GravNN.Trajectories.ReducedGridDist import ReducedGridDist
from GravNN.Trajectories.ReducedRandDist import ReducedRandDist

np.random.seed(1234)
tf.set_random_seed(1234)

class PhysicsInformedNN:
    # Initialize the class
    def __init__(self, x0, x1, x2, a0, a1, a2, layers, PINN=True, activation='tanh'):
        
        X = np.concatenate([x0, x1, x2], 1) # N x 3

        self.lb = X.min(0) # min of each components
        self.ub = X.max(0) # max of each components
        
        self.x0 = x0
        self.x1 = x1
        self.x2 = x2
        
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        
        self.layers = layers

        self.PINN = PINN
        self.activation = activation
        
        # Initialize NN
        self.weights, self.biases = self.initialize_NN(layers)  

        # tf placeholders and graph
        self.sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True,
                                                     log_device_placement=True))
        
        self.x0_tf = tf.placeholder(tf.float32, shape=[None, self.x0.shape[1]])
        self.x1_tf = tf.placeholder(tf.float32, shape=[None, self.x1.shape[1]])
        self.x2_tf = tf.placeholder(tf.float32, shape=[None, self.x2.shape[1]])
        
        self.a0_tf = tf.placeholder(tf.float32, shape=[None, self.a0.shape[1]])
        self.a1_tf = tf.placeholder(tf.float32, shape=[None, self.a1.shape[1]])
        self.a2_tf = tf.placeholder(tf.float32, shape=[None, self.a2.shape[1]])

        
        self.a0_pred, self.a1_pred, self.a2_pred, self.f_a0_pred, self.f_a1_pred, self.f_a2_pred, self.U_pred = self.net_NS(self.x0_tf, self.x1_tf, self.x2_tf)
        
        if self.PINN:
            self.loss = tf.reduce_sum(tf.square(self.a0_tf - self.a0_pred)) + \
                        tf.reduce_sum(tf.square(self.a1_tf - self.a1_pred)) + \
                        tf.reduce_sum(tf.square(self.a2_tf - self.a2_pred)) + \
                        tf.reduce_sum(tf.square(self.a0_tf + self.f_a0_pred)) + \
                        tf.reduce_sum(tf.square(self.a1_tf + self.f_a1_pred)) + \
                        tf.reduce_sum(tf.square(self.a2_tf + self.f_a2_pred))
        else:
            self.loss = tf.reduce_sum(tf.square(self.a0_tf - self.a0_pred)) + \
                        tf.reduce_sum(tf.square(self.a1_tf - self.a1_pred)) + \
                        tf.reduce_sum(tf.square(self.a2_tf - self.a2_pred))
 
        self.optimizer = tf.contrib.opt.ScipyOptimizerInterface(self.loss, 
                                                                method = 'L-BFGS-B', 
                                                                options = {'maxiter': 50000,
                                                                           'maxfun': 50000,
                                                                           'maxcor': 50,
                                                                           'maxls': 50,
                                                                           'ftol' : 1.0 * np.finfo(float).eps})
        
        self.optimizer_Adam = tf.train.AdamOptimizer()
        self.train_op_Adam = self.optimizer_Adam.minimize(self.loss)
        
        init = tf.global_variables_initializer()
        self.sess.run(init)
        
    def initialize_NN(self, layers):        
        weights = []
        biases = []
        num_layers = len(layers) 
        for l in range(0,num_layers-1):
            W = self.xavier_init(size=[layers[l], layers[l+1]])
            b = tf.Variable(tf.zeros([1,layers[l+1]], dtype=tf.float32), dtype=tf.float32)
            weights.append(W)
            biases.append(b)        
        return weights, biases
        
    def xavier_init(self, size):
        in_dim = size[0]
        out_dim = size[1]        
        xavier_stddev = np.sqrt(2/(in_dim + out_dim))
        return tf.Variable(tf.truncated_normal([in_dim, out_dim], stddev=xavier_stddev), dtype=tf.float32)
    
    def neural_net(self, X, weights, biases):
        num_layers = len(weights) + 1
        #eps = 1E-6
        H = 2.0*(X - self.lb)/(self.ub - self.lb) - 1.0
        for l in range(0,num_layers-2):
            W = weights[l]
            b = biases[l]
            if self.activation == 'tanh':
                H = tf.tanh(tf.add(tf.matmul(H, W), b))
            if self.activation == 'relu':
                H = tf.nn.relu(tf.add(tf.matmul(H, W), b))
        W = weights[-1]
        b = biases[-1]
        Y = tf.add(tf.matmul(H, W), b)
        return Y
      
    def net_NS(self, x0, x1, x2):
       
        a_and_U = self.neural_net(tf.concat([x0,x1,x2], 1), self.weights, self.biases)
        a0 = a_and_U[:,0:1]
        a1 = a_and_U[:,1:2]
        a2 = a_and_U[:,2:3]

        U = a_and_U[:,3:4]
        
        if self.PINN: 
            U_x0 = tf.gradients(U, x0)[0]
            U_x1 = tf.gradients(U, x1)[0]
            U_x2 = tf.gradients(U, x2)[0]
        else:
            U_x0 = None
            U_x1 = None
            U_x2 = None

        return a0, a1, a2, U_x0, U_x1, U_x2, U
    
    def callback(self, loss):
        print('Loss:', loss)
    
    def train(self, nIter):
        tf_dict = {self.x0_tf: self.x0, 
                    self.x1_tf: self.x1, 
                    self.x2_tf: self.x2,
                    self.a0_tf: self.a0,
                    self.a1_tf: self.a1,
                    self.a2_tf: self.a2}
        
        start_time = time.time()
        for it in range(nIter):
            self.sess.run(self.train_op_Adam, tf_dict)
            
            # Print
            if it % 10 == 0:
                elapsed = time.time() - start_time
                loss_value = self.sess.run(self.loss, tf_dict)
                print('It: %d, Loss: %.3e, Time: %.2f' % 
                      (it, loss_value, elapsed))
                start_time = time.time()
    
        self.optimizer.minimize(self.sess,
                                feed_dict = tf_dict,
                                fetches = [self.loss],
                                loss_callback = self.callback)
    
    def predict(self, x0_star, x1_star, x2_star):
        tf_dict = {self.x0_tf: x0_star, 
                    self.x1_tf: x1_star, 
                    self.x2_tf: x2_star}

        a0_star = self.sess.run(self.a0_pred, tf_dict)
        a1_star = self.sess.run(self.a1_pred, tf_dict)
        a2_star = self.sess.run(self.a2_pred, tf_dict)
        U_star = self.sess.run(self.U_pred, tf_dict)

        return a0_star, a1_star, a2_star, U_star

    
if __name__ == "__main__": 
    planet = Earth()
    model_file = planet.sh_hf_file
    density_deg = 180
    max_deg = 1000
    save = False

    configurations = {
        "config_entire_map_40000" : {
            'N_train' : [40000],
            'PINN_flag' : [False],
            'epochs' : [10], #[400000],
            'radius_max' : [planet.radius + 10.0],
            'layers' : [[3, 40, 40, 40, 40, 40, 40, 40, 40, 4]],
            'acc_noise' : [0.00],
            'deg_removed' : [2]
        },
    }    

    for key, config in configurations.items():

        radius_min = planet.radius

        activation = 'tanh'
        #activation = 'relu'

        df_file = "continuous_results.data"

        #trajectory = ReducedGridDist(planet, radius_min, degree=density_deg, reduction=0.25)
        #radius _max = None

        # trajectory = ReducedRandDist(planet, [radius_min, config['radius_max'][0]], points=15488*4, degree=density_deg, reduction=0.25)
        # map_trajectory = ReducedGridDist(planet, radius_min, degree=density_deg, reduction=0.25)

        #trajectory = RandomDist(planet, [radius_min, config['radius_max'][0]], points=15488*4)
        trajectory = RandomDist(planet, [radius_min, config['radius_max'][0]], points=259200)

        map_trajectory =  DHGridDist(planet, radius_min, degree=density_deg)

        # trajectory = DHGridDist(planet, radius_min, degree=density_deg)
        # map_trajectory = trajectory
        # radius_max = None

        Call_r0_gm = SphericalHarmonics(model_file, degree=max_deg, trajectory=trajectory)
        accelerations = Call_r0_gm.load()

        Clm_r0_gm = SphericalHarmonics(model_file, degree=config['deg_removed'][0], trajectory=trajectory)
        accelerations_Clm = Clm_r0_gm.load()

        x = Call_r0_gm.positions # position (N x 3)
        a = accelerations - accelerations_Clm
        u = None # potential (N x 1)

        idx_x = np.random.choice(x.shape[0], config['N_train'][0], replace=False) 
        
        # Initial Data
        x0_train = x[:,0].reshape(-1,1)[idx_x] #r
        x1_train = x[:,1].reshape(-1,1)[idx_x] #theta
        x2_train = x[:,2].reshape(-1,1)[idx_x] #phi

        a0_train = a[:,0].reshape(-1,1)[idx_x] #a r
        a1_train = a[:,1].reshape(-1,1)[idx_x] #a theta
        a2_train = a[:,2].reshape(-1,1)[idx_x] #a phi

        # Add Noise if interested
        a0_train = a0_train + config['acc_noise'][0]*np.std(a0_train)*np.random.randn(a0_train.shape[0], a0_train.shape[1])
        a1_train = a1_train + config['acc_noise'][0]*np.std(a1_train)*np.random.randn(a1_train.shape[0], a1_train.shape[1])
        a2_train = a2_train + config['acc_noise'][0]*np.std(a2_train)*np.random.randn(a2_train.shape[0], a2_train.shape[1])

        model = PhysicsInformedNN(x0_train, x1_train, x2_train, 
                                    a0_train, a1_train, a2_train, 
                                    config['layers'][0], config['PINN_flag'][0], activation)

        start = time.time()
        model.train(config['epochs'][0])
        time_delta = np.round(time.time() - start, 2)
                    

        ######################################################################
        ############################# Training Stats #########################
        ######################################################################    

        Call_r0_gm = SphericalHarmonics(model_file, degree=max_deg, trajectory=map_trajectory)
        Call_a = Call_r0_gm.load()
        
        Clm_r0_gm = SphericalHarmonics(model_file, degree=config['deg_removed'][0], trajectory=map_trajectory)
        Clm_a = Clm_r0_gm.load()

        x = Call_r0_gm.positions # position (N x 3)
        a = Call_a - Clm_a
    
        x0 = x[:,0].reshape(-1,1) #r
        x1 = x[:,1].reshape(-1,1) #theta
        x2 = x[:,2].reshape(-1,1) #phi

        a0_pred, a1_pred, a2_pred, U_pred = model.predict(x0, x1, x2)

        acc_pred = np.hstack((a0_pred, a1_pred, a2_pred))

        error = np.abs(np.divide((acc_pred - a), a))*100 # Percent Error for each component
        RSE_Call = np.sqrt(np.square(acc_pred - a))

        params = np.sum([np.prod(v.get_shape().as_list()) for v in tf.trainable_variables()])
        timestamp = pd.Timestamp(time.time(), unit='s').round('s').ctime()
        entries = {
            'timetag' : [timestamp],
            'trajectory' : [trajectory.__class__.__name__],
            'radius_min' : [radius_min],
            'train_time' : [time_delta],
            'degree' : [max_deg],
            'activation' : [activation],

            'rse_mean' : [np.mean(RSE_Call)],
            'rse_std' : [np.std(RSE_Call)],
            'rse_median' : [np.median(RSE_Call)],
            'rse_a0_mean' : [np.mean(RSE_Call[:,0])],
            'rse_a1_mean' : [np.mean(RSE_Call[:,1])],
            'rse_a2_mean' : [np.mean(RSE_Call[:,2])],

            'percent_rel_mean' : [np.mean(error)],
            'percent_rel_std' : [np.std(error)], 
            'percent_rel_median' : [np.median(error)],
            'percent_rel_a0_mean' : [np.mean(error[:,0])], 
            'percent_rel_a1_mean' : [np.mean(error[:,1])], 
            'percent_rel_a2_mean' : [np.mean(error[:,2])],

            'params' : [params]
        }
        entries.update(config)

        ######################################################################
        ############################# Testing Stats ##########################
        ######################################################################    

        grid_true = Grid(trajectory=map_trajectory, accelerations=a)
        grid_pred = Grid(trajectory=map_trajectory, accelerations=acc_pred)
        diff = grid_pred - grid_true

        five_sigma_mask = np.where(grid_true.total > 5*np.std(grid_true.total))
        five_sigma_mask_compliment = np.where(grid_true.total < 5*np.std(grid_true.total))
        five_sig_features = diff.total[five_sigma_mask]
        five_sig_features_comp = diff.total[five_sigma_mask_compliment]

        three_sigma_mask = np.where(grid_true.total > 3*np.std(grid_true.total))
        three_sigma_mask_compliment = np.where(grid_true.total < 3*np.std(grid_true.total))
        three_sig_features = diff.total[three_sigma_mask]
        three_sig_features_comp = diff.total[three_sigma_mask_compliment]

        map_stats = {
            'sigma_3_mean' : [np.average(np.sqrt(np.square(three_sig_features)))],
            'sigma_3_std' : [np.std(np.sqrt(np.square(three_sig_features)))],
            'sigma_3_median' : [np.median(np.sqrt(np.square(three_sig_features)))],

            'sigma_3_c_mean' : [np.average(np.sqrt(np.square(three_sig_features_comp)))],
            'sigma_3_c_std' : [np.std(np.sqrt(np.square(three_sig_features_comp)))],
            'sigma_3_c_median' : [np.median(np.sqrt(np.square(three_sig_features_comp)))],

            'sigma_5_mean' : [np.average(np.sqrt(np.square(five_sig_features)))],
            'sigma_5_std' : [np.std(np.sqrt(np.square(five_sig_features)))],
            'sigma_5_median' : [np.median(np.sqrt(np.square(five_sig_features)))],

            'sigma_5_c_mean' : [np.average(np.sqrt(np.square(three_sig_features_comp)))],
            'sigma_5_c_std' : [np.std(np.sqrt(np.square(three_sig_features_comp)))],
            'sigma_5_c_median' : [np.median(np.sqrt(np.square(three_sig_features_comp)))],

            'max_error' : [np.max(np.sqrt(np.square(diff.total)))]
        }
        entries.update(map_stats)
        df = pd.DataFrame().from_dict(entries).set_index('timetag')

        ######################################################################
        ############################# Plotting ###############################
        ######################################################################    

        mapUnit = 'mGal'
        map_vis = MapVisualization(mapUnit)
        plt.rc('text', usetex=False)

        fig_true, ax = map_vis.plot_grid(grid_true.total, "True Grid [mGal]")
        fig_pred, ax = map_vis.plot_grid(grid_pred.total, "NN Grid [mGal]")
        fig_pert, ax = map_vis.plot_grid(diff.total, "Acceleration Difference [mGal]")

        map_vis.fig_size = (5*4,3.5*4)
        fig, ax = map_vis.newFig()
        vlim = [0, np.max(grid_true.total)*10000.0] 
        plt.subplot(311)
        im = map_vis.new_map(grid_true.total, vlim=vlim, log_scale=False)
        map_vis.add_colorbar(im, '[mGal]', vlim)
        
        plt.subplot(312)
        im = map_vis.new_map(grid_pred.total, vlim=vlim, log_scale=False)
        map_vis.add_colorbar(im, '[mGal]', vlim)
        
        plt.subplot(313)
        im = map_vis.new_map(diff.total, vlim=vlim, log_scale=False)
        map_vis.add_colorbar(im, '[mGal]', vlim)

        ######################################################################
        ############################# Saving #################################
        ######################################################################    

        if save: 
            try: 
                df_all = pd.read_pickle(df_file)
                df_all = df_all.append(df)
                df_all.to_pickle(df_file)
            except: 
                df.to_pickle(df_file)

            directory = os.path.abspath('.') +"/Plots/"+ str(pd.Timestamp(timestamp).to_julian_date()) + "/"
            os.makedirs(directory, exist_ok=True)

            map_vis.save(fig_true, directory + "true.pdf")
            map_vis.save(fig_pred, directory + "pred.pdf")
            map_vis.save(fig_pert, directory + "diff.pdf")
            map_vis.save(fig, directory + "all.pdf")

            # with open(directory + "network.data", 'wb') as f:
            #     pickle.dump(model, f)

          
        #plt.show()
        