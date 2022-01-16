from GravNN.CelestialBodies.Asteroids import Eros
from GravNN.Trajectories import RandomAsteroidDist
from GravNN.Preprocessors import DummyScaler, UniformScaler
from sklearn.preprocessing import MinMaxScaler


def get_default_eros_config():
    data_config = {
        "planet": [Eros()],
        "grav_file": [Eros().obj_200k],
        "distribution": [RandomAsteroidDist],
        "N_dist": [20000],
        "N_train": [2500],
        "N_val": [1500],
        "radius_min": [0],
        "radius_max": [Eros().radius * 3],
        "acc_noise": [0.0],
        "basis": [None],
        "mixed_precision": [False],
        "dtype": ["float32"],
        "analytic_truth": ["poly_stats_"],
        "remove_point_mass": [False],  # remove point mass from polyhedral model
        "x_transformer": [UniformScaler(feature_range=(-1, 1))],
        "u_transformer": [UniformScaler(feature_range=(-1, 1))],
        "a_transformer": [UniformScaler(feature_range=(-1, 1))],
        "scale_by": ["non_dim"],
        "dummy_transformer": [DummyScaler()],
        "override" : [False]
    }
    network_config = {
        'PINN_constraint_fcn' : ['pinn_a'],
        "layers": [[3, 20, 20, 20, 20, 20, 20, 20, 20, 3]],
        "activation": ["gelu"],
        "init_file": [None],
        "epochs": [7500],
        "initializer": ["glorot_normal"],
        "optimizer": ["adam"],
        "batch_size": [131072 // 2],
        "learning_rate": [0.001*2],
        "dropout": [0.0],
        "dtype": ["float32"],
        "skip_normalization": [False],
        "lr_anneal": [False],
        "beta" : [0.0],
        "input_layer": [False],
        "network_type": ["sph_pines_traditional"],
        "custom_input_layer": [None],
        "ref_radius" : [Eros().radius],
        'seed' : [0],
        'transformer_units' : [20],
        'normalization_strategy' : ['uniform'], #'radial, uniform
    }

    config = {}
    config.update(data_config)
    config.update(network_config)
    return config