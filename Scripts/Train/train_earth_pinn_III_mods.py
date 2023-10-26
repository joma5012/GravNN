import os
from pprint import pprint

from GravNN.Networks.Configs import *
from GravNN.Networks.script_utils import save_training
from GravNN.Networks.utils import configure_run_args

os.environ["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"


def run(config):
    from GravNN.Networks.Data import DataSet
    from GravNN.Networks.Model import PINNGravityModel
    from GravNN.Networks.Saver import ModelSaver
    from GravNN.Networks.utils import configure_tensorflow, populate_config_objects

    configure_tensorflow(config)

    # Standardize Configuration
    config = populate_config_objects(config)
    pprint(config)

    # Get data, network, optimizer, and generate model
    data = DataSet(config)
    model = PINNGravityModel(config)
    history = model.train(data)

    saver = ModelSaver(model, history)
    saver.save(df_file=None)

    print(f"Model ID: [{model.config['id']}]")
    return model.config


def run_and_save(df_file, hparams, config):
    args = configure_run_args(config, hparams)
    # return
    configs = [run(args[0][0])]
    save_training(df_file, configs)


def main():
    config = get_default_eros_config()
    config.update(PINN_III())
    config.update(ReduceLrOnPlateauConfig())
    hparams = {
        "N_dist": [6000],
        "N_train": [1000],
        "N_val": [1000],
        "num_units": [20],
        "radius_max": [Eros().radius * 15],
        "loss_fcns": [["mse"]],
        "jit_compile": [True],
        "eager": [False],
        # "jit_compile": [False],
        # "eager": [True],
        "lr_anneal": [False],
        "learning_rate": [0.001],
        "dropout": [0.0],
        "batch_size": [2**18],
        "epochs": [5000],
        "acc_noise": [0.0],
        "preprocessing": [["pines", "r_inv"]],
        "PINN_constraint_fcn": ["pinn_a"],
        "scale_nn_potential": [False],
        "trainable": [False],
        "fuse_models": [False],
        "enforce_bc": [False],
        "uniform_volume": [False],
    }

    df_file = "Data/Dataframes/pinn_III_mods_RMS.data"
    # run_and_save(df_file, hparams, config)

    df_file = "Data/Dataframes/pinn_III_mods_percent.data"
    hparams.update({"loss_fcns": [["percent"]]})
    # run_and_save(df_file, hparams, config)

    df_file = "Data/Dataframes/pinn_III_mods_scaling.data"
    hparams.update({"scale_nn_potential": [True]})
    run_and_save(df_file, hparams, config)

    df_file = "Data/Dataframes/pinn_III_mods_BC.data"
    hparams.update({"enforce_bc": [True], "tanh_k": [0.1]})
    run_and_save(df_file, hparams, config)

    df_file = "Data/Dataframes/pinn_III_mods_fuse.data"
    hparams.update({"fuse_models": [True]})
    run_and_save(df_file, hparams, config)


if __name__ == "__main__":
    main()
