from GravNN.Visualization.VisualizationBase import VisualizationBase
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from GravNN.Networks.Model import load_config_and_model


class HistoryVisualizer(VisualizationBase):
    def __init__(self, model, config, **kwargs):
        super().__init__(**kwargs)
        self.file_directory += os.path.splitext(os.path.basename(__file__))[0] + "/"
        self.model = model
        self.config = config

        self.get_history()

    def get_history(self):
        self.history = self.config['history'][0]
        self.loss = self.history['loss']
        self.val_loss = self.history['val_loss']
        self.percent_error = self.history['percent_mean']
        self.val_percent_error = self.history['val_percent_mean']
        self.epochs = np.arange(0, len(self.loss), 1)

    def plot(self, x, y, label=None, log_x=False, log_y=False):
        if log_x and log_y:
            plt_fcn = plt.loglog
        elif log_x:
            plt_fcn = plt.semilogx
        elif log_y:
            plt_fcn = plt.semilogy
        else:
            plt_fcn = plt.plot
        plt_fcn(x, y, label=label)
        plt.tight_layout()

    def plot_loss(self, **kwargs):
        plt.figure()
        self.plot(self.epochs, self.loss, label='loss', **kwargs)
        self.plot(self.epochs, self.val_loss, label='val loss', **kwargs)
        plt.legend()

if __name__ == "__main__":
    df = pd.read_pickle("Data/Dataframes/test.data")
    model_id = df["id"].values[-1]
    config, model = load_config_and_model(model_id, df)

    vis = HistoryVisualizer(model, config)
    vis.plot_loss(log_y=True)
    plt.show()