import os

import spiceypy as spice

from GravNN.Trajectories.TrajectoryBase import TrajectoryBase


class EphemerisDist(TrajectoryBase):
    def __init__(
        self,
        source,
        target,
        frame,
        start_time,
        end_time,
        sampling_interval,
        meta_kernel="GravNN/Files/Ephmerides/NEAR/near_eros.txt",
        times=None,
        celestial_body=None,
    ):
        """Use python SPICE interface to generate an ephemeris and return the positions of that ephemeris

        Args:
            source (str): object whose trajectory is desired
            target (str): object from which the source trajectory will be defined
            frame (str): the frame in which the trajectory is defined (ex. J2000 or EROS_FIXED)
            start_time (str): trajectory start time in UTC (ex. 'Feb 24, 2000')
            end_time (str): trajectory end time in UTC (ex. 'Feb 24, 2000')
            sampling_interval (int): number of seconds between samples
            meta_kernel (str, optional): SPICE Meta-kernel containing all necessary SPICE binaries and definitions. Defaults to "GravNN/Files/Ephmerides/NEAR/near_eros.txt".
            times (list, optional): Specific times to sample trajectory at (will override start_time and end_time). Defaults to None.
            celestial_body (CelestialBody, optional): Celestial body relevant to the trajectory. Defaults to None.
        """
        self.source = source  # the orbiting object (i.e. spacecraft)
        self.target = target  # the object the source orbits (i.e. asteroid)
        self.frame = frame  # the reference frame of the generated coordinates
        self.start_time = start_time  # UTC i.e. "Jan 07, 2001"
        self.end_time = end_time  # UTC i.e. "Jan 07, 2001"
        self.sampling_interval = sampling_interval  # [s]
        self.times = None
        self.meta_kernel = meta_kernel
        self.celestial_body = celestial_body
        super().__init__()
        pass

    def generate_full_file_directory(self):
        self.trajectory_name = (
            os.path.splitext(os.path.basename(__file__))[0]
            + "/"
            + self.source
            + "_"
            + str(self.target)
            + "_"
            + str(self.frame)
            + "_"
            + str(self.start_time)
            + "_"
            + str(self.end_time)
            + "_"
            + str(self.sampling_interval)
        )

        self.file_directory += self.trajectory_name + "/"
        pass

    def generate(self):
        """Populate spiceypy with necessary kernels, convert UTC times into ET, and
        gather ephemeris positions

        Returns:
            np.array: cartesian position vectors
        """
        spice.furnsh(self.meta_kernel)
        etOne = spice.str2et(self.start_time)
        etTwo = spice.str2et(self.end_time)
        step = int((etTwo - etOne) / self.sampling_interval)
        self.times = [x * (etTwo - etOne) / step + etOne for x in range(step)]
        positions, lightTimes = spice.spkpos(
            self.source,
            self.times,
            self.frame,
            "NONE",
            self.target,
        )
        spice.kclear()
        self.positions = positions * 1000.0  # km -> meters
        return positions


def main():
    import matplotlib.pyplot as plt

    from GravNN.Visualization.VisualizationBase import VisualizationBase

    orbits = [
        "Feb 24, 2000",
        "Mar 03, 2000",
        "Apr 02, 2000",
        "Apr 11, 2000",
        "Apr 22, 2000",
        "Apr 30, 2000",
        "July 07, 2000",
        "July 14, 2000",
        "July 24, 2000",
        "July 31, 2000",
        "Aug 08, 2000",
        "Aug 26, 2000",
        "Sep 05, 2000",
        "Oct 13, 2000",
        "Oct 20, 2000",
        "Oct 25, 2000",
        "Oct 26, 2000",
        "Nov 03, 2000",
        "Dec 07, 2000",
        "Dec 13, 2000",
        "Jan 24, 2001",
        "Jan 28, 2001",
        "Jan 29, 2001",
        "Feb 02, 2001",
        "Feb 06, 2001",
    ]

    # get et values one and two, we could vectorize str2et

    vis = VisualizationBase()

    # https://www.jhuapl.edu/Content/techdigest/pdf/V23-N01/23-01-Holdridge.pdf -- NEAR Orbits
    for i in range(0, len(orbits) - 1):
        [orbits[i], orbits[i + 1]]
        trajectory = EphemerisDist(
            "NEAR",
            "EROS",
            "EROS_FIXED",
            orbits[i],
            orbits[i + 1],
            10 * 60,
        )
        positions = trajectory.positions

        fig, ax = vis.new3DFig()
        ax.plot3D(positions[:, 0], positions[:, 1], positions[:, 2], "gray")

    plt.show()


if __name__ == "__main__":
    main()
