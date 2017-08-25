import time
from genie_python import genie_simulate as gs

# TODO inst prefix
prefix = "NDLT882:MOT:GONIO:"


class AlfDataCollectionSim:
    def __init__(self):
        print "Starting in simulation mode"
        self.api = gs.API()
        self.dae = gs.Dae()
        self.gonio_pv = "TE:NDLT882:SIMPLE:VALUE1:SP"

    def rotate_and_collect_data(self, angle, duration):
        self.api.set_pv_value(self.gonio_pv, angle)
        print "Rotating to " + str(angle) + " deg..."
        time.sleep(3)  # waitfor_move
        print "Beginning Run. Collecting data for " + str(duration) + " s."
        self.dae.begin_run()
        time.sleep(duration)
        self.dae.end_run()
        print "Run complete."

    def matrix_to_rot(self, ubmatrix):
        # TODO translate matrix into instructions for goniometer
        # rot1 = ...
        # genie.set_pv(axis1, rot1) ...
        pass
