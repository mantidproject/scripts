from genie_python import genie as g



class AlfDataCollection:
    def __init__(self):
        self.gonio_pv = "IN:ALF:CS:SB:Rrot"

    def rotate_and_collect_data(self, angle, duration):
        print "Rotating to " + str(angle) + " deg..."
        g.set_pv(self.gonio_pv, angle)
        g.waitfor_move()
        g.change_title(g.get_title + " rot " + str(angle))
        print "Beginning Run " + g.get_title + ". Collecting data for " + str(duration) + " s."
        g.begin()
        g.waitfor_time(seconds=duration)  # waitfor_uamps better?
        g.end()
        print "Run complete."

    def matrix_to_rot(self, ubmatrix):
        # TODO translate matrix into instructions for goniometer
        # rot1 = ...
        # genie.set_pv(axis1, rot1) ...
        pass
