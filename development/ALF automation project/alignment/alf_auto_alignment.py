from alf_data_analysis import *
from alf_data_collection import *
from alf_data_collection_sim import *

##############
# PARAMETERS #
##############
simulate = True

rotation_initial = 0
rotation_step = 5
run_length = 1

lattice_a = 3.82
lattice_b = 3.82
lattice_c = 6.28
lattice_alpha = 90
lattice_beta = 90
lattice_gamma = 90

# TODO create presets based on sample holder material
mask_dspace_min = 0
mask_dspace_max = 2.6

# Peakfinding Parameters
background = 50
resolution_tof = 5000
resolution_phi = 10
resolution_th2 = 10
merge_tolerance = 1

# Set up
print "## SETTING UP"
analyser = AlfDataAnalysis(mask_dspace_min, mask_dspace_max, lattice_a, lattice_b, lattice_c, lattice_alpha,
                           lattice_beta, lattice_gamma, background, resolution_tof, resolution_phi,
                           resolution_th2, merge_tolerance)
if simulate:
    collector = AlfDataCollectionSim()
else:
    collector = AlfDataCollection()

found_ub = False
rotation_current = rotation_initial

print "## STARTING ALIGNMENT"
collector.rotate_and_collect_data(rotation_initial, run_length)

try:
    count = 0
    while not found_ub:
        if analyser.do_next():
            rotation_current += rotation_step
            collector.rotate_and_collect_data(rotation_current, run_length)
            analyser.set_do_next(False)
            count = 0
        count += 1
        if count > 3:
            ans = raw_input("Continue running? (y/n): ")
            if ans.lower() == "n":
                break
            count = 0
        time.sleep(2)

        analyser.process()
        found_ub = analyser.has_found_ub()

    print "UB Matrix found."

except KeyboardInterrupt:
    print "interrupt"
