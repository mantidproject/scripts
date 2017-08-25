import sys
import time

import numpy as np

from mantid.simpleapi import *
from alf_file_event_handler import AlfFileEventHandler

# Static
ws_mask = "ALF_MASK"
all_peaks = "ALF_PEAKS"
# TODO change to actual log folder
run_folder = "C:\\Instrument\\FileWatcherTestBed\\"
crop_index = 1535


class AlfDataAnalysis(object):
    def __init__(self, mask_min, mask_max, lattice_a, lattice_b, lattice_c, lattice_alpha, lattice_beta, lattice_gamma,
                 background, resolution_tof, resolution_phi, resolution_th2, merge_tolerance):
        self.mask_min = mask_min
        self.mask_max = mask_max

        self.lattice_a = lattice_a
        self.lattice_b = lattice_b
        self.lattice_c = lattice_c
        self.lattice_alpha = lattice_alpha
        self.lattice_beta = lattice_beta
        self.lattice_gamma = lattice_gamma

        self.background = background
        self.resolution_tof = resolution_tof
        self.resolution_phi = resolution_phi
        self.resolution_th2 = resolution_th2

        self.merge_tolerance = merge_tolerance

        self.do_next_run = False
        self.found_ub = False

        self.num_last_peaks = 0
        self.last_matrix = np.ndarray(shape=(3, 3), dtype='double')

        if AnalysisDataService.doesExist(all_peaks):
            DeleteWorkspace(all_peaks)

        self.fw = AlfFileEventHandler()

    @staticmethod
    def mask(ws):
        ''' Apply ALF detector mask to WS '''
        if AnalysisDataService.doesExist(ws_mask):
            MaskDetectors(ws, MaskedWorkspace=ws_mask)
        else:
            num_tubes = 24
            tube_size = 64
            # Mask top 3 and bottom 3 spectra per tube
            MaskDetectors(ws, ','.join(
                ['{0}-{1},{2}-{3}'.format(i * tube_size + 1, i * tube_size + 3, i * tube_size + 62, i * tube_size + 64)
                 for
                 i in range(num_tubes)]))
            # Mask first and last tube
            MaskDetectors(ws, '1-65, 1472-1535')
            ExtractMask(ws, OutputWorkspace=ws_mask)

    @staticmethod
    def calc_ub_delta(ub1, ub2):
        ''' Calculate the difference between two UB matrices '''
        ub_delta = np.zeros(shape=(3, 3), dtype='double')
        for i in range(3):
            for j in range(3):
                ub_delta[i][j] = abs(ub1[i][j] - ub2[i][j])
        return ub_delta

    @staticmethod
    def predict(p1, p2):
        predicted = p1 + p2
        # TODO how to convert into angle
        return predicted

    def do_next(self):
        return self.do_next_run

    def set_do_next(self, b):
        self.do_next_run = b

    def has_found_ub(self):
        return self.found_ub

    def pre_process_ws(self, current_ws):
        self.mask(current_ws)
        CropWorkspace(InputWorkspace=current_ws, StartWorkspaceIndex=0, EndWorkspaceIndex=crop_index,
                      OutputWorkspace=current_ws)

        # Set UB Matrix and goniometer information needed for merge
        SetGoniometer(Workspace=current_ws, Axis0="rrot, 0, 1, 0, 1")

        # Filter out Peaks from sample holder
        ConvertUnits(InputWorkspace=current_ws, Target="dSpacing", OutputWorkspace=current_ws)
        MaskBins(InputWorkspace=current_ws, XMin=self.mask_min, XMax=self.mask_max, OutputWorkspace=current_ws)
        ConvertUnits(InputWorkspace=current_ws, Target="TOF", OutputWorkspace=current_ws)

    def process(self):
        new_runs = self.fw.get_updated()  # TODO check if it's already been processed

        if not len(new_runs) == 0:
            for current_run in new_runs:
                path = run_folder + current_run
                current_run = current_run.split(".")[0]
                current_peaks = "PEAKS_" + current_run

                # Load & Preprocess workspace
                Load(path, OutputWorkspace=current_run)
                self.pre_process_ws(current_run)

                # Find peaks on current run
                FindSXPeaks(current_run, PeakFindingStrategy="AllPeaks", AbsoluteBackground=self.background,
                            ResolutionStrategy="AbsoluteResolution", TofResolution=self.resolution_tof,
                            PhiResolution=self.resolution_phi, TwoThetaResolution=self.resolution_th2,
                            OutputWorkspace=current_peaks)

                # Initialise merged Peaks WS
                if not AnalysisDataService.doesExist(all_peaks):
                    RenameWorkspace(InputWorkspace=current_peaks, OutputWorkspace=all_peaks)
                else:
                    CombinePeaksWorkspaces(LHSWorkspace=all_peaks, RHSWorkspace=current_peaks,
                                           CombineMatchingPeaks=True, Tolerance=self.merge_tolerance,
                                           OutputWorkspace=all_peaks)

                all_peaks_ws = AnalysisDataService[all_peaks]
                num_all_peaks = all_peaks_ws.getNumberPeaks()

                # Find UB Matrix and use it to predict peaks
                num_new_peaks = num_all_peaks - self.num_last_peaks
                if all_peaks_ws.getNumberPeaks() > 1 and num_new_peaks > 0:
                    FindUBUsingLatticeParameters(PeaksWorkspace=all_peaks_ws, a=self.lattice_a, b=self.lattice_b,
                                                 c=self.lattice_c, alpha=self.lattice_alpha, beta=self.lattice_beta,
                                                 gamma=self.lattice_gamma, Tolerance=0.2)
                    # PredictPeaks(InputWorkspace=all_peaks_ws, OutputWorkspace="Peaks_Predicted")  # TODO useful?
                    current_matrix = all_peaks_ws.sample().getOrientedLattice().getUB()

                    # Print tracking stuff
                    dif = self.calc_ub_delta(self.last_matrix, current_matrix)
                    print dif
                    avg = 0.0
                    for j in range(3):
                        for k in range(3):
                            avg += dif[j][k]

                    avg /= 9
                    print "New Peaks:\t" + str(num_new_peaks)
                    print "Delta avg:\t" + str(avg)
                    if num_all_peaks > 1:
                        self.last_matrix = AnalysisDataService[all_peaks].sample().getOrientedLattice().getUB()
                        self.num_last_peaks = num_all_peaks

                        if num_all_peaks == 2:
                            # TODO predict next peak
                            pass
                        if num_all_peaks > 2:
                            script_path = sys.path[0]
                            SaveNexus(all_peaks_ws, script_path + "\Out\Peaks_" + current_run + ".nxs")
                            self.found_ub = True

            self.do_next_run = True

        time.sleep(5)

    def predict_peak(self):
        # TODO
        pass