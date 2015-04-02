from mantid.api import DataProcessorAlgorithm

class VesuvioBase(DataProcessorAlgorithm):

    # There seems to be a problem with Python algorithms
    # defining a __init__ method
    _INST = None

    def _load_and_crop_data(self, runs, spectra,
                            ip_file, diff_mode):
        if spectra == "forward":
            spectra = "{0}-{1}".format(*self._INST.forward_spectra)
        elif spectra == "backward":
            spectra = "{0}-{1}".format(*self._INST.backward_spectra)

        if diff_mode == "double":
            diff_mode = "DoubleDifference"
        else:
            diff_mode = "SingleDifference"

        kwargs = {"Filename": runs,
                  "Mode": diff_mode, "InstrumentParFile": ip_file,
                  "SpectrumList": spectra}
        full_range = self._execute_child_alg("LoadVesuvio", **kwargs)
        return self._execute_child_alg("CropWorkspace", InputWorkspace=full_range,
                                       XMin=self._INST.tof_range[0], XMax=self._INST.tof_range[1])

    # ----------------------------------------------------------------------------------------

    def _execute_child_alg(self, name, **kwargs):
        alg = self.createChildAlgorithm(name)
        for name, value in kwargs.iteritems():
            alg.setProperty(name, value)
        alg.execute()
        outputs = list()
        for name in alg.outputProperties():
            outputs.append(alg.getProperty(name).value)
        if len(outputs) == 1:
            return outputs[0]
        else:
            return tuple(outputs)

# -----------------------------------------------------------------------------------------
# Helper to translate from an table workspace to a dictionary. Should be on the workspace
# really ...
# -----------------------------------------------------------------------------------------
class TableWorkspaceDictionaryFacade(object):
    """
    Allows an underlying table workspace to be treated like a read-only dictionary
    """

    def __init__(self, held_object):
        self._table_ws = held_object

    def __getitem__(self, item):
        for row in self._table_ws:
            if row['Name'] == item:
                return row['Value']
        #endfor
        raise KeyError(str(item))

# -----------------------------------------------------------------------------------------
