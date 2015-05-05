# pylint: disable=no-init
from mantid.kernel import *
from mantid.api import *
from vesuvio.algorithms.base import VesuvioBase, TableWorkspaceDictionaryFacade
from vesuvio.fitting import parse_fit_options


#----------------------------------------------------------------------------------------


def create_cuboid_xml(height, width, depth):
    """
    Create the XML string to describe a cuboid of the given dimensions

    @param height Height in metres (Y coordinate)
    @param width Width in metres (X coordinate)
    @param depth Depth in metres (Z coordinate)
    """
    half_height, half_width, half_thick = 0.5*height, 0.5*width, 0.5*depth
    xml_str = \
      " <cuboid id=\"sample-shape\"> " \
      + "<left-front-bottom-point x=\"%f\" y=\"%f\" z=\"%f\" /> " % (half_width,-half_height,half_thick) \
      + "<left-front-top-point x=\"%f\" y=\"%f\" z=\"%f\" /> " % (half_width, half_height, half_thick) \
      + "<left-back-bottom-point x=\"%f\" y=\"%f\" z=\"%f\" /> " % (half_width, -half_height, -half_thick) \
      + "<right-front-bottom-point x=\"%f\" y=\"%f\" z=\"%f\" /> " % (-half_width, -half_height, half_thick) \
      + "</cuboid>"
    return xml_str


#----------------------------------------------------------------------------------------


class MSFlags(dict):

    def __init__(self):
        self['AtomProps'] = None


    def to_algorithm_props(self):
        """
        Converts the options to a dictionary of algorithm
        properties for VesuvioCorrections.
        @return Algorithm property dictionary
        """
        alg_props = dict(self)

        # Remove the AtomProps option (processed later)
        if 'AtomProps' in alg_props:
            alg_props.pop('AtomProps')

        # Convert AtomProps to NumMasses and AtomProperties
        alg_props['NumMasses'] = len(self['AtomProps'])
        alg_props['AtomProperties'] = [prop for atom in self['AtomProps'] for prop in atom]

        return alg_props


#----------------------------------------------------------------------------------------


class VesuvioCorrections(VesuvioBase):

    _input_ws = None
    _output_ws = None

    _save_correction = None
    _save_corrected = None
    _correction_workspaces = None
    _corrected_workspaces = None


    def summary(self):
        return "Apply post fitting steps to vesuvio data"

    def PyInit(self):
        # Inputs
        self.declareProperty(MatrixWorkspaceProperty("InputWorkspace", "", direction=Direction.Input),
                             doc="Input TOF workspace")

        self.declareProperty("WorkspaceIndex", 0,
                             doc="Index of spectrum to calculate corrections for")

        # Gamma background
        self.declareProperty("GammaBackground", True, direction=Direction.Input,
                             doc="If true, correct for the gamma background")

        self.declareProperty(ITableWorkspaceProperty("FitParameters", "", direction=Direction.Input,
                                                     optional=PropertyMode.Optional),
                             doc="Table containing the calculated fit parameters for the data in the workspace")

        float_length_validator = FloatArrayLengthValidator()
        float_length_validator.setLengthMin(1)
        self.declareProperty(FloatArrayProperty("Masses", float_length_validator),
                             doc="Mass values for fitting")

        self.declareProperty("MassProfiles", "", StringMandatoryValidator(),
                             doc="Functions used to approximate mass profile. "
                                 "The format is function=Function1Name,param1=val1,param2=val2;function=Function2Name,param3=val3,param4=val4")

        self.declareProperty("IntensityConstraints", "",
                             doc="A semi-colon separated list of intensity constraints defined as lists e.g "
                                 "[0,1,0,-4];[1,0,-2,0]")

        # Multiple scattering
        self.declareProperty("MultipleScattering", True, direction=Direction.Input,
                             doc="If true, correct for the effects of multiple scattering")

        self.declareProperty("BeamRadius", 2.5,
                             doc="Radius of beam in cm")

        self.declareProperty("SampleHeight", 0.0,
                             doc="Height of sample in cm")

        self.declareProperty("SampleWidth", 0.0,
                             doc="Width of sample in cm")

        self.declareProperty("SampleDepth", 0.0,
                             doc="Depth of sample in cm")

        self.declareProperty("NumMasses", 0,
                             doc="")

        self.declareProperty(FloatArrayProperty("AtomProperties", []),
                             doc="")

        self.declareProperty("SampleDensity", 0.0,
                             doc="Sample density in g/cm^3")

        self.declareProperty("Seed", 123456789,
                             doc="")

        self.declareProperty("NumScatters", 3,
                             doc="")

        self.declareProperty("NumRuns", 10,
                             doc="")

        self.declareProperty("NumEvents", 50000,
                             doc="Number of neutron events")

        self.declareProperty("SmoothNeighbours", 3,
                             doc="")

        self.declareProperty("ScatteringScaleFactor", 1.0,
                             doc="")

        # Outputs
        self.declareProperty(WorkspaceGroupProperty("CorrectionWorkspaces", "",
                                                    direction=Direction.Output,
                                                    optional=PropertyMode.Optional),
                             doc="Workspace group containing correction intensities for each correction")

        self.declareProperty(WorkspaceGroupProperty("CorrectedWorkspaces", "",
                                                    direction=Direction.Output,
                                                    optional=PropertyMode.Optional),
                             doc="Workspace group containing sample corrected with each correction")

        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspace", "", direction=Direction.Output),
                             doc="The name of the output workspace")


    def validateInputs(self):
        errors = dict()

        # The gamma correction requires the table of fitted parameters
        if self.getProperty("GammaBackground").value and self.getProperty("FitParameters").value is None:
            errors["FitParameters"] = "Gamma background correction requires a set of parameters from a fit of the data"

        return errors


    def PyExec(self):
        from mantid.simpleapi import (ExtractSingleSpectrum, GroupWorkspaces)

        self._input_ws = self.getPropertyValue("InputWorkspace")
        spec_idx = self.getProperty("WorkspaceIndex").value
        self._output_ws = self.getPropertyValue("OutputWorkspace")

        self._correction_wsg = self.getPropertyValue("CorrectionWorkspaces")
        self._save_correction = self._correction_wsg != ""
        self._corrected_wsg = self.getPropertyValue("CorrectedWorkspaces")
        self._save_corrected = self._corrected_wsg != ""

        ExtractSingleSpectrum(InputWorkspace=self._input_ws,
                              OutputWorkspace=self._output_ws,
                              WorkspaceIndex=spec_idx)

        self._corrected_workspaces = list()
        self._correction_workspaces = list()

        if self.getProperty("GammaBackground").value:
            self._gamma_correction()

        if self.getProperty("MultipleScattering").value:
            self._ms_correction()

        self.setProperty("OutputWorkspace", self._output_ws)

        if self._save_correction:
            GroupWorkspaces(InputWorkspaces=self._correction_workspaces,
                            OutputWorkspace=self._correction_wsg)
            self.setProperty("CorrectionWorkspaces", self._correction_wsg)

        if self._save_corrected:
            GroupWorkspaces(InputWorkspaces=self._corrected_workspaces,
                            OutputWorkspace=self._corrected_wsg)
            self.setProperty("CorrectedWorkspaces", self._corrected_wsg)


    def _gamma_correction(self):
        from mantid.simpleapi import (CalculateGammaBackground, CloneWorkspace, DeleteWorkspace)

        fit_opts = parse_fit_options(mass_values=self.getProperty("Masses").value,
                                     profile_strs=self.getProperty("MassProfiles").value,
                                     constraints_str=self.getProperty("IntensityConstraints").value)

        params_ws_name = self.getPropertyValue("FitParameters")
        params_dict = TableWorkspaceDictionaryFacade(mtd[params_ws_name])
        func_str = fit_opts.create_function_str(params_dict)

        correction_background_ws = str(self._correction_wsg) + "_GammaBackground"
        corrected_background_ws = str(self._corrected_wsg) + "_GammaBackground"

        CalculateGammaBackground(InputWorkspace=self._output_ws,
                                 ComptonFunction=func_str,
                                 BackgroundWorkspace=correction_background_ws,
                                 CorrectedWorkspace=self._output_ws)

        if self. _save_correction:
            self._correction_workspaces.append(correction_background_ws)
        else:
            DeleteWorkspace(correction_background_ws)

        if self. _save_corrected:
            CloneWorkspace(InputWorkspace=self._output_ws,
                           OutputWorkspace=corrected_background_ws)

            self._corrected_workspaces.append(corrected_background_ws)


    def _ms_correction(self):
        """
        Calculates the contributions from multiple scattering
        on the input data from the set of given options
        """
        from mantid.simpleapi import (CalculateMSVesuvio, CreateSampleShape,
                                      DeleteWorkspace, SmoothData, Minus)

        # Create the sample shape
        # Input dimensions are expected in CM
        CreateSampleShape(InputWorkspace=self._output_ws,
                          ShapeXML=create_cuboid_xml(self.getProperty("SampleHeight").value/100.,
                                                     self.getProperty("SampleWidth").value/100.,
                                                     self.getProperty("SampleDepth").value/100.))

        # Massage options into how algorithm expects them
        total_scatter_correction = str(self._correction_wsg) + "_TotalScattering"
        multi_scatter_correction = str(self._correction_wsg) + "_MultipleScattering"

        # Calculation
        CalculateMSVesuvio(InputWorkspace=self._output_ws,
                           NoOfMasses=self.getProperty("NumMasses").value,
                           SampleDensity=self.getProperty("SampleDensity").value,
                           AtomicProperties=self.getPropertyValue("AtomProperties"),
                           BeamRadius=self.getProperty("BeamRadius").value,
                           TotalScatteringWS=total_scatter_correction,
                           MultipleScatteringWS=multi_scatter_correction)

        # Smooth the output
        smooth_neighbours  = self.getProperty("SmoothNeighbours").value
        SmoothData(InputWorkspace=total_scatter_correction,
                   OutputWorkspace=total_scatter_correction,
                   NPoints=smooth_neighbours)
        SmoothData(InputWorkspace=multi_scatter_correction,
                   OutputWorkspace=multi_scatter_correction,
                   NPoints=smooth_neighbours)

        # Optional scale
        scale_factor = self.getProperty("ScatteringScaleFactor").value
        if scale_factor != 1.0:
            Scale(InputWorkspace=total_scatter_correction,
                  OutputWorkspace=total_scatter_correction,
                  Factor=scale_factor,
                  Operation=Multiply)
            Scale(InputWorkspace=multi_scatter_correction,
                  OutputWorkspace=multi_scatter_correction,
                  Factor=scale_factor,
                  Operation=Multiply)

        if self._save_correction:
            self._correction_workspaces.append(total_scatter_correction)
            self._correction_workspaces.append(multi_scatter_correction)
        else:
            DeleteWorkspace(total_scatter_correction)
            DeleteWorkspace(multi_scatter_correction)

        if self._save_corrected:
            Minus(LHSWorkspace=self._output_ws,
                  RHSWorkspace=total_scatter_correction,
                  OutputWorkspace=self._output_ws)

            total_scatter_corrected = str(self._corrected_wsg) + "_TotalScattering"
            multi_scatter_corrected = str(self._corrected_wsg) + "_MultipleScattering"

            Minus(LHSWorkspace=self._input_ws,
                  RHSWorkspace=total_scatter_correction,
                  OutputWorkspace=total_scatter_corrected)
            Minus(LHSWorkspace=self._input_ws,
                  RHSWorkspace=multi_scatter_correction,
                  OutputWorkspace=multi_scatter_corrected)

            self._corrected_workspaces.append(total_scatter_corrected)
            self._corrected_workspaces.append(multi_scatter_corrected)


# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioCorrections)
