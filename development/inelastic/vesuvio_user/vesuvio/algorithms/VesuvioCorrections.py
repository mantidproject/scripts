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

    _output_ws = None

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
        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspace", "", direction=Direction.Output),
                             doc="The name of the output workspace")


    def validateInputs(self):
        errors = dict()

        # The gamma correction requires the table of fitted parameters
        if self.getProperty("GammaBackground").value and self.getProperty("FitParameters").value is None:
            errors["FitParameters"] = "Gamma background correction requires a set of parameters from a fit of the data"

        return errors


    def PyExec(self):
        from mantid.simpleapi import ExtractSingleSpectrum

        input_ws = self.getPropertyValue("InputWorkspace")
        spec_idx = self.getProperty("WorkspaceIndex").value
        self._output_ws = self.getPropertyValue("OutputWorkspace")

        ExtractSingleSpectrum(InputWorkspace=input_ws,
                              OutputWorkspace=self._output_ws,
                              WorkspaceIndex=spec_idx)

        if self.getProperty("GammaBackground").value:
            self._gamma_correction()

        if self.getProperty("MultipleScattering").value:
            self._ms_correction()

        self.setProperty("OutputWorkspace", self._output_ws)


    def _gamma_correction(self):
        from mantid.simpleapi import (CalculateGammaBackground, DeleteWorkspace)

        fit_opts = parse_fit_options(mass_values=self.getProperty("Masses").value,
                                     profile_strs=self.getProperty("MassProfiles").value,
                                     constraints_str=self.getProperty("IntensityConstraints").value)

        params_ws_name = self.getPropertyValue("FitParameters")
        params_dict = TableWorkspaceDictionaryFacade(mtd[params_ws_name])
        func_str = fit_opts.create_function_str(params_dict)

        num_hist = mtd[self._output_ws].getNumberHistograms()

        CalculateGammaBackground(InputWorkspace=self._output_ws,
                                 ComptonFunction=func_str,
                                 BackgroundWorkspace='__background',
                                 CorrectedWorkspace=self._output_ws)
        DeleteWorkspace('__background')


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
        totalS, multS = str(self._output_ws) + "_totsc", str(self._output_ws) + "_mssc"

        # Calculation
        CalculateMSVesuvio(InputWorkspace=self._output_ws,
                           NoOfMasses=self.getProperty("NumMasses").value,
                           SampleDensity=self.getProperty("SampleDensity").value,
                           AtomicProperties=self.getPropertyValue("AtomProperties"),
                           BeamRadius=self.getProperty("BeamRadius").value,
                           TotalScatteringWS=totalS,
                           MultipleScatteringWS=multS)

        # Smooth the output
        smooth_neighbours  = self.getProperty("SmoothNeighbours").value
        SmoothData(InputWorkspace=totalS,
                   OutputWorkspace=totalS,
                   NPoints=smooth_neighbours)
        SmoothData(InputWorkspace=multS,
                   OutputWorkspace=multS,
                   NPoints=smooth_neighbours)

        # Optional scale
        scale_factor = self.getProperty("ScatteringScaleFactor").value
        if scale_factor != 1.0:
            Scale(InputWorkspace=totalS,
                  OutputWorkspace=totalS,
                  Factor=scale_factor,
                  Operation=Multiply)
            Scale(InputWorkspace=multS,
                  OutputWorkspace=multS,
                  Factor=scale_factor,
                  Operation=Multiply)

        Minus(LHSWorkspace=self._output_ws,
              RHSWorkspace=totalS,
              OutputWorkspace=self._output_ws)

        # Cleanup
        DeleteWorkspace(totalS)
        DeleteWorkspace(multS)


# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioCorrections)
