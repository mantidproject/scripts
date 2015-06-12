# pylint: disable=no-init
from string import Template
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

class VesuvioCorrections(VesuvioBase):

    _input_ws = None
    _output_ws = None
    _correction_workspaces = None
    _linear_fit_table = None

#------------------------------------------------------------------------------

    def summary(self):
        return "Apply post fitting steps to vesuvio data"

#------------------------------------------------------------------------------

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

        self.declareProperty("SampleHeight", 5.0,
                             doc="Height of sample in cm")

        self.declareProperty("SampleWidth", 5.0,
                             doc="Width of sample in cm")

        self.declareProperty("SampleDepth", 5.0,
                             doc="Depth of sample in cm")

        self.declareProperty("SampleDensity", 1.0,
                             doc="Sample density in g/cm^3")

        self.declareProperty("Seed", 123456789,
                             doc="")

        self.declareProperty("NumScatters", 3,
                             doc="")

        self.declareProperty("NumRuns", 10,
                             doc="")

        self.declareProperty("NumEvents", 1000000,
                             doc="Number of neutron events")

        self.declareProperty("SmoothNeighbours", 3,
                             doc="")

        # Outputs
        self.declareProperty(WorkspaceGroupProperty("CorrectionWorkspaces", "",
                                                    direction=Direction.Output,
                                                    optional=PropertyMode.Optional),
                             doc="Workspace group containing correction intensities for each correction")

        self.declareProperty(ITableWorkspaceProperty("LinearFitResult", "",
                                                     direction=Direction.Output,
                                                     optional=PropertyMode.Optional),
                             doc="Table workspace containing the fit parameters used to"
                                 "linearly fit the corrections to the data")

        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspace", "", direction=Direction.Output),
                             doc="The name of the output workspace")

#------------------------------------------------------------------------------

    def validateInputs(self):
        errors = dict()

        if self.getProperty("FitParameters").value is None:
            errors["FitParameters"] = "Corrections require a set of parameters from a fit of the data"

        return errors

#------------------------------------------------------------------------------

    def PyExec(self):
        from mantid.simpleapi import (ExtractSingleSpectrum, GroupWorkspaces,
                                      Scale, Minus, DeleteWorkspace)

        self._input_ws = self.getPropertyValue("InputWorkspace")
        spec_idx = self.getProperty("WorkspaceIndex").value
        self._output_ws = self.getPropertyValue("OutputWorkspace")

        self._correction_wsg = self.getPropertyValue("CorrectionWorkspaces")
        self._linear_fit_table = self.getPropertyValue("LinearFitResult")

        ExtractSingleSpectrum(InputWorkspace=self._input_ws,
                              OutputWorkspace=self._output_ws,
                              WorkspaceIndex=spec_idx)

        self._correction_workspaces = list()

        # Do gamma correction
        if self.getProperty("GammaBackground").value:
            self._correction_workspaces.append(self._gamma_correction())

        # Do multiple scattering correction
        if self.getProperty("MultipleScattering").value:
            self._correction_workspaces.extend(self._ms_correction())

        # Output correction workspaces as a WorkspaceGroup
        if self._correction_wsg != "":
            GroupWorkspaces(InputWorkspaces=self._correction_workspaces,
                            OutputWorkspace=self._correction_wsg)
            self.setProperty("CorrectionWorkspaces", self._correction_wsg)

        # Perform fitting to obtain scale factors for corrections
        if self._linear_fit_table != "":
            # The workspaces to fit for correction scale factors
            fit_corrections = [wks for wks in self._correction_workspaces if 'MultipleScattering' not in wks]

            # Perform fitting of corrections
            params_ws = self._fit_corrections(fit_corrections, self._linear_fit_table)
            self.setProperty("LinearFitResult", params_ws)

            # Scale gamma background
            if self.getProperty("GammaBackground").value:
                gamma_correct_ws = self._get_correction_workspace('GammaBackground')[1]
                gamma_factor = self._get_correction_scale_factor('GammaBackground', fit_corrections, params_ws)
                Scale(InputWorkspace=gamma_correct_ws,
                      OutputWorkspace=gamma_correct_ws,
                      Factor=gamma_factor)

            # Scale multiple scattering
            if self.getProperty("MultipleScattering").value:
                # Use factor of total scattering as this includes single and multiple scattering
                multi_scatter_correct_ws = self._get_correction_workspace('MultipleScattering')[1]
                total_scatter_factor = self._get_correction_scale_factor('TotalScattering', fit_corrections, params_ws)
                Scale(InputWorkspace=multi_scatter_correct_ws,
                      OutputWorkspace=multi_scatter_correct_ws,
                      Factor=total_scatter_factor)

        # Apply corrections
        for correction in self. _correction_workspaces:
            if 'TotalScattering' not in correction:
                Minus(LHSWorkspace=self._output_ws,
                      RHSWorkspace=correction,
                      OutputWorkspace=self._output_ws)

        self.setProperty("OutputWorkspace", self._output_ws)

        # Remove correction workspaces if they are no longer required
        if self._correction_wsg == "":
            for wksp in self._correction_workspaces:
                DeleteWorkspace(wksp)

#------------------------------------------------------------------------------

    def _fit_corrections(self, fit_workspaces, param_table_name):
        func_template = Template("name=TabulatedFunction,Workspace=$ws_name,ties=(Shift=0,XScaling=1)")
        functions = [func_template.substitute(ws_name=wsn) for wsn in fit_workspaces]

        fit = AlgorithmManager.create("Fit")
        fit.initialize()
        fit.setChild(True)
        fit.setProperty("Function", ";".join(functions))
        fit.setProperty("InputWorkspace", self._input_ws)
        fit.setProperty("Output", param_table_name)
        fit.execute()

        return fit.getProperty('OutputParameters').value

#------------------------------------------------------------------------------

    def _get_correction_workspace(self, correction_name, corrections=None):
        if corrections is None:
            corrections = self._correction_workspaces

        for idx, ws_name in enumerate(corrections):
            if correction_name in ws_name:
                return idx, ws_name

        return None, None

#------------------------------------------------------------------------------

    def _get_correction_scale_factor(self, correction_name, corrections, params_ws):
        index = self._get_correction_workspace(correction_name, corrections)[0]
        if index is None:
            raise RuntimeError('No workspace for given correction')

        params_dict = TableWorkspaceDictionaryFacade(params_ws)
        scale_param_name = 'f%d.Scaling' % index

        return params_dict[scale_param_name]

#------------------------------------------------------------------------------

    def _gamma_correction(self):
        from mantid.simpleapi import (CalculateGammaBackground, CloneWorkspace,
                                      DeleteWorkspace)
        correction_background_ws = str(self._correction_wsg) + "_GammaBackground"

        fit_opts = parse_fit_options(mass_values=self.getProperty("Masses").value,
                                     profile_strs=self.getProperty("MassProfiles").value,
                                     constraints_str=self.getProperty("IntensityConstraints").value)
        params_ws_name = self.getPropertyValue("FitParameters")
        params_dict = TableWorkspaceDictionaryFacade(mtd[params_ws_name])
        func_str = fit_opts.create_function_str(params_dict)

        CalculateGammaBackground(InputWorkspace=self._output_ws,
                                 ComptonFunction=func_str,
                                 BackgroundWorkspace=correction_background_ws,
                                 CorrectedWorkspace='__corrected_dummy')
        DeleteWorkspace('__corrected_dummy')

        return correction_background_ws

#------------------------------------------------------------------------------

    def _ms_correction(self):
        """
        Calculates the contributions from multiple scattering
        on the input data from the set of given options
        """
        from mantid.simpleapi import (CalculateMSVesuvio, CreateSampleShape,
                                      DeleteWorkspace, SmoothData, Minus)

        masses = self.getProperty("Masses").value
        params_ws_name = self.getPropertyValue("FitParameters")
        params_dict = TableWorkspaceDictionaryFacade(mtd[params_ws_name])

        atom_props = list()
        for i, mass in enumerate(masses):
            intentisty_prop = 'f%d.Intensity' % i
            c0_prop = 'f%d.C_0' % i

            if intentisty_prop in params_dict:
                intentisy = params_dict[intentisty_prop]
            elif c0_prop in params_dict:
                intentisy = params_dict[c0_prop]
            else:
                continue

            width = params_dict['f%d.Width' % i]

            atom_props.append(mass)
            atom_props.append(intentisy)
            atom_props.append(width)

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
                           NoOfMasses=len(atom_props)/3,
                           SampleDensity=self.getProperty("SampleDensity").value,
                           AtomicProperties=atom_props,
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

        return total_scatter_correction, multi_scatter_correction


# -----------------------------------------------------------------------------------------
AlgorithmFactory.subscribe(VesuvioCorrections)
