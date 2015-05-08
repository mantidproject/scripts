import os
import subprocess

"""
Script to generate a series of images and timepase video for a series of several
reduced workspaces (e.g. a temperature scan).

In order to output a video either ffmpeg or avconv must be installed and
available on the PATH.
On Ubuntu this can be done using "sudo apt-get install libav-tools" and setting
the encoder option in the script to 'avconv'.
"""


########################
#  WORKFLOW FUNCTIONS  #
########################

def subtract_background(workspace_group,
                        workspace=None):
    """
    Subtracts a background from the scan.

    Workspace parameter can wither be an integer denoting an index in the
    workspace group or a string denoting any workspace loaded in Mantid.

    @param workspace_group Name of the workspace group containing the scan
    @param workspace String or integer denoting workspace to subtract
    @return Name of corrected workspace group
    """
    background_ws = None

    if workspace is not None:
        if isinstance(workspace, int):
            background_ws = mtd[workspace_group][workspace]
        elif isinstance(workspace, str):
            background_ws = mtd[workspace]

    if background_ws is not None:
        out_ws = workspace_group + '_bg_sub'
        Minus(LHSWorkspace=workspace_group,
              RHSWorkspace=background_ws,
              OutputWorkspace=out_ws)
        return out_ws

    return workspace_group


def generate_difference(workspace_group):
    """
    Generates difference workspaces for each of the workspaces
    in the original scan workspace group.

    @param workspace_group Name of the workspace group containing the scan
    @return Name of the difference workspace group
    """
    difference_ws = workspace_group + '_diff'

    ws = mtd[workspace_group]
    if ws.size() < 2:
        return workspace_name

    diff_workspaces = []
    for idx in xrange(ws.size() - 1):
        ws1 = ws[idx]
        ws2 = ws[idx + 1]

        ws_name = ws2.name() + '_minus'
        diff_workspaces.append(ws_name)

        Minus(LHSWorkspace=ws1,
              RHSWorkspace=ws2,
              OutputWorkspace=ws_name)

    GroupWorkspaces(InputWorkspaces=diff_workspaces,
                    OutputWorkspace=difference_ws)

    return difference_ws


def generate_video(workspace_group,
                   directory=config['defaultsave.directry'],
                   log_names=[],
                   title_log_name='run_title',
                   colour_scale=None,
                   frame_rate=10,
                   image_filename_format=r'_%d.png',
                   encoder='avconv',
                   **kwargs):
    """
    Generates a timelapse video fram w workspace group.

    If encoder parameter is not set then only a series of images
    will be created.

    @param workspace_group Name of the workspace group
    @param directory Directory in which to save frames and video
    @param log_names Sample log names to add to plot as annotation
    @param colour_scale Range for colour axis scale (None for auto)
    @param frame_rate Frame rate of video
    @param image_filename_format Format for image filenames
    @param encoder Encoder utility (avconv or ffmpeg)
    """

    for i, ws in enumerate(mtd[workspace_group]):
        # Create the plot
        plot = plot2D(ws)
        layer = plot.layer(1)

        # Set Y scale
        if colour_scale is not None:
            layer.setAxisScale(1, *colour_scale)

        run = ws.getRun()

        if len(log_names) > 0:
            # Generate the legend text
            legend_text = ''
            for log in log_names:
                if log in run:
                    entry  = run[log]

                    # Use average value for time series logs
                    if isinstance(entry, FloatTimeSeriesProperty):
                        value = run[log].timeAverageValue()
                    else:
                        value = run[log].value

                    legend_text += '%s: %s\n' % (log, str(value))

            # Add the new legend
            layer.newLegend(legend_text)

        if title_log_name in run:
            plot_title = run[title_log_name].value
            layer.setTitle(plot_title)

        # Get the image filename
        image_filename = os.path.join(directory, image_filename_format % i)

        # Save image of colour fill plot
        plot.exportImage(image_filename)

        # Close
        plot.close()

    if encoder is not None and encoder != '':
        frame_filename = os.path.join(directory, image_filename_format)
        video_filename = os.path.join(directory, workspace_group + '.mp4')

        # Convert frames to timelapse video
        subprocess.call([encoder,
                         '-r', str(frame_rate),
                         '-i', frame_filename,
                         '-c:v', 'mjpeg',
                         '-q:v', '1',
                         video_filename])

        print 'Video file saved to: %s' % (video_filename)


def process_scan(options):
    """
    Contains the workflow for the data processing and video generation.

    @param options Dictionary containing options
    """

    workspace = options['workspace']

    if OPTIONS['replace_bad_values']:
        # Replace any infinate or NaN values
        out_ws_name = workspace + '_clean'
        ReplaceSpecialValues(InputWorkspace=workspace,
                             OutputWorkspace=out_ws_name,
                             NaNValue=0.0,
                             InfinityValue=0.0)
        workspace = out_ws_name

    if OPTIONS['background'] is not None:
        # First subtract a background run
        workspace = subtract_background(workspace,
                                        workspace=OPTIONS['background'])

    if OPTIONS['mode'] == 'difference':
        # Generate difference workspace
        workspace = generate_difference(workspace)

    # Create the images and video
    generate_video(workspace, **OPTIONS)


#############
#  OPTIONS  #
#############

OPTIONS = dict()

# Name of workspace group containing scan
# Workspaces will appear in the output in the sam eorder as they appear in the WorkspaceGroup
OPTIONS['workspace'] = 'MultiFiles'

# Toggle replacement of special values (infinity and NaN) in the input data
OPTIONS['replace_bad_values'] = True

# Remove a background workspace
# Can be either an integer denoting an index in the input workspace group
# or a string denoting the name of a workspace
OPTIONS['background'] = 0

# Mode to plot in, options are:
# normal: plot each of the original workspaces independently
# difference: plots the difference between each pair of workspaces
OPTIONS['mode'] = 'difference'

# Directory to save images and final video in (best to create a new directory for this)
OPTIONS['directory'] = os.path.join(config['defaultsave.directory'], 'tl')

# List of names of log values to add to each plot image
# (time series logs will take the average value)
OPTIONS['log_names'] = ['Stick']

# Name of log entry to use as the title for each plot generated
# If not specified the run title (log name run_title) will be used
#OPTIONS['title_log'] = 'run_title'

# Maximum and minimum values for the colour scale
# Set to None to have this automatically adjust on a per image basis (not recommended as
# the scale will be inconsistent across each of the images)
OPTIONS['colour_scale'] = [-6.0, 6.0]

# Frame rate for the generated video
OPTIONS['frame_rate'] = 5

# Encoder utility to use to create the video
# Options are: ffmpeg, avconv or None
# None will create the series of images but no video
OPTIONS['encoder'] = 'avconv'

process_scan(OPTIONS)
