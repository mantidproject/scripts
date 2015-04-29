import os
import subprocess


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
                   colour_scale=None,
                   frame_rate=10,
                   image_filename_format=r'_%.4d.png',
                   encoder='avconv'):
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

        if len(log_names) > 0:
            # Generate the legend text
            legend_text = ''
            run = ws.getRun()
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

        # Get the image filename
        image_filename = os.path.join(directory, image_filename_format % i)

        # Save image of colour fill plot
        plot.exportImage(image_filename)

        # Close
        plot.close()

    if encoder is not None:
        frame_filename = os.path.join(directory, image_filename_format)
        video_filename = os.path.join(directory, workspace_group + '.mp4')

        # Convert frames to timelapse video
        subprocess.call([encoder,
                         '-r', str(frame_rate),
                         '-i', frame_filename,
                         '-c:v', 'mjpeg',
                         '-q:v', '1',
                         video_filename])


#################
# SETUP OPTIONS #
#################

OPTIONS = dict()

# Name of workspace group containing scan
SCAN_WS = 'osiris_scan'

# Directory to save images and final video in (best to create a new directory for this)
OPTIONS['directory'] = os.path.join(config['defaultsave.directory'], 'tl')

# List of names of log values to add to each plot image
# (time series logs will take the average value)
OPTIONS['log_names'] = ['run_title', 'Stick']

# Maximum and minimum values for the colour scale
OPTIONS['colour_scale'] = [-6.0, 6.0]

# Frame rate for the generated video
OPTIONS['frame_rate'] = 5

# Encoder utility to use to create the video
# Options are: ffmpeg, avconv or None
# None will create the series of images but no video
OPTIONS['encoder'] = None


######################
# PROCESS WORKSPACES #
######################

# Replace any infinate or NaN values
ReplaceSpecialValues(InputWorkspace=SCAN_WS,
                     OutputWorkspace=SCAN_WS,
                     NaNValue=0.0,
                     InfinityValue=0.0)

# First subtract a background run
SCAN_WS = subtract_background(SCAN_WS,
                              workspace=0)

# Generate difference workspace
SCAN_WS = generate_difference(SCAN_WS)

# Create the images and video
generate_video(SCAN_WS, **OPTIONS)
