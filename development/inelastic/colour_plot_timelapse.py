import os
import subprocess

# Directory to save images and final video in (best to create a new directory for this)
# Default to a directory named 'tl' in the default save directory
image_dir = os.path.join(config['defaultsave.directory'], 'tl')
# Name of workspace group to export (workspaces are exported in order displayed in Workspaces dock)
group_workspace = 'NewGroup'
# Frame rate (FPS)
frame_rate = 5
# Scale of Y (colour) axis to be kept common to all plots
colour_scale = [0, 16]
# List of log names to add to the plot2D
log_names = ['run_title']
# Either a name of a workpace to subtract as the background, an index in the workspace group, or None to disable
subtract_background = 0

# Encoder utility to use (avconv works, ffmpeg might)
encoder_utility = 'avconv'
# Image filename pattern
image_format = r'_%d.png'

for i, ws in enumerate(mtd[group_workspace]):
    plot_ws = ws.name()

    if subtract_background is not None:
        background_workspace = None
        if isinstance(subtract_background, str):
            background_workspace = mtd[subtract_background]
        elif isinstance(subtract_background, int):
            background_workspace = mtd[group_workspace][subtract_background]

        if background_workspace is not None:
            plot_ws = '__plot'
            Minus(LHSWorkspace=ws,
                       RHSWOrkspace=background_workspace,
                       OutputWorkspace=plot_ws)

    # Create the plot
    plot = plot2D(plot_ws)

    # Set Y scale
    layer = plot.layer(1)
    layer.setAxisScale(1, *colour_scale)

    if len(log_names) > 0:
        # Generate the legend text
        legend_text = ''
        for log in log_names:
            run = ws.getRun()
            if log in run:
                legend_text += '%s: %s\n' % (log, str(run[log].value))

        # Add the new legend
        layer.newLegend(legend_text)

    # Get the image filename
    image_filename = os.path.join(image_dir, image_format % i)

    # Save image of colour fill plot
    plot.exportImage(image_filename)

    # Close
    plot.close()

    if plot_ws != ws.name():
        DeleteWorkspace(plot_ws)

image_filename_format = os.path.join(image_dir, image_format)
video_filename = os.path.join(image_dir, group_workspace + '.mp4')

# Convert frames to timelapse video
subprocess.call([encoder_utility,
                           '-r', str(frame_rate),
                           '-i', image_filename_format,
                           '-b:v', '1000k',
                           video_filename])
