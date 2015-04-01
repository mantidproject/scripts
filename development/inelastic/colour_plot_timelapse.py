import os
import subprocess

# Directory to save images and final video in (best to create a new directory for this)
# Default to a directory named 'tl' in the default save directory
image_dir = os.path.join(config['defaultsave.directory'], 'tl')
# Name of workspace group to export (workspaces are exported in order displayed in Workspaces dock)
group_workspace = 'elastic_scan_rebin'
# Frame rate (FPS)
frame_rate = 5
# Scale of Y (colour) axis to be kept common to all plots
colour_scale = [0, 5]
# List of log names to add to the plot2D
log_names = ['run_number', 'Temperature']

# Encoder utility to use (avconv works, ffmpeg might)
encoder_utility = 'avconv'
# Image filename pattern
image_format = r'_%d.png'

for i, ws in enumerate(mtd[group_workspace]):
    # Create the plot
    plot = plot2D(ws)

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

image_filename_format = os.path.join(image_dir, image_format)
video_filename = os.path.join(image_dir, group_workspace + '.mp4')

# Convert frames to timelapse video
subprocess.call([encoder_utility,
                           '-r', str(frame_rate),
                           '-i', image_filename_format,
                           '-b:v', '1000k',
                           video_filename])
