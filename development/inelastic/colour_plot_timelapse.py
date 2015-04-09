import os
import subprocess

def generate_video(group_workspace,
                   directory=config['defaultsave.directry'],
                   log_names=[],
                   colour_scale=None,
                   frame_rate=10,
                   image_filename_format=r'_%d.png',
                   encoder='avconv'):
    """
    Generates a timelapse video fram w workspace group.

    @param group_workspace Name of the workspace group
    @param directory Directory in which to save frames and video
    @param log_names Sample log names to add to plot as annotation
    @param colour_scale Range for colour axis scale (None for auto)
    @param frame_rate Frame rate of video
    @param image_filename_format Format for image filenames
    @param encoder Encoder utility (avconv or ffmpeg)
    """

    for i, ws in enumerate(mtd[group_workspace]):
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

    frame_filename = os.path.join(directory, image_filename_format)
    video_filename = os.path.join(directory, group_workspace + '.mp4')

    # Convert frames to timelapse video
    subprocess.call([encoder,
                     '-r', str(frame_rate),
                     '-i', frame_filename,
                     '-b:v', '1000k',
                     video_filename])


# Directory to save images and final video in (best to create a new directory for this)
output_directory = os.path.join(config['defaultsave.directory'], 'tl')

generate_video('osiris_scan',
               directory=output_directory,
               log_names=['run_title', 'Stick'],
               frame_rate=5)
