Using V Curves
==============

Introduction
------------

A "V Curve" is a graph of the size of a star (measured as a half flux radius, or HFD)
versus focus position.  The name refers to the shape of the graph as the star
will shrink as best focus is approached and then grow larger after it is passed.

For a given imaging telescope a V Curve will need to be captured in order to
train the autofocus routine.  Actually several V Curves are normally captured and
then averaged together.

Capturing V Curves With capture_vcurve_script.py
------------------------------------------------

The program "capture_vcurve_script.py" is used for the automated capture of
V Curves.  The program will run a specified number of V Curve captures.

.. code-block:: bash

    usage: capture_vcurve_script.py [-h] [--debugplots] [--simul]
                                    [--focuser_driver FOCUSER_DRIVER]
                                    [--camera_driver CAMERA_DRIVER]
                                    [--exposure_start EXPOSURE_START]
                                    [--exposure_min EXPOSURE_MIN]
                                    [--exposure_max EXPOSURE_MAX]
                                    [--saturation SATURATION]
                                    [--starflux_min STARFLUX_MIN]
                                    [--framesize FRAMESIZE]
                                    [--runoffset RUNOFFSET]
                                    [--hfdcutoff HFDCUTOFF] [--bgthres BGTHRES]
                                    focus_center focus_range focus_nstep focus_dir
                                    nruns

    positional arguments:
      focus_center          Center position of focus run
      focus_range           Range of focus run
      focus_nstep           V Curve number of steps
      focus_dir             IN or OUT
      nruns                 Number of vcurve runs

    optional arguments:
      -h, --help            show this help message and exit
      --debugplots          show debug plots
      --simul               Simulate star
      --focuser_driver FOCUSER_DRIVER
                            Focuser Driver
      --camera_driver CAMERA_DRIVER
                            Camera Driver
      --exposure_start EXPOSURE_START
                            Starting exposure value
      --exposure_min EXPOSURE_MIN
                            Minimum exposure value
      --exposure_max EXPOSURE_MAX
                            Maximum exposure value
      --saturation SATURATION
                            Saturation level for sensor
      --starflux_min STARFLUX_MIN
                            Maximum flux in star
      --framesize FRAMESIZE
                            Size of capture frame, 0=full
      --runoffset RUNOFFSET
                            Shift center of run by this amount
      --hfdcutoff HFDCUTOFF
                            Ignore points with HFD less than this value
      --bgthres BGTHRES     Threshold multiplier for star detection

When run it will create a directory with a named based on the current data and time
and fill it with the images captured at each focus position.  It will also create
a file called "vcurve_fits.json" which contains the fit parameters for each V Curve
captured.  This file is the important output we need.

Averaging V Curves With average_curve_runs.py
---------------------------------------------

Once we have run a series of V Curves we can average them to create a better
estimate of the autofocus parameters needed.  The script "average_vcurve_runs.py"
is used for this purpose.


.. code-block:: bash

    usage: average_vcurve_runs.py [-h] [--debugplots] vcurve_file

    positional arguments:
      vcurve_file   V Curve filename

    optional arguments:
      -h, --help    show this help message and exit
      --debugplots  show debug plots

Basically pass it the path and name of the "vcurve_fits.json" file which is to be averaged.
It will right out a JSON string containing the best fit parameters.
