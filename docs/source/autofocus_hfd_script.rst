Using autofocus_hfd_script.py
=============================

Introduction
------------

The script autofocus_hfd_script.py handles focusing on a bright star near the
center of an image.

Invocation
----------

The invocation of autofocus_auto_star.py is:

.. code-block:: bash

    usage: autofocus_hfd_script.py [-h] [--focus_min FOCUS_MIN]
                                   [--focus_max FOCUS_MAX] [--focus_dir FOCUS_DIR]
                                   [--focus_start FOCUS_START] [--debugplots]
                                   [--debugplotsdelay DEBUGPLOTSDELAY] [--simul]
                                   [--stayopen] [--profile PROFILE]
                                   [--focuser FOCUSER] [--camera CAMERA]
                                   [--exposure_start EXPOSURE_START]
                                   [--exposure_min EXPOSURE_MIN]
                                   [--exposure_max EXPOSURE_MAX]
                                   [--starflux_min STARFLUX_MIN]
                                   [--saturation SATURATION]
                                   [--framesize FRAMESIZE] [--winsize WINSIZE]
                                   [--focusdelay FOCUSDELAY]
                                   [--numaverage NUMAVERAGE]

    optional arguments:
      -h, --help            show this help message and exit
      --focus_min FOCUS_MIN
                            Lowest focus travel allowed
      --focus_max FOCUS_MAX
                            Highest focus travel allowed
      --focus_dir FOCUS_DIR
                            IN or OUT
      --focus_start FOCUS_START
                            Starting focus pos
      --debugplots          show debug plots
      --debugplotsdelay DEBUGPLOTSDELAY
                            Delay (seconds) showing each plot
      --simul               Simulate star
      --stayopen            stay open when done
      --profile PROFILE     Name of equipment profile
      --focuser FOCUSER     Focuser Driver
      --camera CAMERA       Camera Driver
      --exposure_start EXPOSURE_START
                            Starting exposure value
      --exposure_min EXPOSURE_MIN
                            Minimum exposure value
      --exposure_max EXPOSURE_MAX
                            Maximum exposure value
      --starflux_min STARFLUX_MIN
                            Maximum flux in star
      --saturation SATURATION
                            Saturation level for sensor
      --framesize FRAMESIZE
                            Size of capture frame, 0=full
      --winsize WINSIZE     Size of window used to analyze star
      --focusdelay FOCUSDELAY
                            Delay (seconds) after focus moves
      --numaverage NUMAVERAGE
                            Number of images to average


Explanation of arguments
------------------------

The "--focus_min", "--focus_max", and "--focus_dir" arguments define the allowed region
and direction of the focus run.

The hardware drivers can be specified individually with the "--focuser" and "--camera"
arguments, or pulled from an "astroprofile" using the "--profile" argument.

The "--simul" argument will run the program using a simulated star instead of
connecting to a real camera and focuser.  Useful for testing.




