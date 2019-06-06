Using autofocus_auto_star.py
============================

Introduction
------------

The script autofocus_auto_star.py handles finding a focus star, slewing to the
star, autofocusing and then returning to the original position.  It does this by
calling other python utilities which handle most of the actual work.

Invocation
----------

The invocation of autofocus_auto_star.py is:

.. code-block:: bash

    usage: autofocus_auto_star.py [-h] [--profile PROFILE] [--lst LST]
                                  [--onlyside ONLYSIDE] [--lon LON]
                                  [--meridianthres MERIDIANTHRES]
                                  [--maxtries MAXTRIES]
                                  dist mag
    positional arguments:
      dist                  Max distance in degrees
      mag                   Desired mag focus star

    optional arguments:
      -h, --help            show this help message and exit
      --profile PROFILE     Name of astro profile
      --lst LST             Local sidereal time
      --onlyside ONLYSIDE   EAST or WEST side only
      --lon LON             Location longitude
      --meridianthres MERIDIANTHRES
                            How close to meridian is allowed (hh:mm:ss)
      --maxtries MAXTRIES   Number of stars to try before giving up

Explanation of specifying side of pier
--------------------------------------

The "--lon" argument allows the specification of the observing latitude.  Then
script can then compute the local sidereal time.  Optionally the local sidereal
time can be given with the "--lst" argument.

Once the local sidereal time has been determined then the "--onlyside" parameter
can be used to retrict the star to one side of the meridian or the other.  It can
take a value of "EAST" or "WEST" (capitalized!).

The "--meridianthres" argument can be used to create a "keep out" area near the
meridian that excludes choosing a focus star in that area.

Using an astro profile
----------------------

There are no specific settings covered by an astro profile for the autofocus_auto_star.py
script, but several scripts it relies on do.


