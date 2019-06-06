Using find_nearby_stars.py
============================

Introduction
------------

The "script find_nearby_stars.py" handles finding a star within a specified
distance from a RA/DEC position with constraints on brightness.  It is normally
used by the "autofocus_auto_star.py" script but can also be invoked indepently.

Invocation
----------

The invocation of autofocus_auto_star.py is:

.. code-block:: bash

    usage: find_nearby_stars.py [-h] [--minmag MINMAG] [--maxmag MAXMAG]
                                [--verbose] [--outfile OUTFILE] [--force]
                                [--lst LST] [--onlyside ONLYSIDE]
                                [--meridianthres MERIDIANTHRES] [--lon LON]
                                cat ra2000 dec2000 dist

    positional arguments:
      cat                   Catalog to search
      ra2000                RA J2000
      dec2000               DEC J2000
      dist                  Max distance in degrees

    optional arguments:
      -h, --help            show this help message and exit
      --minmag MINMAG
      --maxmag MAXMAG
      --verbose
      --outfile OUTFILE     Output file with candidates
      --force               Overwrite output file
      --lst LST             Local sidereal time
      --onlyside ONLYSIDE   EAST or WEST side only
      --meridianthres MERIDIANTHRES
                            How close to meridian is allowed (hh:mm:ss)
      --lon LON             Location longitude

Program Output
--------------
The program outputs the list of candidate stars to the console and if the argument
"--outfile" is given it will also write CSV output to this file.  The file includes a
header that explains the columns.

Explanation of specifying side of pier
--------------------------------------
The "cat" argument should reference a binary SAO Catalog created with the utilities
in the "find_star" directory.  One such file is in the "data" directory and is
called "SAO_Catalog_m5_p11_filtered.bin" and has stars down to magnitude 11.  It
has been filtered of stars that are close to one another to reduce the chance
of having a another star interfere with the autofocus routine.

The "--lon" argument allows the specification of the observing latitude.  Then
script can then compute the local sidereal time.  Optionally the local sidereal
time can be given with the "--lst" argument.

Once the local sidereal time has been determined then the "--onlyside" parameter
can be used to retrict the star to one side of the meridian or the other.  It can
take a value of "EAST" or "WEST" (capitalized!).

The "--meridianthres" argument can be used to create a "keep out" area near the
meridian that excludes choosing a focus star in that area.


