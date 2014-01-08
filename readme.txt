
--------------------------------------------------------
Note
--------------------------------------------------------

Further development has been moved to the main Boost git repository.


--------------------------------------------------------
General Information
--------------------------------------------------------

This is a partial Boost.Math area for:

* Bernoulli numbers and applications to special functions
* This work has been developed for the 2013 Google Summer of Code

--------------------------------------------------------
Author information
--------------------------------------------------------

Project Developer: Nikhar Agrawal
Project Mentor   : Christopher Kormanyos
Project Advisors : Paul Bristow and John Maddock

--------------------------------------------------------
Project details
--------------------------------------------------------

This project adds several new files to the existing file structure
of Boost.Math. These new files contain the implementations of
certain new functions. These include Bernoulli numbers, extension
of tgamma and lgamma to high precision, and polygamma.

Upon successful completion of this project, the candidate files
in this project may potentially be merged to and/or included
in Boost.Math.

--------------------------------------------------------
Using this project with Boost
--------------------------------------------------------

The files in this project can be used with an existing Boost trunk
or with a distribution of Boost. Do note, however, that at the time
this project is being written, the files require the trunk in order
to compile.

The way to use these files with an existing Boost trunk or distro
is to include the base directory of this project in the "C++ way".
The added include path of this project must be *upstream* of the
include path of the Boost trunk or distro.
