Run with containers
====

Docker
--------

* Prepare input files as mentioned earlier.  

* Install docker on your system. `Docker Installation Instructions <https://docs.docker.com/engine/installation/>`_

* Fetch and use the needed image from dockerhub. (https://hub.docker.com/r/likelet/cirpipe/)

* Run the circPipe. Run the pipeline with the option ``-profile standard,docker``.

Singularity
--------

-  Prepare input files as mentioned earlier.
-  Run the circPipe. Run the pipeline with the option
   ``-profile standard,singularity``.

   .. note:: An image containing all of the software requirements will be automatically fetched and used from singularity hub.



-  Run offline with Singularity Download and transfer the Singularity
   image first:

   .. code:: bash

       singularity pull --name circPipe.simg shub://likelet/cirpipe

   Use ``-with-singularity`` and specify the path to the image file:

   .. code:: bash

       nextflow run /path/to/circPipe -with-singularity circPipe.simg

 .. note:: Remember to pull updated versions of the singularity image when updating the pipeline.

 
