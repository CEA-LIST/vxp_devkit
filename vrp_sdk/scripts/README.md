Build and installation script for the MPFR library
==================================================

About
-----

This script download gmp, mpfr tarballs, extract sources, compile them and install library.


Configuration
-------------

This script compile a RISCV (hard-float and soft-float version) and a x86 linux version of the library.

You can change the configuration of the script by editing variable at the begining of the script.

In addition you can tune two parameters from command line by setting one of this variable :

- To change the installation path of the libraries:

``` 
setenv INSTALL_ROOT /home/560.1.361-EPI/MPFR/aar123/scripts
``` 

- To tune the name of the installation directory of linux version by appending a string describing the linux environment used to compile the library.

``` 
setenv LINUX_VERSION "-sl7.9"
``` 


Execution
---------

``` 
bash ./build_mpfr.sh
``` 
