pch 1.1
=============

User-level changes
------------------
* the additional function **splinex** has been implemented.
* the default number and position of the breaks has been changed.
* a special handling of zero-risk regions has been implemented. This may cause differences with version 1.0.
* information on convergence is reported.

Internal changes
--------------------
* improved convergence of the algorithm
 
Changes in the package structure
----------------
* the NAMESPACE file has been modified to import:  
stats: prcomp  
splines: ns bs  
