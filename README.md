**3denseSTORM** is a MATLAB program for analysis of 2D and 3D STORM data with high molecular density and it is based on the algorithm described in our [paper](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-22-25-31263):

* M. Ovesný, P. Křížek, Z. Švindrych, G. M. Hagen. *High density 3D localization microscopy using sparse support recovery*. Opt. Express 22, 31263-31276 (2014)

Consult the examples to see how to use this program. 3denseSTORM is compatible with [ThunderSTORM](https://github.com/zitmen/thunderstorm/) as it stores results into a CSV file. Also an analytic calibration of defocus can be adopted.

----

My Ph.D. thesis [Computational methods in single molecule localization microscopy](https://www.researchgate.net/publication/311426573_Computational_methods_in_single_molecule_localization_microscopy) now available. It contains detailed description of the algorithm and it's implementation.

----

**Feature to add**

* Mapping of channels in biplane mode through a transformation matrix
* Automated calibration:
    * Homotography estimator for biplane
    * Estimation of real psf from a calibration stack
