# pgee.mixed 0.1.0.9003
* feature/sigma-glm
* sigma is estimated once from glm estimates. 
* Regardless of initial values provided, alpha is estimated from glm estimates.

# pgee.mixed 0.1.0.9002
* fixed bug in CppM and CppM2. (w needs to be squared)

# pgee.mixed 0.1.0.9001

* Added optional ridge component.
* Fixed bugs in examples
* Fixed bug in init\_reset for get\_lambda\_cv\_mixed.

# pgee.mixed 0.1.0.9000

* changed default option for argument `lambda` to handle the default case for `family = "Mixed"`.
* Set `drop=FALSE` in `weighted.scale` to account for 2 column design matrix.

# pgee.mixed 0.1.0

First release.
