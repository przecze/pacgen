PackGen
=======

This library is capable of producing packings of spheres that follow a statistical distribution regarding their diameters filling arbitrary domains defined by triangular meshes.

----------

## Default Containers
### Box
Defined by two points

 - bmin: The minimum corner
 - bmax: The maximum corner

``` C++
 PG::Container* container = new PG::Box(-1.0, -1.0, -1.0,
                                         1.0,  1.0,  1.0);
```

### Cylinder
Defined by two points and a scalar

 - p0: The bottom position
 - p1: The top position
 - rad: The radius

``` C++
 PG::Container* container = new PG::Cylinder(0.0, 0.0, 0.0,
                                             0.0, 5.0, 0.0,
                                             1.0);

```

----------


## Default Radii Generators
### Unifom Distribution
Has 2 parameters:

 - rmin: minimum radius
 -  rmax: maximum radius.
``` C++
PG::NG* ng = new PG::UniformNG(rmin, rmax);
```
### Bernoulli Distribution
Has 3 parameters:

 - rmin: minimum radius.
 - rmax: maximum radius.
 - rmaxprob: the probability to create spheres with the maximum radius.
``` C++
PG::NG* ng = new PG::BernoulliNG(rmin, rmax, rmaxprob);
```
### Truncated Gaussian Distribution
Has 4 parameters:

 - rmin: minimum radius.
 - rmax: maximum radius.
 - mean: mean or expectation of the distribution.
 - stdv: standard deviation.
``` C++
PG::NG* ng = new PG::GaussianNG(rmin, rmax, mean, stdv);
```

