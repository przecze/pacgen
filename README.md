PackGen
=======

This library is capable of producing packings of spheres that follow a statistical distribution regarding their diameters filling arbitrary domains defined by triangular meshes.

![A sphere pack inside the bunny mesh](images/bunny.png?raw=true)

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

### Triangle Mesh
Defined by:

 - v: Array of vertices. (x_1, y_1, z_1, x_2, y_2, z_2, ... , x_nv, y_nv, z_nv)
 - nv: Number of vertices
 - id: Array of triangle indices. (t1_a, t1_b, t1_c, t2_a, t2_b, t2_c, ..., tnid_a, tnid_b, tnid_c)
 - nid: Number of indices.
 - rmax: Maximum radius in the pack.
 - npoints: Number of distance points per dimension. With more resolution the distance computation is more accurate, nevertheless, the memory consumption grows cubically.

``` C++
 void* mesh = PG::Mesh3DCreate(v, nv, id, nid);
 PG::Container* container = new PG::USDF(mesh, rmax, npoints);
```
----------


## Default Radii Generators
### Unifom Distribution
Has 2 parameters:

 - rmin: Minimum radius
 - rmax: Maximum radius.
``` C++
PG::NG* ng = new PG::UniformNG(rmin, rmax);
```
### Bernoulli Distribution
Has 3 parameters:

 - rmin: Minimum radius.
 - rmax: Maximum radius.
 - rmaxprob: The probability to create spheres with the maximum radius.
``` C++
PG::NG* ng = new PG::BernoulliNG(rmin, rmax, rmaxprob);
```
### Truncated Gaussian Distribution
Has 4 parameters:

 - rmin: Minimum radius.
 - rmax: Maximum radius.
 - mean: Mean or expectation of the distribution.
 - stdv: Standard deviation.
``` C++
PG::NG* ng = new PG::GaussianNG(rmin, rmax, mean, stdv);
```

----------
## Demo

### Dragon Pack

[![Dragon Pack](https://i.vimeocdn.com/video/494078915_640.jpg)](https://vimeo.com/109980855 "Dragon Pack")