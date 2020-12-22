Superquadric Collision Detection (SCD)
====
SCD implements the method of (Ruan, 2019) for determining if two objects are in collision. The first object is a convex superquadric and the second object is an ellipsoid. Currently 2D and 3D cases are supported. The method depends on solving for the root of a polynomial. We use Ceres to find this root by minimizing the square of the polynomial. The method returns the point along a line that is the boundary between contact and no contact. By checking the magnitude of the ellipsoid's position relative to the superquadric, we can detect collision. See the paper for more details.
![](/doc/box_ellipse_2d.png)
![](/doc/cube_ellipsoid.png)

## Dependencies
### C++
* Eigen
* Ceres: Polynomial root-finding
* GoogleTest
### MATLAB / Python
* CasADi: Polynomial root-finding with IPOPT
### Notes
* The interface is inspired by that of the [Flexible Collision Library](https://github.com/flexible-collision-library/fcl).
### References

* [S. Ruan, K. L. Poblete, Y. Li, Q. Lin, Q. Ma and G. S. Chirikjian, "Efficient Exact Collision Detection between Ellipsoids and Superquadrics via Closed-form Minkowski Sums," 2019 International Conference on Robotics and Automation (ICRA), Montreal, QC, Canada, 2019, pp. 1765-1771, doi: 10.1109/ICRA.2019.8793496.](https://ieeexplore.ieee.org/abstract/document/8793496)
```
