Superquadric Collision Detection (SCD)
====
SCD implements the method of (Ruan, 2019) for determining if two objects are in collision. The first object is a superquadric and the second object is an ellipsoid.

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
