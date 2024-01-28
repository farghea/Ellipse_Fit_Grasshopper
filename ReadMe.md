
## Features

- **3D Ellipse Fitting in Grasshopper**: Fits ellipses to 3D points within Rhino Grasshopper (which does not have it as a built in block), using Python 2 or Python 3 WITHOUT requiring numpy or scipy.

- **No External Dependencies**: Functions independently with built-in Python libraries, ensuring easy setup and compatibility.

- **Custom Linear Algebra Implementation**: Includes a self-contained set of linear algebra tools necessary for ellipse fitting.

- **Grasshopper Integration**: Integrates with Grasshopper, enhancing computational design with geometric fitting capabilities.

## Methods

### Plane Fitting via Pseudo-Inverse

The process begins by fitting a plane to the provided set of 3D points through the minimization of the L2 normby using the pseudo-inverse of the matrix of points.

### Ellipse Fitting in Local 2D Coordinates

Once the plane is fitted, the set of 3D points is projected onto this plane, transforming the problem into a 2D space. In the local 2D coordinate system, we search for the coefficients of the general ellipse equation:

$$
a x^2 + b xy + c y^2 + d x + e y + g = 0
$$

The coefficients are optimized to minimize the L2 norm, similar to the plane fitting procedure. This optimization is carried out through the computation of the pseudo-inverse of the design matrix formed by the projected points.

### Pseudo-Inverse Calculation via Drazin Inverse Method

The pseudo-inverse is computed iteratively using the Drazin inverse method. This method is particularly chosen for its iterative nature, and simplicity:

$$
A_{i+1} = 2A_i - A_i.A.A_i
$$

where $A_i$ represents the approximation of the pseudo-inverse in the i-th iteration, and $A$ is the matrix for which we seek the pseudo-inverse. The recursion proceeds until the changes in $A_i$ fall below a specified threshold, indicating that the pseudo-inverse has been sufficiently approximated.

### Implementation

The implementation of foregoing concepts is within the `MatrixOperations` and `EllipseFitter` classes. The `MatrixOperations` class contains methods for linear algebra, crucial for the computation of the pseudo-inverse and the transformations required for ellipse fitting. The `EllipseFitter` class leverages these matrix operations to compute the plane and ellipse parameters, providing a seamless interface for fitting ellipses to 3D data points within Grasshopper.

## Requirements

No external libraries are required for the main package. However, the "main_example.py" scripts require `numpy` and `matplotlib` to run.

## How to Use

### For Rhino Grasshopper 
Go into the `rhino_GH_files/` directory to find the files:
- `example.gh_ellipse_fit.gh`: This is your Grasshopper file, all set up.
- `example.3dm`: The Rhino model, for visual folks.
- `Grasshopper_ellipse_fit_code.py`: Just the script. It's already baked into the `.gh` file, so no need to touch this.

The `.gh` file has everything you need.

### Want to Tinker with the Source Code?
If you are in a place without numpy or scipy, check out the `src/` folder:
- `geometry_operations.py`: The heart of the operation with all the custom geometric math.

To see it in action, run `main_example.py`.
