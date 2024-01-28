import copy 
import math 
import ghpythonlib.components as ghc


class MatrixOperations:
    """
    A class for matrix operations: In grasshopper you can't install 
    numpy and other libraries easily, so I had to write this class 
    and half of the linear algebra for this!
    """

    @staticmethod
    def zeros(m, n):
        """
        Generates a matrix of size m x n filled with zeros.
        """
        return [[0.0] * n for _ in range(m)]

    @staticmethod
    def ones(m, n):
        """
        Creates a matrix of size m x n filled with ones.
        """
        return [[1.0] * n for _ in range(m)]

    @staticmethod
    def eye(m, n):
        """
        Produces an m x n identity matrix (1's on the diagonal and 0's elsewhere).
        """
        result = MatrixOperations.zeros(m, n)
        min_index = min(m, n)
        for i in range(min_index):
            result[i][i] = 1.0
        return result

    @staticmethod
    def size(matrix):
        """
        Returns the size of the given matrix as a tuple (rows, columns).
        """
        return len(matrix), len(matrix[0])

    @staticmethod
    def add(mat_a, mat_b):
        """
        Adds two matrices mat_a and mat_b of the same size.
        """
        m, n = MatrixOperations.size(mat_a)
        result = MatrixOperations.zeros(m, n)
        for i in range(m):
            for j in range(n):
                result[i][j] = mat_a[i][j] + mat_b[i][j]
        return result

    @staticmethod
    def transpose(matrix):
        """
        Transposes the given matrix (flips rows and columns).
        """
        m, n = MatrixOperations.size(matrix)
        return [[matrix[j][i] for j in range(m)] for i in range(n)]

    @staticmethod
    def multiply(mat_a, mat_b):
        """
        Multiplies two matrices mat_a and mat_b.
        """
        mA, nA = MatrixOperations.size(mat_a)
        mB, nB = MatrixOperations.size(mat_b)
        result = MatrixOperations.zeros(mA, nB)
        for i in range(mA):
            for j in range(nB):
                for k in range(nA):
                    result[i][j] += mat_a[i][k] * mat_b[k][j]
        return result

    @staticmethod
    def scalar_multiply(matrix, num):
        """
        Multiplies every element in matrix by num.
        """
        m, n = MatrixOperations.size(matrix)
        result = MatrixOperations.ones(m, n)
        for i in range(m):
            for j in range(n):
                result[i][j] *= matrix[i][j] * num
        return result

    @staticmethod
    def sum(matrix):
        """
        Sums up all elements in matrix, returning a matrix of sums for each row.
        """
        m, n = MatrixOperations.size(matrix)
        result = MatrixOperations.zeros(m, 1)
        for i in range(m):
            for j in range(n):
                result[i][0] += matrix[i][j]
        return result

    @staticmethod
    def max(matrix):
        """
        Finds the maximum value in matrix along with its position.
        """
        m, n = MatrixOperations.size(matrix)
        max_value = matrix[0][0]
        max_index = [0, 0]
        for i in range(m):
            for j in range(n):
                if matrix[i][j] > max_value:
                    max_value = matrix[i][j]
                    max_index = [i, j]
        return max_value, max_index

    @staticmethod
    def abs(matrix):
        """
        Applies the absolute value function to each element in matrix.
        """
        m, n = MatrixOperations.size(matrix)
        result = MatrixOperations.zeros(m, n)
        for i in range(m):
            for j in range(n):
                result[i][j] = abs(matrix[i][j])
        return result

    @staticmethod
    def pinv(matrix):
        """
        Computes the pseudo-inverse of matrix.
        - I used an iterative method
        """
        m, n = MatrixOperations.size(matrix)
        A_i_1 = MatrixOperations.transpose(copy.deepcopy(matrix))
        A_i_1 = MatrixOperations.scalar_multiply(A_i_1, 1e-10)
        for _ in range(3000):
            term1 = MatrixOperations.scalar_multiply(A_i_1, 2.0)
            term2 = MatrixOperations.multiply(matrix, A_i_1)
            term3 = MatrixOperations.multiply(A_i_1, term2)
            term4 = MatrixOperations.scalar_multiply(term3, -1.0)
            A_i = MatrixOperations.add(term1, term4)
            A_i_1 = A_i

            B = MatrixOperations.multiply(matrix, A_i)
            my_eye = MatrixOperations.scalar_multiply(MatrixOperations.eye(m, m), -1.0)
            B = MatrixOperations.add(B, my_eye)
            B = MatrixOperations.abs(B)
            max_val, _ = MatrixOperations.max(B)
            if max_val < 1e-10:
                break
        return A_i

    @staticmethod
    def mean(matrix):
        """
        Calculates the mean of each column in matrix.
        """
        m, n = MatrixOperations.size(matrix)
        my_mean = MatrixOperations.transpose(MatrixOperations.sum(MatrixOperations.transpose(matrix)))
        return MatrixOperations.scalar_multiply(my_mean, 1 / m)

    @staticmethod
    def cross(a, b):
        """
        Calculate the cross product of two 3D vectors.

        Args:
            a (list): The first vector, a 3-element list.
            b (list): The second vector, a 3-element list.

        Returns:
            list: The cross product of vectors a and b.
        """
        result = MatrixOperations.zeros(3, 1)
        result[0][0] = a[1][0] * b[2][0] - a[2][0] * b[1][0]
        result[1][0] = -(a[0][0] * b[2][0] - a[2][0] * b[0][0])
        result[2][0] = a[0][0] * b[1][0] - a[1][0] * b[0][0]
        return result
    


# - - - - - - - - - - - - 
class EllipseFitter:
    """
    A class for fitting an ellipse to a set of 3D points.
    """

    def __init__(self, points):
        """
        Initialize the EllipseFitter with a set of points.

        Args:
            points (list): A list of 3D points.
        """
        self.points = points

    def fit_ellipse(self):
        """
        Fits an ellipse to the provided 3D points.

        Returns:
            dict: A dictionary containing ellipse parameters and fitted points.
        """
        mo = MatrixOperations

        # Compute pseudo-inverse of the poin`ts
        pseudo_inv = mo.pinv(self.points)

        # Compute coefficients of the plane
        _, no_of_rows = mo.size(pseudo_inv)
        coef_plane = mo.multiply(pseudo_inv, mo.ones(no_of_rows, 1))
        center = mo.mean(self.points)

        # Orthonormal coordinate system
        e3 = coef_plane
        norm_e3 = math.sqrt(mo.multiply(mo.transpose(e3), e3)[0][0])
        e3 = mo.scalar_multiply(e3, 1 / norm_e3)

        e1 = mo.zeros(3, 1)
        for i in range(3):
            e1[i][0] = self.points[0][i] - center[0][i]
        norm_e1 = math.sqrt(mo.multiply(mo.transpose(e1), e1)[0][0])
        e1 = mo.scalar_multiply(e1, 1 / norm_e1)

        e2 = mo.cross(e3, e1)

        # Transformation matrix
        trans = mo.zeros(3, 3)
        for i in range(3):
            trans[i][0] = e1[i][0]
            trans[i][1] = e2[i][0]
            trans[i][2] = e3[i][0]

        self._transformation_matrix = trans
        # Projection in the new coordinate system
        local_proj = mo.multiply(self.points, trans)
        self.local_proj = local_proj

        # Ellipse fitting
        ellipse_params = self._fit_ellipse_to_projection(local_proj)

        return ellipse_params

    def _fit_ellipse_to_projection(self, proj):
        """
        Fit an ellipse to the projected points.

        Args:
            proj (list): Projected points.

        Returns:
            dict: Parameters of the fitted ellipse.
        """
        mo = MatrixOperations
        m = len(proj)
        x, y = self._extract_coordinates(proj, m)
        mean_xyz = mo.mean(proj)
        x, y = self._center_coordinates(x, y, mean_xyz)

        my_data = self._form_ellipse_equation_data(x, y, m)
        term1 = mo.transpose(mo.sum(mo.transpose(my_data)))
        term2 = mo.multiply(mo.transpose(my_data), my_data)
        term3 = mo.pinv(term2)

        coef_ellipse = mo.multiply(term1, term3)
        ellipse_params = self._compute_ellipse_parameters(coef_ellipse, mean_xyz, proj)

        return ellipse_params
    
    def _extract_coordinates(self, proj, m):
        """Extracts x and y coordinates from the projected points."""
        x = MatrixOperations.zeros(m, 1)
        y = MatrixOperations.zeros(m, 1)
        for i in range(m):
            x[i][0] = proj[i][0]
            y[i][0] = proj[i][1]
        return x, y

    def _center_coordinates(self, x, y, mean_xyz):
        """Centers the x and y coordinates by subtracting their mean."""
        m = len(x)
        mean_x = mean_xyz[0][0]
        mean_y = mean_xyz[0][1]
        for i in range(m):
            x[i][0] -= mean_x
            y[i][0] -= mean_y
        return x, y

    def _form_ellipse_equation_data(self, x, y, m):
        """Forms the data matrix for the ellipse equation."""
        my_data = MatrixOperations.zeros(m, 5)
        for i in range(m):
            my_data[i][0] = x[i][0] * x[i][0]
            my_data[i][1] = x[i][0] * y[i][0]
            my_data[i][2] = y[i][0] * y[i][0]
            my_data[i][3] = x[i][0]
            my_data[i][4] = y[i][0]
        return my_data

    def _compute_ellipse_parameters(self, coef_ellipse, mean_xyz, local_proj):
        """Computes the parameters of the ellipse from its coefficients."""
        # Extracting coefficients
        a0, b0, c0, d0, e0 = coef_ellipse[0]
        mean_x0 = mean_xyz[0][0]
        mean_y0 = mean_xyz[0][1]

        a = copy.deepcopy(a0)
        b = copy.deepcopy(b0)
        c = copy.deepcopy(c0)
        d = copy.deepcopy(d0)
        e = copy.deepcopy(e0)
        
        if min(abs(b/a), abs(b/c)) > 0.001:
            orientation_rad = 0.5 * math.atan(b/(c-a))
            cos_phi = math.cos( orientation_rad )
            sin_phi = math.sin( orientation_rad )

            a = a0*cos_phi**2 - b0*cos_phi*sin_phi + c0*sin_phi**2
            b = 0
            c = a0*sin_phi**2 + b0*cos_phi*sin_phi + c0*cos_phi**2
            d = d0*cos_phi - e0*sin_phi
            e = d0*sin_phi + e0*cos_phi

            mean_x = cos_phi*mean_x0 - sin_phi*mean_y0
            mean_y = sin_phi*mean_x0 + cos_phi*mean_y0
        else:
            orientation_rad = 0
            cos_phi = math.cos( orientation_rad )
            sin_phi = math.sin( orientation_rad )
        
        if (a < 0):
            a = -a
            c = -c
            d = -d
            e = -e

        # Center in the ellipse coordiante system (within the local coordinate system)
        X0 = mean_x - d/2/a
        Y0 = mean_y - e/2/c
        F = 1 + (d*d)/(4*a) + (e*e)/(4*c)
        # Maximum dimeter
        a = ( F/a )**0.5
        # Minimum dimeter
        b = ( F/c )**0.5
        
        long_axis   = 2*max(a,b)
        short_axis  = 2*min(a,b)

        Z0 = local_proj[0][2]

        center0 = self._rotate_center_and_ellipse(orientation_rad, X0, Y0, Z0, a, b, 25)
        
        return {
            "center": center0,
            "axes": [long_axis, short_axis],
            "orientation": orientation_rad
        }
    
    def _rotate_center_and_ellipse(self, orientation_rad, X0, Y0, Z0, a, b, n = 25):
        """
        Constructs ellipse points and rotates them back to the original coordinate system.

        Args:
            orientation_rad (float): orientation of ellipse in radian
            X0, Y0 (float, float): center of ellipse

        Returns:
            list: center of ellipse in local coordinate 
        """
        mo = MatrixOperations
        cos_phi = math.cos(orientation_rad)
        sin_phi = math.sin(orientation_rad)

        # Rotation matrix for the local coordinate axis
        myrotation_matrix = mo.zeros(2, 2)
        myrotation_matrix[0][0], myrotation_matrix[0][1] = cos_phi, sin_phi
        myrotation_matrix[1][0], myrotation_matrix[1][1] = -sin_phi, cos_phi

        # Center of the ellipse in the local coordinate system
        temp = mo.zeros(2, 1)
        temp[0][0] = X0
        temp[1][0] = Y0
        
        p_in = mo.multiply(myrotation_matrix, temp)
        
        # Constructing ellipse 
        theta_degree_list = [360 * i / (n - 1) for i in range(n)]

        ellipse0 = mo.zeros(n, 2)
        for i, theta_deg in enumerate(theta_degree_list):
            theta = theta_deg*math.pi/180
            x = a*math.cos(theta)
            y = b*math.sin(theta)
            ellipse0[i][0] = x
            ellipse0[i][1] = y

        ellipse0 = mo.transpose(mo.multiply(myrotation_matrix, mo.transpose(ellipse0)))
        # Transformed ellipse in the local coordiante yet in 3D
        ellipse_local_coordinate = mo.zeros(n, 3)
        for i in range(n):
            ellipse_local_coordinate[i][0] = ellipse0[i][0] + p_in[0][0]
            ellipse_local_coordinate[i][1] = ellipse0[i][1] + p_in[1][0]
            ellipse_local_coordinate[i][2] = Z0

        # transforming the points in local coordinate system into the global one
        ellipse_global_coordinate = copy.deepcopy(ellipse_local_coordinate)
        inverse_transformation = mo.pinv(self._transformation_matrix)
        ellipse_global_coordinate = mo.multiply(ellipse_global_coordinate, inverse_transformation)
        
        self.ellipse_local_coordinate = ellipse_local_coordinate
        self.ellipse_global_coordinate = ellipse_global_coordinate

        return p_in


# processing input points to MO proper input
points = MatrixOperations.zeros(len(input_points), 3)
for i, pnt in enumerate(input_points):
    for j in range(3):
        points[i][j] = pnt[j]

# fitting the ellipse
ellipse_fitter = EllipseFitter(points)
ellipse_params = ellipse_fitter.fit_ellipse()

# Assigning max/min diameters
max_diameter, min_diameter = ellipse_params['axes']

# Getting points of the fitted ellipse
ellipse_points = []
for pnt in ellipse_fitter.ellipse_global_coordinate:
    gh_pnt = ghc.ConstructPoint(pnt[0], pnt[1], pnt[2])
    ellipse_points.append(gh_pnt)

# Get ellipse curve
ellipse_crv = ghc.Interpolate(ellipse_points, 3, True, 0)['curve']

