#%% For this package you do NOT need any requirement, but for this example you need numpy and matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.geometry_operations import EllipseFitter


if __name__ == '__main__':
    points = [[-13.3935141446573,25.2621310702183,18.5760989374103],[-14.4387801675304,25.2489517822846,17.0836365401485],[-14.647783096337,23.8109853435792,15.9833131513718],[-13.925856833847,22.0608101537935,16.0312135191575],[-12.9163168254492,20.5328639050852,16.6119081571289],[-11.9506126796273,20.6880875388436,18.0706886649447],[-11.8012456010447,22.0866377859798,19.0642609119863],[-12.31827673792,23.898941251315,19.3421829987595]]
    ellipse_fitter = EllipseFitter(points)
    ellipse_params = ellipse_fitter.fit_ellipse()

    print(ellipse_params)

    np_ellipse_local_coordinate = np.array(ellipse_fitter.ellipse_local_coordinate)
    np_local_proj = np.array(ellipse_fitter.local_proj)
    np_points = np.array(points)
    np_ellipse_global = np.array(ellipse_fitter.ellipse_global_coordinate)
    

    fig = plt.figure(figsize=(18/2.54, 8/2.54))
    # Plotting the 2D plot in a "local" coordinate system
    ax1 = fig.add_subplot(1, 2, 1) 
    ax1.plot(np_local_proj[:, 0], np_local_proj[:, 1], '.', markersize=10)
    ax1.plot(np_ellipse_local_coordinate[:, 0], np_ellipse_local_coordinate[:, 1], '-r')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title('Local Coordinate')
    
    # Creating a 3D subplot in the "global" coordinate system
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.scatter(np_points[:, 0], np_points[:, 1], np_points[:, 2], color='b', marker='o', label='Original Points')
    ax2.plot(np_ellipse_global[:, 0], np_ellipse_global[:, 1], np_ellipse_global[:, 2], '-r', label='Fitted Ellipse')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.legend()

    plt.tight_layout()
    plt.show()






 