import numpy as np
import matplotlib.pyplot as plt

# defining plate dimensions:
Lx, Ly = 0.3, 0.4  

# defining griding values of the plate:
Nx, Ny = 31, 41  
dx, dy = Lx / (Nx - 1), Ly / (Ny - 1)
beta = dx / dy
points_count = (Nx-2)*(Ny-2)

# boundry conditions:
T_top = 20.0
T_bottom = 50.0
T_left = 10.0
T_right = 10.0

# defining the main matrix and then giving the first & last row/ fisrt & last column constant values
T = np.zeros((Ny, Nx))
T[0, :] = T_top
T[-1, :] = T_bottom
T[:, 0] = T_left
T[:, -1] = T_right

# assuming a value for all of the indexes (exept first & last column/row) in order to reach answer quicker
T[1:-1, 1:-1] = (T_top + T_bottom + T_left + T_right) / 4

# Relaxation parameters
relaxation_factor = 1.5       
tolerance = 0.0000001       
max_iter = 10000  

# Gauss-Seidel with the use of relaxation factor
for iteration in range(1, max_iter+1):
    old_index_collector = []
    new_index_collector = []
    max_diff = 0.0
    for i in range(1, Ny - 1):
        for j in range(1, Nx - 1):
            T_old = T[i, j]
            T[i,j] = (1 - relaxation_factor) * T[i,j] + (relaxation_factor / (2 * (1 + beta**2))) * (T[i+1, j] + T[i-1, j] + beta**2 * (T[i, j+1] + T[i, j-1]))
            old_index_collector.append(T_old)
            new_index_collector.append(T[i,j])
            
    tolerance_collector = []
    tolerance_flag = True
    for i in range(points_count):
        tolerance_collector.append(abs(old_index_collector[i]-new_index_collector[i]))
    for i in tolerance_collector:
        if i >= tolerance:
            tolerance_flag = False
    
    if tolerance_flag == True:
        print(f"we reach tolerence in {iteration}th iteration")
        print(T.max(),T.min(),T.mean())        
        with open("point-by-point results.txt" , "w") as file:
            file.write("Result:" + "\n")
            for row in T:
                file.write(" ".join([f"{val:.4f}" for val in row]) + "\n")
        x = np.arange(Nx)
        y = np.arange(Ny)
        X, Y = np.meshgrid(x, y)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, T, cmap='hot')
        ax.set_title("dimension of 30cm*40cm - 1cm intervals - point by point Gauss-Seidel")
        ax.set_xlabel("X grid in cm")
        ax.set_ylabel("Y grid in cm")
        ax.set_zlabel("Temperature (°C)")
        plt.show()
        """plt.imshow(T, cmap='hot', origin='lower')
        plt.colorbar(label='Temperature (°C)')
        plt.title('dimension of 30cm*40cm - 1cm intervals - point by point Gauss-Seidel')
        plt.xlabel('X axix in CM')
        plt.ylabel('Y axis in CM')
        plt.show()"""
        break
    if tolerance_flag == False:
        print(f"In {iteration} iteration we didnt reach tolerance")

    