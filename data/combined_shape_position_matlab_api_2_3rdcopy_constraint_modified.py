import matlab.engine

import random
import numpy as np

import math
import matplotlib
import matplotlib.pyplot as plt
import cma
import os
from scipy.io import savemat


# Define the parameters
prefix='pmcopy3'
parameters = {
    'x21': 0.24*2-0.2, 'x22': 0.2*2-0.2, 'x23': 0.14*2-0.2,
    'x24': 0.08*2-0.2, 'x25': 0.04*2-0.2, 'x26': 0.03*2-0.2,
    'x27': 0.03*2-0.2, 'x28': 0.05*2-0.2, 'x29': 0.08*2-0.2,
    'x210': 0.12*2-0.2, 'x211': 0.16*2-0.2, 'x212': 0.2*2-0.2,
    'y21': 0+0.2+0.5, 'y22': 0.02*2+0.2+0.5, 'y23': 0.02*2+0.2+0.5,
    'y24': 0.02*2+0.2+0.5, 'y25': 0.02*2+0.2+0.5, 'y26': 0.01*2+0.2+0.5,
    'y27': 0+0.2+0.5, 'y28': -0.01*2+0.2+0.5, 'y29': -0.01*2+0.2+0.5,
    'y210': -0.01*2+0.2+0.5, 'y211': -0.01*2+0.2+0.5, 'y212': -0.01*2+0.2+0.5,
    'x11': 0.24*2-0.2, 'x12': 0.2*2-0.2, 'x13': 0.14*2-0.2,
    'x14': 0.08*2-0.2, 'x15': 0.04*2-0.2, 'x16': 0.03*2-0.2,
    'x17': 0.03*2-0.2, 'x18': 0.05*2-0.2, 'x19': 0.08*2-0.2,
    'x110': 0.12*2-0.2, 'x111': 0.16*2-0.2, 'x112': 0.2*2-0.2,
    'y11': 0+0.2+0.1, 'y12': 0.02*2+0.2+0.1, 'y13': 0.02*2+0.2+0.1,
    'y14': 0.02*2+0.2+0.1, 'y15': 0.02*2+0.2+0.1, 'y16': 0.01*2+0.2+0.1,
    'y17': 0+0.2+0.1, 'y18': -0.01*2+0.2+0.1, 'y19': -0.01*2+0.2+0.1,
    'y110': -0.01*2+0.2+0.1, 'y111': -0.01*2+0.2+0.1, 'y112': -0.01*2+0.2+0.1,
}
# Save to a .mat file
savemat(f'{prefix}_parameters.mat', parameters)
iteration_counter = 0
with open(f'{prefix}_iteration_count.txt', 'w') as f:
    f.write(str(iteration_counter))

eng = matlab.engine.start_matlab()
matlab_script_path = f'{prefix}_Initilization_Matlab_Comsol.m'
eng.run(matlab_script_path, nargout=0)
matlab_script_path = f'{prefix}_Final_dimen_complete_two_shape_uniform_heat_source.m'
eng.run(matlab_script_path, nargout=0)



num_fins=2
x0=float(eng.workspace['x0'])
fluid_solid_area_width=float(eng.workspace['fluid_solid_area_width'])
fluid_solid_area_height=float(eng.workspace['fluid_solid_area_height'])
x_center_fluid_solid_area=float(eng.workspace['x_center_fluid_solid_area'])
y_center_fluid_solid_area=float(eng.workspace['y_center_fluid_solid_area'])

max_x_fin1=float(eng.workspace['max_x_fin1'])
max_x_fin2=float(eng.workspace['max_x_fin2'])
max_y_fin1=float(eng.workspace['max_y_fin1'])
max_y_fin2=float(eng.workspace['max_y_fin2'])
max_deformation=np.array(eng.workspace['max_deformation'])
max_deform=float(eng.workspace['max_deform'])
x_range=np.array(eng.workspace['x_range'])
y_range=np.array(eng.workspace['y_range'])
x_cen_fin=np.array(eng.workspace['x_cen_fin'])
y_cen_fin=np.array(eng.workspace['y_cen_fin'])

control_pts_locations={}
for i in range(num_fins):
    for j in range(12):
        control_pts_locations[f"x{i+1}{j+1}"]=parameters['x'+str(i+1)+str(j+1)]
        control_pts_locations[f"y{i+1}{j+1}"]=parameters['y'+str(i+1)+str(j+1)]

print(control_pts_locations)
print(y_range[0,0])
print(f"x_cen_fin1, x_cen_fin2, are {x_cen_fin[0]}, {x_cen_fin[1]}, respectively")
print(f"y_cen_fin1, y_cen_fin2,  are {y_cen_fin[0]}, {y_cen_fin[1]}, respectively")
#print(f"max deformation of each fin is , {max_deformation[0]}, {max_deformation[1]}, respectively")
print(f"max deformation of each fin is , {max_deform}")
#max_deformation[0]=0.1
#max_deformation[1]=0.1

T_inlet=float(eng.workspace['T_inlet'])


num_pts_to_move=4
folder_name=f"{prefix}_fin_positioning_shape_CMAES"
if not(os.path.exists(folder_name)):
    os.mkdir(folder_name)
from pathlib import Path
path=Path(folder_name)
print(path)
optimization_class='ave'


def compute_box_penalty(boundary_pts, box_min_y=0.0, box_max_y=1.0, box_min_x=-0.50, box_max_x=0.50):
    penalty_x_lower = np.maximum(0, box_min_x - boundary_pts[:,0])
    penalty_x_upper = np.maximum(0, boundary_pts[:,0] - box_max_x)
    penalty_y_lower = np.maximum(0, box_min_y - boundary_pts[:,1])
    penalty_y_upper = np.maximum(0, boundary_pts[:,1] - box_max_y)
    total_penalty = np.sum(penalty_x_lower + penalty_x_upper + penalty_y_lower + penalty_y_upper)
    return total_penalty

def objective_function_combined(X):
    """
    X layout:
      - Indices [0 : 2*num_fins) -> fin positions (X and Y shifts)
      - Indices [2*num_fins : end) -> shape deformations for all fins
    """
    global iteration_counter
    iteration_counter += 1  # Increment each time function is called
    with open(f'{prefix}_iteration_count.txt', 'w') as f:
        f.write(str(iteration_counter))
    # ----------------------------------------------------------------------
    # 1) Parse the first 2*num_fins entries for fin positions (XX)
    XX=np.zeros(2*num_fins)
    for i in range(num_fins):
        # Map [0..1] -> [-x_range, +x_range], etc.
        XX[i]            =X[i]-0.50 #2 * x_range[i,0] * X[i]            - x_range[i,0]
        XX[num_fins + i] = X[num_fins + i] #2 * y_range[i,0] * X[num_fins + i] - y_range[i,0]

    # ----------------------------------------------------------------------
    # 2) Parse the shape parameters (everything after the position)
    # ----------------------------------------------------------------------
    pos_shape = 2 * num_fins
    shape_params = X[pos_shape:]  # remainder are shape
    # shape_params has length = num_fins * num_pts_to_move * 3
    deformed_matrices = np.reshape(shape_params, (num_fins, num_pts_to_move, 3))

    radius_mid_pts=0.50
    
    for nfin in range(num_fins):
        deformation_each_fin=np.zeros((num_pts_to_move,3))
        deformation_each_fin=deformed_matrices[nfin,:,:]
        boundary_pts=np.zeros((num_pts_to_move,2))
        edgy=np.zeros((num_pts_to_move,))
        for i in range(num_pts_to_move):
            radius =abs(deformation_each_fin[i,0])*max_deform #max(abs(deformation_each_fin[i,0]),0.5)*max_deformation[nfin]
            dangle = (360.0/float(num_pts_to_move))
            angle  = dangle*float(i)+deformation_each_fin[i,1]*dangle/2.0
            x      = radius*math.cos(math.radians(angle))
            y      = radius*math.sin(math.radians(angle))
            edg    = 0.5+0.5*abs(deformation_each_fin[i,2])

            boundary_pts[i,0] = x
            boundary_pts[i,1] = y
            edgy[i]          = edg
            
        center            = np.mean(boundary_pts, axis=0,keepdims=True) 
        boundary_pts -= center

        boundary_pts[:,0]+=x_center_fluid_solid_area+XX[nfin]
        boundary_pts[:,1]+=y_center_fluid_solid_area+XX[num_fins+nfin]

        boundary_pts,edgy= ccw_sort(boundary_pts,edgy) 

        augmented_boundary_pts = np.vstack([boundary_pts, boundary_pts[0,:]])
        
        # Compute list of cartesian angles from one point to the next
        vector = np.diff(augmented_boundary_pts, axis=0)
        angles = np.arctan2(vector[:,1],vector[:,0])
        wrap   = lambda angle: (angle >= 0.0)*angle + (angle < 0.0)*(angle+2*np.pi) 
        angles = wrap(angles)  

        # Create a second list of angles shifted by one point
        # to compute an average of the two at each control point.
        # This helps smoothing the curve around control points
        angles1 = angles
        angles2 = np.roll(angles,1)

        angles  = edgy*angles1 + (1.0-edgy)*angles2 + (np.abs(angles2-angles1) > np.pi)*np.pi
        
        # Add first angle as last angle to close curve
        angles  = np.append(angles, [angles[0]]) 

        for i in range(0,len(augmented_boundary_pts)-1):
            local_ctrl_pts = generate_ctrl_pts(augmented_boundary_pts[i,:],
                                             augmented_boundary_pts[i+1,:],
                                             angles[i],
                                             angles[i+1],
                                             radius_mid_pts)
            #print(local_ctrl_pts)
            if i==0:
                control_pts=local_ctrl_pts
            else:
                control_pts=np.concatenate((control_pts,local_ctrl_pts),axis=0) 
        control_pts_removed=np.delete(control_pts, [4,8,12,15], axis=0)        
        if nfin==0:
            control_pts_total=control_pts_removed
        else:
            control_pts_total=np.concatenate((control_pts_total,control_pts_removed),axis=0) 
            
    penalty = compute_box_penalty(control_pts_total)
    
    penalty_weight=5e3
    if penalty > 0:
    	objective_shape = 1e4 + penalty_weight * penalty  # or np.inf
    else:
    	objective_shape=compute_objective_shape(control_pts_total)
    
    print(objective_shape) 

    return objective_shape  
    
def ccw_sort(points,edgy):
    transfered_pts=points-np.mean(points,axis=0,keepdims=True)
    angles=(np.arctan2(transfered_pts[:,1],transfered_pts[:,0])+0.02)+(np.arctan2(transfered_pts[:,1],transfered_pts[:,0])<-0.02)*2*math.pi
    args=np.argsort(angles)

    
    points=np.array(points)
    sorted_pts=points[args,:]
    sorted_edgy=edgy[args]
    
    return sorted_pts,sorted_edgy   


def generate_ctrl_pts(p1, p2, angle1, angle2, radius):
        dist = compute_distance(p1, p2)
        if (radius == 'random'):
            radius = 0.707*dist*np.random.uniform(low=0.0, high=1.0)
        else:
            radius = 0.707*dist*radius

            # Create array of control pts for cubic Bezier curve
            # First and last points are given, while the two intermediate
            # points are computed from edge points, angles and radius      ##Why 4 points?
            control_pts      = np.zeros((4,2))
            control_pts[0,:] = p1[:]
            control_pts[3,:] = p2[:]
            control_pts[1,:] = p1 + np.array([radius*np.cos(angle1), radius*np.sin(angle1)])
            control_pts[2,:] = p2 + np.array([radius*np.cos(angle2+np.pi), radius*np.sin(angle2+np.pi)])
        return control_pts


# In[5]:


def compute_distance(p1, p2):
    return (np.sqrt(np.inner(p1-p2,p1-p2)))


def compute_objective_shape(control_points):
    try:
        for i in range(num_fins): 
            for j in range(12):
                parameters[f"x{i+1}{j+1}"]=control_points[i*12+j,0]
                parameters[f"y{i+1}{j+1}"]=control_points[i*12+j,1]

    
        savemat(f'{prefix}_parameters.mat', parameters)
        matlab_script_path = f'{prefix}_Final_dimen_complete_two_shape_uniform_heat_source_2.m'
        eng.workspace['iteration'] = iteration_counter
        eng.run(matlab_script_path, nargout=0)
    
        pressure_drop_base_model=-3.0448
        normalized_T_s_ave_base_model=1.4963
        normalized_T_max_base_model=1.5004
    
        T_s_ave=float(eng.workspace['Ts_ave'])
        normalized_T_s_ave=((T_s_ave-T_inlet)/T_inlet)/normalized_T_s_ave_base_model
        
        T_max=float(eng.workspace['T_max'])
        normalized_T_max=((T_max-T_inlet)/T_inlet)/normalized_T_max_base_model
        
        pressure_drop=float(eng.workspace['pressure_drop'])
        normalized_pressure_drop=pressure_drop/pressure_drop_base_model
        
        #objective_both=((T_s_ave-T_inlet)+(T_max-T_inlet))
        objective_both=normalized_T_s_ave+normalized_T_max
        
        #objective_ave=(T_s_ave-T_inlet)
        objective_ave=normalized_T_s_ave
        
        #objective_max=(T_max-T_inlet)
        objective_max=normalized_T_max
        
        #objective_p_max=(-pressure_drop)*(T_s_ave-T_inlet)
        objective_p_max=-pressure_drop
        
        #objective_p_max_T_max=(-pressure_drop)*((T_max-T_inlet)/100)
        objective_p_max_T_max=normalized_pressure_drop*normalized_T_max

        
        folder_path1 = f"{prefix}/"  # Target folder
        #os.makedirs(folder_path1, exist_ok=True)  # Ensure the folder exists
        filename1 = os.path.join(folder_path1, f"control_points_itr_{iteration_counter}.txt")
        with open(filename1, "w") as file:
            for i in range(num_fins):
                file.write(f"Fin {i+1}:\n")  # Header for each fin
                for j in range(12):
                    x_val = control_points[i * 12 + j, 0]
                    y_val = control_points[i * 12 + j, 1]
                    file.write(f"x{i+1}{j+1}: {x_val}, y{i+1}{j+1}: {y_val}\n")
                file.write("\n")  # Separate each fin with a blank line
                
            file.write(f"T_s_ave:{T_s_ave}\n")  # Header for each fin   
            file.write(f"normalized_T_s_ave:{normalized_T_s_ave}\n")  # Header for each fin 
            
            file.write(f"pressure_drop:{pressure_drop}\n")  # Header for each fin
            file.write(f"normalized_pressure_drop:{normalized_pressure_drop}\n")  # Header for each fin
            
            file.write(f"T_max:{T_max}\n")
            file.write(f"normalized_T_max:{normalized_T_max}\n")
            

            T_cons = 700
            penalty=0
            if (T_s_ave>T_cons):
                penalty=100*(T_s_ave-T_cons)
            
            file.write(f"Penalty:{penalty}\n")           
            file.write(f"Cost_function:{objective_p_max+penalty}\n")
      
        
        return (objective_p_max+penalty)
         
    except Exception as e:
        print("COMSOL error occurred:", e)
        return 1e4  # Return a large penalty value to avoid selecting this solution


# Number of decision variables for positions:
num_pos_params = 2 * num_fins  # e.g. 6 for 3 fins

# Number of decision variables for shape:
num_shape_params = num_fins * num_pts_to_move * 3  # e.g. 3*4*3 = 36

# Total dimension:
dim_combined = num_pos_params + num_shape_params  # 6 + 36 = 42 in this example

# Initial guess: just use 0.5 for everything
initial_mean_combined = [0.5]*dim_combined
initial_std_dev_combined = 0.3

# Bounds: for example, 0.0 to 1.0
lower_bounds_combined = [0.0]*dim_combined
upper_bounds_combined = [1.0]*dim_combined

options_combined = {
    'bounds': [lower_bounds_combined, upper_bounds_combined],
    'maxiter': 500,           # or however many
    'verb_disp': 1,
    'verb_filenameprefix': os.path.join(folder_name, f"Method_{optimization_class}_Combined")
}

best_per_iteration_file = os.path.join(folder_name, 'best_per_iteration.txt')
overall_best_file = os.path.join(folder_name, 'overall_best_so_far.txt')
overall_best_value = [float('inf')]
overall_best_cma_iter = [-1]




def cma_callback(es):
    cma_iter = es.countiter            # CMA-ES internal iteration number

    
    # Current population's best value
    curr_gen_fvals = es.fit.fit  # Fitness values of current generation
    curr_gen_best_value = min(curr_gen_fvals)
    
    # Overall best so far
    overall_best_value_now = es.best.f
    
    # Save to per-iteration file (all three: cma_iter, file_iter, curr_gen_best_value, overall_best_value_now)
    with open(best_per_iteration_file, 'a') as f:
        f.write(f"{cma_iter}\t{curr_gen_best_value}\t{overall_best_value_now}\n")
    
    # Update and save overall best file only if improved
    if overall_best_value_now < overall_best_value[0]:
        overall_best_value[0] = overall_best_value_now
        overall_best_cma_iter[0] = cma_iter
        with open(overall_best_file, 'w') as f:
            f.write(f"{overall_best_cma_iter[0]}\t{overall_best_value[0]}\n")
            
            
es_combined = cma.CMAEvolutionStrategy(initial_mean_combined, initial_std_dev_combined, options_combined)

# Optimize:

es_combined.optimize(objective_function_combined, callback=[cma_callback])
result_combined = es_combined.result
optimal_solution_combined = result_combined[0]
print("Optimal X = ", optimal_solution_combined)          





    
