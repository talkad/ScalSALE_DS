

def convert_data(file_path):
    with open(file_path, 'r') as f:
        # data = [float(line) for line in f.readlines()]
        data = [line for line in f.readlines()]

    return data

# 'pressure2', 'mat2', 'dvel_x_dx', 'dvel_x_dy', 'dvel_x_dz', 'dvel_y_dx', 'dvel_y_dy', 'dvel_z_dx', 'dvel_z_dy', 'dvel_z_dz', 'avi_acceleration_x', 'avi_acceleration_y', 'avi_acceleration_z', 'inversed_vertex_mass2', 'acceleration_z', 'velocity_z', 'a_visc2', 'a_visc4', 'inversed_vertex_mass'

# files = [f'material_results/{file}' for file in ['cell_mass', 'vof_total', 'density', 'dp_de', 'dp_drho', 'dt_de', 'pressure', 'dt_drho', 'sie', 'sound_vel', 'temperature', 'vof_new', 'sie_total', 'pressure_sum_total', 'volume_total', 'cell_mas_total', 'dt_de_total', 'temperature_total', 'dp_drho_total', 'density_total', 'sound_total', 'split_debug1', 'split_debug2', 'pressure_total']]

# files = ['total_pressure_result', 'velocity_result', 'total_cell_mass_result', 'total_sie_result',  'sie', 'density', 'temperature'] + files


files = ["pressure"]

for file in files:
    file_path1 = f'/home/talkad/Desktop/ScalSALE_DS_master/src/Scripts/material_results/{file}.txt' # f'/home/talkad/Desktop/ScalSALE_DS_original/src/Scripts/{file}.txt'
    file_path2 = f'/home/talkad/Desktop/ScalSALE_DS_master/src/Scripts/material_results_full/{file}.txt'

    data1, data2 = convert_data(file_path1), convert_data(file_path2)

    for idx, (val1, val2) in enumerate(zip(data1, data2)):
        # print(idx)
        if val1 != val2:
            print(f'idx {idx}:\n    val1 {val1}\n   val2 {val2}')
            print(f"You piece of shit - {file}")
            break

print("MY GUY")

