

def convert_data(file_path):
    with open(file_path, 'r') as f:
        data = [float(line) for line in f.readlines()]

    return data


files = ['split_debug2', 'split_debug1', 'sound_vel', 'temperature', 'density', 'dp_de', 'dt_de', 'dt_drho', 'vof_new', 'sie', 'cell_mass']


for file in files:
    file_path1 = f'/home/talkad_k/Desktop/ScalSALE_DS/src/Scripts/material_results/{file}.txt'
    file_path2 = f'/home/talkad_k/Desktop/ScalSALE_original/ScalSALE_DS/src/Scripts/material_results/{file}.txt'

    data1, data2 = convert_data(file_path1), convert_data(file_path2)

    for val1, val2 in zip(data1, data2):
        if val1 != val2:
            print(f"You piece of shit - {file}")
            break

print("MY GUY")
