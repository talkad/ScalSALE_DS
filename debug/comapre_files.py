

def convert_data(file_path):
    with open(file_path, 'r') as f:
        # data = [float(line) for line in f.readlines()]
        data = [line for line in f.readlines()]

    return data


files = ['split_debug2', 'split_debug1', 'sound_vel', 'temperature', 'density', 'dp_de', 'dt_de', 'dt_drho', 'vof_new', 'sie', 'cell_mass']


for file in files:
    file_path1 = f'/home/talkad/Desktop/ScalSALE_DS_master/src/Scripts/material_results/{file}.txt'
    file_path2 = f'/home/talkad/Desktop/ScalSALE_DS_original/src/Scripts/material_results/{file}.txt'

    data1, data2 = convert_data(file_path1), convert_data(file_path2)

    for idx, (val1, val2) in enumerate(zip(data1, data2)):
        # print(idx)
        if val1 != val2:
            print(f'idx {idx}:\n    val1 {val1}\n   val2 {val2}')
            print(f"You piece of shit - {file}")
            break

print("MY GUY")
