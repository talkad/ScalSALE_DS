

def convert_data(file_path):
    with open(file_path, 'r') as f:
        data = [float(line) for line in f.readlines()]

    return data


file_path1 = '/home/talkad_k/Desktop/ScalSALE_DS/src/Scripts/material_results/pressure.txt'
file_path2 = '/home/talkad_k/Desktop/ScalSALE_original/ScalSALE_DS/src/Scripts/material_results/pressure.txt'

data1, data2 = convert_data(file_path1), convert_data(file_path2)

print(f"data len {len(data1)}")

for val1, val2 in zip(data1, data2):
    if val1 != val2:
        print("You piss of shit")

print("MY GUY")
