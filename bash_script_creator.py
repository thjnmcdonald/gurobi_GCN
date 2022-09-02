location_list = ['2022-08-31_12:00:14-GCN_16_1',
 '2022-08-31_12:01:19-GCN_32_1',
 '2022-08-31_12:02:01-GCN_64_1',
 '2022-08-31_12:38:58-GCN_16_2',
 '2022-08-31_12:04:42-GCN_32_2',
 '2022-08-31_12:45:57-GCN_64_2',
 '2022-08-31_12:06:21-GCN_16_3',
 '2022-08-31_12:57:18-GCN_32_3',
 '2022-08-31_13:02:07-GCN_64_3']

def create_bash_script(location_list):
    print('SECONDS=0')
    mol_len_list = [4,6]
    for j, mol_len in enumerate(mol_len_list):
        for i, location in enumerate(location_list):
            print(f'python GCN_MILP_14_features.py --location {location} --time_lim 10 --mol_len {mol_len} > GCN_outputs/{location}_mol_len_{mol_len}.txt')
            print(f'echo \"progress: {(j)*len(location_list) + i + 1}/{(len(location_list))*(len(mol_len_list))} after: $SECONDS s\"')

create_bash_script(location_list)

# run as
# python bash_script_creator.py > bash.script