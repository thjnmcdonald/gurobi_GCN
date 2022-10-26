location_list = ['2022-08-31_12:05:50-GCN_16_3',
 '2022-08-31_12:51:55-GCN_16_3',
 '2022-08-31_12:06:21-GCN_16_3',
 '2022-08-31_12:50:09-GCN_16_3',
 '2022-08-31_12:51:30-GCN_16_3']

def create_bash_script(location_list):
    print('SECONDS=0')
    print('echo \"started!\"')
    print(f'now=$(date +"%T")')
    print(f'echo \"Started at : $now\"')

    mol_len_list = [4]
    for j, mol_len in enumerate(mol_len_list):
        for i, location in enumerate(location_list):
            print(f'python GCN_MILP_14_features.py --location {location} --time_lim 36000 --mol_len {mol_len} > GCN_outputs/{location[-8:]}_mol_len_{mol_len}.txt')
            print(f'echo \"progress: {(j)*len(location_list) + i + 1}/{(len(location_list))*(len(mol_len_list))} after: $SECONDS s\"')
            print(f'SECONDS=0')
            print(f'now=$(date +"%T")')
            print(f'echo \"Started at : $now\"')

create_bash_script(location_list)

# run as
# python bash_script_creator.py > bash.script
