SECONDS=0
echo "started!"
python GCN_MILP_14_features.py --location 2022-08-31_12:04:42-GCN_32_2 --time_lim 10800 --mol_len 4 > GCN_outputs/2022-08-31_12:04:42-GCN_32_2_mol_len_4.txt
echo "progress: 5/9 after: $SECONDS s"
SECONDS=0
python GCN_MILP_14_features.py --location 2022-08-31_12:45:57-GCN_64_2 --time_lim 10800 --mol_len 4 > GCN_outputs/2022-08-31_12:45:57-GCN_64_2_mol_len_4.txt
echo "progress: 6/9 after: $SECONDS s"
SECONDS=0
python GCN_MILP_14_features.py --location 2022-08-31_12:06:21-GCN_16_3 --time_lim 10800 --mol_len 4 > GCN_outputs/2022-08-31_12:06:21-GCN_16_3_mol_len_4.txt
echo "progress: 7/9 after: $SECONDS s"
SECONDS=0
python GCN_MILP_14_features.py --location 2022-08-31_12:57:18-GCN_32_3 --time_lim 10800 --mol_len 4 > GCN_outputs/2022-08-31_12:57:18-GCN_32_3_mol_len_4.txt
echo "progress: 8/9 after: $SECONDS s"
SECONDS=0
python GCN_MILP_14_features.py --location 2022-08-31_13:02:07-GCN_64_3 --time_lim 10800 --mol_len 4 > GCN_outputs/2022-08-31_13:02:07-GCN_64_3_mol_len_4.txt
echo "progress: 9/9 after: $SECONDS s"
