import pandas as pd
from find_homologs_swissprot import read_msa_build_seq_dict
from sequence_specificity_analysis import analyze_mutations


def create_dicts_by_type_from_msa(msa_path):
    """read mafft msa file and divide the sequences to dicts [record.name: record.seq] by CaV type,
     based on the seq annotation"""
    msa_dict = read_msa_build_seq_dict(msa_path)
    cav1_dict = {}
    cav2_dict = {}
    cav3_dict = {}
    for seq in msa_dict:
        if seq.find('L-type') != -1 or seq.find('Cav1') != -1:
            cav1_dict[seq] = msa_dict[seq]
        elif seq.find('type A') != -1 or seq.find('Cav2') != -1 or seq.find('Cacophony') != -1:
            cav2_dict[seq] = msa_dict[seq]
        elif seq.find('T-type') != -1 or seq.find('subunit T') != -1 or seq.find('Cav3') != -1:
            cav3_dict[seq] = msa_dict[seq]

    return cav1_dict, cav2_dict, cav3_dict


def find_identical_positions_in_specific_type(type_seq_dict, id_threshold=0.9):
    """
     given a list of sequences from a specific type, return a dict with positions where all the seq
     from that type "agree" on (i.e has the same amino acid in that position) excluded gaps
     """
    identical_positions = {}
    last_pos = len(sorted(type_seq_dict.values())[0])
    n = len(type_seq_dict.keys())
    for pos in range(1,last_pos+1):  # for every position in the alignment
        curr_position_dict = {}
        for seq in type_seq_dict.values():
            if seq[pos-1] not in curr_position_dict.keys():
                curr_position_dict[seq[pos-1]] = 0
            curr_position_dict[seq[pos-1]] += 1

        curr_position_dict = sorted(curr_position_dict.items(), key=lambda x: x[1], reverse=True)

        most_frequent_aa = curr_position_dict[0][0]
        most_frequent_aa_cnt = curr_position_dict[0][1]

        if most_frequent_aa_cnt > id_threshold*n and most_frequent_aa != "-":
            identical_positions[pos] = most_frequent_aa

    return identical_positions



def full_compare_identical_positions_between_types(all_types_positions):
    """
     find positions that are identical between all (or most) the sequences from a one specific type,
     even if they are not identical between seq from another type
     """
    compare_positions_dict = {}
    for type_positions in all_types_positions:
        for position in type_positions:
            compare_positions_dict[position] = (all_types_positions[0].get(position), all_types_positions[1].get(position), all_types_positions[2].get(position))

    myKeys = list(compare_positions_dict.keys())
    myKeys.sort()
    sorted_dict = {i: compare_positions_dict[i] for i in myKeys}

    return sorted_dict


def partial_compare_identical_positions_between_types(merge_positions_lst):
    """find only positions that are identical between all the sequences from a specific
     type and different from another type"""
    partial_compare_positions_dict = {}
    for position in merge_positions_lst[0]:
        if position in merge_positions_lst[1] and position in merge_positions_lst[2]:
            if not (merge_positions_lst[0].get(position) == merge_positions_lst[1].get(position) == merge_positions_lst[2].get(position)):
                partial_compare_positions_dict[position] = (merge_positions_lst[0].get(position), merge_positions_lst[1].get(position), merge_positions_lst[2].get(position))

    return partial_compare_positions_dict


def compare_identical_positions_cav1_cav2(merge_positions_lst):
    """find only positions that are identical between all sequences of Cav1 and different from all sequences of Cav2"""
    compare_positions_cav1_cav2 = {}
    for position in merge_positions_lst[0]:
        if position in merge_positions_lst[1]:
            if not (merge_positions_lst[0].get(position) == merge_positions_lst[1].get(position)):
                compare_positions_cav1_cav2[position] = (merge_positions_lst[0].get(position), merge_positions_lst[1].get(position))

    return compare_positions_cav1_cav2


def run_comparison(msa_path, comparison_type="full", taxonomy="", analyze = True):
    cav_index = ['cav1', 'cav2', 'cav3']
    insecta_cav_index = ['cav1 - L type', 'cav2 - A type', 'cav3 - T type']
    index = insecta_cav_index if taxonomy == "insecta" else cav_index

    cav1_dict, cav2_dict, cav3_dict = create_dicts_by_type_from_msa(msa_path)
    num_seq = len(cav1_dict) + len(cav2_dict) + len(cav3_dict)
    print("num of Cav1 seq =", len(cav1_dict), "\nnum of Cav2 seq =", len(cav2_dict),"\nnum of Cav3 seq =", len(cav3_dict))
    cav1_identical_positions = find_identical_positions_in_specific_type(cav1_dict)
    cav2_identical_positions = find_identical_positions_in_specific_type(cav2_dict)
    cav3_identical_positions = find_identical_positions_in_specific_type(cav3_dict)
    all_types_positions_lst = [cav1_identical_positions, cav2_identical_positions, cav3_identical_positions]

    if comparison_type == "partial":
        compare_positions = partial_compare_identical_positions_between_types(all_types_positions_lst)
    elif comparison_type == "full":
        compare_positions = full_compare_identical_positions_between_types(all_types_positions_lst)
    elif comparison_type == "cav1_2":
        compare_positions = compare_identical_positions_cav1_cav2(all_types_positions_lst)
        index = ['cav1', 'cav2']
    else:
        print("need to specify comparison_type (full / partial / cav1_2)")
        return

    compare_positions_df = pd.DataFrame(compare_positions, index=index)
    if analyze:
        analyze_mutations(compare_positions_df, comparison_type, taxonomy, num_seq)

    else:
        compare_positions_df.to_excel(f"{num_seq}seq_{taxonomy}_{comparison_type}_compare_positions.xlsx")
        print(f'{comparison_type}_compare_positions:\n', compare_positions_df)




if __name__ == '__main__':
    # compare 43 cav sequences from SwissProt
    mafft_msa_path = r"C:\Users\sheer\Desktop\projectBenTalLab\mafft\mafft_msa_43CACseq_labels_by_type.fasta"

    # compare 367 insecta cav sequences from UniProt (after apply cd-hit)
    mafft_insecta_msa_path_cd_hit = r"C:\Users\sheer\Desktop\projectBenTalLab\mafft\msa on 367 seq.fasta"

    # compare 1000 insecta cav sequences from UniProt (without cd-hit)
    mafft_insecta_msa_path = r"C:\Users\sheer\Desktop\projectBenTalLab\mafft\msa only insecta sequences from uniprot.fasta"

    #run_comparison(mafft_msa_path, comparison_type="partial", analyze=True)
    #run_comparison(mafft_insecta_msa_path_cd_hit, comparison_type="partial", taxonomy="insecta", analyze=True)
    #run_comparison(mafft_insecta_msa_path_cd_hit, comparison_type="cav1_2", taxonomy="insecta", analyze=True)


