from Bio import SeqIO
import pandas as pd
from find_homologs_swissprot import read_protein_seq, discarded_short_seq, write_dict_to_fasta


"""read data from TSV file"""
def read_tsv_file(tsv_path):
    cav_df = pd.read_csv(tsv_path, sep='\t')
    df = pd.DataFrame(data=cav_df, columns=['Target Name', 'Target Env. Start', 'Target Env. End', 'Description'])
    return df


"""create a mapping dictionary from DataFrame (tsv_df) and update the original dictionary with descriptions"""
def filter_calcium_channels_seq(seq_dict, tsv_df):
    description_mapping = dict(zip(tsv_df['Target Name'], tsv_df['Description']))
    updated_dict = {}
    for key in seq_dict:
        new_key = key + "-" + description_mapping.get(key, 'Description not found')
        updated_dict[new_key] = seq_dict.get(key)

    """after create an updated dictionary with the sequences descriptions, filter seq described as 'calcium' or 'Ca'"""
    calcium_dict = {}
    not_calcium_dict = {}
    cav1_Ltype = {}
    cav2_Atype = {}
    cav3_Ttype = {}
    other_types = {}
    for description in updated_dict:
        if description.find('calcium') != -1 or description.find('Ca') != -1:
            calcium_dict[description] = updated_dict.get(description)
            if description.find('L-type') != -1:
                cav1_Ltype[description] = updated_dict.get(description)
            elif description.find('type A') != -1:
                cav2_Atype[description] = updated_dict.get(description)
            elif description.find('T-type') != -1:
                cav3_Ttype[description] = updated_dict.get(description)
            else:
                other_types[description] = updated_dict.get(description)
        else:
            not_calcium_dict[description] = updated_dict.get(description)

    return calcium_dict, not_calcium_dict, cav1_Ltype, cav2_Atype, cav3_Ttype, other_types


def merge_calcium_dictionaries(dict1, dict2):
    merge_dictionaries = dict2
    for seq in dict1:
        if seq not in dict2:
            merge_dictionaries[seq] = dict1.get(seq)

    return merge_dictionaries


def read_msa_build_dicts_by_type(curr_type_dict, curr_type_msa_path):
    fasta_file = SeqIO.parse(curr_type_msa_path, 'fasta')
    curr_type_seq_dict = {}
    for record in fasta_file:
        curr_type_seq_dict[record.description] = curr_type_dict.get(record.description)

    return curr_type_seq_dict


if __name__ == '__main__':
    """
    path for full length sequences fasta file from Hmmer,
     results for the search of Joel's sequences in UniProt, restricted by taxonomy - Insecta
     """

    fullseq_path_cav1_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\uniprot_insecta\FlyBase_Ca-alpha1T -fullseq_insecta.gz"
    tsv_path_cav1_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\uniprot_insecta\FlyBase_Ca-alpha1T_insecta.tsv"
    print("Cav 1 L-type - Fly - UniProt")
    seq_dict_cav1_fly = read_protein_seq(fullseq_path_cav1_fly)
    cac_seq_dict_cav1_fly = discarded_short_seq(seq_dict_cav1_fly)
    tsv_cav1_df = read_tsv_file(tsv_path_cav1_fly)
    calcium_dict, not_calcium_dict, cav1_Ltype, cav2_Atype, cav3_Ttype, other_types = \
        filter_calcium_channels_seq(cac_seq_dict_cav1_fly, tsv_cav1_df)
    print('calcium_dict', len(calcium_dict))
    print('not_calcium_dict', len(not_calcium_dict.keys()))
    print('cav1_Ltype', len(cav1_Ltype))
    print('cav2_Atype', len(cav2_Atype))
    print('cav3_Ttype', len(cav3_Ttype))
    print('other_types', len(other_types))

    fullseq_path_cav2_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\uniprot_insecta\FlyBase_cac-PJ_green-exon -fullseq_insecta.gz"
    tsv_path_cav2_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\uniprot_insecta\FlyBase_cac-PJ_green-exon -fullseq_insecta.tsv"
    print("\nCav 2 P/Q-type - Fly - UniProt")
    seq_dict_cav2_fly = read_protein_seq(fullseq_path_cav2_fly)
    cac_seq_dict_cav2_fly = discarded_short_seq(seq_dict_cav2_fly)
    tsv_cav2_df = read_tsv_file(tsv_path_cav2_fly)
    filtered_calcium_dict2 = filter_calcium_channels_seq(cac_seq_dict_cav2_fly, tsv_cav2_df)[0]
    print(len(filtered_calcium_dict2))

    fullseq_path_cav3_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\uniprot_insecta\FlyBase_Ca-alpha1T -fullseq_insecta.gz"
    tsv_path_cav3_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\uniprot_insecta\FlyBase_Ca-alpha1T_insecta.tsv"
    print("\nCav 3 T-type - Fly - UniProt")
    seq_dict_cav3_fly = read_protein_seq(fullseq_path_cav3_fly)
    cac_seq_dict_cav3_fly = discarded_short_seq(seq_dict_cav3_fly)
    tsv_cav3_df = read_tsv_file(tsv_path_cav3_fly)
    filtered_calcium_dict3 = filter_calcium_channels_seq(cac_seq_dict_cav3_fly, tsv_cav3_df)[0]
    print(len(filtered_calcium_dict3))

    if calcium_dict == filtered_calcium_dict2 and calcium_dict == filtered_calcium_dict3:
        print('\ncav1_dict, cav2_dict and cav3_dict are all equal (contain the same sequences)')


    #write_dict_to_fasta([cav1_Ltype], r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\insecta_cav1_Ltype_homologs_from_uniprot_hammer.fasta")
    #write_dict_to_fasta([cav2_Atype], r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\insecta_cav2_Atype_homologs_from_uniprot_hammer.fasta")
    #write_dict_to_fasta([cav3_Ttype], r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\insecta_cav3_Ttype_homologs_from_uniprot_hammer.fasta")
    #write_dict_to_fasta([other_types], r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\insecta_other_types_homologs_from_uniprot_hammer.fasta")

    cav1_Ltype_msa_path = r"C:\Users\sheer\Desktop\projectBenTalLab\mafft\mafft job msa 465 seq L type  after cd hit.fasta"
    cav1_seq_dict_after_cdHit = read_msa_build_dicts_by_type(cav1_Ltype, cav1_Ltype_msa_path)
    print('cav1 after cd-hit ', len(cav1_seq_dict_after_cdHit))

    cav2_Atype_msa_path = r"C:\Users\sheer\Desktop\projectBenTalLab\mafft\mafft job mafft job msa 394 seq A type after cd hit.fasta"
    cav2_seq_dict_after_cdHit = read_msa_build_dicts_by_type(cav2_Atype, cav2_Atype_msa_path)
    print('cav2 after cd-hit ', len(cav2_seq_dict_after_cdHit))

    cav3_Ttype_msa_path = r"C:\Users\sheer\Desktop\projectBenTalLab\mafft\msa 101 seq cav3 T type after cd hit.fasta"
    cav3_seq_dict_after_cdHit = read_msa_build_dicts_by_type(cav3_Ttype, cav3_Ttype_msa_path)
    print('cav3 after cd-hit ', len(cav3_seq_dict_after_cdHit))

    merged_dict0 = merge_calcium_dictionaries(cav1_seq_dict_after_cdHit, cav2_seq_dict_after_cdHit)
    merged_dict1 = merge_calcium_dictionaries(cav3_seq_dict_after_cdHit, other_types)
    merged_dict2 = merge_calcium_dictionaries(merged_dict0, merged_dict1)
    print('merged_dict contains', len(merged_dict2), 'sequences')

    write_dict_to_fasta([merged_dict2], r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\insecta_all_types_homologs_after_cdHit.fasta")





