import gzip
import json
from Bio import SeqIO


def build_seq_description_lst(descriptions_df, seq_dict):
    i = 0
    new_seq_id = []
    for key in seq_dict:
        name = descriptions_df['Target Name'].loc[i]
        start = descriptions_df['Target Env. Start'].loc[i]
        end = descriptions_df['Target Env. End'].loc[i]
        description = descriptions_df['Description'].loc[i]
        full_id = name + '/' + str(start) + '-' + str(end)
        if key == full_id:
            if "calcium" in description:
                new_seq_id.append((description, seq_dict[key]))
        i += 1

    return new_seq_id


def read_protein_seq(zip_file_path):
    """reading protein sequences from fa.gz file and build a dictionary {key= rec name, val = seq}"""
    zip_file_path = gzip.open(zip_file_path, "rt")
    fasta_file = SeqIO.parse(zip_file_path, 'fasta')
    seq_dict = {}
    i = 0
    for record in fasta_file:
        seq_dict[record.name] = record.seq

    print("num of sequances =", len(seq_dict))
    return seq_dict


def read_protein_seq_flybase(fasta_path):
    fasta_file = SeqIO.parse(fasta_path, 'fasta')
    seq_dict = {}
    for record in fasta_file:
        start = record.description.find("name")
        end = record.description.find("parent")
        new_label = record.description[start+5:end]
        if new_label.startswith("Ca-alpha1D"):
            new_label = 'Cav1_' + new_label
        elif new_label.startswith("cac"):
            new_label = 'Cav2_' + new_label
        elif new_label.startswith("Ca-alpha1T"):
            new_label = 'Cav3_' + new_label

        seq_dict[new_label] = record.seq

    print("num of sequances =", len(seq_dict))
    return seq_dict


def write_dict_to_fasta(seqs_dicts, output_path, seq_format = "seq_record"):
    output_file = open(output_path, 'w')
    for seq_query in seqs_dicts:
        for seq_id, seq in seq_query.items():
            identifier_line = ">" + seq_id + "\n"
            output_file.write(identifier_line)
            sequence_line = str(seq + "\n")
            output_file.write(sequence_line)

    output_file.close()


def write_list_to_fasta(seqs_lst, output_path):
    output_file = open(output_path, 'w')
    for seq_query in seqs_lst:
        identifier_line = ">" + seq_query[0] + "\n"
        output_file.write(identifier_line)
        sequence_line = str(seq_query[1] + "\n")
        output_file.write(sequence_line)

    output_file.close()
    print("Done writing fasta file!")


def discarded_short_seq(seq_dict):
    """create new dict contains only seq with more than 1500 amino acids/bases"""
    new_seq_dict = {}
    for seq in seq_dict:
        if len(seq_dict[seq]) > 1500:
            new_seq_dict[seq] = seq_dict[seq]
    print(len(new_seq_dict), "seq longer than 1500 bases")

    return new_seq_dict

def filter_CAC_seq(seq_dict):
    """create new dict contains only seq of α1 subunit proteins, encoded by the CACNA1x genes and edit the seq label"""
    new_seq_dict = {}
    by_type_dict = {}
    for seq in seq_dict:
        if seq.startswith("CAC"):
            new_seq_dict[seq[3:]] = seq_dict[seq]
            if seq.startswith("CAC1C") or seq.startswith("CAC1S") or seq.startswith("CAC1D") or seq.startswith("CAC1F") or seq.startswith("CAC1M"):
                new_label = 'Cav1_L-type_'
            elif seq.startswith("CAC1A"):
                new_label = 'Cav2_P/Q-type_'
            elif seq.startswith("CAC1B"):
                new_label = 'Cav2_N-type_'
            elif seq.startswith("CAC1E"):
                new_label = 'Cav2_R-type_'
            elif seq.startswith("CAC1H") or seq.startswith("CAC1G") or seq.startswith("CAC1I"):
                new_label = 'Cav3_T-type_'
            else:               # what to do with this kind of sequences not matching any known CAC gen??
                new_label = seq
            by_type_dict[new_label + seq] = seq_dict[seq]

    print(len(by_type_dict), "seq encoded by the CACNA1x genes")

    return by_type_dict


def read_json_file(json_path):
    """Open and read the JSON file"""
    with open(json_path, 'r') as file:
        data = json.load(file)
    print(data)


def read_msa_build_seq_dict(file_path):
    fasta_file = SeqIO.parse(file_path, 'fasta')
    seq_dict = {}
    for record in fasta_file:
        seq_dict[record.description] = record.seq

    return seq_dict

if __name__ == '__main__':
    #path for full length sequences fasta file from Hmmer, results for the search of Joel's sequences in SwissProt
    fullseq_path_cav1 = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\full_length_fasta\Q13936_CaV1.2 -fullseq.fa.gz"
    fullseq_path_cav2 = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\full_length_fasta\Q00975_CaV2.2 -fullseq.fa.gz"
    fullseq_path_cav3 = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\full_length_fasta\O43497_CaV3.1 -fullseq.fa.gz"
    fullseq_path_cav1_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\full_length_fasta\FlyBase_Ca-alpha1D-PJ -fullseq.fa.gz"
    fullseq_path_cav2_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\full_length_fasta\FlyBase_cac-PJ_green-exon -fullseq.fa.gz"
    fullseq_path_cav3_fly = r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\full_length_fasta\FlyBase_Ca-alpha1T -fullseq.fa.gz"

    print("Cav 1 L-type - Human")
    seq_dict_cav1 = read_protein_seq(fullseq_path_cav1)
    cac_seq_dict_cav1 = filter_CAC_seq(seq_dict_cav1)
    cac_seq_dict_cav1 = discarded_short_seq(cac_seq_dict_cav1)

    print("\nCav 2 N-type - Human")
    seq_dict_cav2 = read_protein_seq(fullseq_path_cav2)
    cac_seq_dict_cav2 = filter_CAC_seq(seq_dict_cav2)
    cac_seq_dict_cav2 = discarded_short_seq(cac_seq_dict_cav2)

    print("\nCav 3 T-type - Human")
    seq_dict_cav3 = read_protein_seq(fullseq_path_cav3)
    cac_seq_dict_cav3 = filter_CAC_seq(seq_dict_cav3)
    cac_seq_dict_cav3 = discarded_short_seq(cac_seq_dict_cav3)

    print("\nCav 1 L-type - Fly")
    seq_dict_cav1_fly = read_protein_seq(fullseq_path_cav1_fly)
    cac_seq_dict_cav1_fly = filter_CAC_seq(seq_dict_cav1_fly)
    cac_seq_dict_cav1_fly = discarded_short_seq(cac_seq_dict_cav1_fly)

    print("\nCav 2 P/Q-type - Fly")
    seq_dict_cav2_fly = read_protein_seq(fullseq_path_cav2_fly)
    cac_seq_dict_cav2_fly = filter_CAC_seq(seq_dict_cav2_fly)
    cac_seq_dict_cav2_fly = discarded_short_seq(cac_seq_dict_cav2_fly)

    print("\nCav 3 T-type - Fly")
    seq_dict_cav3_fly = read_protein_seq(fullseq_path_cav3_fly)
    cac_seq_dict_cav3_fly = filter_CAC_seq(seq_dict_cav3_fly)
    cac_seq_dict_cav3_fly = discarded_short_seq(cac_seq_dict_cav3_fly)


    if cac_seq_dict_cav1 == cac_seq_dict_cav2 and cac_seq_dict_cav1 == cac_seq_dict_cav3 and cac_seq_dict_cav1_fly == cac_seq_dict_cav2_fly and cac_seq_dict_cav1_fly == cac_seq_dict_cav3_fly and cac_seq_dict_cav1_fly == cac_seq_dict_cav1:
        print('\nall Human and Fly cac_seq_dict are equal (contain the same sequences)')


    #write_dict_to_fasta([cac_seq_dict_cav1], r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\CAC_homologs_from_hammer_original_labels.fasta")
    #write_dict_to_fasta([cac_seq_dict_cav1],r"C:\Users\sheer\Desktop\projectBenTalLab\HMMER\CAC_homologs_from_hammer_labels_by_type.fasta")

    print("\nsequences from FlyBase")
    flybase_cav1 = r"C:\Users\sheer\Desktop\projectBenTalLab\FlyBase\FlyBase_cav1.fasta"
    flybase_cav2 = r"C:\Users\sheer\Desktop\projectBenTalLab\FlyBase\FlyBase_cav2.fasta"
    flybase_cav3 = r"C:\Users\sheer\Desktop\projectBenTalLab\FlyBase\FlyBase_cav3.fasta"
    flybase_cav1_dict = discarded_short_seq(read_protein_seq_flybase(flybase_cav1))
    flybase_cav2_dict = discarded_short_seq(read_protein_seq_flybase(flybase_cav2))
    flybase_cav3_dict = discarded_short_seq(read_protein_seq_flybase(flybase_cav3))
    fly_dicts = [flybase_cav1_dict, flybase_cav2_dict, flybase_cav3_dict]
    write_dict_to_fasta(fly_dicts, r"C:\Users\sheer\Desktop\projectBenTalLab\FlyBase\merge_flybase_seq.fasta",
                        seq_format="seq_record")

    merge_flybase_swissprot = [flybase_cav1_dict, flybase_cav2_dict, flybase_cav3_dict, cac_seq_dict_cav1]
    write_dict_to_fasta(merge_flybase_swissprot,
                        r"C:\Users\sheer\Desktop\projectBenTalLab\FlyBase\merge_flybase_and_swissprot_seq.fasta",
                        seq_format="seq_record")









