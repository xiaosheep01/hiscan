# -*- coding: utf-8 -*-
# Development Time: 2023-09-01 10:30:22
# Developer: XiaoYang
import itertools
import os
import pathlib
import platform
import linecache
import re
import numpy
import pandas as pd
from colorama import Fore


def file_type_judge(file_path):
    """
    Determine the type of the input file, apply different reading methods according to different types,
    and return a pandas type data frame.
    :param file_path: input file path
    :return: a pandas dataframe
    """
    if file_path.split(".")[1] == "xlsx":
        content = pd.read_excel(file_path, header=0)
    elif file_path.split(".")[1] == "csv":
        content = pd.read_csv(file_path, sep=",", header=0)
    elif file_path.split(".")[1] == "txt":
        content = pd.read_csv(file_path, sep="\t", header=0)
    else:
        print("%s Warning: Try to read the file in txt format: %s %s" %
              (Fore.YELLOW, Fore.BLUE, file_path))
        content = pd.read_table(file_path, sep="\t", header=0)
        print("%s Note: The file was read successfully: %s %s" %
              (Fore.CYAN, Fore.BLUE, file_path))
    return content


def files_path(dirname, sys_os):
    """
    Gets all file paths of the folder
    :param dirname: path
    :return: files path list
    """
    file_path_list = []
    dir_path_list = []
    for root, dirs, files in os.walk(dirname):
        # root represents the path to the folder being accessed
        # dirs represents the names of all subfolders in this folder, and the list is stored
        # files said evil all child filename, under the folder list is stored
        for item in files:
            file_path_list.append(os.path.join(root, item))

        for item in dirs:
            dir_path_list.append(os.path.join(root, item))

    # Judgment system type
    # MacOS
    if sys_os == "Darwin":
        print("{}Note:Current operating system: {}, remove hidden file path by "
              "default.{}".format(Fore.YELLOW, platform.uname()[0], Fore.RESET))
        for item in file_path_list:
            item_len = len(item.split(os.sep))
            if item.split(os.sep)[item_len-1].startswith("."):
                file_path_list.remove(item)

        # Force out.ds_store files
        file_path_list.remove(os.path.join(dirname, ".DS_Store"))

    # Windows
    elif sys_os == "Windows":
        pass

    # Linux
    elif sys_os == "Linux":
        pass

    return file_path_list


def files_combine(file_path_list):
    content = []

    for item in file_path_list:
        for line in linecache.getlines(item):
            content.append(line)

    return content


def fasta_doc_std(file_list):
    """
    Standardization of fasta documents
    :param file_list:file list by function:files_combine
    :return:Standard fasta content
    """

    seq_name_list = []
    seq_dic = {}
    seq_list = []
    seq_name = ""

    for line in file_list:
        line = line.strip()
        if line.startswith(">"):
            seq_name = line
            seq_name_list.append(seq_name)
            seq_dic[seq_name] = []
        else:
            seq_dic[seq_name].append(line)

    for seq_name in seq_name_list:
        seq_list.append(seq_name + "\n" + "".join(seq_dic[seq_name]))

    return seq_list


def mimic_find(mimic_file_path, seq_list):
    """
    Query whether there is a specified string in the target text, and will result in the form of a list of output
    :param mimic_file_path:mimic file path
    :param seq_list:The path to the file that needs to be queried
    :return:Return match results
    """
    file_suffix = pathlib.Path(mimic_file_path).suffix

    if file_suffix == ".xlsx":
        mimic_list = pd.read_excel(mimic_file_path, header=None)
    elif file_suffix == ".csv":
        mimic_list = pd.read_csv(mimic_file_path, header=None, sep=",")
    elif file_suffix == ".txt":
        mimic_list = pd.read_table(mimic_file_path, header=None, sep="\t")
    else:
        print("{}---Error!Please check that the mimic file suffix is correct!---{}".format(
            Fore.RED, Fore.RESET
        ))

    # Prepare the lists
    NCBI_list = []
    seq_name_list = []
    mimic_name_list = []
    location_list = []
    histone_subunit_list = []
    modification_list = []

    for line in mimic_list.itertuples():

        for item in seq_list:
            annotation_info, seq = item.split("\n")
            NCBI_ID, seq_name = annotation_info.split(" ", 1)
            NCBI_ID = NCBI_ID.replace(">", "")

            # query mimic in sequence
            match_obj = re.finditer(line[1], seq)

            for item in match_obj:
                location = "(" + str(item.span()[0] + 1) + ":" + str(item.end()) + ")"
                # line.1 for mimic; line._2 for histone subunit; line._3 for modification
                # if true then write in lists
                NCBI_list.append(NCBI_ID)                   # NCBI ID of sequence
                seq_name_list.append(seq_name)              # sequence name
                mimic_name_list.append(line[1])             # mimic
                """
                line[1] for mimic;
                line[2] for histone subunits;
                line[3] for modification
                line._fields for the amount of columns, containing index column. Usually [0] for index col.
                """
                # print(len(line._fields))
                if 2 < len(line._fields) < 5:
                    histone_subunit_list.append(line[2])        # histone subunit
                    modification_list.append(line[3])           # modification
                elif len(line._fields) >= 5:
                    print("\033[34m---Note:Excess annotation information (after the third column) is not displayed "
                          "in the result file.---\033[0m")
                    histone_subunit_list.append(line[2])  # histone subunit
                    modification_list.append(line[3])  # modification

                location_list.append(location)              # location in sequence

    if 2 < len(line._fields) < 5:
        # write in xlsx
        df = pd.DataFrame({"NCBI_ID": NCBI_list, "Seq_Name": seq_name_list, "Mimic": mimic_name_list,
                           "Location": location_list, "Histone_Subunit": histone_subunit_list,
                           "Modification": modification_list})
    elif len(line._fields) == 2:
        df = pd.DataFrame({"NCBI_ID": NCBI_list, "Seq_Name": seq_name_list, "Mimic": mimic_name_list,
                           "Location": location_list})

    return df


def mimic_skewness(mimic_df):
    """
    The mimic motif amino acid skewness is only applicable to the case of 5 amino acids!
    :param mimic_df: mimic dataframe
    :return: a dataframe(matrix format) containing the probability of each amino acid
    """
    # mimic should be in col 0
    my_mimic_df = mimic_df.iloc[:, 0]
    my_mimic_df = my_mimic_df.str.split("")

    def remove_blank(x):            # the first to fifth amino acid
        x = x[1:6]
        return x

    my_mimic_df = my_mimic_df.apply(remove_blank)

    def combine(x):
        x = ",".join(x)
        return x

    my_mimic_df = my_mimic_df.apply(combine)

    new_mimic_df = my_mimic_df.str.split(",", expand=True)
    pro_mimic_df = pd.DataFrame()           # mimic probability dataframe
    old_col_name_list = list(new_mimic_df.columns)
    new_col_name_list = ["col_" + str(x) for x in old_col_name_list]
    new_mimic_df.columns = new_col_name_list

    best_mimic_list = []
    best_pro_list = []
    skewness_df = pd.DataFrame()

    for col in list(new_mimic_df.columns):
        pro_dict = {}           # mimic probability dictionary of pei column
        pro_res = new_mimic_df[col].value_counts(normalize=True)

        pro_res_test = pro_res.sort_values(ascending=False)
        best_mimic_list.append(pro_res_test.index.tolist()[0])
        best_pro_list.append(pro_res_test[0])

        index_list = pro_res.index.tolist()

        for index in index_list:
            pro_dict[index] = pro_res.loc[index]

        def pro_match(x):
            x = pro_dict[x]
            return x

        pro_mimic_df[col] = new_mimic_df[col].apply(pro_match)

        def skew_value(x):
            value_tuple = (x, round(pro_dict[x], 5))
            return value_tuple

        skewness_df[col] = new_mimic_df[col].apply(skew_value)

    # calculate best mimic and probability
    best_mimic = "".join(best_mimic_list)
    best_pro = numpy.prod(best_pro_list)
    print()
    print(">>> Highest possible mimic: %s%s" % (Fore.LIGHTYELLOW_EX, best_mimic))
    print(">>> Probability: %s%s" % (Fore.LIGHTYELLOW_EX, best_pro))
    print()

    return skewness_df, new_mimic_df, pro_mimic_df          # return the aa skewness of mimics


def mimic_predict(new_mimic_df, pro_mimic_df):
    """
    To predict possible motifs based on the input motifs, all amino acid permutations and combinations are
    calculated and their probabilities are calculated. It is recommended that no more than 15 motifs be
    entered due to excessive memory consumption.
    :param new_mimic_df: a dataframe of mimic amino acid
    :param pro_mimic_df: a dataframe of mimic amino acid probability
    :return: a dataframe containing all potential motifs combinations with probability
    """
    mimic_num = len(new_mimic_df)
    if 15 <= mimic_num < 20:
        print(Fore.CYAN + "Note: The number of motifs is more than 15, which will take some time to process!")
    elif 20 <= mimic_num < 30:
        print(Fore.YELLOW + "Warning: The number of motifs is more than 20, which will take some time to process!")
    elif mimic_num >= 30:
        print(Fore.YELLOW + "Warning: The number of motifs is more than 30, which will take lots of time to process!")

    # method 1: use itertools package
    new_mimic_res = list(itertools.product(list(new_mimic_df["col_0"]),
                                           list(new_mimic_df["col_1"]),
                                           list(new_mimic_df["col_2"]),
                                           list(new_mimic_df["col_3"]),
                                           list(new_mimic_df["col_4"])))

    new_pro_res = list(itertools.product(
        list(pro_mimic_df["col_0"]),
        list(pro_mimic_df["col_1"]),
        list(pro_mimic_df["col_2"]),
        list(pro_mimic_df["col_3"]),
        list(pro_mimic_df["col_4"])
    ))

    aa_string_list = []
    for i in new_mimic_res:
        aa_string_list.append("".join(i))

    pro_list = []
    for i in new_pro_res:
        pro_list.append(numpy.prod(i))

    result_df = pd.DataFrame()
    result_df["Mimic"] = aa_string_list
    result_df["Probability"] = pro_list
    result_df.sort_values(by="Probability", ascending=False, inplace=True)
    result_df = result_df.drop_duplicates()
    # print(result_df["Probability"].sum())
    return result_df


def obtain_annotation(file_list):
    annotation_list = []
    file_len = len(file_list)

    for num in range(file_len):
        annotation_string = ""

        if file_list[num].startswith("VERSION"):
            NCBI_ID = file_list[num].split(" ", 1)[1]
            annotation_list.append(NCBI_ID.strip())

        elif "ORGANISM" in file_list[num]:
            while not file_list[num + 1].startswith("REFERENCE"):
                annotation_string = annotation_string + file_list[num + 1].strip()
                num = num + 1
            annotation_string = annotation_string.replace(" ", "")
            annotation_string = annotation_string.replace(".", "")
            annotation_list.append(annotation_string)

    """
    "annotation_list" contains the annotation of sequence
    the format is a line for the version name and 
    another line for virus classification.
    """

    # creat classification dict
    realm_dict = {}
    king_dict = {}
    phylum_dict = {}
    class_dict = {}
    order_dict = {}
    family_dict = {}

    for num in range(0, len(annotation_list), 2):

        class_list = annotation_list[num + 1].split(";")
        # Realm info
        if "viria" in annotation_list[num + 1]:
            for item in class_list:
                if item.endswith("viria"):
                    realm_dict[annotation_list[num]] = item
        else:
            realm_dict[annotation_list[num]] = "None"
        # Kingdom info
        if "virae" in annotation_list[num + 1]:
            for item in class_list:
                if item.endswith("virae"):
                    king_dict[annotation_list[num]] = item
        else:
            king_dict[annotation_list[num]] = "None"
        # Phylum info
        if "viricota" in annotation_list[num + 1]:
            for item in class_list:
                if item.endswith("viricota"):
                    phylum_dict[annotation_list[num]] = item
        else:
            phylum_dict[annotation_list[num]] = "None"
        # Class info
        if "viricetes" in annotation_list[num + 1]:
            for item in class_list:
                if item.endswith("viricetes"):
                    class_dict[annotation_list[num]] = item
        else:
            class_dict[annotation_list[num]] = "None"
        # Order info
        if "virales" in annotation_list[num + 1]:
            for item in class_list:
                if item.endswith("virales"):
                    order_dict[annotation_list[num]] = item
        else:
            order_dict[annotation_list[num]] = "None"
        # Family info
        if "viridae" in annotation_list[num + 1]:
            for item in class_list:
                if item.endswith("viridae"):
                    family_dict[annotation_list[num]] = item
        else:
            family_dict[annotation_list[num]] = "None"

    return realm_dict, king_dict, phylum_dict, class_dict, order_dict, family_dict


def add_annotation_to_xlsx(mimic_df, dict_tuples):

    realm_info = dict_tuples[0]
    kingdom_info = dict_tuples[1]
    phylum_info = dict_tuples[2]
    class_info = dict_tuples[3]
    order_info = dict_tuples[4]
    family_info = dict_tuples[5]

    # add annotation
    mimic_df["Realm"] = mimic_df["NCBI_ID"].map(realm_info)
    mimic_df["Kingdom"] = mimic_df["NCBI_ID"].map(kingdom_info)
    mimic_df["Phylum"] = mimic_df["NCBI_ID"].map(phylum_info)
    mimic_df["Class"] = mimic_df["NCBI_ID"].map(class_info)
    mimic_df["Order"] = mimic_df["NCBI_ID"].map(order_info)
    mimic_df["Family"] = mimic_df["NCBI_ID"].map(family_info)

    return mimic_df


def add_host_source(host_ref_file_path, mimic_df):
    host_df = pd.read_excel(host_ref_file_path, sheet_name=0, header=0)
    species_host = host_df[["Species", "Host source"]]
    virus_host = host_df[["Virus name(s)", "Host source"]]
    species_dict = {}
    virus_name_dict = {}

    for index, row in species_host.iterrows():
        key = row.iloc[0]
        value = row.iloc[1]
        species_dict[key] = value

    for index, row in virus_host.iterrows():
        key = row.iloc[0].upper()                # transfer in upper format
        value = row.iloc[1]
        virus_name_dict[key] = value

    # input mimic annotation file
    mimic_file = mimic_df

    # obtain sequence name
    def names(x):
        x = re.findall("\[(.*?)\]", x)[0]
        return x

    mimic_file["Species"] = mimic_file["Seq_Name"].apply(names)

    # match host source
    def host(x):
        if x in species_dict:
            return species_dict[x]
        elif x.upper() in virus_name_dict:          # transfer in upper format
            return virus_name_dict[x.upper()]
        else:
            return "None"

    mimic_file["Host_Source"] = mimic_file["Species"].apply(host)

    return mimic_file


def mimic_stat(result_df):
    mimic_count_res = pd.DataFrame(result_df["Mimic"].value_counts())
    mimic_count_res = mimic_count_res.reset_index()
    mimic_count_res_freq = pd.DataFrame(result_df["Mimic"].value_counts(normalize=True))
    mimic_count_res_freq = mimic_count_res_freq.reset_index()
    result = pd.merge(mimic_count_res, mimic_count_res_freq, on="index")
    result.columns = ["Mimic", "Count", "Frequency"]
    result = result.round({"Frequency": 5})
    print(result.head(5))
    return result


def class_count(result_df, class_type):
    col_list = list(result_df.columns)
    col_list = [item.upper() for item in col_list]
    result_df.columns = col_list
    class_count_res = pd.DataFrame(result_df[class_type].value_counts())
    class_count_res = class_count_res.reset_index()
    class_count_res_freq = pd.DataFrame(result_df[class_type].value_counts(normalize=True))
    class_count_res_freq = class_count_res_freq.reset_index()
    result = pd.merge(class_count_res, class_count_res_freq, on="index")
    result.columns = [class_type, "Count", "Frequency"]
    print(result.head())
    return result


def host_count(result_df):
    host_count_res = pd.DataFrame(result_df["Host_Source"].value_counts())
    host_count_res = host_count_res.reset_index()
    host_count_res_freq = pd.DataFrame(result_df["Host_Source"].value_counts(normalize=True))
    host_count_res_freq = host_count_res_freq.reset_index()
    result = pd.merge(host_count_res, host_count_res_freq, on="index")
    result.columns = ["Host_Source", "Count", "Frequency"]
    print(result.head())
    return result
