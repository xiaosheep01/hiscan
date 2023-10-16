# -*- coding: utf-8 -*-
# Development Time: 2023-09-01 10:30:30
# Developer: XiaoYang

import argparse
import os
import platform
import sys
import time
from colorama import init, Fore
import pandas as pd
import func

example_use = r'''
{}********************** Example of use ***********************{}

  (1) If you have a sequence file(or a folder) and a motif file (without any annotations file):
      hiscan -i your/sequence/path -m your/motif/path -o your/output/path

  (2) If you have a sequence file(or a folder) and a motif file (with NCBI and ICTV annotation):
      hiscan -i your/sequence/path -m your/motif/path -na your/NCBI/file/path -ia your/ICTV/file/path -o your/output/path

  (3) If you already have a motif result to do motif skewness statistics:
      hiscan -i your/motif/result/path -ms -o your/output/path -on your_result_name

  Tip: the input(file/directory) and output(directory) is recommended absolute path.

  Above is just a conceptual example, detailed usage in website: {}https://github.com/xiaosheep01/hiscan

{}***************************  End  ***************************{}

'''.format(Fore.GREEN, Fore.RESET, Fore.BLUE, Fore.GREEN, Fore.RESET)


def starts():
    print(Fore.GREEN + "\n" + "===================================================================")

    print("{}>>> {}Name: Histone Motif Scan (HiScan)".format(Fore.GREEN, Fore.RESET))

    print(
        "{}>>> {}Description: Scanning histone motif in viral protein sequences.".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Version: 1.0 (2023-09-08)".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Author: Yang Xiao".format(Fore.GREEN, Fore.RESET))

    print("{}>>> {}Email: Fredrik1999@163.com".format(Fore.GREEN, Fore.RESET))

    # print("  Citation: XXXX")

    print(Fore.GREEN + "===================================================================" + "\n" + Fore.RESET)

    def parameters():

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            prog="hiscan",
            description="",
            epilog=example_use)

        result_group = parser.add_mutually_exclusive_group()

        result_group.add_argument("-i", type=str, metavar="",
                                  dest="input",
                                  help="the path of the input files")

        result_group.add_argument("-mr", type=str, metavar="",
                                  dest="mimic_result",
                                  help="please input the result file of '-m' parameter")

        result_group.add_argument("-nr", type=str, metavar="",
                                  dest="ncbi_annotation_result",
                                  help="please input the result file of '-na' parameter")

        parser.add_argument("-m", type=str, metavar="",
                            dest="mimic",
                            help="look for the mimic in the sequence file and you need input the mimic file")

        parser.add_argument("-na", type=str, metavar="",
                            dest="ncbi_annotation",
                            help="add NCBI annotation for the result and please input the NCBI annotation file "
                                 "for the query")

        parser.add_argument("-ia", type=str, metavar="",
                            dest="ictv_annotation",
                            help="add ICTV annotation for the result. "
                                 "please input the ICVT annotation file for the query and "
                                 "you could get it from this website: https://ictv.global/vmr")

        count_group = parser.add_mutually_exclusive_group()

        count_group.add_argument("-mc",
                                 dest="mimic_count",
                                 action="store_true",
                                 help="The type and number of mimic in the result file will be counted")

        count_group.add_argument("-ms",
                                 dest="mimic_skewness",
                                 action="store_true",
                                 help="The skewness of each amino acid of input mimic will be calculated")

        count_group.add_argument("-mp",
                                 action="store_true",
                                 dest="mimic_predict",
                                 help="possible motifs will be predicted based on the input motif files")

        count_group.add_argument("-cc",
                                 action="store_true",
                                 dest="classification_count",
                                 help="The results were counted according to the specified classification level, "
                                      "which by default is 'family'. It can be specified with the '-ct' parameter")

        count_group.add_argument("-hc",
                                 action="store_true",
                                 dest="host_count",
                                 help="The type and number of hosts in the result file will be counted")

        parser.add_argument("-ct", type=str, metavar="",
                            dest="class_type",
                            default="family",
                            help="Specify the classification criteria for categorical statistics. "
                                 "The default is' Family 'and you can choose from 'Realm', "
                                 "'Kingdom', 'Phylum', 'Class', 'Order', and 'Family'")

        parser.add_argument("-o", type=str, metavar="",
                            dest="output",
                            default="",
                            help="the path of output file")

        parser.add_argument("-on", type=str, metavar="",
                            dest="output_name",
                            default="Result",
                            help="The name of the output file, default is 'Result'")

        parser.add_argument("-ot", type=str, metavar="",
                            dest="output_type",
                            default="txt",
                            help="the format of output file. you can choose one of formats "
                                 "including 'txt', 'xlsx' or 'csv' as output file format "
                                 "('txt' in default, and 'xlsx' is recommended) ")

        my_args = parser.parse_args(sys.argv[1:])
        return my_args

    my_args = parameters()
    start_time = time.time()
    user_sys_os = platform.uname()[0]
    result = ""
    merged_file = ""
    # MacOS
    if user_sys_os == "Darwin":
        print("{}Note: The current operating system is: {}, remove hidden file path in "
              "default.{}".format(Fore.CYAN, user_sys_os, Fore.RESET))
    # Windows
    elif user_sys_os == "Windows":
        print("{}Note: The current operating system is: {}, if there are any hidden files please manually delete.{}"
              .format(Fore.CYAN, user_sys_os, Fore.RESET))
    # Linux
    elif user_sys_os == "Linux":
        print("{}Note: The current operating system is: {}, if there are any hidden files please manually delete.{}"
              .format(Fore.CYAN, user_sys_os, Fore.RESET))

    """
    Step 1: Document preparation stage.
    In this step, the virus sequence files are merged and standardized.
    This step is necessary!
    """

    if my_args.input:
        if any([my_args.mimic_predict, my_args.mimic_count,
                my_args.classification_count, my_args.host_count,
                my_args.mimic_skewness]):
            pass
        else:
            if os.path.isdir(my_args.input):
                folder_path = func.files_path(my_args.input, user_sys_os)  # The folder path
                merged_file = func.files_combine(folder_path)  # Combined file list
                print("---Files Merged Completed---")
            elif os.path.isfile(my_args.input):
                merged_file = func.files_combine([my_args.input])
            else:
                print(Fore.RED + "---Error! Please make sure that the input is a file or folder, "
                                 "instead of links and so on!---")
                sys.exit()

    """
    Step 2: Find the mimic in the sequence file and write the corresponding information.
    This step is optional.
    """

    if my_args.mimic:
        my_seq = func.fasta_doc_std(merged_file)  # Combined file std
        result = func.mimic_find(my_args.mimic, my_seq)  # return result dataframe
        print("---Mimic Query Completed---")

    # -ms: mimic skewness
    if my_args.mimic_skewness:
        print("---Calculate skewness of mimic amino acid ---")
        input_path = my_args.input
        file_content_df = func.file_type_judge(input_path)
        result = func.mimic_skewness(file_content_df)[0]                # mimic amino acid skewness dataframe
        print("---Skewness of mimic amino acid Completed---")

    # -mp: predict mimic
    if my_args.mimic_predict:
        print("---Predict possible mimic motifs---")
        print(Fore.CYAN + "Note: It is recommended that the number of motifs do not exceed 15; "
                          "otherwise, the memory usage is too large!")
        input_path = my_args.input
        file_content_df = func.file_type_judge(input_path)
        temp_result = func.mimic_skewness(file_content_df)
        new_mimic_df = temp_result[1]           # mimic amino acid dataframe
        mimic_pro_df = temp_result[2]           # mimic amino acid probability dataframe
        result = func.mimic_predict(new_mimic_df, mimic_pro_df)
        print("---Possible mimic motifs Completed---")

    if my_args.mimic_result:

        if my_args.mimic_result.split(".")[1] == "xlsx":
            result = pd.read_excel(my_args.mimic_result, header=0)
        elif my_args.mimic_result.split(".")[1] == "txt":
            result = pd.read_table(my_args.mimic_result, sep="\t", header=0)
        elif my_args.mimic_result.split(".")[1] == "csv":
            result = pd.read_csv(my_args.mimic_result, sep=",", header=0)
        else:
            print("%s Warning: Try to read the file in txt format: %s %s" %
                  (Fore.YELLOW, Fore.BLUE, my_args.mimic_result))
            result = pd.read_table(my_args.mimic_result, sep="\t", header=0)
            print("%s Note: The file was read successfully: %s %s" %
                  (Fore.CYAN, Fore.BLUE, my_args.mimic_result))

    """
    Step 3:Create a sequence annotation file.
    This step is optional.
    """

    if my_args.ncbi_annotation:  # add ncbi annotation
        if os.path.isdir(my_args.ncbi_annotation):
            annotation_content = func.files_combine(func.files_path(my_args.ncbi_annotation, user_sys_os))
            dict_tuple = func.obtain_annotation(annotation_content)
            result = func.add_annotation_to_xlsx(result, dict_tuple)
            print("---Annotation Supplement Completed---")
        elif os.path.isfile(my_args.ncbi_annotation):
            annotation_content = func.files_combine([my_args.ncbi_annotation])
            dict_tuple = func.obtain_annotation(annotation_content)
            result = func.add_annotation_to_xlsx(result, dict_tuple)
            print("---NCBI Annotation Supplement Completed---")
        else:
            print(Fore.RED + "---Error! Please check the annotations file path is correct!---")
            sys.exit()

    if my_args.ncbi_annotation_result:
        if my_args.ncbi_annotation_result.split(".")[1] == "xlsx":
            result = pd.read_excel(my_args.ncbi_annotation_result, header=0)
        elif my_args.ncbi_annotation_result.split(".")[1] == "txt":
            result = pd.read_table(my_args.ncbi_annotation_result, sep="\t", header=0)
        elif my_args.ncbi_annotation_result.split(".")[1] == "csv":
            result = pd.read_csv(my_args.ncbi_annotation_result, sep=",", header=0)
        else:
            print("%s Warning: Try to read the file in txt format: %s %s" %
                  (Fore.YELLOW, Fore.BLUE, my_args.ncbi_annotation_result))
            result = pd.read_csv(my_args.ncbi_annotation_result, sep="\t", header=0)
            print("%s Note: The file was read successfully: %s %s" %
                  (Fore.CYAN, Fore.BLUE, my_args.ncbi_annotation_result))

    if my_args.ictv_annotation:  # add ictv annotation

        result = func.add_host_source(my_args.ictv_annotation, result)
        print("---ICTV Annotation Supplement Completed---")

    """
    Step 4:Statistics analysis.
    This step is optional.
    """

    # -mc: mimic count
    if my_args.mimic_count:
        file_path = my_args.input
        file_content_df = func.file_type_judge(file_path)

        col_list = list(file_content_df.columns)
        if "Mimic" not in col_list:
            print(Fore.RED + "---Error! Please change the column name of mimic to 'Mimic'!---")
            sys.exit()
        elif "Mimic" in col_list:
            print("---The mimic statistics(partial) are as follows---")
            result = func.mimic_stat(file_content_df)
            print("---Mimic count completed---")

    # -cc: classification count
    if my_args.classification_count:
        file_path = my_args.input
        file_content_df = func.file_type_judge(file_path)

        if my_args.class_type.upper() == "FAMILY":
            class_type = "FAMILY"
            print("---The classification statistics(partial) are as follows---")
            result = func.class_count(file_content_df, class_type)
            print("---The classification count completed---")
        elif my_args.class_type.upper() in ("REALM", "KINGDOM", "PHYLUM", "CLASS", "ORDER"):
            class_type = my_args.class_type.upper()
            print("---The classification statistics(partial) are as follows---")
            result = func.class_count(file_content_df, class_type)
            print("---The classification count completed---")
        else:
            print(Fore.RED + "ERROR! Please check whether the category name in the form or the entered "
                  "category name meets the requirements!")
            print(Fore.RED + "ERROR! The category name you input is '%s', please change to one of the following: "
                             "'realm', 'kingdom', 'phylum', 'class', 'order', 'family'." % my_args.class_type)
            sys.exit()

    # -hc: host count
    if my_args.host_count:
        file_path = my_args.input
        file_content_df = func.file_type_judge(file_path)

        col_list = list(file_content_df.columns)
        if "Host_Source" in col_list:
            print("---The host statistics(partial) are as follows---")
            result = func.host_count(file_content_df)
            print("---Hosts count completed---")
        elif "Host_Source" not in col_list:
            print(Fore.RED + """ERROR! Please check if the name of the host information column is 'Host_Source', 
            if not, please change it to 'Host_Source' .""")
            sys.exit()

    # Output final result
    if my_args.output:  # output path
        my_args.output = os.path.realpath(my_args.output)
        if my_args.output == "":  # blank path
            print(Fore.RED + "---Error! The output path is not specified!---")
            sys.exit()

        elif os.path.exists(my_args.output):
            if os.path.isdir(my_args.output):
                # txt file in default
                if my_args.output_type.upper() == "TXT":
                    if my_args.output_name == "Result":
                        result_path = my_args.output + os.sep + "Result.txt"
                    else:
                        result_path = my_args.output + os.sep + my_args.output_name + ".txt"
                    result.to_csv(result_path, index=False, sep="\t")
                    print("---Results stored in the path: %s %s %s---" % (Fore.BLUE, result_path, Fore.RESET))
                # xlsx file
                elif my_args.output_type.upper() == "XLSX":
                    if my_args.output_name == "Result":
                        result_path = my_args.output + os.sep + "Result.xlsx"
                    else:
                        result_path = my_args.output + os.sep + my_args.output_name + ".xlsx"
                    result.to_excel(result_path, index=False)
                    print("---Results stored in the path: %s %s %s---" % (Fore.BLUE, result_path, Fore.RESET))
                # csv file
                elif my_args.output_type.upper() == "CSV":
                    if my_args.output_name == "Result":
                        result_path = my_args.output + os.sep + "Result.csv"
                    else:
                        result_path = my_args.output + os.sep + my_args.output_name + ".csv"
                    result.to_csv(result_path, index=False, sep=",")
                    print("---Results stored in the path: %s %s %s---" % (Fore.BLUE, result_path, Fore.RESET))
                else:
                    print("{}Warning: This format '{}{}{}' is not currently supported for writing, "
                          "the default format has been used.".format(Fore.YELLOW, Fore.RED, my_args.output_type,
                                                                     Fore.YELLOW))
                    print("{}Warning: Results are saved in 'txt' format, recommend use notepad to open after "
                          "modifying the file suffix to 'txt'.".format(Fore.YELLOW))
                    if my_args.output_name == "Result":
                        result_path = my_args.output + os.sep + "Result.txt"
                    else:
                        result_path = my_args.output + os.sep + my_args.output_name + ".txt"
                    result.to_csv(result_path, index=False, sep="\t")
                    print("---Results stored in the path: %s %s %s---" % (Fore.BLUE, result_path, Fore.RESET))

            elif os.path.isfile(my_args.output):  # a specified file as output
                print("{}Error! The output path is a file path, please pass in a folder path.".format(Fore.RED))
                sys.exit()
        # invalid path
        elif not os.path.exists(my_args.output):
            # create a directory
            os.makedirs(my_args.output)
            print("%sWarning: The current path does not exist, automatically create the path: %s %s" %
                  (Fore.YELLOW, Fore.BLUE, os.path.realpath(my_args.output)))

            # txt file
            if my_args.output_type.upper() == "TXT":
                if my_args.output_name == "Result":
                    result_path = my_args.output + os.sep + "Result.txt"
                else:
                    result_path = my_args.output + os.sep + my_args.output_name + ".txt"
                result.to_csv(result_path, index=False, sep="\t")
                print("---Results stored in the path: %s %s %s---" %
                      (Fore.BLUE, os.path.realpath(result_path), Fore.RESET))
            # xlsx file
            elif my_args.output_type.upper() == "XLSX":
                if my_args.output_name == "Result":
                    result_path = my_args.output + os.sep + "Result.xlsx"
                else:
                    result_path = my_args.output + os.sep + my_args.output_name + ".xlsx"
                result.to_excel(result_path, index=False)
                print("---Results stored in the path: %s %s %s---" %
                      (Fore.BLUE, os.path.realpath(result_path), Fore.RESET))
            # csv file
            elif my_args.output_type.upper() == "CSV":
                if my_args.output_name == "Result":
                    result_path = my_args.output + os.sep + "Result.csv"
                else:
                    result_path = my_args.output + os.sep + my_args.output_name + ".csv"
                result.to_csv(result_path, index=False, sep=",")
                print("---Results stored in the path: %s %s %s---" %
                      (Fore.BLUE, os.path.realpath(result_path), Fore.RESET))
            else:
                print("{}Warning: This output format '{}{}{}' is not currently supported for writing in, "
                      "the default format has been used.{}".format(Fore.YELLOW, Fore.RED, my_args.output_type,
                                                                   Fore.YELLOW, Fore.RESET))
                print("{}Warning: Result is saved in 'txt' format, recommend use notepad to open it.{}"
                      .format(Fore.YELLOW, Fore.RESET))
                if my_args.output_name == "Result":
                    result_path = my_args.output + os.sep + "Result.txt"
                else:
                    result_path = my_args.output + os.sep + my_args.output_name + ".txt"
                result.to_csv(result_path, index=False, sep="\t")
                print("---Results stored in the path: %s %s %s---" %
                      (Fore.BLUE, os.path.realpath(result_path), Fore.RESET))

    # elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = float(elapsed_time % 60)
    print("---Take %s%02d:%02d:%02.2f%s time in total---" % (Fore.LIGHTYELLOW_EX, hours, minutes, seconds, Fore.RESET))


if __name__ == "__main__":
    init(wrap=True, autoreset=True)
    starts()
