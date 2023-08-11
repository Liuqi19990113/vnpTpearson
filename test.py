files_list = [i for i in range(0,240)]
cut_number = 8
file_number_bin = int(len(files_list)/cut_number)
new_files_lists = [files_list[i:i+file_number_bin] for i in range(0, len(files_list), file_number_bin)]
print(new_files_lists)