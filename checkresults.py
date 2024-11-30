# get all txt files in folder 'result'
import os

result_folder = 'C:\\Users\\cybcat\\OneDrive\\文档\\GitHub\\PFTool\\results'
txt_files = [os.path.join(result_folder, f) for f in os.listdir(result_folder) if f.endswith('.txt')]

# read each txt file and print the content
for txt_file in txt_files:
    with open(txt_file, 'r') as f:
        content = f.readlines()
        if len(content) == 3:
            # read a stringed list from the third line
            lst = eval(content[2])
            if type(lst) == list:
                if sum(abs(_) for _ in lst[0]) < 1e5:
                    print(txt_file, 'is a valid result')

dic = dict()
for txt_file in txt_files:
    with open(txt_file, 'r') as f:
        content = f.readlines()
        if len(content) == 3:
            # read a stringed list from the third line
            lst = eval(content[1])
            if type(lst) == list and content[2] == "None":
                if sum(abs(_) for _ in lst[0]) < 1e5 and min(abs(_) for _ in lst[0]) < 1e3 and sum(abs(_) for _ in lst[1]) < 1e5 and min(abs(_) for _ in lst[1]) < 1e3:
                    pair = ('-' in txt_file, txt_file.replace('_', '__')[50:53])
                    if pair not in dic:
                        dic[pair] = [txt_file]
                    else:
                        dic[pair].append(txt_file)


for _ in dic:
    if len(dic[_]) <= 2:
        print(dic[_], 'is a valid result')