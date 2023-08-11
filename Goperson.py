import sys
import os
import centrality
import dataget
import cn2andcn4
import numpy as np
import personcov

#Variables definition
centrality_interval = np.array([0, 5, 10, 20, 30, 40, 50])
#centrality_interval = np.linspace(0,50,100)

def mainprocess(files_list,cut_number):
    new_files_lists = np.array_split(files_list,cut_number)
    print([len(l) for l in new_files_lists])
    i = 0
    cov_v2pt_std,cov_v2pt_2sub,cov_v2pt_3sub = np.zeros((cut_number,len(centrality_interval))),np.zeros((cut_number,len(centrality_interval))),np.zeros((cut_number,len(centrality_interval)))
    cov_v3pt_std,cov_v3pt_2sub,cov_v3pt_3sub = np.zeros((cut_number,len(centrality_interval))),np.zeros((cut_number,len(centrality_interval))),np.zeros((cut_number,len(centrality_interval)))
    per_v2pt_std,per_v2pt_2sub,per_v2pt_3sub = np.zeros((cut_number,len(centrality_interval))),np.zeros((cut_number,len(centrality_interval))),np.zeros((cut_number,len(centrality_interval)))
    per_v3pt_std,per_v3pt_2sub,per_v3pt_3sub = np.zeros((cut_number,len(centrality_interval))),np.zeros((cut_number,len(centrality_interval))),np.zeros((cut_number,len(centrality_interval)))
    for files_lists in new_files_lists:
        j = 0
        mudic = centrality.mult_dic(files_lists)
        centrality_results = centrality.centrality_sort(mudic,centrality_interval)
        for one_cent_bin in centrality_results:
            print('In file group {} and centrality group {}'.format(i+1, j+1))
            para_list_for_cumulants,para_for_person = dataget.readdata(one_cent_bin)
            para_list_for_c22andc24 = [list_item+[2] for list_item in para_list_for_cumulants]
            para_list_for_c32andc34 = [list_item+[3] for list_item in para_list_for_cumulants]
            para_list_for_v2pt_person = [list_item+[2] for list_item in para_for_person]
            para_list_for_v3pt_person = [list_item+[3] for list_item in para_for_person]
            c24, c22sub = cn2andcn4.runflow(para_list_for_c22andc24)
            c34, c32sub = cn2andcn4.runflow(para_list_for_c32andc34)
            cov_v2pt_std[i][j], cov_v2pt_2sub[i][j], cov_v2pt_3sub[i][j], per_v2pt_std[i][j], per_v2pt_2sub[i][j], per_v2pt_3sub[i][j] = personcov.person_go(para_list_for_v2pt_person,c24,c22sub)
            cov_v3pt_std[i][j], cov_v3pt_2sub[i][j], cov_v3pt_3sub[i][j], per_v3pt_std[i][j], per_v3pt_2sub[i][j], per_v3pt_3sub[i][j] = personcov.person_go(para_list_for_v3pt_person,c34,c32sub)
            j+=1
        i+=1
    cov_v2pt_std_fin = np.column_stack((np.mean(cov_v2pt_std, axis=0), np.std(cov_v2pt_std, axis=0)))
    cov_v2pt_2sub_fin = np.column_stack((np.mean(cov_v2pt_2sub, axis=0), np.std(cov_v2pt_2sub, axis=0)))
    cov_v2pt_3sub_fin = np.column_stack((np.mean(cov_v2pt_3sub, axis=0), np.std(cov_v2pt_3sub, axis=0)))
    cov_v3pt_std_fin = np.column_stack((np.mean(cov_v3pt_std, axis=0), np.std(cov_v3pt_std, axis=0)))
    cov_v3pt_2sub_fin = np.column_stack((np.mean(cov_v3pt_2sub, axis=0), np.std(cov_v3pt_2sub, axis=0)))
    cov_v3pt_3sub_fin = np.column_stack((np.mean(cov_v3pt_3sub, axis=0), np.std(cov_v3pt_3sub, axis=0)))
    per_v2pt_std_fin = np.column_stack((np.nanmean(per_v2pt_std, axis=0), np.nanstd(per_v2pt_std, axis=0)))
    per_v2pt_2sub_fin = np.column_stack((np.nanmean(per_v2pt_2sub, axis=0), np.nanstd(per_v2pt_2sub, axis=0)))
    per_v2pt_3sub_fin = np.column_stack((np.nanmean(per_v2pt_3sub, axis=0), np.nanstd(per_v2pt_3sub, axis=0)))
    per_v3pt_std_fin = np.column_stack((np.nanmean(per_v3pt_std, axis=0), np.nanstd(per_v3pt_std, axis=0)))
    per_v3pt_2sub_fin = np.column_stack((np.nanmean(per_v3pt_2sub, axis=0), np.nanstd(per_v3pt_2sub, axis=0)))
    per_v3pt_3sub_fin = np.column_stack((np.nanmean(per_v3pt_3sub, axis=0), np.nanstd(per_v3pt_3sub, axis=0)))



    return cov_v2pt_std_fin,cov_v2pt_2sub_fin,cov_v2pt_3sub_fin,cov_v3pt_std_fin,cov_v3pt_2sub_fin,cov_v3pt_3sub_fin\
            ,per_v2pt_std_fin,per_v2pt_2sub_fin,per_v2pt_3sub_fin,per_v3pt_std_fin,per_v3pt_2sub_fin,per_v3pt_3sub_fin


#go!
if __name__ == "__main__":
    print("go!")
    files_dir_list = sys.argv
    del files_dir_list[0]
    whichtype=(files_dir_list[0].split('_'))[3] 
    whereami = os.getcwd()
    os.chdir(whereami)
    files_list = []
    for dict in files_dir_list:
        abs_dict_path = os.path.abspath(dict)
        hdf5_list = os.listdir(dict)
        for hdf5_file in hdf5_list:
            files_list.append(os.path.join(abs_dict_path,hdf5_file))
    print("Target files: {}".format(files_list))
    cov_v2pt_std,cov_v2pt_2sub,cov_v2pt_3sub,\
    cov_v3pt_std,cov_v3pt_2sub,cov_v3pt_3sub,\
    per_v2pt_std,per_v2pt_2sub,per_v2pt_3sub,\
    per_v3pt_std,per_v3pt_2sub,per_v3pt_3sub = mainprocess(files_list,10)
    print("cov_v2pt_3sub = {}".format(cov_v2pt_3sub))
    print("cov_v3pt_3sub = {}".format(cov_v3pt_3sub))
    print("per_v2pt_3sub = {}".format(per_v2pt_3sub))
    print("per_v3pt_3sub = {}".format(per_v3pt_3sub))
    results_name = whichtype+'_vnpt.npz'
    np.savez(results_name, cov_v2pt_std = cov_v2pt_std, cov_v2pt_2sub=cov_v2pt_2sub, cov_v2pt_3sub=cov_v2pt_3sub, \
             cov_v3pt_std=cov_v3pt_std,cov_v3pt_2sub=cov_v3pt_2sub, cov_v3pt_3sub=cov_v3pt_3sub, \
            per_v2pt_std=per_v2pt_std, per_v2pt_2sub=per_v2pt_2sub,per_v2pt_3sub=per_v2pt_3sub, \
            per_v3pt_std=per_v3pt_std, per_v3pt_2sub= per_v3pt_2sub, per_v3pt_3sub = per_v3pt_3sub)
