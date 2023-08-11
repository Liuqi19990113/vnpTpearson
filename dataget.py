import numpy as np
import h5py
eta_std = [-2.5, 2.5]
eta_sub1 = [-2.5,-0.75]
eta_sub2 = [0.75, 2.5]
eta_sub3 = [-0.5,0.5]
y_cut = [-0.1, 0.1]
pt_cut = [0.5, 5]

def readdata(hdf5file_event_list: list) -> list:
    para_for_cumu = []
    para_for_person = []
    for hdf5file_event in hdf5file_event_list:
        hdf5_path, event_name = hdf5file_event.split('*')
        with h5py.File(hdf5_path, 'r') as f:
            nsample = f[event_name]['sample'][-1]
            sample = f[event_name]['sample']
            phi = f[event_name]['phi']
            charge = f[event_name]['charge']
            eta = f[event_name]['eta']
            pt = f[event_name]['pT']
            y = f[event_name]['y']
            for i in range(1, nsample+1): 
                phi_std = phi[(sample == i) & (charge != 0) & (eta > eta_std[0]) & 
                            (eta < eta_std[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]  #You can change the condition.
                phi_sub1 = phi[(sample == i) & (charge != 0) & (eta > eta_sub1[0]) & 
                            (eta < eta_sub1[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]  #You can change the condition.
                phi_sub2 = phi[(sample == i) & (charge != 0) & (eta > eta_sub2[0]) & 
                            (eta < eta_sub2[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]  #You can change the condition.
                phi_sub3 = phi[(sample == i) & (charge != 0) & (eta > eta_sub3[0]) & 
                            (eta < eta_sub3[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]  #You can change the condition.
                if (len(phi_std) <5 or len(phi_sub1) <5 or len(phi_sub2) <5 or len(phi_sub3) <5):
                    continue
#                pt_all = pt[(sample == i) & (charge != 0) & (y > y_cut[0]) & (y < y_cut[1])]
                pt_in_std = pt[(sample == i) & (charge != 0)  & (eta > eta_std[0]) & (eta < eta_std[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]
                pt_in_sub1 = pt[(sample == i) & (charge != 0)  & (eta > eta_sub1[0]) & (eta < eta_sub1[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]
                pt_in_sub2 = pt[(sample == i) & (charge != 0)  & (eta > eta_sub2[0]) & (eta < eta_sub2[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]
                pt_in_sub3 = pt[(sample == i) & (charge != 0)  & (eta > eta_sub3[0]) & (eta < eta_sub3[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]
                para_for_cumu.append([phi_std,phi_sub1,phi_sub2])
                para_for_person.append([phi_std,phi_sub1,phi_sub2,phi_sub3,pt_in_std,pt_in_sub1,pt_in_sub2,pt_in_sub3])
    return para_for_cumu, para_for_person
#    return para_list_for_cumulants


