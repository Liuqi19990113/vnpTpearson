import numpy as np
from multiprocessing.pool import Pool
from numpy import sqrt

def qn1(phi, n):
    phi = (np.array(phi)).astype(complex)
    numerator = (np.exp(1j*n*phi)).sum()
    denominator = len(phi)
#    print(denominator)
    return numerator/denominator


def on2(phi, pt, n, ptmean):
#    print(phi)
    phi = (np.array(phi)).astype(complex)
#    print(phi)
    pt = np.array(pt)
    numerator = (np.exp(1j*n*phi)*(pt-ptmean)).sum()
    denominator = len(phi)
#    print(denominator)
    return numerator/denominator


def pmk(pt, m, ptmean):
    pt = np.array(pt)
    numerator = (np.power(pt-ptmean, m)).sum()
    denominator = len(pt)
#    print(denominator)
    return numerator/denominator


def tauk(weight, k):
    weight = np.array(weight)
    numerator = (np.power(weight, k+1)).sum()
    denominator = np.power(weight.sum(), k+1)
    return numerator/denominator


def event_cov_std(phi_std, pt_std, n, ptmean):
    event_qn1 = qn1(phi_std, n)
    tau1 = tauk(np.ones(len(phi_std)), 1)
    tau2 = tauk(np.ones(len(phi_std)), 2)
    p11 = pmk(pt_std, 1, ptmean)
    p13 = pmk(pt_std, 1, ptmean)
    event_on2 = on2(phi_std, pt_std, n, ptmean)
    event_std_num = (np.abs(event_qn1)*np.abs(event_qn1)-tau1)*p11-2*tau1*(event_on2*np.conjugate(event_qn1)).real+2*tau2*p13
    event_std_dom = 1 - 3*tau1 + 2*tau2
    event_std = event_std_num/event_std_dom
    return event_std, event_std_dom


def event_cov_2sub(phi_sub1, phi_sub2, pt_sub1, pt_sub2, n, ptmean):
    event_2sub_num_part1 = (qn1(phi_sub1, n)*pmk(pt_sub1, 1, ptmean)-tauk(np.ones(len(phi_sub1)), 1)*on2(phi_sub1,pt_sub1,n,ptmean))*np.conjugate(qn1(phi_sub2, n))
    event_2sub_num_part2 = (qn1(phi_sub2, n)*pmk(pt_sub2, 1, ptmean)-tauk(np.ones(len(phi_sub2)), 1)*on2(phi_sub2,pt_sub2,n,ptmean))*np.conjugate(qn1(phi_sub1, n))
    event_2sub_num = (event_2sub_num_part1 + event_2sub_num_part2).real
    event_2sub_dom = 2 - tauk(np.ones(len(phi_sub1)), 1) - tauk(np.ones(len(phi_sub2)), 1)
#    print(len(phi_sub1))
    event_2sub = event_2sub_num/event_2sub_dom
    return event_2sub, event_2sub_dom


def event_cov_3sub(phi_sub1, phi_sub2, pt_sub3, n, ptmean):
    event_3sub = (qn1(phi_sub1, n)*np.conjugate(qn1(phi_sub2, n))).real*pmk(pt_sub3, 1, ptmean)
    return event_3sub, 1


def event_varpt(pt,ptmean):
    p11 = pmk(pt, 1,ptmean)
    p22 = pmk(pt, 2,ptmean)
    tau1 = tauk(np.ones(len(pt)), 1)
    return (p11*p11-tau1*p22)/(1-tau1), 1-tau1

def meanpT(pt_list):
    return np.mean(pt_list)


def varvn(cn4std, cn22sub):
    return cn4std + cn22sub*cn22sub


def event_cov_process(para_list):
    phi_std = para_list[0]
    phi_sub1 = para_list[1]
    phi_sub2 = para_list[2]
    phi_sub3 = para_list[3]
    pt_std = para_list[4]
    pt_sub1 = para_list[5]
    pt_sub2 = para_list[6]
    pt_sub3 = para_list[7]
    n = para_list[8]
    ptmean = para_list[9]
    event_std, event_std_w = event_cov_std(phi_std, pt_std, n,ptmean)
#    print(event_std)
    event_2sub, event_2sub_w = event_cov_2sub(phi_sub1, phi_sub2, pt_sub1, pt_sub2, n,ptmean)
    event_3sub, event_3sub_w = event_cov_3sub(phi_sub1, phi_sub2, pt_sub3, n,ptmean)
    event_var_pt, event_var_pt_w = event_varpt(pt_std,ptmean)
    return event_std, event_std_w, event_2sub, event_2sub_w, event_3sub, event_3sub_w, event_var_pt, event_var_pt_w

def weight_mean(value_list,weight_list):
    return (value_list*weight_list).sum()/weight_list.sum()


def person_go(all_event_list, cn4std, cn22sub):
    all_pt_list = [ll[4] for ll in all_event_list]
    with Pool() as pool:
        all_pt_mean = np.array(pool.map(meanpT,all_pt_list))
    all_pt_weight = np.array([len(event_pt) for event_pt in all_pt_list])
    pt_weight_mean = weight_mean(all_pt_mean,all_pt_weight)
    new_all_event_list = [listitem+[pt_weight_mean] for listitem in all_event_list]
#    print(new_all_event_list)
    with Pool() as pool:
        all_cov_list = pool.map(event_cov_process, new_all_event_list)
    cov_std_list = np.array([event_cov[0] for event_cov in all_cov_list])
    cov_std_w_list = np.array([event_cov[1] for event_cov in all_cov_list])
#    print(cov_std_w_list)
    cov_2sub_list = np.array([event_cov[2] for event_cov in all_cov_list])
    cov_2sub_w_list = np.array([event_cov[3] for event_cov in all_cov_list])
    cov_3sub_list = np.array([event_cov[4] for event_cov in all_cov_list])
    cov_3sub_w_list = np.array([event_cov[5] for event_cov in all_cov_list])
    var_pt_list = np.array([event_cov[6] for event_cov in all_cov_list])
#    print(var_pt_list)
    var_pt_w_list = np.array([event_cov[7] for event_cov in all_cov_list])
    fin_cov_std = weight_mean(cov_std_list,cov_std_w_list)
    fin_cov_2sub = weight_mean(cov_2sub_list,cov_2sub_w_list)
    fin_cov_3sub = weight_mean(cov_3sub_list,cov_3sub_w_list)
    print('This cov_std = {}'.format(fin_cov_std))
    print('This cov_2sub = {}'.format(fin_cov_2sub))
    print('This cov_3sub = {}'.format(fin_cov_3sub))
    fin_varpt = weight_mean(var_pt_list,var_pt_w_list)
    fin_varvn = varvn(cn4std, cn22sub)
    print('This varpt and varvn = {}, {}'.format(fin_varpt,fin_varvn))
    fin_per_std = fin_cov_std/(sqrt(fin_varpt)*sqrt(fin_varvn))
    fin_per_2sub = fin_cov_2sub/(sqrt(fin_varpt)*sqrt(fin_varvn))
    fin_per_3sub = fin_cov_3sub/(sqrt(fin_varpt)*sqrt(fin_varvn))
    print('This person_std = {}'.format(fin_per_std))
    print('This person_2sub = {}'.format(fin_per_2sub))
    print('This person_3sub = {}'.format(fin_per_3sub))
    return fin_cov_std, fin_cov_2sub, fin_cov_3sub, fin_per_std, fin_per_2sub, fin_per_3sub

    



