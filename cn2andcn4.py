import numpy as np
from multiprocessing.pool import Pool

def qn(phi, n) -> complex:
    phi = (np.array([phi])).astype(complex)
    qn_vector = np.exp(1j*n*phi).sum()
    return qn_vector

def eve_weight(mult, n) -> float:
    weight = 1.
    for i in range(mult - n + 1, mult + 1):
        i_float = float(i)
        weight = i_float*weight
    return weight

def eve_weight_etagap(mult1, mult2) -> float:
    return mult1*mult2

def single_ave2(qn, mult) -> float:
    qn_conj = np.conjugate(qn)
    ave2_single_event = (((qn*qn_conj).real - mult)/(eve_weight(mult, 2)))
    return ave2_single_event

def single_ave2_with_eta_gap(qn_1,qn_2,mult_1,mult_2) -> float:
    qn_2_conj = np.conjugate(qn_2)
    return float((qn_1*qn_2_conj).real)/float((mult_1*mult_2))

def single_ave4(qn,q2n,mult) -> float:
    qn_conj = np.conjugate(qn)
    q2n_conj = np.conjugate(q2n)
    abs_qn_fp = (qn*qn*qn_conj*qn_conj).real
    abs_q2n_sq = (q2n_conj*q2n).real
    single_ave4_numerator = (abs_qn_fp + abs_q2n_sq - 2*(q2n*qn_conj*qn_conj).real - 
                4*(mult - 2)*(qn*qn_conj).real + 2*mult*(mult - 3))
    single_ave4_denominator = eve_weight(mult, 4)
    ave4_single_event = single_ave4_numerator/single_ave4_denominator
    return ave4_single_event

def all_ave(weight_list: list, ave_list: list) -> float:
    '''This function is used to calculate <<2>> and <<4>>
    '''
    all_event_ave_numerator = (weight_list*ave_list).sum()
    all_event_ave_denominator = weight_list.sum()
    ave_all_events = all_event_ave_numerator/all_event_ave_denominator
    return ave_all_events

def cn2(ave2_array: list, weight2_array: list) -> tuple:
    '''This function is used to calculate cn2.
    '''
    ave2_array = np.array([ave2_array])
    weight2_array = np.array([weight2_array])
    cn2 = all_ave(weight2_array, ave2_array)
    return cn2

def cn4(ave2_array: list, ave4_array: list, weight2_array: list, weight4_array: list) ->tuple:
    '''This function is used to calculate intvn4 and its error.
    '''
    ave2_array = np.array([ave2_array])
    weight2_array = np.array([weight2_array])
    ave4_array = np.array([ave4_array])
    weight4_array = np.array([weight4_array])
    all_ave2 = all_ave(weight2_array, ave2_array)
    all_ave4 = all_ave(weight4_array, ave4_array)
    cn4 = all_ave4 - 2*all_ave2*all_ave2
    return cn4

def cn2tsub(ave2_etagap_array: list, weight2_etagap_array: list) -> tuple:
    '''This function is used to calculate intvn2 with etagap and its error.
    '''
    ave2_etagap_array = np.array([ave2_etagap_array])
    weight2_etagap_array = np.array([weight2_etagap_array])
    cn2_etagap = all_ave(weight2_etagap_array, ave2_etagap_array)
    return cn2_etagap

def sigleeveprocesser(para_list):
    '''return <2> <4> <2>withsub for an event'''
    phi_std, phi_sub1, phi_sub2 = para_list[0], para_list[1], para_list[2]
    n = para_list[3]
    mult_std, mult_1, mult_2 = len(phi_std), len(phi_sub1), len(phi_sub2)
    qn_std, q2n_std, qn_1, qn_2 = qn(phi_std, n), qn(phi_std, 2*n), qn(phi_sub1, n), qn(phi_sub2, n)
    sinave2 = single_ave2(qn_std, mult_std)
    weight2 = eve_weight(mult_std, 2)
    sinave4 = single_ave4(qn_std, q2n_std, mult_std)
    weight4 = eve_weight(mult_std, 4)
    ave2etagap = single_ave2_with_eta_gap(qn_1, qn_2, mult_1, mult_2)
    weight2withgap = eve_weight_etagap(mult_1, mult_2)
    return sinave2, sinave4, ave2etagap, weight2, weight4, weight2withgap

def runflow(eberesults):
    with Pool() as pool:
        single_eve_list = pool.map(sigleeveprocesser, eberesults)
    ave2_array = [eve_results[0] for eve_results in single_eve_list]
    ave4_array = [eve_results[1] for eve_results in single_eve_list]
    ave2etagap_array = [eve_results[2] for eve_results in single_eve_list]
    weight2_array = [eve_results[3] for eve_results in single_eve_list]
    weight4_array = [eve_results[4] for eve_results in single_eve_list]
    weight2etagap_array = [eve_results[5] for eve_results in single_eve_list]
#    resultcn2 = cn2(ave2_array, weight2_array)
    resultcn4 = cn4(ave2_array,ave4_array,weight2_array,weight4_array)
    resultcn2sub = cn2tsub(ave2etagap_array,weight2etagap_array)
#    return resultcn2, resultcn4, resultcn2sub
    return resultcn4, resultcn2sub




