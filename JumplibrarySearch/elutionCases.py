import pandas as pd
import numpy as np


def get_max_int_rt(rt_list, rt_int_dict):
    #first reverse dictionary to rt_int_dict
#     dict1_rev = dict((v,k) for k,v in int_rt_dict.items())
    max_int = 0
    reqd_rt = -1
    if rt_list != -1:
        for index, x in enumerate(rt_list):
    #         intensity = dict1_rev[x]
            intensity = rt_int_dict[x]
            if intensity > max_int:
                max_int = intensity
                reqd_rt = rt_list[index]
    return reqd_rt


def getMaxIntCluster(list_of_clusters, rt_int_dict):
    max_int_list = []
    final_rt_list = []
    final_rt_val = 0
    for lst in list_of_clusters:
        max_int = 0
        for val in lst:
            intensity = rt_int_dict[val]
            if intensity > max_int:
                max_int = intensity
                final_rt_val = val
        max_int_list.append(max_int)
        final_rt_list.append(final_rt_val)
    
    final_rt_max_int_dict = int_rt_dict(max_int_list, final_rt_list)
    
#     selected_rt = final_rt_max_int_dict[np.max(list(final_rt_max_int_dict.keys()))]
    
    return final_rt_max_int_dict



#lst1 = rt_list; lst2 = intensity_list
def int_rt_dict(lst1, lst2):
    dict1 = dict(zip(lst1, lst2))
    return dict1






def get_multipeakcluster(row, eps):
#RT_Clust_eps_1 -- for final_RT_singleton == -1, see if all RT_Clust_eps_1 has 2 or more peaks (Case 2), else Case3
    singletonOnly = row.final_RT_singleton
    list_of_cluster = row.RT_Clust_eps_1
    if len(max_int_rt_dict.keys()) == 1:
        #found the singleton cluster
        return list(max_int_rt_dict.values())[0]
    else:
        return -1


'''
CASE1 all final_RT_singleton == -1 are multiclusters

3 more cases of multiclusters
CASE2 = multipeak clusters [MC] --> The clusters have multiple peaks
CASE4 = single peak clusters/ sigleton cluster [SC]
CASE3 = [MC]+[SC]
'''

#get Case 2, 3, 4
#get multipeak cluster (Case 2)



# evaluate the RT cluster peaks
def evalute_rt_cluster(row,column):
    result = "-1"
    list_of_clusters = row[column]
    clust_Type = [] # if 0 then it refers to single peak else 1 -- multiple peaks
    new_list = []
    if len(list_of_clusters) == 1:
        new_list = list_of_clusters
        clust_Type.append(-1)
        result = "Case1"
    else:
        for lst in list_of_clusters:
            if len(lst) > 1:
                clust_Type.append(1)
                new_list.append(lst)
            else:
                clust_Type.append(0)
    
    if clust_Type.count(1) == len(list_of_clusters):
        #CASE2 = multipeak clusters [MC] --> The clusters have multiple peaks
        result = "Case2" 
    
    
    if clust_Type.count(0) == len(list_of_clusters):
        #CASE4 = single peak clusters/ sigleton cluster [SC]
        result = "Case4" 
    
    if result == "-1":
        result = "Case3" 
    
    
    if new_list == []:
        return pd.Series([-1, result])
    else:
        return pd.Series([new_list, result])





#evaluate case2 first

def inferRT_Case2(row, clusterExplore, column1, column2): #clusterExplore which clustType you want to check#column1 is which you want to retain if the condition does not meet # column 2 is what you want to analyze
#     final_RT_singleton = row["final_RT_singleton"]
    final_RT_singleton = row[column1]
    rt_int_dict = row["rt_int_dict"]
    clustType = row["clusterType"]
    rt_clusters = row[column2]
#     rt_clusters = row["RT_peaks_evaluate_eps{}".format(eps)]
    
    if clustType == clusterExplore:
        select_first_rt_cluster =  rt_clusters[0]
        rt = get_max_int_rt(select_first_rt_cluster, rt_int_dict)
        return rt
    else:
        return final_RT_singleton



#evaluate case4 

def inferRT_Case4(row, clusterExplore, column1, column2): #clusterExplore which clustType you want to check#column1 is which you want to retain if the condition does not meet # column 2 is what you want to analyze
#     final_RT_singleton = row["final_RT_singleton"]
    final_RT_singleton = row[column1]
    max_int_rt_dict = row["max_int_rt_dict"]
    clustType = row["clusterType"]
    rt_clusters = row[column2]
#     rt_clusters = row["RT_peaks_evaluate_eps{}".format(eps)]
    
    if clustType == clusterExplore:
#         print (max_int_rt_dict)
        rt = max_int_rt_dict[np.max(list(max_int_rt_dict.keys()))]
        return rt
    else:
        return final_RT_singleton


# evaluate Case 3 for further mini Case2 and mini Case1
# evaluate the RT cluster peaks
def evalute_rt_cluster_case3(row,column):
    case = row["clusterType"]
    result = "-1"
    list_of_clusters = row[column]
    
    if case == "Case3":
        clust_Type = [] # if 0 then it refers to single peak else 1 -- multiple peaks
        new_list = []
        if len(list_of_clusters) == 1:
            new_list = list_of_clusters
            clust_Type.append(-1)
            result = "subCase1"
        else:
            for lst in list_of_clusters:
                if len(lst) > 1:
                    clust_Type.append(1)
                    new_list.append(lst)
                else:
                    clust_Type.append(0)

        if clust_Type.count(1) == len(list_of_clusters):
            #CASE2 = multipeak clusters [MC] --> The clusters have multiple peaks
            result = "subCase2" 


        if clust_Type.count(0) == len(list_of_clusters):
            #CASE4 = single peak clusters/ sigleton cluster [SC]
            result = "subCase4" 

        if result == "-1":
            result = "subCase3" 

    else:
        result = "No sub case"
    return result



def inferRT_case3_subcase2(row, eps):
    final_RT_singleton = row["final_RT_case4"]
    rt_int_dict = row["rt_int_dict"]
    subclustType = row["subClusterTypeCase3"]
    rt_clusters = row["RT_peaks_evaluate_eps{}".format(eps)]
    
    if subclustType == "subCase2":
        select_first_rt_cluster =  rt_clusters[0]
        rt = get_max_int_rt(select_first_rt_cluster, rt_int_dict)
        return rt
    else:
        return final_RT_singleton


def inferRT_case3_subcase1(row, eps):
    final_RT_singleton = row["final_RT_case3_subcase2"]
    rt_int_dict = row["rt_int_dict"]
    subclustType = row["subClusterTypeCase3"]
    rt_clusters = row["RT_peaks_evaluate_eps{}".format(eps)]
    
    if subclustType == "subCase1":
        final_rt_max_int_dict = getMaxIntCluster(rt_clusters, rt_int_dict)
        rt = final_rt_max_int_dict[max(final_rt_max_int_dict.keys())]
        return rt
    else:
        return final_RT_singleton

############################################ WEIGHTED RT FUNCTIONS ###########################

def inferRT_case3_subcase2_wtrt(row, eps):
    final_RT_singleton = row["final_RT_case4"]
    rt_int_dict = row["rt_int_dict"]
    subclustType = row["subClusterTypeCase3"]
    rt_clusters = row["RT_peaks_evaluate_eps{}".format(eps)]
    
    if subclustType == "subCase2":
        select_first_rt_cluster =  [rt_clusters[0]]
        rt = weighted_average_each_cluster(select_first_rt_cluster, rt_int_dict)
        return rt[0]
    else:
        return final_RT_singleton



def inferRT_Case2_wtrt(row, clusterExplore, column1, column2): #clusterExplore which clustType you want to check#column1 is which you want to retain if the condition does not meet # column 2 is what you want to analyze
#     final_RT_singleton = row["final_RT_singleton"]
    final_RT_singleton = row[column1]
    rt_int_dict = row["rt_int_dict"]
    clustType = row["clusterType"]
    rt_clusters = row[column2]
#     rt_clusters = row["RT_peaks_evaluate_eps{}".format(eps)]
    
    if clustType == clusterExplore:
        select_first_rt_cluster =  [rt_clusters[0]]
        rt = weighted_average_each_cluster(select_first_rt_cluster, rt_int_dict)
        return rt[0]
    else:
        return final_RT_singleton


def select_singleton_cluster_wtrt(row):
    rt_list = row.weighted_rt_list
    if len(rt_list) == 1:
        #found the singleton cluster
        return rt_list[0]
    else:
        return -1



def weighted_average_each_cluster(list_of_clusters, rt_int_dict):

    weighted_rt_list = []
    for lst in list_of_clusters:
        num = 0
        den = 0
        for val in lst:
            intensity = rt_int_dict[val]
            num = num+(val*intensity)
            den += intensity
        weightedRT = num/den
        weighted_rt_list.append(weightedRT)
    return weighted_rt_list