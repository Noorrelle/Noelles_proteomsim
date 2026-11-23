from pyteomics import achrom

def predict_lc_retention_times(peptides):
    """
    Calculate the predicted retention time for each peptide.

    Iterates over a list of peptides and calculates the retention time
    using the AChrom algorithm with the Guo pH 7.0 retention coefficients.

    Parameters
    ----------
    peptides : list, str
        A list of peptide sequences for which retention times are calculated.

    Returns
    -------
    dict
        A dictionary mapping each peptide sequence to its predicted
        retention time.
    """
    
    retentiontimes_dic = {}
    if type(peptides)==str:
        retentiontimes_dic = round(achrom.calculate_RT(peptides, achrom.RCs_guo_ph7_0),2)
    else:
        for pep in peptides: 
            retentiontimes_dic[pep]= round(achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0),2)
    return(retentiontimes_dic)

import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    Plot the distribution of peptide retention times as a histogram.

    Creates a histogram to visualize the frequency distribution of
    retention times, with a configurable number of bins.

    Parameters
    ----------
    retention_times : list or array
        A collection of peptide retention times (in minutes).
    resolution : int, optional
        Number of bins used to divide the retention time range.
        Default is 30.

    Returns
    -------
    None
        Displays a histogram with retention time (min) on the x-axis
        and frequency on the y-axis.
    """
    plt.hist(retention_times, resolution, color="#B53389")
    plt.xlabel('Retention time in min')
    plt.ylabel('Frequency') 

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Select peptides whose retention times fall within a specified range.

    Parameters
    ----------
    peptide_rt_map : dict
        Dictionary mapping peptide sequences to their retention times.
    lower_ret_time : float
        Lower boundary of the retention time window (in minutes).
    upper_ret_time : float
        Upper boundary of the retention time window (in minutes).

    Returns
    -------
    dict
        Dictionary containing only peptides with retention times
        within the specified window (inclusive).
    """
    FilteredRetentiontimes_dic = {}
    if len(peptide_rt_map) == 1:
        ret = list(peptide_rt_map.values())[0]
        if ret >= lower_ret_time and ret <= upper_ret_time:
            FilteredRetentiontimes_dic=peptide_rt_map
    else:
        for pep,ret in peptide_rt_map.items():
            if ret >= lower_ret_time and ret <= upper_ret_time:
                FilteredRetentiontimes_dic[pep]=ret
    return(FilteredRetentiontimes_dic)  
            
   