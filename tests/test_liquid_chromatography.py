from Noelles_proteomsim import (
    predict_lc_retention_times,
    select_retention_time_window,
)
from pyteomics import achrom

def test_predict_lc_retention_times():
    peptides = 'GSH'
    expected1 = round(achrom.calculate_RT("GSH", achrom.RCs_guo_ph7_0, raise_no_mod=False),2)

    actual = predict_lc_retention_times(peptides)
    
    assert actual == expected1

test_predict_lc_retention_times()

def test_select_retention_time_window():
    peptide_rt_map = {"A":5, "B":10, "C":12, "D":9}
    peptide_rt_map2 = {"A":5}
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=5, upper_ret_time=11)
    selected2 = select_retention_time_window(peptide_rt_map2, lower_ret_time=5, upper_ret_time=11)

    assert len(selected) == 3
    assert len(selected2) == 1

test_select_retention_time_window()
#more room to test, what happens if none of them is in and also tets which peptides are inside 