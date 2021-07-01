from glob import glob
import itertools
import math
import numpy as np
import os
import pandas as pd
import pickle
import re
import statistics


#region global vars
lipidclass_dict =  {
    'Fatty acyls': ['FA', 'NAGly', 'NAGlySer', 'NAOrn', 'NAE', 'CAR', 'FAHFA'],
    'Glycerolipids': ['DG', 'EtherDG', 'DGDG', 'EtherDGDG', 'MGDG', 'EtherMGDG', 
                      'SQDG', 'EtherSMGDG', 'MG', 'ADGGA', 'DGCC', 'DGGA', 
                      'DGTS/A', 'LDGCC', 'LDGTS/A', 'EtherTG', 'TG'],
    'Glycerophospholipids': ['LPA', 'PA', 'EtherLPC', 'EtherPC', 'LPC', 'PC', 
                             'EtherLPE', 'EtherPE', 'EtherPE(P)', 'PlasmPE', 
                             'LNAPE', 'LPE', 'PE', 'BMP', 'EtherLPG', 'EtherPG', 
                             'HBMP', 'LPG', 'PG', 'CL', 'DLCL', 'MLCL',
                             'Ac2PIM1', 'Ac2PIM2', 'Ac3PIM2', 'Ac4PIM2', 
                             'EtherPI', 'LPI', 'PI', 'EtherPS', 'LNAPS', 'LPS', 
                             'PS', 'PEtOH', 'PMeOH', 'EtherOxPE', 'OxPC', 'OxPE', 
                             'OxPG', 'OxPI', 'OxPS'],
    'Prenol lipids': ['VAE', 'CoQ', 'Vitamine E'],
    'Saccharolipids': ['LipidA'], 
    'Sphingolipids': ['GM3', 'SHexCer', 'SHexCer+O', 'Cer_ADS', 'Cer_AP', 
                      'Cer_AS', 'Cer_BDS', 'Cer_BS', 'Cer_HDS', 'Cer_HS', 
                      'Cer_EBDS', 'Cer_EODS', 'Cer_EOS', 'Cer_NDS', 'Cer_NP', 
                      'Cer_NS', 'CerP', 'AHexCer', 
                      'HexCer_ADS', 'HexCer_AP', 'HexCer_AS', 'HexCer_BDS', 
                      'HexCer_BS', 'HexCer_HDS', 'HexCer_HS', 'HexCer_EOS', 
                      'HexCer_NDS', 'HexCer_NP', 'HexCer_NS', 
                      'Hex2Cer', 'Hex3Cer', 'ASM', 'PE-Cer', 'PE-Cer+O', 
                      'PI-Cer', 'SM', 'SM+O', 
                      'PhytoSph', 'SL', 'SL+O', 'DHSph', 'Sph'],
    'Sterol lipids': ['CASulfate', 'CA', 'DCAE', 'GDCAE', 'GLCAE', 'TDCAE', 
                      'TLCAE', 'AHexCAS', 'AHexCS', 'AHexSIS', 'AHexBRS', 
                      'AHexSTS', 'Vitamine D', 'SSulfate', 'BRSE', 'CASE', 'CE', 
                      'Cholesterol', 'SHex', 'SISE', 'STSE', 'SPE', 'BAHex', 
                      'BASulfate', 'SPEHex', 'SPGHex', 'BRSLPHex', 'BRSPHex', 
                      'CASLPHex', 'CASPHex', 'SISLPHex', 'SISPHex', 'STSLPHex', 
                      'STSPHex']
                    }
ex_mass = {'C': 12, '13C': 13.0033548378, 'H': 1.007825, 'D': 2.01410178, 
           'N': 14.003074, 'O': 15.994915, 'P': 30.973762, 'S': 31.972071,
           'e': 0.00054858, 'H+': 1.00727642, 'H2O': 18.010565, 'CO2': 43.98983}
ref_oad_ratio = {'OAD01': 0.10, 'OAD02': 0.25, 'OAD03': 0.50, 'OAD04': 0.01, 
                 'OAD05': 0.01, 'OAD06': 0.05, 'OAD07': 0.10, 'OAD08': 0.06, 
                 'OAD09': 0.06, 'OAD10': 0.10, 'OAD11': 0.02, 'OAD12': 0.02, 
                 'OAD13': 0.04, 'OAD14': 0.05, 'OAD15': 0.20, 'OAD16': 0.40, 
                 'OAD17': 0.03, 'OAD18': 0.10, 'OAD19': 0.20, 'OAD20': 0.01}
#endregion

def generate_ref_oad_nl_and_type(each_comb, ontology, deuterium):
    ref_oad_d = {}
    c_ms, c13_ms, o_ms = ex_mass['C'], ex_mass['13C'], ex_mass['O']
    h_ms, d_ms = ex_mass['H'], ex_mass['D']
    h2o_ms = ex_mass['H2O']
    dHs = deuterium*(d_ms - h_ms)
    c13 = c13_ms - c_ms
    mass_shift_counter = 0
    for pos in each_comb:
        ms_sht = 2*ex_mass['H']*mass_shift_counter
        #region Pre
        """ -C*(X-1)-2*H*(X-1)+H+O,             ex) -95.12246   in n-9  """
        pre_mnO_AH = c_ms*(pos-1)+2*h_ms*(pos-1)-h_ms-o_ms-ms_sht+dHs
        """ -C*(X-1)-2*H*(X-1)+O,               ex) -96.130285  in n-9  """
        pre_mnO = c_ms*(pos-1)+2*h_ms*(pos-1)-o_ms-ms_sht+dHs
        """ -C*(X-1)-2*H*(X-1)-H+O,             ex) -97.13811   in n-9  """
        pre_mnO_LH = c_ms*(pos-1)+2*h_ms*(pos-1)+h_ms-o_ms-ms_sht+dHs
        """ -C*(X-1)-2*H*(X-1)-2H+O,            ex) -98.145935  in n-9  """
        pre_mnO_L2H = c_ms*(pos-1)+2*h_ms*(pos-1)+2*h_ms-o_ms-ms_sht+dHs
        """ -C*(X-1)-2*H*(X-1)+H+O-H2O,         ex) -113.13303  in n-9  """
        pre_mnO_AH_H2O = pre_mnO_AH + h2o_ms
        """ -C*(X-1)-2*H*(X-1)+O-H2O,           ex) -114.14085  in n-9  """
        pre_mnO_H2O = pre_mnO + h2o_ms
        """ -C*(X-1)-2*H*(X-1)-H+O-H2O,         ex) -115.14868  in n-9  """
        pre_mnO_LH_H2O = pre_mnO_LH + h2o_ms
        #endregion
        #region Exa
        """ -C*(X)-2*H*(X)+O,                   ex) -110.145935 in n-9  """
        exa_mnO = c_ms*(pos)+2*h_ms*(pos)-o_ms-ms_sht+dHs
        """ -C*(X-1)-2*H*(X-1),                 ex) -124.12520  in n-9  """
        exa = c_ms*(pos)+2*h_ms*(pos-1)-ms_sht+dHs
        """ -C*(X)-2*H*(X)-H,                   ex) -125.13303  in n-9  """
        exa_LH = c_ms*(pos)+2*h_ms*(pos)-h_ms-ms_sht+dHs
        """ -C*(X)-2*H*(X),                     ex) -126.14085  in n-9  """
        exa_L2H = c_ms*(pos)+2*h_ms*(pos)-ms_sht+dHs
        #endregion
        #region Post
        """ -C*(X+1)-2*H*(X)-H+O,               ex) -123.15376  in n-9  """
        post_mnO_LH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-o_ms-ms_sht+dHs
        """ -C*(X+1)-2*H*(X)-2H+O,              ex) -124.16159  in n-9  """
        post_mnO_L2H = c_ms*(pos+1)+2*h_ms*(pos+1)-o_ms-ms_sht+dHs
        """ -C*(X+1)-2*H*(X)+H,                 ex) -137.133025 in n-9  """
        post_AH = c_ms*(pos+1)+2*h_ms*(pos)-h_ms-ms_sht+dHs
        """ -C*(X+1)-2*H*(X),                   ex) -138.14085  in n-9  """
        post = c_ms*(pos+1)+2*h_ms*(pos)-ms_sht+dHs
        """ -C*(X+1)-2*H*(X)-H,                 ex) -139.148675 in n-9  """
        post_LH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-ms_sht+dHs
        """ -C*(X+1)-2*H*(X)-2H,                ex) -140.1565   in n-9  """
        post_L2H = c_ms*(pos+1)+2*h_ms*(pos)+2*h_ms-ms_sht+dHs
        #endregion

        ref_oad_d[f'n-{pos}/dis@n-{pos-1}/+O+H/OAD01'] = pre_mnO_AH
        ref_oad_d[f'n-{pos}/dis@n-{pos-1}/+O/OAD02'] = pre_mnO
        ref_oad_d[f'n-{pos}/dis@n-{pos-1}/+O-H/OAD03'] = pre_mnO_LH
        ref_oad_d[f'n-{pos}/dis@n-{pos-1}/+O-2H/OAD04'] = pre_mnO_L2H
        ref_oad_d[f'n-{pos}/dis@n-{pos-1}/+O+H-H2O/OAD05'] = pre_mnO_AH_H2O
        ref_oad_d[f'n-{pos}/dis@n-{pos-1}/+O-H2O/OAD06'] = pre_mnO_H2O
        ref_oad_d[f'n-{pos}/dis@n-{pos-1}/+O-H-H2O/OAD07'] = pre_mnO_LH_H2O

        ref_oad_d[f'n-{pos}/dis@n-{pos}/+O/OAD08'] = exa_mnO
        ref_oad_d[f'n-{pos}/dis@n-{pos}/none/OAD09'] = exa
        ref_oad_d[f'n-{pos}/dis@n-{pos}/-H/OAD10'] = exa_LH
        ref_oad_d[f'n-{pos}/dis@n-{pos}/-2H/OAD11'] = exa_L2H
        ref_oad_d[f'n-{pos}/dis@n-{pos+1}/+O-H/OAD12'] = post_mnO_LH
        ref_oad_d[f'n-{pos}/dis@n-{pos+1}/+O-2H/OAD13'] = post_mnO_L2H
        ref_oad_d[f'n-{pos}/dis@n-{pos+1}/+H/OAD14'] = post_AH
        ref_oad_d[f'n-{pos}/dis@n-{pos+1}/none/OAD15'] = post
        ref_oad_d[f'n-{pos}/dis@n-{pos+1}/-H/OAD16'] = post_LH
        ref_oad_d[f'n-{pos}/dis@n-{pos+1}/-2H/OAD17'] = post_L2H
        if ontology in lipidclass_dict['Sphingolipids']:
            ref_oad_d[f'n-{pos}/dis@n-{pos+1}/-H2O/OAD18'] \
            = post + h2o_ms
            ref_oad_d[f'n-{pos}/dis@n-{pos+1}/-H-H2O/OAD19'] \
            = post_LH + h2o_ms
            ref_oad_d[f'n-{pos}/dis@n-{pos+1}/-2H-H2O/OAD20'] \
            = post_L2H + h2o_ms
        mass_shift_counter += 1
    return ref_oad_d


def query_essential_diagnostic_ions(df, ref_oad_dict, db_in_sphingobase,
    c_num, db_num, for_range, tolerance, must_nl_cut_off_dict, structure_dict):
    tol = tolerance
    frag = 'frag m/z'
    diagnostic_nl_dict = {}
    ontology = structure_dict['Ontology']
    ref_mz = structure_dict['Ref precursor Mz']
    # diagnostic_1_type = must_nl_cut_off_dict['diagnostic_1'][0]
    diagnostic_1_cutoff = must_nl_cut_off_dict['diagnostic_1'][1]
    # diagnostic_2_type = must_nl_cut_off_dict['diagnostic_2'][0]
    diagnostic_2_cutoff = must_nl_cut_off_dict['diagnostic_2'][1]
    cut_1_df = df[df['Ratio(%)'] >= diagnostic_1_cutoff]
    cut_2_df = df[df['Ratio(%)'] >= diagnostic_2_cutoff]
    sph_df = df[df['Ratio(%)'] >= must_nl_cut_off_dict['sphingobase']]
    if ontology not in lipidclass_dict['Sphingolipids']:
        for pos, tag_and_nl in ref_oad_dict.items():
            each_pos_bool = []
            tags_of_dgn_1 = [f'n-{each}/dis@n-{each-1}/+O-H/OAD03' for each in pos]
            tags_of_dgn_2 = [f'n-{each}/dis@n-{each+1}/-H/OAD16' for each in pos]
            dgn_1_nls = [tag_and_nl[tag] for tag in tags_of_dgn_1]
            dgn_2_nls = [tag_and_nl[tag] for tag in tags_of_dgn_2]
            head_dgn_1_mzs = [ref_mz - nl - tol for nl in dgn_1_nls]
            tail_dgn_1_mzs = [ref_mz - nl + tol for nl in dgn_1_nls]
            head_dgn_2_mzs = [ref_mz - nl - tol for nl in dgn_2_nls]
            tail_dgn_2_mzs = [ref_mz - nl + tol for nl in dgn_2_nls]
            for i in for_range:
                head_1_mz = head_dgn_1_mzs[i]
                tail_1_mz = tail_dgn_1_mzs[i]
                head_2_mz = head_dgn_2_mzs[i]
                tail_2_mz = tail_dgn_2_mzs[i]
                dgn_1_df = cut_1_df[(cut_1_df[frag] > head_1_mz)
                                  & (cut_1_df[frag] < tail_1_mz)]
                dgn_2_df = cut_2_df[(cut_2_df[frag] > head_2_mz)
                                  & (cut_2_df[frag] < tail_2_mz)]
                len_of_dgn_1_df = len(dgn_1_df)
                len_of_dgn_2_df = len(dgn_2_df)
                if len_of_dgn_1_df > 0 and len_of_dgn_2_df > 0 :
                    each_pos_bool.append(True)
                else:
                    each_pos_bool.append(False)
            diagnostic_nl_dict[pos] = each_pos_bool
    else:
        if db_in_sphingobase:
            last_loop = db_num -1
            for pos, tag_and_nl in ref_oad_dict.items():
                each_pos_bool = []
                tags_of_dgn_1 = [
                    f'n-{each}/dis@n-{each-1}/+O-H/OAD03' for each in pos]
                tags_of_dgn_2 = [
                    f'n-{each}/dis@n-{each+1}/-H/OAD16' for each in pos]
                tags_of_dgn_3 = [
                    f'n-{each}/dis@n-{each-1}/+O-H-H2O/OAD07' for each in pos]
                tags_of_dgn_4 = [
                    f'n-{each}/dis@n-{each+1}/-H-H2O/OAD19' for each in pos]
                dgn_1_nls = [tag_and_nl[tag] for tag in tags_of_dgn_1]
                dgn_2_nls = [tag_and_nl[tag] for tag in tags_of_dgn_2]
                dgn_3_nls = [tag_and_nl[tag] for tag in tags_of_dgn_3]
                dgn_4_nls = [tag_and_nl[tag] for tag in tags_of_dgn_4]
                head_dgn_1_mzs = [ref_mz - nl - tol for nl in dgn_1_nls]
                tail_dgn_1_mzs = [ref_mz - nl + tol for nl in dgn_1_nls]
                head_dgn_2_mzs = [ref_mz - nl - tol for nl in dgn_2_nls]
                tail_dgn_2_mzs = [ref_mz - nl + tol for nl in dgn_2_nls]
                head_dgn_3_mzs = [ref_mz - nl - tol for nl in dgn_3_nls]
                tail_dgn_3_mzs = [ref_mz - nl + tol for nl in dgn_3_nls]
                head_dgn_4_mzs = [ref_mz - nl - tol for nl in dgn_4_nls]
                tail_dgn_4_mzs = [ref_mz - nl + tol for nl in dgn_4_nls]
                db_not_in_end_sph = c_num < 0 or c_num - pos[-1] != 4
                for i in for_range:
                    if i < last_loop or db_not_in_end_sph:
                        head_1_mz = head_dgn_1_mzs[i]
                        tail_1_mz = tail_dgn_1_mzs[i]
                        head_2_mz = head_dgn_2_mzs[i]
                        tail_2_mz = tail_dgn_2_mzs[i]
                        head_3_mz = head_dgn_3_mzs[i]
                        tail_3_mz = tail_dgn_3_mzs[i]
                        head_4_mz = head_dgn_4_mzs[i]
                        tail_4_mz = tail_dgn_4_mzs[i]
                        dgn_1_df = cut_1_df[(cut_1_df[frag] > head_1_mz)
                                           &(cut_1_df[frag] < tail_1_mz)]
                        dgn_2_df = cut_2_df[(cut_2_df[frag] > head_2_mz)
                                           &(cut_2_df[frag] < tail_2_mz)]
                        dgn_3_df = cut_1_df[(cut_1_df[frag] > head_3_mz)
                                           &(cut_1_df[frag] < tail_3_mz)]
                        dgn_4_df = cut_2_df[(cut_2_df[frag] > head_4_mz)
                                           &(cut_2_df[frag] < tail_4_mz)]
                        len_of_dgn_1_df = len(dgn_1_df)
                        len_of_dgn_2_df = len(dgn_2_df)
                        len_of_dgn_3_df = len(dgn_3_df)
                        len_of_dgn_4_df = len(dgn_4_df)
                        if ((len_of_dgn_1_df > 0 and len_of_dgn_2_df > 0)
                            or (len_of_dgn_3_df > 0 and len_of_dgn_4_df > 0)):
                            each_pos_bool.append(True)
                        else:
                            each_pos_bool.append(False)
                    else:
                        #region Tag list
                        p = pos[-1]
                        post_nl = f'n-{p}/dis@n-{p+1}/none/OAD15'
                        post_LH_nl = f'n-{p}/dis@n-{p+1}/-H/OAD16'
                        post_L2H_nl = f'n-{p}/dis@n-{p+1}/-2H/OAD17'
                        post_H2O_nl = f'n-{p}/dis@n-{p+1}/-H2O/OAD18'
                        post_LH_H2O_nl = f'n-{p}/dis@n-{p+1}/-H-H2O/OAD19'
                        post_L2H_H2O_nl = f'n-{p}/dis@n-{p+1}/-2H-H2O/OAD20'
                        #endregion
                        #region frag m/z list
                        h_post_mz = ref_mz - tag_and_nl[post_nl] - tol
                        t_post_mz = ref_mz - tag_and_nl[post_nl] + tol
                        h_post_LH_mz = ref_mz - tag_and_nl[post_LH_nl] - tol
                        t_post_LH_mz = ref_mz - tag_and_nl[post_LH_nl] + tol
                        h_post_L2H_mz = ref_mz - tag_and_nl[post_L2H_nl] - tol
                        t_post_L2H_mz = ref_mz - tag_and_nl[post_L2H_nl] + tol
                        h_post_H2O_mz = ref_mz - tag_and_nl[post_H2O_nl] - tol
                        t_post_H2O_mz = ref_mz - tag_and_nl[post_H2O_nl] + tol
                        h_post_LH_H2O_mz = ref_mz - tag_and_nl[post_LH_H2O_nl] - tol
                        t_post_LH_H2O_mz = ref_mz - tag_and_nl[post_LH_H2O_nl] + tol
                        h_post_L2H_H2O_mz = ref_mz - tag_and_nl[post_L2H_H2O_nl] - tol
                        t_post_L2H_H2O_mz = ref_mz - tag_and_nl[post_L2H_H2O_nl] + tol
                        #endregion
                        post_bool = any((sph_df[frag] > h_post_mz)
                                       &(sph_df[frag] < t_post_mz))
                        post_LH_bool = any((sph_df[frag] > h_post_LH_mz)
                                          &(sph_df[frag] < t_post_LH_mz))
                        post_L2H_bool = any((sph_df[frag] > h_post_L2H_mz)
                                           &(sph_df[frag] < t_post_L2H_mz))
                        post_H2O_bool = any((sph_df[frag] > h_post_H2O_mz)
                                           &(sph_df[frag] < t_post_H2O_mz))
                        post_LH_H2O_bool = any((sph_df[frag] > h_post_LH_H2O_mz)
                                              &(sph_df[frag] < t_post_LH_H2O_mz))
                        post_L2H_H2O_bool = any((sph_df[frag] > h_post_L2H_H2O_mz)
                                               &(sph_df[frag] < t_post_L2H_H2O_mz))
                        posts = [post_bool, post_LH_bool, 
                                 post_L2H_bool].count(True)
                        post_h2os = [post_H2O_bool, post_LH_H2O_bool, 
                                     post_L2H_H2O_bool].count(True)
                        Total_bool = posts >= 2 or post_h2os >= 2
                        if Total_bool:
                            each_pos_bool.append(True)
                        else:
                            each_pos_bool.append(False)
                diagnostic_nl_dict[pos] = each_pos_bool
        else:
            for pos, tag_and_nl in ref_oad_dict.items():
                each_pos_bool = []
                tags_of_dgn_1 = [
                    f'n-{each}/dis@n-{each-1}/+O-H/OAD03' for each in pos]
                tags_of_dgn_2 = [
                    f'n-{each}/dis@n-{each+1}/-H/OAD16' for each in pos]
                tags_of_dgn_3 = [
                    f'n-{each}/dis@n-{each-1}/+O-H-H2O/OAD07' for each in pos]
                tags_of_dgn_4 = [
                    f'n-{each}/dis@n-{each+1}/-H-H2O/OAD19' for each in pos]
                dgn_1_nls = [tag_and_nl[tag] for tag in tags_of_dgn_1]
                dgn_2_nls = [tag_and_nl[tag] for tag in tags_of_dgn_2]
                dgn_3_nls = [tag_and_nl[tag] for tag in tags_of_dgn_3]
                dgn_4_nls = [tag_and_nl[tag] for tag in tags_of_dgn_4]
                head_dgn_1_mzs = [ref_mz - nl - tol for nl in dgn_1_nls]
                tail_dgn_1_mzs = [ref_mz - nl + tol for nl in dgn_1_nls]
                head_dgn_2_mzs = [ref_mz - nl - tol for nl in dgn_2_nls]
                tail_dgn_2_mzs = [ref_mz - nl + tol for nl in dgn_2_nls]
                head_dgn_3_mzs = [ref_mz - nl - tol for nl in dgn_3_nls]
                tail_dgn_3_mzs = [ref_mz - nl + tol for nl in dgn_3_nls]
                head_dgn_4_mzs = [ref_mz - nl - tol for nl in dgn_4_nls]
                tail_dgn_4_mzs = [ref_mz - nl + tol for nl in dgn_4_nls]
                for i in for_range:
                    head_1_mz = head_dgn_1_mzs[i]
                    tail_1_mz = tail_dgn_1_mzs[i]
                    head_2_mz = head_dgn_2_mzs[i]
                    tail_2_mz = tail_dgn_2_mzs[i]
                    head_3_mz = head_dgn_3_mzs[i]
                    tail_3_mz = tail_dgn_3_mzs[i]
                    head_4_mz = head_dgn_4_mzs[i]
                    tail_4_mz = tail_dgn_4_mzs[i]
                    dgn_1_df = cut_1_df[(cut_1_df[frag] > head_1_mz)
                                       &(cut_1_df[frag] < tail_1_mz)]
                    dgn_2_df = cut_2_df[(cut_2_df[frag] > head_2_mz)
                                       &(cut_2_df[frag] < tail_2_mz)]
                    dgn_3_df = cut_1_df[(cut_1_df[frag] > head_3_mz)
                                       &(cut_1_df[frag] < tail_3_mz)]
                    dgn_4_df = cut_2_df[(cut_2_df[frag] > head_4_mz)
                                       &(cut_2_df[frag] < tail_4_mz)]
                    len_of_dgn_1_df = len(dgn_1_df)
                    len_of_dgn_2_df = len(dgn_2_df)
                    len_of_dgn_3_df = len(dgn_3_df)
                    len_of_dgn_4_df = len(dgn_4_df)
                    if ((len_of_dgn_1_df > 0 and len_of_dgn_2_df > 0)
                        or (len_of_dgn_3_df > 0 and len_of_dgn_4_df > 0)):
                        each_pos_bool.append(True)
                    else:
                        each_pos_bool.append(False)
                diagnostic_nl_dict[pos] = each_pos_bool
    return diagnostic_nl_dict