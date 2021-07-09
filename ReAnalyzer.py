import math


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
    ref_oad_d, tag_and_nl = {}, {}
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

        tag_and_nl[f'n-{pos}/dis@n-{pos-1}/+O+H/OAD01'] = pre_mnO_AH
        tag_and_nl[f'n-{pos}/dis@n-{pos-1}/+O/OAD02'] = pre_mnO
        tag_and_nl[f'n-{pos}/dis@n-{pos-1}/+O-H/OAD03'] = pre_mnO_LH
        tag_and_nl[f'n-{pos}/dis@n-{pos-1}/+O-2H/OAD04'] = pre_mnO_L2H
        tag_and_nl[f'n-{pos}/dis@n-{pos-1}/+O+H-H2O/OAD05'] = pre_mnO_AH_H2O
        tag_and_nl[f'n-{pos}/dis@n-{pos-1}/+O-H2O/OAD06'] = pre_mnO_H2O
        tag_and_nl[f'n-{pos}/dis@n-{pos-1}/+O-H-H2O/OAD07'] = pre_mnO_LH_H2O

        tag_and_nl[f'n-{pos}/dis@n-{pos}/+O/OAD08'] = exa_mnO
        tag_and_nl[f'n-{pos}/dis@n-{pos}/none/OAD09'] = exa
        tag_and_nl[f'n-{pos}/dis@n-{pos}/-H/OAD10'] = exa_LH
        tag_and_nl[f'n-{pos}/dis@n-{pos}/-2H/OAD11'] = exa_L2H
        tag_and_nl[f'n-{pos}/dis@n-{pos+1}/+O-H/OAD12'] = post_mnO_LH
        tag_and_nl[f'n-{pos}/dis@n-{pos+1}/+O-2H/OAD13'] = post_mnO_L2H
        tag_and_nl[f'n-{pos}/dis@n-{pos+1}/+H/OAD14'] = post_AH
        tag_and_nl[f'n-{pos}/dis@n-{pos+1}/none/OAD15'] = post
        tag_and_nl[f'n-{pos}/dis@n-{pos+1}/-H/OAD16'] = post_LH
        tag_and_nl[f'n-{pos}/dis@n-{pos+1}/-2H/OAD17'] = post_L2H
        if ontology in lipidclass_dict['Sphingolipids']:
            tag_and_nl[f'n-{pos}/dis@n-{pos+1}/-H2O/OAD18'] \
            = post + h2o_ms
            tag_and_nl[f'n-{pos}/dis@n-{pos+1}/-H-H2O/OAD19'] \
            = post_LH + h2o_ms
            tag_and_nl[f'n-{pos}/dis@n-{pos+1}/-2H-H2O/OAD20'] \
            = post_L2H + h2o_ms
        mass_shift_counter += 1
    ref_oad_d[each_comb] = tag_and_nl
    return ref_oad_d

def query_essential_diagnostic_ions(df, ref_oad_dict, db_in_SPB,
    c_num, db_num, tolerance, must_nl_cut_off_dict, structure_dict):
    tol = tolerance
    for_range = range(0, db_num)
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
            tags_of_dgn_1 = [f'n-{each}/dis@n-{each-1}/+O-H/OAD03' 
                             for each in pos]
            tags_of_dgn_2 = [f'n-{each}/dis@n-{each+1}/-H/OAD16' 
                             for each in pos]
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
        if db_in_SPB:
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

def calc_presence_ratios_and_score(ref_oad_dict, cut_df, ref_precursor_mz, 
    ms_tol_ppm, sph_set):
    counter = 0
    score_dict = {}
    tol = math_floor(ms_tol_ppm*ref_precursor_mz/(1000*1000), 6)
    def get_n_discription(txt):
        if ',)' in txt:
            edit = 'n-' + txt.replace('(', '').replace(',)', '')
        else:
            edit = 'n-' + txt.replace('(', '').replace(')', '').replace(' ', '')
        return edit
    for positions, tag_and_nl in ref_oad_dict.items():
        presence_counter, ratio_sum = 0, 0
        next_to_3oh = sph_set[0] and ((sph_set[1]-positions[-1]) == 4)
        each_score_d = {'Positions': '', 'N-description': '',
                        'Score': 0, 'Ratio sum': 0, 'Presence': 0,
                        'Notice': 'ReAnalysis', 'Measured peaks': [],
                        'Ref peaks': [], 'Peaks dict': {}}
        # Measured peaks: [[Measured m/z, Measured ratio, ppm], [...]]
        # Ref peaks: [[OAD type, Ref m/z, Ref NL, Ref ratio], [...]]
        # Peaks dict = {'n-9/dis@n-8/+O/OAD03': [Ref m/z, Ref delta, 
        #                                        Measured m/z, 
        #                                        Measured ratio, ppm]}
        peaks_dict = {}
        ref_peaks, measured_peaks = [], []
        for tag, ref_nl in tag_and_nl.items():
            db = int(tag.split('/')[0].replace('n-', ''))
            ref_mz = math_floor((ref_precursor_mz - ref_nl), 4)
            if next_to_3oh and db == positions[-1]:
                tag_num = int(tag.split('OAD')[-1])
                if tag_num >= 15:
                    mz, ratio, ppm = query_matched_ion_by_ppm(cut_df, ref_mz, tol)
                else: continue
            else:
                mz, ratio, ppm = query_matched_ion_by_ppm(cut_df, ref_mz, tol)
            if mz > 0:
                presence_counter += 1
                ratio_sum += ratio
            peaks_dict[tag] = [ref_mz, ref_nl, mz, ratio, ppm]
        #region Rel int ratio
        # base_peaks = []
        # for key, v in peaks_dict.items():
        #     oad_type = key.split('/')[-1]
        #     if oad_type == 'OAD03': base_peaks.append(v[3])
        # if next_to_3oh and len(positions) == 1:
        #     for key, v in peaks_dict.items():
        #         #append -> [OAD type, Ref m/z, Ref NL, Ref ratio]
        #         oad_type = key.split('/')[-1]
        #         ref_peaks.append([key, v[0], v[1], ref_oad_ratio[oad_type]])
        #         if v[2] > 0: #[Measured m/z, Measured ratio, ppm]
        #             measured_peaks.append([v[2], v[3], v[4]])
        # else:
        #     base_ratio = statistics.mean(base_peaks)
        #     for key, v in peaks_dict.items():
        #         #append -> [OAD type, Ref m/z, Ref NL, Ref ratio]
        #         oad_type = key.split('/')[-1]
        #         ref_peaks.append([key, v[0], v[1], 
        #                         get_rel_ratio(base_ratio, oad_type)])
        #         if v[2] > 0: #[Measured m/z, Measured ratio, ppm]
        #             measured_peaks.append([v[2], v[3], v[4]])
        #endregion
        for key, v in peaks_dict.items():
            #append -> [OAD type, Ref m/z, Ref NL, Ref ratio]
            ref_ratio = get_ref_ratio_via_db_position(positions, key)
            ref_peaks.append([key, v[0], v[1], ref_ratio])
            if v[2] > 0: #[Measured m/z, Measured ratio, ppm]
                measured_peaks.append([v[2], v[3], v[4]])
        each_score_d['Positions'] = positions
        each_score_d['N-description'] = get_n_discription(str(positions))
        each_score_d['Measured peaks'] = measured_peaks
        each_score_d['Ref peaks'] = ref_peaks
        each_score_d['Peaks dict'] = peaks_dict
        presence = math_floor(presence_counter/len(tag_and_nl)*100, 4)
        ratio_sum = math_floor(ratio_sum, 4)
        acts, refs = construct_ratio_vec(peaks_dict, ref_peaks)
        # acts, refs = constrcut_ranking_vec(peaks_dict, ref_peaks)
        each_score_d['Score'] = get_msms_similarity_score(acts, refs)
        each_score_d['Ratio sum'] = ratio_sum
        each_score_d['Presence'] = presence
        score_dict[counter] = each_score_d
        counter += 1
    return score_dict

def get_ref_ratio_via_db_position(positions, key):
    pre = ['OAD01', 'OAD02', 'OAD03', 'OAD04', 'OAD05', 'OAD06', 'OAD07']
    post = ['OAD15', 'OAD16', 'OAD17']
    db_num, oad_type = len(positions), key.split('/')[-1]
    db_delta = [positions[i+1]-positions[i] for i in range(db_num-1)]
    if 2 in db_delta: # delta == 2 -> Conjugated C=C, ex) n-7,9
        conj_pos = db_delta.index(2)
        conj_start = f'n-{positions[conj_pos]}'
        conj_end = f'n-{positions[conj_pos+1]}'
        if ((key.startswith(conj_start) and oad_type in post) or
            (key.startswith(conj_end) and oad_type in pre)):
            ref_ratio = 0.015
            return ref_ratio
    first_db, last_db = f'n-{positions[0]}', f'n-{positions[-1]}'
    half_pre = (db_num >= 3 and key.startswith(first_db) and oad_type in pre)
    half_post = (db_num >= 3 and key.startswith(last_db) and oad_type in post)
    ref_ratio = ref_oad_ratio[oad_type]
    if half_pre or half_post:
        ref_ratio = math_floor(ref_ratio/2, 3)
    return ref_ratio

def query_matched_ion_by_ppm(cut_df, ref_mz, tol):
    calc_ppm = lambda mz, ref: math_floor(abs((mz-ref)/ref*1000*1000), 2)
    start, end = ref_mz-tol, ref_mz+tol
    ex_df = cut_df[(cut_df['frag m/z'] >= start) & (cut_df['frag m/z'] <= end)]
    if len(ex_df) == 1:
        mz, ratio = ex_df['frag m/z'].values[0], ex_df['Ratio(%)'].values[0]
        ppm = calc_ppm(mz, ref_mz)
    elif len(ex_df) >= 2:
        mzs, ratios = ex_df['frag m/z'].values, ex_df['Ratio(%)'].values
        ppms = [calc_ppm(mz, ref_mz) for mz in mzs]
        li = [[mz, ratio, ppm] for mz, ratio, ppm in zip(mzs, ratios, ppms)]
        s_li = sorted(li, key=lambda x: x[2])
        mz, ratio, ppm = s_li[0][0], s_li[0][1], s_li[0][2]
    else:
        mz, ratio, ppm = 0, 0, 0
    return mz, ratio, ppm

def construct_ratio_vec(peaks_dict, ref_peaks):
    # Peaks dict = {'n-9/dis@n-8/+O/OAD03': 
    #               [Ref m/z, Ref delta, Measured m/z, Measured ratio, ppm]}
    # ref_peaks = [[OAD type, Ref m/z, Ref NL, Ref ratio], ...]
    sort_peaks = dict(sorted(peaks_dict.items(), key=lambda x: x[1][0]))
    ratio_dict = {v[0]: {'Measured': 0, 'Ref': 0} for v in sort_peaks.values()}
    for key, v in sort_peaks.items():
        # oad, ref_mz = key.split('/')[-1], v[0]
        ref_mz = v[0]
        ratio_dict[ref_mz]['Measured'] = math_floor(v[3], 4)
        ref_ratio = [li[3] for li in ref_peaks if li[0] == key][0]
        stacked_ratio = math_floor(ref_ratio+ratio_dict[ref_mz]['Ref'], 2)
        ratio_dict[ref_mz]['Ref'] = stacked_ratio
    acts = [v['Measured'] for v in ratio_dict.values()]
    refs = [v['Ref'] for v in ratio_dict.values()]
    return acts, refs

def get_msms_similarity_score(acts, refs):
    # acts = [v[3] for v in peaks_dict.values()]
    # nozero_acts = [v[3] for v in peaks_dict.values() if v[3] > 0]
    # refs = [ref_oad_ratio[tag.split('/')[-1]] for tag in peaks_dict.keys()]
    nozero_acts = [v for v in acts if v > 0]
    try:
        act_min_digit = len(str(min(nozero_acts)).split('.')[1])
    except:
        act_min_digit = 0
    ref_min_digit = len(str(min(refs)).split('.')[1])
    digit = act_min_digit if act_min_digit >= ref_min_digit else ref_min_digit
    acts = [int(v*10**digit) for v in acts]
    refs = [int(v*10**digit) for v in refs]
    if sum(acts) > 0 and sum(refs):
        # rev_dotp = calc_dot_product(acts, refs)
        # score = math_floor(rev_dotp, 2)
        sim = calc_similarity_score(acts, refs)
        score = math_floor(sim, 2)
        # score = math_floor(rev_dotp+sim, 2)
    else: score = 0
    return score

def calc_dot_product(acts, refs):
    dotp =  lambda acts, refs: sum(act*ref for act, ref in zip(acts, refs))**2
    scalar = lambda vec: sum(v**2 for v in vec)
    dotp_sml = lambda acts, refs: dotp(acts, refs)/(scalar(acts)*scalar(refs))
    score = dotp_sml(acts, refs)*1000
    return score

def calc_similarity_score(acts, refs):
    sqrt_dotp =  lambda acts, refs: sum(math.sqrt(act*ref) for act, ref in zip(acts, refs))
    scalar = lambda vec: sum(vec)
    sim = lambda acts, refs: sqrt_dotp(acts, refs)/math.sqrt(scalar(acts)*scalar(refs))
    score = sim(acts, refs)*1000
    return score

def set_oad_graph_dict_value(oad_dict, lipid_info, auto_magnification):
    #region Data structure
    # dict = {0:              {'Positions': '', 'N-description': '', 'Score': float, 
    #                          'Ratio sum': float, 'Presence': float, 'Notice': '',
    #                          'Measured peaks': [[Measured m/z, Measured ratio, ppm], [...]],
    #                          'Ref peaks': [[OAD type, Ref m/z, Ref NL, Ref ratio], [...]],
    #                          'Peaks dict': {'n-9/dis@n-8/+O/': [Ref m/z, Ref delta, Measured m/z, Measured ratio, ppm]}
    #                         },
    #         1:              {....},
    #         Determined db:  {'Positions': '', 'N-description': '', 'Score': float,
    #                          'Ratio sum': float, 'Presence': float, 'Notice': '',
    #                          'Measured peaks': [[Measured m/z, Measured ratio, ppm], [...]],
    #                          'Ref peaks': [[OAD type, Ref m/z, Ref NL, Ref ratio], [...]]
    #                         }
    #         }
    #endregion
    graph_dict = {'MS2 Mz': 0, 'Ref precursor Mz': 0, 'Ontology': '', 
                  'x-range': [], 'Magnification': 1, 'Bar_width': 1}
    d_len = len(oad_dict)
    xtick_num = 10
    if lipid_info['MS2 Mz'] > 0: precursor_mz = lipid_info['MS2 Mz']
    else: precursor_mz = lipid_info['Precursor Mz']
    ref_precursor_mz = lipid_info['Ref precursor Mz']
    ontology = lipid_info['Ontology']
    exp_peaks = 'Measured peaks'
    ref_peaks = 'Ref peaks'
    #region x-ragne
    measured_oad_ions_3, measured_oad_ions_2, measured_oad_ions_1 = [], [], []
    min_mzs = []
    if d_len == 4:
        for rank, d in oad_dict['Moiety-1'].items():
            if d['N-description'] == 'Unresolved':
                continue
            measured_oad_ions_1 = [li[1] for li in d[exp_peaks]]
            ref_mz_1 = [li[1] for li in d[ref_peaks]]
            min_mzs.append(min(ref_mz_1))
    if d_len == 5:
        for rank, d in oad_dict['Moiety-2'].items():
            if d['N-description'] == 'Unresolved':
                continue
            measured_oad_ions_2 = [li[1] for li in d[exp_peaks]]
            ref_mz_2 = [li[1] for li in d[ref_peaks]]
            min_mzs.append(min(ref_mz_2))
        for rank, d in oad_dict['Moiety-1'].items():
            if d['N-description'] == 'Unresolved':
                continue
            measured_oad_ions_1 = [li[1] for li in d[exp_peaks]]
            ref_mz_1 = [li[1] for li in d[ref_peaks]]
            min_mzs.append(min(ref_mz_1))
    if d_len == 6:
        for rank, d in oad_dict['Moiety-3'].items():
            if d['N-description'] == 'Unresolved':
                continue
            measured_oad_ions_3 = [li[1] for li in d[exp_peaks]]
            ref_mz_3 = [li[1] for li in d[ref_peaks]]
            min_mzs.append(min(ref_mz_3))
        for rank, d in oad_dict['Moiety-2'].items():
            if d['N-description'] == 'Unresolved':
                continue
            ref_mz_2 = [li[1] for li in d[ref_peaks]]
            measured_oad_ions_2 = [li[1] for li in d[exp_peaks]]
            min_mzs.append(min(ref_mz_2))
        for rank, d in oad_dict['Moiety-1'].items():
            if d['N-description'] == 'Unresolved':
                continue
            measured_oad_ions_1 = [li[1] for li in d[exp_peaks]]
            ref_mz_1 = [li[1] for li in d[ref_peaks]]
            min_mzs.append(min(ref_mz_1))
    if min_mzs:
        x_min = min(min_mzs) -50
    else:
        x_min = precursor_mz - 250
    x_max = math.ceil(precursor_mz/xtick_num)*xtick_num+5
    #endregion
    bar_width = math.floor(10*(x_max-x_min)/400)/10
    #region Magnification
    max_ratio_list = []
    max_criteria = 30
    if measured_oad_ions_3: max_ratio_list.append(max(measured_oad_ions_3))
    if measured_oad_ions_2: max_ratio_list.append(max(measured_oad_ions_2))
    if measured_oad_ions_1: max_ratio_list.append(max(measured_oad_ions_1))
    if auto_magnification == 10:
        if max_ratio_list:
            final_magnification = math.floor(max_criteria/max(max_ratio_list))
            if final_magnification > 500:
                final_magnification = 500
        else:
            final_magnification = 30
    else:
        final_magnification = auto_magnification
    #endregion
    graph_dict['MS2 Mz'] = precursor_mz
    graph_dict['Ref precursor Mz'] = ref_precursor_mz
    graph_dict['Ontology'], graph_dict['x-range'] = ontology, [x_min, x_max]
    graph_dict['Magnification'] = final_magnification
    graph_dict['Bar_width'] = bar_width
    return graph_dict


#region Utilities
def math_floor(num, digit):
    floored = math.floor(num*10**digit)/(10**digit)
    return floored

#endregion









