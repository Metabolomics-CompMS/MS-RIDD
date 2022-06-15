import itertools
import math
import pandas as pd
import pickle
import re
import pprint
import time


#region lipidomics infomartions
no_acyl_list = ['CoQ', 'Vitamine E', 'CASulfate', 'CA', 'Vitamine D', 
                'SSulfate', 'Cholesterol', 'SHex', 'SPE', 'BAHex', 
                'BASulfate', 'SPEHex', 'SPGHex']
mono_acyl_list = ['FA', 'NAE', 'CAR', 'MG', 'LDGCC', 'LDGTS/A', 
                  'LPA', 'EtherLPC', 'LPC', 'EtherLPE', 'LPE', 'EtherLPG', 'LPG', 'LPI', 'LPS', 
                  'VAE', 'PhytoSph', 'DHSph', 'Sph', 'DCAE', 'GDCAE', 'GLCAE', 'TDCAE', 'TLCAE', 
                  'AHexCAS', 'AHexCS', 'AHexSIS', 'AHexBRS', 'AHexSTS', 
                  'BRSE', 'CASE', 'CE', 'SISE', 'STSE']
tri_acyl_list = ['ADGGA', 'EtherTG', 'TG', 'HBMP', 'MLCL', 
                 'Cer_EBDS', 'Cer_EODS', 'Cer_EOS', 'AHexCer', 'HexCer_EOS', 'ASM']
triacyls_sphingolipids = ['AHexCer', 'ASM', 'Cer_EBDS', 'Cer_EODS', 'Cer_EOS', 'HexCer_EOS']
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
adduct_dict_neg = {'[M-H]-': -1.00727642, '[M+CH3COO]-': 59.01385292}
adduct_dict_pos = {'[M]+': -0.00054858, '[M+H]+': 1.00727642, '[M+NH4]+': 18.03382555, '[M+Na]+': 22.9892207 }
msdial_std_columns = ['ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name', 
                      'Adduct type', 'Post curation result', 'Fill %', 'MS/MS assigned', 
                      'Reference RT', 'Reference m/z', 'Formula', 'Ontology', 
                      'INCHIKEY', 'SMILES', 'Annotation tag (VS1.0)', 'RT matched', 
                      'm/z matched', 'MS/MS matched', 'Comment', 
                      'Manually modified for quantification', 'Isotope tracking parent ID', 
                      'Isotope tracking weight number', 'Total score', 
                      'RT similarity', 'Dot product', 'Reverse dot product', 
                      'Fragment presence %', 'S/N average', 'Spectrum reference file name', 
                      'Spectrum reference file name', 'MS1 isotopic spectrum', 'MS/MS spectrum']
exclude_subclass = ['Others', 'CoQ', 'CL', 'OxTG', 'FAHFATG', 
                    'Vitamin_D', 'Vitamin_E', 'Vitamin D', 'Vitamin E']
ref_oad_ratio = {
    'OAD01': 0.10, 'OAD02': 0.25, 'OAD03': 0.50, 'OAD04': 0.01, 
    'OAD05': 0.01, 'OAD06': 0.05, 'OAD07': 0.10, 'OAD08': 0.06, 
    'OAD09': 0.06, 'OAD10': 0.10, 'OAD11': 0.02, 'OAD12': 0.02, 
    'OAD13': 0.04, 'OAD14': 0.05, 'OAD15': 0.20, 'OAD16': 0.40, 
    'OAD17': 0.03, 'OAD18': 0.10, 'OAD19': 0.20, 'OAD20': 0.01
}

rel_oad_ratio = {
    'OAD01': 0.20, 'OAD02': 0.50, 'OAD03': 1.00, 'OAD04': 0.02, 
    'OAD05': 0.02, 'OAD06': 0.10, 'OAD07': 0.20, 'OAD08': 0.12, 
    'OAD09': 0.20, 'OAD10': 0.40, 'OAD11': 0.04, 'OAD12': 0.04, 
    'OAD13': 0.08, 'OAD14': 0.10, 'OAD15': 0.40, 'OAD16': 0.80, 
    'OAD17': 0.06, 'OAD18': 0.20, 'OAD19': 0.40, 'OAD20': 0.02
}

# rel_oad_ratio = {'OAD01': 0.15, 'OAD02': 0.46, 'OAD03': 1.00, 'OAD04': 0.06, 
#                  'OAD05': 0.05, 'OAD06': 0.15, 'OAD07': 0.30, 'OAD08': 0.11, 
#                  'OAD09': 0.13, 'OAD10': 0.31, 'OAD11': 0.06, 'OAD12': 0.04, 
#                  'OAD13': 0.11, 'OAD14': 0.09, 'OAD15': 0.35, 'OAD16': 0.71, 
#                  'OAD17': 0.07, 'OAD18': 0.18, 'OAD19': 0.32, 'OAD20': 0.04}
#endregion

class SingleAnalyzer(object):
    """ Definition of class that generates data object to analyze, visualize,
        modify, and save OAD-MS/MS spectra and annotations for single data 
        format.
    """
    def __init__(self, tabel_format, directry, prefix, input_data, 
                 ms_tolerance_ppm, must_nl_cut_off_dict, cut_off_ratio, 
                 file_name, sec_rep, sec_bar, each_rep, each_bar, timer):
        """
        Args:
            tabel_format
            directry
            prefix
            input_data
            ms_tolerance_ppm
            must_nl_cut_off_dict (dict): 
                types of essential ions and relative intensity threshold
            cut_off_ratio
            file_name
            sec_rep
            sec_bar
            each_rep
            each_bar
            
        """
        sec_rep.set("Data Pre-processing")
        sec_bar["maximum"] = 9
        #region Data preprocessing
        if tabel_format == 'Alignment':
            raw_table = pd.read_csv(input_data, skiprows=[0,1,2,3], sep='\t')
            raw_table = raw_table.rename(columns={
                'Alignment ID': 'ID', 'Average Rt(min)': 'RT(min)', 
                'Average Mz': 'Precursor m/z'})
        elif tabel_format == 'PeakList':
            raw_table = pd.read_csv(input_data, sep='\t')
            raw_table = raw_table.rename(columns={
                'PeakID': 'ID', 'Title': 'Metabolite name', 'RT (min)': 'RT(min)', 
                'Adduct': 'Adduct type', 'InChIKey': 'INCHIKEY',
                'MSMS spectrum': 'MS/MS spectrum'})
        # elif tabel_format == 'Merged text':
        #     raw_table = pd.read_csv(input_data, sep='\t')
        new_table = get_annotated_df(raw_table)
        columns = new_table.columns.values
        if 'Height' not in columns:
            new_table['Height'] = 0
            samples = [col for col in columns 
                       if col not in msdial_std_columns and 'Blank' not in col]
            for row, df in new_table.iterrows():
                new_table.loc[row:row, ['Height']] = df[samples].mean()
        if 'Data from' not in columns:
            new_table['Data from'] = ''
        self.target_table = new_table
        self.target_table.fillna({'Comment': ''})
        comments = self.target_table['Comment'].values.tolist()
        comments = ['' if isinstance(v, float) else v for v in comments]
        ex_list = ['isotope of', 'adduct linked to', 'Unit',
                   'similar chromatogram', 'found in']
        user_comment_list = [refine_comments(txt, ex_list) for txt in comments]
        self.target_table['User comment'] = user_comment_list
        #endregion
        sec_rep.set("Extracting MS/MS")
        sec_bar.step(1)
        #region Extracting MS/MS (Modified ver)
        self.msms_dict = {}
        target_id_list = self.target_table['ID'].values.tolist()
        total = set_each_prgbar(each_bar, target_id_list)
        for i, target_id in enumerate(target_id_list, start=1):
            each_rep.set("{}/{}".format(i, total))
            df = self.target_table[self.target_table['ID'] == target_id]
            pair_list = str(df['MS/MS spectrum'].values[0]).split(' ')
            fragment_list = [pair for pair in pair_list if pair != '']
            mz_list = [float(v.split(':')[0]) for v in fragment_list]
            intensity_list = [int(v.split(':')[1]) for v in fragment_list]
            msms_df = pd.DataFrame({'frag m/z': mz_list, 'intensity': intensity_list})
            self.msms_dict[target_id] = msms_df
            each_bar.step(1)
        #endregion
        sec_rep.set("Constructing Structure Database")
        sec_bar.step(1)
        #region Construcing lipid structural info dict
        self.lipid_structural_info_dict = {}
        self.target_table['Precise m/z'] = 0
        self.target_table['Precise m/z type'] = ''
        total = set_each_prgbar(each_bar, self.target_table)
        for i, (idx, one_df) in enumerate(self.target_table.iterrows(), start=1):
            each_rep.set("{}/{}".format(i, total))
            table_id = one_df['ID']
            current_structural_dict = extract_lipid_structural_info(one_df)
            msms_df = self.msms_dict[table_id]
            ref_mz = current_structural_dict['Ref precursor Mz']
            ms1_mz = one_df['Precursor m/z']
            ms2_mz = 0
            ref_front = ref_mz - 0.01
            ref_tail = ref_mz + 0.01
            ms2_df = msms_df[(msms_df['frag m/z']>=ref_front)
                            &(msms_df['frag m/z']<=ref_tail)]
            if len(ms2_df) > 0:
                ms2_mz = ms2_df['frag m/z'].values[0]
            ms1_ppm = (ms1_mz-ref_mz)/ref_mz*1000*1000
            ms2_ppm = (ms2_mz-ref_mz)/ref_mz*1000*1000
            ms1_ppm = abs(math_floor(ms1_ppm, 2))
            ms2_ppm = abs(math_floor(ms2_ppm, 2))
            mz_type = ''
            determined_mz = 0
            if ms1_ppm <= 10:
                mz_type, determined_mz = 'MS1', math_floor(ms1_mz, 4)
            elif ms2_ppm <= 10:
                mz_type, determined_mz = 'MS2', math_floor(ms2_mz, 4)
            else:
                ref_O_front = ref_mz + ex_mass['O'] - 0.01
                ref_O_tail = ref_mz + ex_mass['O'] + 0.01
                ms2_O_df = msms_df[(msms_df['frag m/z']>=ref_O_front)
                                 & (msms_df['frag m/z']<=ref_O_tail)]
                len_ms2_O_df = len(ms2_O_df)
                if len_ms2_O_df > 0:
                    ms2_O_mz = ms2_O_df['frag m/z'].values[0] - ex_mass['O']
                    ms2_O_ppm = (ms2_O_mz-ref_mz)/ref_mz*1000*1000
                    ms2_O_ppm = abs(math_floor(ms2_O_ppm, 2))
                    if ms2_O_ppm <= 10:
                        mz_type, determined_mz = 'MS2+O', math_floor(ms2_O_mz, 4)
                else:
                    mz_type, determined_mz = 'MS1>10ppm', math_floor(ms1_mz, 4)
            self.target_table.loc[idx:idx, ['Precise m/z type', 'Precise m/z']] \
            = mz_type, determined_mz
            current_structural_dict['Precise precursor Mz'] \
            = [mz_type, determined_mz]
            current_structural_dict['MS2 Mz'] = math_floor(ms2_mz, 4)
            self.lipid_structural_info_dict[table_id] = current_structural_dict
            each_bar.step(1)
        #endregion
        sec_rep.set("MS/MS Data processing")
        sec_bar.step(1)
        #region Calcurating Ratio(%) and Delta (Modified ver)
        digit, up = 4, 1
        total = set_each_prgbar(each_bar, self.msms_dict)
        for i, (idx, mz_int_df) in enumerate(self.msms_dict.items(), start=1):
            each_rep.set("{}/{}".format(i, total))
            max_int = mz_int_df['intensity'].max()
            ref_mz = self.lipid_structural_info_dict[idx]['Ref precursor Mz']
            det_mz = self.lipid_structural_info_dict[idx]['Precise precursor Mz'][1]
            mzs = mz_int_df['frag m/z'].values.tolist()
            ints = mz_int_df['intensity'].values.tolist()
            deltas = [math_floor((det_mz - mz), digit) for mz in mzs]
            ratios = [math_floor((v/max_int)*100, digit) for v in ints]
            new_msms_df = pd.DataFrame({'frag m/z': mzs, 'intensity': ints,
                                        'Delta': deltas, 'Ratio(%)': ratios})
            new_msms_df = new_msms_df[new_msms_df['frag m/z'] <= (ref_mz + up)]
            new_msms_df = new_msms_df.iloc[::-1].reset_index(drop=True)
            self.msms_dict[idx] = new_msms_df
            each_bar.step(1)
        #endregion
        sec_rep.set("Searching CID fragment ions")
        sec_bar.step(1)
        #region CID fragment ions search
        self.cid_result_dict = {}
        total = set_each_prgbar(each_bar, self.lipid_structural_info_dict)
        for i, (idx, structure_dict) in enumerate(
            self.lipid_structural_info_dict.items(), start=1):
            each_rep.set("{}/{}".format(i, total))
            msms_df = self.msms_dict[idx]
            self.cid_result_dict[idx] = search_cid_fragment_ions(
                structure_dict, msms_df, ms_tolerance_ppm)
            each_bar.step(1)
        #endregion
        sec_rep.set("Updating Structure Database")
        sec_bar.step(1)
        #region Plasmalogen search
        plasm_candidate_df = self.target_table[self.target_table['Ontology'].str.contains('Ether')]
        total = set_each_prgbar(each_bar, plasm_candidate_df)
        for i, (row, df) in enumerate(plasm_candidate_df.iterrows(), start=1):
            each_rep.set("{}/{}".format(i, total))
            table_id = df['ID']
            subclass_result_dict = self.cid_result_dict[table_id]['Lipid subclass']
            if subclass_result_dict:
                plasm_ions = [key for key, v in subclass_result_dict.items() \
                              if 'Plasmalogen' in key and v[1]>0]
                if plasm_ions:
                    structure_dict = self.lipid_structural_info_dict[table_id]
                    db_1 = structure_dict['Each moiety info']['db-1']
                    if db_1 > 0:
                        name = df['Metabolite name']
                        ontology = df['Ontology'].replace('Ether', 'Plasm')
                        self.target_table.loc[row:row, 'Ontology'] = ontology
                        structure_dict['Ontology'] = ontology
                        structure_dict['Each moiety info']['db-1'] = db_1 -1
                        if (db_1 -1 == 0):
                            pre_v = structure_dict['Unsaturated moiety']
                            structure_dict['Unsaturated moiety'] = pre_v -1
                        self.lipid_structural_info_dict[table_id] = structure_dict
            each_bar.step(1)
        #endregion
        sec_rep.set("Analyzing OAD-MS/MS")
        sec_bar.step(1)
        #region OAD analysis
        self.oad_result_dict = {}
        self.determined_db_pos_dict = {}
        total = set_each_prgbar(each_bar, self.lipid_structural_info_dict)
        for i, (table_id, info_dict) in enumerate(
            self.lipid_structural_info_dict.items(), start=1):
            name = self.target_table[
                self.target_table['ID'] == table_id]['Metabolite name'].values[0]
            name = name.split('|')[1] if '|' in name else name
            each_rep.set("{}/{} : {}".format(i, total, name))
            unsaturated_moieties_num = info_dict['Unsaturated moiety']
            ontology = info_dict['Ontology']
            solved_moiety = len(info_dict['Each moiety info'])
            is_moiety_solved = True
            if ontology in triacyls_sphingolipids and solved_moiety < 6:
                is_moiety_solved = False
            elif ontology not in mono_acyl_list and solved_moiety < 4:
                is_moiety_solved = False
            if unsaturated_moieties_num > 0 and is_moiety_solved:
                msms_df = self.msms_dict[table_id]
                self.oad_result_dict[table_id] = determine_db_positions(
                    unsaturated_moieties_num,  info_dict, msms_df, 
                    ms_tolerance_ppm, must_nl_cut_off_dict, cut_off_ratio)
            else:
                self.oad_result_dict[table_id] = {
                    'Resolved level': 'None', 
                    'Validated num': 0, 
                    'Each bools': [False], 
                    'Moiety-1': {}
                }
            each_bar.step(1)
        #endregion
        sec_rep.set("Reflecting Analysis Results")
        sec_bar.step(1)
        #region Update OAD analysis result into metabolite name & Generating Graph dict
        self.target_table['OAD result name'] = ''
        self.target_table['Solved level'] = 'None'
        self.graph_dict = {}
        total = set_each_prgbar(each_bar, self.oad_result_dict)
        for i, (table_id, result_dict) in enumerate(self.oad_result_dict.items(), start=1):
            each_rep.set("{}/{}".format(i, total))
            idx = self.target_table[self.target_table['ID'] == table_id].index[0]
            lipid_info = self.lipid_structural_info_dict[table_id]
            name = self.target_table['Metabolite name'][idx]
            oad_result_name = name.split('|')[1] if '|' in name else name
            if any(result_dict['Each bools']):
                oad_result_name = determine_oad_metabolite_name_N_description(
                    oad_result_dict=result_dict, structure_dict=lipid_info,
                    metabolite_name=oad_result_name)
            self.target_table.loc[idx:idx, ['OAD result name', 'Solved level']] \
                = oad_result_name, result_dict['Resolved level']
            self.graph_dict[table_id] = {'OAD': {}}
            self.graph_dict[table_id]['OAD'] = set_oad_graph_dict_value(
                                                result_dict, lipid_info)
            each_bar.step(1)
        #endregion
        sec_rep.set("Finalizing")
        sec_bar.step(1)
        each_rep.set("")
        #region Generating analysis result files
        self.target_table = self.target_table.sort_values(
            ['OAD result name', 'Precursor m/z'])
        dataframe_path = f'{directry}/{prefix}_analysis_table.pkl'
        msms_path = f'{directry}/{prefix}_extracted_msms.pickle'
        cid_reslut_path = f'{directry}/{prefix}_cid_result.pickle'
        oad_result_path = f'{directry}/{prefix}_oad_result.pickle'
        structure_info_path = f'{directry}/{prefix}_structure_info.pickle'
        graph_info_path = f'{directry}/{prefix}_graph_info.pickle'
        self.target_table.to_pickle(dataframe_path)
        with open(msms_path, 'wb') as output_msms_file:
            pickle.dump(self.msms_dict, output_msms_file)
        with open(cid_reslut_path, 'wb') as output_cid_file:
            pickle.dump(self.cid_result_dict, output_cid_file)
        with open(oad_result_path, 'wb') as output_oad_file:
            pickle.dump(self.oad_result_dict, output_oad_file)
        with open(structure_info_path, 'wb') as output_lipif_file:
            pickle.dump(self.lipid_structural_info_dict, output_lipif_file)
        with open(graph_info_path, 'wb') as output_graph_file:
            pickle.dump(self.graph_dict, output_graph_file)
        #endregion
        sec_rep.set("Analysis Finished")
        sec_bar.step(1)

class BatchAnalyzer(object):
    """ Definition of class that generates data object to analyze, visualize,
        modify, and save OAD-MS/MS spectra and annotations for batch data 
        format.
    """
    def __init__(self, directry, prefix, alignment_path, peakid_path, 
        peaklists_dict, ms_tolerance_ppm, must_nl_cut_off_dict, cut_off_ratio, 
        normalized, sec_rep, sec_bar, each_rep, each_bar, timer):
        """
        Args:
            
        """
        sec_rep.set("Data Pre-processing")
        sec_bar["maximum"] = 9
        #region Data preprocessing
        #Alignment table
        raw_df_1 = pd.read_csv(alignment_path, skiprows=[0,1,2,3], sep='\t')
        raw_df_1 = raw_df_1.rename(columns={
                'Average Rt(min)': 'RT(min)', 'Average Mz': 'Precursor m/z'})
        align_df = get_annotated_df(raw_df_1)
        #PeakID table
        raw_df_2 = pd.read_csv(peakid_path, skiprows=[0,1,2,3], sep='\t')
        raw_df_2 = raw_df_2.rename(columns={
                'Average Rt(min)': 'RT(min)', 'Average Mz': 'Precursor m/z'})
        peakid_df = get_annotated_df(raw_df_2)
        #PeakLists
        peaklists = [ [name, pd.read_csv(path, sep='\t')]
            for name, path in peaklists_dict.items()]
        peaklists_df_d = {li[0]: li[1] for li in peaklists}
        #region Generating target_table
        ion_v_col = 'pmol/mg tissue' if normalized else 'Height'
        new_cols = ['ID', 'Alignment ID', 'PeakID', 'RT(min)', 
                    'Precursor m/z', 'Metabolite name', 'Adduct type', 
                    ion_v_col, 'Data from', 'Reference RT', 
                    'Reference m/z', 'Formula', 'Ontology','INCHIKEY', 
                    'SMILES', 'Comment']    
        ids, align_ids, peakids, rts, mzs, names, adducts, ion_values, datafrom, \
        ref_rts, ref_mzs, formulas, ontologies, inchikeys, smiles, comments \
        = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
        count = 0
        for idx, one_df in align_df.iterrows():
            align_id = one_df['Alignment ID']
            name, adduct = one_df['Metabolite name'], one_df['Adduct type']
            inchikey, smile = one_df['INCHIKEY'], one_df['SMILES']
            formula, ontology = one_df['Formula'], one_df['Ontology']
            ref_rt , ref_mz = one_df['Reference RT'], one_df['Reference m/z']
            comment = one_df['Comment']
            for li in peaklists:
                sample, peak_df = li[0], li[1]
                ex_peaklist = peakid_df[peakid_df['Alignment ID']==align_id]
                peak_id = ex_peaklist[sample].values[0]
                if peak_id >= 0: 
                    ex_peak_df = peak_df[peak_df['PeakID']==peak_id]
                    msms = str(ex_peak_df['MSMS spectrum'].values[0])
                    if msms != 'nan':
                        # height = ex_peak_df['Height'].values[0]
                        ion_value = math_floor(one_df[sample], 3)
                        rt = ex_peak_df['RT (min)'].values[0]
                        mz = ex_peak_df['Precursor m/z'].values[0]
                        ids.append(count), align_ids.append(align_id), 
                        peakids.append(peak_id), rts.append(rt), mzs.append(mz)
                        names.append(name), adducts.append(adduct)
                        ion_values.append(ion_value), datafrom.append(sample)
                        ref_rts.append(ref_rt), ref_mzs.append(ref_mz)
                        formulas.append(formula), ontologies.append(ontology)
                        inchikeys.append(inchikey), smiles.append(smile)
                        comments.append(comment)
                        count += 1
        new_table = pd.DataFrame({'ID': ids, 'Alignment ID': align_ids,
            'PeakID': peakids, 'RT(min)': rts, 'Precursor m/z': mzs,
            'Metabolite name': names, 'Adduct type': adducts, 
            ion_v_col: ion_values, 
            'Data from': datafrom, 'Reference RT': ref_rts, 
            'Reference m/z': ref_mzs, 'Formula': formulas, 
            'Ontology': ontologies, 'INCHIKEY': inchikeys, 'SMILES':smiles, 
            'Comment': comments})
        #endregion
        self.target_table = new_table
        self.target_table.fillna({'Comment': ''})
        comments = self.target_table['Comment'].values.tolist()
        comments = ['' if isinstance(v, float) else v for v in comments]
        ex_list = ['isotope of', 'adduct linked to', 'Unit',
                   'similar chromatogram', 'found in']
        user_comment_list = [refine_comments(txt, ex_list) for txt in comments]
        self.target_table['User comment'] = user_comment_list
        #endregion
        sec_rep.set("Extracting MS/MS")
        sec_bar.step(1)
        #region Extracting MS/MS (Modified ver)
        self.msms_dict = {}
        target_id_list = self.target_table['ID'].values.tolist()
        total = set_each_prgbar(each_bar, target_id_list)
        for i, (row, one_df) in enumerate(self.target_table.iterrows(), start=1):
            each_rep.set("{}/{}".format(i, total))
            peak_id, sample = one_df['PeakID'], one_df['Data from']
            target_id = one_df['ID']
            ex_df = peaklists_df_d[sample]
            ex_df = ex_df[ex_df['PeakID']==peak_id]
            pair_list = str(ex_df['MSMS spectrum'].values[0]).split(' ')
            fragment_list = [pair for pair in pair_list if pair != '']
            mz_list = [float(v.split(':')[0]) for v in fragment_list]
            intensity_list = [int(v.split(':')[1]) for v in fragment_list]
            msms_df = pd.DataFrame({'frag m/z': mz_list, 
                                    'intensity': intensity_list})
            self.msms_dict[target_id] = msms_df
            each_bar.step(1)
        #endregion
        sec_rep.set("Constructing Structure Database")
        sec_bar.step(1)
        #region Construcing lipid structural info dict
        self.lipid_structural_info_dict = {}
        self.target_table['Precise m/z'] = 0
        self.target_table['Precise m/z type'] = ''
        total = set_each_prgbar(each_bar, self.target_table)
        for i, (idx, one_df) in enumerate(self.target_table.iterrows(), start=1):
            each_rep.set("{}/{}".format(i, total))
            table_id = one_df['ID']
            current_structural_dict = extract_lipid_structural_info(one_df)
            msms_df = self.msms_dict[table_id]
            ref_mz = current_structural_dict['Ref precursor Mz']
            ms1_mz = one_df['Precursor m/z']
            ms2_mz = 0
            ref_front = ref_mz - 0.01
            ref_tail = ref_mz + 0.01
            ms2_df = msms_df[(msms_df['frag m/z']>=ref_front)
                            &(msms_df['frag m/z']<=ref_tail)]
            if len(ms2_df) > 0:
                ms2_mz = ms2_df['frag m/z'].values[0]
            ms1_ppm = (ms1_mz-ref_mz)/ref_mz*1000*1000
            ms2_ppm = (ms2_mz-ref_mz)/ref_mz*1000*1000
            ms1_ppm = abs(math_floor(ms1_ppm, 2))
            ms2_ppm = abs(math_floor(ms2_ppm, 2))
            mz_type = ''
            determined_mz = 0
            if ms1_ppm <= 10:
                mz_type, determined_mz = 'MS1', math_floor(ms1_mz, 4)
            elif ms2_ppm <= 10:
                mz_type, determined_mz = 'MS2', math_floor(ms2_mz, 4)
            else:
                ref_O_front = ref_mz + ex_mass['O'] - 0.01
                ref_O_tail = ref_mz + ex_mass['O'] + 0.01
                ms2_O_df = msms_df[(msms_df['frag m/z']>=ref_O_front)
                                  &(msms_df['frag m/z']<=ref_O_tail)]
                len_ms2_O_df = len(ms2_O_df)
                if len_ms2_O_df > 0:
                    ms2_O_mz = ms2_O_df['frag m/z'].values[0] - ex_mass['O']
                    ms2_O_ppm = (ms2_O_mz-ref_mz)/ref_mz*1000*1000
                    ms2_O_ppm = abs(math_floor(ms2_O_ppm, 2))
                    if ms2_O_ppm <= 10:
                        mz_type, determined_mz = 'MS2+O', math_floor(ms2_O_mz, 4)
                else:
                    mz_type, determined_mz = 'MS1>10ppm', math_floor(ms1_mz, 4)
            self.target_table.loc[idx:idx, ['Precise m/z type', 'Precise m/z']] \
            = mz_type, determined_mz
            current_structural_dict['Precise precursor Mz'] \
            = [mz_type, determined_mz]
            current_structural_dict['MS2 Mz'] = math_floor(ms2_mz, 4)
            self.lipid_structural_info_dict[table_id] = current_structural_dict
            each_bar.step(1)
        #endregion
        sec_rep.set("MS/MS Data processing")
        sec_bar.step(1)
        #region Calcurating Ratio(%) and Delta (Modified ver)
        digit, up = 4, 1
        total = set_each_prgbar(each_bar, self.msms_dict)
        for i, (idx, mz_int_df) in enumerate(self.msms_dict.items(), start=1):
            each_rep.set("{}/{}".format(i, total))
            max_int = mz_int_df['intensity'].max()
            ref_mz = self.lipid_structural_info_dict[idx]['Ref precursor Mz']
            det_mz = self.lipid_structural_info_dict[idx]['Precise precursor Mz'][1]
            mzs = mz_int_df['frag m/z'].values.tolist()
            ints = mz_int_df['intensity'].values.tolist()
            deltas = [math_floor((det_mz - mz), digit) for mz in mzs]
            ratios = [math_floor((v/max_int)*100, digit) for v in ints]
            new_msms_df = pd.DataFrame({'frag m/z': mzs, 'intensity': ints,
                                        'Delta': deltas, 'Ratio(%)': ratios})
            new_msms_df = new_msms_df[new_msms_df['frag m/z'] <= (ref_mz + up)]
            new_msms_df = new_msms_df.iloc[::-1].reset_index(drop=True)
            self.msms_dict[idx] = new_msms_df
            each_bar.step(1)
        #endregion
        sec_rep.set("Searching CID fragment ions")
        sec_bar.step(1)
        #region CID fragment ions search
        self.cid_result_dict = {}
        total = set_each_prgbar(each_bar, self.lipid_structural_info_dict)
        for i, (idx, structure_dict) in enumerate(self.lipid_structural_info_dict.items(), start=1):
            each_rep.set("{}/{}".format(i, total))
            msms_df = self.msms_dict[idx]
            try:
                self.cid_result_dict[idx] = search_cid_fragment_ions(
                    structure_dict, msms_df, ms_tolerance_ppm)
            except:
                print(f'ID: {idx}')
                pprint.pprint(structure_dict)
            each_bar.step(1)
        #endregion
        sec_rep.set("Updating Structure Database")
        sec_bar.step(1)
        #region Plasmalogen search
        plasm_candidate_df = self.target_table[self.target_table['Ontology'].str.contains('Ether')]
        total = set_each_prgbar(each_bar, plasm_candidate_df)
        for i, (row, df) in enumerate(plasm_candidate_df.iterrows(), start=1):
            each_rep.set("{}/{}".format(i, total))
            table_id = df['ID']
            subclass_result_dict = self.cid_result_dict[table_id]['Lipid subclass']
            if subclass_result_dict:
                plasm_ions = [key for key, v in subclass_result_dict.items() \
                              if 'Plasmalogen' in key and v[1]>0]
                if plasm_ions:
                    structure_dict = self.lipid_structural_info_dict[table_id]
                    db_1 = structure_dict['Each moiety info']['db-1']
                    if db_1 > 0:
                        name = df['Metabolite name']
                        ontology = df['Ontology'].replace('Ether', 'Plasm')
                        self.target_table.loc[row:row, 'Ontology'] = ontology
                        structure_dict['Ontology'] = ontology
                        structure_dict['Each moiety info']['db-1'] = db_1 -1
                        if (db_1 -1 == 0):
                            pre_v = structure_dict['Unsaturated moiety']
                            structure_dict['Unsaturated moiety'] = pre_v -1
                        self.lipid_structural_info_dict[table_id] = structure_dict
            each_bar.step(1)
        #endregion
        sec_rep.set("Analyzing OAD-MS/MS")
        sec_bar.step(1)
        #region OAD analysis
        self.oad_result_dict = {}
        self.determined_db_pos_dict = {}
        total = set_each_prgbar(each_bar, self.lipid_structural_info_dict)
        for i, (table_id, info_dict) in enumerate(
            self.lipid_structural_info_dict.items(), start=1):
            name = self.target_table[self.target_table['ID'] == table_id]['Metabolite name'].values[0]
            name = name.split('|')[1] if '|' in name else name
            each_rep.set("{}/{} : {}".format(i, total, name))
            unsaturated_moieties_num = info_dict['Unsaturated moiety']
            ontology = info_dict['Ontology']
            solved_moiety = len(info_dict['Each moiety info'])
            is_moiety_solved = True
            if ontology in triacyls_sphingolipids and solved_moiety < 6:
                is_moiety_solved = False
            elif ontology not in mono_acyl_list and solved_moiety < 4:
                is_moiety_solved = False
            if unsaturated_moieties_num > 0 and is_moiety_solved:
                msms_df = self.msms_dict[table_id]
                try:
                    self.oad_result_dict[table_id] = determine_db_positions(
                        unsaturated_moieties_num,  info_dict, msms_df, 
                        ms_tolerance_ppm, must_nl_cut_off_dict, cut_off_ratio)
                except:
                    print(name)
                    pprint.pprint(info_dict)
            else:
                self.oad_result_dict[table_id] = {'Resolved level': 'None', 
                    'Validated num': 0, 'Each bools': [False], 'Moiety-1': {}}
            each_bar.step(1)
        #endregion
        sec_rep.set("Reflecting Analysis Results")
        sec_bar.step(1)
        #region Update OAD analysis result into metabolite name & Generating Graph dict
        self.target_table['OAD result name'] = ''
        self.target_table['Solved level'] = 'None'
        self.graph_dict = {}
        total = set_each_prgbar(each_bar, self.oad_result_dict)
        for i, (table_id, result_dict) in enumerate(self.oad_result_dict.items(), start=1):
            each_rep.set("{}/{}".format(i, total))
            idx = self.target_table[self.target_table['ID'] == table_id].index[0]
            lipid_info = self.lipid_structural_info_dict[table_id]
            name = self.target_table['Metabolite name'][idx]
            oad_result_name = name.split('|')[1] if '|' in name else name
            if any(result_dict['Each bools']):
                oad_result_name = determine_oad_metabolite_name_N_description(
                    oad_result_dict=result_dict, structure_dict=lipid_info,
                    metabolite_name=oad_result_name)
            self.target_table.loc[idx:idx, ['OAD result name', 'Solved level']] \
                = oad_result_name, result_dict['Resolved level']
            self.graph_dict[table_id] = {'OAD': {}}
            self.graph_dict[table_id]['OAD'] \
                = set_oad_graph_dict_value(result_dict, lipid_info)
            each_bar.step(1)
        #endregion
        sec_rep.set("Finalizing")
        sec_bar.step(1)
        each_rep.set("")
        #region Generating analysis result files
        self.target_table = self.target_table.sort_values(
            ['OAD result name', 'Precursor m/z'])
        dataframe_path = f'{directry}/{prefix}_analysis_table.pkl'
        msms_path = f'{directry}/{prefix}_extracted_msms.pickle'
        cid_reslut_path = f'{directry}/{prefix}_cid_result.pickle'
        oad_result_path = f'{directry}/{prefix}_oad_result.pickle'
        structure_info_path = f'{directry}/{prefix}_structure_info.pickle'
        graph_info_path = f'{directry}/{prefix}_graph_info.pickle'
        self.target_table.to_pickle(dataframe_path)
        with open(msms_path, 'wb') as output_msms_file:
            pickle.dump(self.msms_dict, output_msms_file)
        with open(cid_reslut_path, 'wb') as output_cid_file:
            pickle.dump(self.cid_result_dict, output_cid_file)
        with open(oad_result_path, 'wb') as output_oad_file:
            pickle.dump(self.oad_result_dict, output_oad_file)
        with open(structure_info_path, 'wb') as output_lipif_file:
            pickle.dump(self.lipid_structural_info_dict, output_lipif_file)
        with open(graph_info_path, 'wb') as output_graph_file:
            pickle.dump(self.graph_dict, output_graph_file)
        #endregion
        sec_rep.set("Analysis Finished")
        sec_bar.step(1)

#region Utilities
def get_annotated_df(raw_df):
    """ Extract annotated metabolites in given dataframe

    Args:
        raw_df (Dataframe): Dataframe as MS-DIAL format.
    
    Returns:
        df (Dataframe): Returns the dataframe only with annotated metabolites
    """
    df = raw_df.dropna(subset=['Metabolite name'])
    df = df.dropna(subset=['MS/MS spectrum'])
    df = df[df['Metabolite name'] != 'Unknown']
    df = df[~df['Metabolite name'].str.startswith('w/o')]
    df = df[~df['Metabolite name'].str.startswith('RIKEN')]
    df = df[~df['Metabolite name'].str.startswith('Unsettled')]
    for subclass in exclude_subclass:
        df = df[df['Ontology'] != subclass]
    return df

def refine_comments(txt, excluding_list):
    """ Exclude the unwanted words in given string

    Args:
        txt (str): target string
        excluding_list (list): list of unwated words

    Returns:
        txt (str): Returns the string without unwanted words
    """
    if txt == '': return txt
    else:
        if '; ' in txt: txt = txt.split('; ')[0]
        for tag in excluding_list:
            if tag in txt: return ''
        return txt

def set_each_prgbar(prgbar, li):
    """ Set 0 and total value to progress bar

    Args:
        prgbar (object): Progressbar object
        li (list): list of target metabolites

    Returns:
        int: Returns total values (int)
    """
    prgbar["value"] = 0
    prgbar["maximum"] = len(li)
    return len(li)

#region Not used function
def each_process_timer(label, start):
    hour, minute, second = 0, 0, 0
    while start:
        second += 1
        time.sleep(1)
        if second == 60:
            minute += 1
            second = 0
            if minute == 60:
                minute = 0
                hour += 1
        if minute >= 10 and second >= 10:
            label.set("0{h}:{m}:{s}".format(h=str(hour),
            m=str(minute), s=str(second)))
        elif minute < 10 and second >= 10:
            label.set("0{h}:0{m}:{s}".format(h=str(hour),
            m=str(minute), s=str(second)))
        elif minute >= 10 and second < 10:
            label.set("0{h}:{m}:0{s}".format(h=str(hour),
            m=str(minute), s=str(second)))
        else:
            label.set("0{h}:0{m}:0{s}".format(h=str(hour),
            m=str(minute), s=str(second)))
    if not start:
        hour, minute, second = 0, 0, 0
#endregion

def math_floor(num, digit):
    floored = math.floor(num*10**digit)/(10**digit)
    return floored
#endregion

#region Graph dict
def set_oad_graph_dict_value(oad_dict, lipid_info):
    """ Generate dictionary for visualizing graphs of OAD-MS/MS spectra

    Args:
        oad_dict (dict): dict contains below
            Resolved level (str): 'All', 'Partial', 'None'
            Validated num (int): 0-3
            Each bools (list[bool]): Whether each moiety was resolved
            Moiety-1~3 (dict): Contains 'Score', 'Presence', and 'Ratio sum'
        lipid_info (dict): dict containing belows
            'Status', 'Adduct', 'Precursor Mz', 'MS2 Mz', 'Precise precursor Mz',
            'Ref precursor Mz', 'RT(min)', 'Reference RT', 'Ontology', 'Brutto', 
            'Valid moiety num', 'Each moiety info', 'Unsaturated moiety', 
            'Unsaturated sphingobase', 'Deuterium', 'SMILES', 'Formula', 
            'Atom dict', 'Oxidized', 'NL type'

    Returns:
        graph_dict (dict):  Returns dict containing MS2 Mz, Ref precursor Mz,
                            Ontology, x-range, Magnification, Bar_width
    """
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
            measured_oad_ions_1 = [li[1] for li in d[exp_peaks] if li[1] > 0]
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
    final_magnification= 10
    max_criteria = 30
    if measured_oad_ions_3: max_ratio_list.append(max(measured_oad_ions_3))
    if measured_oad_ions_2: max_ratio_list.append(max(measured_oad_ions_2))
    if measured_oad_ions_1: max_ratio_list.append(max(measured_oad_ions_1))
    if max_ratio_list:
        final_magnification = math.floor(max_criteria/max(max_ratio_list))
        if final_magnification > 500:
            final_magnification = 500
    #endregion
    graph_dict['MS2 Mz'] = precursor_mz
    graph_dict['Ref precursor Mz'] = ref_precursor_mz
    graph_dict['Ontology'], graph_dict['x-range'] = ontology, [x_min, x_max]
    graph_dict['Magnification'] = final_magnification
    graph_dict['Bar_width'] = bar_width
    return graph_dict
#endregion

#region lipid structural info dict
def extract_lipid_structural_info(df):
    """ Generate lipid_info_d (dict) as basic information for OAD-MS/MS analysis

    Args:
        df (Dataframe): input dataframe of one metabolite as MS-DIAL format

    Returns:
        lipid_info_d (dict):  Returns dict containing belows
            'Status', 'Adduct', 'Precursor Mz', 'MS2 Mz', 'Precise precursor Mz',
            'Ref precursor Mz', 'RT(min)', 'Reference RT', 'Ontology', 'Brutto', 
            'Valid moiety num', 'Each moiety info', 'Unsaturated moiety', 
            'Unsaturated sphingobase', 'Deuterium', 'SMILES', 'Formula', 
            'Atom dict', 'Oxidized', 'NL type'
    """
    lipid_info_d = {
        'Status': '', 'Adduct': '', 'Precursor Mz': 0, 'MS2 Mz': 0,
        'Precise precursor Mz': '', 'Ref precursor Mz': 0, 
        'RT(min)': 0, 'Reference RT': 0, 'Ontology': '',
        'Brutto': '', 'Valid moiety num': 0, 'Each moiety info': {}, 
        'Unsaturated moiety': 0, 'Unsaturated sphingobase': False, 
        'Deuterium': 0,
        'SMILES': '', 'Formula': '', 'Atom dict': '', 'Oxidized': 0, 
        'NL type': []}
    metabolite_name = str(df['Metabolite name'])
    lipid_info_d['Status'] = 'Identified'
    lipid_info_d['Adduct'] = str(df['Adduct type'])
    lipid_info_d['Precursor Mz'] = math_floor(float(df['Precursor m/z']), 4)
    lipid_info_d['Ref precursor Mz'] = math_floor(float(df['Reference m/z']), 4)
    lipid_info_d['RT(min)'] = math_floor(float(df['RT(min)']), 3)
    lipid_info_d['Reference RT'] = math_floor(float(df['Reference RT']), 3)
    ontology = str(df['Ontology'])
    unsaturated_sphingobase = False
    if ontology != 'EtherPE':
        lipid_info_d['Ontology'] = ontology
    elif (ontology == 'EtherPE') and 'P-' in metabolite_name:
        lipid_info_d['Ontology'] = 'PlasmPE'
    else:
        lipid_info_d['Ontology'] = ontology
    lipid_info_d['SMILES'] = str(df['SMILES'])
    lipid_info_d['Formula'] = str(df['Formula'])
    is_deuteriumed = re.findall(r'\(d\d\)', metabolite_name)
    if is_deuteriumed:
        lipid_info_d['Deuterium'] = int(re.findall(r'\d', is_deuteriumed[0])[0])
            
    if ('|' in metabolite_name):
        split_name = metabolite_name.split('|')
        moiety_list = re.findall(r'\d+\:\d+', split_name[1])
        lipid_info_d['Valid moiety num'] = len(moiety_list)
        moiety_dict, unsaturated_moiety_num, unsaturated_sphingobase \
        = get_chain_and_db_value(moiety_list, ontology)
        lipid_info_d['Unsaturated moiety'] = unsaturated_moiety_num
        lipid_info_d['Each moiety info'] = moiety_dict
    else:
        if ontology in no_acyl_list:
            pass
        elif ontology in mono_acyl_list:
            moiety_list = re.findall(r'\d+\:\d+', metabolite_name)
            if len(moiety_list) > 1:
                if ontology in lipidclass_dict['Sterol lipids']:
                    moiety = moiety_list[1]
                    moiety_list = [moiety]
                else:   
                    moiety_list.remove('0:0')
            lipid_info_d['Valid moiety num'] = len(moiety_list)
            moiety_dict, unsaturated_moiety_num, unsaturated_sphingobase \
            = get_chain_and_db_value(moiety_list, ontology)
            lipid_info_d['Unsaturated moiety'] = unsaturated_moiety_num
            lipid_info_d['Each moiety info'] = moiety_dict
        else:
            split_name = metabolite_name.split()
            moiety_list = re.findall(r'\d+\:\d+', split_name[1])
            lipid_info_d['Valid moiety num'] = len(moiety_list)
            moiety_dict, unsaturated_moiety_num, unsaturated_sphingobase \
            = get_chain_and_db_value(moiety_list, ontology)
            lipid_info_d['Unsaturated moiety'] = unsaturated_moiety_num
            lipid_info_d['Each moiety info'] = moiety_dict
    try:
        chain_sum = 0
        db_sum = 0
        for i, v in enumerate(moiety_dict.values()):
            if i%2 == 0:
                chain_sum += int(v)
            else:
                db_sum += int(v)
        lipid_info_d['Brutto'] = str(chain_sum) + ':' + str(db_sum)
    except Exception:
        pass
    lipid_info_d['Unsaturated sphingobase'] = unsaturated_sphingobase
    return lipid_info_d

def get_chain_and_db_value(moiety_list, ontology):
    """ Extract numbers of carbons and double bonds in moieties

    Args:
        moiety_list (list): list of moieties (str)
        ontology (str): Lipid subclass

    Returns:
        moiety_dict (dict): Returns dict contains numbers of carbons and double
                            bond in each moiety
        unsaturated_moiety_num (int):   number of unsaturated moieties
        unsaturated_sphingobase (bool): Whether there is unsaturation in 
                                        sphingoid base
    """
    unsaturated_moiety_num = 0
    unsaturated_sphingobase = False
    if len(moiety_list) == 1:
        chain_and_db = moiety_list[0].split(':')
        moiety_dict = {'chain-1': int(chain_and_db[0]), 
                        'db-1': int(chain_and_db[1])}
        if moiety_dict['db-1'] > 0:
            unsaturated_moiety_num = 1
    elif len(moiety_list) == 2:
        chain_and_db_1 = moiety_list[0].split(':')
        chain_and_db_2 = moiety_list[1].split(':')
        moiety_dict = {'chain-1': int(chain_and_db_1[0]), 
                    'db-1': int(chain_and_db_1[1]),
                    'chain-2': int(chain_and_db_2[0]), 
                    'db-2': int(chain_and_db_2[1])}
        if ontology != 'AHexCer':
            moiety_dict = {'chain-1': int(chain_and_db_1[0]), 
                            'db-1': int(chain_and_db_1[1]),
                            'chain-2': int(chain_and_db_2[0]), 
                            'db-2': int(chain_and_db_2[1])}
        else:
            moiety_dict = {'chain-1': int(chain_and_db_2[0]), 
                            'db-1': int(chain_and_db_2[1]),
                            'chain-2': int(chain_and_db_1[0]), 
                            'db-2': int(chain_and_db_1[1])}
        if moiety_dict['db-1'] == 0 and moiety_dict['db-2'] == 0:
            unsaturated_moiety_num = 0
        elif moiety_dict['db-1'] == 0 or moiety_dict['db-2'] == 0:
            unsaturated_moiety_num = 1
        else:
            unsaturated_moiety_num = 2
        if (ontology in lipidclass_dict['Sphingolipids']) and (moiety_dict['db-1'] > 0):
            unsaturated_sphingobase = True
    elif len(moiety_list) == 3:
        chain_and_db_1 = moiety_list[0].split(':')
        chain_and_db_2 = moiety_list[1].split(':')
        chain_and_db_3 = moiety_list[2].split(':')
        if ontology != 'AHexCer':
            moiety_dict = {'chain-1': int(chain_and_db_1[0]), 
                            'db-1': int(chain_and_db_1[1]),
                            'chain-2': int(chain_and_db_2[0]), 
                            'db-2': int(chain_and_db_2[1]), 
                            'chain-3': int(chain_and_db_3[0]), 
                            'db-3': int(chain_and_db_3[1])}
        else:
            moiety_dict = {'chain-1': int(chain_and_db_2[0]), 
                            'db-1': int(chain_and_db_2[1]),
                            'chain-2': int(chain_and_db_3[0]), 
                            'db-2': int(chain_and_db_3[1]), 
                            'chain-3': int(chain_and_db_1[0]), 
                            'db-3': int(chain_and_db_1[1])}
        if moiety_dict['db-1'] == 0 and moiety_dict['db-2'] == 0 and moiety_dict['db-3'] == 0:
            unsaturated_moiety_num = 0
        elif ((moiety_dict['db-1'] == 0 and moiety_dict['db-2'] == 0) 
            or (moiety_dict['db-1'] == 0 and moiety_dict['db-3'] == 0)
            or (moiety_dict['db-2'] == 0 and moiety_dict['db-3'] == 0)):
            unsaturated_moiety_num = 1
        elif moiety_dict['db-1'] == 0 or moiety_dict['db-2'] == 0 or moiety_dict['db-3'] == 0:
            unsaturated_moiety_num = 2
        else:
            unsaturated_moiety_num = 3
        if (ontology in lipidclass_dict['Sphingolipids']) and (moiety_dict['db-1'] > 0):
            unsaturated_sphingobase = True
    return moiety_dict, unsaturated_moiety_num, unsaturated_sphingobase
#endregion

#region OAD fragmentation analysis
def determine_db_positions(unsaturated_moieties_num, structure_dict, 
    msms_df, ms_tolerance_ppm, must_nl_cut_off_dict, cut_off_ratio):
    """ Determine C=C positions in an inputed metabolite

    Args:
        unsaturated_moieties_num (int): number of unsaturated moieties
        structure_dict (dict): dict generated by extract_lipid_structural_info
        msms_df (Dataframe): Extracted MS/MS spectrum
        ms_tolerance_ppm (int): mass tolerance
        must_nl_cut_off_dict: 
            types of essential ions and relative intensity threshold
        cut_off_ratio (float): relative intensity threshold

    Returns:
        total_result_dict (dict): Returns dict contains below
            Resolved level (str): 'All', 'Partial', 'None'
            Validated num (int): 0-3
            Each bools (list[bool]): Whether each moiety was resolved
            Moiety-1~3 (dict): Containing 'Score', 'Presence', and 'Ratio sum'
    """
    #region parameter settings
    ref_precursor_mz = structure_dict['Ref precursor Mz']
    ontology = structure_dict['Ontology']
    ms_tolerance = math_floor(ms_tolerance_ppm*ref_precursor_mz/(1000*1000), 6)
    db_in_sphingobase = structure_dict['Unsaturated sphingobase']
    cut_off_df = msms_df[msms_df['Ratio(%)'] >= cut_off_ratio]
    score_cutoff = 0
    ref_oad_dict = {}
    ref_oad_dict = generate_reference_oad_dict(
        structure_dict, unsaturated_moieties_num, msms_df, 
        ms_tolerance, must_nl_cut_off_dict
    )
    unresolved_d = {
        0: {'Positions': '', 'N-description': 'Unresolved',
            'Score': '###', 'Ratio sum': '###', 'Presence': '###',
            'Notice': 'Unresolved', 
            'Measured peaks': [[0, 0, 0]], 
            'Ref peaks': [['', 0, 0, 0]], 
            'Peaks dict': {'none': [0, 0, 0, 0, 0]}
            }
    }
    #endregion

    if unsaturated_moieties_num == 1:
        #region settings
        ref_oad_dict_1 = ref_oad_dict['1']
        total_result_dict = {}
        db_num_1 = structure_dict['Each moiety info']['db-1']
        c_num_1 = structure_dict['Each moiety info']['chain-1']
        if db_num_1 == 0:
            db_num_1 = structure_dict['Each moiety info']['db-2']
            if db_num_1 == 0:
                db_num_1 = structure_dict['Each moiety info']['db-3']
            c_num_1 = 0
        #endregion
        """ Checking diagnostic ions """
        #region
        diagnostic_ions_result_dict_1 = {}
        for_range_1 = range(0, db_num_1)
        diagnostic_ions_result_dict_1 = query_essential_diagnostic_ions(
            df=msms_df, ref_oad_dict=ref_oad_dict_1, 
            db_in_sphingobase=db_in_sphingobase, c_num=c_num_1, db_num=db_num_1, 
            for_range=for_range_1, tolerance=ms_tolerance, 
            must_nl_cut_off_dict=must_nl_cut_off_dict, 
            structure_dict=structure_dict
        )
        #endregion
        """ Calculating presence rate and ratio sum via reference """
        #region
        sph_set = [db_in_sphingobase, c_num_1]
        score_info_dict = calc_presence_ratios_and_score(
            ref_oad_dict_1, diagnostic_ions_result_dict_1, cut_off_df, 
            ref_precursor_mz, ms_tolerance_ppm, sph_set
        )
        #endregion
        """ Sorting the result by MS/MS simularity """
        #region
        score_info_dict = score_thresholding(score_info_dict, score_cutoff)
        unresolved_d_1 = unresolved_d.copy()
        if score_info_dict:
            sorted_dict = get_score_sorted_dict(score_info_dict)
            sorted_dict[len(sorted_dict)] = unresolved_d_1[0]
            revolved_level = 'All'
            validated_num = 1
            bool_list = [True]
        else:
            sorted_dict = unresolved_d_1
            revolved_level = 'None'
            validated_num = 0
            bool_list = [False]
        sorted_dict = add_determined_db_info(sorted_dict)
        total_result_dict = {
            'Resolved level': revolved_level, 
            'Validated num': validated_num,
            'Each bools': bool_list, 
            'Moiety-1': sorted_dict
        }
        return total_result_dict
        #endregion

    if unsaturated_moieties_num == 2:
        #region settings
        ref_oad_dict_1 = ref_oad_dict['1']
        ref_oad_dict_2 = ref_oad_dict['2']
        total_result_dict = {}
        valid_moiety_num = structure_dict['Valid moiety num']
        if valid_moiety_num == 2:
            db_num_1 = structure_dict['Each moiety info']['db-1']
            c_num_1 = structure_dict['Each moiety info']['chain-1']
            db_num_2 = structure_dict['Each moiety info']['db-2']
        elif valid_moiety_num == 3:
            db_num_1 = structure_dict['Each moiety info']['db-2']
            c_num_1 = 0
            db_num_2 = structure_dict['Each moiety info']['db-3']
            if db_num_1 == 0:
                db_num_1 = structure_dict['Each moiety info']['db-1']
                c_num_1 = structure_dict['Each moiety info']['chain-1']
                db_num_2 = structure_dict['Each moiety info']['db-3']
            if db_num_2 == 0:
                db_num_1 = structure_dict['Each moiety info']['db-1']
                c_num_1 = structure_dict['Each moiety info']['chain-1']
                db_num_2 = structure_dict['Each moiety info']['db-2']
        #endregion
        """ Checking diagnostic ions """
        #region
        # Acyl-1
        diagnostic_ions_result_dict_1 = {}
        for_range_1 = range(0, db_num_1)
        diagnostic_ions_result_dict_1 = query_essential_diagnostic_ions(
            df=msms_df, ref_oad_dict=ref_oad_dict_1, 
            db_in_sphingobase=db_in_sphingobase, c_num=c_num_1, db_num=db_num_1, 
            for_range=for_range_1, tolerance=ms_tolerance, 
            must_nl_cut_off_dict=must_nl_cut_off_dict, 
            structure_dict=structure_dict)
        # Acyl-2
        diagnostic_ions_result_dict_2 = {}
        for_range_2 = range(0, db_num_2)
        diagnostic_ions_result_dict_2 = query_essential_diagnostic_ions(
            df=msms_df, ref_oad_dict=ref_oad_dict_2, 
            db_in_sphingobase=db_in_sphingobase, c_num=-1, db_num=db_num_2, 
            for_range=for_range_2, tolerance=ms_tolerance, 
            must_nl_cut_off_dict=must_nl_cut_off_dict, 
            structure_dict=structure_dict)
        #endregion
        """ Calculating presence rate and ratio sum via reference """
        #region
        # Acyl-2
        score_info_dict_2 = calc_presence_ratios_and_score(
            ref_oad_dict_2, diagnostic_ions_result_dict_2, cut_off_df, 
            ref_precursor_mz, ms_tolerance_ppm, [False, c_num_1])
        # Acyl-1
        sph_set = [db_in_sphingobase, c_num_1]
        score_info_dict_1 = calc_presence_ratios_and_score(
            ref_oad_dict_1, diagnostic_ions_result_dict_1, cut_off_df, 
            ref_precursor_mz, ms_tolerance_ppm, sph_set)
        #endregion
        """ Sorting the result by MS/MS simularity """
        #region
        score_info_dict_1 = score_thresholding(score_info_dict_1, score_cutoff)
        score_info_dict_2 = score_thresholding(score_info_dict_2, score_cutoff)
        unresolved_d_1 = unresolved_d.copy()
        unresolved_d_2 = unresolved_d.copy()
        if score_info_dict_2 and score_info_dict_1:
            sorted_dict_2 = get_score_sorted_dict(score_info_dict_2)
            sorted_dict_1 = get_score_sorted_dict(score_info_dict_1)
            sorted_dict_1, sorted_dict_2 = check_same_dbs_position_in_two(
                dict_1=sorted_dict_1, dict_2=sorted_dict_2, tag='2')
            sorted_dict_2[len(sorted_dict_2)] = unresolved_d_2[0]
            sorted_dict_1[len(sorted_dict_1)] = unresolved_d_1[0]
            revolved_level = 'All'
            validated_num = 2
            bool_list = [True, True]
        elif score_info_dict_2:
            sorted_dict_2 = get_score_sorted_dict(score_info_dict_2)
            sorted_dict_1 = unresolved_d_1
            sorted_dict_2[len(sorted_dict_2)] = unresolved_d_2[0]
            revolved_level = 'Partial'
            validated_num = 1
            bool_list = [False, True]
        elif score_info_dict_1:
            sorted_dict_2 = unresolved_d_2
            sorted_dict_1 = get_score_sorted_dict(score_info_dict_1)
            sorted_dict_1[len(sorted_dict_1)] = unresolved_d_1[0]
            revolved_level = 'Partial'
            validated_num = 1
            bool_list = [True, False]
        else:
            sorted_dict_2 = unresolved_d_2
            sorted_dict_1 = unresolved_d_1
            revolved_level = 'None'
            validated_num = 0
            bool_list = [False, False]
        sorted_dict_1 = add_determined_db_info(sorted_dict_1)
        sorted_dict_2 = add_determined_db_info(sorted_dict_2)
        total_result_dict = {'Resolved level': revolved_level,
                             'Validated num': validated_num,
                             'Each bools': bool_list,
                             'Moiety-1': sorted_dict_1,
                             'Moiety-2': sorted_dict_2}
        return total_result_dict
        #endregion

    if unsaturated_moieties_num == 3:
        #region settings
        ref_oad_dict_1 = ref_oad_dict['1']
        ref_oad_dict_2 = ref_oad_dict['2']
        ref_oad_dict_3 = ref_oad_dict['3']
        total_result_dict = {}
        db_num_1 = structure_dict['Each moiety info']['db-1']
        c_num_1 = structure_dict['Each moiety info']['chain-1']
        db_num_2 = structure_dict['Each moiety info']['db-2']
        db_num_3 = structure_dict['Each moiety info']['db-3']
        #endregion
        """ Checking diagnostic ions """
        #region 
        # Acyl-1
        diagnostic_ions_result_dict_1 = {}
        for_range_1 = range(0, db_num_1)
        diagnostic_ions_result_dict_1 = query_essential_diagnostic_ions(
            df=msms_df, ref_oad_dict=ref_oad_dict_1, 
            db_in_sphingobase=db_in_sphingobase, c_num=c_num_1, db_num=db_num_1, 
            for_range=for_range_1, tolerance=ms_tolerance, 
            must_nl_cut_off_dict=must_nl_cut_off_dict, 
            structure_dict=structure_dict)
        # Acyl-2
        diagnostic_ions_result_dict_2 = {}
        for_range_2 = range(0, db_num_2)
        diagnostic_ions_result_dict_2 = query_essential_diagnostic_ions(
            df=msms_df, ref_oad_dict=ref_oad_dict_2, 
            db_in_sphingobase=db_in_sphingobase, c_num=-1, db_num=db_num_2, 
            for_range=for_range_2, tolerance=ms_tolerance, 
            must_nl_cut_off_dict=must_nl_cut_off_dict, 
            structure_dict=structure_dict)
        # Acyl-3
        diagnostic_ions_result_dict_3 = {}
        for_range_3 = range(0, db_num_3)
        diagnostic_ions_result_dict_3 = query_essential_diagnostic_ions(
            df=msms_df, ref_oad_dict=ref_oad_dict_3, 
            db_in_sphingobase=db_in_sphingobase, c_num=-1, db_num=db_num_3, 
            for_range=for_range_3, tolerance=ms_tolerance, 
            must_nl_cut_off_dict=must_nl_cut_off_dict, 
            structure_dict=structure_dict)
        #endregion
        """ Calculating presence rate and ratio sum via reference """
        #region
        # Acyl-3
        score_info_dict_3 = calc_presence_ratios_and_score(
            ref_oad_dict_3, diagnostic_ions_result_dict_3, cut_off_df, 
            ref_precursor_mz, ms_tolerance_ppm, [False, 0])
        # Acyl-2
        score_info_dict_2 = calc_presence_ratios_and_score(
            ref_oad_dict_2, diagnostic_ions_result_dict_2, cut_off_df, 
            ref_precursor_mz, ms_tolerance_ppm, [False, 0])
        # Acyl-1
        sph_set = [db_in_sphingobase, c_num_1]
        score_info_dict_1 = calc_presence_ratios_and_score(
            ref_oad_dict_1, diagnostic_ions_result_dict_1, cut_off_df, 
            ref_precursor_mz, ms_tolerance_ppm, sph_set)
        #endregion
        """ Sorting the result by MS/MS simularity """
        #region
        score_info_dict_1 = score_thresholding(score_info_dict_1, score_cutoff)
        score_info_dict_2 = score_thresholding(score_info_dict_2, score_cutoff)
        score_info_dict_3 = score_thresholding(score_info_dict_3, score_cutoff)
        unresolved_d_1 = unresolved_d.copy()
        unresolved_d_2 = unresolved_d.copy()
        unresolved_d_3 = unresolved_d.copy()
        if score_info_dict_1 and score_info_dict_2 and score_info_dict_3:
            sorted_dict_3 = get_score_sorted_dict(score_info_dict_3)
            sorted_dict_2 = get_score_sorted_dict(score_info_dict_2)
            sorted_dict_1 = get_score_sorted_dict(score_info_dict_1)
            sorted_dict_1, sorted_dict_2, sorted_dict_3 \
            = check_same_dbs_position_in_three(dict_1=sorted_dict_1, 
                                               dict_2=sorted_dict_2, 
                                               dict_3=sorted_dict_3)
            sorted_dict_3[len(sorted_dict_3)] = unresolved_d_3[0]
            sorted_dict_2[len(sorted_dict_2)] = unresolved_d_2[0]
            sorted_dict_1[len(sorted_dict_1)] = unresolved_d_1[0]
            revolved_level = 'All'
            validated_num = 3
            bool_list = [True, True, True]
        elif score_info_dict_2 and score_info_dict_3:
            sorted_dict_3 = get_score_sorted_dict(score_info_dict_3)
            sorted_dict_2 = get_score_sorted_dict(score_info_dict_2)
            sorted_dict_1 = unresolved_d_1
            sorted_dict_2, sorted_dict_3 = check_same_dbs_position_in_two(
                dict_1=sorted_dict_2, dict_2=sorted_dict_3, tag='3')
            sorted_dict_3[len(sorted_dict_3)] = unresolved_d_3[0]
            sorted_dict_2[len(sorted_dict_2)] = unresolved_d_2[0]
            revolved_level = 'Partial'
            validated_num = 2
            bool_list = [False, True, True]
        elif score_info_dict_1 and score_info_dict_3:
            sorted_dict_3 = get_score_sorted_dict(score_info_dict_3)
            sorted_dict_2 = unresolved_d_2
            sorted_dict_1 = get_score_sorted_dict(score_info_dict_1)
            sorted_dict_1, sorted_dict_3 = check_same_dbs_position_in_two(
                dict_1=sorted_dict_1, dict_2=sorted_dict_3, tag='3')
            sorted_dict_3[len(sorted_dict_3)] = unresolved_d_3[0]
            sorted_dict_1[len(sorted_dict_1)] = unresolved_d_1[0]
            bool_list = [True, False, True]
            revolved_level = 'Partial'
            validated_num = 2
        elif score_info_dict_1 and score_info_dict_2:
            sorted_dict_3 = unresolved_d_3
            sorted_dict_2 = get_score_sorted_dict(score_info_dict_2)
            sorted_dict_1 = get_score_sorted_dict(score_info_dict_1)
            sorted_dict_1, sorted_dict_2 = check_same_dbs_position_in_two(
                dict_1=sorted_dict_1, dict_2=sorted_dict_2, tag='2')
            sorted_dict_2[len(sorted_dict_2)] = unresolved_d_2[0]
            sorted_dict_1[len(sorted_dict_1)] = unresolved_d_1[0]
            bool_list = [True, True, False]
            revolved_level = 'Partial'
            validated_num = 2
        elif score_info_dict_3:
            sorted_dict_3 = get_score_sorted_dict(score_info_dict_3)
            sorted_dict_3[len(sorted_dict_3)] = unresolved_d_3[0]
            sorted_dict_2 = unresolved_d_2
            sorted_dict_1 = unresolved_d_1
            bool_list = [False, False, True]
            revolved_level = 'Partial'
            validated_num = 1
        elif score_info_dict_2:
            sorted_dict_3 = unresolved_d_3
            sorted_dict_2 = get_score_sorted_dict(score_info_dict_2)
            sorted_dict_2[len(sorted_dict_2)] = unresolved_d_2[0]
            sorted_dict_1 = unresolved_d_1
            bool_list = [False, True, False]
            revolved_level = 'Partial'
            validated_num = 1
        elif score_info_dict_1:
            sorted_dict_3 = unresolved_d_3
            sorted_dict_2 = unresolved_d_2
            sorted_dict_1 = get_score_sorted_dict(score_info_dict_1)
            sorted_dict_1[len(sorted_dict_1)] = unresolved_d_1[0]
            bool_list = [True, False, False]
            revolved_level = 'Partial'
            validated_num = 1
        else:
            sorted_dict_3 = unresolved_d_3
            sorted_dict_2 = unresolved_d_2
            sorted_dict_1 = unresolved_d_1
            bool_list = [False, False, False]
            revolved_level = 'None'
            validated_num = 0
        sorted_dict_1 = add_determined_db_info(sorted_dict_1)
        sorted_dict_2 = add_determined_db_info(sorted_dict_2)
        sorted_dict_3 = add_determined_db_info(sorted_dict_3)
        total_result_dict = {'Resolved level': revolved_level,
                             'Validated num': validated_num,
                             'Each bools': bool_list,
                             'Moiety-1': sorted_dict_1,
                             'Moiety-2': sorted_dict_2,
                             'Moiety-3': sorted_dict_3}
        return total_result_dict
        #endregion

def generate_reference_oad_dict(structure_dict, unsaturated_moieties_num, 
    msms_df, ms_tolerance, must_nl_cut_off_dict):
    """ Generate in-silico OAD fragment ions according to OAD fragments rules

    Args:
        structure_dict (dict): dict generated by extract_lipid_structural_info
        unsaturated_moieties_num (int): number of unsaturated moieties
        msms_df (Dataframe): Extracted MS/MS spectrum
        ms_tolerance (float): mass tolerance
        must_nl_cut_off_dict: 
            types of essential ions and relative intensity threshold

    Returns:
        ref_oad_dict (dict): {
            '1~3 (Moiety-No.)': {
                C=C position (tuple): {NL type tag: NL value}
            }
        }
    """
    moieties_info_length = len(structure_dict['Each moiety info'])
    db_in_sphingobase = structure_dict['Unsaturated sphingobase']
    ontology = structure_dict['Ontology']
    deuterium = structure_dict['Deuterium']
    db_combs_d = {}
    ref_oad_dict = {}
    if unsaturated_moieties_num == 1:
        if moieties_info_length == 2:
            chain_num = structure_dict['Each moiety info']['chain-1'] - 1
            db_num = structure_dict['Each moiety info']['db-1']
        elif moieties_info_length == 4:
            if structure_dict['Each moiety info']['db-1'] > 0:
                chain_num = structure_dict['Each moiety info']['chain-1'] - 1
                db_num = structure_dict['Each moiety info']['db-1']
            elif structure_dict['Each moiety info']['db-2'] > 0:
                chain_num = structure_dict['Each moiety info']['chain-2'] - 1
                db_num = structure_dict['Each moiety info']['db-2']
        elif moieties_info_length == 6:
            if structure_dict['Each moiety info']['db-1'] > 0:
                chain_num = structure_dict['Each moiety info']['chain-1'] - 1
                db_num = structure_dict['Each moiety info']['db-1']
            elif structure_dict['Each moiety info']['db-2'] > 0:
                chain_num = structure_dict['Each moiety info']['chain-2'] - 1
                db_num = structure_dict['Each moiety info']['db-2']
            elif structure_dict['Each moiety info']['db-3'] > 0:
                chain_num = structure_dict['Each moiety info']['chain-3'] - 1
                db_num = structure_dict['Each moiety info']['db-3']
        db_combs_d['1'] = calculate_db_pairs(
            chain_num, db_num, structure_dict, msms_df, ms_tolerance,
            db_in_sphingobase, must_nl_cut_off_dict, deuterium
        )
        ref_oad_dict['1'] = generate_ref_oad_nl_and_type(
            db_combs_d['1'], ontology, deuterium
        )
    elif unsaturated_moieties_num == 2:
        if moieties_info_length == 4:
            chain_num_1 = structure_dict['Each moiety info']['chain-1'] - 1
            db_num_1 = structure_dict['Each moiety info']['db-1']
            chain_num_2 = structure_dict['Each moiety info']['chain-2'] - 1
            db_num_2 = structure_dict['Each moiety info']['db-2']
            if ontology == 'SM': deuterium = 0
            if ontology == 'Cer_NS': deuterium = structure_dict['Deuterium']
            db_combs_d['1'] = calculate_db_pairs(chain_num_1, db_num_1, 
                structure_dict, msms_df, ms_tolerance, 
                db_in_sphingobase, must_nl_cut_off_dict, deuterium)
            db_in_sphingobase = False
            if ontology == 'SM': deuterium = structure_dict['Deuterium']
            if ontology == 'Cer_NS': deuterium = 0
            db_combs_d['2'] = calculate_db_pairs(chain_num_2, db_num_2, 
                structure_dict, msms_df, ms_tolerance,
                db_in_sphingobase, must_nl_cut_off_dict, deuterium)
        elif moieties_info_length == 6:
            if structure_dict['Each moiety info']['db-1'] > 0 and structure_dict['Each moiety info']['db-2'] > 0:
                chain_num_1 = structure_dict['Each moiety info']['chain-1'] - 1
                db_num_1 = structure_dict['Each moiety info']['db-1']
                chain_num_2 = structure_dict['Each moiety info']['chain-2'] - 1
                db_num_2 = structure_dict['Each moiety info']['db-2']
                db_combs_d['1'] = calculate_db_pairs(chain_num_1, db_num_1, 
                    structure_dict, msms_df, ms_tolerance,
                    db_in_sphingobase, must_nl_cut_off_dict, deuterium)
                db_in_sphingobase = False
                db_combs_d['2'] = calculate_db_pairs(chain_num_2, db_num_2, 
                    structure_dict, msms_df, ms_tolerance,
                    db_in_sphingobase, must_nl_cut_off_dict, deuterium)
            elif structure_dict['Each moiety info']['db-1'] > 0 and structure_dict['Each moiety info']['db-3'] > 0:
                chain_num_1 = structure_dict['Each moiety info']['chain-1'] - 1
                db_num_1 = structure_dict['Each moiety info']['db-1']
                chain_num_2 = structure_dict['Each moiety info']['chain-3'] - 1
                db_num_2 = structure_dict['Each moiety info']['db-3']
                db_combs_d['1'] = calculate_db_pairs(chain_num_1, db_num_1, 
                    structure_dict, msms_df, ms_tolerance,
                    db_in_sphingobase, must_nl_cut_off_dict, deuterium)
                db_in_sphingobase = False
                db_combs_d['2'] = calculate_db_pairs(chain_num_2, db_num_2, 
                    structure_dict, msms_df, ms_tolerance,
                    db_in_sphingobase, must_nl_cut_off_dict, deuterium)
            elif structure_dict['Each moiety info']['db-2'] > 0 and structure_dict['Each moiety info']['db-3'] > 0:
                chain_num_1 = structure_dict['Each moiety info']['chain-2'] - 1
                db_num_1 = structure_dict['Each moiety info']['db-2']
                chain_num_2 = structure_dict['Each moiety info']['chain-3'] - 1
                db_num_2 = structure_dict['Each moiety info']['db-3']
                db_combs_d['1'] = calculate_db_pairs(chain_num_1, db_num_1, 
                    structure_dict, msms_df, ms_tolerance,
                    db_in_sphingobase, must_nl_cut_off_dict, deuterium)
                db_combs_d['2'] = calculate_db_pairs(chain_num_2, db_num_2, 
                    structure_dict, msms_df, ms_tolerance,
                    db_in_sphingobase, must_nl_cut_off_dict, deuterium)
        if ontology == 'SM': deuterium = 0
        if ontology == 'Cer_NS': deuterium = structure_dict['Deuterium']
        ref_oad_dict['1'] = generate_ref_oad_nl_and_type(
            db_combs_d['1'], ontology, deuterium
        )
        if ontology == 'SM': deuterium = structure_dict['Deuterium']
        if ontology == 'Cer_NS': deuterium = 0
        ref_oad_dict['2'] = generate_ref_oad_nl_and_type(
            db_combs_d['2'], ontology, deuterium
        )
    elif unsaturated_moieties_num == 3:
        chain_num_1 = structure_dict['Each moiety info']['chain-1'] - 1
        db_num_1 = structure_dict['Each moiety info']['db-1']
        chain_num_2 = structure_dict['Each moiety info']['chain-2'] - 1
        db_num_2 = structure_dict['Each moiety info']['db-2']
        chain_num_3 = structure_dict['Each moiety info']['chain-3'] - 1
        db_num_3 = structure_dict['Each moiety info']['db-3']
        db_combs_d['1'] = calculate_db_pairs(chain_num_1, db_num_1, 
            structure_dict, msms_df, ms_tolerance,
            db_in_sphingobase, must_nl_cut_off_dict, deuterium)
        db_in_sphingobase = False
        db_combs_d['2'] = calculate_db_pairs(chain_num_2, db_num_2, 
            structure_dict, msms_df, ms_tolerance,
            db_in_sphingobase, must_nl_cut_off_dict, deuterium)
        db_combs_d['3'] = calculate_db_pairs(chain_num_3, db_num_3, 
            structure_dict, msms_df, ms_tolerance,
            db_in_sphingobase, must_nl_cut_off_dict, deuterium)
        ref_oad_dict['1'] = generate_ref_oad_nl_and_type(
            db_combs_d['1'], ontology, deuterium)
        ref_oad_dict['2'] = generate_ref_oad_nl_and_type(
            db_combs_d['2'], ontology, deuterium)
        ref_oad_dict['3'] = generate_ref_oad_nl_and_type(
            db_combs_d['3'], ontology, deuterium)
    return ref_oad_dict

#Diagnostic ions
def calculate_db_pairs(chain_num, db_num, structure_dict, msms_df, 
    ms_tolerance, db_in_sphingobase, must_nl_cut_off_dict, deuterium):
    """ Returns C=C positions theoretically possible

    Args:
        chain_num (int): number of carbons in a moiety
        db_num (int): number of double bonds in a moiety
        structure_dict (dict): dict generated by extract_lipid_structural_info
        msms_df (Dataframe): Extracted MS/MS spectrum
        ms_tolerance (float): mass tolerance
        db_in_sphingobase (bool): Whether there is unsaturation in sphingoid base
        must_nl_cut_off_dict: 
            types of essential ions and relative intensity threshold
        deuterium (int): number of deuteriums in a moiety

    Returns:
        refined_db_combs (list[tuple]): list of C=C positions candidates
    """
    delete_cand_h, delete_cand_t = [], []
    ontology = structure_dict['Ontology']
    # class_based_cutoff = get_class_specific_cutoff(ontology, ion=1)
    # cut_off_df = msms_df[msms_df['Ratio(%)'] >= class_based_cutoff]
    cut_off_df = msms_df[
        msms_df['Ratio(%)'] >= must_nl_cut_off_dict['diagnostic_1'][1]
    ]
    essential_ion_type = must_nl_cut_off_dict['diagnostic_1'][0]
    ref_mz = structure_dict['Ref precursor Mz']
    h2o_ms = ex_mass['H2O']
    dHs = deuterium*(ex_mass['D'] - ex_mass['H'])

    """ Generating essential OAD neutral loss """
    if db_in_sphingobase: # Return all possible double bond position list
        chain_range = range(3, chain_num-2)
        db_combs = list(itertools.combinations(chain_range, db_num))    
        i = 0
        remove_tuples_index_set = set()
        for v in db_combs:
            counter = 1
            while counter < db_num:
                if v[counter] - v[counter-1] == 1:
                    remove_tuples_index_set.add(i)
                    break
                counter += 1
            i += 1
        comb_list_len = len(db_combs)
        refined_db_combs = [
            db_combs[i] for i in range(0, comb_list_len) 
            if i not in remove_tuples_index_set
        ]
        return refined_db_combs
    else: # Generating ref OAD ions list
        essential_nl_dict_head = {}
        for position in range(3, chain_num): # 1st C=C
            mass_shift = 0
            ref_oad_ion = simulate_diagnostic_1_ions(
                position, essential_ion_type, mass_shift, dHs
            )
            essential_nl_dict_head[position] = ref_oad_ion
        if db_num >= 2: # multiple C=C
            essential_nl_dict_tail = {}
            mass_shift_counter = db_num - 1
            mass_shift = 2*ex_mass['H']*mass_shift_counter
            for position in range(3+db_num-1, chain_num):
                ref_oad_ion = simulate_diagnostic_1_ions(
                    position, essential_ion_type, mass_shift, dHs
                )
                essential_nl_dict_tail[position] = ref_oad_ion

    #region Selecting unnecessary candidates
    if ontology not in lipidclass_dict['Sphingolipids']:
        for pos, nl in essential_nl_dict_head.items(): # 1st C=C
            head_mz = ref_mz - nl - ms_tolerance
            tail_mz = ref_mz - nl + ms_tolerance
            extracted_df = cut_off_df[
                (cut_off_df['frag m/z'] >= head_mz)
               &(cut_off_df['frag m/z'] <= tail_mz)
            ]
            len_of_extracted_df = len(extracted_df)
            if len_of_extracted_df == 0:
                delete_cand_h.append(pos)
        if db_num >= 2: # multiple C=C
            for pos, nl in essential_nl_dict_tail.items():
                head_mz = ref_mz - nl - ms_tolerance
                tail_mz = ref_mz - nl + ms_tolerance
                extracted_df = cut_off_df[
                    (cut_off_df['frag m/z'] >= head_mz)
                   &(cut_off_df['frag m/z'] <= tail_mz)
                ]
                len_of_extracted_df = len(extracted_df)
                if len_of_extracted_df == 0:
                    delete_cand_t.append(pos)
    else:
        for pos, nl in essential_nl_dict_head.items(): # 1st C=C
            head_mz = ref_mz - nl - ms_tolerance
            tail_mz = ref_mz - nl + ms_tolerance
            head_H2O_mz = ref_mz - nl - ms_tolerance - h2o_ms
            tail_H2O_mz = ref_mz - nl + ms_tolerance - h2o_ms
            extracted_df = cut_off_df[
                (cut_off_df['frag m/z'] >= head_mz)
               &(cut_off_df['frag m/z'] <= tail_mz)
            ]
            extracted_H2O_df = cut_off_df[
                (cut_off_df['frag m/z'] >= head_H2O_mz)
               &(cut_off_df['frag m/z'] <= tail_H2O_mz)
            ]
            len_of_extracted_df = len(extracted_df)
            len_of_extracted_H2O_df = len(extracted_H2O_df)
            if (len_of_extracted_df == 0) and (len_of_extracted_H2O_df == 0):
                delete_cand_h.append(pos)
        if db_num >= 2: # multiple C=C
            for pos, nl in essential_nl_dict_tail.items():
                head_mz = ref_mz - nl - ms_tolerance
                tail_mz = ref_mz - nl + ms_tolerance
                head_H2O_mz = ref_mz - nl - ms_tolerance - h2o_ms
                tail_H2O_mz = ref_mz - nl + ms_tolerance - h2o_ms
                extracted_df = cut_off_df[
                    (cut_off_df['frag m/z'] >= head_mz)
                   &(cut_off_df['frag m/z'] <= tail_mz)
                ]
                extracted_H2O_df = cut_off_df[
                    (cut_off_df['frag m/z'] >= head_H2O_mz)
                   &(cut_off_df['frag m/z'] <= tail_H2O_mz)
                ]
                len_of_extracted_df = len(extracted_df)
                len_of_extracted_H2O_df = len(extracted_H2O_df)
                if (len_of_extracted_df == 0) and (len_of_extracted_H2O_df == 0):
                    delete_cand_t.append(pos)
    #endregion
    #region Return possible double bond position list
    chain_range = range(3, chain_num)
    db_combs = list(itertools.combinations(chain_range, db_num))    
    db_combs = [
        each_comb for each_comb in db_combs if each_comb[0] not in delete_cand_h 
        and each_comb[-1] not in delete_cand_t
    ]
    i = 0
    remove_tuples_index_set = set()
    for v in db_combs:
        counter = 1
        while counter < db_num:
            if v[counter] - v[counter-1] == 1:
                remove_tuples_index_set.add(i)
                break
            counter += 1
        i += 1
    comb_list_len = len(db_combs)
    refined_db_combs = [
        db_combs[i] for i in range(0, comb_list_len) 
        if i not in remove_tuples_index_set
    ]
    #endregion
    return refined_db_combs

#Diagnostic ions
def simulate_diagnostic_1_ions(pos, ion_type, ms_sht, dHs):
    """ Simulate neutral loss value of one essential ion

    Args:
        pos (int): n-terminal position of C=C
        ion_type (str): 'OAD02', 'OAD03', or 'OAD04'
        ms_sht (int): mass shift value
        dHs (float): m/z of deuterium

    Returns: neutral loss value (float) of essential ion
    """
    c_ms, h_ms, o_ms = ex_mass['C'], ex_mass['H'], ex_mass['O']
    if ion_type == 'OAD03':
        """ -C*(X-1)-2*H*(X-1)-H+O,     ex) -97.13811   in n-9  """
        pre_mnO_LH = c_ms*(pos-1)+2*h_ms*(pos-1)+h_ms-o_ms - ms_sht + dHs
        return pre_mnO_LH
    elif ion_type == 'OAD02':
        """ -C*(X-1)-2*H*(X-1)+O,       ex) -96.130285  in n-9  """
        pre_mnO = c_ms*(pos-1)+2*h_ms*(pos-1)-o_ms - ms_sht + dHs
        return pre_mnO
    elif ion_type == 'OAD04':
        """ -C*(X-1)-2*H*(X-1)-2H+O,    ex) -98.145935  in n-9  """
        pre_mnO_L2H = c_ms*(pos-1)+2*h_ms*(pos-1)+2*h_ms-o_ms - ms_sht + dHs
        return pre_mnO_L2H

#region Not used function
def get_class_specific_cutoff(ontology, ion):
    base_cutoff = 0.01
    if ontology == 'PC':
        return base_cutoff if ion == 1 else base_cutoff
    else:
        return base_cutoff
#endregion

#NL registration
def generate_ref_oad_nl_and_type(db_positions_comb_list, ontology, deuterium):
    """ Generates in-silico OAD fragment ions tagged with neutral loss type codes
    Args:
        db_positions_comb_list (list[tuple]): list of C=C positions
        ontology (str): Lipid subclass
        deuterium (int): Number of deuterium

    Returns: 
        ref_oad_dict (dict): {
            C=C position (tuple): {NL type tag: NL value}
        }
    """
    ref_oad_dict = {}
    c_ms, c13_ms, o_ms = ex_mass['C'], ex_mass['13C'], ex_mass['O']
    h_ms, d_ms = ex_mass['H'], ex_mass['D']
    h2o_ms = ex_mass['H2O']
    dHs = deuterium*(d_ms - h_ms)
    c13 = c13_ms - c_ms
    for each_comb in db_positions_comb_list:
        # peak_list = []
        # each_combs_d = {'Positions': '', 'Delta': ''}
        # delta_and_tag_dict = {}
        each_combs_d = {}
        mass_shift_counter = 0
        for pos in each_comb:
            ms_sht = 2*ex_mass['H']*mass_shift_counter
            #region Pre
            # """ -C*(X-1)-2*H*(X-1)+3*O,     ex) -64.140455  in n-9  """
            # pre_triO = c_ms*(pos-1)+2*h_ms*(pos-1)-3*o_ms-ms_sht+dHs
            # """ -C*(X-1)-2*H*(X-1)-H+3*O,   ex) -65.14828   in n-9  """
            # pre_triO_LH = c_ms*(pos-1)+2*h_ms*(pos-1)+h_ms-3*o_ms-ms_sht+dHs

            # """ -C*(X-3)-13C*2-2*H*(X-1)+H+O        ex) -95.13140   in n-9  """
            # pre_mnO_AH = c_ms*(pos-1)+2*h_ms*(pos-1)+h_ms-o_ms-ms_sht+dHs-2*c13
            # """ -C*(X-1)-13C-2*H*(X-1)+O,           ex) -96.134755  in n-9  """
            # pre_mnO = c_ms*(pos-1)+2*h_ms*(pos-1)+h_ms-o_ms-ms_sht+dHs-c13
            """ -C*(X-1)-2*H*(X-1)+H+O,             ex) -95.12246   in n-9  """
            pre_mnO_AH = c_ms*(pos-1)+2*h_ms*(pos-1)-h_ms-o_ms-ms_sht+dHs
            """ -C*(X-1)-2*H*(X-1)+O,               ex) -96.130285  in n-9  """
            pre_mnO = c_ms*(pos-1)+2*h_ms*(pos-1)-o_ms-ms_sht+dHs
            """ -C*(X-1)-2*H*(X-1)-H+O,             ex) -97.13811   in n-9  """
            pre_mnO_LH = c_ms*(pos-1)+2*h_ms*(pos-1)+h_ms-o_ms-ms_sht+dHs
            """ -C*(X-1)-2*H*(X-1)-2H+O,            ex) -98.145935  in n-9  """
            pre_mnO_L2H = c_ms*(pos-1)+2*h_ms*(pos-1)+2*h_ms-o_ms-ms_sht+dHs
            # """ -C*(X-3)-2*13C-2*H*(X-1)+H+O-H2O,   ex) -113.14196  in n-9  """
            # pre_mnO_AH_H2O = pre_mnO_AH + h2o_ms
            # """ -C*(X-1)-13C-2*H*(X-1)+O-H2O,       ex) -114.14532  in n-9  """
            # pre_mnO_H2O = pre_mnO + h2o_ms
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
            # """ -C*(X-1)-13C-2*H*(X-1),             ex) -124.12967  in n-9  """
            # exa = c_ms*(pos)+2*h_ms*(pos)-h_ms-ms_sht+dHs-c13
            """ -C*(X-1)-2*H*(X-1),                 ex) -124.12520  in n-9  """
            exa = c_ms*(pos)+2*h_ms*(pos-1)-ms_sht+dHs
            """ -C*(X)-2*H*(X)-H,                   ex) -125.13303  in n-9  """
            exa_LH = c_ms*(pos)+2*h_ms*(pos)-h_ms-ms_sht+dHs
            """ -C*(X)-2*H*(X),                     ex) -126.14085  in n-9  """
            exa_L2H = c_ms*(pos)+2*h_ms*(pos)-ms_sht+dHs
            #endregion
            #region Post
            # """ -C*(X+1)-2*H*(X)+H+2*O,     ex) -105.143195 in n-9  """
            # post_diO_AH = c_ms*(pos+1)+2*h_ms*(pos)-h_ms-2*o_ms-ms_sht+dHs
            # """ -C*(X+1)-2*H*(X)+2*O,       ex) -106.15102  in n-9  """
            # post_diO = c_ms*(pos+1)+2*h_ms*(pos)-2*o_ms-ms_sht+dHs
            # """ -C*(X+1)-2*H*(X)-H+2*O,     ex) -107.158845 in n-9  """
            # post_diO_LH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-2*o_ms-ms_sht+dHs

            # """ -C*(X+1)-13C-2*H*(X)-H+O,           ex) -123.15376  in n-9  """
            # post_mnO_LH = c_ms*(pos+1)+2*h_ms*(pos+1)-o_ms-ms_sht+dHs-c13
            """ -C*(X+1)-2*H*(X)-H+O,               ex) -123.15376  in n-9  """
            post_mnO_LH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-o_ms-ms_sht+dHs
            """ -C*(X+1)-2*H*(X)-2H+O,              ex) -124.16159  in n-9  """
            post_mnO_L2H = c_ms*(pos+1)+2*h_ms*(pos+1)-o_ms-ms_sht+dHs
            # """ -C*(X-1)-2*13C-2*H*(X)+H,           ex) -137.14196  in n-9  """
            # post_AH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-ms_sht+dHs-2*c13
            # """ -C*(X)-13C-2*H*(X),                 ex) -138.14532  in n-9  """
            # post = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-ms_sht+dHs-c13
            """ -C*(X+1)-2*H*(X)+H,                 ex) -137.133025 in n-9  """
            post_AH = c_ms*(pos+1)+2*h_ms*(pos)-h_ms-ms_sht+dHs
            """ -C*(X+1)-2*H*(X),                   ex) -138.14085  in n-9  """
            post = c_ms*(pos+1)+2*h_ms*(pos)-ms_sht+dHs
            """ -C*(X+1)-2*H*(X)-H,                 ex) -139.148675 in n-9  """
            post_LH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-ms_sht+dHs
            """ -C*(X+1)-2*H*(X)-2H,                ex) -140.1565   in n-9  """
            post_L2H = c_ms*(pos+1)+2*h_ms*(pos)+2*h_ms-ms_sht+dHs
            #endregion

            # each_combs_d[f'n-{pos}/dis@n-{pos-1}/+3O/OAD01'] = pre_triO
            # each_combs_d[f'n-{pos}/dis@n-{pos-1}/+3O-H/OAD02'] = pre_triO_LH
            # each_combs_d[f'n-{pos}/dis@n-{pos-1}/+O+H/'] = pre_mnO_AH
            each_combs_d[f'n-{pos}/dis@n-{pos-1}/+O+H/OAD01'] = pre_mnO_AH
            each_combs_d[f'n-{pos}/dis@n-{pos-1}/+O/OAD02'] = pre_mnO
            each_combs_d[f'n-{pos}/dis@n-{pos-1}/+O-H/OAD03'] = pre_mnO_LH
            each_combs_d[f'n-{pos}/dis@n-{pos-1}/+O-2H/OAD04'] = pre_mnO_L2H
            each_combs_d[f'n-{pos}/dis@n-{pos-1}/+O+H-H2O/OAD05'] = pre_mnO_AH_H2O
            each_combs_d[f'n-{pos}/dis@n-{pos-1}/+O-H2O/OAD06'] = pre_mnO_H2O
            each_combs_d[f'n-{pos}/dis@n-{pos-1}/+O-H-H2O/OAD07'] = pre_mnO_LH_H2O

            each_combs_d[f'n-{pos}/dis@n-{pos}/+O/OAD08'] = exa_mnO
            each_combs_d[f'n-{pos}/dis@n-{pos}/none/OAD09'] = exa
            each_combs_d[f'n-{pos}/dis@n-{pos}/-H/OAD10'] = exa_LH
            each_combs_d[f'n-{pos}/dis@n-{pos}/-2H/OAD11'] = exa_L2H
            # each_combs_d[f'n-{pos}/dis@n-{pos+1}/+2O+H/'] = post_diO_AH
            # each_combs_d[f'n-{pos}/dis@n-{pos+1}/+2O/OAD07'] = post_diO
            # each_combs_d[f'n-{pos}/dis@n-{pos+1}/+2O-H/OAD08'] = post_diO_LH
            # each_combs_d[f'n-{pos}/dis@n-{pos+1}/+H/'] = post_AH
            each_combs_d[f'n-{pos}/dis@n-{pos+1}/+O-H/OAD12'] = post_mnO_LH
            each_combs_d[f'n-{pos}/dis@n-{pos+1}/+O-2H/OAD13'] = post_mnO_L2H
            each_combs_d[f'n-{pos}/dis@n-{pos+1}/+H/OAD14'] = post_AH
            each_combs_d[f'n-{pos}/dis@n-{pos+1}/none/OAD15'] = post
            each_combs_d[f'n-{pos}/dis@n-{pos+1}/-H/OAD16'] = post_LH
            each_combs_d[f'n-{pos}/dis@n-{pos+1}/-2H/OAD17'] = post_L2H
            if ontology in lipidclass_dict['Sphingolipids']:
                each_combs_d[f'n-{pos}/dis@n-{pos+1}/-H2O/OAD18'] \
                = post + h2o_ms
                each_combs_d[f'n-{pos}/dis@n-{pos+1}/-H-H2O/OAD19'] \
                = post_LH + h2o_ms
                each_combs_d[f'n-{pos}/dis@n-{pos+1}/-2H-H2O/OAD20'] \
                = post_L2H + h2o_ms

            mass_shift_counter += 1

        ref_oad_dict[each_comb] = each_combs_d       
    return ref_oad_dict

#region Not used function
#NL registration
def simulate_oad_delta_with_tag(max_score_dict, ontology):
    db_positions_list = max_score_dict['Positions']
    mass_shift_counter = 0
    tag_d = {}
    c_ms, h_ms, o_ms = ex_mass['C'], ex_mass['H'], ex_mass['O']
    h2o_ms = ex_mass['H2O']
    for pos in db_positions_list:
        ms_sht = 2*ex_mass['H']*mass_shift_counter
        #region Pre
        # """ -C*(X-1)-2*H*(X-1)+3*O,     ex) -64.140455  in n-9  """
        # pre_triO = c_ms*(pos-1)+2*h_ms*(pos-1)-3*o_ms - ms_sht
        # """ -C*(X-1)-2*H*(X-1)-H+3*O,   ex) -65.14828   in n-9  """
        # pre_triO_LH = c_ms*(pos-1)+2*h_ms*(pos-1)+h_ms-3*o_ms - ms_sht
        """ -C*(X-1)-2*H*(X-1)+H+O,     ex) -95.12246   in n-9  """
        pre_mnO_AH = c_ms*(pos-1)+2*h_ms*(pos-1)-h_ms-o_ms - ms_sht
        """ -C*(X-1)-2*H*(X-1)+O,       ex) -96.130285  in n-9  """
        pre_mnO = c_ms*(pos-1)+2*h_ms*(pos-1)-o_ms - ms_sht
        """ -C*(X-1)-2*H*(X-1)-H+O,     ex) -97.13811   in n-9  """
        pre_mnO_LH = c_ms*(pos-1)+2*h_ms*(pos-1)+h_ms-o_ms - ms_sht
        """ -C*(X-1)-2*H*(X-1)-2H+O,    ex) -98.145935  in n-9  """
        pre_mnO_L2H = c_ms*(pos-1)+2*h_ms*(pos-1)+2*h_ms-o_ms - ms_sht
        """ -C*(X-1)-2*H*(X-1)+H+O-H2O, ex) -113.13303  in n-9  """
        pre_mnO_AH_H2O = pre_mnO_AH + h2o_ms
        """ -C*(X-1)-2*H*(X-1)+O-H2O,   ex) -114.14085  in n-9  """
        pre_mnO_H2O = pre_mnO + h2o_ms
        """ -C*(X-1)-2*H*(X-1)-H+O-H2O, ex) -115.14868  in n-9  """
        pre_mnO_LH_H2O = pre_mnO_LH + h2o_ms
        #endregion
        #region Exa
        """ -C*(X)-2*H*(X)+O,           ex) -110.145935 in n-9  """
        exa_mnO = c_ms*(pos)+2*h_ms*(pos)-o_ms - ms_sht
        """ -C*(X)-2*H*(X-1),           ex) -124.12520 in n-9  """
        exa = c_ms*(pos)+2*h_ms*(pos-1) - ms_sht
        """ -C*(X)-2*H*(X)-H,           ex) -125.13303 in n-9  """
        exa_LH = c_ms*(pos)+2*h_ms*(pos)-h_ms - ms_sht
        """ -C*(X)-2*H*(X),             ex) -126.14085 in n-9  """
        exa_L2H = c_ms*(pos)+2*h_ms*(pos) - ms_sht
        #endregion
        #region Post
        # """ -C*(X+1)-2*H*(X)+H+2*O,     ex) -105.143195 in n-9  """
        # post_diO_AH = c_ms*(pos+1)+2*h_ms*(pos)-h_ms-2*o_ms - ms_sht
        # """ -C*(X+1)-2*H*(X)+2*O,       ex) -106.15102  in n-9  """
        # post_diO = c_ms*(pos+1)+2*h_ms*(pos)-2*o_ms - ms_sht
        # """ -C*(X+1)-2*H*(X)-H+2*O,     ex) -107.158845 in n-9  """
        # post_diO_LH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-2*o_ms - ms_sht
        """ -C*(X+1)-2*H*(X)-H+O,       ex) -123.15376  in n-9  """
        post_mnO_LH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms-o_ms - ms_sht
        """ -C*(X+1)-2*H*(X)-2H+O,      ex) -124.16159  in n-9  """
        post_mnO_L2H = c_ms*(pos+1)+2*h_ms*(pos+1)-o_ms - ms_sht
        """ -C*(X+1)-2*H*(X)+H,         ex) -137.133025 in n-9  """
        post_AH = c_ms*(pos+1)+2*h_ms*(pos)-h_ms - ms_sht
        """ -C*(X+1)-2*H*(X),           ex) -138.14085  in n-9  """
        post = c_ms*(pos+1)+2*h_ms*(pos) - ms_sht
        """ -C*(X+1)-2*H*(X)-H,         ex) -139.148675 in n-9  """
        post_LH = c_ms*(pos+1)+2*h_ms*(pos)+h_ms - ms_sht
        """ -C*(X+1)-2*H*(X)-2H,        ex) -140.1565   in n-9  """
        post_L2H = c_ms*(pos+1)+2*h_ms*(pos)+2*h_ms - ms_sht
        #endregion

        # tag_d[f'n-{pos}/dis@n-{pos-1}/+3O/OAD01'] = pre_triO
        # tag_d[f'n-{pos}/dis@n-{pos-1}/+3O-H/OAD02'] = pre_triO_LH
        # tag_d[f'n-{pos}/dis@n-{pos-1}/+O+H/'] = pre_mnO_AH
        tag_d[f'n-{pos}/dis@n-{pos-1}/+O+H/OAD01'] = pre_mnO_AH
        tag_d[f'n-{pos}/dis@n-{pos-1}/+O/OAD02'] = pre_mnO
        tag_d[f'n-{pos}/dis@n-{pos-1}/+O-H/OAD03'] = pre_mnO_LH
        tag_d[f'n-{pos}/dis@n-{pos-1}/+O-2H/OAD04'] = pre_mnO_L2H
        tag_d[f'n-{pos}/dis@n-{pos-1}/+O+H-H2O/OAD05'] = pre_mnO_AH_H2O
        tag_d[f'n-{pos}/dis@n-{pos-1}/+O-H2O/OAD06'] = pre_mnO_H2O
        tag_d[f'n-{pos}/dis@n-{pos-1}/+O-H-H2O/OAD07'] = pre_mnO_LH_H2O

        tag_d[f'n-{pos}/dis@n-{pos}/+O/OAD08'] = exa_mnO
        tag_d[f'n-{pos}/dis@n-{pos}/none/OAD09'] = exa
        tag_d[f'n-{pos}/dis@n-{pos}/-H/OAD10'] = exa_LH
        tag_d[f'n-{pos}/dis@n-{pos}/-2H/OAD11'] = exa_L2H
        # tag_d[f'n-{pos}/dis@n-{pos+1}/+2O+H/'] = post_diO_AH
        # tag_d[f'n-{pos}/dis@n-{pos+1}/+2O/OAD07'] = post_diO
        # tag_d[f'n-{pos}/dis@n-{pos+1}/+2O-H/OAD08'] = post_diO_LH
        tag_d[f'n-{pos}/dis@n-{pos+1}/+O-H/OAD12'] = post_mnO_LH
        tag_d[f'n-{pos}/dis@n-{pos+1}/+O-2H/OAD13'] = post_mnO_L2H
        # tag_d[f'n-{pos}/dis@n-{pos+1}/+H/']  = post_AH
        tag_d[f'n-{pos}/dis@n-{pos+1}/+H/OAD14'] = post_AH
        tag_d[f'n-{pos}/dis@n-{pos+1}/none/OAD15'] = post
        tag_d[f'n-{pos}/dis@n-{pos+1}/-H/OAD16'] = post_LH
        tag_d[f'n-{pos}/dis@n-{pos+1}/-2H/OAD17'] = post_L2H
        if ontology in lipidclass_dict['Sphingolipids']:
            tag_d[f'n-{pos}/dis@n-{pos+1}/-H2O/OAD18'] \
            = post + h2o_ms
            tag_d[f'n-{pos}/dis@n-{pos+1}/-H-H2O/OAD19'] \
            = post_LH + h2o_ms
            tag_d[f'n-{pos}/dis@n-{pos+1}/-2H-H2O/OAD20'] \
            = post_L2H + h2o_ms   
        mass_shift_counter += 1
    return tag_d
#endregion

#region Not used function
#NL registration
def get_ref_oad_ratio(oad_type):
    if oad_type == 'OAD01': return 0.1
    elif oad_type == 'OAD02': return 0.25
    elif oad_type == 'OAD03': return 0.5
    elif oad_type == 'OAD04': return 0.01
    elif oad_type == 'OAD05': return 0.01
    elif oad_type == 'OAD06': return 0.05
    elif oad_type == 'OAD07': return 0.1
    elif oad_type == 'OAD08': return 0.06
    elif oad_type == 'OAD09': return 0.1
    elif oad_type == 'OAD10': return 0.2
    elif oad_type == 'OAD11': return 0.02
    elif oad_type == 'OAD12': return 0.02
    elif oad_type == 'OAD13': return 0.04
    elif oad_type == 'OAD14': return 0.05
    elif oad_type == 'OAD15': return 0.2
    elif oad_type == 'OAD16': return 0.4
    elif oad_type == 'OAD17': return 0.03
    elif oad_type == 'OAD18': return 0.1
    elif oad_type == 'OAD19': return 0.2
    elif oad_type == 'OAD20': return 0.01
#endregion

#NL registration&Diagnostic ions
def query_essential_diagnostic_ions(df, ref_oad_dict, db_in_sphingobase,
    c_num, db_num, for_range, tolerance, must_nl_cut_off_dict, structure_dict):
    """ Checking whether essential ions for each C=C positions were detected

    Args:
        df (Dataframe): Extracted MS/MS spectra
        ref_oad_dict (dict): {C=C position (tuple): {NL type tag: NL value}}
        db_in_sphingobase (bool): Whether there is unsaturation in sphingoid base
        c_num (int): number of carbons in a moiety
        db_num (int): number of double bonds in a moiety
        for_range: range(c_num)
        tolerance (float): mass tolerance
        must_nl_cut_off_dict: 
            types of essential ions and relative intensity threshold
        structure_dict (dict): dict generated by extract_lipid_structural_info

    Returns:
        diagnostic_nl_dict (dict): {
            C=C positions: whether resolved or not (list[bool])
        }
    """
    tol = tolerance
    frag = 'frag m/z'
    diagnostic_nl_dict = {}
    ontology = structure_dict['Ontology']
    ref_mz = structure_dict['Ref precursor Mz']
    diagnostic_1_type = must_nl_cut_off_dict['diagnostic_1'][0]
    diagnostic_2_type = must_nl_cut_off_dict['diagnostic_2'][0]
    diagnostic_1_cutoff = must_nl_cut_off_dict['diagnostic_1'][1]
    diagnostic_2_cutoff = must_nl_cut_off_dict['diagnostic_2'][1]
    # diagnostic_1_cutoff = get_class_specific_cutoff(ontology=ontology, ion=1)
    # diagnostic_2_cutoff = get_class_specific_cutoff(ontology=ontology, ion=2)
    post_tag1 = '+O-H/OAD03' if diagnostic_1_type == 'OAD03' else '+O/OAD02'
    post_tag2 = '-H/OAD16' if diagnostic_2_type == 'OAD16' else 'none/OAD15'
    cut_1_df = df[df['Ratio(%)'] >= diagnostic_1_cutoff]
    cut_2_df = df[df['Ratio(%)'] >= diagnostic_2_cutoff]
    sph_df = df[df['Ratio(%)'] >= must_nl_cut_off_dict['sphingobase']]
    if ontology not in lipidclass_dict['Sphingolipids']:
        for pos, tag_and_nl in ref_oad_dict.items():
            each_pos_bool = []
            tags_of_dgn_1 \
            = [f'n-{each}/dis@n-{each-1}/{post_tag1}' for each in pos]
            tags_of_dgn_2 \
            = [f'n-{each}/dis@n-{each+1}/{post_tag2}' for each in pos]
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
        post_tag3 = '+O-H-H2O/OAD07' if diagnostic_1_type == 'OAD03' else '+O-H2O/OAD06'
        post_tag4 = '-H-H2O/OAD19' if diagnostic_2_type == 'OAD16' else '-H2O/OAD18'
        if db_in_sphingobase:
            last_loop = db_num -1
            for pos, tag_and_nl in ref_oad_dict.items():
                each_pos_bool = []
                tags_of_dgn_1 \
                = [f'n-{each}/dis@n-{each-1}/{post_tag1}' for each in pos]
                tags_of_dgn_2 \
                = [f'n-{each}/dis@n-{each+1}/{post_tag2}' for each in pos]
                tags_of_dgn_3 \
                = [f'n-{each}/dis@n-{each-1}/{post_tag3}' for each in pos]
                tags_of_dgn_4 \
                = [f'n-{each}/dis@n-{each+1}/{post_tag4}' for each in pos]
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
                        h_post_LH_H2O_mz \
                        = ref_mz - tag_and_nl[post_LH_H2O_nl] - tol
                        t_post_LH_H2O_mz \
                        = ref_mz - tag_and_nl[post_LH_H2O_nl] + tol
                        h_post_L2H_H2O_mz \
                        = ref_mz - tag_and_nl[post_L2H_H2O_nl] - tol
                        t_post_L2H_H2O_mz \
                        = ref_mz - tag_and_nl[post_L2H_H2O_nl] + tol
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
                tags_of_dgn_1 = [f'n-{each}/dis@n-{each-1}/{post_tag1}' for each in pos]
                tags_of_dgn_2 = [f'n-{each}/dis@n-{each+1}/{post_tag2}' for each in pos]
                tags_of_dgn_3 = [f'n-{each}/dis@n-{each-1}/{post_tag3}' for each in pos]
                tags_of_dgn_4 = [f'n-{each}/dis@n-{each+1}/{post_tag4}' for each in pos]
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

#Diagnostic ions
def calc_presence_ratios_and_score(ref_oad_dict, diagnostic_ions_result_dict, 
    cut_df, ref_precursor_mz, ms_tolerance_ppm, sph_set):
    """ Construct data (dict) for ranking C=C positions candidates

    Args:
        ref_oad_dict (dict): {C=C position (tuple): {NL type tag: NL value}}
        diagnostic_ions_result_dict: {
            C=C positions: whether resolved or not (list[bool])
        }
        cut_df (Dataframe): Filtered (>= cut off ratio) OAD-MS/MS spectrum
        ref_precursor_mz (float): Reference m/z of precursor ion
        ms_tolerance_ppm (int): mass tolerance threshold based on ppm
        sph_set (list): [db_in_sphingobase, carbon number in sphingoid base]

    Returns:
        each_score_d (dict): dict contains belows
            'Positions', 'N-description', 'Score', 'Ratio sum', 'Presence',
            'Notice', 'Measured peaks', 'Ref peaks', 'Peaks dict'
    """
    counter = 0
    score_dict = {}
    tol = math_floor(ms_tolerance_ppm*ref_precursor_mz/(1000*1000), 6)
    def get_n_discription(txt):
        if ',)' in txt:
            edit = 'n-' + txt.replace('(', '').replace(',)', '')
        else:
            edit = 'n-' + txt.replace('(', '').replace(')', '').replace(' ', '')
        return edit
    for positions, tag_and_nl in ref_oad_dict.items():
        if all(diagnostic_ions_result_dict[positions]):
            presence_counter, ratio_sum = 0, 0
            next_to_3oh = sph_set[0] and ((sph_set[1]-positions[-1]) == 4)
            each_score_d = {
                'Positions': '', 'N-description': '',
                'Score': 0, 'Ratio sum': 0, 'Presence': 0,
                'Notice': '', 'Measured peaks': [],
                'Ref peaks': [], 'Peaks dict': {}
            }
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
                        mz, ratio, ppm = query_matched_ion_by_ppm(
                            cut_df, ref_mz, tol
                        )
                    else: continue
                else:
                    mz, ratio, ppm = query_matched_ion_by_ppm(
                        cut_df, ref_mz, tol
                    )
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

#Function which can implement isotope ion correction
def get_ref_ratio_via_db_position(positions, key):
    """ Returns reference ralative intensity of each OAD fragment ion

    Args:
        positions (tuple): C=C positions
        key: tag of OAD NL type

    Returns:
        ref_ratio (float): heulistric relative intensity via NL types
    """
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
    """ Query OAD fragment ions within given mass tolerance (ppm)

    Args:
        cut_df (Dataframe): Filtered (>= cut off ratio) OAD-MS/MS spectrum
        ref_mz (float): Reference m/z of fragment ion
        tol (float): mass tolerance

    Returns:
        mz (float): m/z of detected ion
        ratio (float): Relative intensity of detected ion
        ppm (float): Calculated mass difference
    """
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
    """ Construct relative intensity vector (list[float]) of 
        measured and reference ions

    Args:
        peaks_dict (dict): {
            tag of NL types: [
                Ref m/z, Ref delta, Measured m/z, Measured ratio, ppm
            ]
        }
        ref_peaks (list[list]): [[OAD type, Ref m/z, Ref NL, Ref ratio], ...]

    Returns:
        acts (list[float]): relative intensity vector of measured ions
        refs (list[float]): relative intensity vector of reference ions
    """
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

#region Not used function
def constrcut_ranking_vec(peaks_dict, ref_peaks):
    # peaks_dict = {'n-9/dis@n-8/+O/OAD03': 
    #               [Ref m/z, Ref delta, Measured m/z, Measured ratio, ppm]}
    # ref_peaks = [[OAD type, Ref m/z, Ref NL, Ref ratio], ...]
    oad_types, ref_mzs, ref_ratio, mzs, exp_ratio = [], [], [], [], []
    for key, v in peaks_dict.items():
        oad_types.append(key), ref_mzs.append(v[0])
        mzs.append(v[2]), exp_ratio.append(v[3])
        for i, li in enumerate(ref_peaks):
            if key == li[0]: break
        ref_ratio.append(ref_peaks[i][3])
    df = pd.DataFrame({'OAD type': oad_types, 'Ref m/z': ref_mzs, 
                       'Ref ratio': ref_ratio,'Exp m/z': mzs,
                       'Exp ratio': exp_ratio})
    rank_method = 'min'
    r_df = df.rank(method=rank_method)
    rank_list = list(set(r_df['Ref ratio'].values))
    min_rank, max_rank = rank_list[0], rank_list[-1]
    ref_ranks, exp_ranks = r_df['Ref ratio'].values, r_df['Exp ratio'].values
    reconst_exp_ranks = []
    #region Average or Maximum
    if rank_method == 'average' or rank_method == 'max':
        for ref, exp in zip(ref_ranks, exp_ranks):
            rerank = 0.0
            if ref == min_rank and exp <= ref:
                rerank = exp
            else:
                idx = rank_list.index(ref)
                if rank_list[idx-1] < exp <= rank_list[idx]:
                    rerank = exp
            reconst_exp_ranks.append(rerank)
    #endregion
    #region Minimum
    if rank_method == 'min':
        for ref, exp in zip(ref_ranks, exp_ranks):
            rerank = 0.0
            idx = rank_list.index(ref)
            if ref == max_rank:
                if rank_list[idx-1] < exp <= rank_list[idx]:
                    rerank = exp
            else:
                if rank_list[idx] <= exp < rank_list[idx+1]:
                    rerank = exp
            reconst_exp_ranks.append(rerank)
    #endregion
    return reconst_exp_ranks, ref_ranks
#endregion

#region Not used function
def reconstruct_db_positions(score_dict):
    no1_pos = score_dict[0]['Positions']
    each_pos_cands = {i+1: {} for i in range(len(no1_pos))}
    # each_pos_cands = {1: {6: 0.3, 7: 0.1},
    #                   2: {9, 0.4, 10: 0.04},
    #                   3: {12, 0.2, 13: 0.09},
    #                   4: {15: 0.3, 16: 0.03}}
    for num, each_score_d in score_dict.items():
        dbs = each_score_d['Positions']
        peaks = each_score_d['Peaks dict']
        for j, pos in enumerate(dbs, start=1):
            oad3 = f'n-{pos}/dis@n-{pos-1}/+O-H/OAD03'
            oad7 = f'n-{pos}/dis@n-{pos-1}/+O-H-H2O/OAD07'
            oad15 = f'n-{pos}/dis@n-{pos+1}/none/OAD15'
            oad16 = f'n-{pos}/dis@n-{pos+1}/-H/OAD16'
            oad17 = f'n-{pos}/dis@n-{pos+1}/-2H/OAD17'
            oad18 = f'n-{pos}/dis@n-{pos+1}/-H2O/OAD18'
            oad19 = f'n-{pos}/dis@n-{pos+1}/-H-H2O/OAD19'
            oad20 = f'n-{pos}/dis@n-{pos+1}/-2H-H2O/OAD20'
            # if (oad3 in peaks) and (oad16 in peaks):
            #     each_pos_cands[j][pos] = peaks[oad3][3] + peaks[oad16][3]
            # elif (oad7 in peaks) and (oad19 in peaks):
            #     each_pos_cands[j][pos] = peaks[oad7][3] + peaks[oad19][3]
            # else:
            #     sph_bool1 = [k in peaks for k in [oad15, oad16, oad17]]
            #     sph_bool2 = [k in peaks for k in [oad18, oad19, oad20]]
            #     if (sph_bool1.count(True) >= 2) and 
            head = f'n-{pos}/'
            each_pos_cands[j][pos] = sum(
                [v[3] for k, v in peaks.items() if k.startswith(head)])
    intense_pos = []
    for pos_and_ratio in each_pos_cands.values():
        sorted_pos = sorted(pos_and_ratio.items(), key=lambda x:x[1])
        intense_pos.append(sorted_pos[0][0])
    new_no1_pos = tuple(intense_pos)
    try:
        keys = [k for k, v in score_dict.items() if v['Positions'] == new_no1_pos]
        key = keys[0]
    except:
        for v in score_dict.values():
            print(f"Positions: {v['Positions']}")
        print(f'New No1. :{new_no1_pos}')
        # print(keys)
        raise KeyError
    score_dict[key]['Score'] += 999
    return score_dict
#endregion


def score_thresholding(score_info_dict, score_cutoff):
    """ Filtering candidates by given score

    Args:
        score_info_dict (dict): {
            'Positions', 'N-description', 'Score', 'Ratio sum', 'Presence',
            'Notice', 'Measured peaks', 'Ref peaks', 'Peaks dict'
        }
        score_cutoff (float): threshold of ms/ms similarity score

    Returns:
        cutoff_dict (dict): score_info_dict thresholded by score cut off
    """
    cutoff_dict = {rank: v for rank, v in score_info_dict.items() 
                   if v['Score'] >= score_cutoff}
    cutoff_dict = {i: v for i, v in enumerate(cutoff_dict.values())}
    return cutoff_dict

def get_presence_and_ratio_sorted_dict(score_info_dict):
    """ Ranking candidates by presence and sum of OAD fragment ions

    Args:
        score_info_dict (dict): {
            'Positions', 'N-description', 'Score', 'Ratio sum', 'Presence',
            'Notice', 'Measured peaks', 'Ref peaks', 'Peaks dict'
        }

    Returns:
        sorted_dict (dict): score_info_dict sorted 
                            by presence and sum of OAD fragment ions
    """
    column_list = ['Index', 'Score', 'Presence', 'Ratio sum']
    rows = range(len(score_info_dict))
    scoring_df = pd.DataFrame(columns=column_list, index=rows)
    for idx, v in score_info_dict.items():
        scoring_df.loc[idx:idx, column_list] = [idx, v['Score'], 
                                                v['Presence'], v['Ratio sum']]
    sorted_df = scoring_df.sort_values(['Presence', 'Ratio sum'], 
                                        ascending=[False, False])
    sorted_index = sorted_df['Index'].to_list()
    sorted_dict = {n:score_info_dict[idx] for n, idx in enumerate(sorted_index)}
    return sorted_dict

def get_score_sorted_dict(score_info_dict):
    """ Ranking candidates by MS/MS similarity score, 
        and presence and sum of OAD fragment ions

    Args:
        score_info_dict (dict): {
            'Positions', 'N-description', 'Score', 'Ratio sum', 'Presence',
            'Notice', 'Measured peaks', 'Ref peaks', 'Peaks dict'
        }

    Returns:
        sorted_dict (dict): score_info_dict sorted by score
    """
    column_list = ['Index', 'Score', 'Presence', 'Ratio sum']
    rows = range(len(score_info_dict))
    scoring_df = pd.DataFrame(columns=column_list, index=rows)
    for idx, v in score_info_dict.items():
        scoring_df.loc[idx:idx, column_list] = [idx, v['Score'], 
                                                v['Presence'], v['Ratio sum']]
    sorted_df = scoring_df.sort_values(['Score', 'Presence', 'Ratio sum'], 
                                        ascending=[False, False, False])
    sorted_index = sorted_df['Index'].to_list()
    sorted_dict = {n:score_info_dict[idx] for n, idx in enumerate(sorted_index)}
    return sorted_dict

#region Not used function
def check_same_db_position_in_two(dict_1, dict_2, tag):
    first_dbs_in_moiety_2 = [v['Positions'][0] for v in dict_2.values()]
    first_dbs_in_moiety_1 = [v['Positions'][0] for v in dict_1.values()]
    first_db_in_moiety_2 = first_dbs_in_moiety_2[0]
    first_db_in_moiety_1 = first_dbs_in_moiety_1[0]
    if first_db_in_moiety_1 != first_db_in_moiety_2:
        return dict_1, dict_2
    else:
        equal_pos_idx_list = []
        unequal_pos_idx_list = []
        refined_dict_1 = {}
        end = len(dict_1)
        for key, pos in enumerate(first_dbs_in_moiety_1):
            if pos == first_db_in_moiety_2:
                equal_pos_idx_list.append(key)
            else:
                unequal_pos_idx_list.append(key)
        for idx in equal_pos_idx_list:
            dict_1[idx]['Notice'] = 'Same as Moiety-{}'.format(tag)
        for new_idx, old_idx in enumerate(unequal_pos_idx_list):
            refined_dict_1[new_idx] = dict_1[old_idx]
        start = len(unequal_pos_idx_list)
        for new_idx, old_idx in zip(range(start, end), equal_pos_idx_list):
            refined_dict_1[new_idx] = dict_1[old_idx]
        return refined_dict_1, dict_2
#endregion

def check_same_dbs_position_in_two(dict_1, dict_2, tag):
    """ Checking 2 candidates between moiety-1 and moiety-2
        and rearrange C=C positions to escape duplicative position

    Args:
        dict_1 and dict_2 (dict): {
            'Positions', 'N-description', 'Score', 'Ratio sum', 'Presence',
            'Notice', 'Measured peaks', 'Ref peaks', 'Peaks dict'
        }
        tag (str): target moiety number ('1', '2', '3')

    Returns:
        dict_1 and dict_2 (dict):
            No.1 candidate of C=C positions will be re-ranked 
            if those of another moiety were same.
    """
    first_db_in_moiety_2 = next(iter(dict_2.values()))['Positions']
    first_db_in_moiety_1 = next(iter(dict_1.values()))['Positions']
    if len(first_db_in_moiety_2) >= len(first_db_in_moiety_1):
        dup_db = tuple(first_db_in_moiety_2[i] 
                    for i in range(len(first_db_in_moiety_1)))
        if first_db_in_moiety_1 != dup_db:
            return dict_1, dict_2
        else:
            equal_pos_idx_list, unequal_pos_idx_list = [], []
            refined_dict_1 = {}
            end = len(dict_1)
            dbs_in_moiety_1 = [v['Positions'] for v in dict_1.values()]
            for key, dbs in enumerate(dbs_in_moiety_1):
                dup_bools = []
                for pos, dup in zip(dbs, dup_db):
                    if pos == dup: dup_bools.append(True)
                    else: dup_bools.append(False)
                if any(dup_bools): equal_pos_idx_list.append(key)
                else: unequal_pos_idx_list.append(key)
            if len(unequal_pos_idx_list) == 0:
                equal_pos_idx_list, unequal_pos_idx_list = [], []
                for key, dbs in enumerate(dbs_in_moiety_1):
                    if dbs[0] == dup_db[0]: equal_pos_idx_list.append(key)
                    else: unequal_pos_idx_list.append(key)
            for idx in equal_pos_idx_list:
                dict_1[idx]['Notice'] = f'Same as Moiety-{tag}'
            for new_idx, old_idx in enumerate(unequal_pos_idx_list):
                refined_dict_1[new_idx] = dict_1[old_idx]
            start = len(unequal_pos_idx_list)
            for new_idx, old_idx in zip(range(start, end), equal_pos_idx_list):
                refined_dict_1[new_idx] = dict_1[old_idx]
        return refined_dict_1, dict_2
    else:
        dup_db = tuple(first_db_in_moiety_1[i] 
                    for i in range(len(first_db_in_moiety_2)))
        if first_db_in_moiety_2 != dup_db:
            return dict_1, dict_2
        else:
            equal_pos_idx_list, unequal_pos_idx_list = [], []
            refined_dict_2 = {}
            end = len(dict_2)
            dbs_in_moiety_2 = [v['Positions'] for v in dict_2.values()]
            for key, dbs in enumerate(dbs_in_moiety_2):
                dup_bools = []
                for pos, dup in zip(dbs, dup_db):
                    if pos == dup: dup_bools.append(True)
                    else: dup_bools.append(False)
                if any(dup_bools): equal_pos_idx_list.append(key)
                else: unequal_pos_idx_list.append(key)
            if len(unequal_pos_idx_list) == 0:
                equal_pos_idx_list, unequal_pos_idx_list = [], []
                for key, dbs in enumerate(dbs_in_moiety_2):
                    if dbs[0] == dup_db[0]: equal_pos_idx_list.append(key)
                    else: unequal_pos_idx_list.append(key)
            for idx in equal_pos_idx_list:
                dict_2[idx]['Notice'] = f'Same as Moiety-{tag}'
            for new_idx, old_idx in enumerate(unequal_pos_idx_list):
                refined_dict_2[new_idx] = dict_2[old_idx]
            start = len(unequal_pos_idx_list)
            for new_idx, old_idx in zip(range(start, end), equal_pos_idx_list):
                refined_dict_2[new_idx] = dict_2[old_idx]
        return dict_1, refined_dict_2

#region Not used function
def check_same_db_position_in_three(dict_1, dict_2, dict_3):
    def no_duplicates(li):
        return len(li) == len(set(li))
    first_db_list = [dict_1[0]['Positions'][0], 
                     dict_2[0]['Positions'][0], 
                     dict_3[0]['Positions'][0]]
    if no_duplicates(first_db_list):
        return dict_1, dict_2, dict_3
    else:
        refined_dict_2, dict_3 = check_same_db_position_in_two(
                                dict_1=dict_2, dict_2=dict_3, tag='3')
        first_db_in_moiety_2 = refined_dict_2[0]['Positions'][0]
        first_db_in_moiety_3 = dict_3[0]['Positions'][0]
        if first_db_in_moiety_2 == first_db_in_moiety_3:
            refined_dict_1, refined_dict_2 = check_same_db_position_in_two(
                                dict_1=dict_1, dict_2=refined_dict_2, tag='2&3')
        else:
            first_dbs_in_moiety_1 = [v['Positions'][0] for v in dict_1.values()]
            check_list = []
            refined_dict_1 = {}
            for pos in first_dbs_in_moiety_1:
                li = [pos, first_db_in_moiety_2, first_db_in_moiety_3]
                check_list.append(no_duplicates(li))
            if any(check_list):
                unequal_list = [i for i, t_or_f in enumerate(check_list) if t_or_f]
                equal_list = [i for i, t_or_f in enumerate(check_list) if not t_or_f]
                end = len(first_dbs_in_moiety_1)
                for new_idx, old_idx in enumerate(unequal_list):
                    refined_dict_1[new_idx] = dict_1[old_idx]
                start = len(refined_dict_1)
                for new_idx, old_idx in zip(range(start, end), equal_list):
                    first_db = dict_1[old_idx]['Positions'][0]
                    if first_db == first_db_in_moiety_3:
                        dict_1[old_idx]['Notice'] = 'Same as Moiety-3'
                    else:
                        dict_1[old_idx]['Notice'] = 'Same as Moiety-2'
                    refined_dict_1[new_idx] = dict_1[old_idx]
            else:
                for idx, pos in enumerate(first_dbs_in_moiety_1):
                    if pos == first_db_in_moiety_3:
                        dict_1[idx]['Notice'] = 'Same as Moiety-3'
                    else:
                        dict_1[idx]['Notice'] = 'Same as Moiety-2'
                refined_dict_1 = dict_1
        return refined_dict_1, refined_dict_2, dict_3
#endregion

def check_same_dbs_position_in_three(dict_1, dict_2, dict_3):
    """ Checking 3 candidates amoung moiety-1, moiety-2 and moiety-3
        and rearrange C=C positions to escape duplicative position

    Args:
        dict_1, dict_2 and dict_3 (dict): {
            'Positions', 'N-description', 'Score', 'Ratio sum', 'Presence',
            'Notice', 'Measured peaks', 'Ref peaks', 'Peaks dict'
        }

    Returns:
        dict_1, dict_2 and dict_3 (dict):
            No.1 candidate of C=C positions will be re-ranked 
            if those of another moiety were same.
    """
    def no_duplicates(li):
        return len(li) == len(set(li))
    refined_dict_2, dict_3 = check_same_dbs_position_in_two(
                            dict_1=dict_2, dict_2=dict_3, tag='3')
    first_dbs_in_moiety_2 = refined_dict_2[0]['Positions']
    first_dbs_in_moiety_3 = dict_3[0]['Positions']
    determined_db_3 = tuple(first_dbs_in_moiety_3[i] 
                            for i in range(len(first_dbs_in_moiety_2)))
    if first_dbs_in_moiety_2 == determined_db_3:
        refined_dict_1, refined_dict_2 = check_same_dbs_position_in_two(
                            dict_1=dict_1, dict_2=refined_dict_2, tag='2&3')
    else:
        dbs_in_moiety_1 = [v['Positions'] for v in dict_1.values()]
        # first_dbs_in_moiety_1 = dict_1[0]['Positions']
        determined_db_2 = tuple(first_dbs_in_moiety_2[i] 
                     for i in range(len(dbs_in_moiety_1[0])))
        determined_db_3 = tuple(first_dbs_in_moiety_3[i] 
                     for i in range(len(dbs_in_moiety_1[0])))
        check_list = []
        refined_dict_1 = {}
        for pos in dbs_in_moiety_1:
            li = [pos, determined_db_2, determined_db_3]
            check_list.append(no_duplicates(li))
        if any(check_list):
            unequal_list = [i for i, t_or_f in enumerate(check_list) if t_or_f]
            equal_list = [i for i, t_or_f in enumerate(check_list) if not t_or_f]
            end = len(dbs_in_moiety_1)
            for new_idx, old_idx in enumerate(unequal_list):
                refined_dict_1[new_idx] = dict_1[old_idx]
            start = len(refined_dict_1)
            for new_idx, old_idx in zip(range(start, end), equal_list):
                db = dict_1[old_idx]['Positions']
                if db == determined_db_3:
                    dict_1[old_idx]['Notice'] = 'Same as Moiety-3'
                else:
                    dict_1[old_idx]['Notice'] = 'Same as Moiety-2'
                refined_dict_1[new_idx] = dict_1[old_idx]
        else:
            for idx, pos in enumerate(dbs_in_moiety_1):
                if pos == determined_db_3:
                    dict_1[idx]['Notice'] = 'Same as Moiety-3'
                else:
                    dict_1[idx]['Notice'] = 'Same as Moiety-2'
            refined_dict_1 = dict_1
    return refined_dict_1, refined_dict_2, dict_3

def add_determined_db_info(oad_dict):
    """ Add determined C=C positions to 'Determined db' key
    """
    if not oad_dict:
        return oad_dict
    # ref_peaks = []
    # measured_peaks = []
    # for key, v in oad_dict[0]['Peaks dict'].items(): #[OAD type, Ref m/z, Ref NL, Ref ratio]
    #     oad_type = key.split('/')[-1]
    #     ref_peaks.append([key, v[0], v[1], get_ref_oad_ratio(oad_type)])
    #     if v[2] > 0: #[Measured m/z, Measured ratio, ppm]
    #         measured_peaks.append([v[2], v[3], v[4]])
    # oad_dict['Determined db'] = {'Positions': oad_dict[0]['Positions'],
    #                              'N-description': oad_dict[0]['N-description'],
    #                              'Score': oad_dict[0]['Score'],
    #                              'Ratio sum': oad_dict[0]['Ratio sum'],
    #                              'Presence': oad_dict[0]['Presence'],
    #                              'Notice': oad_dict[0]['Notice'],
    #                              'Peaks dict': oad_dict[0]['Peaks dict'],
    #                              'Measured peaks': measured_peaks,
    #                              'Ref peaks': ref_peaks}
    oad_dict['Determined db'] = oad_dict[0]
    return oad_dict

def determine_oad_metabolite_name_N_description(oad_result_dict, structure_dict, 
    metabolite_name):
    """ Adds C=C positions to metabolite name as n-terminal desctription

    Args:
        oad_result_dict (dict): {
            Resolved level (str): 'All', 'Partial', 'None'
            Validated num (int): 0-3
            Each bools (list[bool]): Whether each moiety was resolved
            Moiety-1~3 (dict): Contains 'Score', 'Presence', and 'Ratio sum'
        }
        structure_dict (dict): {
            'Status', 'Adduct', 'Precursor Mz', 'MS2 Mz', 'Precise precursor Mz',
            'Ref precursor Mz', 'RT(min)', 'Reference RT', 'Ontology', 'Brutto', 
            'Valid moiety num', 'Each moiety info', 'Unsaturated moiety', 
            'Unsaturated sphingobase', 'Deuterium', 'SMILES', 'Formula', 
            'Atom dict', 'Oxidized', 'NL type'
        }
        metabolite_name (str): inputed metabolite name

    Returns:
        oad_name (str):
            Metabolite name whose C=C positions are added and described 
            by n-terminal method
    """
    moieties = structure_dict['Each moiety info']
    unsaturated_moieties_num = structure_dict['Unsaturated moiety']
    if 'Plasm' in structure_dict['Ontology']:
        target_str = 'O-' + str(moieties['chain-1']) + ':' + str(moieties['db-1']+1)
        replace_str = 'P-' + str(moieties['chain-1']) + ':' + str(moieties['db-1'])
        metabolite_name = metabolite_name.replace(target_str, replace_str)
    def refine_db_position_str(raw_txt):
        refined_str = raw_txt
        if '[' in raw_txt:
            refined_str = raw_txt.replace('[', '').replace(']', '').replace(' ', '')
        elif '(' in raw_txt:
            refined_str = raw_txt.replace('(', '').replace(',)', '').replace(')', '')
        refined_str = refined_str.replace(' ', '')
        return refined_str
    if unsaturated_moieties_num == 1:
        db_positions = str(
            oad_result_dict['Moiety-1']['Determined db']['Positions'])
        db_positions = refine_db_position_str(db_positions)
        if len(moieties) == 2:
            target_str = str(moieties['chain-1']) + ':' + str(moieties['db-1'])
        elif len(moieties) == 4:
            if moieties['db-1'] > 0:
                target_str = str(moieties['chain-1']) + ':' + str(moieties['db-1'])
            else:
                target_str = str(moieties['chain-2']) + ':' + str(moieties['db-2'])
        elif len(moieties) == 6:
            if moieties['db-1'] > 0:
                target_str = str(moieties['chain-1']) + ':' + str(moieties['db-1'])
            elif moieties['db-2'] > 0:
                target_str = str(moieties['chain-2']) + ':' + str(moieties['db-2'])
            elif moieties['db-3'] > 0:
                target_str = str(moieties['chain-3']) + ':' + str(moieties['db-3'])
        db_sol_moiety = target_str + '(n-' + db_positions + ')'
        oad_name = metabolite_name.replace(target_str, db_sol_moiety)
    elif unsaturated_moieties_num == 2:
        try:
            db_positions_2 = str(oad_result_dict['Moiety-2']['Determined db']['Positions'])
        except KeyError:
            db_positions_2 = ''
        try:
            db_positions_1 = str(oad_result_dict['Moiety-1']['Determined db']['Positions'])
        except KeyError:
            db_positions_1 = ''
        db_positions_2 = refine_db_position_str(db_positions_2)
        db_positions_1 = refine_db_position_str(db_positions_1)
        if len(moieties) == 4:
            target_str_1 = str(moieties['chain-1']) + ':' + str(moieties['db-1'])
            target_str_2 = str(moieties['chain-2']) + ':' + str(moieties['db-2'])
        elif len(moieties) == 6:
            if moieties['db-1'] > 0 and moieties['db-2'] > 0:
                target_str_1 = str(moieties['chain-1']) + ':' + str(moieties['db-1'])
                target_str_2 = str(moieties['chain-2']) + ':' + str(moieties['db-2'])
            elif moieties['db-1'] > 0 and moieties['db-3'] > 0:
                target_str_1 = str(moieties['chain-1']) + ':' + str(moieties['db-1'])
                target_str_2 = str(moieties['chain-3']) + ':' + str(moieties['db-3'])
            elif moieties['db-2'] > 0 and moieties['db-3'] > 0:
                target_str_1 = str(moieties['chain-2']) + ':' + str(moieties['db-2'])
                target_str_2 = str(moieties['chain-3']) + ':' + str(moieties['db-3'])
        db_sol_1 = target_str_1 + '(n-' + db_positions_1 + ')'
        db_sol_2 = target_str_2 + '(n-' + db_positions_2 + ')'
        rep_name = metabolite_name.replace(target_str_1, 'target1', 1)
        rep_name = rep_name.replace(target_str_2, 'target2', 1)
        oad_name = rep_name.replace('target1', db_sol_1).replace('target2', db_sol_2)
    elif unsaturated_moieties_num == 3:
        try:
            db_positions_3 = str(oad_result_dict['Moiety-3']['Determined db']['Positions'])
        except KeyError:
            db_positions_3 = ''
        try:
            db_positions_2 = str(oad_result_dict['Moiety-2']['Determined db']['Positions'])
        except KeyError:
            db_positions_2 = ''
        try:
            db_positions_1 = str(oad_result_dict['Moiety-1']['Determined db']['Positions'])
        except KeyError:
            db_positions_1 = ''
        db_positions_3 = refine_db_position_str(db_positions_3)
        db_positions_2 = refine_db_position_str(db_positions_2)
        db_positions_1 = refine_db_position_str(db_positions_1)
        target_str_1 = str(moieties['chain-1']) + ':' + str(moieties['db-1'])
        target_str_2 = str(moieties['chain-2']) + ':' + str(moieties['db-2'])
        target_str_3 = str(moieties['chain-3']) + ':' + str(moieties['db-3'])
        db_sol_1 = target_str_1 + '(n-' + db_positions_1 + ')'
        db_sol_2 = target_str_2 + '(n-' + db_positions_2 + ')'
        db_sol_3 = target_str_3 + '(n-' + db_positions_3 + ')'
        rep_name = metabolite_name.replace(target_str_1, 'target1', 1)
        rep_name = rep_name.replace(target_str_2, 'target2', 1)
        rep_name = rep_name.replace(target_str_3, 'target3', 1)
        oad_name = rep_name.replace('target1', db_sol_1).replace('target2', db_sol_2).replace('target3', db_sol_3)
    oad_name = oad_name.replace('(n-)', '')
    return oad_name

#region Not used function
def get_oad_similarity_score(act_presence_rate, act_ratio_sum):
    # ref_ratio_sum = sum(list(ref_oad_nl_ratio_dict.values()))
    oad_similarity_score = round(act_presence_rate*act_ratio_sum, 1)
    return oad_similarity_score
#endregion

#region Not used function
def get_rel_ratio(base_ratio, oad_type):
    ratio = math_floor(base_ratio*rel_oad_ratio[oad_type], 2)
    return ratio
#endregion

def get_msms_similarity_score(acts, refs):
    """ Calculates MS/MS similarity of reverse dot-product

    Args:
        acts (list[float]): relative intensity vector of measured ions
        refs (list[float]): relative intensity vector of reference ions

    Returns:
        score (float): reverse dot-product between measured and reference ions
    """
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

#region Not used function
def calc_dot_product(acts, refs):
    dotp =  lambda acts, refs: sum(act*ref for act, ref in zip(acts, refs))**2
    scalar = lambda vec: sum(v**2 for v in vec)
    dotp_sml = lambda acts, refs: dotp(acts, refs)/(scalar(acts)*scalar(refs))
    score = dotp_sml(acts, refs)*1000
    return score
#endregion


def calc_similarity_score(acts, refs):
    """ Calculates squared reverse dot-product

    Args:
        acts (list[float]): relative intensity vector of measured ions
        refs (list[float]): relative intensity vector of reference ions

    Returns:
        score (float): dot-product between measured and reference ions
    """
    sqrt_dotp =  lambda acts, refs: sum(math.sqrt(act*ref) for act, ref in zip(acts, refs))
    scalar = lambda vec: sum(vec)
    sim = lambda acts, refs: sqrt_dotp(acts, refs)/math.sqrt(scalar(acts)*scalar(refs))
    score = sim(acts, refs)*1000
    return score
#endregion

#region CID fragmentation analysis
def simulate_acyl_loss(structure_dict):
    """ Calculates neutral loss value of acyl-chains and those of the combinations

    Args:
        structure_dict (dict): {
            'Status', 'Adduct', 'Precursor Mz', 'MS2 Mz', 'Precise precursor Mz',
            'Ref precursor Mz', 'RT(min)', 'Reference RT', 'Ontology', 'Brutto', 
            'Valid moiety num', 'Each moiety info', 'Unsaturated moiety', 
            'Unsaturated sphingobase', 'Deuterium', 'SMILES', 'Formula', 
            'Atom dict', 'Oxidized', 'NL type'
        }

    Returns:
        acyl_loss_dict (dict): {
            'acyl-1 loss', 'acyl-2 loss', 'acyl-3 loss', '2acyls loss',
            'acyl-1/2 loss', 'acyl-2/3 loss', 'acyl-1/3 loss'
        }
    """
    valid_moiety_num = structure_dict['Valid moiety num']
    each_moiety_dict = structure_dict['Each moiety info']
    ref_precursor_mz = structure_dict['Ref precursor Mz']
    acyl_loss_dict = {}
    if valid_moiety_num == 1:
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        chain_1_loss = ref_precursor_mz - (ex_mass['C']*chain_1_num + (2*(chain_1_num-1-db_1_num))*ex_mass['H'] + ex_mass['O'])
        acyl_loss_dict['acyl-1 loss'] = chain_1_loss
    elif valid_moiety_num == 2:
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        chain_2_num = each_moiety_dict['chain-2']
        db_2_num = each_moiety_dict['db-2']
        chain_1_loss = ref_precursor_mz - (ex_mass['C']*chain_1_num + (2*(chain_1_num-1-db_1_num))*ex_mass['H'] + ex_mass['O'])
        chain_2_loss = ref_precursor_mz - (ex_mass['C']*chain_2_num + (2*(chain_2_num-1-db_2_num))*ex_mass['H'] + ex_mass['O'])
        both_acyls_loss = (ref_precursor_mz - (ex_mass['C']*(chain_1_num+chain_2_num) + 
                            2*((chain_1_num+chain_2_num)-2-(db_1_num+db_2_num))*ex_mass['H'] + 2*ex_mass['O']))
        acyl_loss_dict['acyl-1 loss'] = chain_1_loss
        acyl_loss_dict['acyl-2 loss'] = chain_2_loss
        acyl_loss_dict['2acyls loss'] = both_acyls_loss
    elif valid_moiety_num == 3:
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        chain_2_num = each_moiety_dict['chain-2']
        db_2_num = each_moiety_dict['db-2']
        chain_3_num = each_moiety_dict['chain-3']
        db_3_num = each_moiety_dict['db-3']
        chain_1_loss = ref_precursor_mz - (ex_mass['C']*chain_1_num + (2*(chain_1_num-1-db_1_num))*ex_mass['H'] + ex_mass['O'])
        chain_2_loss = ref_precursor_mz - (ex_mass['C']*chain_2_num + (2*(chain_2_num-1-db_2_num))*ex_mass['H'] + ex_mass['O'])
        chain_3_loss = ref_precursor_mz - (ex_mass['C']*chain_3_num + (2*(chain_3_num-1-db_3_num))*ex_mass['H'] + ex_mass['O'])
        acyls_1_2_loss = (ref_precursor_mz - (ex_mass['C']*(chain_1_num+chain_2_num) + 
                            2*((chain_1_num+chain_2_num)-2-(db_1_num+db_2_num))*ex_mass['H'] + 2*ex_mass['O']))
        acyls_1_3_loss = (ref_precursor_mz - (ex_mass['C']*(chain_1_num+chain_3_num) + 
                            2*((chain_1_num+chain_3_num)-2-(db_1_num+db_3_num))*ex_mass['H'] + 2*ex_mass['O']))
        acyls_2_3_loss = (ref_precursor_mz - (ex_mass['C']*(chain_2_num+chain_3_num) + 
                            2*((chain_2_num+chain_3_num)-2-(db_2_num+db_3_num))*ex_mass['H'] + 2*ex_mass['O']))
        acyl_loss_dict['acyl-1 loss'] = chain_1_loss
        acyl_loss_dict['acyl-2 loss'] = chain_2_loss
        acyl_loss_dict['acyl-3 loss'] = chain_3_loss
        acyl_loss_dict['acyls-1/2 loss'] = acyls_1_2_loss
        acyl_loss_dict['acyls-1/3 loss'] = acyls_1_3_loss
        acyl_loss_dict['acyls-2/3 loss'] = acyls_2_3_loss
    return acyl_loss_dict

def simulate_free_fa(structure_dict):
    """ Calculates m/z of free fatty acids

    Args:
        structure_dict (dict): {
            'Status', 'Adduct', 'Precursor Mz', 'MS2 Mz', 'Precise precursor Mz',
            'Ref precursor Mz', 'RT(min)', 'Reference RT', 'Ontology', 'Brutto', 
            'Valid moiety num', 'Each moiety info', 'Unsaturated moiety', 
            'Unsaturated sphingobase', 'Deuterium', 'SMILES', 'Formula', 
            'Atom dict', 'Oxidized', 'NL type'
        }

    Returns:
        free_fa_dict (dict): {'fa-1', 'fa-2', 'fa-3'}
    """
    valid_moiety_num = structure_dict['Valid moiety num']
    each_moiety_dict = structure_dict['Each moiety info']
    free_fa_dict = {}
    if valid_moiety_num == 1:
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        fa_1 = ex_mass['C']*chain_1_num + 2*(chain_1_num-db_1_num)*ex_mass['H'] + 2*ex_mass['O'] - ex_mass['H+']
        free_fa_dict['fa-1'] = fa_1
    elif valid_moiety_num == 2:
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        chain_2_num = each_moiety_dict['chain-2']
        db_2_num = each_moiety_dict['db-2']
        fa_1 = ex_mass['C']*chain_1_num + 2*(chain_1_num-db_1_num)*ex_mass['H'] + 2*ex_mass['O'] - ex_mass['H+']
        free_fa_dict['fa-1'] = fa_1
        fa_2 = ex_mass['C']*chain_2_num + 2*(chain_2_num-db_2_num)*ex_mass['H'] + 2*ex_mass['O'] - ex_mass['H+']
        free_fa_dict['fa-2'] = fa_2
    elif valid_moiety_num == 3:
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        chain_2_num = each_moiety_dict['chain-2']
        db_2_num = each_moiety_dict['db-2']
        chain_3_num = each_moiety_dict['chain-3']
        db_3_num = each_moiety_dict['db-3']
        fa_1 = ex_mass['C']*chain_1_num + 2*(chain_1_num-db_1_num)*ex_mass['H'] + 2*ex_mass['O'] - ex_mass['H+']
        free_fa_dict['fa-1'] = fa_1
        fa_2 = ex_mass['C']*chain_2_num + 2*(chain_2_num-db_2_num)*ex_mass['H'] + 2*ex_mass['O'] - ex_mass['H+']
        free_fa_dict['fa-2'] = fa_2
        fa_3 = ex_mass['C']*chain_3_num + 2*(chain_3_num-db_3_num)*ex_mass['H'] + 2*ex_mass['O'] - ex_mass['H+']
        free_fa_dict['fa-3'] = fa_3
    return free_fa_dict

def simulate_sphingo_base(structure_dict):
    """ Calculates m/z of sphingoid bases

    Args:
        structure_dict (dict): {
            'Status', 'Adduct', 'Precursor Mz', 'MS2 Mz', 'Precise precursor Mz',
            'Ref precursor Mz', 'RT(min)', 'Reference RT', 'Ontology', 'Brutto', 
            'Valid moiety num', 'Each moiety info', 'Unsaturated moiety', 
            'Unsaturated sphingobase', 'Deuterium', 'SMILES', 'Formula', 
            'Atom dict', 'Oxidized', 'NL type'
        }

    Returns:
        sphingo_base_dict (dict): {
            'Sph_base', 'Sph-H2O', 'PhytoSph_base', 'PhytoSph-H2O'
        }
    """
    each_moiety_dict = structure_dict['Each moiety info']
    sphingo_base_dict = {}
    chain_1_num = each_moiety_dict['chain-1']
    db_1_num = each_moiety_dict['db-1']
    sphingo_base_dict['Sph_base'] = chain_1_num*ex_mass['C'] + (2*(chain_1_num-db_1_num)+3)*ex_mass['H'] + ex_mass['O']*2 + ex_mass['N']
    sphingo_base_dict['Sph-H2O'] = sphingo_base_dict['Sph_base'] - ex_mass['H2O']
    sphingo_base_dict['PhytoSph_base'] = chain_1_num*ex_mass['C'] + (2*(chain_1_num-db_1_num)+3)*ex_mass['H'] + ex_mass['O']*3 + ex_mass['N']
    sphingo_base_dict['PhytoSph-H2O'] = sphingo_base_dict['PhytoSph_base'] - ex_mass['H2O']
    return sphingo_base_dict

def simulate_sphingo_base_for_acylhex(structure_dict):
    """ Calculates m/z of sphingoid bases for acyl-hexosyl-ceramides

    Args:
        structure_dict (dict): {
            'Status', 'Adduct', 'Precursor Mz', 'MS2 Mz', 'Precise precursor Mz',
            'Ref precursor Mz', 'RT(min)', 'Reference RT', 'Ontology', 'Brutto', 
            'Valid moiety num', 'Each moiety info', 'Unsaturated moiety', 
            'Unsaturated sphingobase', 'Deuterium', 'SMILES', 'Formula', 
            'Atom dict', 'Oxidized', 'NL type'
        }

    Returns:
        sphingo_base_dict (dict): {
            'Sph_base', 'Sph-H2O', 'PhytoSph_base', 'PhytoSph-H2O'
        }
    """
    each_moiety_dict = structure_dict['Each moiety info']
    sphingo_base_for_acylhex_dict = {}
    chain_2_num = each_moiety_dict['chain-2']
    db_2_num = each_moiety_dict['db-2']
    sphingo_base_for_acylhex_dict['Sph_base'] = chain_2_num*ex_mass['C'] + (2*(chain_2_num-db_2_num)+3)*ex_mass['H'] + ex_mass['O']*2 + ex_mass['N']
    sphingo_base_for_acylhex_dict['Sph-H2O'] = sphingo_base_for_acylhex_dict['Sph_base'] - ex_mass['H2O']
    sphingo_base_for_acylhex_dict['PhytoSph_base'] = chain_2_num*ex_mass['C'] + (2*(chain_2_num-db_2_num)+3)*ex_mass['H'] + ex_mass['O']*3 + ex_mass['N']
    sphingo_base_for_acylhex_dict['PhytoSph-H2O'] = sphingo_base_for_acylhex_dict['PhytoSph_base'] - ex_mass['H2O']
    return sphingo_base_for_acylhex_dict

def search_cid_fragment_ions(structure_dict, dataframe, ms_tolerance_ppm):
    """ Query reference fragment ions based on 
        collision induced dissociation (CID) library

    Args:
        structure_dict (dict): {
            'Status', 'Adduct', 'Precursor Mz', 'MS2 Mz', 'Precise precursor Mz',
            'Ref precursor Mz', 'RT(min)', 'Reference RT', 'Ontology', 'Brutto', 
            'Valid moiety num', 'Each moiety info', 'Unsaturated moiety', 
            'Unsaturated sphingobase', 'Deuterium', 'SMILES', 'Formula', 
            'Atom dict', 'Oxidized', 'NL type'
        }
        dataframe (Dataframe): Extracted MS/MS spectrum
        ms_tolerance_ppm (int): mass tolerance 

    Returns:
        cid_fragmentation_dict (dict): 
            dict contains diagnostic ions based on CID fragmentation
    """
    current_lipidclass = structure_dict['Ontology']
    adduct_type = structure_dict['Adduct']
    ref_precursor_mz = structure_dict['Ref precursor Mz']
    valid_moiety_num = structure_dict['Valid moiety num']
    
    each_moiety_dict = structure_dict['Each moiety info']
    acyl_loss_dict = simulate_acyl_loss(structure_dict)
    free_fa_dict = simulate_free_fa(structure_dict)
    sphingo_base_dict = simulate_sphingo_base(structure_dict)
    if len(each_moiety_dict) >= 4:
        sphingo_base_for_acylhex_dict = simulate_sphingo_base_for_acylhex(structure_dict)
    diagnostic_ions_dict = {}
    #region CID fragmetation library
    if current_lipidclass == 'FA' and adduct_type == '[M-H]-':
        pass

    if current_lipidclass == 'NAGly' and (adduct_type == '[M+H]+' or adduct_type == '[M+NH4]+'):
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        diagnostic_ions_dict['Glycine'] = 76.03930542
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
        if valid_moiety_num == 2:
            chain_2_num = each_moiety_dict['chain-2']
            db_2_num = each_moiety_dict['db-2']
            diagnostic_ions_dict['acyl-2-H2O'] = ex_mass['C']*chain_2_num + 2*ex_mass['H']*(chain_2_num - 1 - db_2_num) + ex_mass['O'] - ex_mass['H'] - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2'] = ex_mass['C']*chain_2_num + 2*ex_mass['H']*(chain_2_num - 1 - db_2_num) + ex_mass['O'] - ex_mass['H']

    if current_lipidclass == 'NAGly' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['Glycine'] = 74.02475258
        diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
        diagnostic_ions_dict['acyl-1+H2O+CO2 loss'] = acyl_loss_dict['acyl-1 loss'] - (ex_mass['H2O'] + ex_mass['CO2'])

    if current_lipidclass == 'NAGlySer' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['GlySer-O'] = 145.0618666
        diagnostic_ions_dict['Serine'] = 106.0498704
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
        if valid_moiety_num == 2:
            chain_2_num = each_moiety_dict['chain-2']
            db_2_num = each_moiety_dict['db-2']
            diagnostic_ions_dict['acyl-2-Gly'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O'] - (106.0499 - ex_mass['H'])
            diagnostic_ions_dict['acyl-2'] = ex_mass['C']*chain_2_num + 2*ex_mass['H']*(chain_2_num - 1 - db_2_num) + ex_mass['O'] - ex_mass['H']

    if current_lipidclass == 'NAGlySer' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
        diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
        diagnostic_ions_dict['acyl-1+H2O+CH2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O'])
        diagnostic_ions_dict['acyl-1+2H2O+CO2 loss'] = acyl_loss_dict['acyl-1 loss'] - 2*ex_mass['H2O'] - ex_mass['CO2']

    if current_lipidclass == 'NAOrn' and adduct_type == '[M+H]+':
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        diagnostic_ions_dict['Ornithine C4H8N+'] = 70.063
        diagnostic_ions_dict['Ornithine -H2O'] = 115.086
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
        diagnostic_ions_dict['acyl-1+2H2O loss'] = acyl_loss_dict['acyl-1 loss'] - 2*ex_mass['H2O']
        if valid_moiety_num == 2:
            chain_2_num = each_moiety_dict['chain-2']
            db_2_num = each_moiety_dict['db-2']

    if current_lipidclass == 'NAE' and adduct_type == '[M+H]+':
        pass

    if current_lipidclass == 'CAR' and adduct_type == '[M]+':
        diagnostic_ions_dict['C4H5O2+'] = 85.028406

    if current_lipidclass == 'FAHFA' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
        diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']

    if current_lipidclass == 'DG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['NH3+H2O loss'] = ref_precursor_mz - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 DMAG'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']

    if current_lipidclass == 'EtherDG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'MGDG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C6H10O5+NH3+H2O loss'] = ref_precursor_mz - 179.079374 - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 DMAG'] = acyl_loss_dict['acyl-2 loss'] - 179.079374 - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - 179.079374 - ex_mass['H2O']

    if current_lipidclass == 'MGDG' and adduct_type == '[M+CH3COO]-':
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'EtherMGDG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C6H10O5+NH3+H2O loss'] = ref_precursor_mz - 179.079374 - ex_mass['H2O']
        diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - 179.079374 - ex_mass['H2O']
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'EtherMGDG' and adduct_type == '[M+CH3COO]-':
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']
        diagnostic_ions_dict['acyl-2 loss'] = acyl_loss_dict['acyl-1 loss'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'DGDG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C12H20O10+NH3+H2O loss'] = ref_precursor_mz - 341.132199 - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 DMAG'] = acyl_loss_dict['acyl-2 loss'] - 341.132199 - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - 341.132199 - ex_mass['H2O']

    if current_lipidclass == 'DGDG' and adduct_type == '[M+CH3COO]-':
        if valid_moiety_num == 2:
            diagnostic_ions_dict['2acyls loss'] = acyl_loss_dict['2acyls loss'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+'] - ex_mass['H2O']
            diagnostic_ions_dict['2acyls+H2O loss'] = acyl_loss_dict['2acyls loss'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+'] - 2*ex_mass['H2O']
            if valid_moiety_num == 2:
                diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
                diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'SQDG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C6H10O7S loss'] = ref_precursor_mz - 226.014723836 - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 DMAG'] = acyl_loss_dict['acyl-2 loss'] - 226.014723836 - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+'] - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - 226.014723836 - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+'] - ex_mass['H2O']

    if current_lipidclass == 'SQDG' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['Sulfoquinovosyl-H2O'] = 225.007447384
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']
            diagnostic_ions_dict['acyl-1 loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2 loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O']

    if current_lipidclass == 'MG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['NH3 loss'] = ref_precursor_mz - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
        diagnostic_ions_dict['NH3+H2O loss'] = ref_precursor_mz - adduct_dict_pos['[M+NH4]+'] - ex_mass['H2O'] + ex_mass['H+']

    if current_lipidclass == 'ADGGA' and adduct_type == '[M+NH4]+':
        chain_3_num = each_moiety_dict['chain-3']
        db_3_num = each_moiety_dict['db-3']
        acyl_3 = ex_mass['C']*chain_3_num + 2*ex_mass['H']*(chain_3_num - db_3_num) + ex_mass['O']
        grucronic_acid = ex_mass['C']*6 + ex_mass['H']*6 + ex_mass['O']*6
        diagnostic_ions_dict['AcylGrucroniAcid'] = acyl_3 + grucronic_acid + ex_mass['H+']
        diagnostic_ions_dict['acyl-1 DMAG'] = (acyl_loss_dict['acyl-2 loss'] - 
            acyl_3 - grucronic_acid - adduct_dict_pos['[M+NH4]+'] - ex_mass['H2O'] + ex_mass['H+'])
        diagnostic_ions_dict['acyl-2 DMAG'] = (acyl_loss_dict['acyl-1 loss'] - 
            acyl_3 - grucronic_acid - adduct_dict_pos['[M+NH4]+'] - ex_mass['H2O'] + ex_mass['H+'])

    if current_lipidclass == 'ADGGA' and adduct_type == '[M-H]-':
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'DGCC' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C6H14NO2+'] = 132.101905115
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 loss'] = acyl_loss_dict['acyl-1 loss']
            diagnostic_ions_dict['acyl-2 loss'] = acyl_loss_dict['acyl-2 loss']
            diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O']

    if current_lipidclass == 'DGTS/A' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C7H14NO2+'] = 144.101905115
        diagnostic_ions_dict['C10H22NO5+'] = 236.149249233
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 loss'] = acyl_loss_dict['acyl-1 loss']
            diagnostic_ions_dict['acyl-2 loss'] = acyl_loss_dict['acyl-2 loss']
            diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O']

    if current_lipidclass == 'DGGA' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['GrucroniAcid+NH3 loss'] = ref_precursor_mz - 194.042652622 - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 DMAG'] = (acyl_loss_dict['acyl-2 loss'] - 
                194.042652622 - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+'])
            diagnostic_ions_dict['acyl-2 DMAG'] = (acyl_loss_dict['acyl-1 loss'] - 
                194.042652622 - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+'])

    if current_lipidclass == 'DGGA' and adduct_type == '[M-H]-':
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'LDGCC' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C6H14NO2+'] = 132.101905115

    if current_lipidclass == 'LDGTS/A' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C7H14NO2+'] = 144.101905115
        diagnostic_ions_dict['C10H22NO5+'] = 236.149249233

    if current_lipidclass == 'TG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - adduct_dict_pos['[M+NH4]+'] - ex_mass['H2O'] + ex_mass['H+']
        diagnostic_ions_dict['acyl-2+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - adduct_dict_pos['[M+NH4]+'] - ex_mass['H2O'] + ex_mass['H+']
        diagnostic_ions_dict['acyl-3+H2O loss'] = acyl_loss_dict['acyl-3 loss'] - adduct_dict_pos['[M+NH4]+'] - ex_mass['H2O'] + ex_mass['H+']

    if current_lipidclass == 'TG' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
        diagnostic_ions_dict['acyl-2+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O']
        diagnostic_ions_dict['acyl-3+H2O loss'] = acyl_loss_dict['acyl-3 loss'] - ex_mass['H2O']

    if current_lipidclass == 'EtherTG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['ether-link loss'] = acyl_loss_dict['acyl-1 loss'] - adduct_dict_pos['[M+NH4]+'] - 4*ex_mass['H'] + ex_mass['H+']
        diagnostic_ions_dict['acyl-2+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - adduct_dict_pos['[M+NH4]+'] - ex_mass['H2O'] + ex_mass['H+']
        diagnostic_ions_dict['acyl-3+H2O loss'] = acyl_loss_dict['acyl-3 loss'] - adduct_dict_pos['[M+NH4]+'] - ex_mass['H2O'] + ex_mass['H+']
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'EtherTG' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['ether-link loss'] = acyl_loss_dict['acyl-1 loss'] - 4*ex_mass['H']
        diagnostic_ions_dict['acyl-2+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O']
        diagnostic_ions_dict['acyl-3+H2O loss'] = acyl_loss_dict['acyl-3 loss'] - ex_mass['H2O']
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'PA' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'LPA' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394

    if current_lipidclass == 'PC' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C5H15NO4P+'] = 184.073320905
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 loss'] = acyl_loss_dict['acyl-1 loss']
            diagnostic_ions_dict['acyl-2 loss'] = acyl_loss_dict['acyl-2 loss']
            diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O']
    
    if current_lipidclass == 'PC' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['C2H5O4PNa+'] = 146.98176634
        diagnostic_ions_dict['C3H9N loss'] = ref_precursor_mz - (ex_mass['C']*3 + ex_mass['H']*9 + ex_mass['N'])
        diagnostic_ions_dict['C5H14NO4P loss'] = ref_precursor_mz - 183.06604
    
    if current_lipidclass == 'PC' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['CH3 loss'] = ref_precursor_mz - adduct_dict_neg['[M+CH3COO]-'] - 15.0234738
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'EtherPC' and adduct_type == '[M+H]+':
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02
        diagnostic_ions_dict['C5H15NO4P+'] = 184.073320905

    if current_lipidclass == 'EtherPC' and adduct_type == '[M+CH3COO]-':
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['CH3 loss'] = ref_precursor_mz - adduct_dict_neg['[M+CH3COO]-'] - 15.0234738
        diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02    
    
    if current_lipidclass == 'LPC' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C5H15NO4P+'] = 184.073320905
        diagnostic_ions_dict['C5H14NO+'] = 104.106990495

    if current_lipidclass == 'LPC' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['C3H9N loss'] = ref_precursor_mz - (ex_mass['C']*3 + ex_mass['H']*9 + ex_mass['N'])

    if current_lipidclass == 'LPC' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['CH3 loss'] = ref_precursor_mz - adduct_dict_neg['[M+CH3COO]-'] - 15.0234738

    if current_lipidclass == 'EtherLPC' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C5H15NO4P+'] = 184.073320905
        diagnostic_ions_dict['C2H6O4P+'] = 124.999821612
        diagnostic_ions_dict['C5H14NO+'] = 104.106990495
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'EtherLPC' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['CH3 loss'] = ref_precursor_mz - adduct_dict_neg['[M+CH3COO]-'] - 15.0234738
        diagnostic_ions_dict['C4H11NO4P-'] = 168.043117937
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'PE' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C2H8NO4P loss'] = ref_precursor_mz - 141.019094261
        if valid_moiety_num == 2:
            chain_1_num = each_moiety_dict['chain-1']
            db_1_num = each_moiety_dict['db-1']
            chain_2_num = each_moiety_dict['chain-2']
            db_2_num = each_moiety_dict['db-2']
            diagnostic_ions_dict['acyl-1'] = ex_mass['C']*chain_1_num + 2*ex_mass['H']*(chain_1_num - db_1_num) + ex_mass['O'] - ex_mass['H']
            diagnostic_ions_dict['acyl-2'] = ex_mass['C']*chain_2_num + 2*ex_mass['H']*(chain_2_num - db_2_num) + ex_mass['O'] - ex_mass['H']
        
    if current_lipidclass == 'PE' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['C2H5N loss'] = ref_precursor_mz - (ex_mass['C']*2 + ex_mass['H']*5 + ex_mass['N'])
        diagnostic_ions_dict['C2H8NO4P loss'] = ref_precursor_mz - (ex_mass['C']*2 + ex_mass['H']*8 + ex_mass['N'] + ex_mass['O']*4 + ex_mass['P'])

    if current_lipidclass == 'PE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C5H11NO5P-'] = 196.038032559
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'EtherPE' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C2H8NO4P loss'] = ref_precursor_mz - 141.019094261
        position = each_moiety_dict['chain-1'] - 1
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['ester-link+C2H8NO3P'] = (ex_mass['C']*chain_1_num 
            + 2*ex_mass['H']*(chain_1_num - db_1_num + 1) + ex_mass['O'] + 124.01635)
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02    
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - (141.019094261 - ex_mass['O'] + ex_mass['H']*2)

    if current_lipidclass == 'EtherPE' and adduct_type == '[M-H]-':
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-2 loss'] = acyl_loss_dict['acyl-2 loss']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02    

    if current_lipidclass == 'EtherPE(P)' and adduct_type == '[M+H]+':
        position = each_moiety_dict['chain-1'] - 1
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['C2H8NO4P loss'] = ref_precursor_mz - 141.019094261
        diagnostic_ions_dict['ester-link+C2H8NO3P'] = (ex_mass['C']*chain_1_num 
            + 2*ex_mass['H']*(chain_1_num - db_1_num + 1) + ex_mass['O'] + 124.01635)
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02    
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - (141.019094261 - ex_mass['O'] + ex_mass['H']*2)

    if current_lipidclass == 'PlasmPE' and adduct_type == '[M+H]+':
        position = each_moiety_dict['chain-1'] - 1
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['C2H8NO4P loss'] = ref_precursor_mz - 141.019094261
        diagnostic_ions_dict['ester-link+C2H8NO3P'] = (ex_mass['C']*chain_1_num + 2*ex_mass['H']*(chain_1_num - db_1_num + 1) 
                                                        + ex_mass['O'] + 124.01635)
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02  
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - (141.019094261 - ex_mass['O'] + ex_mass['H']*2)

    if current_lipidclass == 'NAPE' and adduct_type == '[M+H]+':
        pass

    if current_lipidclass == 'LPE' and (adduct_type == '[M+H]+' or adduct_type == '[M+Na]+'):
        diagnostic_ions_dict['C2H8NO4P loss'] = ref_precursor_mz - 141.019094261

    if current_lipidclass == 'LPE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']

    if current_lipidclass == 'EtherLPE' and adduct_type == '[M+H]+':
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['C3H8NO4P loss'] = ref_precursor_mz - (ex_mass['C']*3+ex_mass['H']*8+ex_mass['N']+ex_mass['O']*4+ex_mass['P']) - ex_mass['H']
        diagnostic_ions_dict['C3H10NO5P loss'] = ref_precursor_mz - (ex_mass['C']*3+ex_mass['H']*10+ex_mass['N']+ex_mass['O']*5+ex_mass['P']) - ex_mass['H']
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'EtherLPE' and adduct_type == '[M-H]-':
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['C2H7NO4P-'] = 140.0118196
        diagnostic_ions_dict['C2H6N+H2O loss'] = ref_precursor_mz - (ex_mass['C']*2+ex_mass['H']*6+ex_mass['N']+ex_mass['H2O']) + ex_mass['H+']
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02    

    if current_lipidclass == 'LNAPE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
        diagnostic_ions_dict['acyl-1 loss'] = acyl_loss_dict['acyl-1 loss']
        diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']

    if current_lipidclass == 'PG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C3H8O6P loss'] = ref_precursor_mz - 171.0053034 - adduct_dict_pos['[M+NH4]+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 DMAG'] = acyl_loss_dict['acyl-2 loss'] - 171.0053034 - adduct_dict_pos['[M+NH4]+']
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - 171.0053034 - adduct_dict_pos['[M+NH4]+']

    if current_lipidclass == 'PG' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['C3H8O6P loss'] = ref_precursor_mz - 171.0053034 - adduct_dict_pos['[M+Na]+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 DMAG'] = acyl_loss_dict['acyl-2 loss'] - 171.0053034 - adduct_dict_pos['[M+Na]+']
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - 171.0053034 - adduct_dict_pos['[M+Na]+']

    if current_lipidclass == 'PG' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'EtherPG' and adduct_type == '[M-H]-':
        chain_1_num = each_moiety_dict['chain-1']
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        diagnostic_ions_dict['acyl-1+O'] =  ex_mass['C']*chain_1_num + 2*ex_mass['H']*(chain_1_num + 1 - db_1_num) + 2*ex_mass['O'] - ex_mass['H+']
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'LPG' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C3H8O6P loss'] = ref_precursor_mz - 171.0053034 - adduct_dict_pos['[M+NH4]+']
    if current_lipidclass == 'LPG' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C3H8O6P loss'] = ref_precursor_mz - 171.0053034 - ex_mass['H']

    if current_lipidclass == 'LPG' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394

    if current_lipidclass == 'EtherLPG' and adduct_type == '[M-H]-':
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        diagnostic_ions_dict['ether-link'] = ex_mass['C']*chain_1_num + 2*ex_mass['H']*(chain_1_num - db_1_num + 1) + ex_mass['O'] - ex_mass['H+']
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        position = each_moiety_dict['chain-1'] - 1
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'BMP' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C3H8O6P loss'] = ref_precursor_mz - 171.0053034 - adduct_dict_pos['[M+NH4]+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 DMAG'] = acyl_loss_dict['acyl-2 loss'] - 171.0053034 - adduct_dict_pos['[M+NH4]+']
            diagnostic_ions_dict['acyl-2 DMAG'] = acyl_loss_dict['acyl-1 loss'] - 171.0053034 - adduct_dict_pos['[M+NH4]+']
    
    if current_lipidclass == 'HBMP' and adduct_type == '[M+NH4]+':
        chain_1_num = each_moiety_dict['chain-1']
        db_1_num = each_moiety_dict['db-1']
        chain_2_num = each_moiety_dict['chain-2']
        db_2_num = each_moiety_dict['db-2']
        chain_3_num = each_moiety_dict['chain-3']
        db_3_num = each_moiety_dict['db-3']
        glycerol = ex_mass['C']*3 + ex_mass['H']*5 + ex_mass['O']*2
        diagnostic_ions_dict['acyl-1+H3PO4 loss'] = (acyl_loss_dict['acyl-1 loss'] - 
            (ex_mass['H']*3+ex_mass['P']+ex_mass['O']*4 + glycerol + adduct_dict_pos['[M+NH4]+']))
        diagnostic_ions_dict['acyl-1'] = ex_mass['C']*chain_1_num + 2*ex_mass['H']*(chain_1_num - db_1_num) + ex_mass['O'] + glycerol
        diagnostic_ions_dict['acyl-2'] = ex_mass['C']*chain_2_num + 2*ex_mass['H']*(chain_2_num - db_2_num) + ex_mass['O'] + glycerol
        diagnostic_ions_dict['acyl-3'] = ex_mass['C']*chain_3_num + 2*ex_mass['H']*(chain_3_num - db_3_num) + ex_mass['O'] + glycerol

    if current_lipidclass == 'HBMP' and adduct_type == '[M-H]-':
        if valid_moiety_num == 3:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']
            diagnostic_ions_dict['fa-3'] = free_fa_dict['fa-3']
            diagnostic_ions_dict['acyl-1+C3H6O4P'] = free_fa_dict['fa-1'] + (ex_mass['C']*3+ex_mass['H']*6+ex_mass['O']*4+ex_mass['P']) - ex_mass['H+']

    if current_lipidclass == 'CL' and adduct_type == '[M+NH4]+':
        pass
    if current_lipidclass == 'CL' and adduct_type == '[M-H]-':
        pass

    if current_lipidclass == 'MLCL' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        if valid_moiety_num == 3:
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']
            diagnostic_ions_dict['fa-3'] = free_fa_dict['fa-3']
            glycerol_phosphate = ex_mass['C']*3+ex_mass['H']*6+ex_mass['O']*5+ex_mass['P']
            diagnostic_ions_dict['acyl-1+C3H6O5P'] = acyl_loss_dict['acyls-2/3 loss'] - glycerol_phosphate
            diagnostic_ions_dict['acyl-2-3+C3H6O5P'] = acyl_loss_dict['acyl-1 loss'] - glycerol_phosphate - ex_mass['H+']

    if current_lipidclass == 'DLCL' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'PI' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C6H13O9P+NH3 loss'] = ref_precursor_mz - 277.056268097

    if current_lipidclass == 'PI' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['C6H13O9P loss'] = ref_precursor_mz - 260.029723
        diagnostic_ions_dict['C6H13O9P+Na loss'] = ref_precursor_mz - 260.029723 - adduct_dict_pos['[M+Na]+'] + ex_mass['H+']
        diagnostic_ions_dict['C6H13O9PNa+'] = 283.0189444

    if current_lipidclass == 'PI' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C6H10O8P-'] = 241.01188
        diagnostic_ions_dict['C9H14O9P-'] = 297.03809
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'EtherPI' and adduct_type == '[M-H]-':
        chain_1_num = each_moiety_dict['chain-1']
        position = each_moiety_dict['chain-1'] - 1
        db_1_num = each_moiety_dict['db-1']
        mass_shift = 0
        if (db_1_num - 1) > 0:
            mass_shift = 2*ex_mass['H']*(db_1_num - 1)
        ref_oad_01 = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O'] - mass_shift
        ref_oad_01_rad = (position-1)*ex_mass['C']+2*(position-1)*ex_mass['H']-ex_mass['O']+ex_mass['H'] - mass_shift
        ref_oad_02 = (position-1)*ex_mass['C']+2*position*ex_mass['H']-ex_mass['O'] - mass_shift
        diagnostic_ions_dict['C6H10O8P-'] = 241.01188
        diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']
        diagnostic_ions_dict['ether-link+C3H5O4P'] = acyl_loss_dict['acyl-2 loss'] - 180.06339
        diagnostic_ions_dict['Plasmalogen OAD1'] = ref_precursor_mz - ref_oad_01
        diagnostic_ions_dict['Plasmalogen OAD1rad'] = ref_precursor_mz - ref_oad_01_rad
        diagnostic_ions_dict['Plasmalogen OAD2'] = ref_precursor_mz - ref_oad_02

    if current_lipidclass == 'LPI' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['C6H13O9P+NH3 loss'] = ref_precursor_mz - 277.056268097
    if current_lipidclass == 'LPI' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C6H13O9P loss'] = ref_precursor_mz - 260.029723 - ex_mass['H']

    if current_lipidclass == 'LPI' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C6H10O8P-'] = 241.01188
        diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
        diagnostic_ions_dict['acyl-1 loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']

    if current_lipidclass == 'PS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['PhosphoSerine loss'] = ref_precursor_mz - 185.00892
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-1 loss'] = acyl_loss_dict['acyl-1 loss']
            diagnostic_ions_dict['acyl-2 loss'] = acyl_loss_dict['acyl-2 loss']
            diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
            diagnostic_ions_dict['acyl-2+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O']

    if current_lipidclass == 'PS' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['PhosphoSerine loss'] = ref_precursor_mz - 185.00892

    if current_lipidclass == 'PS' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['Serine loss'] = ref_precursor_mz - 87.03203
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'LPS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['PhosphoSerine loss'] = ref_precursor_mz - 185.00892

    if current_lipidclass == 'LPS' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['Serine loss'] = ref_precursor_mz - 87.03203
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394

    if current_lipidclass == 'LNAPS' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        diagnostic_ions_dict['N-acyl+Serine loss'] = acyl_loss_dict['acyl-2 loss'] - 87.03203

    if current_lipidclass == 'PMeOH' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['HPO4CH3-'] = 110.985269187
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'PEtOH' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['HPO4C2H5-'] = 125.000919251
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['fa-2'] = free_fa_dict['fa-2']

    if current_lipidclass == 'OxPC' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['CH3 loss'] = ref_precursor_mz - adduct_dict_neg['[M+CH3COO]-'] - 15.0234738
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']

    if current_lipidclass == 'OxPE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C5H11NO5P-'] = 196.038032559
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']

    if current_lipidclass == 'EtherOxPE' and adduct_type == '[M-H]-':
        pass

    if current_lipidclass == 'OxPG' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C3H6O5P-'] = 152.995833394
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']

    if current_lipidclass == 'OxPI' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C6H10O8P-'] = 241.01188
        diagnostic_ions_dict['C9H14O9P-'] = 297.03809
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']

    if current_lipidclass == 'OxPS' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['Serine loss'] = ref_precursor_mz - 87.03203
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa-1'] = free_fa_dict['fa-1']

    if current_lipidclass == 'VAE' and (adduct_type == '[M+H]+' or adduct_type == '[M+Na]+'):
        diagnostic_ions_dict['acyl-1 loss'] = 269.2263771

    if current_lipidclass == 'CoQ' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C10H13O4+'] = 197.08081642

    if current_lipidclass == 'Vitamine E' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        diagnostic_ions_dict['C10H11O2-'] = 163.0753564

    if current_lipidclass == 'GM3' and (adduct_type == '[M+H]+' or adduct_type == '[M+NH4]+'):
        if adduct_type == '[M+H]+':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
            diagnostic_ions_dict['C23H35NO17 loss'] = ref_precursor_mz - (617.21671342+ex_mass['H2O']-ex_mass['H']*2)
        elif adduct_type == '[M+NH4]+':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
            diagnostic_ions_dict['C23H35NO17 loss'] = ref_precursor_mz - (617.21671342+ex_mass['H2O']-ex_mass['H']*2) - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
        diagnostic_ions_dict['C11H18NO8+'] = 291.095416511 + ex_mass['H+']
        diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
        diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']

    if current_lipidclass == 'GM3' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C11H16NO8-'] = 291.095416511 - ex_mass['H+']

    if current_lipidclass == 'SHexCer' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['H2SO4 loss'] = ref_precursor_mz - (ex_mass['H']*2+ex_mass['S']+ex_mass['O']*4)
        diagnostic_ions_dict['H2SO4+C6H10O5 loss'] = ref_precursor_mz - (244.02528852 + ex_mass['H2O'] - ex_mass['H']) + ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'SHexCer' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['HSO4-'] = (ex_mass['H']*2 + ex_mass['S'] + ex_mass['O']*4) - ex_mass['H+']

    if current_lipidclass == 'SHexCer+O' and adduct_type == '[M+H]+':
        pass
    if current_lipidclass == 'SHexCer+O' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['HSO4-'] = (ex_mass['H']*2 + ex_mass['S'] + ex_mass['O']*4) - ex_mass['H+']
    
    if current_lipidclass == 'Cer_EOS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
        if valid_moiety_num >= 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'Cer_EOS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if valid_moiety_num == 3:
            diagnostic_ions_dict['N-acyl fa'] = free_fa_dict['fa-3']
            diagnostic_ions_dict['O-acyl loss'] = acyl_loss_dict['acyl-3 loss']
            diagnostic_ions_dict['O-acyl amide'] = free_fa_dict['fa-2'] + ex_mass['N'] + ex_mass['H']

    if current_lipidclass == 'Cer_EBDS' and adduct_type == '[M+CH3COO]-':
        if valid_moiety_num == 3:
            diagnostic_ions_dict['N-acyl fa'] = free_fa_dict['fa-3']
            diagnostic_ions_dict['N-acyl loss'] = acyl_loss_dict['acyl-3 loss']
            diagnostic_ions_dict['O-acyl fa'] = free_fa_dict['fa-2'] + ex_mass['N'] + ex_mass['H'] - (ex_mass['C']*3+ex_mass['H']*5+ex_mass['O'])
    
    if current_lipidclass == 'CerP' and adduct_type == '[M+H]+':
        pass
    if current_lipidclass == 'Cer_AP' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
        diagnostic_ions_dict['2H2O loss'] = ref_precursor_mz - ex_mass['H2O']*2
        if valid_moiety_num == 2:
            diagnostic_ions_dict['PhytoSph'] = sphingo_base_dict['PhytoSph_base'] + ex_mass['H+']
            diagnostic_ions_dict['PhytoSph-H2O'] = sphingo_base_dict['PhytoSph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['PhytoSph-3H2O'] = sphingo_base_dict['PhytoSph-H2O'] - ex_mass['H2O']*2 + ex_mass['H+']

    if current_lipidclass == 'Cer_AP' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa+NC(C)CO']= free_fa_dict['fa-2'] + 70.028745 + ex_mass['H+']
            diagnostic_ions_dict['fa+O']= free_fa_dict['fa-2'] + ex_mass['O']
            diagnostic_ions_dict['fa-CO-3H']= free_fa_dict['fa-2'] - (ex_mass['C']+ex_mass['O']+ex_mass['H']*2)

    if current_lipidclass == 'Cer_NP' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
            diagnostic_ions_dict['2H2O loss'] = ref_precursor_mz - ex_mass['H2O']*2
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['2H2O loss'] = ref_precursor_mz - ex_mass['H2O']*2 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa+NC(C)CO']= free_fa_dict['fa-2'] + 70.028745 + ex_mass['H+'] - ex_mass['O']
            diagnostic_ions_dict['N-acyl fa']= free_fa_dict['fa-2'] + ex_mass['N'] + ex_mass['H+'] - ex_mass['O']
            diagnostic_ions_dict['PhytoSph-CH5NO'] = sphingo_base_dict['PhytoSph_base'] - 50.061143

    if current_lipidclass == 'Cer_ADS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
            diagnostic_ions_dict['CH4O2 loss'] = ref_precursor_mz - 48.02114
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['CH4O2 loss'] = ref_precursor_mz - 48.02114 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph_base'] = sphingo_base_dict['Sph_base'] - ex_mass['H+']
            diagnostic_ions_dict['fa-3H'] = free_fa_dict['fa-2'] - ex_mass['H']*2
            diagnostic_ions_dict['fa-CO-3H'] = free_fa_dict['fa-2'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O'])

    if current_lipidclass == 'Cer_BDS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph+C2H2O'] = sphingo_base_dict['Sph_base'] + (ex_mass['C']*2+ex_mass['H']*2+ex_mass['O']) - ex_mass['H+']
            diagnostic_ions_dict['Sph+C-2H'] = sphingo_base_dict['Sph_base'] + ex_mass['C'] -ex_mass['H']*2 - ex_mass['H+']
            diagnostic_ions_dict['Sph-NCCO-3H'] = sphingo_base_dict['Sph_base'] -(ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']) - ex_mass['H+']

    if current_lipidclass == 'Cer_NDS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'Cer_NDS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['CH2O loss'] = ref_precursor_mz - (ex_mass['C']+ex_mass['H']*2+ex_mass['O'])
            diagnostic_ions_dict['CH4O2 loss'] = ref_precursor_mz - 48.02114
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['CH2O loss'] = ref_precursor_mz - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['CH4O2 loss'] = ref_precursor_mz - 48.02114 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-3H'] = free_fa_dict['fa-2'] - (ex_mass['O']+ex_mass['H']*2)
            diagnostic_ions_dict['acyl+C2H3N'] = free_fa_dict['fa-2'] + (ex_mass['C']*2+ex_mass['H']*3+ex_mass['N'])
            diagnostic_ions_dict['Sph-NCCO-3H'] = sphingo_base_dict['Sph_base'] -(ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']) - ex_mass['H+']

    if current_lipidclass == 'Cer_AS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
            diagnostic_ions_dict['CH4O2 loss'] = ref_precursor_mz - 48.02114
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['CH4O2 loss'] = ref_precursor_mz - 48.02114 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl-CH2O'] = free_fa_dict['fa-2'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O'])
            diagnostic_ions_dict['Sph-NCCO-H'] = sphingo_base_dict['Sph_base'] -(ex_mass['C']*2+ex_mass['H']*5+ex_mass['N']+ex_mass['O']) - ex_mass['H+']
            diagnostic_ions_dict['Sph-NCCO-3H'] = sphingo_base_dict['Sph_base'] -(ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']) - ex_mass['H+']

    if current_lipidclass == 'Cer_BS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-NCCO-H'] = sphingo_base_dict['Sph_base'] - (ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']) - ex_mass['H+']
            diagnostic_ions_dict['Sph+C2H2O'] = sphingo_base_dict['Sph_base'] + (ex_mass['C']*2+ex_mass['H']*2+ex_mass['O']) - ex_mass['H+']

    if current_lipidclass == 'Cer_NS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'Cer_NS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
            diagnostic_ions_dict['CH2O loss'] = ref_precursor_mz - (ex_mass['C']+ex_mass['H']*2+ex_mass['O'])
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['CH2O loss'] = ref_precursor_mz - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl+C2H3N'] = free_fa_dict['fa-2'] + (ex_mass['C']*2+ex_mass['H']*2+ex_mass['N']+ex_mass['H+']) - ex_mass['O']
            diagnostic_ions_dict['Sph-NCCO-3H'] = sphingo_base_dict['Sph_base'] -(ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']) - ex_mass['H+']
            diagnostic_ions_dict['acyl-3H'] = free_fa_dict['fa-2'] - (ex_mass['O']+ex_mass['H']*2)

    if current_lipidclass == 'Cer_HS' and adduct_type == '[M+H]+':
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'Cer_HS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
            diagnostic_ions_dict['CH4O2 loss'] = ref_precursor_mz - 48.02114
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['CH4O2 loss'] = ref_precursor_mz - 48.02114 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl+C2H3N'] = free_fa_dict['fa-2'] + (ex_mass['C']*2+ex_mass['H']*2+ex_mass['N']+ex_mass['H+']) - ex_mass['O']
            diagnostic_ions_dict['Sph-NCCO-3H'] = sphingo_base_dict['Sph_base'] -(ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']) - ex_mass['H+']

    if current_lipidclass == 'Cer_HDS' and adduct_type == '[M+H]+':
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'Cer_HDS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        pass

    if current_lipidclass == 'Hex2Cer' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        diagnostic_ions_dict['2Hex loss'] = diagnostic_ions_dict['C6H10O5 loss'] - 162.052833

    if current_lipidclass == 'Hex3Cer' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        diagnostic_ions_dict['2Hex loss'] = diagnostic_ions_dict['Hex loss'] - 162.052833
        diagnostic_ions_dict['3Hex loss'] = diagnostic_ions_dict['2Hex loss'] - 162.052833

    if current_lipidclass == 'HexCer_EOS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833
        diagnostic_ions_dict['Hex+H2O loss'] = ref_precursor_mz - 162.052833 - ex_mass['H2O']
        if valid_moiety_num >= 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'HexCer_EOS' and adduct_type == '[M+CH3COO]-':
        if valid_moiety_num == 3:
            diagnostic_ions_dict['N-acyl fa'] = free_fa_dict['fa-3']
            diagnostic_ions_dict['O-acyl loss'] = acyl_loss_dict['acyl-3 loss'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['O-acyl+Hex loss'] = acyl_loss_dict['acyl-3 loss'] - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['O-acyl amide'] = free_fa_dict['fa-2'] + ex_mass['N'] + ex_mass['H']

    if current_lipidclass == 'AHexCer' and adduct_type == '[M+CH3COO]-':
        if valid_moiety_num >= 2:
            diagnostic_ions_dict['hex-acyl fa'] = free_fa_dict['fa-1']
            diagnostic_ions_dict['acyl-1 loss'] = acyl_loss_dict['acyl-1 loss'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['acyl-1+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O'] - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
            diagnostic_ions_dict['acyl-1+Hex loss'] = acyl_loss_dict['acyl-1 loss'] - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        # diagnostic_ions_dict['Sph-C2H7NO loss'] = (ref_precursor_mz - (sphingo_base_dict['Sph_base'] - 
        #                                             (ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']))
        #                                              - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+'])
        # diagnostic_ions_dict['Sph-C2H7NO+H2O loss'] = (ref_precursor_mz - (sphingo_base_dict['Sph_base'] - 
        #                                             (ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O'])) -ex_mass['H2O']
        #                                              - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+'])

    if current_lipidclass == 'HexCer_AP' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833
        diagnostic_ions_dict['Hex+H2O loss'] = ref_precursor_mz - 162.052833 - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['PhytoSph'] = sphingo_base_dict['PhytoSph_base'] + ex_mass['H+']
            diagnostic_ions_dict['PhytoSph-H2O'] = sphingo_base_dict['PhytoSph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['PhytoSph-3H2O'] = sphingo_base_dict['PhytoSph-H2O'] - ex_mass['H2O']*2 + ex_mass['H+']

    if current_lipidclass == 'HexCer_AP' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['fa+O']= free_fa_dict['fa-2'] + ex_mass['O']
            diagnostic_ions_dict['fa-CO-3H']= free_fa_dict['fa-2'] - (ex_mass['C']+ex_mass['O']+ex_mass['H']*2)
            diagnostic_ions_dict['fa+C3H5NO']= free_fa_dict['fa-2'] + (ex_mass['C']*3+ex_mass['H']*6+ex_mass['N']+ex_mass['O']) - ex_mass['H+']

    if current_lipidclass == 'HexCer_NDS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833 - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']

    if current_lipidclass == 'HexCer_NDS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl+C2H3N'] = free_fa_dict['fa-2'] + (ex_mass['C']*2+ex_mass['H']*2+ex_mass['N']+ex_mass['H+']) - ex_mass['O']
            diagnostic_ions_dict['acyl-3H'] = free_fa_dict['fa-2'] - (ex_mass['O']+ex_mass['H']*2)
            diagnostic_ions_dict['Sph-NCCO-3H'] = sphingo_base_dict['Sph_base'] -(ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']) - ex_mass['H+']

    if current_lipidclass == 'HexCer_NS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
        diagnostic_ions_dict['Hex+H2O loss'] = ref_precursor_mz - 162.052833 - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'HexCer_NS' and (adduct_type == '[M+CH3COO]-' or adduct_type == '[M-H]-'):
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl+C2H3N'] = free_fa_dict['fa-2'] + (ex_mass['C']*2+ex_mass['H']*2+ex_mass['N']+ex_mass['H+']) - ex_mass['O']
            diagnostic_ions_dict['acyl-3H'] = free_fa_dict['fa-2'] - (ex_mass['O']+ex_mass['H']*2)
            diagnostic_ions_dict['Sph-NCCO-3H'] = sphingo_base_dict['Sph_base'] -(ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']+ex_mass['O']) - ex_mass['H+']

    if current_lipidclass == 'HexCer_HS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
        diagnostic_ions_dict['Hex+H2O loss'] = ref_precursor_mz - 162.052833 - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'HexCer_HS' and adduct_type == '[M+CH3COO]-':
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']

    if current_lipidclass == 'HexCer_HDS' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['H2O loss'] = ref_precursor_mz - ex_mass['H2O']
        diagnostic_ions_dict['Hex+H2O loss'] = ref_precursor_mz - 162.052833 - ex_mass['H2O']
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'HexCer_HDS' and adduct_type == '[M+CH3COO]-':
        if adduct_type == '[M-H]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833
        elif adduct_type == '[M+CH3COO]-':
            diagnostic_ions_dict['Hex loss'] = ref_precursor_mz - 162.052833 - adduct_dict_neg['[M+CH3COO]-'] - ex_mass['H+']

    if current_lipidclass == 'SM' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C5H15NO4P+'] = 184.073320905
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']

    if current_lipidclass == 'SM' and adduct_type == '[M+Na]+':
        diagnostic_ions_dict['C3H9N loss'] = ref_precursor_mz - (ex_mass['C']*3 + ex_mass['H']*9 + ex_mass['N'])
        diagnostic_ions_dict['C5H14NO4P loss'] = ref_precursor_mz - 183.06604
        diagnostic_ions_dict['C5H16NO5P loss'] = ref_precursor_mz - 183.06604 -ex_mass['H2O'] - adduct_dict_pos['[M+Na]+'] +ex_mass['H+']

    if current_lipidclass == 'SM' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['CH3 loss'] = ref_precursor_mz - adduct_dict_neg['[M+CH3COO]-'] - 15.0234738
        diagnostic_ions_dict['C4H11NO4P-'] = 168.043117937
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl+CH3 loss'] = acyl_loss_dict['acyl-2 loss']- 15.0234738 - adduct_dict_neg['[M+CH3COO]-']

    if current_lipidclass == 'ASM' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['CH3 loss'] = ref_precursor_mz - adduct_dict_neg['[M+CH3COO]-'] - 15.0234738

    if current_lipidclass == 'SM+O' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C5H15NO4P+'] = 184.073320905

    if current_lipidclass == 'SM+O' and adduct_type == '[M+CH3COO]-':
        diagnostic_ions_dict['CH3 loss'] = ref_precursor_mz - adduct_dict_neg['[M+CH3COO]-'] - 15.0234738
        diagnostic_ions_dict['C4H11NO4P-'] = 168.043117937

    if current_lipidclass == 'PE-Cer' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C2H7NO4P-'] = 140.0118196
        if valid_moiety_num == 2:
            diagnostic_ions_dict['Sph-C2H7NO loss'] = (ref_precursor_mz - (sphingo_base_dict['Sph_base'] - 
                                                        (ex_mass['C']*2+ex_mass['H']*9+ex_mass['N']+ex_mass['O'])))

    if current_lipidclass == 'PE-Cer+O' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C2H7NO4P-'] = 140.0118196
        if valid_moiety_num == 2:
            diagnostic_ions_dict['PhytoSph-C2H7NO loss'] = (ref_precursor_mz - (sphingo_base_dict['PhytoSph_base'] - 
                                                        (ex_mass['C']*2+ex_mass['H']*9+ex_mass['N']+ex_mass['O'])))

    if current_lipidclass == 'PI-Cer+O' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['C6H13O9P loss'] = ref_precursor_mz - 260.029722

    if current_lipidclass == 'PI-Cer+O' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['C6H10O8P-'] = 241.011877388
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl+O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['O']

    if current_lipidclass == 'SL' and (adduct_type == '[M+NH4]+' or adduct_type == '[M+H]+'):
        if adduct_type == '[M+H]+' and valid_moiety_num == 2:
            diagnostic_ions_dict['acyl loss'] = acyl_loss_dict['acyl-2 loss']
            diagnostic_ions_dict['acyl+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O']
            diagnostic_ions_dict['Sph-C2H8N loss'] = (ref_precursor_mz - (sphingo_base_dict['Sph-H2O'] - 
                                                        (ex_mass['C']*2+ex_mass['H']*7+ex_mass['N'])))
            diagnostic_ions_dict['Sph-C2H8N+H2O loss'] = (ref_precursor_mz - (sphingo_base_dict['Sph-H2O'] - 
                                                        (ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']) + ex_mass['H2O']))
        elif adduct_type == '[M+NH4]+' and valid_moiety_num == 2:
            diagnostic_ions_dict['acyl loss'] = acyl_loss_dict['acyl-2 loss'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
            diagnostic_ions_dict['acyl+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
            diagnostic_ions_dict['Sph-C2H8N loss'] = (ref_precursor_mz - (sphingo_base_dict['Sph-H2O'] - 
                                                        (ex_mass['C']*2+ex_mass['H']*7+ex_mass['N'])) - 
                                                        adduct_dict_pos['[M+NH4]+'] + ex_mass['H+'])
            diagnostic_ions_dict['Sph-C2H8N+H2O loss'] = (ref_precursor_mz - (sphingo_base_dict['Sph-H2O'] - 
                                                        (ex_mass['C']*2+ex_mass['H']*7+ex_mass['N']) + ex_mass['H2O']) - 
                                                        adduct_dict_pos['[M+NH4]+'] + ex_mass['H+'])
        diagnostic_ions_dict['C2H6NO3S+'] = 124.0062904

    if current_lipidclass == 'SL' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['SO3-'] = 79.95736258
        diagnostic_ions_dict['HSO3-'] = 80.96409042
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl loss'] = acyl_loss_dict['acyl-2 loss']

    if current_lipidclass == 'SL+O' and (adduct_type == '[M+NH4]+' or adduct_type == '[M+H]+'):
        diagnostic_ions_dict['C2H6NO3S+'] = 124.0062904
        if adduct_type == '[M+H]+' and valid_moiety_num == 2:
            diagnostic_ions_dict['acyl+O+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['O'] - ex_mass['H2O']
            diagnostic_ions_dict['acyl+O']= free_fa_dict['fa-2'] - ex_mass['H'] + ex_mass['H+'] 
        elif adduct_type == '[M+NH4]+' and valid_moiety_num == 2:
            diagnostic_ions_dict['acyl+O+H2O loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['O'] - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
            diagnostic_ions_dict['acyl+O']= free_fa_dict['fa-2'] - ex_mass['H'] + ex_mass['H+'] 

    if current_lipidclass == 'SL+O' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['SO3-'] = 79.95736258
        if valid_moiety_num == 2:
            diagnostic_ions_dict['acyl loss'] = acyl_loss_dict['acyl-2 loss'] - ex_mass['O']

    if current_lipidclass == 'PhytoSph' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['PhytoSph-H2O'] = sphingo_base_dict['PhytoSph-H2O'] + ex_mass['H+']
        diagnostic_ions_dict['PhytoSph-2H2O'] = sphingo_base_dict['PhytoSph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
        diagnostic_ions_dict['PhytoSph-3H2O'] = sphingo_base_dict['PhytoSph-H2O'] - ex_mass['H2O']*2 + ex_mass['H+']
        diagnostic_ions_dict['PhytoSph-H2O-CH2O'] = sphingo_base_dict['PhytoSph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'DHSph' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
        diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
        diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'Sph' and adduct_type == '[M+H]+':
        diagnostic_ions_dict['Sph-H2O'] = sphingo_base_dict['Sph-H2O'] + ex_mass['H+']
        diagnostic_ions_dict['Sph-2H2O'] = sphingo_base_dict['Sph-H2O'] - ex_mass['H2O'] + ex_mass['H+']
        diagnostic_ions_dict['Sph-H2O-CH2O'] = sphingo_base_dict['Sph-H2O'] - (ex_mass['C']+ex_mass['H']*2+ex_mass['O']) + ex_mass['H+']

    if current_lipidclass == 'DCAE' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['acyl+2H2O loss'] = acyl_loss_dict['acyl-1 loss'] - 2*ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']
    
    if current_lipidclass == 'DCAE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['acyl+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O']
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'GDCAE' and adduct_type == '[M+NH4]+':
        pass
    if current_lipidclass == 'GDCAE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'GLCAE' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['acyl+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']

    if current_lipidclass == 'GLCAE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'TDCAE' and adduct_type == '[M+NH4]+':
        pass
    if current_lipidclass == 'TDCAE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'TLCAE' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['acyl+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']

    if current_lipidclass == 'TLCAE' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'SSulfate' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['HSO4-'] = 96.96010327

    if current_lipidclass == 'AHexCS' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['Sterol'] = 369.351578

    if current_lipidclass == 'AHexCS' and (adduct_type == '[M-H]-' or adduct_type == '[M+CH3COO]-'):
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']
    
    if current_lipidclass == 'AHexBRS' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['BR-Sterol'] = 381.3515764

    if current_lipidclass == 'AHexBRS' and (adduct_type == '[M-H]-' or adduct_type == '[M+CH3COO]-'):
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'AHexCAS' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['CA-Sterol'] = 383.3672264

    if current_lipidclass == 'AHexCAS' and (adduct_type == '[M-H]-' or adduct_type == '[M+CH3COO]-'):
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'AHexSIS' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['SI-Sterol'] = 397.3828764

    if current_lipidclass == 'AHexSIS' and (adduct_type == '[M-H]-' or adduct_type == '[M+CH3COO]-'):
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'AHexSTS' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['ST-Sterol'] = 395.3672264

    if current_lipidclass == 'AHexSTS' and (adduct_type == '[M-H]-' or adduct_type == '[M+CH3COO]-'):
        diagnostic_ions_dict['fa']= free_fa_dict['fa-1']

    if current_lipidclass == 'Vitamine D' and (adduct_type == '[M+H]+' or adduct_type == '[M+Na]+'):
        pass

    if current_lipidclass == 'Cholesterol' and adduct_type == '[M-H2O+H]-':
        diagnostic_ions_dict['Sterol'] = 369.351578

    if current_lipidclass == 'SHex' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['acyl+H2O loss'] = acyl_loss_dict['acyl-1 loss'] - ex_mass['H2O'] - adduct_dict_pos['[M+NH4]+'] + ex_mass['H+']

    if current_lipidclass == 'SHex' and adduct_type == '[M-H]-':
        diagnostic_ions_dict['Hexose'] = 179.0561136

    if current_lipidclass == 'BRSE' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['BR-Sterol'] = 381.3515764

    if current_lipidclass == 'CASE' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['CA-Sterol'] = 383.3672264

    if current_lipidclass == 'CE' and (adduct_type == '[M+NH4]+' or adduct_type == '[M+Na]+'):
        diagnostic_ions_dict['Sterol'] = 369.351578

    if current_lipidclass == 'SISE' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['SI-Sterol'] = 395.3672264

    if current_lipidclass == 'STSE' and adduct_type == '[M+NH4]+':
        diagnostic_ions_dict['ST-Sterol'] = 395.3672264
    #endregion
    
    # diagnostic_ions_dict = {key: math_floor(v, 4) for key, v in diagnostic_ions_dict.items()}
    tuples_in_list = sorted(diagnostic_ions_dict.items(), key=lambda x: x[1], reverse=True)
    diagnostic_ions_dict = {t[0]: math_floor(t[1], 4) for t in tuples_in_list}
    cid_fragmentation_dict = {}
    ms_tolerance = math_floor(ms_tolerance_ppm*ref_precursor_mz/(1000*1000), 4)
    lipidsubclass_ions = []
    moiety_ions = []
    for key in diagnostic_ions_dict:
        if ('acyl' in key) or ('fa'in key) or ('Sph' in key) or ('ester' in key):
            moiety_ions.append(key)
        else:
            lipidsubclass_ions.append(key)
    cid_fragmentation_dict['Lipid subclass'] = categorize_ions(
        df=dataframe, diagnostic_ions=diagnostic_ions_dict, keys=lipidsubclass_ions)
    cid_fragmentation_dict['Moiety'] = categorize_ions(
        df=dataframe, diagnostic_ions=diagnostic_ions_dict, keys=moiety_ions)
    return cid_fragmentation_dict

def categorize_ions(df, diagnostic_ions, keys):
    """ Query diagnostic ions in measured OAD-MS/MS spectrum

    Args:
        df (Dataframe): Extracted OAD-MS/MS spectrum
        diagnostic_ions (flaot): m/z of reference ions
        keys (str): tag of CID fragment ion

    Returns:
        input_dict (dict): {tag: [ref_mz, measured_mz, intensity, ratio, ppm]}
    """
    input_dict = {}
    presence_count = 0
    if keys:
        for key in keys:
            ref_mz = math_floor(diagnostic_ions[key], 4)
            ref_front = ref_mz - 0.01
            ref_tail = ref_mz + 0.01
            extract_df = df[(df['frag m/z']>=ref_front) & (df['frag m/z']<=ref_tail)]
            extract_df = extract_df.reset_index(drop=True)
            if len(extract_df) == 1:
                measured_mz = extract_df['frag m/z'][0]
                intensity = extract_df['intensity'][0]
                ratio = extract_df['Ratio(%)'][0]
                ppm = (measured_mz-ref_mz)/ref_mz*1000*1000
                ppm = math_floor(ppm, 2)
                input_dict[key] = [ref_mz, measured_mz, intensity, ratio, ppm]
                presence_count += 1
            elif len(extract_df) > 1:
                sorted_df = extract_df.sort_values('intensity', ascending=False).reset_index(drop=True)
                measured_mz = sorted_df['frag m/z'][0]
                intensity = extract_df['intensity'][0]
                ratio = extract_df['Ratio(%)'][0]
                ppm = (measured_mz-ref_mz)/ref_mz*1000*1000
                ppm = math_floor(ppm, 2)
                input_dict[key] = [ref_mz, measured_mz, intensity, ratio, ppm]
                presence_count += 1
            else:
                input_dict[key] = [ref_mz, 0, 0, 0, 0]
        input_dict['Presence'] = math_floor((presence_count/len(keys)*100), 4)
    return input_dict    
#endregion


 