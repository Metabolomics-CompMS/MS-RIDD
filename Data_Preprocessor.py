import datetime
import math
import pandas as pd
import re


#region Lipid subclass information
fatty_aclys = ['FA', 'NAGly', 'NAGlySer', 'NAOrn', 'NAE', 'CAR', 'FAHFA']
glycerolipids = ['DG', 'EtherDG', 'DGDG', 'EtherDGDG', 'MGDG', 'EtherMGDG', 
                 'SQDG', 'EtherSMGDG' 'MG', 'ADGGA', 'DGCC', 'DGGA', 'DGTS/A', 
                 'LDGCC', 'LDGTS/A', 'EtherTG', 'TG', 'OxTG', 'FAHFATG']
glycerophospholipids = ['LPA', 'PA', 'EtherLPC', 'EtherPC', 'LPC', 'PC', 
                        'EtherLPE', 'EtherPE', 'EtherPE(P)', 'PlasmPE', 'LNAPE',
                        'LPE', 'PE', 'BMP', 'EtherLPG', 'EtherPG', 'HBMP', 
                        'LPG', 'PG', 'CL', 'DLCL', 'MLCL', 'Ac2PIM1', 'Ac2PIM2',
                        'Ac3PIM2', 'Ac4PIM2', 'EtherPI', 'LPI', 'PI', 'EtherPS',
                        'LNAPS', 'LPS', 'PS', 'PEtOH', 'PMeOH', 'EtherOxPE', 
                        'OxPC', 'OxPE', 'OxPG', 'OxPI', 'OxPS']
prenol_lipids = ['VAE', 'CoQ', 'Vitamin E', 'Vitamin_E']
saccharolipids = ['LipidA']
sphingolipids = ['GM3', 'SHexCer', 'SHexCer+O', 
                 'Cer_ADS', 'Cer_AP', 'Cer_AS', 'Cer_BDS', 'Cer_BS', 'Cer_HDS', 
                 'Cer_HS', 'Cer_EBDS', 'Cer_EODS', 'Cer_EOS', 'Cer_NDS', 
                 'Cer_NP', 'Cer_NS', 'CerP', 'AHexCer', 'HexCer_ADS', 
                 'HexCer_AP', 'HexCer_AS', 'HexCer_BDS', 'HexCer_BS', 
                 'HexCer_HDS', 'HexCer_HS', 'HexCer_EOS', 'HexCer_NDS', 
                 'HexCer_NP', 'HexCer_NS', 'Hex2Cer', 'Hex3Cer', 'ASM', 
                 'PE_Cer', 'PE_Cer+O', 'PI_Cer', 'SM', 'SM+O',
                 'PhytoSph', 'SL', 'SL+O', 'DHSph', 'Sph']
sterol_lipids = ['CASulfate', 'BileAcid', 
                 'DCAE', 'GDCAE', 'GLCAE', 'TDCAE', 'TLCAE', 
                 'AHexCAS', 'AHexCS', 'AHexSIS', 'AHexBRS', 'AHexSTS', 
                 'Vitamin D', 'Vitamin_D', 
                 'SSulfate', 'BRSE', 'CASE', 'CE', 'Cholesterol', 
                 'SHex', 'SISE', 'STSE', 'SPE', 'BAHex', 'BASulfate', 'SPEHex', 
                 'SPGHex', 'BRSLPHex', 'BRSPHex', 'CASLPHex', 'CASPHex', 
                 'SISLPHex', 'SISPHex', 'STSLPHex', 'STSPHex']
lipidclass_dict =  {'Fatty acyls': fatty_aclys,
                    'Glycerolipids': glycerolipids,
                    'Glycerophospholipids': glycerophospholipids,
                    'Prenol lipids': prenol_lipids,
                    'Saccharolipids': saccharolipids, 
                    'Sphingolipids': sphingolipids,
                    'Sterol lipids': sterol_lipids
                    }
#endregion
#region Acyl chains information
no_acyl_list = ['Others', 'CoQ', 'CA', 
                'Vitamin D', 'Vitamin_D', 'Vitamin E', 'Vitamin_E', 
                'Cholesterol', 'SSulfate', 'CASulfate', 'BASulfate', 
                'SHex', 'SPE', 'BAHex', 'SPEHex', 'SPGHex']
mono_acyl_list = ['FA', 'NAE', 'CAR', 
                  'MG', 'LDGCC', 'LDGTS/A', 
                  'LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPS', 
                  'EtherLPC', 'EtherLPE', 'EtherLPG', 
                  'PhytoSph', 'DHSph', 'Sph', 'VAE', 
                  'DCAE', 'GDCAE', 'GLCAE', 'TDCAE', 'TLCAE', 
                  'AHexCAS', 'AHexCS', 'AHexSIS', 'AHexBRS', 'AHexSTS', 
                  'CE', 'BRSE', 'CASE', 'SISE', 'STSE']
tri_acyl_list = ['ADGGA', 'TG', 'EtherTG', 'HBMP', 'MLCL', 
                 'Cer_EBDS', 'Cer_EODS', 'Cer_EOS', 
                 'AHexCer', 'HexCer_EOS', 'ASM']
triacyls_sphingolipids = ['ASM', 'AHexCer',  'HexCer_EOS'
                          'Cer_EBDS', 'Cer_EODS', 'Cer_EOS',]
#endregion


class DataPreprocessor(object):
    def __init__(self, path1, fmt1, path2, fmt2):
        dt_now = datetime.datetime.now()
        idx = str(dt_now).find('.') +1
        stamp = str(dt_now)[:idx].replace('-', '').replace(' ', '_').replace('.', 's')
        timestamp = stamp.replace(':', 'h', 1).replace(':', 'm', 1)
        directory = path2.rsplit('/', 1)[0]
        self.timestamp = timestamp
        self.directory = directory
        self.path1, self.path2 = path1, path2
        self.fmt1, self.fmt2 = fmt1, fmt2
        
    def merge_bipolarity_cid_data(self, output='pos'):
        neg_path, neg_format = self.path1, self.fmt1
        pos_path, pos_format = self.path2, self.fmt2
        neg = self.check_input_format(neg_path, neg_format)
        pos = self.check_input_format(pos_path, pos_format)
        neg = self.add_moiety_info(neg)
        pos = self.add_moiety_info(pos)
        neg, pos = self.complement_moiety_info(neg, pos)
        raw_neg_path = f'{self.directory}/{self.timestamp}_merged_neg_data.txt'
        raw_pos_path = f'{self.directory}/{self.timestamp}_merged_pos_data.txt'
        neg.to_csv(raw_neg_path, index=False, sep='\t')
        pos.to_csv(raw_pos_path, index=False, sep='\t')
        pos = self.extract_annotated_molecules(pos)
        neg = self.extract_annotated_molecules(neg)
        list_neg_path = f'{self.directory}/{self.timestamp}_annotation_list_neg.txt'
        list_pos_path = f'{self.directory}/{self.timestamp}_annotation_list_pos.txt'
        neg.to_csv(list_neg_path, index=False, sep='\t')
        pos.to_csv(list_pos_path, index=False, sep='\t')
        if output == 'pos':
            pos['CCS'] = 0
            txt_path = f'{self.directory}/{self.timestamp}_txt_library_pos.txt'
            pos = self.rename_for_txt_library(pos)
            pos.to_csv(txt_path, index=False, sep='\t',
                       columns=['Name', 'MZ', 'RT', 'Adduct', 'InChIKey',
                                'Formula', 'SMILES', 'Ontology', 'CCS'])
        elif output == 'neg':
            neg['CCS'] = 0
            txt_path = f'{self.directory}/{self.timestamp}_txt_library_neg.txt'
            neg = self.rename_for_txt_library(neg)
            neg.to_csv(txt_path, index=False, sep='\t',
                       columns=['Name', 'MZ', 'RT', 'Adduct', 'InChIKey',
                                'Formula', 'SMILES', 'Ontology', 'CCS'])

    def merge_cid_and_oad_data(self):
        cid_path, cid_format = self.path1, self.fmt1
        oad_path, oad_format = self.path2, self.fmt2
        cid_raw_table = self.check_input_format(cid_path, cid_format)
        oad = self.check_input_format(oad_path, oad_format)
        cid_table = self.extract_annotated_molecules(cid_raw_table)
        mass_tolerance = 0.01
        rt_tolerance = 0.5
        column_list = ['Metabolite name', 'Adduct type', 'Reference RT', 
                       'Reference m/z', 'Formula', 'Ontology', 'INCHIKEY', 
                       'SMILES']
        for row, df in cid_table.iterrows():
            front_mz = df['Average Mz'] - mass_tolerance
            tail_mz = df['Average Mz'] + mass_tolerance
            front_rt = df['Average Rt(min)'] - rt_tolerance
            tail_rt = df['Average Rt(min)'] + rt_tolerance
            extracted_oad_df = oad[
                (oad['Average Mz'] > front_mz)&(oad['Average Mz'] < tail_mz)
                &(oad['Average Rt(min)'] > front_rt)&(oad['Average Rt(min)'] < tail_rt)]
            df_len = len(extracted_oad_df)
            if df_len == 1:
                target_row = extracted_oad_df.index[0]
                for colomn in column_list:
                    oad.loc[target_row:target_row, colomn] = df[colomn]
            elif df_len >= 2:
                target_rows = list(extracted_oad_df.index)
                for target_row in target_rows:
                    for colomn in column_list:
                        oad.loc[target_row:target_row, colomn] = df[colomn]
        # oad = oad.rename(columns={'Alignment ID': 'ID', 
        #                           'Average Rt(min)': 'RT(min)', 
        #                           'Average Mz': 'Precursor m/z'})
        save_path = f'{self.directory}/{self.timestamp}_merged_OAD_data.txt'
        oad.to_csv(save_path, index=False, sep='\t')

    def rename_for_txt_library(self, df):
        df = df.rename(columns={'Metabolite name': 'Name',
                                'Reference m/z': 'MZ',
                                'Average Rt(min)': 'RT',
                                'Adduct type': 'Adduct',
                                'INCHIKEY': 'InChIKey'})
        return df

    def check_input_format(self, path, fmt):
        if fmt == 'Alignment':
            raw_table = pd.read_csv(path, skiprows=[0,1,2,3], sep='\t')
            raw_table = raw_table.rename(columns={'Alignment ID': 'ID'})
        elif fmt == 'PeakList':
            raw_table = pd.read_csv(path, sep='\t')
            raw_table = raw_table.rename(columns={'PeakID': 'ID', 
                                                  'Title': 'Metabolite name', 
                                                  'RT (min)': 'Average Rt(min)', 
                                                  'Precursor m/z': 'Average Mz', 
                                                  'Adduct': 'Adduct type', 
                                                  'InChIKey': 'INCHIKEY', 
                                                  'MSMS spectrum': 'MS/MS spectrum'})
        return raw_table
    
    def extract_annotated_molecules(self, df):
        ex_df = df.dropna(subset=['Metabolite name'])
        ex_df = ex_df[ex_df['Metabolite name'] != 'Unknown']
        ex_df = ex_df[~ex_df['Metabolite name'].str.startswith('w/o')]
        ex_df = ex_df[~ex_df['Metabolite name'].str.startswith('RIKEN')]
        ex_df = ex_df[~ex_df['Metabolite name'].str.startswith('Unsettled')]
        ex_df = ex_df.fillna({'Reference RT': 0})
        return ex_df

    def exclude_IS(self, df):
        df = df[~df['Metabolite name'].str.contains('\(d\d+\)')]
        return df

    def add_moiety_info(self, df):
        new_cols = ['chains solved', 'Brutto', 'Total chain', 'Total db', 
                    'chain-1', 'db-1', 'chain-2', 'db-2', 
                    'chain-3', 'db-3', 'chain-4', 'db-4']
        for col in new_cols:
            df[col] = 0 if col != 'chains solved' else False
        annotated_df = self.extract_annotated_molecules(df)
        names = list(annotated_df['Metabolite name'].values)
        ontologies = list(annotated_df['Ontology'].values)
        idxs, cols = list(annotated_df.index), list(annotated_df.columns)
        col_pos = [cols.index(new) for new in new_cols]
        def get_chain_and_db(moiety):
            chain_and_db = moiety.split(':')
            chain, db = int(chain_and_db[0]), int(chain_and_db[1])
            return chain, db
        new_values = []
        exception_cls = ['CL', 'AHexCer', 'ASM', 
                         'CerEBDS', 'Cer_EOS', 'HexCer_EOS']
        for name, ontology in zip(names, ontologies):
            solved, brutto, total_chain, total_db = False, 0, 0, 0
            chian_1, db_1, chian_2, db_2, chian_3, db_3, chian_4, db_4 \
                = 0, 0, 0, 0, 0, 0, 0, 0
            if '|' in name:
                solved = True
                split_name = name.split('|')
                brutto = re.findall(r'\d+\:\d+', split_name[0])[0]
                total_chain, total_db = get_chain_and_db(brutto)
                moieties = re.findall(r'\d+\:\d+', split_name[1])
                if len(moieties) == 2:
                    chian_1, db_1 = get_chain_and_db(moieties[0])
                    chian_2, db_2 = get_chain_and_db(moieties[1])
                    if ontology in exception_cls:
                        solved = False
                    if ontology == 'AHexCer':
                        #AHexCer 60:2;3O|AHexCer (O-18:1)42:1;3O
                        chian_1, chian_2 = chian_2, chian_1
                        db_1, db_2 = db_2, db_1
                elif len(moieties) == 3:
                    chian_1, db_1 = get_chain_and_db(moieties[0])
                    chian_2, db_2 = get_chain_and_db(moieties[1])
                    chian_3, db_3 = get_chain_and_db(moieties[2])
                elif len(moieties) == 4:
                    chian_1, db_1 = get_chain_and_db(moieties[0])
                    chian_2, db_2 = get_chain_and_db(moieties[1])
                    chian_3, db_3 = get_chain_and_db(moieties[2])
                    chian_4, db_4 = get_chain_and_db(moieties[3])
            else:
                if ontology in no_acyl_list:
                    pass
                elif ontology in mono_acyl_list:
                    solved = True
                    if ontology in sterol_lipids and '/' in name:
                        brutto = re.findall(r'\d+\:\d+', name)[1]
                    else:
                        brutto = re.findall(r'\d+\:\d+', name)[0]
                    total_chain, total_db = get_chain_and_db(brutto)
                    chian_1, db_1 = get_chain_and_db(brutto)
                elif '(d' in name:
                    solved = True
                    moieties = re.findall(r'\d+\:\d+', name)
                    if len(moieties) == 2:
                        chian_1, db_1 = get_chain_and_db(moieties[0])
                        chian_2, db_2 = get_chain_and_db(moieties[1])
                    if len(moieties) == 3:
                        chian_1, db_1 = get_chain_and_db(moieties[0])
                        chian_2, db_2 = get_chain_and_db(moieties[1])
                        chian_3, db_3 = get_chain_and_db(moieties[2])
                    total_chain = chian_1 + chian_2 + chian_3 + chian_4
                    total_db = db_1 + db_2 + db_3 + db_4
                    brutto = str(total_chain) + ':' + str(total_db)
                else:
                    brutto = re.findall(r'\d+\:\d+', name)[0]
                    total_chain, total_db = get_chain_and_db(brutto)
            all_info = [solved, brutto, total_chain, total_db, 
                        chian_1, db_1, chian_2, db_2, 
                        chian_3, db_3, chian_4, db_4]
            new_values.append(all_info)
            # df.loc[idx:idx, new_cols] = all_info
        df.iloc[idxs, col_pos] = new_values
        return df

    def complement_moiety_info(self, raw_neg, raw_pos):
        neg = self.extract_annotated_molecules(raw_neg)
        pos = self.extract_annotated_molecules(raw_pos)
        raw_neg['chains complement'], raw_pos['chains complement'] = '', ''
        raw_neg['comple delta RT'], raw_pos['comple delta RT'] = 0, 0
        neg, pos = self.exclude_IS(neg), self.exclude_IS(pos)
        # temporary change Brutto of PlasmPE
        pos = self.temp_change_brutto_of_plasm(pos)
        unsolved_neg = neg[(neg['chains solved'] == False)
                          |(neg['Ontology'] == 'EtherPE')]
        unsolved_pos = pos[pos['chains solved'] == False]
        comple_cols = ['Metabolite name', 'chains solved', 
                       'chain-1', 'db-1', 'chain-2', 'db-2',
                       'chain-3', 'db-3', 'chain-4', 'db-4',
                       'chains complement', 'comple delta RT']
        plasm_cols = ['Ontology', 'Brutto', 'Total db']
        def get_abs_delta(rt, f_rt):
            return self.math_floor((abs(rt-f_rt)), 3)
        #region Neg
        for row, one in unsolved_neg.iterrows():
            brutto = one['Brutto']
            ontology = one['Ontology']
            rt = one['Average Rt(min)']
            if ontology == 'EtherPE':
                find = pos[((pos['Ontology']==ontology)|(pos['Ontology']=='PlasmPE'))
                          &(pos['Brutto']==brutto)
                          &(pos['chains solved']==True)]
            else:
                find = pos[(pos['Ontology']==ontology)
                          &(pos['Brutto']==brutto)
                          &(pos['chains solved']==True)]
            if len(find) > 1:
                idxs = list(find.index)
                f_rts = list(find['Average Rt(min)'].values)
                rt_d = {idx: get_abs_delta(rt, f_rt)
                        for idx, f_rt in zip(idxs, f_rts)}
                cand_idx = sorted(rt_d.items(), key=lambda x:x[1])[0][0]
                cols = list(find.columns)
                find = find.loc[cand_idx:cand_idx, cols]
            if len(find) == 1:
                new_ont = find['Ontology'].values[0]
                if new_ont == 'PlasmPE':
                    total_chain = find['Total chain'].values[0]
                    total_db = find['Total db'].values[0]
                    brutto = str(total_chain) + ':' + str(total_db)
                    raw_neg.loc[row:row, plasm_cols] = new_ont, brutto, total_db
                name = find['Metabolite name']
                chain_1, db_1 = find['chain-1'], find['db-1']
                chain_2, db_2 = find['chain-2'], find['db-2']
                chain_3, db_3 = find['chain-3'], find['db-3']
                chain_4, db_4 = find['chain-4'], find['db-4']
                pos_id = find['ID']
                delta = get_abs_delta(rt, find['Average Rt(min)'])
                raw_neg.loc[row:row, comple_cols] = name, True, \
                    chain_1,db_1, chain_2,db_2, chain_3,db_3, chain_4,db_4, \
                    pos_id, delta
        #endregion
        #region Pos
        for row, one in unsolved_pos.iterrows():
            brutto = one['Brutto']
            ontology = one['Ontology']
            rt = one['Average Rt(min)']
            find = neg[(neg['Ontology']==ontology)
                      &(neg['Brutto']==brutto)
                      &(neg['chains solved']==True)]
            if len(find) > 1:
                idxs = list(find.index)
                f_rts = list(find['Average Rt(min)'].values)
                rt_d = {idx: get_abs_delta(rt, f_rt)
                        for idx, f_rt in zip(idxs, f_rts)}
                cand_idx = sorted(rt_d.items(), key=lambda x:x[1])[0][0]
                cols = list(find.columns)
                find = find.loc[cand_idx:cand_idx, cols]
            if len(find) == 1:
                name = find['Metabolite name']
                chain_1, db_1 = find['chain-1'], find['db-1']
                chain_2, db_2 = find['chain-2'], find['db-2']
                chain_3, db_3 = find['chain-3'], find['db-3']
                chain_4, db_4 = find['chain-4'], find['db-4']
                neg_id = find['ID']
                delta = get_abs_delta(rt, find['Average Rt(min)'])
                raw_pos.loc[row:row, comple_cols] = name, True, \
                    chain_1,db_1, chain_2,db_2, chain_3,db_3, chain_4,db_4, \
                    neg_id, delta
        #endregion
        #complement triacyl ceramides
        raw_neg, raw_pos = self.complement_triacyl_ceramides(
                            raw_neg=raw_neg, raw_pos=raw_pos)
        #return result
        return raw_neg, raw_pos

    def temp_change_brutto_of_plasm(self, df):
        plasm_df = df[df['Ontology'] == 'PlasmPE']
        idxs = list(plasm_df.index)
        for idx, chain, db in zip(
            idxs, plasm_df['Total chain'], plasm_df['Total db']):
            temp_brutto = str(chain) + ':' + str(db+1)
            df.loc[idx:idx, ['Brutto']] = temp_brutto
        return df

    def complement_triacyl_ceramides(self, raw_neg, raw_pos):
        neg = self.extract_annotated_molecules(raw_neg)
        pos = self.extract_annotated_molecules(raw_pos)
        unsolved_neg = neg[neg['chains solved'] == False]
        unsolved_pos = pos[pos['chains solved'] == False]
        f_neg = unsolved_neg[unsolved_neg['Metabolite name'].str.contains('\|')]
        f_pos = unsolved_pos[unsolved_pos['Metabolite name'].str.contains('\|')]
        cand_neg = f_neg[(f_neg['Ontology']=='Cer_EBDS')
                        |(f_neg['Ontology']=='Cer_EODS')
                        |(f_neg['Ontology']=='Cer_EOS')
                        |(f_neg['Ontology']=='HexCer_EOS')]
        cand_pos = f_pos[(f_pos['Ontology']=='Cer_EBDS')
                        |(f_pos['Ontology']=='Cer_EODS')
                        |(f_pos['Ontology']=='Cer_EOS')
                        |(f_pos['Ontology']=='HexCer_EOS')]
        comple_cols = ['Metabolite name', 'chains solved', 'chain-1', 'db-1', 
                       'chain-2', 'db-2', 'chain-3', 'db-3', 
                       'chains complement', 'comple delta RT']
        def get_abs_delta(rt, f_rt):
            return self.math_floor((abs(rt-f_rt)), 3)
        # Neg -> FA, Pos -> Sphingobase
        for row, one in cand_pos.iterrows():
            ontology = one['Ontology']
            brutto = one['Brutto']
            rt = one['Average Rt(min)']
            find = cand_neg[(cand_neg['Ontology'] == ontology)
                           &(cand_neg['Brutto'] == brutto)]
            if len(find) == 0:
                continue
            elif len(find) > 1:
                idxs = list(find.index)
                f_rts = list(find['Average Rt(min)'].values)
                rt_d = {idx: get_abs_delta(rt, f_rt)
                        for idx, f_rt in zip(idxs, f_rts)}
                cand_idx = sorted(rt_d.items(), key=lambda x:x[1])[0][0]
                cols = list(find.columns)
                find = find.loc[cand_idx:cand_idx, cols]
            # HexCer 67:4;4O|HexCer 16:1;2O/51:3;2O
            name = one['Metabolite name'].values[0]
            head_name, tail_name = name.split('|')[0], name.split('|')[1]
            cls_type = tail_name.split(' ', 1)[0] # HexCer
            fa_chain, fa_db = find['chain-2'].values[0], find['db-2'].values[0]
            sph_chain, sph_db = one['chain-1'].values[0], one['db-1'].values[0]
            unsol_c, unsol_db = one['chain-2'].values[0], one['db-2'].values[0]
            new_chain_2, new_db_2 = (unsol_c-fa_chain), (unsol_db-fa_db-1)
            name = '{head}|{cls} {sph}:{db1};2O/{n_acyl}:{db2}(FA {fa})'.format(
            head=head_name, cls=cls_type, sph=sph_chain, db1=sph_db,
            n_acyl=new_chain_2, bd2=new_db_2, fa=str(fa_chain)+':'+str(fa_db))
            f_rt = find['Average Rt(min)'].values[0]
            delta_rt = get_abs_delta(rt, f_rt)
            #input -> Pos
            id_neg = find['Alignment ID'].values[0]
            raw_pos.loc[row:row, comple_cols] = name, True, \
                sph_chain, sph_db, new_chain_2, new_db_2, fa_chain, fa_db, \
                id_neg, delta_rt
            #input -> Neg
            id_pos = one['Alignment ID'].values[0]
            idx_neg = list(find.index)[0]
            raw_neg.loc[idx_neg:idx_neg, comple_cols] = name, True, \
                sph_chain, sph_db, new_chain_2, new_db_2, fa_chain, fa_db, \
                id_pos, delta_rt
        return raw_neg, raw_pos

    def math_floor(self, num, digit):
        floored = math.floor(num*10**digit)/(10**digit)
        return floored