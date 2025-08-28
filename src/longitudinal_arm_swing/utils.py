import math
import pandas as pd
import re
import sys
import warnings

from collections import defaultdict
from pathlib import Path
from scipy.stats import wilcoxon, ranksums, norm
from tabulate import tabulate

sys.path.append(str(Path("src/").resolve()))

from longitudinal_arm_swing.constants import *

def defaultdict_to_dict(d):
    if isinstance(d, defaultdict):
        d = {k: defaultdict_to_dict(v) for k, v in d.items()}
    elif isinstance(d, dict):
        d = {k: defaultdict_to_dict(v) for k, v in d.items()}
    return d


def split_into_two_lines(text):
    words = text.split()
    mid = len(words) // 2
    return ' '.join(words[:mid]) + '\n' + ' '.join(words[mid:])


def extract_med_info(visit_data, med_colname):
    """Allocate participants into groups based on medication status."""
    return {
        'med': [id for id in visit_data.loc[visit_data[med_colname].isin([1, '1']), 'id'].tolist()],
        'no_med': [id for id in visit_data.loc[visit_data[med_colname].isin([0, '0']), 'id'].tolist()]
    }

def get_med_ids(med_info_ids, group):
    """Extract the participants who have started medication before the study."""
    return med_info_ids['visits'][1][group]['med']

def get_no_med_ids(med_info_ids, group, start_med_week_df):
    """Extract the participants who have never started medication."""
    return list(set([
        x for x in med_info_ids['visits'][1][group]['no_med'] 
        if (
            x in med_info_ids['visits'][2][group]['no_med'] and 
            x in med_info_ids['visits'][3][group]['no_med']
        ) or (
            x not in start_med_week_df['ID'].tolist() and (
                x not in med_info_ids['visits'][2][group]['med'] or 
                x not in med_info_ids['visits'][3][group]['med']
            )
        ) or (
            x in start_med_week_df['ID'].tolist() and start_med_week_df.loc[start_med_week_df['ID'] == x, 'StartWeek'].values[0] > 104  # Starting after study
        )]))

def get_start_med_ids(med_info_ids, group, start_med_week_df):
    """Extract the participants who started levodopa medication between visit 1 and visit 3."""
    return list(set([
        x for x in med_info_ids['visits'][1][group]['no_med'] 
        if (
            x in med_info_ids['visits'][2][group]['med'] or 
            x in med_info_ids['visits'][3][group]['med'] 
        ) or (
            x in start_med_week_df['ID'].tolist() and start_med_week_df.loc[start_med_week_df['ID'] == x, 'StartWeek'].values[0] <= 104  # Starting during study
        )
    ]))

def set_data_types(df, numeric_cols):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        
        # Replace USER_MISSING with NaN
        df = df.replace(to_replace=r'.*USER_MISSING*.', value=np.nan, regex=True)

        # Convert numeric columns
        df.loc[:, numeric_cols] = df.loc[:, numeric_cols].apply(pd.to_numeric, errors='coerce')

        return df

def add_updrs_columns(visit_df, dataset):
    # Create masks for 'right' and 'left'
    right_mask = visit_df['WatchSide'] == 'right'
    left_mask = visit_df['WatchSide'] == 'left'

    # Initialize Series for new scores
    rigidity_off_ws = pd.Series(index=visit_df.index, dtype=int)
    rigidity_off_non_ws = pd.Series(index=visit_df.index, dtype=int)
    bradykinesia_off_ws = pd.Series(index=visit_df.index, dtype=int)
    bradykinesia_off_non_ws = pd.Series(index=visit_df.index, dtype=int)
    hypokinesia_upper_off_ws = pd.Series(index=visit_df.index, dtype=int)
    hypokinesia_upper_off_non_ws = pd.Series(index=visit_df.index, dtype=int)

    # Compute score for watch and non-watch side
    rigidity_off_ws[right_mask] = visit_df.loc[right_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['rigidity']['right_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    bradykinesia_off_ws[right_mask] = visit_df.loc[right_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['bradykinesia']['watch_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    rigidity_off_non_ws[right_mask] = visit_df.loc[right_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['rigidity']['left_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)   
    bradykinesia_off_non_ws[right_mask] = visit_df.loc[right_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['bradykinesia']['non_watch_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    hypokinesia_upper_off_ws[right_mask] = visit_df.loc[right_mask, ['Up3OfRigRue', 'Up3OfFiTaYesDev', 'Up3OfProSYesDev', 'Up3OfHaMoYesDev']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    hypokinesia_upper_off_non_ws[right_mask] = visit_df.loc[right_mask, ['Up3OfRigLue', 'Up3OfFiTaNonDev', 'Up3OfProSNonDev', 'Up3OfHaMoNonDev']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)

    # Compute for left watch side
    rigidity_off_ws[left_mask] = visit_df.loc[left_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['rigidity']['left_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    bradykinesia_off_ws[left_mask] = visit_df.loc[left_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['bradykinesia']['watch_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    rigidity_off_non_ws[left_mask] = visit_df.loc[left_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['rigidity']['right_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    bradykinesia_off_non_ws[left_mask] = visit_df.loc[left_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['bradykinesia']['non_watch_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    hypokinesia_upper_off_ws[left_mask] = visit_df.loc[left_mask, ['Up3OfRigLue', 'Up3OfFiTaYesDev', 'Up3OfProSYesDev', 'Up3OfHaMoYesDev']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    hypokinesia_upper_off_non_ws[left_mask] = visit_df.loc[left_mask, ['Up3OfRigRue', 'Up3OfFiTaNonDev', 'Up3OfProSNonDev', 'Up3OfHaMoNonDev']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    
    # Compute all new columns at once
    new_columns = {
        f'updrs_3_rigidity_off_ws': rigidity_off_ws,
        f'updrs_3_rigidity_off_non_ws': rigidity_off_non_ws,
        f'updrs_3_bradykinesia_off_ws': bradykinesia_off_ws,
        f'updrs_3_bradykinesia_off_non_ws': bradykinesia_off_non_ws,
        f'updrs_3_hypokinesia_off_ws': rigidity_off_ws + bradykinesia_off_ws,
        f'updrs_3_hypokinesia_off_non_ws': rigidity_off_non_ws + bradykinesia_off_non_ws,
        'updrs_3_hypokinesia_upper_off_ws': hypokinesia_upper_off_ws,
        'updrs_3_hypokinesia_upper_off_non_ws': hypokinesia_upper_off_non_ws
    }

    if dataset == 'ppp': 
        rigidity_on_ws = pd.Series(index=visit_df.index, dtype=int)
        rigidity_on_non_ws = pd.Series(index=visit_df.index, dtype=int)
        bradykinesia_on_ws = pd.Series(index=visit_df.index, dtype=int)
        bradykinesia_on_non_ws = pd.Series(index=visit_df.index, dtype=int)
        hypokinesia_upper_on_ws = pd.Series(index=visit_df.index, dtype=int)
        hypokinesia_upper_on_non_ws = pd.Series(index=visit_df.index, dtype=int)

        rigidity_on_ws[right_mask] = visit_df.loc[right_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['rigidity']['right_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        rigidity_on_non_ws[right_mask] = visit_df.loc[right_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['rigidity']['left_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        bradykinesia_on_ws[right_mask] = visit_df.loc[right_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['bradykinesia']['watch_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        bradykinesia_on_non_ws[right_mask] = visit_df.loc[right_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['bradykinesia']['non_watch_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)  
        hypokinesia_upper_on_ws[right_mask] = visit_df.loc[right_mask, ['Up3OnRigRue', 'Up3OnFiTaYesDev', 'Up3OnProSYesDev', 'Up3OnHaMoYesDev']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        hypokinesia_upper_on_non_ws[right_mask] = visit_df.loc[right_mask, ['Up3OnRigLue', 'Up3OnFiTaNonDev', 'Up3OnProSNonDev', 'Up3OnHaMoNonDev']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)       

        rigidity_on_ws[left_mask] = visit_df.loc[left_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['rigidity']['left_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        rigidity_on_non_ws[left_mask] = visit_df.loc[left_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['rigidity']['right_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        bradykinesia_on_ws[left_mask] = visit_df.loc[left_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['bradykinesia']['watch_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        bradykinesia_on_non_ws[left_mask] = visit_df.loc[left_mask, UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['bradykinesia']['non_watch_side']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        hypokinesia_upper_on_ws[left_mask] = visit_df.loc[left_mask, ['Up3OnRigLue', 'Up3OnFiTaYesDev', 'Up3OnProSYesDev', 'Up3OnHaMoYesDev']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
        hypokinesia_upper_on_non_ws[left_mask] = visit_df.loc[left_mask, ['Up3OnRigRue', 'Up3OnFiTaNonDev', 'Up3OnProSNonDev', 'Up3OnHaMoNonDev']].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)

        new_columns['updrs_3_rigidity_on_ws'] = rigidity_on_ws
        new_columns['updrs_3_rigidity_on_non_ws'] = rigidity_on_non_ws
        new_columns['updrs_3_bradykinesia_on_ws'] = bradykinesia_on_ws
        new_columns['updrs_3_bradykinesia_on_non_ws'] = bradykinesia_on_non_ws
        new_columns['updrs_3_hypokinesia_on_ws'] = rigidity_on_ws + bradykinesia_on_ws
        new_columns['updrs_3_hypokinesia_on_non_ws'] = rigidity_on_non_ws + bradykinesia_on_non_ws
        new_columns['updrs_3_hypokinesia_upper_on_ws'] = hypokinesia_upper_on_ws
        new_columns['updrs_3_hypokinesia_upper_on_non_ws'] = hypokinesia_upper_on_non_ws

    new_columns['updrs_1_total'] = visit_df[UPDRS_PART_1_COLS].apply(
        lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    new_columns['updrs_2_total'] = visit_df[UPDRS_PART_2_COLS].apply(
        lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    new_columns['updrs_3_off_total'] = visit_df[UPDRS_PART_3_OFF_COLS].apply(
        lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    
    if dataset == 'ppp':
        new_columns['updrs_3_on_total'] = visit_df[UPDRS_PART_3_ON_COLS].apply(
            lambda row: row.sum() if row.isna().sum() == 0 else float('nan'), axis=1)
    else:
        new_columns['updrs_3_on_total'] = np.nan

    return pd.concat([visit_df, pd.DataFrame(new_columns)], axis=1)

def determine_affected_side_watch_side_mapping(df, score_category):
    mas_mask = df[f'updrs_3_{score_category}_off_ws'] > df[f'updrs_3_{score_category}_off_non_ws']
    las_mask = df[f'updrs_3_{score_category}_off_ws'] < df[f'updrs_3_{score_category}_off_non_ws']
    unknown_mask = df[f'updrs_3_{score_category}_off_ws'] == df[f'updrs_3_{score_category}_off_non_ws']

    mas_ids = list(df[mas_mask]['id'].values)
    las_ids = list(df[las_mask]['id'].values)
    rest_ids = list(df[unknown_mask]['id'].values)

    side_ids = {
        'mas': mas_ids,
        'las': las_ids,
        'rest': rest_ids,
    }

    mas_ids_subjective = [
        subject for subject in side_ids['rest'] 
        if (
            df.loc[df['id']==subject, 'MostAffSide'].values[0] in AFFECTED_SIDE_CATEGORIES['right'] and 
            df.loc[df['id']==subject, 'WatchSide'].values[0]=='right'
        ) or (
            df.loc[df['id']==subject,'MostAffSide'].values[0] in AFFECTED_SIDE_CATEGORIES['left'] and 
            df.loc[df['id']==subject, 'WatchSide'].values[0]=='left'
        )
    ]

    las_ids_subjective = [
        subject for subject in side_ids['rest'] 
        if (
            df.loc[df['id']==subject,'MostAffSide'].values[0] in AFFECTED_SIDE_CATEGORIES['right'] and 
            df.loc[df['id']==subject, 'WatchSide'].values[0]=='left'
        ) or (
            df.loc[df['id']==subject,'MostAffSide'].values[0] in AFFECTED_SIDE_CATEGORIES['left'] and 
            df.loc[df['id']==subject, 'WatchSide'].values[0]=='right'
        )
    ]

    rest_ids = [
        subject for subject in side_ids['rest'] 
        if subject not in mas_ids_subjective and subject not in las_ids_subjective
    ]

    side_ids['mas'] = side_ids['mas'] + mas_ids_subjective
    side_ids['las'] = side_ids['las'] + las_ids_subjective
    side_ids['rest'] = rest_ids

    return side_ids

def generate_ids_bins(df, dataset, start_ids, watch_side_dict, dysk_start_exclude):

    # Identify participants with unknown or switched watch sides
    watch_side_unknown, watch_side_switched = [], []
    if dataset in PD_DATASETS:
        watch_side_switched_per_visit = {'1-2': [], '2-3': []}
    else:
        watch_side_switched_per_visit = {'1-2': []}
    for subject in watch_side_dict:
        watch_side_subject = watch_side_dict[subject]['final']
        if pd.isna(watch_side_subject):
            watch_side_unknown.append(subject)
        elif watch_side_subject == 'switched':
            watch_side_switched.append(subject)

        if watch_side_dict[subject]['1'] in ['left', 'right'] and watch_side_dict[subject]['2'] in ['left', 'right'] and watch_side_dict[subject]['1'] != watch_side_dict[subject]['2']:
            watch_side_switched_per_visit['1-2'].append(subject)
        elif dataset in PD_DATASETS and watch_side_dict[subject]['1'] in ['left', 'right'] and pd.isna(watch_side_dict[subject]['2']) and watch_side_dict[subject]['3'] in ['left', 'right'] and watch_side_dict[subject]['1'] != watch_side_dict[subject]['3']:
            watch_side_switched_per_visit['1-2'].append(subject)

        if dataset in PD_DATASETS:
            if watch_side_dict[subject]['2'] in ['left', 'right'] and watch_side_dict[subject]['3'] in ['left', 'right'] and watch_side_dict[subject]['2'] != watch_side_dict[subject]['3']:
                watch_side_switched_per_visit['2-3'].append(subject)
            elif watch_side_dict[subject]['3'] in ['left', 'right'] and pd.isna(watch_side_dict[subject]['2']) and watch_side_dict[subject]['1'] in ['left', 'right'] and watch_side_dict[subject]['1'] != watch_side_dict[subject]['3']:
                watch_side_switched_per_visit['2-3'].append(subject)

    # Initialize dictionary to store per exclusion criteria the participants
    ids_bins = {}
    for visit_nr in VISITS_PER_DATASET[dataset]:
        df_visit = df.loc[df['visit'] == visit_nr]
        visit_ids = df_visit['id'].unique().tolist()
        walking_aid_ids = df_visit.loc[df_visit['walking_aid'] == True, 'id'].unique().tolist()
        dysk_ids = df_visit.loc[df_visit[f'at_least_{dysk_start_exclude}_dyskinesia'] == True, 'id'].unique().tolist()

        # Bin participants by clinical conditions
        ids_bins[f'no_clinical_data_visit_{visit_nr}'] = [
            x for x in start_ids if x.replace('_no_med', '').replace('_med', '') not in visit_ids 
        ]
        ids_bins[f'walking_aid_visit_{visit_nr}'] = walking_aid_ids
        ids_bins[f'at_least_{dysk_start_exclude}_dyskinesia_visit_{visit_nr}'] = dysk_ids

    # Add watch side exclusions
    ids_bins['watch_side_unknown'] = watch_side_unknown
    ids_bins['watch_side_switched'] = watch_side_switched
    ids_bins['watch_side_switched_visits_1_2'] = watch_side_switched_per_visit['1-2']

    if dataset in PD_DATASETS:
        ids_bins['watch_side_switched_visits_2_3'] = watch_side_switched_per_visit['2-3']

    return ids_bins

def determine_excluded_ids_med_split(start_ids, med_info_ids, dataset, base_exclude_by_category_dict, start_med=False):
        excluded_ids_split_by_med = {}
        if dataset != 'controls':  
            med_groups = ['med', 'no_med']     
            excluded_ids = []
            ids_remaining = []
            for med_status in med_groups:
                med_status_ids = med_info_ids['groups'][med_status][dataset]
                lng_ids_med_status = [id for id in start_ids if id in med_status_ids]

                if start_med and med_status == 'med':
                    base_exclude_by_category_dict_med_state = base_exclude_by_category_dict['post_med'].copy()
                elif start_med and med_status == 'no_med':
                    base_exclude_by_category_dict_med_state = base_exclude_by_category_dict['pre_med'].copy()
                else:
                    base_exclude_by_category_dict_med_state = base_exclude_by_category_dict

                excluded_ids_split_by_med[med_status] = {
                    k: sorted([
                        subject for subject in v 
                        if subject in lng_ids_med_status
                    ]) 
                    for k, v in base_exclude_by_category_dict_med_state.items()
                }

                # Flatten the excluded IDs and add to the total excluded list
                excluded_ids.extend([y for z in excluded_ids_split_by_med[med_status].values() for y in z])

                # Filter remaining IDs that are not excluded
                ids_remaining.extend([
                    x for x in lng_ids_med_status if 
                    all(x not in v for v in excluded_ids_split_by_med[med_status].values())
                ])
            
            ids_remaining = sorted(set(ids_remaining))

        else:
            excluded_ids_split_by_med = base_exclude_by_category_dict
            excluded_ids = [y for z in excluded_ids_split_by_med.values() for y in z]

            ids_remaining = sorted([x for x in start_ids if x not in excluded_ids])

        return excluded_ids_split_by_med, ids_remaining


def allocate_start_meds(med_info_ids):
    for dataset in PD_DATASETS:
        for subject in med_info_ids['groups']['start_med'][dataset]:
            med_info_ids['groups']['med'][dataset].append(subject)
            med_info_ids['groups']['no_med'][dataset].append(subject)

    return med_info_ids

def strip_med_suffix_id(id):
    return re.sub(r'_(no_)?med$', '', id)

def calculate_pvalue(x, y, test):
    if type(x) == np.ndarray:
        x = x.tolist()

    if type(y) == np.ndarray:
        y = y.tolist()
        
    if test == 'wilcoxon':
        return wilcoxon(x, y)
    elif test == 'ranksums':
        return ranksums(x, y)
    else:
        raise ValueError('Invalid test')
    
def plot_significance(ax, x_min, x_max, pvalue, y_min_significance, text_size, color='k'):
    if pvalue < 0.0001:
        asterisks = '****'
    elif pvalue < 0.001:
        asterisks = '***'
    elif pvalue < 0.01:
        asterisks = '**'
    elif pvalue < 0.05:
        asterisks = '*'
    else:
        asterisks = 'ns'
    
    bottom, top = ax.get_ylim()
    y_range = top - bottom

    bar_height = y_min_significance
    bar_tips = bar_height - y_range * 0.01
    text_height = bar_height + y_range * 0.01

    ax.plot([x_min, x_min, x_max, x_max], [bar_tips, bar_height, bar_height, bar_tips], lw=1, c=color)
    ax.text((x_max + x_min)/2, text_height, asterisks, ha='center', va='bottom', c=color, size=text_size)

def make_json_serializable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.generic, np.number)):
        return obj.item()
    elif isinstance(obj, dict):
        return {k: make_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [make_json_serializable(i) for i in obj]
    else:
        return obj
    
def calculate_ci(x, corr, ci=0.95):
    if ci == 0.95:
        multiplier = 1.96
    else:
        raise ValueError(f"Currently, only ci == 0.95 is implemented.")
    
    std_error = 1.0 / math.sqrt(len(x) - 3)
    delta = multiplier * std_error
    lower = math.tanh(math.atanh(corr) - delta)
    upper = math.tanh(math.atanh(corr) + delta)

    return lower, upper


def bootstrap_ci(data, n_boot=1000, ci=95):
    boot_means = []
    for _ in range(n_boot):
        sample = np.random.choice(data, size=len(data), replace=True)
        boot_means.append(np.nanmean(sample))
    lower = np.nanpercentile(boot_means, (100 - ci) / 2)
    upper = np.nanpercentile(boot_means, 100 - (100 - ci) / 2)
    return lower, upper


def bootstrap_ci_diffs(data1, data2, n_boot=1000, ci=95):
    diffs = []
    for _ in range(n_boot):
        sample1 = np.random.choice(data1, size=len(data1), replace=True)
        sample2 = np.random.choice(data2, size=len(data2), replace=True)
        
        diffs.append(np.median(sample1) - np.median(sample2))

    lower = np.percentile(diffs, (100 - ci) / 2)
    upper = np.percentile(diffs, 100 - (100 - ci) / 2)
    return lower, upper


def bca_ci(data1, data2, n_boot=1000, ci=95):
    alpha = (100 - ci) / 100
    data1 = np.array(data1)
    data2 = np.array(data2)

    # Observed statistic
    stat_obs = np.median(data1) - np.median(data2)

    # Bootstrap resampling
    boot_stats = []
    for _ in range(n_boot):
        sample1 = np.random.choice(data1, size=len(data1), replace=True)
        sample2 = np.random.choice(data2, size=len(data2), replace=True)
        boot_stats.append(np.median(sample1) - np.median(sample2))
    boot_stats = np.array(boot_stats)

    # Bias correction (z0)
    z0 = norm.ppf((boot_stats < stat_obs).mean())

    # Jackknife resampling to compute acceleration (a)
    n1, n2 = len(data1), len(data2)
    jack_stats = []

    for i in range(n1):
        jack_sample1 = np.delete(data1, i)
        stat = np.median(jack_sample1) - np.median(data2)
        jack_stats.append(stat)

    for i in range(n2):
        jack_sample2 = np.delete(data2, i)
        stat = np.median(data1) - np.median(jack_sample2)
        jack_stats.append(stat)

    jack_stats = np.array(jack_stats)
    jack_mean = np.mean(jack_stats)
    numer = np.sum((jack_mean - jack_stats)**3)
    denom = 6.0 * (np.sum((jack_mean - jack_stats)**2))**1.5
    a = numer / denom if denom != 0 else 0.0

    # Compute adjusted percentiles
    z_alpha1 = norm.ppf(alpha / 2)
    z_alpha2 = norm.ppf(1 - alpha / 2)

    pct1 = norm.cdf(z0 + (z0 + z_alpha1) / (1 - a * (z0 + z_alpha1)))
    pct2 = norm.cdf(z0 + (z0 + z_alpha2) / (1 - a * (z0 + z_alpha2)))

    lower = np.percentile(boot_stats, 100 * pct1)
    upper = np.percentile(boot_stats, 100 * pct2)

    return lower, upper
