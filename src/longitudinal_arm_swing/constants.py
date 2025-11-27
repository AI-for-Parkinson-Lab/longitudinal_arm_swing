import numpy as np
import os
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

BASE_PATH = os.getenv('BASE_PATH')
PATH_PREPROCESSED_DATA = Path(BASE_PATH, 'preprocessed_data')
PATH_CLINICAL_DATA = Path(BASE_PATH, 'clinical')
PATH_IDS = Path(BASE_PATH, 'ids')
PATH_STATS = Path(BASE_PATH, 'stats')
PATH_FIGURES = Path(BASE_PATH, 'figures')

LEDD_PPP_FILENAME = 'ledd.csv'

WATCH_SIDE_MAPPING = {
    1: 'right', 2: 'left',
}

AFFECTED_SIDE_MAPPING = {
    1: 'unilateral: right', 2: 'unilateral: left', 3: 'bilateral: right more than left',
    4: 'bilateral: left more than right', 5: 'bilateral: left as much as right', 6: 'nor right nor left'
}

AFFECTED_SIDE_RENAMING = {
    'mas': 'Most affected side',
    'las': 'Least affected side',
}

MEASURES_RENAMING = {
    'range_of_motion': 'Range of motion',
    'peak_angular_velocity': 'Peak angular velocity',
}

MEASURES_ABBREV = {
    'range_of_motion': 'RoM',
    'peak_angular_velocity': 'PAV'
}

UPDRS_SCORE_RENAMING = {
    'updrs_3_hypokinesia_off_ws': 'MDS-UPDRS Part III\nSubscore OFF',
    'updrs_3_hypokinesia_on_ws': 'MDS-UPDRS Part III\nSubscore ON',
    'updrs_2_total': 'MDS-UPDRS Part II Total',
    'Up3OfGait': 'MDS-UPDRS Part III\nGait OFF',
    'Up3OnGait': 'MDS-UPDRS Part III\nGait ON'
}

SEGMENT_DURATION_RENAMING = {
    '0_10': r'$<$ 10 s',
    '0_20': r'$<$ 20 s',
    '10_20': '10-20 s',
    '20_inf': r'$\geq$ 20 s',
    '10_30': '10-30 s',
    '30_inf': r'$\geq$ 30 s',
    '30_60': '30-60 s',
    '60_120': '60-120 s',
    '120_inf': r'$\geq$ 120 s',
}

UPDRS_HYPOKINESIA_SIDE_MAPPING = {}
for med_stage in ['Of', 'On']:
    UPDRS_HYPOKINESIA_SIDE_MAPPING[med_stage] = {
        'rigidity': {
            'right_side': [f'Up3{med_stage}RigRue', f'Up3{med_stage}RigRle'],
            'left_side': [f'Up3{med_stage}RigLle', f'Up3{med_stage}RigLue'],
        },
        'bradykinesia': {
            'watch_side': [f'Up3{med_stage}LAgiYesDev', f'Up3{med_stage}FiTaYesDev', f'Up3{med_stage}ToTaYesDev', f'Up3{med_stage}ProSYesDev', f'Up3{med_stage}HaMoYesDev'],
            'non_watch_side': [f'Up3{med_stage}HaMoNonDev', f'Up3{med_stage}LAgiNonDev', f'Up3{med_stage}ToTaNonDev', f'Up3{med_stage}FiTaNonDev', f'Up3{med_stage}ProSNonDev']
        },
        'other': [f'Up3{med_stage}Gait', f'Up3{med_stage}Facial', f'Up3{med_stage}RigNec', f'Up3{med_stage}Speech', f'Up3{med_stage}Arise', f'Up3{med_stage}StaPos']
    }
    
ANALYSIS_RENAMING = {
    'cs': 'Cross-sectional',
    'l1tf': 'L1 trend filtering',
    'srm': 'Standardized Response Mean',
    'regr': 'Regression'
}

GROUP_RENAMING = {
    'pd_med': 'PD medicated',
    'pd_no_med': 'PD unmedicated',
    'controls': 'Controls',
    'no_med': 'Unmedicated',
    'start_med': 'Started medication',
    'med': 'Medicated',
}

AGGREGATION_RENAMING = {
    'median': 'Median',
    '95p': r'95$^{\mathrm{th}}$ percentile',
    'mean_cov': 'CoV',
    'median_cov': 'CoV',
    'mean_std': 'SD',
    'median_std': 'SD',
}

AFFECTED_SIDE_CATEGORIES = {
    'left': ['unilateral: left', 'bilateral: left more than right'],
    'right': ['unilateral: right', 'bilateral: right more than left']
}

GENERAL_COLS = ['id', 'WatchSide', 'PrefHand']
GENERAL_COLS_VISIT_1 = GENERAL_COLS + ['Age', 'Gender']
GENERAL_COLS_VISIT_2 = GENERAL_COLS

GENERAL_COLS_PD = GENERAL_COLS + [
    'MotComOffStateTime', 'MotComDailyOffState', 'MotComDailyAwake',
    'MotComDysKinTime', 'MotComDailyDysKin', 'Up3OfHoeYah'
]
GENERAL_COLS_PD_VISIT_1 = GENERAL_COLS_PD + ['Age', 'Gender', 'MostAffSide', 'DiagParkYear', 'DiagParkMonth']
GENERAL_COLS_PD_VISITS_23 = GENERAL_COLS_PD

GENERAL_COLS_CONTROLS_VISIT_1 = GENERAL_COLS_VISIT_1 + ['AssessYear', 'AssessMonth']
GENERAL_COLS_CONTROLS_VISIT_2 = GENERAL_COLS_VISIT_2
GENERAL_COLS_PPP_VISIT_1 = GENERAL_COLS_PD_VISIT_1 + ['AssessYear', 'AssessWeekNum', 'YearsSinceDiagFloat']
GENERAL_COLS_DENOVO_VISIT_1 = GENERAL_COLS_PD_VISIT_1 + ['AssessYear', 'AssessMonth', 'YearsSinceDiagFloat']

UPDRS_PART_1_COLS = [f'Updrs2It0{x}' for x in range(7,10)] + [f'Updrs2It{x}' for x in range(10, 14)]
UPDRS_PART_2_COLS = [f'Updrs2It{x}' for x in range(14,27) if x != 23]
UPDRS_PART_3_COLS = {
    'bradykinesia': {
        'off': UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['bradykinesia']['watch_side'] + UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['bradykinesia']['non_watch_side'],
        'on': UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['bradykinesia']['watch_side'] + UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['bradykinesia']['non_watch_side']
    },
    'rigidity': {
        'off': UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['rigidity']['right_side'] + UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['rigidity']['left_side'],
        'on': UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['rigidity']['right_side'] + UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['rigidity']['left_side']
    },
    'rest': {
        'off': UPDRS_HYPOKINESIA_SIDE_MAPPING['Of']['other'],
        'on': UPDRS_HYPOKINESIA_SIDE_MAPPING['On']['other']
    }
}
PDQ_MOBILITY_COLS = [f'Pdq39It{i:02}' for i in range(1, 11)]
PDQ_ADL_COLS = [f'Pdq39It{i:02}' for i in range(11, 17)]
PDQ_EMO_COLS = [f'Pdq39It{i:02}' for i in range(17, 23)]
PDQ_STIGMA_COLS = [f'Pdq39It{i:02}' for i in range(23, 27)]
PDQ_SOCSUP_COLS = [f'Pdq39It{i:02}' for i in range(27, 30) if i != 28] + ['Pdq39It28a', 'Pdq39It28b']
PDQ_COG_COLS = [f'Pdq39It{i:02}' for i in range(30, 34)]
PDQ_COMM_COLS = [f'Pdq39It{i:02}' for i in range(34, 37)]
PDQ_BODDIS_COLS = [f'Pdq39It{i:02}' for i in range(37, 40)]
PDQ_ALL_COLS = PDQ_MOBILITY_COLS + PDQ_ADL_COLS + PDQ_EMO_COLS + PDQ_STIGMA_COLS + PDQ_SOCSUP_COLS + PDQ_COG_COLS + PDQ_COMM_COLS + PDQ_BODDIS_COLS

PDQ_BRADY_RIG_COLS = ['Pdq39It11', 'Pdq39It12', 'Pdq39It13', 'Pdq39It14', 'Pdq39It15', 'Pdq39It16']
PDQ_BALANCE_COLS = ['Pdq39It09']
PDQ_GAIT_COLS = ['Pdq39It04', 'Pdq39It05', 'Pdq39It06']

UPDRS_PART_3_OFF_COLS = UPDRS_PART_3_COLS['bradykinesia']['off'] + UPDRS_PART_3_COLS['rigidity']['off'] + UPDRS_PART_3_COLS['rest']['off']
UPDRS_PART_3_ON_COLS = UPDRS_PART_3_COLS['bradykinesia']['on'] + UPDRS_PART_3_COLS['rigidity']['on'] + UPDRS_PART_3_COLS['rest']['on']

PPP_UPDRS_COLS = UPDRS_PART_1_COLS + UPDRS_PART_2_COLS + UPDRS_PART_3_OFF_COLS + UPDRS_PART_3_ON_COLS + PDQ_ALL_COLS + ['ParkinMedUser']
DENOVO_UPDRS_COLS = UPDRS_PART_1_COLS + UPDRS_PART_2_COLS + UPDRS_PART_3_OFF_COLS + PDQ_ALL_COLS + ['Up3OfParkMedic']

UPDRS_COLS_PER_DATASET = {
    'ppp': PPP_UPDRS_COLS,
    'denovo': DENOVO_UPDRS_COLS
}

NUMERIC_COLS_VISITS = ['PrefHand', 'MotComDailyOffState', 'MotComDailyAwake', 'MotComDailyDysKin', 'Up3OfHoeYah']
NUMERIC_COLS_VISIT_1 = NUMERIC_COLS_VISITS + ['Age', 'Gender', 'DiagParkYear', 'DiagParkMonth', 'AssessYear']
NUMERIC_COLS_VISIT_1_PER_DATASET = {
    'ppp': NUMERIC_COLS_VISIT_1 + ['AssessWeekNum', 'MonthSinceDiag'],
    'denovo': NUMERIC_COLS_VISIT_1 + ['AssessMonth']
}

STORE_CLINICAL_COLS_PER_DATASET = {
    'ppp': {
        1: GENERAL_COLS_PPP_VISIT_1 + PPP_UPDRS_COLS,
        2: GENERAL_COLS_PD_VISITS_23 + PPP_UPDRS_COLS,
        3: GENERAL_COLS_PD_VISITS_23 + PPP_UPDRS_COLS
    },
    'denovo':{
        1: GENERAL_COLS_DENOVO_VISIT_1 + DENOVO_UPDRS_COLS,
        2: GENERAL_COLS_PD_VISITS_23 + DENOVO_UPDRS_COLS,
        3: GENERAL_COLS_PD_VISITS_23 + DENOVO_UPDRS_COLS
    },
    'controls': {
        1: GENERAL_COLS_CONTROLS_VISIT_1,
        2: GENERAL_COLS_CONTROLS_VISIT_2
    }
}

PD_DATASETS = ['ppp', 'denovo']
DATASETS = PD_DATASETS + ['controls']
MED_STATUSES = ['med', 'no_med', 'start_med']
DATASETS_MEDS = [f"{x} {med_status}" for x in PD_DATASETS for med_status in MED_STATUSES]
DATASETS_MEDS.append('controls')

WEEKS_PER_DATASET = {
    'ppp': [1] + list(np.arange(2, 106, 2)),
    'denovo': [1] + list(np.arange(2, 106, 2)),
    'controls': [1] + list(np.arange(2, 54, 2)),
}

valid_week_edges = {
    'ppp': {
        'strict': {
            'start': 2,
            'end': 100
        },
        'loose': {
            'start': [2, 4, 6],
            'end': [96, 98, 100, 102, 104]
        }
    },
    'denovo': {
        'strict': {
            'start': 2,
            'end': 100
        },
        'loose': {
            'start': [2, 4, 6],
            'end': [96, 98, 100, 102, 104]
        }
    },
    'controls': {
        'strict': {
            'start': 2,
            'end': 52
        },
        'loose': {
            'start': [2, 4, 6],
            'end': [48, 50, 52]
        }
    }
}

VISITS_PER_DATASET = {'ppp': [1, 2, 3], 'denovo': [1, 2, 3], 'controls': [1, 2]}

AT_LEAST_SLIGHT_DYSKINESIA_VALS = ['1: Heel licht', '2: Licht', '3: Matig', '4: Ernstig']
AT_LEAST_SIGNIFICANT_DYSKINESIA_VALS = ['2: Licht', '3: Matig', '4: Ernstig']

PATHS_PER_DATASET = {}
for dataset in DATASETS:
    PATHS_PER_DATASET[dataset] = {
        'aggregations': PATH_PREPROCESSED_DATA / dataset / '4.aggregation',
        'clinical': PATH_CLINICAL_DATA / dataset,
        'ids': PATH_IDS / dataset
    }

VISIT_FILENAMES_PER_DATASET = {dataset: {visit_nr: f'visit_{visit_nr}.csv' for visit_nr in VISITS_PER_DATASET[dataset]} for dataset in DATASETS}

PROCESSED_IDS_FILENAME = 'processed_ids_per_week.json'
VISIT_WEEKS_FILENAME = 'visit_weeks.json'
PPP_START_MED_FILENAME = 'start_week_ppp.csv'
DENOVO_START_MED_WEEK_FILENAME = 'start_week_denovo.csv'

STATS_FILENAMES_PER_DATASET = {}
for dataset in DATASETS:
    STATS_FILENAMES_PER_DATASET[dataset] = f'{dataset}_participant_stats.json'

# IDs
# Participants with alternative diagnosis in study: remove from all analyses
IDS_ALTERNATIVE_DIAGNOSIS = ['POMU04E638F5EC3A95C0', 'POMU065F753B3E97FF42', 'POMU10F33F607BAF7942', 'POMU4861269DC8E779A6', 'POMU584ADDDA24B571C1', 
                             'POMU785C84D79EA74864', 'POMU7C7AE0CD7DBDF0E2', 'POMU8AAE6ADB5A1BE3C6', 'POMUA20B1A38DFC93AEB', 'POMUA2EA522320DEB7B4',
                             'POMUC27BD4A6AF046175', 'POMUCEA19060522828C9', 'POMUDC96BAA36834CD4C', 'POMUECAA325B2EB33E7B', 'POMUF2A9725453BD4C5B']

# Participants with LEDD missing in visit 1 but taking medication: remove from regression analysis
IDS_LEDD_MISSING_VISIT_1 = ['POMU066326B8F70E150E', 'POMU0A109E0D97672361', 'POMU2E891447BF53DB23', 'POMU428FEF5AA8B909DC', 'POMU5A02B606856FA12F',
                            'POMU6080FFB910C10DC3', 'POMU74E11AFD8EE3BB6E', 'POMU82D03EF523FB9F4D''POMUA249715D7C6FDE8F']

# Exclude for all analyses
IDS_START_MED_TIME_UNKNOWN = ['POMU2E891447BF53DB23', 'POMU74E11AFD8EE3BB6E', 'POMU40D5BDD6FE082EFC', 'POMUEE612759679A830D',
                              'POMU8612F20403B70AF7', 'POMU961541776845B387', 'POMUA64681A2AFA522D6', 'POMUF43C40289F727F8E']
IDS_MED_INFO_PARTICIPANT_MISSING = ['POMU066326B8F70E150E', 'POMU0A109E0D97672361', 'POMU428FEF5AA8B909DC', 
                                    'POMU6080FFB910C10DC3', 'POMUA249715D7C6FDE8F', 'POMUF43C40289F727F8E']

# Exclude for regression
IDS_USE_ANTICHOLINERGIC_MEDS = ['POMU5EFACB489706C45A', 'POMUB2C4A660F7DE131B']
IDS_LEDD_DOUBTFUL = ['POMU6CA359B35463F82A', 'POMUEC313C797213BF9C', 'POMUF6514A4F9B462AB1']
IDS_DIAGNOSIS_DATE_DOUBTFUL = ['POMU600C11F136E6FB4D', 'POMU26A3A790D7E1B9CA']

ANALYSIS_STEPS = ['cs', 'l1tf', 'srm', 'srm', 'regr']
CS_ANALYSIS_WEEKS = [1, 2]

# Colorblind-safe colors
COLOR_GROUP_1 = '#E69F00'  # Orange
COLOR_GROUP_2 = '#56B4E9'  # Sky blue
COLOR_GROUP_NEUTRAL = "#686868" # Dark gray, neutral

PLOT_TITLESIZE = 18
PLOT_LABELSIZE = 16
PLOT_TICKSIZE = 14
PLOT_LEGEND_TITLESIZE = 16
PLOT_LEGEND_FONTSIZE = 14
PLOT_TEXTSIZE = 14
