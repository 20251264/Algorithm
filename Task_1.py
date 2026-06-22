{\rtf1\ansi\ansicpg949\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red0\green0\blue0;}
{\*\expandedcolortbl;;\csgray\c0;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww13140\viewh12980\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs22 \cf2 \CocoaLigature0 # -*- coding: utf-8 -*-\
\
import os\
import glob\
import time  # \uc0\u49884 \u44036  \u52769 \u51221  \u46972 \u51060 \u48652 \u47084 \u47532  \u52628 \u44032 \
import pandas as pd\
import numpy as np\
from rdkit import Chem\
from rdkit.Chem import AllChem\
from rdkit import DataStructs\
from rdkit import RDLogger\
\
RDLogger.DisableLog('rdApp.*')\
\
# \uc0\u51204 \u52404  \u49884 \u51089  \u49884 \u44036  \u44592 \u47197 \
start_total = time.time()\
\
\pard\pardeftab720\partightenfactor0

\fs24 \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 \outl0\strokewidth0 \strokec3 # 1. Positive (\uc0\u51032 \u50557 \u54408 ) \u45936 \u51060 \u53552  \u51204 \u49688  \u47196 \u46300  \u48143  \u54609 \u44144 \u54532 \u47536 \u53944  \u48320 \u54872 
\fs22 \cf2 \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
print("==> [Step 1] Loading Positive (Drugs) Dataset...")\
start_step1 = time.time()\
positive_df = pd.read_csv('drugs.csv')\
pos_smiles = positive_df['SMILES'].dropna().tolist()\
\
pos_fps = []\
for s in pos_smiles:\
    mol = Chem.MolFromSmiles(s)\
    if mol:\
        pos_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))\
\
\pard\pardeftab720\partightenfactor0

\fs24 \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 \outl0\strokewidth0 \strokec3 # Positive \uc0\u45236 \u48512  \u49345 \u54840  \u50976 \u49324 \u46020  \u44228 \u49328 
\fs22 \cf2 \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
pos_max_sims = []\
for i in range(len(pos_fps)):\
    sims = DataStructs.BulkTanimotoSimilarity(pos_fps[i], pos_fps)\
    sims_filtered = [s for j, s in enumerate(sims) if i != j]\
    if sims_filtered:\
        pos_max_sims.append(max(sims_filtered))\
end_step1 = time.time()\
print("    Finished Step 1. (Time taken: \{:.2f\} seconds)".format(end_step1 - start_step1))\
\
\pard\pardeftab720\partightenfactor0

\fs24 \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 \outl0\strokewidth0 \strokec3 # 2. ZINC DB \uc0\u49692 \u52264 \u51201  \u52376 \u47532  (Batch Processing\u51004 \u47196  \u47700 \u47784 \u47532  \u51208 \u50557 )
\fs22 \cf2 \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
print("==> [Step 2] Screening ZINC Database (64 Cores MPI Running)...")\
start_step2 = time.time()\
zinc_files = glob.glob('zinc_db/*.smi')\
zinc_max_sims = []\
zinc_valid_smiles = []\
\
for file_path in zinc_files:\
    df = pd.read_csv(file_path, sep=r'\\s+', header=None, names=['SMILES', 'ID'])\
    smiles_list = df['SMILES'].dropna().tolist()\
    \
    for s in smiles_list:\
        mol = Chem.MolFromSmiles(s)\
        if mol:\
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)\
            sims = DataStructs.BulkTanimotoSimilarity(fp, pos_fps)\
            max_sim = max(sims)\
            \
            zinc_max_sims.append(max_sim)\
            zinc_valid_smiles.append(s)\
\
zinc_arr = np.array(zinc_max_sims)\
pos_arr = np.array(pos_max_sims)\
end_step2 = time.time()\
print("    Finished Step 2. (Time taken: \{:.2f\} seconds)".format(end_step2 - start_step2))\
\
\pard\pardeftab720\partightenfactor0

\fs24 \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 \outl0\strokewidth0 \strokec3 # 3. \uc0\u45936 \u51060 \u53552  \u44592 \u48152  \u46041 \u51201  \u51076 \u44228 \u44050  \u44228 \u49328 
\fs22 \cf2 \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
zinc_median = np.median(zinc_arr)\
pos_median = np.median(pos_arr)\
optimal_threshold = round(float((zinc_median + pos_median) / 2.0), 2)\
\
\pard\pardeftab720\partightenfactor0

\fs24 \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 \outl0\strokewidth0 \strokec3 # 4. \uc0\u44228 \u49328 \u46108  \u51076 \u44228 \u44050 \u51004 \u47196  \u52572 \u51333  Negative \u45936 \u51060 \u53552  \u54596 \u53552 \u47553  \u48143  \u51200 \u51109 
\fs22 \cf2 \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
final_negative_smiles = []\
for max_sim, smi in zip(zinc_max_sims, zinc_valid_smiles):\
    if max_sim < optimal_threshold:\
        final_negative_smiles.append(smi)\
\
negative_df = pd.DataFrame(\{'SMILES': final_negative_smiles\})\
negative_df.to_csv('final_negative_set.csv', index=False)\
success_ratio = (len(final_negative_smiles) / float(len(zinc_max_sims))) * 100.0\
\
# \uc0\u51204 \u52404  \u51333 \u47308  \u49884 \u44036  \u44228 \u49328 \
end_total = time.time()\
total_duration_sec = end_total - start_total\
total_duration_min = total_duration_sec / 60.0\
\
# 5. 
\fs24 \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 \outl0\strokewidth0 \strokec3 \uc0\u44208 \u44284  \u52636 \u47141 
\fs22 \cf2 \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 \
print("\\n" + "="*60)\
print(" Task 1. \uc0\u45936 \u51060 \u53552  \u44592 \u48152  \u46041 \u51201  \u51076 \u44228 \u44050  \u49440 \u51221  \u48143  \u51204 \u49688  \u51312 \u49324  \u44208 \u44284 ")\
print("="*60)\
print("- \uc0\u52509  \u48516 \u49437 \u46108  \u51032 \u50557 \u54408 (Positive) \u48516 \u51088  \u49688  : \{\} \u44060 ".format(len(pos_fps)))\
print("- \uc0\u52509  \u48516 \u49437 \u46108  ZINC(Negative \u54980 \u48372 ) \u48516 \u51088  \u49688 : \{\} \u44060 ".format(len(zinc_max_sims)))\
print("-" * 60)\
print("[\uc0\u48516 \u54252  \u48516 \u47532  \u48143  \u44228 \u49328 \u46108  \u51076 \u44228 \u44050  \u51648 \u54364 ]")\
print("- ZINC \uc0\u52572 \u45824  \u50976 \u49324 \u46020  \u48516 \u54252  \u51473 \u49900 \u44050  (\u48393 \u50864 \u47532 ) : \{:.4f\}".format(zinc_median))\
print("- Drugs \uc0\u52572 \u45824  \u50976 \u49324 \u46020  \u48516 \u54252  \u51473 \u49900 \u44050  (\u48393 \u50864 \u47532 ): \{:.4f\}".format(pos_median))\
print("\uc0\u55357 \u56393  \u45936 \u51060 \u53552  \u44592 \u48152  \u52572 \u51201  \u46020 \u52636  \u51076 \u44228 \u44050  (Threshold) : \{\}".format(optimal_threshold))\
print("-" * 60)\
print("[\uc0\u52572 \u51333  \u49440 \u48324  \u44208 \u44284 ]")\
print("- \uc0\u51312 \u44148  \u47564 \u51313  \u50976 \u54952  \u51020 \u49457 (Negative) \u48516 \u51088  \u49688  : \{\} \u44060 ".format(len(final_negative_smiles)))\
print("- \uc0\u44228 \u49328 \u46108  \u51076 \u44228 \u44050  \u48120 \u47564  \u54596 \u53552  \u53685 \u44284 \u50984    : \{:.2f\} %".format(success_ratio))\
print("- \uc0\u51032 \u50557 \u54408 (Positive) \u45936 \u51060 \u53552 \u51032  \u54616 \u50948  5% \u44221 \u44228 \u44050  \u49688 \u52824 : \{:.4f\}".format(np.percentile(pos_arr, 5)))\
print("-" * 60)\
print("[\uc0\u9201 \u65039  MPI \u48337 \u47148  \u50672 \u49328  \u49548 \u50836  \u49884 \u44036  \u49457 \u45733  \u47532 \u54252 \u53944 ]")\
print("- \uc0\u51032 \u50557 \u54408 (Drugs) \u54609 \u44144 \u54532 \u47536 \u53944  \u50672 \u49328  \u49884 \u44036    : \{:.2f\} \u52488 ".format(end_step1 - start_step1))\
print("- ZINC \uc0\u45824 \u44508 \u47784  DB \u51204 \u49688  \u51312 \u49324  \u49828 \u53356 \u47532 \u45789  \u49884 \u44036 : \{:.2f\} \u52488 ".format(end_step2 - start_step2))\
print("\uc0\u55357 \u56393  64\u44060  \u53076 \u50612  \u52509  \u49692 \u49688  \u50672 \u49328  \u49548 \u50836  \u49884 \u44036      : \{:.2f\} \u52488  (\{:.2f\} \u48516 )".format(total_duration_sec, total_duration_min))\
print("="*60 + "\\n")}