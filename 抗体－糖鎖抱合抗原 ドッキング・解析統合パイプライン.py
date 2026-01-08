import json
import os
import pandas as pd
from rdkit import Chem
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.SeqUtils import seq1

# ANARCIがインストールされている環境を想定
try:
    from anarci import anarci
except ImportError:
    anarci = None

class AntibodyDockingWorkflow:
    def __init__(self, job_name, h_chain="H", l_chain="L", ant_prot_chain="A", ant_gly_chain="B"):
        self.job_name = job_name
        self.h_chain = h_chain
        self.l_chain = l_chain
        self.ant_prot_chain = ant_prot_chain
        self.ant_gly_chain = ant_gly_chain

    # --- PHASE 1: ドッキング準備 (JSON生成) ---

    def _get_af3_atom_name(self, smiles, smarts_pattern):
        """RDKitを使用してSMILESからAF3形式の原子名を特定"""
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts(smarts_pattern)
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            raise ValueError(f"Pattern '{smarts_pattern}' not found in SMILES.")
        
        target_idx = matches[0][0]
        atom = mol.GetAtomWithIdx(target_idx)
        symbol = atom.GetSymbol()
        count = sum(1 for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() == symbol and i <= target_idx)
        return f"{symbol}{count}"

    def generate_docking_json(self, prot_seq, gly_smiles, h_seq, l_seq, bond_res_idx, terminal_smarts):
        """抗体と糖鎖抱合抗原の複合体予測用JSONを作成"""
        bond_atom_ligand = self._get_af3_atom_name(gly_smiles, terminal_smarts)
        
        data = {
            "name": f"{self.job_name}_Docking",
            "modelSeeds": [1],
            "sequences": [
                {"protein": {"id": self.ant_prot_chain, "sequence": prot_seq}},
                {"protein": {"id": self.h_chain, "sequence": h_seq}},
                {"protein": {"id": self.l_chain, "sequence": l_seq}},
                {"ligand": {"id": self.ant_gly_chain, "smiles": gly_smiles}}
            ],
            "bondedAtomPairs": [
                {
                    "at1": {"resChainId": self.ant_prot_chain, "resIdx": bond_res_idx, "atomName": "NZ"},
                    "at2": {"resChainId": self.ant_gly_chain, "resIdx": 1, "atomName": bond_atom_ligand}
                }
            ]
        }
        
        path = f"{self.job_name}_docking_input.json"
        with open(path, "w") as f:
            json.dump(data, f, indent=4)
        print(f"✅ [PHASE 1] Docking JSON created: {path}")

    # --- PHASE 2: 解析 (パラトープ抽出 & CDRマッピング) ---

    def _get_cdr_labels(self, sequence, scheme="imgt"):
        """ANARCIを用いて残基番号とCDRラベルの対応表を作成"""
        if anarci is None: return {}
        results = anarci([("seq", sequence)], scheme=scheme)
        numbering = results[0][0]
        if numbering is None: return {}

        labels = {}
        for pos_data, aa in numbering[0]:
            num = pos_data[0]
            label = "FW"
            if 27 <= num <= 38: label = "CDR1"
            elif 56 <= num <= 65: label = "CDR2"
            elif 105 <= num <= 117: label = "CDR3"
            labels[num] = label
        return labels

    def analyze_paratope(self, cif_file, distance_cutoff=4.5):
        """複合体構造からパラトープを抽出しレポート化"""
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("complex", cif_file)[0]
        
        # 原子群の定義
        antigen_atoms = [a for c in [self.ant_prot_chain, self.ant_gly_chain] for a in structure[c].get_atoms()]
        antibody_chains = {self.h_chain: "Heavy", self.l_chain: "Light"}
        
        ns = NeighborSearch([a for c in antibody_chains.keys() for a in structure[c].get_atoms()])
        contacts = set()
        for a_atom in antigen_atoms:
            for res in ns.search(a_atom.coord, distance_cutoff, level="R"):
                contacts.add(res)

        # データ集計
        report = []
        for chain_id, chain_type in antibody_chains.items():
            seq = "".join([seq1(r.get_resname()) for r in structure[chain_id] if "CA" in r])
            cdr_map = self._get_cdr_labels(seq)
            
            for res in sorted(list(contacts), key=lambda x: x.id[1]):
                if res.get_parent().id == chain_id:
                    res_num = res.id[1]
                    report.append({
                        "Chain": chain_type,
                        "ResNum": res_num,
                        "ResName": res.get_resname(),
                        "CDR_Region": cdr_map.get(res_num, "FW/Unknown"),
                        "Target": "Antigen"
                    })
        
        df = pd.DataFrame(report)
        print(f"✅ [PHASE 2] Paratope analysis complete for: {cif_file}")
        return df

# --- 実務運用シナリオ ---

# 1. パイプライン初期化
pipeline = AntibodyDockingWorkflow("AntiTACA_V1")

# 2. 【準備】JSON生成 (AlphaFold3に投入する前に実行)
pipeline.generate_docking_json(
    prot_seq="MKTII...", # 抗原タンパク
    gly_smiles="CC(=O)N...", # 糖鎖
    h_seq="EVQLV...", # 抗体H鎖
    l_seq="DIQMT...", # 抗体L鎖
    bond_res_idx=105,
    terminal_smarts="C(=O)N"
)

# --- (ここでAlphaFold3で構造予測を実行) ---

# 3. 【解析】結果ファイル(CIF)からパラトープを自動抽出
# result_df = pipeline.analyze_paratope("AntiTACA_V1_Docking_model_0.cif")
# print(result_df)
# result_df.to_csv("paratope_report.csv", index=False)