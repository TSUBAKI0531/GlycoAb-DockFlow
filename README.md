# GlycoAb-DockFlow
### 抗体－糖鎖抱合抗原 ドッキング・解析統合パイプライン

本リポジトリは、糖鎖抱合ワクチン（Glycoconjugate Vaccine）の開発において、AlphaFold 3を活用した抗原－抗体複合体の構造予測と、その相互作用解析を自動化するPythonツールキットです。

## 背景
糖鎖抱合抗原は、キャリアタンパク質、リンカー、糖鎖抗原（TACA等）が複雑に組み合わさっており、抗体との相互作用部位（パラトープ）の特定が困難です。本ツールは、これら一連のワークフローを計算機科学的に支援します。

## 主な機能
1. **Automated AF3 Input Generation**: 
   - 複雑な糖鎖SMILESから、リンカーの結合原子をSMARTSパターンを用いて自動特定。
   - 共有結合指定を含むAlphaFold 3用ドッキングJSONを自動生成。
2. **Paratope Extraction & CDR Mapping**: 
   - 予測された複合体構造（CIF）から、糖鎖およびタンパク質に接触する抗体残基を自動抽出。
   - ANARCIを用いたIMGTナンバリングに基づき、接触残基をCDR領域へ自動マッピング。
3. **Quantitative Ranking**: 
   - SASA（溶媒露出面積）等の指標に基づき、複数の予測モデルから最適な候補を選別。

## インストール
```bash
conda install -c bioconda anarci
pip install rdkit-pypi biopython pandas freesasa py3dmol

## 使い方
1. AntibodyDockingWorkflow を用いてドッキング用JSONを作成。
2. AlphaFold 3 で計算実行。
3. analyze_paratope() を用いて、結合部位の解析レポート（CSV）を出力。

## 著者
[TSUBAKI0531] (Ph.D. in Agriculture)
