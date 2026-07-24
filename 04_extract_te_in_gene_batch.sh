#!/usr/bin/env bash
# Batch extract TE intervals overlapping gene windows (gene body +/- 2 kb).
# Requires: gene BEDs under ./up2k_gene_down2k_beds/ and TE BEDs from step 02.
# Adjust paths if your local layout differs.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SCRIPT_DIR}/03_extract_te_in_gene.py"
python "$PY" --gene-bed ./up2k_gene_down2k_beds/ECD04_up2k_down2k.bed --te-bed Brapa_ECD04_EDTA.bed --output Brapa_ECD04_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/C4012_v1_up2k_down2k.bed --te-bed Bcarinata_C4012_v1_EDTA.bed --output Bcarinata_C4012_v1_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/zd_1_up2k_down2k.bed --te-bed Bcarinata_zd_1_EDTA.bed --output Bcarinata_zd_1_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/ASM1677196v1_up2k_down2k.bed --te-bed Bcarinata_ASM1677196v1_EDTA.bed --output Bcarinata_ASM1677196v1_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/ASM4058406v1_up2k_down2k.bed --te-bed Bcarinata_ASM4058406v1_EDTA.bed --output Bcarinata_ASM4058406v1_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/10167_v1_up2k_down2k.bed --te-bed Bcarinata_10167_v1_EDTA.bed --output Bcarinata_10167_v1_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/HJ_up2k_down2k.bed --te-bed Bjuncea_HJ_EDTA.bed --output Bjuncea_HJ_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/varuna.v0_up2k_down2k.bed --te-bed Bjuncea_varuna.v0_EDTA.bed --output Bjuncea_varuna.v0_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/tumida.v1.5_up2k_down2k.bed --te-bed Bjuncea_tumida.v1.5_EDTA.bed --output Bjuncea_tumida.v1.5_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/W1701_up2k_down2k.bed --te-bed Boleracea_W1701_EDTA.bed --output Boleracea_W1701_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/PL021_up2k_down2k.bed --te-bed Boleracea_PL021_EDTA.bed --output Boleracea_PL021_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/D101_up2k_down2k.bed --te-bed Boleracea_D101_EDTA.bed --output Boleracea_D101_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/HDEM_up2k_down2k.bed --te-bed Boleracea_HDEM_EDTA.bed --output Boleracea_HDEM_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/No.06-9-28_up2k_down2k.bed --te-bed Boleracea_No.06-9-28_EDTA.bed --output Boleracea_No.06-9-28_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/M249_up2k_down2k.bed --te-bed Boleracea_M249_EDTA.bed --output Boleracea_M249_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/07-DH-33_up2k_down2k.bed --te-bed Boleracea_07-DH-33_EDTA.bed --output Boleracea_07-DH-33_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/JS_up2k_down2k.bed --te-bed Boleracea_JS_EDTA.bed --output Boleracea_JS_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/C2_up2k_down2k.bed --te-bed Bnigra_C2_EDTA.bed --output Bnigra_C2_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Ni100_up2k_down2k.bed --te-bed Bnigra_Ni100_EDTA.bed --output Bnigra_Ni100_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Sangam_up2k_down2k.bed --te-bed Bnigra_Sangam_EDTA.bed --output Bnigra_Sangam_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/CN115125_up2k_down2k.bed --te-bed Bnigra_CN115125_EDTA.bed --output Bnigra_CN115125_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/ZS11_up2k_down2k.bed --te-bed Bnapus_ZS11_EDTA.bed --output Bnapus_ZS11_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Shengli_up2k_down2k.bed --te-bed Bnapus_Shengli_EDTA.bed --output Bnapus_Shengli_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Xiaoyun_up2k_down2k.bed --te-bed Bnapus_Xiaoyun_EDTA.bed --output Bnapus_Xiaoyun_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Darmor_v10_up2k_down2k.bed --te-bed Bnapus_Darmor_v10_EDTA.bed --output Bnapus_Darmor_v10_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Tapidor_up2k_down2k.bed --te-bed Bnapus_Tapidor_EDTA.bed --output Bnapus_Tapidor_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/QuintaA_up2k_down2k.bed --te-bed Bnapus_QuintaA_EDTA.bed --output Bnapus_QuintaA_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Gangan_up2k_down2k.bed --te-bed Bnapus_Gangan_EDTA.bed --output Bnapus_Gangan_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Express617_up2k_down2k.bed --te-bed Bnapus_Express617_EDTA.bed --output Bnapus_Express617_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Westar_up2k_down2k.bed --te-bed Bnapus_Westar_EDTA.bed --output Bnapus_Westar_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Ningyou7_up2k_down2k.bed --te-bed Bnapus_Ningyou7_EDTA.bed --output Bnapus_Ningyou7_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/No2127_up2k_down2k.bed --te-bed Bnapus_No2127_EDTA.bed --output Bnapus_No2127_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Zheyou7_up2k_down2k.bed --te-bed Bnapus_Zheyou7_EDTA.bed --output Bnapus_Zheyou7_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Z1_up2k_down2k.bed --te-bed Brapa_Z1_EDTA.bed --output Brapa_Z1_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/Chiifu_up2k_down2k.bed --te-bed Brapa_Chiifu_EDTA.bed --output Brapa_Chiifu_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/HSCT_up2k_down2k.bed --te-bed Brapa_HSCT_EDTA.bed --output Brapa_HSCT_EDTA_gene.bed
python "$PY" --gene-bed ./up2k_gene_down2k_beds/AJ_up2k_down2k.bed --te-bed Brapa_AJ_EDTA.bed --output Brapa_AJ_EDTA_gene.bed
