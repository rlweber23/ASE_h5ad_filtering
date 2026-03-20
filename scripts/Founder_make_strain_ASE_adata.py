#!/usr/bin/env python3
"""
Build per-strain allele-specific AnnData files for a given tissue.

Example
-------
python make_ase_adatas.py -t Heart
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


# ----------------------------- CONFIG -----------------------------

TISSUE_TO_PLATE = {
    "CortexHippocampus": ["igvf_003","igvf_011","igvf_012"],
    "Heart": ["igvf_003","igvf_004","igvf_012"],
    "Liver": ["igvf_004","igvf_005"],
    "Gastrocnemius": ["igvf_010","igvf_011","igvf_012"],
    "Adrenal": ["igvf_008","igvf_008b","igvf_009","igvf_012"],
    "DiencephalonPituitary": ["igvf_005","igvf_007","igvf_012"],
    "Kidney": ["igvf_009","igvf_010"],
    "GonadsFemale": ["igvf_008","igvf_008b","igvf_009"],
    "GonadsMale": ["igvf_008","igvf_008b","igvf_009"]
}
ALT_TO_F1 = {
    "A_J": "B6AF1J",
    "NOD_ShiLtJ": "B6NODF1J",
    "NZO_HlLtJ": "B6NZOF1J",
    "WSB_EiJ": "B6WSBF1J",
    "129S1_SvImJ": "B6129S1F1J",
    "CAST_EiJ": "B6CASTF1J",
    "PWK_PhJ": "B6PWKF1J",
}

strain_short_dict = {
    "A_J": "AJ",
    "NOD_ShiLtJ": "NODJ",
    "NZO_HlLtJ": "NZOJ",
    "WSB_EiJ": "WSBJ",
    "129S1_SvImJ": "129S1J",
    "CAST_EiJ": "CASTJ",
    "PWK_PhJ": "PWKJ",
}

ALT_STRAINS: List[str] = [
    "129S1_SvImJ",
    "A_J",
    "CAST_EiJ",
    "NOD_ShiLtJ",
    "NZO_HlLtJ",
    "PWK_PhJ",
    "WSB_EiJ",
]


# --------------------------- HELPERS ------------------------------

def get_bc1(text: str) -> str:
    return text[-8:]


def get_bc2(text: str) -> str:
    return text[8:16]


def get_bc3(text: str) -> str:
    return text[:8]


def load_barcode_maps(barcode_map_csv: Path) -> Dict[str, Dict[str, str]]:
    """Return {'bc1': {sequence->well}, 'bc2': {...}, 'bc3': {...}}."""
    df = pd.read_csv(barcode_map_csv)
    required_cols = {"bc", "bc_sequence", "bc_well"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"barcode_map_df is missing columns: {missing}")
    return (
        df.groupby("bc")
          .apply(lambda g: dict(zip(g["bc_sequence"], g["bc_well"])))
          .to_dict()
    )


def load_biotype_and_chrom_maps(genes_with_biotype_tsv: Path,
                                gene_id_chr_tsv: Path) -> tuple[set[str], Dict[str, str]]:
    """Return (protein_coding_gene_ids, gene_id -> chromosome)"""

    g36 = (
        pd.read_csv(genes_with_biotype_tsv, sep="\t", header=None)
        .rename(columns={0: "gene_id", 1: "gene_name", 2: "biotype"})
    )
    g36["gene_id_stripped"] = g36["gene_id"].str.split(".").str[0]
    protein_coding = set(g36.loc[g36["biotype"] == "protein_coding", "gene_id_stripped"])

    id_loc = pd.read_csv(gene_id_chr_tsv, sep="\t", header=None)
    # assume columns: [chrom, gene_id_with_version, ...]
    id_loc["gene_id_stripped"] = id_loc[1].str.split(".").str[0]
    chrom_map = dict(zip(id_loc["gene_id_stripped"], id_loc[0]))
    return protein_coding, chrom_map


def find_h5ads(directory: Path) -> List[Path]:
    if not directory.exists():
        return []
    return sorted([p for p in directory.iterdir() if p.suffix == ".h5ad"])


def ensure_columns(adata: ad.AnnData, cols: List[str], where: str = "obs"):
    """Raise a clear error if required columns are missing."""
    frame = adata.obs if where == "obs" else adata.var
    missing = [c for c in cols if c not in frame.columns]
    if missing:
        raise KeyError(f"Missing columns in adata.{where}: {missing}")


def process_tissue(
    tissue: str,
    base_path: Path,
    barcode_map_csv: Path,
    genes_with_biotype_tsv: Path,
    gene_id_chr_tsv: Path,
    annotation_dir: Path,
    out_dir: Path,
):
    
    
    if tissue not in TISSUE_TO_PLATE:
        raise KeyError(f"Unknown tissue '{tissue}'. Valid keys: {sorted(TISSUE_TO_PLATE)}")
    
    barcode_maps = load_barcode_maps(barcode_map_csv)
    pc_genes, chrom_map = load_biotype_and_chrom_maps(genes_with_biotype_tsv, gene_id_chr_tsv)
    
    
    # Load annotations once
    anno_csv = annotation_dir / f"{tissue}_Founder_processed_annotated_obs.csv"
    # if not anno_csv.exists():
        # raise FileNotFoundError(f"Annotation file not found: {anno_csv}")
    
    annotations = pd.read_csv(anno_csv)
    
    if "cellID" not in annotations.columns:
        raise KeyError("Annotation file must include a 'cellID' column.")
    subtype_dict = dict(zip(annotations["cellID"], annotations.get("subtype")))
    # Keep just the columns we need (presence-checked)
    keep_cols = [
        "cellID", "Mouse_Tissue_ID", "Sex", "Genotype",
        "leiden", "general_celltype", "general_CL_ID", "celltype", "CL_ID"
    ]
    keep_cols = [c for c in keep_cols if c in annotations.columns]
    sub_anno = annotations[keep_cols].set_index("cellID")
    
    
    
    plate = TISSUE_TO_PLATE[tissue]
    
    plate_paths = []
    for plate_sub in plate:
        tissue_root_sub = base_path / plate_sub
    
        plate_paths.append(tissue_root_sub)
    
    
    for strain in ALT_STRAINS:
        adatas = []
        
        for tissue_root in plate_paths:
            print(f"\n[INFO] === {strain} ===")
            strain_dir = tissue_root / strain
            h5ad_files = find_h5ads(strain_dir)
            if not h5ad_files:
                print(f"[WARN] No .h5ad files found in {strain_dir}. Skipping.")
                continue
            
            for f in h5ad_files:
                print(f"[INFO] Reading: {tissue_root.stem} {f.name}")
                basename = os.path.basename(f)
                adata = sc.read_h5ad(str(f))
                adata.obs.reset_index(inplace = True)
                #adata.obs.reset_index(drop = True, inplace = True)
                #adata.obs.rename(columns = {'bc':'barcode'}, inplace = True)
                # We need 'barcode' in obs
                ensure_columns(adata, ["barcode"], where="obs")
    
                # parse barcodes
                adata.obs["bc1_sequence"] = adata.obs["barcode"].astype(str).apply(get_bc1)
                adata.obs["bc2_sequence"] = adata.obs["barcode"].astype(str).apply(get_bc2)
                adata.obs["bc3_sequence"] = adata.obs["barcode"].astype(str).apply(get_bc3)
    
                # map to wells
                adata.obs["bc1_well"] = adata.obs["bc1_sequence"].map(barcode_maps.get("bc1", {}))
                adata.obs["bc2_well"] = adata.obs["bc2_sequence"].map(barcode_maps.get("bc2", {}))
                adata.obs["bc3_well"] = adata.obs["bc3_sequence"].map(barcode_maps.get("bc3", {}))
    
                # ensure strings (avoid 'nan' object issues)
                for col in ["bc1_well", "bc2_well", "bc3_well"]:
                    adata.obs[col] = adata.obs[col].astype(str)
    
                subpool = basename.split('.')[0]
                adata.obs["subpool"] = subpool
    
                # plate_id = TISSUE_TO_PLATE[tissue]
                plate_id = tissue_root.stem
                adata.obs["cellID"] = (
                    adata.obs["bc1_well"] + "_" +
                    adata.obs["bc2_well"] + "_" +
                    adata.obs["bc3_well"] + "_" +
                    adata.obs["subpool"] + "_" +
                    plate_id
                )
    
                adatas.append(adata)  
    
                
            if not adatas:
                print(f"[WARN] No valid AnnData loaded for {strain}. Skipping.")
                continue
            
            
            merged = ad.concat(adatas, join="outer", axis=0, merge="same")
            print(f"[INFO] Merged shape for {strain}: {merged.shape}")
    
            
            if "gene_id" not in merged.var.columns:
                print("[WARN] 'gene_id' not found in .var; creating from var_names.")
                merged.var["gene_id"] = merged.var_names
            
            merged.var["gene_id"] = merged.var["gene_id"].astype(str).str.split(".").str[0]
            merged.var["gene_name_unmodified"] = merged.var_names.str.replace(f"{strain}_", "", regex=False)
            merged.var["gene_id_unmodified"] = merged.var["gene_id"].str.replace(f"{strain}_", "", regex=False)
            
            
            merged.var["allele"] = np.where(
                merged.var_names.str.startswith(strain),
                strain,
                "B6J"
            )
            merged.var["B6J_allele"] = merged.var["allele"] == "B6J"
            merged.var["ALT_allele"] = merged.var["allele"] == strain
            
            # merged = merged[:, merged.var["gene_id_unmodified"].isin(pc_genes)].copy()
            merged.var["chromosome"] = merged.var["gene_id_unmodified"].map(chrom_map)
    
            chromosomes_to_exclude = ['chrY','chrM','chr1_GL456210v1_random','chr1_GL456212v1_random',
                                      'chr1_GL456221v1_random','chr1_GL456211v1_random','chr5_GL456354v1_random','chr5_JH584296v1_random',
                                      'chr5_JH584297v1_random','chr5_JH584298v1_random','chr5_JH584299v1_random',
                                      'chr7_GL456219v1_random','chrUn_JH584304v1','chrY_JH584303v1_random']
            
            merged = merged[:, ~merged.var["chromosome"].isin(chromosomes_to_exclude)].copy()
            
            merged.obs.index = merged.obs["cellID"].astype(str)
            
            
            merged.obs["subtype"] = merged.obs.index.map(subtype_dict)
            print(f'merged size: {merged.shape}')
            print(f'{merged.obs.head()}')
            print('===============================================================================')
            
            
            merged_filt = merged[~merged.obs["subtype"].isna()].copy()
            print(f'merged_filt size: {merged_filt.shape}')
            
            
            merged_filt.obs = merged_filt.obs.join(sub_anno, how="left")
            # Remove low quality
            if "subtype" in merged_filt.obs:
                merged_filt = merged_filt[merged_filt.obs["subtype"] != "low quality"].copy()
            
            strain_short = strain_short_dict[strain]
            merged_filt_strain = merged_filt[merged_filt.obs['Genotype'].isin(['B6J',strain_short])]
            
            out_tissue_dir = out_dir / tissue
            out_tissue_dir.mkdir(parents=True, exist_ok=True)
            out_fp = out_tissue_dir / f"{strain}_{tissue}_Founder_ASE.h5ad"
            print(f"[INFO] Writing: {out_fp}  Shape: {merged_filt_strain.shape}")
            merged_filt_strain.write_h5ad(out_fp)
    
    
           
    



def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create tissue- and strain-specific ASE AnnData files."
    )
    p.add_argument("-t", "--tissue", required=True, help=f"Tissue key: {sorted(TISSUE_TO_PLATE)}")

    # Optional path overrides (defaults match your environment)
    p.add_argument(
        "--base-path",
        default="/share/crsp/lab/seyedam/share/igvf_allele_specific_mapping/ASE_kb/test/cellbender/",
        help="Base directory containing plate folders (igvf_0xx).",
    )
    p.add_argument(
        "--barcode-map",
        default="/share/crsp/lab/seyedam/weberrl/F1_IGVF/ref_files/barcode_map_df.csv",
        help="CSV with columns: bc, bc_sequence, bc_well.",
    )
    p.add_argument(
        "--genes-with-biotype",
        default="/share/crsp/lab/seyedam/weberrl/genome_files/gencode_36/igvf_gtf_gene_biotypes.tsv",
        help="TSV with gene_id, gene_name, biotype (no header).",
    )
    p.add_argument(
        "--gene-id-chr-map",
        default="/share/crsp/lab/seyedam/share/igvf_references/gencode_M36_gene_ID_chromosome_map.tsv",
        help="TSV with chromosome and gene_id(with version).",
    )
    p.add_argument(
        "--annotation-dir",
        default="/share/crsp/lab/seyedam/weberrl/F1_IGVF/standard_processing/anndatas/annotated_processed/",
        help="Directory with {tissue}_F1_processed_annotated_obs.csv files.",
    )
    p.add_argument(
        "--out-dir",
        default="/share/crsp/lab/seyedam/weberrl/F1_IGVF/allele_specific_processing/g36/tissue_ASE_adatas_processed",
        help="Output directory (per-strain subfolders will be created).",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    process_tissue(
        tissue=args.tissue,
        base_path=Path(args.base_path),
        barcode_map_csv=Path(args.barcode_map),
        genes_with_biotype_tsv=Path(args.genes_with_biotype),
        gene_id_chr_tsv=Path(args.gene_id_chr_map),
        annotation_dir=Path(args.annotation_dir),
        out_dir=Path(args.out_dir),
    )





































