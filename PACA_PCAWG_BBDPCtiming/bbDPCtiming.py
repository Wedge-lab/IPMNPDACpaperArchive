import sys
import os
import logging
import pandas as pd
from pathlib import Path
from natsort import natsort_keygen

logging.basicConfig(level=logging.INFO,stream=sys.stdout, format="%(levelname)s: %(message)s") #INFO logs go to stdout rather than stderr by default
logger = logging.getLogger(__name__)

# Constants
DPC_SUFFIX = "_DPoutput_2000iters_1000burnin_seed123"
ITERS_TAG  = "2000iters_1000burnin"
CCF_SUBCLONE_THRESHOLD = 0.95

DPC_CCF_COLS = ["cluster", "location_CCF"]

SSDPI_COLS = [
    "chr", "end", "subclonal.CN",
    "nMaj1", "nMin1", "frac1",
    "nMaj2", "nMin2", "frac2",
    "no.chrs.bearing.mut",
]

def _build_paths(sample_id: str, ssDPI_path: str, ssDPCO_path: str, out_path: str) -> dict:
    """Centralise all path construction """
    dpc_dir = Path(ssDPCO_path) / f"{sample_id}{DPC_SUFFIX}"   
    return {
        "ssdpi"     : Path(ssDPI_path) / f"{sample_id}_ssDPI.txt",
        "cluster"   : dpc_dir / f"{sample_id}_{ITERS_TAG}_bestClusterInfo.txt",
        "bed"       : dpc_dir / f"{sample_id}_{ITERS_TAG}_bestConsensusAssignments.bed",
        "merge_out" : Path(out_path) / f"{sample_id}_ssdpi_dpcOutMerge.csv",
        "timing_out": Path(out_path) / f"{sample_id}_bbDPCtimg.csv",
    }


def _classify_timing(row: pd.Series) -> str:
    """
    Assign a timing label to one mutation row.

    Labels
    ------
    unSp         CCF ≤ 0  (missing / artefacts: unspecified)
    subclone     0 < CCF < 0.95
    cloneEarly   clonal + amplified + mutation on >1 copy
    possiblyLate clonal + amplified + minor allele = 1 + mutation on ≤1 copy
    cloneLate    clonal + amplified + other
    cloneNA      clonal, not amplified
    """
    ccf = row.location_CCF

    if ccf <= 0:
        return "unSp"
    if 0 < ccf < CCF_SUBCLONE_THRESHOLD:
        return "subclone"

    amplified = row.nMaj1 > 1 or row.nMaj2 > 1
    if not amplified:
        return "cloneNA"
    if row.noChrsBearingMut > 1:
        return "cloneEarly"
    if (row.nMin1 == 1 or row.nMin2 == 1) and row.noChrsBearingMut <= 1:
        return "possiblyLate"
    return "cloneLate"


def bbDPC_timing(sample_id: str, ssDPI_path: str, ssDPCO_path: str, out_path: str) -> None:
    """
    Merge DPC clustering output with ssDPI data and write timing-labelled CSVs.

    Parameters
    ----------
    sample_id   : Sample identifier.
    ssDPI_path  : Directory containing the ssDPI input file.
    ssDPCO_path : Directory containing DPC output sub-folders.
    out_path    : Directory for output CSVs.
    """
    paths = _build_paths(sample_id, ssDPI_path, ssDPCO_path, out_path)
    logger.info("Processing sample: %s", sample_id)

    # --- 1. Load inputs ---
    ssdpi       = pd.read_csv(paths["ssdpi"],   sep="\t", usecols=SSDPI_COLS)
    best_cluster = pd.read_csv(paths["cluster"], sep="\t")
    dpc_out      = pd.read_csv(paths["bed"],     sep="\t")

    # --- 2. Prepare ssDPI ---
    ssdpi = ssdpi.copy()   # avoid SettingWithCopyWarning on the column we add next
    ssdpi["chr_pos"] = ssdpi["chr"].astype(str) + "_" + ssdpi["end"].astype(str)

    # --- 3. Prepare cluster CCF table ---
    dpc_out = dpc_out[["chr", "end", "cluster"]].fillna(0)

    best_cluster = best_cluster.rename(
        columns={"cluster.no": "cluster", "location": "location_CCF"}
    )[["cluster", "location_CCF"]]

    # --- 4. Attach CCF to every DPC locus ---
    dpc_ccf = (
        dpc_out
        .merge(best_cluster, on="cluster", how='left')
        .assign(chr_pos=lambda df: df["chr"].astype(str) + "_" + df["end"].astype(str))
        .sort_values("chr_pos", key=natsort_keygen())
        [["chr_pos", "cluster", "location_CCF"]]
    )

    # --- 5. Fill in ssDPI loci absent from DPC output, Merge & save intermediate ---
    missing_loci = set(ssdpi["chr_pos"]) - set(dpc_ccf["chr_pos"])
    if missing_loci:
        logger.info("  %d ssDPI loci have no CNA entry; assigning cluster=-1", len(missing_loci))
    merged = (
        ssdpi
            .merge(dpc_ccf, on='chr_pos', how='left')
            .rename(columns={'no.chrs.bearing.mut':'noChrsBearingMut'})
    )
    merged[DPC_CCF_COLS] = merged[DPC_CCF_COLS].fillna(-1)
    merged.to_csv(paths["merge_out"], index=False)
    logger.info("  Saved merge file → %s", paths["merge_out"])

    # --- 6. Classify timing & save final output ---
    merged["bbDPCtiming"] = merged.apply(_classify_timing, axis=1)
    (
        merged
        .sort_values("chr_pos", key=natsort_keygen())
        .to_csv(paths["timing_out"], index=False)
    )
    logger.info("  Saved timing file → %s", paths["timing_out"])


def main() -> None:
    if len(sys.argv) != 5:
        print("Usage: python3 bbdpctiming.py <sampleID> <ssdpcInputPath> <ssdpcOutputPath> <outputPath>")
        sys.exit(1)
    _, sample_id, ssDPI_path, ssDPCO_path, out_path = sys.argv
    bbDPC_timing(sample_id, ssDPI_path, ssDPCO_path, out_path)


if __name__ == "__main__":
    main()
    
