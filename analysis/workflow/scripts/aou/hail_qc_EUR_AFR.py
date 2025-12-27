#!/usr/env/bin python
"""
My own hail QC script for filtering phased GTs
Adapted from https://github.com/CAST-genomics/cast-workflows/blob/0e9f9a6e/gwas/aou/hail_runner_EUR_AFR.py
"""

import os
import argparse
from pathlib import Path

import hail as hl
import pandas as pd

SMALLNUM = 10e-400


class HailRunner:

    import hail as hl

    def __init__(
        self,
        gts: Path,
        output: Path,
        sample_files_dir: Path = None,
        sample_call_rate: float = 0.9,
        variant_call_rate: float = 0.9,
        MAF: float = 0.01,
        HWE: float = 1e-100,
        GQ: float = 20,
    ):
        self.gts = gts
        self.output = output
        self.samples_file_dir = sample_files_dir
        self.sample_call_rate = sample_call_rate
        self.variant_call_rate = variant_call_rate
        self.MAF = MAF
        self.HWE = HWE
        self.GQ = GQ
        self.gwas = None
        self.data = None
        self.method = "hail"
        self.setup()

    def run(self):
        # Set up hail
        hl.init(default_reference="GRCh38")

        # Load genotypes
        mt = hl.import_vcf(self.gts)

        if self.samples_file_dir is None:
            os.system(
                "gsutil -u ${GOOGLE_PROJECT} cp ${WORKSPACE_BUCKET}/samples/EUR_WHITE.csv ."
            )
            os.system(
                "gsutil -u ${GOOGLE_PROJECT} cp ${WORKSPACE_BUCKET}/samples/AFR_BLACK.csv ."
            )
            samps_dir = Path(".")
        else:
            samps_dir = self.samples_file_dir
        eur_sample_ids = pd.read_csv(samps_dir / "EUR_WHITE.csv")["person_id"].astype(str)
        afr_sample_ids = pd.read_csv(samps_dir / "AFR_BLACK.csv")["person_id"].astype(str)
        eur_tbl = hl.Table.from_pandas(pd.DataFrame(eur_sample_ids), key="person_id")
        afr_tbl = hl.Table.from_pandas(pd.DataFrame(afr_sample_ids), key="person_id")
        ids = pd.DataFrame(self.ptcovar["person_id"])
        sample_tbl = hl.Table.from_pandas(ids, key="person_id")

        data = mt.filter_cols(
            hl.is_defined(eur_tbl[mt.s])
            | hl.is_defined(afr_tbl[mt.s])
            | hl.is_defined(sample_tbl[mt.s])
        )

        # filter multiallelics
        data = data.filter_rows(hl.len(data.alleles) == 2)

        # Genotype QC
        data = data.annotate_entries(FT=hl.coalesce(data.FT, "PASS"))
        data = data.filter_entries(data.FT == "PASS")
        data = data.filter_entries(data.GQ >= self.GQ)  # 20

        # Run variant_qc separately for each group
        data = data.annotate_cols(
            eur_cohort=hl.is_defined(eur_tbl[data.s]),
            afr_cohort=hl.is_defined(afr_tbl[data.s]),
            sample_cohort=hl.is_defined(sample_tbl[data.s]),
        )

        eur_qc = hl.variant_qc(data.filter_cols(data.eur_cohort))
        afr_qc = hl.variant_qc(data.filter_cols(data.afr_cohort))
        sample_qc = hl.variant_qc(data.filter_cols(data.sample_cohort))

        # Annotate per-group QC back to the main MT
        eur_rows = eur_qc.rows()
        afr_rows = afr_qc.rows()
        sample_rows = sample_qc.rows()

        data = data.annotate_rows(
            eur_AF=eur_rows[data.row_key].variant_qc.AF,
            eur_HWE=eur_rows[data.row_key].variant_qc.p_value_hwe,
            afr_AF=afr_rows[data.row_key].variant_qc.AF,
            afr_HWE=afr_rows[data.row_key].variant_qc.p_value_hwe,
            sample_AF=sample_rows[data.row_key].variant_qc.AF,
            sample_HWE=sample_rows[data.row_key].variant_qc.p_value_hwe,
        )

        # Variant filter: keep if either group passes all thresholds
        data = data.filter_rows(
            ((hl.min(data.eur_AF) >= self.MAF) & (data.eur_HWE > self.HWE))
            | ((hl.min(data.afr_AF) >= self.MAF) & (data.afr_HWE > self.HWE))
            | ((hl.min(data.sample_AF) >= self.MAF) & (data.sample_HWE > self.HWE))
        )

        # now filter samples to given cohort
        # ids = pd.DataFrame(self.ptcovar['person_id'])
        # sample_tbl = hl.Table.from_pandas(ids, key="person_id")
        data = data.filter_cols(hl.is_defined(sample_tbl[data.s]))

        # sample QC
        data = hl.sample_qc(data)
        data = data.filter_cols(
            data.sample_qc.call_rate >= self.sample_call_rate, keep=True
        )  # 0.9
        # remaining variant QC
        data = hl.variant_qc(data)
        data = data.filter_rows(
            data.variant_qc.call_rate >= self.variant_call_rate, keep=True
        )  # 0.9

        # Keep track of data
        self.data = data

        # write data out
        hl.export_vcf(self.data, self.output)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("output", help="path to output file")
    parser.add_argument("genotypes", help="path to genotypes .vcf.bgz file")
    parser.add_argument(
        "--samples-file-dir",
        help="Path to EUR_WHITE.csv and AFR_BLACK.csv files",
        default=None,
    )
    parser.add_argument(
        "--sample-call-rate",
        help="Apply minimum sample call rate QC",
        type=float,
        default=0.90,
    )
    parser.add_argument(
        "--variant-call-rate",
        help="Apply minimum variant call rate QC",
        type=float,
        default=0.90,
    )
    parser.add_argument(
        "--MAF", help="Apply minor allele frequency QC", type=float, default=0.01
    )
    parser.add_argument(
        "--HWE", help="Apply HWE p-value cutoff QC", type=float, default=1e-100
    )
    parser.add_argument(
        "--GQ", help="Apply minimun genotype score QC", type=int, default=20
    )
    args = parser.parse_args()

    runner = HailRunner(
        Path(args.output),
        Path(args.genotypes),
        Path(args.samples_file_dir) if args.samples_file_dir is not None else None,
        sample_call_rate=args.sample_call_rate,
        variant_call_rate=args.variant_call_rate,
        MAF=args.MAF,
        HWE=args.HWE,
        GQ=args.GQ,
    )
    runner.run()


if __name__ == "__main__":
    main()
