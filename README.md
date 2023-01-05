# Format Coveter

This is a take home interview assignment. Please write a Python program that converts the provided input files to the specified output files.
Please apply best practices to organize and document your Python program. Please feel free to use third-party libraries (e.g. pandas).

## Usage

Please write a Python program that can be invoked as a bash command. The program should take one positional argument (`input_manifest`) and one optional argument (`output_directory`).
If the optional argument `output_directory` is not provided, the default should be the current directory.

```bash
format_conveter [--output_directory OUTPUT_DIRECTORY] input_manifest
```

## Input

### Manifest

| Sample ID  | Subject ID | MRD Type | Specimen Type | Visit Code | sampleDir                    | mrdDir           |
| ---------- | ---------- | -------- | ------------- | ---------- | ---------------------------- | ---------------- |
| 86755_94   | NHL_2_4    | Baseline | Plasma        | C1D1       | /path/to/main_dir/86755_94   | /path/to/mrd_dir |
| NHL_2.1    | NHL_2_4    | FollowUp | Plasma        | C2D1       | /path/to/main_dir/NHL_2.1    | /path/to/mrd_dir |
| 126713_131 | NHL_2_4    | Normal   | PBMC          | C1D1       | /path/to/main_dir/126713_131 | /path/to/mrd_dir |
| 158376_62  | NHL_6_8    | Baseline | Plasma        | C1D1       | /path/to/main_dir/158376_62  | /path/to/mrd_dir |
| NHL_6.1    | NHL_6_8    | FollowUp | Plasma        | C2D1       | /path/to/main_dir/NHL_6.1    | /path/to/mrd_dir |
| 165520_66  | NHL_6_8    | Normal   | PBMC          | C1D1       | /path/to/main_dir/165520_66  | /path/to/mrd_dir |

### sampleDir

The output from the main pipeline. A set of input files are provided under `data/input/main_dir`.
Each sample has one sub-folder. The naming convention is `[Sample ID]/qcMetrics/[Sample ID]_qc_metrics_summary.csv`.
The file contains the QC metrics of each sample.

```bash
main_dir/
├── 126713_131
│   └── qcMetrics
│       └── 126713_131_qc_metrics_summary.csv
├── 158376_62
│   └── qcMetrics
│       └── 158376_62_qc_metrics_summary.csv
├── 165520_66
│   └── qcMetrics
│       └── 165520_66_qc_metrics_summary.csv
├── 86755_94
│   └── qcMetrics
│       └── 86755_94_qc_metrics_summary.csv
├── NHL_2.1
│   └── qcMetrics
│       └── NHL_2.1_qc_metrics_summary.csv
└── NHL_6.1
    └── qcMetrics
        └── NHL_6.1_qc_metrics_summary.csv
```

### mrdDir

The output from the MRD pipeline. A set of input files are provided under `data/input/mrd_dir`.
Normal (PBMC) samples are excluded at this stage.

- Under `reporters`, each plasma sample has one sub-folder. The naming convention is `reporters/[Sample ID]/[Sample ID]_*.vcf`.
  The files contain SNVs detected in each sample.
- Under `monitor`, each subject has one sub-folder. The naming convention is `monitor/[Subject ID]/[Sample ID]/snv_out_withsubject.txt`.
  The file contains the ctDNA level of each plasma sample (e.g. mean AF).

```bash
mrd_dir/
├── monitor
│   ├── NHL_2_4
│   │   ├── 86755_94
│   │   │   └── snv_out_withsubject.txt
│   │   └── NHL_2.1
│   │       └── snv_out_withsubject.txt
│   └── NHL_6_8
│       ├── 158376_62
│       │   └── snv_out_withsubject.txt
│       └── NHL_6.1
│           └── snv_out_withsubject.txt
└── reporters
    ├── 158376_62
    │   ├── 158376_62_fixed.vcf
    │   ├── 158376_62_germline_annotated.vcf
    │   └── 158376_62_reporters.vcf
    ├── 86755_94
    │   ├── 86755_94_fixed.vcf
    │   ├── 86755_94_germline_annotated.vcf
    │   └── 86755_94_reporters.vcf
    ├── NHL_2.1
    │   ├── NHL_2.1_fixed.vcf
    │   └── NHL_2.1_germline_annotated.vcf
    └── NHL_6.1
        ├── NHL_6.1_fixed.vcf
        └── NHL_6.1_germline_annotated.vcf
```

## Output

Your Python program should generate three output files. A set of examples are provided under `data/output`.

1. QC Metrics

   An comma-delimited aggregated file (`qc_metrics.csv`) based on four `[Sample ID]_qc_metrics_summary.csv` from `Plasma` samples.
   The columns are described below. The output file should contain QC metrics from four `Plasma` samples, each sample should have three rows (one row for each measurement).

   | Column Name | Description                                         | Data Source                  |
   | ----------- | --------------------------------------------------- | ---------------------------- |
   | PATNUM      |                                                     | `Subject ID` in the manifest |
   | VISIT       |                                                     | `Visit Code` in the manifest |
   | ACCSNM      |                                                     | `Sample ID` in the manifest  |
   | PGRUNID     |                                                     | `MRD Type` in the manifest   |
   | PGTESTCD    | Genomics test or examination short name             | see Appendix A               |
   | PGTEST      | Genomics test or examination long name              | see Appendix A               |
   | BMSRESN     | Numeric result or finding in agreed units           | see Appendix A               |
   | BMSORESU    | Unit of measurement associated with the test result | see Appendix A               |

   Appendix A
   | PGTESTCD | PGTEST                             | BMSRESN                                                                                     | BMSORESU |
   | -------- | ---------------------------------- | ------------------------------------------------------------------------------------------- | -------- |
   | DNAIPMS  | DNA Input Mass                     | `InputMass` column in `[Sample ID]_qc_metrics_summary.csv`                                  | ng       |
   | INSPVL   | Volume Measurement, Input Specimen | `PlasmaVolume` column in `[Sample ID]_qc_metrics_summary.csv`                               | mL       |
   | SEQOTRT  | Sequencing On-Target Rate          | `ONTARGET_SORTEDBAM_CAPTURE: On-target rate` column in `[Sample ID]_qc_metrics_summary.csv` | %        |

2. MRD Summary

   An comma-delimited aggregated file (`mrd_summary.csv`) based on all four `snv_out_withsubject.txt`.
   The columns are described below. The output file should contain ctDNA level from all four samples.
   Each `Baseline` plasma sample should have one row that contains the `DCMMNCNC` (mean MMPM (mutant molecule per mL of plasma)) of the sample.
   Each `FollowUp` plasma sample should have two rows (one row of `DCMMNCNC` and one row of `DCMMNLCB`).

   | Column Name | Description                                         | Data Source                  |
   | ----------- | --------------------------------------------------- | ---------------------------- |
   | PATNUM      |                                                     | `Subject ID` in the manifest |
   | VISIT       |                                                     | `Visit Code` in the manifest |
   | ACCSNM      |                                                     | `Sample ID` in the manifest  |
   | PGRUNID     |                                                     | `MRD Type` in the manifest   |
   | PGTESTCD    | Genomics test or examination short name             | see Appendix B               |
   | PGTEST      | Genomics test or examination long name              | see Appendix B               |
   | GNFRESN     | Numeric result or finding in agreed units           | see Appendix B               |
   | GNFORESU    | Unit of measurement associated with the test result | see Appendix B               |

   Appendix B
   | PGTESTCD | PGTEST                                   | GNFRESN                                                                           | GNFORESU  |
   | -------- | ---------------------------------------- | --------------------------------------------------------------------------------- | --------- |
   | DCMMNCNC | Cancer DNA Molecules, Mean Concentration | `MMPM_MEAN` column in `snv_out_withsubject.txt`                                   | copies/mL |
   | DCMMNLCB | Cancer DNA Molec, Mean Conc Log Chng Bsl | Calculated as: log10(`DCMMNCNC` of the `Baseline` / `DCMMNCNC` of the `FollowUp`) |           |

3. SNV

   An comma-delimited aggregated file (`snv.csv`) based on all four `[Sample ID]_germline_annotated.vcf`.
   The columns are described below. The output file should contain SNVs from all four samples.

   | Column Name | Description                                                                   | Data Source                  |
   | ----------- | ----------------------------------------------------------------------------- | ---------------------------- |
   | PATNUM      |                                                                               | `Subject ID` in the manifest |
   | VISIT       |                                                                               | `Visit Code` in the manifest |
   | ACCSNM      |                                                                               | `Sample ID` in the manifest  |
   | PGRUNID     |                                                                               | `MRD Type` in the manifest   |
   | PGTESTCD    | Genomics test or examination short name                                       | see Appendix C               |
   | PGTEST      | Genomics test or examination long name                                        | see Appendix C               |
   | GNFRESN     | Numeric result or finding in agreed units                                     | see Appendix C               |
   | GNFORESU    | Unit of measurement associated with the test result                           | see Appendix C               |
   | PFGRPID     | ID used to group all tests referring to the same set of pharmacogenomics test | see Appendix C               |

   Appendix C
   | PGTESTCD | PGTEST           | GNFRESN                                                               | GNFORESU | PFGRPID                                                             |
   | -------- | ---------------- | --------------------------------------------------------------------- | -------- | ------------------------------------------------------------------- |
   | UQDPTH   | Read Depth       | `DP` in the `INFO` column of the `[Sample ID]_germline_annotated.vcf` |          | [Sample ID_Chromosome:Position_Reference Allele_Alternative Allele] |
   | ALFRQ    | Allele Frequency | `AF` in the `INFO` column of the `[Sample ID]_germline_annotated.vcf` | %        | [Sample ID_Chromosome:Position_Reference Allele_Alternative Allele] |
