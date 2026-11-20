#!/bin/bash

jobid_pre=$(sbatch /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/001_vascular_reintegrate_preprocess.sh | awk '{print $4}')

jobid_harm=$(sbatch --dependency=afterok:$jobid_pre /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/001_vascular_reintegrate_harmony.sh | awk '{print $4}')
jobid_scANVI=$(sbatch --dependency=afterok:$jobid_pre /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/001_vascular_reintegrate_scANVI.sh | awk '{print $4}')
jobid_scPoli=$(sbatch --dependency=afterok:$jobid_pre /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/001_vascular_reintegrate_scPoli.sh | awk '{print $4}')

jobid_bench=$(sbatch --dependency=afterok:$jobid_harm,$jobid_scANVI,$jobid_scPoli --array=1-4 /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/001_vascular_integration_evaluation.sh | awk '{print $4}')

sbatch --dependency=afterok:$jobid_bench /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/001_vascular_compile_benchmarking.sh