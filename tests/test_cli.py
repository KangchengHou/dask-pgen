import os
import subprocess
import pandas as pd
import numpy as np
import tempfile
import dapgen
from .test_score import python_score


def test_score():
    """
    Test cli `dapgen score`
    """

    tmp_dir = tempfile.TemporaryDirectory()
    tmp_dir_name = tmp_dir.name
    # make up score file
    # build df_weight
    plink_path = dapgen.get_test_data("plink2.merged.pgen")
    np.random.seed(1234)
    df_weight = dapgen.read_pvar(plink_path.replace(".pgen", ".pvar"))
    weight_cols = ["WEIGHT1", "WEIGHT2", "WEIGHT3"]
    for col in weight_cols:
        df_weight[col] = np.random.normal(size=len(df_weight))
    df_weight.to_csv(os.path.join(tmp_dir_name, "weight.tsv"), index=False, sep="\t")

    for center in [True, False]:
        cmds = [
            "dapgen score",
            f"--plink {plink_path}",
            f"--weights {os.path.join(tmp_dir_name, 'weight.tsv')}",
            f"--out {os.path.join(tmp_dir_name, 'score.tsv')}",
        ]
        if center:
            cmds.append("--center True")
        subprocess.check_call(" ".join(cmds), shell=True)

        # # full observation
        df_score = pd.read_csv(
            os.path.join(tmp_dir_name, "score.tsv"), index_col=0, sep="\t"
        )
        true_score, true_df_snp = python_score(
            plink_path, df_weight, weight_cols, center=center
        )
        assert np.allclose(df_score.values, true_score)

    tmp_dir.cleanup()

    return
