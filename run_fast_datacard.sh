#!/bin/sh
# ls /vols/cms/VBFHinv/20210819-UL-NewEGamma/
for dir in test_df_MTR_2017_2020v1 test_df_MTR_2018_2020v1 test_df_VTR_2017_2020v1 test_df_VTR_2018_2020v1;
do
    cd $dir
    fast_datacard SR_cfg.yaml
    fast_datacard Zmumu_cfg.yaml
    fast_datacard Zee_cfg.yaml
    fast_datacard Wmunu_cfg.yaml
    fast_datacard Wenu_cfg.yaml
    cd ../
done
