# PPG_quality_determinants

This repository provides code to enable reproduction of the analysis described in:

Charlton PH et al., 'Determinants of photoplethysmography signal quality at the wrist' [under review], https://doi.org/10.36227/techrxiv.172954491.17588920/v1

Take the following steps to reproduce the analysis:
1. You will need the Aurora-BP Dataset. See [here](https://github.com/microsoft/aurorabp-sample-data) for details of how to access it.
   - Save the dataset to a folder on your computer called `raw_data`, and note down the path of this folder (e.g. `/Users/petercharlton/Documents/Data/Aurora/raw_data/`).
2. Collate the dataset using `collate_aurora_data.m`. To do so:
   - Adjust the `root_folder` (on line 15) to be the path of the folder containing the `raw_data` folder (e.g. `/Users/petercharlton/Documents/Data/Aurora/`).
   - Adjust `up.pt_data_filename` to be the path of the `participants.tsv` file (e.g. `/Users/petercharlton/Documents/Data/Aurora/raw_data/participants.tsv`)
   - This will produce many files, corresponding to the different activities performed by participants in the oscillometric and auscultatory subsets. Only some of these files are used in the following analysis.
3. Download the `ppg-quality` toolbox from [here](https://ppg-quality.readthedocs.io/).
   - Add it to the Matlab path.
4. Run `aurora_ppg_quality_assessment2` to perform the analysis. To do so:
   - Set `up.root_data_path` (in the `setup_up` function) to be the same as the `root_folder` specified in step 2 above (e.g. `/Users/petercharlton/Documents/Data/Aurora/`).
   - Set `up.plots_path` (in the `setup_up` function) to a folder where you're happy for the results plots to be saved.
   - Note that Fig. 1 and 4 of the paper won't be reproduced exactly (they were generated using an old version of the code).
5. If you wish, you can reproduce Fig. 2 of the paper ('The AC and DC components of the PPG signal') using `ppg_components_plot.m`.
