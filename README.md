# scWT
The scripts for scWT analysis:

First run the shell script in pre-processing to handle the raw data.

Then scWT_fig_xx.R could be uesd for figure generating.


scWT_fig_01.R is used for fig1 generating.

scWT_fig_02.R is used for fig2 generating.

scWT_fig_03.R is used for fig3 generating.

```bash
├── pre-processing
│   ├── run_scWT_bulk.sh
│   ├── run_scWT_drop.sh
│   └── sub-script
│       ├── BCFilter.py
│       ├── calcCells_other.R
│       ├── calcCells.R
│       ├── calcCov.py
│       ├── calcHybrid.py
│       ├── calcReads.py
│       ├── calcTrans.py
│       └── calcTrans.R
├── scWT_fig_01.R
├── scWT_fig_02.R
└── scWT_fig_03.R
```
