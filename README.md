# SMACostEff
Cost-effectiveness model for the treatment of SMA with Zolgensma and Spinraza in The Netherlands

Before running the model using script "Run/ZolgCostEffNL.R", run the script package_install.R such that all required packages will be installed.

## Directory overview
Below are the contents specified of all directories. For now, we only included the files relevant for the basecase and scenario analyses. We have not yet explored the files of the DSA and PSA.

### Run
Includes all scripts to run the model in R:
- ZolgCostEffNL.R: Runs the basecase, only BSC patients can relapse
- ZolgCostEffNL_relapse.R: Runs the scenario analyses, all patients can relapse, different probabilities for relapse
- microDSA_range.R: DSA
- microPSA_FULL_multicore.R: PSA

### Scripts
Includes all R scripts with parameters
#### Transitions
- p_1D.R: SMA1 --> Death (Treatment) 
- P_1D_BSC.R: SMA 1 --> Death (BSC) 
- p_2D.R: SMA2 --> Death 
- lifetable_NL.R: SMA 3 --> Death
- P_10.R: SMA 1 --> SMA 0 (Relapse, Progression-free survival - Overall survival = Event curve)
#### Input parameters for different treatments and BSC
- Zolginput_NL.R: Input parameters for Zolgensma (OA)
- spininput_NL.R: Input parameters for Spinraza (N)
- BSC_inputNL.R: Input parameters for BSC

### Files 
Includes all scripts and files relevant for the survival models/functions
