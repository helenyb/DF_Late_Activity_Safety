# DF_Late_Activity_Safety

R Code README for “Dose Finding Studies for Therapies with Late-Onset Safety and Efficacy Outcomes”

This folder contains 5 function files:

1. data_generation_TTE.R
Contains functions for data generation required for all methods. This includes functions to find parameters for data generation.

2. compare_Yuan_v4.R
Function to implement one simulation of trial conducted according to Yuan & Yin design

3. CRM_2_combo_v11.R
Function to implement one simulation of trial conducted according to Joint TITE-CRM design.

4. CRM_2_combo_v12.R
Function to implement one simulation of trial conducted according to Joint TITE-CRM design, with added option for no time-to-event.

5. LiuJohnson_compare_v2.R
Function to implement one simulation of trial conducted according to Liu & Johnson design.

And 9 implementation files:

1. TTE_simsv4.R
Simulations for Joint TITE-CRM.

2. TTE_simsv4_noSTOP.R
Simulations for Joint TITE-CRM (minimal stopping rules)

3. TTE_simsv4_noTITE.R
Simulations for Joint CRM.

4. TTE_simsv4_noTITEnoSTOP.R
Simulations for Joint CRM. (minimal stopping rules)

5. TTE_simsLJ_v4.R
Simulations for Liu & Johnson design.

6. TTE_simsLJ_v4_noSTOP.R
Simulations for Liu & Johnson design. (minimal stopping rules)

7. TTE_simsYY_v4.R
Simulations for Yin & Yuan design.

8. TTE_simsYY_v4_noSTOP.R
Simulations for Yin & Yuan design. (minimal stopping rules)

9. eff_patterns.R
Simulations for alternative efficacy time trends.
