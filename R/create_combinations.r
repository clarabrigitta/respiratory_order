# create different model combinations

# imm_period / imm_sd define a lognormal prior on the duration of immunity:
#   dlnorm(meanlog = log(imm_period), sdlog = imm_sd)
#
# lb/ub hold the bounds for all 13 fitted parameters, in the order used by
# `idx` in fit_model_rcpp.r / fit_hpc.r:
#   1-5  detection rates, 6-10 susceptibility, 11 immunity duration,
#   12   log10(imports),  13   p_inf
#
# lb[11]/ub[11] are hard bounds on the immunity duration. NOTE: they are still
# NOT used by the prior, which truncates at [0, Inf] instead. 

create_combinations <- function(){

  combinations <- list(list(name = "fluA",
                            inc_period = 2,
                            inf_period = 5,
                            # 95%: 274-1945 d. Protection lost over 2-7 yrs via
                            # antigenic drift (Smith 2004; Bedford 2014);
                            # Woolthuis 2017 assume 3.5 yrs. Weakly identified
                            # here: almost no influenza in the 2020-22 window.
                            imm_period = 730,
                            imm_sd = 0.5,
                            lb = c(rep(0, 5), rep(0, 5), 180, -3, 0),
                            ub = c(rep(1, 5), rep(1, 5), 2555, 3, 1)),
                       list(name = "fluB",
                            inc_period = 2,
                            inf_period = 6,
                            # 95%: 411-2918 d. As fluA but longer: B lineages
                            # drift more slowly than A(H3N2) (Bedford 2014).
                            imm_period = 1095,
                            imm_sd = 0.5,
                            lb = c(rep(0, 5), rep(0, 5), 180, -3, 0),
                            ub = c(rep(1, 5), rep(1, 5), 3650, 3, 1)),
                       list(name = "RSV",
                            inc_period = 4.98,
                            inf_period = 6.16,
                            # 95%: 116-457 d. Lang 2022 review: 19 RSV models
                            # use 183-203 d; one uses ~359 d (the old value
                            # here). Hall 1991 rechallenge; Falsey 2006.
                            imm_period = 230,
                            imm_sd = 0.35,
                            lb = c(rep(0, 5), rep(0, 5), 60, -3, 0),
                            ub = c(rep(1, 5), rep(1, 5), 730, 3, 1)),
                       list(name = "hCOV",
                            inc_period = 3,
                            inf_period = 3.5,
                            # 95%: 159-626 d. Kissler 2020 estimate 45 wks;
                            # Edridge 2020 find reinfection from 6 months and
                            # common by 12 months. Best-evidenced of the eight.
                            imm_period = 315,
                            imm_sd = 0.35,
                            lb = c(rep(0, 5), rep(0, 5), 60, -3, 0),
                            ub = c(rep(1, 5), rep(1, 5), 730, 3, 1)),
                       list(name = "AdV",
                            inc_period = 6,
                            inf_period = 5.5,
                            # 95%: 101-719 d. Type-specific antibodies last
                            # years, but "Adenovirus" pools many types, so
                            # protection against the group is shorter. Weak.
                            imm_period = 270,
                            imm_sd = 0.5,
                            lb = c(rep(0, 5), rep(0, 5), 30, -3, 0),
                            ub = c(rep(1, 5), rep(1, 5), 1095, 3, 1)),
                       list(name = "RV",
                            inc_period = 2,
                            inf_period = 11,
                            # 95%: 34-240 d. ~160 types; repeat infections are
                            # mostly with a different type, same-species
                            # reinfection at ~55-68 d in a child cohort.
                            imm_period = 90,
                            imm_sd = 0.5,
                            lb = c(rep(0, 5), rep(0, 5), 14, -3, 0),
                            ub = c(rep(1, 5), rep(1, 5), 730, 3, 1)),
                       list(name = "hMPV",
                            inc_period = 4,
                            inf_period = 10.5,
                            # 95%: 114-548 d. Little direct human evidence;
                            # borrowed from RSV (same family), with lifelong
                            # reinfection and little protection at 58 wks in
                            # an animal challenge study.
                            imm_period = 250,
                            imm_sd = 0.4,
                            lb = c(rep(0, 5), rep(0, 5), 60, -3, 0),
                            ub = c(rep(1, 5), rep(1, 5), 730, 3, 1)),
                       list(name = "PIV",
                            inc_period = 2.6,
                            inf_period = 10,
                            # 95%: 91-438 d. Glezen 1984 (Houston Family
                            # Study): >=2/3 of children infected with PIV3 in
                            # each of their first two years. Data pool all
                            # four types, so group-level protection is shorter.
                            imm_period = 200,
                            imm_sd = 0.4,
                            lb = c(rep(0, 5), rep(0, 5), 30, -3, 0),
                            ub = c(rep(1, 5), rep(1, 5), 730, 3, 1)))

  return(combinations)

}
