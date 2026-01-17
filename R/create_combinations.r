# create different model combinations

create_combinations <- function(){
  
  combinations <- list(list(name = "fluA",
                            p_inf = 0.0972,
                            inc_period = 2,
                            inf_period = 5,
                            imm_period = 360), 
                       list(name = "fluB",
                            p_inf = 0.0972,
                            inc_period = 2,
                            inf_period = 6,
                            imm_period = 360),
                       list(name = "RSV",
                            p_inf = 0.0972,
                            inc_period = 4.98,
                            inf_period = 6.16,
                            imm_period = 358.9), 
                       list(name = "hCOV",
                            p_inf = 0.0972,
                            inc_period = 3,
                            inf_period = 3.5,
                            imm_period = 360), 
                       list(name = "AdV",
                            p_inf = 0.0972,
                            inc_period = 6,
                            inf_period = 5.5,
                            imm_period = 360), 
                       list(name = "RV",
                            p_inf = 0.0972,
                            inc_period = 2,
                            inf_period = 11,
                            imm_period = 360), 
                       list(name = "hMPV",
                            p_inf = 0.0972,
                            inc_period = 4,
                            inf_period = 10.5,
                            imm_period = 360))
                       
  return(combinations)
  
}
