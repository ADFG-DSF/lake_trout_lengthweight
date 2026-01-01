{
  {
    script1start <- Sys.time()
    print(script1start)
    source("R/2_laketrout_lwmodels.R", print.eval = FALSE)
    aa_script1time <- Sys.time() - script1start
    print(aa_script1time)
  }

  {
    script2start <- Sys.time()
    print(script2start)
    source("R/2pt1_lwmodels_loocv.R", print.eval = FALSE)
    aa_script2time <- Sys.time() - script2start
    print(aa_script2time)
  }
}
