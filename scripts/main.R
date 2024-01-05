# gono model history matching
# with emulation workflow
#
# see tidy models for inspiration
# https://www.tmwr.org/workflows.html


# create analysis object
gono_hmer <- 
  hmer_analysis() |> 
  add_input_ranges(names, data) |> 
  add_targets(names, data) |> 
  add_model(gono_model)

# run initial wave
gono_hmer@run_wave
gono_hmer <- run_wave(gono_hmer)


# run first wave
gono_hmer@run_wave
gono_hmer <- run_wave(gono_hmer)
