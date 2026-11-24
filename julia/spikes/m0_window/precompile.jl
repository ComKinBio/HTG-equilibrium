# Workload recorded while building the sysimage: a full start-up, every component
# driven by the self-test, a screenshot, a memory report and a clean shut-down.
ENV["GENIE_ENV"] = "prod"
ENV["M0_SELFTEST"] = "1"
ENV["M0_AUTOCLOSE"] = "1"
ENV["M0_IDLE"] = "1"
ENV["M0_SHOT"] = tempname() * ".png"
using M0App
M0App.main()
rm(ENV["M0_SHOT"]; force = true)
