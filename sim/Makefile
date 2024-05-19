run_comp:
	vcs -sverilog -f tb.f -top tb -R -kdb +lint=PCWM +lint=TFIPC-L -timescale=1ns/1ps -l log -P $(VERDI_HOME)/share/PLI/VCS/LINUX64/verdi.tab $(VERDI_HOME)/share/PLI/VCS/LINUX64/pli.a -R  -debug_access -debug_region+lib+cell  +fsdb+glitch=0 +fsdb+delta
clean:
	@rm -rf csrc simv.daidir simv verdiLog novas.rc verdi_config_file work novas_dump.log log tb.fsdb ucli.key
