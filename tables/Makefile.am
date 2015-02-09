EXTRA_DIST = bufrtabb_13.csv bufrtabb_12.csv bufrtabb_11.csv bufrtabb_6.csv bufrtabb_2.csv \
	  bufrtabd_13.csv bufrtabd_12.csv bufrtabd_11.csv bufrtabd_6.csv bufrtabd_2.csv \
	  localtabb_65535_6.csv localtabb_255_1.csv localtabb_85_1.csv localtabb_85_10.csv \
          localtabb_85_2.csv localtabb_85_0.csv localtabb_85_12.csv localtabb_247_7.csv \
          localtabb_247_8.csv localtabb_85_14.csv bufrtabb_14.csv \
	  localtabd_65535_6.csv localtabd_255_1.csv localtabd_85_1.csv localtabd_85_10.csv \
          localtabd_85_2.csv localtabd_85_0.csv localtabd_85_12.csv localtabd_247_7.csv \
          localtabd_247_8.csv localtabd_85_14.csv bufrtabd_14.csv \
          localtabb_247_9.csv localtabd_247_9.csv bufrtabb_15.csv bufrtabd_15.csv \
          bufrtabb_16.csv bufrtabd_16.csv 

install-data-local:
	test -z "$(prefix)/share/bbufr/tables/" || $(MKDIR_P) "$(prefix)/share/bbufr/tables/"
	@for file in $(EXTRA_DIST); do \
		$(install_sh_DATA) $$file "$(prefix)/share/bbufr/tables/" || exit $$?; \
	done
          