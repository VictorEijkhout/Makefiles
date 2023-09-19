################################################################
####
#### Makefile for software bookkeeping
####
################################################################

info ::
	@echo "make clean : non-recursive"
	@echo "make allclean : recursive"

.PHONY: clean allclean
clean ::
	@rm -f *~ *.log
allclean :: clean
	@for d in * ; do \
	  if [ -d "$${d}" ] ; then \
	    echo " .. cleaning $$d" \
	     && ( cd "$${d}" && make --no-print-directory clean ) \
	  ; fi \
	done

