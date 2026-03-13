################################################################
####
#### Makefile for software bookkeeping
####
################################################################

info ::
	@echo "make clean : non-recursive"
	@echo "make allclean : recursive"

.PHONY: clean localclean
localclean : 
	@rm -f *~ *.log install*.o* install.slurm
clean :: localclean
	@for d in * ; do \
	  if [ -d "$${d}" ] ; then \
	    echo " .. cleaning $$d" \
	     && ( cd "$${d}" \
	         && if [ -f Configuration ] ; then \
	                mpm.py clean \
	            ; else \
	                make --no-print-directory clean \
	            ; fi \
	         ) \
	  ; fi \
	done

