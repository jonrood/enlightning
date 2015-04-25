SAMRAI        = /Users/someuser/SAMRAI-v3.9.1/install
SRCDIR        = .
SUBDIR        = .
OBJECT        = $(SAMRAI)

RHS           = rhs_hybrid_high_llf.m4

default: main	

include $(OBJECT)/config/Makefile.config

CXX_OBJS      = main.o enlightning.o
F77_OBJS      = rhs.o runge_kutta.o tag_cells.o update_state.o

main:	$(CXX_OBJS) $(F77_OBJS) $(LIBSAMRAIDEPEND)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXX_OBJS) $(F77_OBJS)	\
	$(LIBSAMRAI) $(LDLIBS) -o main

clean-check:
	$(SAMCLEAN)

clean:	clean-check
	$(RM) *.f *.o main

reset:	
	$(RM) enlightning.viz/* enlightning.record/mic-* \
        enlightning.restart/* enlightning_log.txt

mic:	
	cd enlightning.record && make && cd ..

cleanmic:
	$(RM) enlightning.record/mic

cleanall:	clean
	$(MAKE) reset
	$(MAKE) cleanmic

all:	main
	$(MAKE) mic

FORTRAN		= $(SRCDIR)/fortran
M4DIRS		= -DFORTDIR=$(FORTRAN) $(SAMRAI_M4_FLAGS)

rhs.o:	$(FORTRAN)/$(RHS)
	$(M4) $(M4DIRS) $(FORTRAN)/$(RHS) > rhs.f
	$(F77) $(FFLAGS) -c rhs.f -o $@

tag_cells.o:	$(FORTRAN)/tag_cells.m4
	$(M4) $(M4DIRS) $(FORTRAN)/tag_cells.m4 > tag_cells.f
	$(F77) $(FFLAGS) -c tag_cells.f -o $@

runge_kutta.o:	$(FORTRAN)/runge_kutta.m4
	$(M4) $(M4DIRS) $(FORTRAN)/runge_kutta.m4 > runge_kutta.f
	$(F77) $(FFLAGS) -c runge_kutta.f -o $@

update_state.o:	$(FORTRAN)/update_state.m4
	$(M4) $(M4DIRS) $(FORTRAN)/update_state.m4 > update_state.f
	$(F77) $(FFLAGS) -c update_state.f -o $@

