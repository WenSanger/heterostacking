#Define the Fortran compiler
FC = gfortran

# Define the object files
OBJS = solve_mn.o twist_angle.o distance_func.o

# Define the target
solve_mn: $(OBJS)
	$(FC) -o $@ $^

# Rule for compiling .f90 files to .o files
%.o: %.f90
	$(FC) -c $<

# Clean up
clean:
	rm -f solve_mn $(OBJS)

