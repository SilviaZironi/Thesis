# FC = gfortran-mp-12 -w -O2
FC = gfortran -w -O2

OPT = -llapack -lblas

RUNO    = m_nrtype.o m_nrutil.o m_constants.o m_nr_sub.o cdbonn.o n4lo500_2017.o ichiral.o m_triton_sub.o m_interpolation.o wigner.o main.o \

all: run

run: $(RUNO)
	${FC} -o run $(RUNO) ${OPT}

#run:    $(RUNO)
#		${FC} ${OPT} -o run $(RUNO) $(FLAGS)

m_nrtype.o          : m_nrtype.f90
					${FC} -c m_nrtype.f90

m_nrutil.o          : m_nrutil.f90
					${FC} -c m_nrutil.f90

m_constants.o       : m_constants.f90
					${FC} -c m_constants.f90

m_nr_sub.o          : m_nr_sub.f90
					${FC} -c m_nr_sub.f90

cdbonn.o            : cdbonn.f
					${FC} -c cdbonn.f

n4lo500_2017.o      : n4lo500_2017.f
					${FC} -c n4lo500_2017.f

ichiral.o           : ichiral.f90
					${FC} -c ichiral.f90

m_triton_sub.o      : m_triton_sub.f90
					${FC} -c m_triton_sub.f90

m_interpolation.o   : m_interpolation.f90
					${FC} -c m_interpolation.f90

wigner.o            : wigner.f90
					${FC} -c wigner.f90

main.o              : main.f90
					${FC} -c main.f90

# For a fresh build 
clean:
	rm -f run *.o *.mod *.MOD
