FC = gfortran

OPT = -O3 -llapack -lblas

RUNO    = m_nrtype.o m_nrutil.o m_constants.o m_nr_sub.o m_trimod_sub.o m_interpolation.o main.o \

all: run

run: $(RUNO)
	${FC} ${OPT} -o run $(RUNO)

m_nrtype.o          : m_nrtype.f90
					${FC} -c m_nrtype.f90

m_nrutil.o          : m_nrutil.f90
					${FC} -c m_nrutil.f90

m_constants.o       : m_constants.f90
					${FC} -c m_constants.f90

m_nr_sub.o          : m_nr_sub.f90
					${FC} -c m_nr_sub.f90

m_trimod_sub.o      : m_trimod_sub.f90
					${FC} -c m_trimod_sub.f90

m_interpolation.o   : m_interpolation.f90
					${FC} -c m_interpolation.f90

main.o              : main.f90
					${FC} -c main.f90

# For a fresh build 
clean:
	rm -f *.o *.mod *.MOD
