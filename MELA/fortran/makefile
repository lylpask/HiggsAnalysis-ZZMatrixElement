ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

CMSROOT = ./
#MODULES=$(CMSROOT)modules/
#OBJECTS=$(CMSROOT)objects/
MODULES=$(CMSROOT)
OBJECTS=$(CMSROOT)
ROOFITINCLUDE = 
RM  = /bin/rm
INCLUDE = -I$(CMSROOT) $(ROOFITINCLUDE)
Comp = gfort

# Linking LHAPDF
# link pdfs via LHA library ('Yes' or 'No')
UseLHAPDF=No
# directory which contains libLHAPDF.a, libLHAPDF.la, libLHAPDF.so
MyLHADir=/home/schulze/lib/LHAPDF-6.1.5/lib/
# remember to export
#          LD_LIBRARY_PATH=/.../LHAPDF-x.y.z/lib/:${LD_LIBRARY_PATH}
#          LHAPDF_DATA_PATH=/.../LHAPDF-x.y.z/share/LHAPDF/:${LHAPDF_DATA_PATH}
ifeq ($(UseLHAPDF),Yes)
   LHAPDFflags = -L$(MyLHADir) -lLHAPDF -DuseLHAPDF=1
else
   LHAPDFflags = -DuseLHAPDF=0
endif


ifeq ($(Comp),ifort)
	fcomp = ifort -fpp -O2 -vec-report0 -Dcompiler=1 -lifcore -fPIC
endif
ifeq ($(Comp),gfort)
	fcomp = gfortran -O0 -ffree-line-length-none -Dcompiler=2 $(LHAPDFflags) -fno-automatic -fno-f2c  -fPIC -g 
endif

ifeq ($(Comp),gfort)
	ccomp = gcc
endif
ifeq ($(Comp),ifort)
	ccomp = icc
endif

ifeq ($(Comp),gfort)
	clib =  -lm -lgfortran
endif
ifeq ($(Comp),ifort)
	clib =  -lm -lirc
endif


LIB = libjhugenmela.so

all:    libs
	g++ -Wl,-soname,$(LIB) -shared -o $(LIB) *.o
#	mv *.o $(OBJECTS)
#	mv *.mod $(MODULES)

libs: NNPDFDriver.o mod_Parameters.o mod_Misc.o mod_Kinematics.o mod_JHUGen.o mod_Higgs_MatEl.o mod_Graviton_MatEl.o mod_Zprime_MatEl.o mod_HiggsJJ_MatEl.o mod_HiggsJ_MatEl.o mod_VHiggs_MatEl.o mod_TTBH_MatEl.o mod_TH_MatEl.o

NNPDFDriver.o: ./pdfs/NNPDFDriver.f
	@echo " "
	@echo " compiling NNPDFDriver.f with "$(Comp)
	$(fcomp) -c ./pdfs/NNPDFDriver.f
	
	
vegas.o: ./vegas.f
	@echo " "
	@echo " compiling vegas.f with "$(Comp)
	$(fcomp) -c ./vegas.f


mod_Parameters.o: mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Parameters.F90 with "$(Comp)
	$(fcomp) -c mod_Parameters.F90

mod_Misc.o: mod_Misc.F90
	@echo " "
	@echo " compiling mod_Misc.F90 with "$(Comp)
	$(fcomp) -c mod_Misc.F90

mod_Kinematics.o: mod_Kinematics.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Kinematics.F90 with "$(Comp)
	$(fcomp) -c mod_Kinematics.F90

mod_JHUGen.o: mod_JHUGen.F90 mod_Kinematics.F90 mod_Parameters.F90 mod_Misc.F90
	@echo " "
	@echo " compiling mod_JHUGen.F90 with "$(Comp)
	$(fcomp) -c mod_JHUGen.F90


mod_Higgs_MatEl.o: mod_Higgs_MatEl.F90 mod_Misc.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Higgs_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_Higgs_MatEl.F90


mod_Zprime_MatEl.o: mod_Zprime_MatEl.F90 mod_Misc.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Zprime_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_Zprime_MatEl.F90

mod_Graviton_MatEl.o: mod_Graviton_MatEl.F90 mod_Misc.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Graviton_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_Graviton_MatEl.F90

mod_HiggsJ_MatEl.o: mod_HiggsJ_MatEl.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_HiggsJ_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_HiggsJ_MatEl.F90

mod_HiggsJJ_MatEl.o: mod_HiggsJJ_MatEl.F90 mod_Misc.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_HiggsJJ_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_HiggsJJ_MatEl.F90

mod_VHiggs_MatEl.o: mod_VHiggs_MatEl.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_VHiggs_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_VHiggs_MatEl.F90 

mod_TTBH_MatEl.o: mod_TTBH_MatEl.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_TTBH_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_TTBH_MatEl.F90  

mod_TH_MatEl.o: mod_TH_MatEl.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_TH_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_TH_MatEl.F90  

clean:
	@echo " deleting object files"
	rm -f *.so $(OBJECTS)*.o $(MODULES)*.mod

endif

# supresses command calls
.SILENT:
