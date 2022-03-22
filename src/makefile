# Specify which compiler to use:
CXX= g++

# Use CXXFLAGS to specify additional compiler options to use:
CXXFLAGS = -Wall -std=c++0x

# Define variables:
EXE_FILE= cpscdexe 
OBJ_FILES= main.o slipSystem.o deformationGradient.o SCDWrapper.o Bundle.o CascadeDamage.o Damage.o Object.o OneLine.o rvgs.o cpdf.o
HEADER_FILES= slipSystem.h deformationGradient.h constant.h SCDWrapper.h Bundle.h CascadeDamage.h Damage.h Object.h OneLine.h rvgs.h cpdf.h gnuplot-iostream.h gnuplot_i.h

# Build executable:
exe : $(EXE_FILE)
$(EXE_FILE) : $(OBJ_FILES) 
	$(CXX) -o $(EXE_FILE) $(OBJ_FILES) -g -lm -lboost_iostreams -lboost_system -lboost_filesystem
	# chmod ugo+x $(EXE_FILE)	# turn on execute permissionss

# Rules:
.cpp.o:	$(HEADER_FILES)
	$(CXX) $(CXXFLAGS) -c  $*.cpp 

clean:
	rm $(EXE_FILE) $(OBJ_FILES)
