CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)

OS_NAME:=$(shell uname -s | tr A-Z a-z)
ifeq ($(OS_NAME),darwin)
STDINCDIR := -I/opt/local/include
STDLIBDIR := -L/opt/local/lib
else
STDINCDIR :=-I/media/data/cmorgoth/git/DarkMatterAna/
STDLIBDIR := 
endif

CPPFLAGS := -Wl,--no-as-needed $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)

CPPFLAGS += -g

TARGET = DM_ANA_Err
TARGET1 = DataAnaEff
#SRC = src/testV3.cc src/DM_StackPlotsV3.cc src/DM_1DRatioV3.cc src/DM_2DRatioV3.cc src/DM_BaseV3.cc src/DM_DY_HTBinsV3.cc src/DM_TT_LSLHV3.cc src/DM_DataV3.cc src/DM_ZJetsNuNuV3.cc src/DM_WJetsHTBinsV3.cc src/DM_METPlotsV3.cc src/DM_KinePlotsV3.cc

#SRC = src/error_studyV3.cc src/DM_StackPlotsV3.cc src/DM_1DRatioV3.cc src/DM_2DRatioV3.cc src/DM_BaseV3.cc src/DM_DY_HTBinsV3.cc src/DM_TT_LSLHV3.cc src/DM_DataV3.cc src/DM_ZJetsNuNuV3.cc src/DM_WJetsHTBinsV3.cc src/DM_METPlotsV3.cc src/DM_KinePlotsV3.cc

SRC = src/testV3_MRCategories.cc src/DM_StackPlotsV3.cc src/DM_1DRatioV3.cc src/DM_2DRatioV3.cc src/DM_BaseV3.cc src/DM_DY_HTBinsV3.cc src/DM_TT_LSLHV3.cc src/DM_DataV3.cc src/DM_ZJetsNuNuV3.cc src/DM_WJetsHTBinsV3.cc src/DM_METPlotsV3.cc src/DM_KinePlotsV3.cc

SRC1 = simple_data_yield_study.cc

OBJ = $(SRC:.cc=.o)

OBJ1 = $(SRC1:.cc=.o)

all : $(TARGET) $(TARGET1)

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET1) : $(OBJ1)
	$(LD) $(CPPFLAGS) -o $(TARGET1) $(OBJ1) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET) $(TARGET1) *~

