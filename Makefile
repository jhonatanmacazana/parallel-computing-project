TARGET    := $(notdir $(CURDIR))
BUILD_S     := build/seq
BUILD_P     := build/par
OUT_S       := out/seq
OUT_P       := out/par
SOURCES_S   := src/sequential
SOURCES_P   := src/parallel
INCLUDES    := include


#----------------------------------------------------------------------
# Languages Standard
#----------------------------------------------------------------------
C_STANDARD   := -std=c11
CXX_STANDARD := -std=c++11

#----------------------------------------------------------------------
# Defined Symbols
#----------------------------------------------------------------------
DEFINES += -DDEBUG
DEFINES += -DENABLE_LOGGING

#----------------------------------------------------------------------
# Sources & Files
#----------------------------------------------------------------------
OUTPUT_S    := $(CURDIR)/$(OUT_S)/$(TARGET)
OUTPUT_P    := $(CURDIR)/$(OUT_P)/$(TARGET)
SYMBOLS_S   := $(CURDIR)/$(BUILD_S)/$(TARGET).out
SYMBOLS_P   := $(CURDIR)/$(BUILD_P)/$(TARGET).out

CSOURCES   := $(foreach dir,$(SOURCES),$(wildcard $(dir)/*.c))
CXXSOURCES := $(foreach dir,$(SOURCES),$(wildcard $(dir)/*.cpp))

OBJS       := $(patsubst %,$(BUILD)/%.o,$(basename $(CSOURCES)) $(basename $(CXXSOURCES)) $(basename $(ASMSOURCES)))
DEPS       := $(patsubst %,$(BUILD)/%.d,$(basename $(CSOURCES)) $(basename $(CXXSOURCES)) $(basename $(ASMSOURCES)))


CXXFLAGS = -O3 -Wall 
###-fopenmp

.PHONY: all clean help

FLAGS := -lm -fopenmp

#----------------------------------------------------------------------
# Targets
#----------------------------------------------------------------------

all: seq

seq: ${SOURCES_S}/main.cpp
	@mkdir -p $(OUT_S)
	g++ $(FLAGS) -g ${SOURCES_S}/main.cpp -o ${OUT_S}/${TARGET}

seq-debug: ${SOURCES_S}/main.cpp
	@mkdir -p $(OUT_S)
	g++ $(FLAGS) -DDEBUG -g ${SOURCES_S}/main.cpp -o ${OUT_S}/${TARGET}

seq-output: ${SOURCES_S}/main.cpp
	@mkdir -p $(OUT_S)
	g++ $(FLAGS) -DEXPORT -g ${SOURCES_S}/main.cpp -o ${OUT_S}/${TARGET}

par: ${SOURCES_P}/main.cpp
	@mkdir -p $(OUT_P)
	g++ $(FLAGS) -g ${SOURCES_P}/main.cpp -o ${OUT_P}/${TARGET}

par-debug: ${SOURCES_P}/main.cpp
	@mkdir -p $(OUT_P)
	g++ $(FLAGS) -DDEBUG -g ${SOURCES_P}/main.cpp -o ${OUT_P}/${TARGET}

par-output: ${SOURCES_P}/main.cpp
	@mkdir -p $(OUT_P)
	g++ $(FLAGS) -DEXPORT -g ${SOURCES_P}/main.cpp -o ${OUT_P}/${TARGET}


run-seq: ${OUTPUT_S}
	@$<

run-par: ${OUTPUT_P}
	@$<


clean:
	@$(RM) -r $(BUILD) $(OUT)

help:
	@echo available targets: all {seq,par}-{test,debug,output} run-{seq,par} clean