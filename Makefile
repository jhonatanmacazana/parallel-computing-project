#======================================================================
# Makefile for mixed C/C++ sources (igoticecream 2019)
#======================================================================

#----------------------------------------------------------------------
# Project Structure
#----------------------------------------------------------------------
# TARGET    is the name of the output
# BUILD     is the directory where object files & intermediate files will be placed
# OUT       is the directory where target files will be placed
# SOURCES   is a list of directories containing source code
# INCLUDES  is a list of directories containing header files
# LIBRARIES is a list of directories containing libraries, this must be the top level containing include and lib
# RESOURCES is a list of directories containing other files
# DATA      is a list of directories containing data files
#----------------------------------------------------------------------
TARGET    := $(notdir $(CURDIR))
BUILD     := build
OUT       := out
SOURCES   := src
INCLUDES  := include
LIBRARIES := /usr/local
RESOURCES :=
DATA      :=

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
OUTPUT     := $(CURDIR)/$(OUT)/$(TARGET)
SYMBOLS    := $(CURDIR)/$(BUILD)/$(TARGET).out

CSOURCES   := $(foreach dir,$(SOURCES),$(wildcard $(dir)/*.c))
CXXSOURCES := $(foreach dir,$(SOURCES),$(wildcard $(dir)/*.cpp))
ASMSOURCES := $(foreach dir,$(SOURCES),$(wildcard $(dir)/*.s))

OBJS       := $(patsubst %,$(BUILD)/%.o,$(basename $(CSOURCES)) $(basename $(CXXSOURCES)) $(basename $(ASMSOURCES)))
DEPS       := $(patsubst %,$(BUILD)/%.d,$(basename $(CSOURCES)) $(basename $(CXXSOURCES)) $(basename $(ASMSOURCES)))

INCLUDE    += $(addprefix -I,$(foreach dir,$(INCLUDES), $(wildcard $(dir))))
INCLUDE    += $(addprefix -I,$(foreach dir,$(LIBRARIES),$(wildcard $(dir)/include)))

LDLIBS     += $(addprefix -L,$(foreach dir,$(LIBRARIES),$(wildcard $(dir)/lib)))
LDLIBS     += $(addprefix -L,$(foreach dir,$(LIBRARIES),$(wildcard $(dir)/library)))

#----------------------------------------------------------------------
# Compiler & Linker
#----------------------------------------------------------------------
CC  := gcc-9
CXX := g++-9
AS  := gcc-9
AR  := gcc-ar-9 # or /usr/local/opt/binutils/bin/ar
NM  := gcc-nm-9 # or /usr/local/opt/binutils/bin/nm

#----------------------------------------------------------------------
# Compiler & Linker Flags
#----------------------------------------------------------------------
# CPPFLAGS  C and C++ Compiler Flags
# CFLAGS    C Compiler Flags
# CXXFLAGS  C++ Compiler Flags
# ASFLAGS   ASM Compiler Flags
# LDFLAGS   Linker Flags
#----------------------------------------------------------------------
CPPFLAGS  += $(DEFINES) $(INCLUDE)
CPPFLAGS  += -ggdb
CPPFLAGS  += -g3
CPPFLAGS  += -Og
#CPPFLAGS += -O3
CPPFLAGS  += -fmessage-length=0
CPPFLAGS  += -fsigned-char
CPPFLAGS  += -fPIC

CFLAGS    += $(CPPFLAGS)
CFLAGS    += -Wall
CFLAGS    += -Wextra
CFLAGS    += -Werror
CFLAGS    += -pedantic
CFLAGS    += -pedantic-errors
CFLAGS    += -Wfatal-errors
CFLAGS    += -Wpacked
CFLAGS    += -Winline
CFLAGS    += -Wfloat-equal
CFLAGS    += -Wconversion
CFLAGS    += -Wpointer-arith
CFLAGS    += -Wdisabled-optimization
CFLAGS    += -Wunknown-pragmas
CFLAGS    += -Wno-unused-parameter
CFLAGS    += -Wno-unused-function
CFLAGS    += -Wno-unused-variable

CXXFLAGS  += $(CFLAGS)
CXXFLAGS  += -Weffc++
CXXFLAGS  += -Wfloat-equal
CXXFLAGS  += -Wsign-promo
CXXFLAGS  += -Wmissing-declarations
CXXFLAGS  += -Woverloaded-virtual
CXXFLAGS  += -Wmissing-format-attribute
CXXFLAGS  += -Wold-style-cast
CXXFLAGS  += -Wshadow
CXXFLAGS  += -Wctor-dtor-privacy

ASFLAGS   += $(CPPFLAGS)
ASFLAGS   += -x assembler-with-cpp

# LDFLAGS   += -lmbedtls
# LDFLAGS   += -lmbedx509
# LDFLAGS   += -lmbedcrypto

#----------------------------------------------------------------------
# Compiler & Linker Commands
#----------------------------------------------------------------------
# LINK.o      link object files to binary
# COMPILE.c   compile C source files
# COMPILE.cpp compile C++ source files
#----------------------------------------------------------------------
ifeq ($(strip $(CXXSOURCES)),)
    LD := $(CC)
else
    LD := $(CXX)
endif

DEPFLAGS   += -MT $@
DEPFLAGS   += -MMD
DEPFLAGS   += -MP
DEPFLAGS   += -MF $(BUILD)/$*.d

LINK.o      = $(LD) $(LDFLAGS) $(LDLIBS) $^ -o $@
COMPILE.c   = $(CC) $(C_STANDARD) $(CFLAGS) $(DEPFLAGS) -c $< -o $@
COMPILE.cpp = $(CXX) $(CXX_STANDARD) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@
COMPILE.s   = $(AS) $(ASFLAGS) $(DEPFLAGS) -c $< -o $@
SYMBOLS.out = $(NM) -CSn $@ > $(SYMBOLS)

#----------------------------------------------------------------------
# Special Built-in Target
#----------------------------------------------------------------------
# .SUFFIXES     disable built-in wildcard rules
# .INTERMEDIATE make will treat targets as intermediate files, and delete them
# .PRECIOUS     make will not be deleted after it is no longer needed. Keep objects to speed up recompilation
# .PHONY        make will run this targets unconditionally, regardless of whether a file with that name exists or what its last-modification time is
#----------------------------------------------------------------------
.SUFFIXES:
.INTERMEDIATE:
.PRECIOUS: $(OBJS) $(DEPS)
.PHONY: all clean help

#----------------------------------------------------------------------
# Targets
#----------------------------------------------------------------------
all: $(OUTPUT)

$(OUTPUT): $(OBJS)
	@mkdir -p $(@D)
	$(LINK.o)
	$(SYMBOLS.out)

$(BUILD)/%.o: %.c
	@mkdir -p $(@D)
	$(COMPILE.c)

$(BUILD)/%.o: %.cpp
	@mkdir -p $(@D)
	$(COMPILE.cpp)

$(BUILD)/%.o: %.s
	@mkdir -p $(@D)
	$(COMPILE.s)

run: $(OUTPUT)
	@$<

clean:
	@$(RM) -r $(BUILD) $(OUT)

help:
	@echo available targets: all run clean

-include $(DEPS)