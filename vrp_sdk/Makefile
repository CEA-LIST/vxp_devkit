# Copyright 2023 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##
#  @file   Makefile
#  @author Cesar Fuguet-Tortolero
##

HDRS      =
SRCS      =
CFLAGS    =

ifdef NOVRPMPFR
include Makefile.mpfr
include Makefile.novrp
else
ifdef NOVRP
include Makefile.novrp
else
ifdef RISCVBARECEA
include Makefile.riscv-bare-cea
else
include Makefile.vrp
endif
endif
endif

###
#  set compilation flags
###
ifdef DEBUG_VRP_SDK
CFLAGS   += -DDEBUG=$(DEBUG)
endif
ifdef VBLAS_ENABLE_HWPF
CFLAGS   += -DVBLAS_ENABLE_HWPF=$(VBLAS_ENABLE_HWPF)
endif


INCLUDES += $(STDLIB_INCLUDES)
INCLUDES += -Iinclude
INCLUDES += -I$(MATRIX_SDK_HOME)/include

ifdef MPFR
INCLUDES += -I$(MPFR_PATH)/include
endif

CFLAGS	+= $(INCLUDES) -ggdb

CXXFLAGS +=  $(INCLUDES) -ggdb

###
#  set object files
###
OBJS      = $(subst .s,.o,$(filter %.s,$(SRCS))) \
            $(subst .S,.o,$(filter %.S,$(SRCS))) \
            $(subst .cpp,.o,$(filter %.cpp,$(SRCS))) \
            $(subst .c,.o,$(filter %.c,$(SRCS)))

###
#  rules
###
.PHONY: all
all: $(TARGET).a $(TARGET).size

$(OBJS): $(SRCS) $(HDRS)

$(TARGET).a: $(OBJS)
	$(ECHO) "Build static library: $@"
	$(Q)$(AR) rcs $@ $^

$(TARGET).size: $(TARGET).a
	$(Q)$(SZ) $^ > $@

.PHONY: clean
clean:
	$(ECHO) "Clean all object files"
	$(Q)$(RM) $(OBJS)
	$(Q)$(RM) $(TARGET).a $(TARGET).size

%.o: %.s
	$(ECHO) "Build object: $@"
	$(Q)$(AS) $(ASFLAGS) -o $@ $<

%.o: %.cpp
	$(ECHO) "Build object: $@"
	$(Q)$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.S
	$(ECHO) "Build object: $@"
	$(Q)$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.c
	$(ECHO) "Build object: $@"
	$(Q)$(CC) $(CFLAGS) -c -o $@ $<
