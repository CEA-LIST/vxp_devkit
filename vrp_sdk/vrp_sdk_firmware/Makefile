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
# 
# 
# Authors       : Jerome Fereyre
# Creation Date : August, 2023
# Description   : 

VRP_SDK_FIRMWARE		?= cg
TARGET		        	= vrp_solver_${VRP_SDK_FIRMWARE}
BUILD_DIR=$(shell readlink -f ./build)
OBJS=${BUILD_DIR}/${TARGET}_firmware.o $(BUILD_DIR)/utils/cpu_tool.o

CFLAGS=$(shell pkg-config --cflags vrp_sdk_vrp_${BSP}) -I./utils -ggdb -O0 -Wall
CXXFLAGS=$(CFLAGS)
LDFLAGS= $(shell pkg-config --libs vrp_sdk_vrp_${BSP})

all: ${TARGET}

clean: 
	-rm -Rf $(BUILD_DIR) $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) ${CXXFLAGS} -o $@ $^  $(LDFLAGS) 

$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/utils/%.o: utils/%.c
	mkdir -p $(BUILD_DIR)/utils
	$(CC) $(CFLAGS) -c -o $@ $<

bin: $(TARGET)
	$(OC) -S -O binary $(TARGET) $(TARGET).x.bin

mem: bin
	$(BIN2MEM) --in $(TARGET).x.bin --out $(TARGET).x.mem --outw 64

dump: $(TARGET)
	$(OD) -D $< > $<.x.dump	