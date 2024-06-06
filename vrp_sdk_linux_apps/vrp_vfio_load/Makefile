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
# Authors       : Riccardo ALIDORI
# Creation Date : August, 2022
# Description   : 
APP = vrpinit

# Add any other object files to this list below
SRCS=$(wildcard *.c)
APP_OBJS = $(SRCS:.c=.o)
CFLAGS += -ggdb -O0
LDFLAGS += -ldl -lpthread

all: build

build: $(APP)

$(APP): $(APP_OBJS)
	$(CC) -o $@ $(APP_OBJS) $(LDFLAGS) $(LDLIBS)
clean:
	rm -f $(APP) *.o
