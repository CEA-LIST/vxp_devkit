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
# Author : Riccardo ALIDORI
#

rm -f /usr/share/vrp/_done

# VRP dyamic linking script

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>/usr/share/vrp/dyn_link.log 2>&1
# Everything below will go to the file 'dyn_link.log':

rm /usr/share/vrp/vrp_fw.x*

echo "VRP Malloc Address: 0x$1"
echo "__ADDR_RAM_UNCACHED__: 0x$1"
echo "__ADDR_RAM_CACHED__:  0x$2"
echo "__RAM_SIZE__: $3M"
echo "Opening file: $4"

# keep always the original Linker script
cp /usr/share/vrp/scripts/vrp_linkcmds_template /usr/share/vrp/scripts/vrp_linkcmds
sed -i -e s/__ADDR_RAM_UNCACHED__/0x$1/g /usr/share/vrp/scripts/vrp_linkcmds
sed -i -e s/__ADDR_RAM_CACHED__/0x$2/g /usr/share/vrp/scripts/vrp_linkcmds
sed -i -e s/__RAM_SIZE__/$3M/g /usr/share/vrp/scripts/vrp_linkcmds

# Call the RISCV Linker
riscv64-unknown-elf-g++ -march=rv64imac -mabi=lp64 -mcmodel=medany -mrelax -nostdlib -static -Wl,--gc-sections -Wl,--print-memory-usage -Wl,--relax -T /usr/share/vrp/scripts/vrp_linkcmds -u _printf_float -o /usr/share/vrp/vrp_fw.x $4 -lc -lgcc
riscv64-unknown-elf-objcopy -S -O binary /usr/share/vrp/vrp_fw.x  /usr/share/vrp/vrp_fw.x.bin
riscv64-unknown-elf-objdump -D /usr/share/vrp/vrp_fw.x > /usr/share/vrp/vrp_fw.x.dump

truncate --no-create -s $3M  /usr/share/vrp/vrp_fw.x.bin

sync -f /usr/share/vrp/vrp_fw.x.bin

# This is a flag for the kernel module, which indicates that the compilation is done
touch /usr/share/vrp/_done

exit 0
