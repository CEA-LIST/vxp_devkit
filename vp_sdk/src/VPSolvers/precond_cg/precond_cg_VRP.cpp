/**
* Copyright 2023 CEA Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
* 
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
* 
*     http://www.apache.org/licenses/LICENSE-2.0
* 
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
**/
/**
 * Authors       : Jerome Fereyre
 * Creation Date : August, 2023
 * Description   : 
 **/

#include "precond_cg_kernel.hpp"
#include "VRPSDK/perfcounters/cpu.h"
#include "VRPSDK/vblas_perfmonitor.h"
#include "VPSDK/VBLASConfig.hpp"
#include <time.h>

namespace VPFloatPackage::Solver {

	int precond_cg(int precision, int transpose, int n, double * x, matrix_t A, matrix_t iM, double * b, double tolerance, uint16_t exponent_size, int32_t stride_size, char * log_buffer, uint64_t log_buffer_size) {
			int l_iteration_count;
			
			VBLASPERFMONITOR_INITIALIZE;

			VPFloatPackage::VBLAS::VBLAS_Init();

			clock_t t0, t1;
			uint64_t instr0, instr1;
			uint64_t dmiss0, dmiss1;
			uint64_t imiss0, imiss1;
			dmiss0 = cpu_dmiss();
			imiss0 = cpu_imiss();
			instr0 = cpu_instructions();

			VPFloatArray Xv(x, n);
			VPFloatArray Bv(b, n);

			std::cout << "transpose : " << transpose << std::endl;

			dmiss0 = cpu_dmiss();
			imiss0 = cpu_imiss();
			instr0 = cpu_instructions();
			t0 = clock();

			l_iteration_count = precond_cg_vp(precision, transpose, n, Xv, A, iM, Bv, tolerance, exponent_size, stride_size);
			t1 = clock();
			dmiss1 = cpu_dmiss();
			imiss1 = cpu_imiss();
			instr1 = cpu_instructions();

			double ipc = ((double)(instr1-instr0)/(t1-t0));

			VPFloatPackage::VBLAS::VBLAS_Destroy();

			VBLASPERFMONITOR_DISPLAY;

			std::cout << "PCG: instructions: " << instr1 - instr0 << " / dmiss =" << dmiss1 - dmiss0 << "/ imiss =" << imiss1 - imiss0 << " / ipc = " << ipc << " / elapsed_time=" << (((double)(t1-t0))/CORE_REFCLK) << std::endl;

			return l_iteration_count;			
	}

}