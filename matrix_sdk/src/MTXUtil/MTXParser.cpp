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

#include "MTXUtil/MTXParser.hpp"

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstddef>

/******************************************************************************
 * Internal structure definitions
 *****************************************************************************/
// Structure defined to map memory to a complex number
typedef struct zcomplex {
  double real;
  double imag;
} zcomplex_t;

// Structure describing matrix cell meta information
typedef struct cell_meta {
  int col_index;
  double * value;
  struct cell_meta * next;
} cell_meta_t;

// Structure describing matrix row meta information
typedef struct row_meta {
  int row_index;
  int nb_elements_in_row;
  struct row_meta * next;
  cell_meta_t * first_cell;
  cell_meta_t * last_cell;
} row_meta_t;

// Structure describing matrix meta information
typedef struct mtx_meta {
  int non_empty_row_count;
  int nb_rows;
  int nb_cols;
  int nb_non_zero_elements; 
  bool complex;
  row_meta_t * row_metas;
  row_meta_t * last_row_meta;
} mtx_meta_t;

/******************************************************************************
 * Internal functions definitions
 *****************************************************************************/
/*
 * Initialization function for the matrix meta information structure
 */
void mtx_meta_init(mtx_meta_t * a_mtx_meta);

/*
 * Function release memory used by matrix meta information structure and other linked structures.
 */
void mtx_meta_clear(mtx_meta_t * a_mtx_meta);

/*
 * Function used to linked a cell meta information structure to a row
 */
bool row_meta_add_cell(mtx_meta_t * a_mtx_meta, row_meta_t * a_row_meta, cell_meta_t * a_cell_meta, int a_col_index, double * a_real);

/*
 * Function that returns the structure containing meta information for a desired row of the matrix.
 * If needed the row meta information structure is created and initialized. This occures on first row retrieve.
 */
row_meta_t * mtx_meta_get_row_index(mtx_meta_t * a_mtx_meta, int a_new_index);

/*
 * Debug function used to display matrix meta information content
 */
void display_mtx_meta(mtx_meta_t * a_mtx_meta);

/*
 * Function used to generate the CSR matrix structure from matrix meta information.
 */
crs_t * build_crs_matrix_from_meta(mtx_meta_t * a_mtx_meta);


/******************************************************************************
 * Internal functions implementation
 *****************************************************************************/

void mtx_meta_init(mtx_meta_t * a_mtx_meta) {
    a_mtx_meta->nb_rows = 0;
    a_mtx_meta->nb_cols = 0;
    a_mtx_meta->nb_non_zero_elements = 0;
    a_mtx_meta->complex = false;
    a_mtx_meta->non_empty_row_count = 0;
    a_mtx_meta->row_metas = NULL;
    a_mtx_meta->last_row_meta = NULL;
}

void mtx_meta_clear(mtx_meta_t * a_mtx_meta) {
    row_meta_t * l_current_row_meta = a_mtx_meta->row_metas;
    row_meta_t * l_next_row_meta = NULL;

    while ( l_current_row_meta != NULL) {
        l_next_row_meta = l_current_row_meta->next;
        free(l_current_row_meta);
        l_current_row_meta = l_next_row_meta;
    }
}

void display_mtx_meta(mtx_meta_t * a_mtx_meta) {
  row_meta_t * l_current_row_meta = a_mtx_meta->row_metas;

  while ( l_current_row_meta != NULL ) {
    printf("%d => count : %d\n",
        l_current_row_meta->row_index,
        l_current_row_meta->nb_elements_in_row);
    l_current_row_meta = l_current_row_meta->next;
  }

  printf("real row count : %d\n", a_mtx_meta->non_empty_row_count);
}

bool row_meta_add_cell(mtx_meta_t * a_mtx_meta, row_meta_t * a_row_meta, cell_meta_t * a_cell_meta, int a_col_index, double * a_value) {
    if ( a_row_meta->first_cell == NULL ) {                                     // Insert first element in the list
        a_cell_meta->col_index = a_col_index;
        a_cell_meta->value = a_value;
        a_cell_meta->next = NULL;
        a_row_meta->first_cell = a_cell_meta;
        a_row_meta->last_cell = a_cell_meta;
    } else if ( a_col_index > a_row_meta->last_cell->col_index ) {   // Insert at end
        a_cell_meta->col_index = a_col_index;
        a_cell_meta->value = a_value;
        a_row_meta->last_cell->next = a_cell_meta;
        a_cell_meta->next = NULL;
        a_row_meta->last_cell = a_cell_meta;
    } else if ( a_col_index < a_row_meta->first_cell->col_index ) {  // Insert at begining
        a_cell_meta->col_index = a_col_index;
        a_cell_meta->value = a_value;
        a_cell_meta->next = a_row_meta->first_cell;
        a_row_meta->first_cell = a_cell_meta;
    }  else {                                                                    // Insert in list
        cell_meta_t * l_current_cell = a_row_meta->first_cell;
        cell_meta_t * l_previous_cell = NULL;

        // Walk through the list to find the right place
        while ( l_current_cell->col_index < a_col_index ) {
            l_previous_cell = l_current_cell;
            l_current_cell = l_current_cell->next;
        }

        if ( l_current_cell->col_index == a_col_index ) {
            printf("Warning: cell(%d,%d) already registered. Ignore this occurence\n", a_row_meta->row_index, l_current_cell->col_index);
            return false;
        } else {
            a_cell_meta->col_index = a_col_index;
            a_cell_meta->value = a_value;
            l_previous_cell->next = a_cell_meta;
            a_cell_meta->next = l_current_cell;
        }
    }
    


    a_row_meta->nb_elements_in_row++;
    a_mtx_meta->nb_non_zero_elements++;

    return true;
}

row_meta_t * create_new_row_meta(int a_new_index) {
    row_meta_t * l_new_row_meta = (row_meta_t *) malloc(sizeof(row_meta_t));

    l_new_row_meta->row_index = a_new_index;
    l_new_row_meta->nb_elements_in_row = 0;
    l_new_row_meta->next = NULL;
    l_new_row_meta->first_cell = NULL;
    l_new_row_meta->last_cell = NULL;

    return l_new_row_meta;
}


row_meta_t * mtx_meta_get_row_index(mtx_meta_t * a_mtx_meta, int a_new_index) {
    row_meta_t * l_row_meta = NULL;

    if ( a_mtx_meta->row_metas == NULL || a_new_index > a_mtx_meta->last_row_meta->row_index ) { 
        // This row index is greater than all the previous ones
        // Add a new row_meta
        l_row_meta = create_new_row_meta(a_new_index);

        if ( a_mtx_meta->last_row_meta != NULL ) {
            a_mtx_meta->last_row_meta->next = l_row_meta;
        }
        a_mtx_meta->last_row_meta = l_row_meta;

        if ( a_mtx_meta->row_metas == NULL ) {
            a_mtx_meta->row_metas = l_row_meta;
        }

        a_mtx_meta->non_empty_row_count++; // We register a new row

    } else {
        row_meta_t * l_current_row_meta = a_mtx_meta->row_metas;
        row_meta_t * l_previous_row_meta = NULL;

        // Iterate to find the right place for insertion in the row metadata list
        while ( l_current_row_meta != NULL && l_current_row_meta->row_index < a_new_index ) {        // Current element describes a row with smaller index than the searched one
            l_previous_row_meta = l_current_row_meta;
            l_current_row_meta = l_current_row_meta->next;
        }

        if ( l_current_row_meta == NULL || l_current_row_meta->row_index > a_new_index ) {
            l_row_meta = create_new_row_meta(a_new_index);

            l_row_meta->next = l_current_row_meta;

            if ( l_previous_row_meta != NULL ) {
                l_previous_row_meta->next = l_row_meta;
            }
            a_mtx_meta->non_empty_row_count++; // We register a new row
        }

        if ( l_current_row_meta->row_index == a_new_index ) { // We have the right element
            l_row_meta = l_current_row_meta;
        }
    }

    return l_row_meta;
}

void MTXParser::displayCSR(crs_t * SPM) {
    int l_row_index, l_col_index;

    printf("SPM->M = %d\n", SPM->M);
    printf("SPM->N = %d\n", SPM->N);
    printf("SPM->NZNUM = %d\n", SPM->NZNUM);

    printf("SPM->ROW_PTR = ");
    for(l_row_index = 0; l_row_index < SPM->M; l_row_index++) {
        printf("%d,", SPM->ROW_PTR[l_row_index]);
    }
    printf("%d\n", SPM->ROW_PTR[l_row_index]);
    
    printf("SPM->COL_IND = ");
    for(l_col_index = 0; l_col_index < SPM->NZNUM; l_col_index++) {
        printf("%d,", SPM->COL_IND[l_col_index]);
    }
    printf("\n");

    for(l_row_index = 0; l_row_index < SPM->M; l_row_index++) {
        int l_row_index_1_based = l_row_index + 1;
        for (l_col_index = SPM->ROW_PTR[l_row_index]; l_col_index < SPM->ROW_PTR[l_row_index+1]; l_col_index++) {
            if ( SPM->VALUE_TYPE == COMPLEX ) {
                printf("%d %d %lf %lf\n", l_row_index_1_based, SPM->COL_IND[l_col_index-1], SPM->VAL[(l_col_index-1)*2], SPM->VAL[(l_col_index-1)*2+1]);
            } else {
                printf("%d %d %lf\n", l_row_index_1_based, SPM->COL_IND[l_col_index-1], SPM->VAL[l_col_index-1]);
            }
        }
    }

}

crs_t * build_crs_matrix_from_meta(mtx_meta_t * a_mtx_meta) {
    crs_t * l_crs = (crs_t *)malloc(sizeof(crs_t));
    row_meta_t * l_current_row_meta = a_mtx_meta->row_metas;
    cell_meta_t * l_current_cell_meta  = NULL;
    int l_row_ptr_index = 0;
    int l_cell_index_in_row = 0;
    int l_current_csr_col_val_start_index = 0;

    if ( a_mtx_meta->complex ) {
        std::cout << "Set VALUE_TYPE to COMPLEX" << std::endl;
        l_crs->VALUE_TYPE = COMPLEX;
    } else {
        std::cout << "Set VALUE_TYPE to REAL" << std::endl;
        l_crs->VALUE_TYPE = REAL;
    }
    l_crs->M = a_mtx_meta->nb_rows;
    l_crs->N = a_mtx_meta->nb_cols;
    l_crs->NZNUM = a_mtx_meta->nb_non_zero_elements;
    l_crs->ROW_PTR = (int *) malloc((l_crs->M+1) * sizeof(int));
    l_crs->COL_IND = (int *) malloc(l_crs->NZNUM * sizeof(int));

    if ( l_crs->VALUE_TYPE == COMPLEX ){
        l_crs->VAL = (double *) malloc(l_crs->NZNUM * sizeof(double) * 2);
    } else {
        l_crs->VAL = (double *) malloc(l_crs->NZNUM * sizeof(double));
    }

    l_crs->ROW_PTR[l_row_ptr_index++] = 1;

    // Walk through the list of rows and fill the CSR accordingly
    while (l_row_ptr_index < l_crs->M+1) {
        // If needed add entries for empty lines in ROW_PTR
        if ( l_current_row_meta == NULL || l_row_ptr_index < l_current_row_meta->row_index ) {
            printf("%d empty line : %d\n", l_crs->ROW_PTR[l_row_ptr_index-1], l_row_ptr_index);
            l_crs->ROW_PTR[l_row_ptr_index] = l_crs->ROW_PTR[l_row_ptr_index-1];
        } else {
            // After empty line process next valid line
            l_crs->ROW_PTR[l_row_ptr_index] = l_crs->ROW_PTR[l_row_ptr_index - 1] + l_current_row_meta->nb_elements_in_row;

            // Walk through the list of values an store them to the right place in CSR
            l_current_cell_meta = l_current_row_meta->first_cell;
            l_cell_index_in_row = 0;

            while ( l_current_cell_meta != NULL ) {
                l_crs->COL_IND[l_current_csr_col_val_start_index + l_cell_index_in_row ] = l_current_cell_meta->col_index;

                if ( l_crs->VALUE_TYPE == COMPLEX ) {
                    zcomplex_t * l_complex = (zcomplex_t *)(l_current_cell_meta->value);
                    l_crs->VAL[(l_current_csr_col_val_start_index * 2) + ( l_cell_index_in_row * 2 )] = l_complex->real;
                    l_crs->VAL[(l_current_csr_col_val_start_index * 2) + ( l_cell_index_in_row * 2 ) + 1] = l_complex->imag;
                } else {
                    l_crs->VAL[l_current_csr_col_val_start_index  + l_cell_index_in_row] = *(l_current_cell_meta->value);
                }
                
                l_current_cell_meta = l_current_cell_meta->next;
                l_cell_index_in_row++;
            }

            l_current_csr_col_val_start_index += l_current_row_meta->nb_elements_in_row;
            l_current_row_meta = l_current_row_meta->next;
        }
        l_row_ptr_index++;
    }

    // If needed add entries for empty lines in ROW_PTR
    for(l_row_ptr_index; l_row_ptr_index < l_crs->M+1; l_row_ptr_index++ ) {
        l_crs->ROW_PTR[l_row_ptr_index] = l_crs->ROW_PTR[l_row_ptr_index-1];
    }

    return l_crs;
}


crs_t  * MTXParser::parseFileToCRS(std::string a_mtx_file_path) {
    FILE * l_mtx_file;
    MM_typecode l_matcode;
    
    crs_t * l_crs = NULL; 
    int l_nb_fields;
    int l_mtx_element_counter = 0 ;
    int l_row_index_1_based, l_col_index_1_based;
    mtx_meta_t l_mtx_meta;
    row_meta_t * l_current_row_meta = NULL, * l_current_symmetric_row_meta = NULL;
    cell_meta_t * l_cell_meta_pool = NULL;
    int l_cell_meta_index = 0;

    double * l_value_pool = NULL;
    int l_nb_cells = 0, l_nb_cells_in_mtx_file=0;
    int l_value_index = 0;
    bool l_row_added = false;

    mtx_meta_init(&l_mtx_meta);

    if ((l_mtx_file = fopen(a_mtx_file_path.c_str(), "r")) == NULL) {
        fprintf(stderr,"could not find this file %s\n", a_mtx_file_path.c_str());
        return NULL;
    }

    if (mm_read_banner(l_mtx_file, &l_matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        return NULL;
    }

    if (mm_read_mtx_crd_size(l_mtx_file, &l_mtx_meta.nb_rows, &l_mtx_meta.nb_cols, &l_nb_cells_in_mtx_file)  !=0 ) {
        return NULL;
    }

    l_nb_cells = l_nb_cells_in_mtx_file;
    
    if ( mm_is_symmetric(l_matcode) ) {
        fprintf(stderr, " symmetric ");
        l_nb_cells *= 2; 
    }

    // Allocate all cell metadata strucutures needed
    l_cell_meta_pool = (cell_meta_t *) malloc(sizeof(cell_meta_t) * l_nb_cells);

    // Allocate memory to store all values
    if (mm_is_complex(l_matcode) && mm_is_matrix(l_matcode) && mm_is_sparse(l_matcode) ) {        
        fprintf(stderr, " complex ");
        l_mtx_meta.complex = true;
        l_value_pool = (double *)malloc(sizeof(zcomplex_t) * l_nb_cells);    
    } else {
        l_value_pool = (double *)malloc(sizeof(double) * l_nb_cells);
    }

    fprintf(stderr,"matrix %d x %d with %d non-zeros\n",l_mtx_meta.nb_rows, l_mtx_meta.nb_cols, l_nb_cells_in_mtx_file);

    // Parse mtx file
    for (l_mtx_element_counter=0; l_mtx_element_counter < l_nb_cells_in_mtx_file; l_mtx_element_counter++)
    {
        if ( l_mtx_meta.complex ) { // Processing file describing a complex matrix
            if (l_mtx_meta.nb_cols == 1 ) { // Processing file describing a vector
                l_nb_fields = fscanf(l_mtx_file, "%lg %lg\n", &l_value_pool[l_value_index], &l_value_pool[l_value_index+1]) ;
                l_row_index_1_based=l_mtx_element_counter+1;
                l_col_index_1_based=1;
            } else {
                l_nb_fields = fscanf(l_mtx_file, "%d %d %lg %lg\n", &l_row_index_1_based, &l_col_index_1_based, &l_value_pool[l_value_index], &l_value_pool[l_value_index+1]) ;
            }
        } else { // Processing file describing a real matrix
            if (l_mtx_meta.nb_cols == 1 ) { // Processing file describing a vector
                l_nb_fields = fscanf(l_mtx_file, "%lg\n", &l_value_pool[l_value_index]);
                l_row_index_1_based=l_mtx_element_counter+1;
                l_col_index_1_based=1;
            } else { 
                l_nb_fields = fscanf(l_mtx_file, "%d %d %lg\n", &l_row_index_1_based, &l_col_index_1_based, &l_value_pool[l_value_index]) ;
            }
        }

        if ( l_nb_fields == 0) {
            printf("Fail parsing MTX file line.\n");
            break;
        }

        if ( l_current_row_meta == NULL || l_current_row_meta->row_index != l_row_index_1_based ) {
            l_current_row_meta = mtx_meta_get_row_index(&l_mtx_meta, l_row_index_1_based);
        }

        l_row_added = row_meta_add_cell(&l_mtx_meta, l_current_row_meta, &l_cell_meta_pool[l_cell_meta_index++], l_col_index_1_based, &l_value_pool[l_value_index]);
        // l_mtx_element_counter++;

        // If matrix is symmetric register sysmmetric elements
        if ( l_row_added && mm_is_symmetric(l_matcode) && ( l_row_index_1_based != l_col_index_1_based ) ) {
            // je rajoute le symetrique a la fin
            if ( l_current_symmetric_row_meta == NULL || l_current_symmetric_row_meta->row_index != l_col_index_1_based ) {
                l_current_symmetric_row_meta = mtx_meta_get_row_index(&l_mtx_meta, l_col_index_1_based);
            }    

            // Parsing complex element line
            l_row_added = row_meta_add_cell(&l_mtx_meta, l_current_symmetric_row_meta, &l_cell_meta_pool[l_cell_meta_index++], l_row_index_1_based, &l_value_pool[l_value_index]);
            // l_mtx_element_counter++; fix
        }

        if ( l_mtx_meta.complex ) {
            // Parsing complex element line
            l_value_index += 2;
        } else {
            // Parsing real element line
            l_value_index += 1;
        }
    }

    fclose(l_mtx_file);

    l_crs = build_crs_matrix_from_meta(&l_mtx_meta);

    free(l_cell_meta_pool);
    free(l_value_pool);
   
    mtx_meta_clear(&l_mtx_meta);

    return l_crs;
}