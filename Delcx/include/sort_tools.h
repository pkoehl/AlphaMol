/*********************************************************************************
 *      A set of tools for sorting
 *********************************************************************************/

#ifndef SORT_TOOLS_H
#define SORT_TOOLS_H

  #include <vector>
  #include "Vertex.h"

/*********************************************************************************
  Class 
 *********************************************************************************/

  class SortingTools {

	public:
		void valsort2(int a, int b, int *ia, int *ib, int *iswap);
		void valsort3(int a, int b, int c, int *ia, int *ib, int *ic, int *iswap);
		void valsort4(int a, int b, int c, int d, int *ia, int *ib, int *ic, int *id, int *iswap);
		void valsort5(int a, int b, int c, int d, int e, 
			int *ia, int *ib, int *ic, int *id, int *ie, int *iswap);
		void isort_indx(int *list, int *idx, int *nswap, int n);
		void isort_swap(int *list, int *nswap, int n);
		void isort4_swap(int *a, int *b, int *c, int *d, int *nswap);
		void sort4_sign(int *list, int *idx, int *nswap, int n);

  };

/*********************************************************************************
 valsort2: sort two integers and count number of swaps
 *********************************************************************************/

  void SortingTools::valsort2(int a, int b, int *ia, int *ib, int *iswap)
  {
	*iswap = 1;
	if( a > b) {
		*ia = b;
		*ib = a;
		*iswap = -*iswap;
	} else {
		*ia = a;
		*ib = b;
	}
  }

/*********************************************************************************
 valsort3: sort three integers and count number of swaps
 *********************************************************************************/

  void SortingTools::valsort3(int a, int b, int c, int *ia, int *ib, int *ic, int *iswap)
  {
	valsort2(a, b, ia, ib, iswap);

	*ic = c;

	int temp;
	if(*ib > *ic) {
		temp = *ib;
		*ib = *ic;
		*ic = temp;
		*iswap = -*iswap;
		if(*ia > *ib) {
			temp = *ia;
			*ia = *ib;
			*ib = temp;
			*iswap = -*iswap;
		}
	}
  }

/*********************************************************************************
 valsort4: sort four integers and count number of swaps
 *********************************************************************************/

  void SortingTools::valsort4(int a, int b, int c, int d, int *ia, int *ib, int *ic, int *id, int *iswap)
  {
	valsort3(a, b, c, ia, ib, ic, iswap);

	*id = d;

	int temp;
	if(*ic > *id) {
		temp = *ic;
		*ic = *id;
		*id = temp;
		*iswap = -*iswap;
		if(*ib > *ic) {
			temp = *ib;
			*ib = *ic;
			*ic = temp;
			*iswap = -*iswap;
			if(*ia > *ib) {
				temp = *ia;
				*ia = *ib;
				*ib = temp;
				*iswap = -*iswap;
			}
		}
	}
   }


/*********************************************************************************
 valsort5: sort five integers and count number of swaps
 *********************************************************************************/

  void SortingTools::valsort5(int a, int b, int c, int d, int e, 
	int *ia, int *ib, int *ic, int *id, int *ie, int *iswap)
  {
	valsort4(a, b, c, d, ia, ib, ic, id, iswap);

	*ie = e;
	int temp;

	if(*id > *ie) {
		temp = *id;
		*id = *ie;
		*ie = temp;
		*iswap = -*iswap;
		if(*ic > *id) {
			temp = *ic;
			*ic = *id;
			*id = temp;
			*iswap = -*iswap;
			if(*ib > *ic) {
				temp = *ib;
				*ib = *ic;
				*ic = temp;
				*iswap = -*iswap;
				if(*ia > *ib) {
					temp = *ia;
					*ia = *ib;
					*ib = temp;
					*iswap = -*iswap;
				}
			}
		}
	}

  }

/*********************************************************************************
 isort_indx: sort a list of integers and define index array
 *********************************************************************************/

  void SortingTools::isort_indx(int *list, int *idx, int *nswap, int n)
  {
	for(int i = 0; i < n; i++) idx[i] = i;

	*nswap = 0;

	int a;

	for(int i = 0; i < n-1; i++) {
		for(int j = i+1; j < n; j++) {
			if(list[i] > list[j]) {
				a = list[i];
				list[i] = list[j];
				list[j] = a;
				a = idx[i];
				idx[i] = idx[j];
				idx[j] = a;
				(*nswap)++;
			}
		}
	}

  }

/*********************************************************************************
 isort_swap: sort a list of integers and count number of swaps
 *********************************************************************************/

  void SortingTools::isort_swap(int *list, int *nswap, int n)
  {
	*nswap = 0;

	int a;

	for(int i = 0; i < n-1; i++) {
		for(int j = i+1; j < n; j++) {
			if(list[i] > list[j]) {
				a = list[i];
				list[i] = list[j];
				list[j] = a;
				(*nswap)++;
			}
		}
	}
  }

/*********************************************************************************
 isort4_swap: sort a list of integers and count number of swaps
 *********************************************************************************/

  void SortingTools::isort4_swap(int *a, int *b, int *c, int *d, int *nswap)
  {

	int temp;
	if(*a > *b) {
		temp = *a;
		*a = *b;
		*b = temp;
		*nswap = 1;
	}

	if(*a > *c) {
		temp = *a;
		*a = *c;
		*c = temp;
		(*nswap)++;
	}

	if(*a > *d) {
		temp = *a;
		*a = *d;
		*d = temp;
		(*nswap)++;
	}

	if(*b > *c) {
		temp = *b;
		*b = *c;
		*c = temp;
		(*nswap)++;
	}

	if(*b > *d) {
		temp = *b;
		*b = *d;
		*d = temp;
		(*nswap)++;
	}

	if(*c > *d) {
		temp = *c;
		*c = *d;
		*d = temp;
		(*nswap)++;
	}
  }

/*********************************************************************************
 sort4_sign: 
 sorts the list of 4 numbers, and computes the signature of the permutation
 *********************************************************************************/

  void SortingTools::sort4_sign(int *list, int *idx, int *nswap, int n)
  {
	for(int i = 0; i < n; i++) idx[i] = i;

	*nswap = 1;

	int a;

	for(int i = 0; i < n-1; i++) {
		for(int j = i+1; j < n; j++) {
			if(list[i] > list[j]) {
				a = list[i];
				list[i] = list[j];
				list[j] = a;
				a = idx[i];
				idx[i] = idx[j];
				idx[j] = a;
				*nswap = -*nswap;
			}
		}
	}

  }

#endif
