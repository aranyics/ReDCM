/* multifit_nlin/gsl_multifit_nlin_ext.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MULTIFIT_NLIN_EXT_H__
#define __GSL_MULTIFIT_NLIN_EXT_H__

#include <stdlib.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


int gsl_multifit_fdfsolver_dif_df_simple(const gsl_vector *x, gsl_multifit_function_fdf *fdf,
                                  const gsl_vector *f, gsl_matrix *J);
int gsl_multifit_fdfsolver_dif_fdf_simple(const gsl_vector *x, gsl_multifit_function_fdf *fdf,
                                   gsl_vector *f, gsl_matrix *J);

__END_DECLS

#endif /* __GSL_MULTIFIT_NLIN_EXT_H__ */
