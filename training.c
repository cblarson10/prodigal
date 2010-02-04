/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2009 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#include "training.h"

/* Reads a training file to use for gene prediction */

int read_training_file(char *fn, struct _training *tinf) {
  size_t rv;
  FILE *fh;
  fh = fopen(fn, "rb");
  if(fh == NULL) return -1;
  rv = fread(tinf, sizeof(struct _training), 1, fh);
  fclose(fh);
  if(rv != 1) return -1;
  return 0;
}

/* Writes a training file to use for a later run of gene prediction */

int write_training_file(char *fn, struct _training *tinf) {
  size_t rv;
  FILE *fh;
  fh = fopen(fn, "wb");
  if(fh == NULL) return -1;
  rv = fwrite(tinf, sizeof(struct _training), 1, fh);
  fclose(fh);
  if(rv != 1) return -1;
  return 0;
}