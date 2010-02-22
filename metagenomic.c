/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2010 University of Tennessee / UT-Battelle

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

#include "metagenomic.h"


/*******************************************************************************
  Initialize the metagenomic bins with the 30 precalculated training files
  from the model organisms that best represent all of microbial Genbank.  
*******************************************************************************/
void initialize_metagenomic_bins(struct _metagenomic_bin *meta) {
  initialize_metagenome_0(meta[0].tinf);
  initialize_metagenome_1(meta[1].tinf);
  initialize_metagenome_2(meta[2].tinf);
  initialize_metagenome_3(meta[3].tinf);
  initialize_metagenome_4(meta[4].tinf);
  initialize_metagenome_5(meta[5].tinf);
  initialize_metagenome_6(meta[6].tinf);
  initialize_metagenome_7(meta[7].tinf);
  initialize_metagenome_8(meta[8].tinf);
  initialize_metagenome_9(meta[9].tinf);
  initialize_metagenome_10(meta[10].tinf);
  initialize_metagenome_11(meta[11].tinf);
  initialize_metagenome_12(meta[12].tinf);
  initialize_metagenome_13(meta[13].tinf);
  initialize_metagenome_14(meta[14].tinf);
  initialize_metagenome_15(meta[15].tinf);
  initialize_metagenome_16(meta[16].tinf);
  initialize_metagenome_17(meta[17].tinf);
  initialize_metagenome_18(meta[18].tinf);
  initialize_metagenome_19(meta[19].tinf);
  initialize_metagenome_20(meta[20].tinf);
  initialize_metagenome_21(meta[21].tinf);
  initialize_metagenome_22(meta[22].tinf);
  initialize_metagenome_23(meta[23].tinf);
  initialize_metagenome_24(meta[24].tinf);
  initialize_metagenome_25(meta[25].tinf);
  initialize_metagenome_26(meta[26].tinf);
  initialize_metagenome_27(meta[27].tinf);
  initialize_metagenome_28(meta[28].tinf);
  initialize_metagenome_29(meta[29].tinf);
}

/*******************************************************************************
  Calculate the coding score of a sequence interval in all 6 frames and        
  return the highest score.                                                    
*******************************************************************************/
double score_sample(unsigned char *seq, unsigned char *rseq, int slen, int
                    begin, int end, struct _training *tinf) {
  int i, j;
  double max = 0.0, cur = 0.0;
 
  for(i = 0; i < 3; i++) {
    cur = 0.0;
    for(j = begin+i; j <= end-5 && j <= slen-5; j+=3)
      cur += tinf->gene_dc[mer_ndx(6, seq, j)];
    if(cur > max) max = cur; 
    cur = 0.0;
    for(j = begin+i; j <= end-5 && j <= slen-5; j+=3)
      cur += tinf->gene_dc[mer_ndx(6, rseq, slen-j-6)];
    if(cur > max) max = cur; 
  } 
  return max;
}

/*******************************************************************************
  Given the 30 metagenomic training files and a sequence, randomly sample
  sequence fragments and score them using each of the files.  Sort the bins
  from best fit to worst fit.
*******************************************************************************/
void determine_top_bins(unsigned char *seq, unsigned char *rseq, int slen,
                        double gc, struct _metagenomic_bin *meta) {
  int i, j;
  double nsamp = 0.0, rnd = 0.0;

  srand(time(NULL));

  for(i = 0; i < 30; i++) {
    meta[i].index = i;
    meta[i].weight = 0.0;
    meta[i].gc = fabs(gc - meta[i].tinf->gc);
  }
  if(slen < 2*SAMPLE_LEN) {
    for(i = 0; i < 30; i++) 
      meta[i].weight = dmax(0.0, score_sample(seq, rseq, slen, 0, slen-1, 
                            meta[i].tinf)); 
  }
  else {
    nsamp = ((double)slen*5.0/SAMPLE_LEN);
    if(nsamp > MAX_SAMPLE) nsamp = MAX_SAMPLE;
    for(i = 0; i < nsamp; i++) {
      rnd = (int)(((double)rand())/((double)RAND_MAX) * (slen-SAMPLE_LEN-1));
      for(j = 0; j < 30; j++) {
        meta[j].weight += dmax(0.0, score_sample(seq, rseq, slen, i*SAMPLE_LEN,
                               (i+1)*SAMPLE_LEN-1, meta[j].tinf));
      }
    }
  }

  qsort(meta, 30, sizeof(struct _metagenomic_bin), &compare_meta_bins); 
}

/* Sorting routine for metagenomic bins */

int compare_meta_bins(const void *v1, const void *v2) {
  struct _metagenomic_bin *n1, *n2;
  n1 = (struct _metagenomic_bin *)v1;
  n2 = (struct _metagenomic_bin *)v2;
  if(n1->weight < n2->weight) return 1;
  if(n1->weight > n2->weight) return -1;
  if(n1->gc < n2->gc) return -1;
  if(n1->gc > n2->gc) return 1;
  return 0;
}

