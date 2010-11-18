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
  Initialize the metagenomic bins with the precalculated training files
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
  initialize_metagenome_30(meta[30].tinf);
  initialize_metagenome_31(meta[31].tinf);
  initialize_metagenome_32(meta[32].tinf);
  initialize_metagenome_33(meta[33].tinf);
  initialize_metagenome_34(meta[34].tinf);
  initialize_metagenome_35(meta[35].tinf);
  initialize_metagenome_36(meta[36].tinf);
  initialize_metagenome_37(meta[37].tinf);
  initialize_metagenome_38(meta[38].tinf);
  sprintf(meta[0].desc, "%d|%s|%s|%.1f|%d|%d", 0,
          "Mycoplasma mycoides subsp. mycoides SC str. Gladysdale MU clone",
          "B", 23.95, meta[0].tinf->trans_table, meta [0].tinf->uses_sd);
  sprintf(meta[1].desc, "%d|%s|%s|%.1f|%d|%d", 1,
          "Mycoplasma gallisepticum str. R(high), complete genome.",
          "B", 31.47, meta[1].tinf->trans_table, meta [1].tinf->uses_sd);
  sprintf(meta[2].desc, "%d|%s|%s|%.1f|%d|%d", 2,
          "Mycoplasma conjunctivae HRC/581 chromosome, complete genome.",
          "B", 28.49, meta[2].tinf->trans_table, meta [2].tinf->uses_sd);
  sprintf(meta[3].desc, "%d|%s|%s|%.1f|%d|%d", 3,
          "Mycoplasma pneumoniae FH, complete genome.",
          "B", 40.00, meta[3].tinf->trans_table, meta [3].tinf->uses_sd);
  sprintf(meta[4].desc, "%d|%s|%s|%.1f|%d|%d", 4,
          "Mycoplasma genitalium G37, complete genome.",
          "B", 31.69, meta[4].tinf->trans_table, meta [4].tinf->uses_sd);
  sprintf(meta[5].desc, "%d|%s|%s|%.1f|%d|%d", 5,
          "Erwinia pyrifoliae DSM 12163 complete genome, culture collection",
          "B", 53.41, meta[5].tinf->trans_table, meta [5].tinf->uses_sd);
  sprintf(meta[6].desc, "%d|%s|%s|%.1f|%d|%d", 6,
          "Coxiella burnetii CbuK_Q154 chromosome, complete genome.",
          "B", 42.67, meta[6].tinf->trans_table, meta [6].tinf->uses_sd);
  sprintf(meta[7].desc, "%d|%s|%s|%.1f|%d|%d", 7,
          "Petrotoga mobilis SJ95, complete genome.",
          "B", 34.12, meta[7].tinf->trans_table, meta [7].tinf->uses_sd);
  sprintf(meta[8].desc, "%d|%s|%s|%.1f|%d|%d", 8,
          "Streptococcus pyogenes MGAS2096 chromosome, complete genome.",
          "B", 38.73, meta[8].tinf->trans_table, meta [8].tinf->uses_sd);
  sprintf(meta[9].desc, "%d|%s|%s|%.1f|%d|%d", 9,
          "Nitrobacter winogradskyi Nb-255, complete genome.",
          "B", 62.05, meta[9].tinf->trans_table, meta [9].tinf->uses_sd);
  sprintf(meta[10].desc, "%d|%s|%s|%.1f|%d|%d", 10,
          "Burkholderia sp. 383 chromosome 2, complete sequence.",
          "B", 66.73, meta[10].tinf->trans_table, meta [10].tinf->uses_sd);
  sprintf(meta[11].desc, "%d|%s|%s|%.1f|%d|%d", 11,
          "Anaplasma marginale str. St. Maries, complete genome.",
          "B", 49.76, meta[11].tinf->trans_table, meta [11].tinf->uses_sd);
  sprintf(meta[12].desc, "%d|%s|%s|%.1f|%d|%d", 12,
          "Candidatus Methanosphaerula palustris E1-9c, complete genome.",
          "A", 55.35, meta[12].tinf->trans_table, meta [12].tinf->uses_sd);
  sprintf(meta[13].desc, "%d|%s|%s|%.1f|%d|%d", 13,
          "Yersinia pseudotuberculosis YPIII chromosome, complete genome.",
          "B", 47.53, meta[13].tinf->trans_table, meta [13].tinf->uses_sd);
  sprintf(meta[14].desc, "%d|%s|%s|%.1f|%d|%d", 14,
          "Rickettsia africae ESF-5 chromosome, complete genome.",
          "B", 32.40, meta[14].tinf->trans_table, meta [14].tinf->uses_sd);
  sprintf(meta[15].desc, "%d|%s|%s|%.1f|%d|%d", 15,
          "Candidatus Amoebophilus asiaticus 5a2 chromosome, complete genome.",
          "B", 35.05, meta[15].tinf->trans_table, meta [15].tinf->uses_sd);
  sprintf(meta[16].desc, "%d|%s|%s|%.1f|%d|%d", 16,
          "Methanothermobacter thermautotrophicus str. Delta H chromosome,",
          "A", 49.54, meta[16].tinf->trans_table, meta [16].tinf->uses_sd);
  sprintf(meta[17].desc, "%d|%s|%s|%.1f|%d|%d", 17,
          "Methanococcus maripaludis C7, complete genome.",
          "A", 33.28, meta[17].tinf->trans_table, meta [17].tinf->uses_sd);
  sprintf(meta[18].desc, "%d|%s|%s|%.1f|%d|%d", 18,
          "Sulfolobus islandicus Y.G.57.14 chromosome, complete genome.",
          "A", 35.39, meta[18].tinf->trans_table, meta [18].tinf->uses_sd);
  sprintf(meta[19].desc, "%d|%s|%s|%.1f|%d|%d", 19,
          "Streptomyces griseus subsp. griseus NBRC 13350, complete genome.",
          "B", 72.23, meta[19].tinf->trans_table, meta [19].tinf->uses_sd);
  sprintf(meta[20].desc, "%d|%s|%s|%.1f|%d|%d", 20,
          "Aeropyrum pernix K1, complete genome.",
          "A", 56.31, meta[20].tinf->trans_table, meta [20].tinf->uses_sd);
  sprintf(meta[21].desc, "%d|%s|%s|%.1f|%d|%d", 21,
          "Chlorobium phaeobacteroides BS1, complete genome.",
          "B", 48.93, meta[21].tinf->trans_table, meta [21].tinf->uses_sd);
  sprintf(meta[22].desc, "%d|%s|%s|%.1f|%d|%d", 22,
          "Clostridium botulinum F str. 230613, complete genome.",
          "B", 28.30, meta[22].tinf->trans_table, meta [22].tinf->uses_sd);
  sprintf(meta[23].desc, "%d|%s|%s|%.1f|%d|%d", 23,
          "Candidatus Desulforudis audaxviator MP104C, complete genome.",
          "B", 60.85, meta[23].tinf->trans_table, meta [23].tinf->uses_sd);
  sprintf(meta[24].desc, "%d|%s|%s|%.1f|%d|%d", 24,
          "Synechococcus sp. CC9311, complete genome.",
          "B", 52.45, meta[24].tinf->trans_table, meta [24].tinf->uses_sd);
  sprintf(meta[25].desc, "%d|%s|%s|%.1f|%d|%d", 25,
          "Synechococcus sp. WH 7803, complete genome.",
          "B", 60.24, meta[25].tinf->trans_table, meta [25].tinf->uses_sd);
  sprintf(meta[26].desc, "%d|%s|%s|%.1f|%d|%d", 26,
          "Halomicrobium mukohataei DSM 12286, complete genome.",
          "A", 65.63, meta[26].tinf->trans_table, meta [26].tinf->uses_sd);
  sprintf(meta[27].desc, "%d|%s|%s|%.1f|%d|%d", 27,
          "Natrialba magadii ATCC 43099, complete genome.",
          "A", 61.42, meta[27].tinf->trans_table, meta [27].tinf->uses_sd);
  sprintf(meta[28].desc, "%d|%s|%s|%.1f|%d|%d", 28,
          "Desulfurococcus kamchatkensis 1221n chromosome, complete genome.",
          "A", 45.34, meta[28].tinf->trans_table, meta [28].tinf->uses_sd);
  sprintf(meta[29].desc, "%d|%s|%s|%.1f|%d|%d", 29,
          "Methanosarcina acetivorans C2A chromosome, complete genome.",
          "A", 42.68, meta[29].tinf->trans_table, meta [29].tinf->uses_sd);
  sprintf(meta[30].desc, "%d|%s|%s|%.1f|%d|%d", 30,
          "Bacteroides fragilis NCTC 9343 chromosome, complete genome.",
          "B", 43.19, meta[30].tinf->trans_table, meta [30].tinf->uses_sd);
  sprintf(meta[31].desc, "%d|%s|%s|%.1f|%d|%d", 31,
          "Chlorobaculum parvum NCIB 8327, complete genome.",
          "B", 55.80, meta[31].tinf->trans_table, meta [31].tinf->uses_sd);
  sprintf(meta[32].desc, "%d|%s|%s|%.1f|%d|%d", 32,
          "Staphylothermus hellenicus DSM 12710, complete genome.",
          "A", 36.84, meta[32].tinf->trans_table, meta [32].tinf->uses_sd);
  sprintf(meta[33].desc, "%d|%s|%s|%.1f|%d|%d", 33,
          "Methanococcus voltae A3, complete genome.",
          "A", 28.59, meta[33].tinf->trans_table, meta [33].tinf->uses_sd);
  sprintf(meta[34].desc, "%d|%s|%s|%.1f|%d|%d", 34,
          "Archaeoglobus fulgidus DSM 4304, complete genome.",
          "A", 48.58, meta[34].tinf->trans_table, meta [34].tinf->uses_sd);
  sprintf(meta[35].desc, "%d|%s|%s|%.1f|%d|%d", 35,
          "Sulfolobus tokodaii str. 7 chromosome, complete genome.",
          "A", 32.79, meta[35].tinf->trans_table, meta [35].tinf->uses_sd);
  sprintf(meta[36].desc, "%d|%s|%s|%.1f|%d|%d", 36,
          "Pyrobaculum aerophilum str. IM2 chromosome, complete genome.",
          "A", 51.36, meta[36].tinf->trans_table, meta [36].tinf->uses_sd);
  sprintf(meta[37].desc, "%d|%s|%s|%.1f|%d|%d", 37,
          "Candidatus Phytoplasma mali, complete genome.",
          "B", 21.39, meta[37].tinf->trans_table, meta [37].tinf->uses_sd);
  sprintf(meta[38].desc, "%d|%s|%s|%.1f|%d|%d", 38,
          "Thermus thermophilus HB27, complete genome.",
          "B", 69.44, meta[38].tinf->trans_table, meta [38].tinf->uses_sd);
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
  Given the metagenomic training files and a sequence, randomly sample
  sequence fragments and score them using each of the files.  Sort the bins
  from best fit to worst fit.
*******************************************************************************/
void determine_top_bins(unsigned char *seq, unsigned char *rseq, int slen,
                        double gc, struct _metagenomic_bin *meta) {
  int i, j;
  double nsamp = 0.0, rnd = 0.0;

  srand(time(NULL));

  for(i = 0; i < NUM_META; i++) {
    meta[i].index = i;
    meta[i].weight = 0.0;
    meta[i].gc = fabs(gc - meta[i].tinf->gc);
  }

  nsamp = ((double)slen/SAMPLE_LEN);
  if(nsamp < MAX_SAMPLE) {
    for(i = 0; i < nsamp; i++) 
      for(j = 0; j < NUM_META; j++) {
        meta[j].weight += dmax(0.0, score_sample(seq, rseq, slen, i*SAMPLE_LEN,
                              (i+1)*SAMPLE_LEN-1, meta[j].tinf)); 
      }
  }
  else {
    for(i = 0; i < MAX_SAMPLE; i++) {
      rnd = (int)(((double)rand())/((double)RAND_MAX) * (slen-SAMPLE_LEN-1));
      for(j = 0; j < NUM_META; j++) {
        meta[j].weight += dmax(0.0, score_sample(seq, rseq, slen, rnd, 
                               rnd+SAMPLE_LEN-1, meta[j].tinf));
      }
    }
  }

  qsort(meta, NUM_META, sizeof(struct _metagenomic_bin), &compare_meta_bins);
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

