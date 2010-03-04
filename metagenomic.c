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
  sprintf(meta[0].desc, "%d|%s|%s|%s|%.1f|%d|%d", 0, 
          "Micrococcus luteus NCTC 2665", "B", "Actinobacteria", 72.9,
          meta[0].tinf->trans_table, meta[0].tinf->uses_sd);
  sprintf(meta[1].desc, "%d|%s|%s|%s|%.1f|%d|%d", 1, 
          "Thiobacillus denitrificans ATCC 25259", "B", "Betaproteobacteria",
          66.1, meta[1].tinf->trans_table, meta[1].tinf->uses_sd);
  sprintf(meta[2].desc, "%d|%s|%s|%s|%.1f|%d|%d", 2, 
          "Bradyrhizobium sp. ORS278", "B", "Alphaproteobacteria", 65.5,
          meta[2].tinf->trans_table, meta[2].tinf->uses_sd);
  sprintf(meta[3].desc, "%d|%s|%s|%s|%.1f|%d|%d", 3, 
          "Mycobacterium abscessus", "B", "Actinobacteria", 64.1,
          meta[3].tinf->trans_table, meta[3].tinf->uses_sd);
  sprintf(meta[4].desc, "%d|%s|%s|%s|%.1f|%d|%d", 4, 
          "Halorhabdus utahensis DSM 12940", "A", "Euryarchaeota", 62.9,
          meta[4].tinf->trans_table, meta[4].tinf->uses_sd);
  sprintf(meta[5].desc, "%d|%s|%s|%s|%.1f|%d|%d", 5, 
          "Hyphomonas neptunium ATCC 15444", "B", "Alphaproteobacteria", 61.9,
          meta[5].tinf->trans_table, meta[5].tinf->uses_sd);
  sprintf(meta[6].desc, "%d|%s|%s|%s|%.1f|%d|%d", 6, 
          "Pseudomonas putida W619", "B", "Gammaproteobacteria", 61.4,
          meta[6].tinf->trans_table, meta[6].tinf->uses_sd);
  sprintf(meta[7].desc, "%d|%s|%s|%s|%.1f|%d|%d", 7, 
          "Synechococcus sp. WH 7803", "B", "Cyanobacteria", 60.2,
          meta[7].tinf->trans_table, meta[7].tinf->uses_sd);
  sprintf(meta[8].desc, "%d|%s|%s|%s|%.1f|%d|%d", 8, 
          "Mycobacterium leprae Br4923", "B", "Actinobacteria", 57.8,
          meta[8].tinf->trans_table, meta[8].tinf->uses_sd);
  sprintf(meta[9].desc, "%d|%s|%s|%s|%.1f|%d|%d", 9, 
          "Pelobacter carbinolicus DSM 2380", "B", "Deltaproteobacteria", 55.1,
          meta[9].tinf->trans_table, meta[9].tinf->uses_sd);
  sprintf(meta[10].desc, "%d|%s|%s|%s|%.1f|%d|%d", 10, 
          "Pyrobaculum arsenaticum DSM 13514", "A", "Crenarchaeota", 55.1,
          meta[10].tinf->trans_table, meta[10].tinf->uses_sd);
  sprintf(meta[11].desc, "%d|%s|%s|%s|%.1f|%d|%d", 11, 
          "Corynebacterium glutamicum ATCC 13032", "B", "Actinobacteria", 53.8,
          meta[11].tinf->trans_table, meta[11].tinf->uses_sd);
  sprintf(meta[12].desc, "%d|%s|%s|%s|%.1f|%d|%d", 12, 
          "Synechococcus sp. CC9311", "B", "Cyanobacteria", 52.4,
          meta[12].tinf->trans_table, meta[12].tinf->uses_sd);
  sprintf(meta[13].desc, "%d|%s|%s|%s|%.1f|%d|%d", 13, 
          "Pectobacterium atrosepticum SCRI1043", "B","Gammaproteobacteria",
          51.0, meta[13].tinf->trans_table, meta[13].tinf->uses_sd);
  sprintf(meta[14].desc, "%d|%s|%s|%s|%.1f|%d|%d", 14, 
          "Prosthecochloris aestuarii DSM 271", "B", "Bacteroidetes/Chlorobi",
          50.1, meta[14].tinf->trans_table, meta[14].tinf->uses_sd);
  sprintf(meta[15].desc, "%d|%s|%s|%s|%.1f|%d|%d", 15, 
          "Anaplasma marginale str. Florida", "B", "Alphaproteobacteria", 49.8,
          meta[15].tinf->trans_table, meta[15].tinf->uses_sd);
  sprintf(meta[16].desc, "%d|%s|%s|%s|%.1f|%d|%d", 16, 
          "Nitrosomonas eutropha C91", "B", "Betaproteobacteria", 48.5,
          meta[16].tinf->trans_table, meta[16].tinf->uses_sd);
  sprintf(meta[17].desc, "%d|%s|%s|%s|%.1f|%d|%d", 17, 
          "Tropheryma whipplei TW08/27", "B", "Actinobacteria", 46.3,
          meta[17].tinf->trans_table, meta[17].tinf->uses_sd);
  sprintf(meta[18].desc, "%d|%s|%s|%s|%.1f|%d|%d", 18, 
          "Pyrococcus abyssi GE5", "A", "Euryarchaeota", 44.7,
          meta[18].tinf->trans_table, meta[18].tinf->uses_sd);
  sprintf(meta[19].desc, "%d|%s|%s|%s|%.1f|%d|%d", 19, 
          "Sulfurovum sp. NBC37-1", "B", "Epsilonproteobacteria", 43.9,
          meta[19].tinf->trans_table, meta[19].tinf->uses_sd);
  sprintf(meta[20].desc, "%d|%s|%s|%s|%.1f|%d|%d", 20, 
          "Bacteroides fragilis YCH46", "B", "Bacteroidetes/Chlorobi", 43.2,
          meta[0].tinf->trans_table, meta[20].tinf->uses_sd);
  sprintf(meta[21].desc, "%d|%s|%s|%s|%.1f|%d|%d", 21, 
          "Marinomonas sp. MWYL1", "B", "Gammaproteobacteria", 42.6,
          meta[21].tinf->trans_table, meta[21].tinf->uses_sd);
  sprintf(meta[22].desc, "%d|%s|%s|%s|%.1f|%d|%d", 22, 
          "Legionella pneumophila str. Corby", "B", "Gammaproteobacteria", 
          38.5, meta[22].tinf->trans_table, meta[22].tinf->uses_sd);
  sprintf(meta[23].desc, "%d|%s|%s|%s|%.1f|%d|%d", 23, 
          "Oceanobacillus iheyensis HTE831", "B", "Firmicutes", 35.7,
          meta[23].tinf->trans_table, meta[23].tinf->uses_sd);
  sprintf(meta[24].desc, "%d|%s|%s|%s|%.1f|%d|%d", 24, 
          "Wolbachia endosymbiont of Drosophila melanogaster", "B", 
          "Alphaproteobacteria", 35.2,
          meta[24].tinf->trans_table, meta[24].tinf->uses_sd);
  sprintf(meta[25].desc, "%d|%s|%s|%s|%.1f|%d|%d", 25, 
          "Sulfolobus islandicus L.S.2.15", "A", "Crenarchaeota", 35.1,
          meta[25].tinf->trans_table, meta[25].tinf->uses_sd);
  sprintf(meta[26].desc, "%d|%s|%s|%s|%.1f|%d|%d", 26, 
          "Francisella philomiragia subsp. philomiragia ATCC 25017", "B",
          "Gammaproteobacteria", 32.6,
          meta[26].tinf->trans_table, meta[26].tinf->uses_sd);
  sprintf(meta[27].desc, "%d|%s|%s|%s|%.1f|%d|%d", 27, 
          "Rickettsia felis URRWXCal2", "B", "Alphaproteobacteria", 32.5,
          meta[27].tinf->trans_table, meta[27].tinf->uses_sd);
  sprintf(meta[28].desc, "%d|%s|%s|%s|%.1f|%d|%d", 28, 
          "Clostridium acetobutylicum ATCC 824", "B", "Firmicutes", 30.9,
          meta[28].tinf->trans_table, meta[28].tinf->uses_sd);
  sprintf(meta[29].desc, "%d|%s|%s|%s|%.1f|%d|%d", 29, 
          "Mycoplasma arthritidis 158L3-1", "B", "Firmicutes", 30.7,
          meta[29].tinf->trans_table, meta[29].tinf->uses_sd);

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

