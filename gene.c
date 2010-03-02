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

#include "gene.h"

/* Copies genes from the dynamic programming to a final array */

int add_genes(struct _gene *glist, struct _node *nod, int dbeg) {
  int path, pos, ctr;

  if(dbeg == -1) return 0;
  path = dbeg; ctr = 0; pos = -1;  
  while(nod[path].traceb != -1) path = nod[path].traceb;

  while(path != -1) {
    if(nod[path].elim == 1) { path = nod[path].tracef; continue; }
    if(nod[path].strand == 1 && nod[path].type != STOP) {
      glist[ctr].begin = nod[path].ndx+1;
      glist[ctr].start_ndx = path;
    }
    if(nod[path].strand == -1 && nod[path].type == STOP) {
      glist[ctr].begin = nod[path].ndx-1;
      glist[ctr].stop_ndx = path;
    }
    if(nod[path].strand == 1 && nod[path].type == STOP) {
      glist[ctr].end = nod[path].ndx+3;
      glist[ctr].stop_ndx = path;
      ctr++;
    }
    if(nod[path].strand == -1 && nod[path].type != STOP) {
      glist[ctr].end = nod[path].ndx+1;
      glist[ctr].start_ndx = path;
      ctr++;
    }
    path = nod[path].tracef;
    if(ctr == MAX_GENES) {
      fprintf(stderr, "warning, max # of genes exceeded, truncating...\n");
      return ctr;
    }
  }
  return ctr;
}

/*******************************************************************************
  This routine attempts to solve the problem of extremely close starts.  If two
  potential starts are 5 amino acids or less away from each other, this routine
  sets their coding equal to each other and lets the RBS/operon/ATG-GTG-TTG
  start features determine which start to use, under the assumption that 5 or
  less words of coding is too weak a signal to use to select the proper start.

  In addition, we try to correct TTG (or whatever start codon is rare) starts
  that have an RBS score, an upstream score, and a coding score all superior
  to whatever start we initially chose.

  This routine was tested on numerous genomes and found to increase overall
  performance. 
*******************************************************************************/
void tweak_final_starts(struct _gene *genes, int ng, struct _node *nod,
                        int nn, struct _training *tinf) {
  int i, j, ndx, mndx, maxndx[2];
  double sc, igm, tigm, maxsc[2], maxigm[2];

  for(i = 0; i < ng; i++) {
    ndx = genes[i].start_ndx;
    sc = nod[ndx].sscore + nod[ndx].cscore;
    igm = 0.0;
    if(i > 0 && nod[ndx].strand == 1 && nod[genes[i-1].start_ndx].strand == 1)
      igm = intergenic_mod(&nod[genes[i-1].stop_ndx], &nod[ndx], tinf);
    if(i > 0 && nod[ndx].strand == 1 && nod[genes[i-1].start_ndx].strand == -1)
      igm = intergenic_mod(&nod[genes[i-1].start_ndx], &nod[ndx], tinf);
    if(i < ng-1 && nod[ndx].strand == -1 && nod[genes[i+1].start_ndx].strand 
       == 1)
      igm = intergenic_mod(&nod[ndx], &nod[genes[i+1].start_ndx], tinf);
    if(i < ng-1 && nod[ndx].strand == -1 && nod[genes[i+1].start_ndx].strand 
       == -1)
      igm = intergenic_mod(&nod[ndx], &nod[genes[i+1].stop_ndx], tinf);

    /* Search upstream and downstream for the #2 and #3 scoring starts */
    maxndx[0] = -1; maxndx[1] = -1; maxsc[0] = 0; maxsc[1] = 0;
    maxigm[0] = 0; maxigm[1] = 0;
    for(j = ndx-100; j < ndx+100; j++) {
      if(j < 0 || j >= nn || j == ndx) continue;
      if(nod[j].type == STOP || nod[j].stop_val != nod[ndx].stop_val)
        continue;

      tigm = 0.0;
      if(i > 0 && nod[j].strand == 1 && nod[genes[i-1].start_ndx].strand == 1)
      {
        if(nod[genes[i-1].stop_ndx].ndx - nod[j].ndx > MAX_SAM_OVLP) continue;
        tigm = intergenic_mod(&nod[genes[i-1].stop_ndx], &nod[j], tinf);
      }
      if(i > 0 && nod[j].strand == 1 && nod[genes[i-1].start_ndx].strand == -1)
      {
        if(nod[genes[i-1].start_ndx].ndx - nod[j].ndx >= 0) continue;
        tigm = intergenic_mod(&nod[genes[i-1].start_ndx], &nod[j], tinf);
      }
      if(i < ng-1 && nod[j].strand == -1 && nod[genes[i+1].start_ndx].strand 
         == 1) {
        if(nod[j].ndx - nod[genes[i+1].start_ndx].ndx >= 0) continue;
        tigm = intergenic_mod(&nod[j], &nod[genes[i+1].start_ndx], tinf);
      }
      if(i < ng-1 && nod[j].strand == -1 && nod[genes[i+1].start_ndx].strand 
         == -1) {
        if(nod[j].ndx - nod[genes[i+1].stop_ndx].ndx > MAX_SAM_OVLP) continue;
        tigm = intergenic_mod(&nod[j], &nod[genes[i+1].stop_ndx], tinf);
      }
 
      if(maxndx[0] == -1) {
        maxndx[0] = j;
        maxsc[0] = nod[j].cscore + nod[j].sscore;
        maxigm[0] = tigm;
      }
      else if(nod[j].cscore + nod[j].sscore + tigm > maxsc[0]) {
        maxndx[1] = maxndx[0];
        maxsc[1] = maxsc[0];
        maxigm[1] = maxigm[0];
        maxndx[0] = j;
        maxsc[0] = nod[j].cscore + nod[j].sscore;
        maxigm[0] = tigm;
      }
      else if(maxndx[1] == -1 || nod[j].cscore + nod[j].sscore + tigm > 
              maxsc[1]) { 
        maxndx[1] = j;
        maxsc[1] = nod[j].cscore + nod[j].sscore;
        maxigm[1] = tigm;
      }
    }

    /* Change the start if it's a TTG with better coding/RBS/upstream score */
    /* Also change the start if it's <=15bp but has better coding/RBS       */
    for(j = 0; j < 2; j++) {
      mndx = maxndx[j];
      if(mndx == -1) continue;

      /* Start of less common type but with better coding, rbs, and */
      /* upstream.  Must be 18 or more bases away from original.    */
      if(nod[mndx].tscore < nod[ndx].tscore && maxsc[j]-nod[mndx].tscore >= 
         sc-nod[ndx].tscore+tinf->st_wt && nod[mndx].rscore > nod[ndx].rscore
         && nod[mndx].uscore > nod[ndx].uscore && nod[mndx].cscore > 
         nod[ndx].cscore && abs(nod[mndx].ndx-nod[ndx].ndx) > 15) {
        maxsc[j] += nod[ndx].tscore-nod[mndx].tscore;
      }

      /* Close starts.  Ignore coding and see if start has better rbs */
      /* and type. */
      else if(abs(nod[mndx].ndx-nod[ndx].ndx) <= 15 && nod[mndx].rscore+
              nod[mndx].tscore > nod[ndx].rscore+nod[ndx].tscore) {
        if(nod[ndx].cscore > nod[mndx].cscore) maxsc[j] += nod[ndx].cscore -
                                               nod[mndx].cscore;
        if(nod[ndx].uscore > nod[mndx].uscore) maxsc[j] += nod[ndx].uscore - 
                                               nod[mndx].uscore;
        if(igm > maxigm[j]) maxsc[j] += igm - maxigm[j]; 
      }
  
      else maxsc[j] = -1000.0;
    }

    /* Change the gene coordinates to the new maximum. */
    mndx = -1;
    for(j = 0; j < 2; j++) {
      if(maxndx[j] == -1) continue;
      if(mndx == -1 && maxsc[j]+maxigm[j] > sc+igm)
        mndx = j;
      else if(mndx >= 0 && maxsc[j]+maxigm[j] > maxsc[mndx]+maxigm[mndx])
        mndx = j; 
    }
    if(mndx != -1 && nod[maxndx[mndx]].strand == 1) {
      genes[i].start_ndx = maxndx[mndx];
      genes[i].begin = nod[maxndx[mndx]].ndx+1;
    } 
    else if(mndx != -1 && nod[maxndx[mndx]].strand == -1) {
      genes[i].start_ndx = maxndx[mndx];
      genes[i].end = nod[maxndx[mndx]].ndx+1;
    } 
  }
}

/* Print the genes.  'Flag' indicates which format to use. */
void print_genes(FILE *fp, struct _gene *genes, int ng, struct _node *nod, 
                 int slen, int flag, int sctr, int is_meta, char *mdesc) {
  int i;
  char left[50], right[50];

  strcpy(left, "");
  strcpy(right, "");

  if(flag == 0) fprintf(fp, "DEFINITION           Prodigal_Seq_%d (%d bp)\n", 
                        sctr, slen);
  else if(flag == 1) fprintf(fp, "sequence_prodigal=%d|%d\n", sctr, slen);
  else if(flag == 2) fprintf(fp, "# Prodigal_Seq_%d (%d bp)\n", sctr, slen);

  if(is_meta == 1 && flag == 0) fprintf(fp, "SOURCE               %s", mdesc);
  else if(is_meta == 1 && flag == 1) {
    for(i = 0; i < strlen(mdesc); i++) if(mdesc[i] == '\t') mdesc[i] = '|';
    fprintf(fp, "classification_prodigal=%s", mdesc);
  }
  else if(flag == 2) fprintf(fp, "# Classification: %s", mdesc);
  else if(flag == 3) fprintf(fp, "##classification: %s", mdesc);

  for(i = 0; i < ng; i++) {
    if(nod[genes[i].start_ndx].strand == 1) {
      if(nod[genes[i].start_ndx].edge == 1) sprintf(left, "<%d", 
         genes[i].begin);
      else sprintf(left, "%d", genes[i].begin);
      if(nod[genes[i].stop_ndx].edge == 1) sprintf(right, ">%d", 
         genes[i].end);
      else sprintf(right, "%d", genes[i].end);
      if(flag == 0) fprintf(fp, "     CDS             %s..%s\n", left, right);
      if(flag == 1)
        fprintf(fp, "gene_prodigal=%d|1|f|y|y|3|0|%d|%d|%d|%d|-1|-1|1.0\n", i+1,
                genes[i].begin, genes[i].end, genes[i].begin, genes[i].end);
      if(flag == 2) fprintf(fp, ">%d_%d_%d_+\n", i+1, genes[i].begin, 
                            genes[i].end);
      if(flag == 3)
        fprintf(fp, "Prod_Seq_%d\tProdigal\tCDS\t%d\t%d\t.\t+\t0\n", sctr, 
                genes[i].begin, genes[i].end);
    }
    else {
      if(nod[genes[i].stop_ndx].edge == 1) sprintf(left, "<%d", 
         genes[i].begin);
      else sprintf(left, "%d", genes[i].begin);
      if(nod[genes[i].start_ndx].edge == 1) sprintf(right, ">%d", 
         genes[i].end);
      else sprintf(right, "%d", genes[i].end);
      if(flag == 0) fprintf(fp, "     CDS             complement(%s..%s)\n", 
                            left, right);
      if(flag == 1)
        fprintf(fp, "gene_prodigal=%d|1|r|y|y|3|0|%d|%d|%d|%d|-1|-1|1.0\n", i+1,
               slen+1-genes[i].end, slen+1-genes[i].begin,
               slen+1-genes[i].end, slen+1-genes[i].begin);
      if(flag == 2) fprintf(fp, ">%d_%d_%d_-\n", i+1, genes[i].begin, 
                            genes[i].end);
      if(flag == 3)
        fprintf(fp, "Prod_Seq_%d\tProdigal\tCDS\t%d\t%d\t.\t-\t0\n", sctr, 
                genes[i].begin, genes[i].end);
    }
  }

  if(flag == 0) fprintf(fp, "//\n");
}

/* Print the gene translations */
void write_translations(FILE *fh, struct _gene *genes, int ng, struct 
                        _node *nod, unsigned char *seq, unsigned char *rseq, 
                        int slen, struct _training *tinf, int sctr) {
  int i, j;

  fprintf(fh, "# Prodigal Sequence %d\n", sctr);
  for(i = 0; i < ng; i++) {
    if(nod[genes[i].start_ndx].strand == 1) {
      fprintf(fh, ">Prodigal Gene %d # %d # %d # 1\n", i+1, genes[i].begin,
              genes[i].end);
      for(j = genes[i].begin; j < genes[i].end-3; j+=3) {
        fprintf(fh, "%c", amino(seq, j-1, tinf, j==genes[i].begin?1:0 &&
                (1-nod[genes[i].start_ndx].edge)));
        if((j-genes[i].begin)%180 == 177) fprintf(fh, "\n");
      }
      if((j-genes[i].begin)%180 != 0) fprintf(fh, "\n");
    }
    else {
      fprintf(fh, ">Prodigal Gene %d # %d # %d # -1\n", i+1, genes[i].begin,
              genes[i].end);
      for(j = slen+1-genes[i].end; j < slen+1-genes[i].begin-3; j+=3) {
        fprintf(fh, "%c", amino(rseq, j-1, tinf, j==slen+1-genes[i].end?1:0 &&
                (1-nod[genes[i].start_ndx].edge)));
        if((j-slen-1+genes[i].end)%180 == 177) fprintf(fh, "\n");
      }
      if((j-slen-1+genes[i].end)%180 != 0) fprintf(fh, "\n");
    }
  }
}

/* Print the gene nucleotide sequences */
void write_nucleotide_seqs(FILE *fh, struct _gene *genes, int ng, struct 
                           _node *nod, unsigned char *seq, unsigned char *rseq,
                           unsigned char *useq, int slen, struct _training 
                           *tinf, int sctr) {
  int i, j;

  fprintf(fh, "# Prodigal Sequence %d\n", sctr);
  for(i = 0; i < ng; i++) {
    if(nod[genes[i].start_ndx].strand == 1) {
      fprintf(fh, ">Prodigal Gene %d # %d # %d # 1\n", i+1, genes[i].begin,
              genes[i].end);
      for(j = genes[i].begin-1; j < genes[i].end; j++) {
        if(is_a(seq, j) == 1) fprintf(fh, "A");
        else if(is_t(seq, j) == 1) fprintf(fh, "T");
        else if(is_g(seq, j) == 1) fprintf(fh, "G");
        else if(is_c(seq, j) == 1 && is_n(useq, j) == 0) fprintf(fh, "C");
        else fprintf(fh, "N");
        if((j-genes[i].begin+1)%70 == 69) fprintf(fh, "\n");
      }
      if((j-genes[i].begin+1)%70 != 0) fprintf(fh, "\n");
    }
    else {
      fprintf(fh, ">Prodigal Gene %d # %d # %d # -1\n", i+1, genes[i].begin,
              genes[i].end);
      for(j = slen-genes[i].end; j < slen+1-genes[i].begin; j++) {
        if(is_a(rseq, j) == 1) fprintf(fh, "A");
        else if(is_t(rseq, j) == 1) fprintf(fh, "T");
        else if(is_g(rseq, j) == 1) fprintf(fh, "G");
        else if(is_c(rseq, j) == 1 && is_n(useq, slen-1-j) == 0) 
          fprintf(fh, "C");
        else fprintf(fh, "N");
        if((j-slen+genes[i].end)%70 == 69) fprintf(fh, "\n");
      }
      if((j-slen+genes[i].end)%70 != 0) fprintf(fh, "\n");
    }
  }
}
