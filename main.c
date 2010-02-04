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

#include "sequence.h"
#include "node.h"
#include "dprog.h"
#include "gene.h"

void version();
void usage(char *);
void help();

int main(int argc, char *argv[]) {

  int rv, slen, nn, ng, i, ipath, *gc_frame, do_training, output;
  int closed, do_mask, nmask, force_nonsd, user_tt;
  unsigned char *seq, *rseq;
  char *train_file, *start_file, *trans_file;
  struct _node *nodes;
  struct _gene *genes;
  struct _training tinf;
  mask mlist[MAX_MASKS];

  /* Allocate memory and initialize variables */
  seq = (unsigned char *)malloc(MAX_SEQ*sizeof(unsigned char));
  rseq = (unsigned char *)malloc(MAX_SEQ*sizeof(unsigned char));
  nodes = (struct _node *)malloc(MAX_NODES*sizeof(struct _node));
  genes = (struct _gene *)malloc(MAX_GENES*sizeof(struct _gene));
  if(seq == NULL || rseq == NULL || nodes == NULL || genes == NULL) {
    fprintf(stderr, "Malloc failed on sequence/orfs\n\n"); exit(1);
  }
  memset(seq, 0, MAX_SEQ*sizeof(unsigned char));
  memset(rseq, 0, MAX_SEQ*sizeof(unsigned char));
  memset(nodes, 0, MAX_NODES*sizeof(struct _node));
  memset(genes, 0, MAX_GENES*sizeof(struct _gene));
  memset(&tinf, 0, sizeof(struct _training));
  nn = 0; slen = 0; ipath = 0; ng = 0; user_tt = 0;
  train_file = NULL; do_training = 0;
  start_file = NULL; trans_file = NULL;
  output = 0; closed = 0; nmask = 0; do_mask = 0; force_nonsd = 0;

  /***************************************************************************
    Set the start score weight.  Changing this number can dramatically
    affect the performance of the program.  Some genomes want it high (6+),
    and some prefer it low (2.5-3).  Attempts were made to determine this 
    weight dynamically, but none were successful.  Therefore, we just 
    manually set the weight to an average value that seems to work decently 
    for 99% of genomes.  This problem may be revisited in future versions.
  ***************************************************************************/
  tinf.st_wt = 4.35;
  tinf.trans_table = 11;

  /* Parse the command line arguments */
  for(i = 1; i < argc; i++) {
    if(i == argc-1 && (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "-T") == 0
       || strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "-A") == 0 ||
       strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "-g") == 0 ||
       strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "-O") == 0 ||
       strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "-S") == 0))
      usage("-g, -s, -t, and -o options require parameters.");
    else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "-C") == 0)
      closed = 1;
    else if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "-M") == 0)
      do_mask = 1;
    else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "-N") == 0)
      force_nonsd = 1;
    else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-H") == 0) help();
    else if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "-V") == 0) version();
    else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "-T") == 0) {
      train_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "-A") == 0) {
      trans_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "-S") == 0) {
      start_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "-G") == 0) {
      tinf.trans_table = atoi(argv[i+1]);
      if(tinf.trans_table < 1 || tinf.trans_table > 23 || tinf.trans_table == 7
         || tinf.trans_table == 8 || (tinf.trans_table >= 17 && tinf.trans_table
         <= 20))
        usage("Invalid translation table specified.");
      user_tt = tinf.trans_table;
      i++;
    }
    else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "-O") == 0) {
      if(strncmp(argv[i+1], "0", 1) == 0 || strcmp(argv[i+1], "gbk") == 0 ||
         strcmp(argv[i+1], "GBK") == 0)
        output = 0;
      else if(strncmp(argv[i+1], "1", 1) == 0 || strcmp(argv[i+1], "gca") == 0
              || strcmp(argv[i+1], "GCA") == 0)
        output = 1;
      else if(strncmp(argv[i+1], "2", 1) == 0 || strcmp(argv[i+1], "sco") == 0
              || strcmp(argv[i+1], "SCO") == 0)
        output = 2;
      else if(strncmp(argv[i+1], "3", 1) == 0 || strcmp(argv[i+1], "gff") == 0
              || strcmp(argv[i+1], "GFF") == 0)
        output = 3;
      else usage("Invalid output format specified.");
      i++;
    }
    else usage("Unknown option.");
  }

  /* Print header */
  fprintf(stderr, "-------------------------------------\n");
  fprintf(stderr, "PRODIGAL v1.20 [September, 2009]     \n");
  fprintf(stderr, "Univ of Tenn / Oak Ridge National Lab\n");
  fprintf(stderr, "Doug Hyatt, Loren Hauser, et al.     \n");
  fprintf(stderr, "-------------------------------------\n");

  /* Read in the training file (if specified) */
  if(train_file != NULL) {
    rv = read_training_file(train_file, &tinf);
    if(rv == -1) do_training = 1;
    else {
      fprintf(stderr, "Reading in training data from file %s...", train_file);
      if(force_nonsd == 1) { 
        fprintf(stderr, "error, cannot force non-SD finder with a training");
        fprintf(stderr, " file already created!\n"); exit(1);
      }
      if(user_tt > 0 && user_tt != tinf.trans_table) { 
        fprintf(stderr, "warning, user-specified translation table does not ");
        fprintf(stderr, "match the one in the specified training file! \n");
      }
      fflush(stderr);
      if(rv != 0) { fprintf(stderr, "error!\n"); exit(1); }
      fprintf(stderr, "done!\n"); fflush(stderr);
    }
  }

  /***************************************************************************
    Read in the sequence and create a reverse complement sequence.  Also 
    gather a dicodon background for the sequence.
  ***************************************************************************/
  fprintf(stderr, "Reading in the sequence..."); fflush(stderr);
  slen = read_seq(stdin, seq, &(tinf.gc), do_mask, mlist, &nmask, do_training);
  if(slen == 0) {
    fprintf(stderr, "Sequence read failed (file must be Fasta, Genbank, or");
    fprintf(stderr, "  EMBL format).\n\n");
    exit(1);
  }
  if((do_training == 1 || (do_training == 0 && train_file == NULL)) && slen <
      1000) {
    fprintf(stderr, "Sequence must be 1000 characters (only ");
    fprintf(stderr, "%d read).\n\n", slen);
    exit(1);
  }
  rcom_seq(seq, rseq, slen);
  if(do_training == 0) {
    fprintf(stderr, "%d bp read, %.2f pct GC\n", slen, tinf.gc*100.0);
    fflush(stderr);
  }
  else {
    fprintf(stderr, "%d bp seq created, %.2f pct GC\n", slen, tinf.gc*100.0);
    fflush(stderr);
  }

  /***************************************************************************
    Find all the potential starts and stops, sort them, and create a 
    comprehensive list of nodes for dynamic programming.
  ***************************************************************************/
  fprintf(stderr, "Locating all potential starts and stops..."); fflush(stderr);
  nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
  if(nn == 0) {
    fprintf(stderr, "0 nodes, exiting...\n");
    fflush(stderr);
    exit(0);
  }
  qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
  fprintf(stderr, "%d nodes\n", nn); fflush(stderr);

  if(do_training == 1 || (do_training == 0 && train_file == NULL)) {

    /***********************************************************************
      Scan all the ORFS looking for a potential GC bias in a particular 
      codon position.  This information will be used to acquire a good 
      initial set of genes.
    ***********************************************************************/
    fprintf(stderr, "Looking for GC bias in different frames...");
    fflush(stderr);
    gc_frame = calc_most_gc_frame(seq, slen);
    if(gc_frame == NULL) {
      fprintf(stderr, "Malloc failed on gc frame plot\n\n");
      exit(1);
    }
    record_gc_bias(gc_frame, nodes, nn, &tinf);
    fprintf(stderr, "frame bias scores: %.2f %.2f %.2f\n", tinf.bias[0],
            tinf.bias[1], tinf.bias[2]); fflush(stderr);

    /***********************************************************************
      Do an initial dynamic programming routine with just the GC frame bias
      used as a scoring function.  This will get an initial set of genes to
      train on. 
    ***********************************************************************/
    fprintf(stderr, "Building initial set of genes to train from...");
    fflush(stderr);
    record_overlapping_starts(nodes, nn, &tinf, 0);
    ipath = dprog(nodes, nn, &tinf, 0);
    fprintf(stderr, "done!\n"); fflush(stderr);

    /***********************************************************************
      Gather dicodon statistics for the training set.  Score the entire set
      of nodes.                               
    ***********************************************************************/
    fprintf(stderr, "Creating coding model and scoring nodes...");
    fflush(stderr);
    calc_dicodon_gene(&tinf, seq, rseq, slen, nodes, ipath);
    raw_coding_score(seq, rseq, slen, nodes, nn, &tinf);
    fprintf(stderr, "done!\n"); fflush(stderr);

    /***********************************************************************
      Determine if this organism uses Shine-Dalgarno or not and score the 
      nodes appropriately.
    ***********************************************************************/
    fprintf(stderr, "Examining upstream regions and training starts...");
    fflush(stderr);
    rbs_score(seq, rseq, slen, nodes, nn, &tinf);
    train_starts_sd(seq, rseq, slen, nodes, nn, &tinf);
    determine_sd_usage(&tinf);
    if(force_nonsd == 1) tinf.uses_sd = 0;
    if(tinf.uses_sd == 0) train_starts_nonsd(seq, rseq, slen, nodes, nn, &tinf);
    fprintf(stderr, "done!\n"); fflush(stderr);
  }

  if(do_training == 1) {
    fprintf(stderr, "Writing data to training file %s...", train_file);
    fflush(stderr);
    rv = write_training_file(train_file, &tinf);
    if(rv != 0) { fprintf(stderr, "error!\n"); exit(1); }
    else { fprintf(stderr, "done!\n"); exit(0); }
  }

  /***************************************************************************
    Second dynamic programming, using the dicodon statistics as the scoring
    function.                                
  ***************************************************************************/
  fprintf(stderr, "Building set of genes from statistics..."); fflush(stderr);

  score_nodes(seq, rseq, slen, nodes, nn, &tinf);
  if(start_file != NULL) write_start_file(start_file, nodes, nn, &tinf);
  record_overlapping_starts(nodes, nn, &tinf, 1);
  ipath = dprog(nodes, nn, &tinf, 1);
  eliminate_bad_genes(nodes, ipath, &tinf);
  ng = add_genes(genes, nodes, ipath);
  tweak_final_starts(genes, ng, nodes, nn, &tinf);
  fprintf(stderr, "done!\n"); fflush(stderr);

  /* Output the genes */
  print_genes(genes, ng, nodes, slen, output);
  if(trans_file != NULL)
    write_translations(trans_file, genes, ng, nodes, seq, rseq, slen, &tinf);
  exit(0);
}

void version() {
  fprintf(stderr, "\nProdigal V1.20: September, 2009\n\n");
  exit(0);
}

void usage(char *msg) {
  fprintf(stderr, "\n%s\n", msg);
  fprintf(stderr, "\nUsage:  prodigal [-a <trans_file>] [-c] [-g <tr_table>] ");
  fprintf(stderr, "[-h] [-m] [-n]\n");
  fprintf(stderr, "        [-o <output type>] [-s <start_file>] ");
  fprintf(stderr, "[-t <training_file>] [-v]\n");
  fprintf(stderr, "        < <input_file> > <output_file>\n\n");
  fprintf(stderr, "Do 'prodigal -h' for more information.\n\n");
  exit(1);
}

void help() {
  fprintf(stderr, "\nUsage:  prodigal [-a <trans_file>] [-c] [-g <tr_table>]");
  fprintf(stderr, " [-h] [-m] [-n]\n");
  fprintf(stderr, "        [-o <output type>] [-s <start_file>] ");
  fprintf(stderr, "[-t <training_file>] [-v]\n");
  fprintf(stderr, "        < <input_file> > <output_file>\n\n");
  fprintf(stderr, "         -h:  Print help menu.\n");
  fprintf(stderr, "         -a:  Write protein translations to the selected ");
  fprintf(stderr, "file.\n");
  fprintf(stderr, "         -c:  Closed ends.  Do not allow genes to run off ");
  fprintf(stderr, "edges.\n");
  fprintf(stderr, "         -g:  Specify a translation table to use (default");
  fprintf(stderr, " 11).\n");
  fprintf(stderr, "         -m:  Treat runs of n's as masked sequence and do");
  fprintf(stderr, " not build genes across \n");
  fprintf(stderr, "              them.\n");
  fprintf(stderr, "         -n:  Bypass the Shine-Dalgarno trainer and force");
  fprintf(stderr, " the program to scan\n");
  fprintf(stderr, "              for motifs.\n");
  fprintf(stderr, "         -o:  Select output format (gbk, gff, or sco)\n");
  fprintf(stderr, "         -s:  Write all potential genes (with scores) to");
  fprintf(stderr, " the selected file.\n");
  fprintf(stderr, "         -t:  Write a training file (if none exists); ");
  fprintf(stderr, "otherwise, read and use\n");
  fprintf(stderr, "              the specified training file.\n");
  fprintf(stderr, "         -v:  Print version number and exit.\n");
  fprintf(stderr, "\n         Single sequence :  prodigal < fasta.seq > ");
  fprintf(stderr, "my.genes\n");
  fprintf(stderr, "\n         Training on draft: cat fasta_seqs.* | prodigal");
  fprintf(stderr, "  -t my.trn\n");
  fprintf(stderr, "         Running on draft:  prodigal -t my.trn < ");
  fprintf(stderr, "fasta_seqs.1\n");
  fprintf(stderr, "                            prodigal -t my.trn < ");
  fprintf(stderr, "fasta_seqs.2 etc.\n");
  fprintf(stderr, "\n");
  exit(0);
}
