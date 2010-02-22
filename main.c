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

#include "sequence.h"
#include "metagenomic.h"
#include "node.h"
#include "dprog.h"
#include "gene.h"
  #include <sys/time.h>
#define MIN_SINGLE_GENOME 1000

void version();
void usage(char *);
void help();

int main(int argc, char *argv[]) {

  int rv, slen, nn, ng, i, ipath, *gc_frame, do_training, output, max_phase;
  int closed, do_mask, nmask, force_nonsd, user_tt, is_meta, num_seq;
  double max_score;
  unsigned char *seq, *rseq;
  char *train_file, *start_file, *trans_file; 
  char *input_file, *output_file;
  FILE *input_ptr, *output_ptr, *start_ptr, *trans_ptr;
  struct _node *nodes;
  struct _gene *genes;
  struct _training tinf;
  struct _metagenomic_bin meta[30];
  mask mlist[MAX_MASKS];

  /* Allocate memory and initialize variables */
  seq = (unsigned char *)malloc(MAX_SEQ*sizeof(unsigned char));
  rseq = (unsigned char *)malloc(MAX_SEQ*sizeof(unsigned char));
  nodes = (struct _node *)malloc(MAX_NODES*sizeof(struct _node));
  genes = (struct _gene *)malloc(MAX_GENES*sizeof(struct _gene));
  if(seq == NULL || rseq == NULL || nodes == NULL || genes == NULL) {
    fprintf(stderr, "\nError: Malloc failed on sequence/orfs\n\n"); exit(1);
  }
  memset(seq, 0, MAX_SEQ*sizeof(unsigned char));
  memset(rseq, 0, MAX_SEQ*sizeof(unsigned char));
  memset(nodes, 0, MAX_NODES*sizeof(struct _node));
  memset(genes, 0, MAX_GENES*sizeof(struct _gene));
  memset(&tinf, 0, sizeof(struct _training));
  for(i = 0; i < 30; i++) {
    memset(&meta[i], 0, sizeof(struct _metagenomic_bin));
    meta[i].tinf = (struct _training *)malloc(sizeof(struct _training));
    if(meta[i].tinf == NULL) {
      fprintf(stderr, "\nError: Malloc failed on training structure.\n\n"); 
      exit(1);
    }
    memset(meta[i].tinf, 0, sizeof(struct _training));
  }
  nn = 0; slen = 0; ipath = 0; ng = 0; nmask = 0;
  user_tt = 0; is_meta = 0; num_seq = 0; 
  max_phase = 0; max_score = -100.0;
  train_file = NULL; do_training = 0;
  start_file = NULL; trans_file = NULL;
  start_ptr = stdout; trans_ptr = stdout;
  input_file = NULL; output_file = NULL;
  input_ptr = stdin; output_ptr = stdout;
  output = 0; closed = 0; do_mask = 0; force_nonsd = 0;

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
       strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-F") == 0 ||
       strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "-S") == 0 ||
       strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "-I") == 0 ||
       strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "-O") == 0 ||
       strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-P") == 0))
      usage("-a/-f/-g/-i/-o/-p/-s options require parameters.");
    else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "-C") == 0)
      closed = 1;
    else if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "-M") == 0)
      do_mask = 1;
    else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "-N") == 0)
      force_nonsd = 1;
    else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-H") == 0) help();
    else if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "-V") == 0) version();
    else if(strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "-A") == 0) {
      trans_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "-I") == 0) {
      input_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "-O") == 0) {
      output_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "-S") == 0) {
      start_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "-T") == 0) {
      train_file = argv[i+1];
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
    else if(strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-P") == 0) {
      if(strncmp(argv[i+1], "0", 1) == 0 || strcmp(argv[i+1], "single") == 0 ||
         strcmp(argv[i+1], "SINGLE") == 0) {
        is_meta = 0;
      }
      else if(strncmp(argv[i+1], "0", 1) == 0 || strcmp(argv[i+1], "meta") == 0 ||
         strcmp(argv[i+1], "META") == 0) {
        is_meta = 1;
      }
      else usage("Invalid meta/single genome type specified.");
      i++;
    }
    else if(strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-F") == 0) {
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
  fprintf(stderr, "PRODIGAL v2.00 [March, 2010]         \n");
  fprintf(stderr, "Univ of Tenn / Oak Ridge National Lab\n");
  fprintf(stderr, "Doug Hyatt, Loren Hauser, et al.     \n");
  fprintf(stderr, "-------------------------------------\n");

  /* Read in the training file (if specified) */
  if(train_file != NULL) {
    if(is_meta == 1) {
      fprintf(stderr, "\nError: cannot specify metagenomic sequence with a");
      fprintf(stderr, " training file.\n");
      exit(2);
    } 
    rv = read_training_file(train_file, &tinf);
    if(rv == 1) do_training = 1;
    else {
      if(force_nonsd == 1) { 
        fprintf(stderr, "\nError: cannot force non-SD finder with a training");
        fprintf(stderr, " file already created!\n"); exit(3);
      }
      fprintf(stderr, "Reading in training data from file %s...", train_file);
      if(user_tt > 0 && user_tt != tinf.trans_table) { 
        fprintf(stderr, "Warning: user-specified translation table does not ");
        fprintf(stderr, "match the one in the specified training file! \n");
        fflush(stderr);
      }
      if(rv == -1) { 
        fprintf(stderr, "\n\nError: training file did not read correctly!\n"); 
        exit(4); 
      }
      fprintf(stderr, "done!\n"); 
      fprintf(stderr, "-------------------------------------\n");
      fflush(stderr);
    }
  }

  /* Check i/o files (if specified) and prepare them for reading/writing */
  if(input_file != NULL) {
    input_ptr = fopen(input_file, "r");
    if(input_ptr == NULL) {
      fprintf(stderr, "\nError opening input file %s.\n\n", input_file);
      exit(5);
    }
  }
  if(output_file != NULL) {
    output_ptr = fopen(output_file, "w");
    if(output_ptr == NULL) {
      fprintf(stderr, "\nError opening output file %s.\n\n", output_file);
      exit(6);
    }
  }
  if(start_file != NULL) {
    start_ptr = fopen(start_file, "w");
    if(start_ptr == NULL) {
      fprintf(stderr, "\nError opening start file %s.\n\n", start_file);
      exit(7);
    }
  }
  if(trans_file != NULL) {
    trans_ptr = fopen(trans_file, "w");
    if(trans_ptr == NULL) {
      fprintf(stderr, "\nError opening translation file %s.\n\n", trans_file);
      exit(8);
    }
  }

  /***************************************************************************
    Single Genome Training:  Read in the sequence(s) and perform the
    training on them.
  ***************************************************************************/
  if(is_meta == 0 && (do_training == 1 || (do_training == 0 && train_file == 
     NULL))) {
    fprintf(stderr, "Request:  Single Genome, Phase:  Training\n");
    fprintf(stderr, "Reading in the sequence(s) to train..."); fflush(stderr);
    slen = read_seq_single(input_ptr, seq, &(tinf.gc), do_mask, mlist, 
                           &nmask);
    if(slen == 0) {
      fprintf(stderr, "\n\nSequence read failed (file must be Fasta, ");
      fprintf(stderr, "Genbank, or EMBL format).\n\n");
      exit(9);
    }
    if(slen < MIN_SINGLE_GENOME) {
      fprintf(stderr, "\n\nSequence must be %d", MIN_SINGLE_GENOME);
      fprintf(stderr, " characters (only %d read).\n(Consider", slen);
      fprintf(stderr, " running with the -p meta option.)\n\n");
      exit(10);
    }
    rcom_seq(seq, rseq, slen);
    fprintf(stderr, "%d bp seq created, %.2f pct GC\n", slen, tinf.gc*100.0);
    fflush(stderr);

    /***********************************************************************
      Find all the potential starts and stops, sort them, and create a 
      comprehensive list of nodes for dynamic programming.
    ***********************************************************************/
    fprintf(stderr, "Locating all potential starts and stops..."); 
    fflush(stderr);
    nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
    if(nn == 0) {
      fprintf(stderr, "0 nodes, exiting...\n");
      fflush(stderr);
      exit(0);
    }
    qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
    fprintf(stderr, "%d nodes\n", nn); fflush(stderr);

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
      exit(11);
    }
    record_gc_bias(gc_frame, nodes, nn, &tinf);
    fprintf(stderr, "frame bias scores: %.2f %.2f %.2f\n", tinf.bias[0],
            tinf.bias[1], tinf.bias[2]); fflush(stderr);
    free(gc_frame);

    /***********************************************************************
      Do an initial dynamic programming routine with just the GC frame
      bias used as a scoring function.  This will get an initial set of 
      genes to train on. 
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

    /* If training specified, write the training file and exit. */
    if(do_training == 1) {
      fprintf(stderr, "Writing data to training file %s...", train_file);
      fflush(stderr);
      rv = write_training_file(train_file, &tinf);
      if(rv != 0) { 
        fprintf(stderr, "\nError: could not write training file!\n"); 
        exit(12); 
      }
      else { fprintf(stderr, "done!\n"); exit(0); }
    }

    /* Rewind input file */    
    fprintf(stderr, "-------------------------------------\n");
    if(fseek(input_ptr, 0, SEEK_SET) == -1) {
      fprintf(stderr, "\nError: could not rewind input file.\n"); exit(13);
    }

    /* Reset all the sequence/dynamic programming variables */
    memset(seq, 0, slen*sizeof(unsigned char));
    memset(rseq, 0, slen*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0; nmask = 0;
  }

  /* Initialize the training files for a metagenomic request */
  else if(is_meta == 1) {
    fprintf(stderr, "Request:  Metagenomic, Phase:  Training\n");
    fprintf(stderr, "Initializing training files...");
    initialize_metagenomic_bins(meta);
    fprintf(stderr, "done!\n");
    fprintf(stderr, "-------------------------------------\n");
    fflush(stderr);
  }

  /* Print out header for gene finding phase */
  if(is_meta == 1)
    fprintf(stderr, "Request:  Metagenomic, Phase:  Gene Finding\n");
  else fprintf(stderr, "Request:  Single Genome, Phase:  Gene Finding\n");
  fflush(stderr);

  /* Read and process each sequence in the file in succession */
  while((slen = next_seq_multi(input_ptr, seq, &num_seq, &(tinf.gc), do_mask, 
                               mlist, &nmask)) != -1) {
    rcom_seq(seq, rseq, slen);
    if(slen == 0) {
      fprintf(stderr, "\nSequence read failed (file must be Fasta, ");
      fprintf(stderr, "Genbank, or EMBL format).\n\n");
      exit(14);
    }

    fprintf(stderr, "Finding genes in sequence #%d (%d bp)...", num_seq, slen);
    fflush(stderr);

    if(is_meta == 0) { /* Single Genome Version */

      /***********************************************************************
        Find all the potential starts and stops, sort them, and create a 
        comprehensive list of nodes for dynamic programming.
      ***********************************************************************/
      nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
      qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

      /***********************************************************************
        Second dynamic programming, using the dicodon statistics as the
        scoring function.                                
      ***********************************************************************/
      score_nodes(seq, rseq, slen, nodes, nn, &tinf);
      if(start_ptr != stdout) 
        write_start_file(start_ptr, nodes, nn, &tinf, num_seq);
      record_overlapping_starts(nodes, nn, &tinf, 1);
      ipath = dprog(nodes, nn, &tinf, 1);
      eliminate_bad_genes(nodes, ipath, &tinf);
      ng = add_genes(genes, nodes, ipath);
      tweak_final_starts(genes, ng, nodes, nn, &tinf);
      fprintf(stderr, "done!\n"); fflush(stderr);

    }

    else { /* Metagenomic Version */

      determine_top_bins(seq, rseq, slen, tinf.gc, meta);

      /***********************************************************************
        We do the dynamic programming multiple times on metagenomic fragments,
        once for each of the top model organisms in our training set.
      ***********************************************************************/
      for(i = 0; i < NUM_BIN; i++) {
        if(i == 0 || meta[i].tinf->trans_table != 
           meta[i-1].tinf->trans_table) {
          memset(nodes, 0, nn*sizeof(struct _node));
          nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, 
                         meta[i].tinf);
          qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
        }
        reset_node_scores(nodes, nn);
        score_nodes(seq, rseq, slen, nodes, nn, meta[i].tinf);
        record_overlapping_starts(nodes, nn, meta[i].tinf, 1);
        ipath = dprog(nodes, nn, meta[i].tinf, 1);
        if(i == 0 || nodes[ipath].score > max_score) {
          max_phase = i;
          max_score = nodes[ipath].score;
          eliminate_bad_genes(nodes, ipath, meta[i].tinf);
          ng = add_genes(genes, nodes, ipath);
          tweak_final_starts(genes, ng, nodes, nn, meta[i].tinf);
        }
      }    

      fprintf(stderr, "done!\n"); fflush(stderr);

      /* Recover the nodes for the best of the runs */
      if(start_ptr != stdout || meta[max_phase].tinf->trans_table !=
         meta[i-1].tinf->trans_table) {
        memset(nodes, 0, nn*sizeof(struct _node));
        nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask,
                       meta[max_phase].tinf);
        qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
        score_nodes(seq, rseq, slen, nodes, nn, meta[max_phase].tinf);
        if(start_ptr != stdout) 
          write_start_file(start_ptr, nodes, nn, meta[max_phase].tinf, 
                           num_seq);
      }
    }

    /* Output the genes */
    print_genes(output_ptr, genes, ng, nodes, slen, output, num_seq);
    fflush(stdout);
    if(trans_ptr != stdout) {
      if(is_meta == 0)
        write_translations(trans_ptr, genes, ng, nodes, seq, rseq, slen, &tinf, 
                           num_seq);
      else
        write_translations(trans_ptr, genes, ng, nodes, seq, rseq, slen,
                           meta[max_phase].tinf, num_seq);
    }

    /* Reset all the sequence/dynamic programming variables */
    memset(seq, 0, slen*sizeof(unsigned char));
    memset(rseq, 0, slen*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0; nmask = 0;
  }

  /* Free all memory */
  if(seq != NULL) free(seq);
  if(rseq != NULL) free(rseq);
  if(nodes != NULL) free(nodes);
  if(genes != NULL) free(genes);
  for(i = 0; i < 30; i++) if(meta[i].tinf != NULL) free(meta[i].tinf);

  /* Close all the filehandles and exit */
  if(input_ptr != stdin) fclose(input_ptr);
  if(output_ptr != stdout) fclose(output_ptr);
  if(start_ptr != stdout) fclose(start_ptr);
  if(trans_ptr != stdout) fclose(trans_ptr);
  exit(0);
}

void version() {
  fprintf(stderr, "\nProdigal V2.00: March, 2010\n\n");
  exit(0);
}

void usage(char *msg) {
  fprintf(stderr, "\n%s\n", msg);
  fprintf(stderr, "\nUsage:  prodigal [-a trans_file] [-c] [-f output_type] [-g tr_table]");
  fprintf(stderr, " [-h]\n");
  fprintf(stderr, "                 [-i input_file] [-m] [-n] [-o output_file] [-p ");
  fprintf(stderr, "mode]\n");
  fprintf(stderr, "                 [-s start_file] [-t training_file] [-v]\n");
  fprintf(stderr, "\nDo 'prodigal -h' for more information.\n\n");
  exit(15);
}

void help() {
  fprintf(stderr, "\nUsage:  prodigal [-a trans_file] [-c] [-f output_type] [-g tr_table]");
  fprintf(stderr, " [-h]\n");
  fprintf(stderr, "                 [-i input_file] [-m] [-n] [-o output_file] [-p ");
  fprintf(stderr, "mode]\n");
  fprintf(stderr, "                 [-s start_file] [-t training_file] [-v]\n");
  fprintf(stderr, "\n         -a:  Write protein translations to the selected ");
  fprintf(stderr, "file.\n");
  fprintf(stderr, "         -c:  Closed ends.  Do not allow genes to run off ");
  fprintf(stderr, "edges.\n");
  fprintf(stderr, "         -f:  Select output format (gbk, gff, or sco).  ");
  fprintf(stderr, "Default is gbk.\n");
  fprintf(stderr, "         -g:  Specify a translation table to use (default");
  fprintf(stderr, " 11).\n");
  fprintf(stderr, "         -h:  Print help menu and exit.\n");
  fprintf(stderr, "         -i:  Specify input file (default reads from ");
  fprintf(stderr, "stdin).\n");
  fprintf(stderr, "         -m:  Treat runs of n's as masked sequence and do");
  fprintf(stderr, " not build genes across \n");
  fprintf(stderr, "              them.\n");
  fprintf(stderr, "         -n:  Bypass the Shine-Dalgarno trainer and force");
  fprintf(stderr, " the program to scan\n");
  fprintf(stderr, "              for motifs.\n");
  fprintf(stderr, "         -o:  Specify output file (default writes to ");
  fprintf(stderr, "stdout).\n");
  fprintf(stderr, "         -p:  Select procedure (single or meta).  Default");
  fprintf(stderr, " is single.\n");
  fprintf(stderr, "         -s:  Write all potential genes (with scores) to");
  fprintf(stderr, " the selected file.\n");
  fprintf(stderr, "         -t:  Write a training file (if none exists); ");
  fprintf(stderr, "otherwise, read and use\n");
  fprintf(stderr, "              the specified training file.\n");
  fprintf(stderr, "         -v:  Print version number and exit.\n");
  fprintf(stderr, "\n         Single FASTA sequence, finished genome:\n");
  fprintf(stderr, "           prodigal -i fasta.seq -o my.genes\n");
  fprintf(stderr, "         Single FASTA sequence, metagenomic:\n");
  fprintf(stderr, "           prodigal -p meta -i fasta.seq -o my.genes\n");
  fprintf(stderr, "         Multiple FASTA sequence, single genome:\n");
  fprintf(stderr, "           prodigal -i fasta.seqs -o my.genes\n");
  fprintf(stderr, "         Multiple FASTA sequence, metagenomic:\n");
  fprintf(stderr, "           prodigal -p meta -i fasta.seqs -o my.genes\n");
  fprintf(stderr, "\n         Training on separate files: cat fasta_seqs.* ");
  fprintf(stderr, " | prodigal -t my.trn\n");
  fprintf(stderr, "         Running on separate files:  prodigal -t my.trn");
  fprintf(stderr, " -i fasta_seqs.1\n");
  fprintf(stderr, "                                     prodigal -t my.trn");
  fprintf(stderr, " -i fasta_seqs.2 etc.\n");
  fprintf(stderr, "\n");
  exit(0);
}
