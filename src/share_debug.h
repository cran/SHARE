

#define long int



typedef struct node *Link1;
typedef struct node {long snp_pos; Link1 next;} Node;
typedef struct hap_model { long nhaps; long nsnps; double dev; long hapbase; long *snp_set; long *haplo_vec; double *hfreq; double *coef; } hmodel;




/********** function prototypes **********************/
 
void xshare(long *indx_subj,long *nsubj,long *nobs,long *subj_rep,long *csctl,long *nloci,long *nfold,long *nhap,long *hap1code,long *hap2code,long *uhap,double *happrob,double *wgt,long  *maxsnps, double *deviance, double *tol, long *verbose, long *phase, long *Minherit); 


void finalsubset(long *indx_subj,long *nsubj,long *nobs,long *subj_rep,long *csctl,long *nloci,long *nhap,long *hap1code,long *hap2code,long *uhap,double *happrob,double *wgt,long *maxsnps, long *bestsize, long *output_snp_set, long *out_haplo_vec, double *out_hap_freq, double *tol, long *phase, double *varstore, double *coef, long *out_nhaps, long *Minherit, long *verbose);

void shrink_phase_infer(long nsnp,long *nrecord,long *nhap,long *hap1code,long *hap2code,long *uhap,double *happrob,double *wgt,long nsubj,long *subj_rep, double tol, long *verbose,long ntestdata,long *test_hap1code,long *test_hap2code, FILE *file);  

void varest(hmodel *best_model,long *nsubj,long *nobs,long *subj_rep,long *csctl,long *nloci,long *nhap,long *hap1code,long *hap2code,long *uhap,double *happrob,double *wgt,long *snp_set, long nsnp_set, double cutoff, long Minherit, long *verbose,FILE *file, long *phase, double *varstore);

/*
  void extern_cross_val(long *nrecords,long *y,long *nsubj,long *subj_reps, long *haplotype_vec,long *nsnps,long *nhaps,long *hap1code,long *hap2code, double *freq, long *verbose, long *phase, double *dev); */

long **hap_shrink(long nrecords,long nsnps,long nhaps,long *fhaplotype_vec, double *hfreq,long *hap1code,long *hap2code,long *snp_set,long nsnp_set, long *out_nhaps,long *out_hap1code,long *out_hap2code, double *out_hfreq, long *verbose, FILE *file);


double haplo_dist(long *haplotype1, long *haplotype2,double *locus_freq,long nlocus);

void haplo_cluster(long nrecords,long *hap1code,long *hap2code,long **haplotype,double *freq,double min_freq,long nlocus,long *nhap,long *verbose, FILE *file);


void iwls_bin(long n,double **x,long ncov,long *y,double *weights,double *mustart,long maxit,double tol,double *coef,double *deviance,long *conv, long *verbose, FILE *file);  


hmodel *hap_shrink_reg(long *y,long  *haplotype_vec,long  nsnps,long  nhaps,long  *hap1code,long  *hap2code,double *weight,double *hfreq, long nrecords,long *nreps,long *snp_set,long nsnp_set,double cutoff,double *mustart,long maxit, long Minherit, double tol, long *verbose, FILE *file);  

double cross_val(hmodel *result_model,long nrecords,long *y,long nsubj,long *subj_reps, long *haplotype_vec,long nsnps,long nhaps,long *hap1code,long *hap2code, double *freq, long Minherit, long *verbose, long *phase, FILE *file);


void stepwise_search_alpha(long *indx_subj,long *nsubj,long *nobs, long *subj_rep,double *wgt, long *csctl,long *nloci,long *nhap,long *hap1code,long *hap2code,long *uhap,double *happrob,long *bestsize, long *Minherit, double *deviance, long  *maxsnp,double *tol,double *alpha,long *verbose,long *phase, long *output_snp_set, long *out_haplo_vec, double *out_hap_freq, double *varstore, double *coef, long *out_nhaps);


long *insertionSort(long *numbers, long array_size, long *order);

double **double_vec_to_mat(double *Yvec, long nrow, long ncol);

long **long_vec_to_mat(long *Yvec, long nrow, long ncol);

double **double_matrix(long nrow, long ncol);

long **long_matrix(long nrow, long ncol);

double *double_vec(long n);

long *long_vec(long n);

long *long_mat_to_vec(long **Ymat, long nrow, long ncol);

static void errmsg(char *string);

void dqr_(double *, long *, long *, double *, long *, double *, double *, double *, double *,long *,long *, double *, double *);

long max_long(long x, long y);
long min_long(long x, long y);

double max_double(double x, double y);
double min_double(double x, double y);

void print_matrix_long(long **m, long nrow, long ncol, FILE *file);
void print_matrix_double(double **m, long nrow, long ncol, FILE *file);
void print_vector_double(double *m, long n, FILE *file);
void print_vector_long(long *m, long n, FILE *file);
void print_list_long(Node m, FILE *file);
void print_hmodel(hmodel *m, FILE *file);


hmodel *new_hmodel(long nsnps,long nhaps);

void copy_hmodel(hmodel *current, hmodel *best);

void Free_hmodel(hmodel *x);

void insertSort(long *numbers, long array_size);
