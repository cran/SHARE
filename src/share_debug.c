
/* This code is to execute the cross validation + model searching step in the SHARE algorithm */
/* Author: James Dai */

// 2010-07-15 : Valerie Obenchain
// 			- Added non-genetic covariate handling



# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <R.h>
# include "share_debug.h"  
/*# include "fortran.h"   */
/* static hmodel *final_model; */

   
/* this function takes input from a phased genotype dataset with case-control status (after random permutation, output the the prediction deviances of a ladder of models up to size maxsnps  */

void xshare(long *indx_subj,long *nsubj,long *nobs,long *subj_rep,long *csctl,long *nloci,long *nfold,long *nhap,long *hap1code,long *hap2code,long *uhap,double *happrob,double *wgt,long  *maxsnps,
long *nngcov, double *ngcov_vec,
double *deviance,double *tol,long *verbose, long *phase, long *Minherit)
{

  /**   
	long *indx_subj              vector of subject ids, length nobs  
	long *nsubj,                 total number of distinct subjects  
	long *nobs,                  total number of observations   
	long *subj_rep,              vector of number of possible pairs for each subject, length nsubj  
	long *csctl,                 vector of case-control status, length nobs  
	long *nloci,                 total number of SNPs considered  
	long *nfold,                 the number of cross-validation fold  
	long *nhap,                  number of haplotypes phased using all data , to start with 
	long *hap1code,              the vector of codes for haplotype 1, length nobs  
	long *hap2code,              the vector of codes for haplotype 2, length nobs  
	long *uhap,                  vector of unique haplotype  
	double *happrob,             the vector of haplotype probabilities  
	double *wgt,                 the posterior probabilities of each pair of haplotypes, length nobs     
	long  *maxsnps,              maximal number of SNPs to be considered in the selected set  
	long *nngcov,                the number of non-genetic covariates  
	double *ngcov_vec,           vector of non-genetic covariates  
	double *deviance,            the vector of deviances from cross-validations  
	double *tol                  the convergence parameter
	long *verbose		     		  the indicator for print-out when debugging
	long *phase		     			  the indicator whether phase is known
	long *Minherit               the integer indexing the mode of inheritance: 1-additive; 2-dominant; 3-recessive
  **/

     long ncross, i, j, k, count, count1, count2, nextra, *cvlabel,*cvlabel2, *nobs_train, train_sample_size, maxit=1000, test_nsub, full_nhap;
     long *train_hap1code, *train_hap2code, *train_csctl, *test_size, *train_subj_rep, *test_hap1code, *test_hap2code, *test_csctl, *test_subj_rep, *snp_set, *current_snp_set;
     double *train_wgt, *test_wgt, **train_ngcov, **test_ngcov, obssum, wgtsum, meany, temp_phi, min_phi, *mustart, **dev_mat, min_freq, **ngcov_mat;
     hmodel *temp, *base_model, *best_model; 
     Node SNP_out;
     Link1 current_snp= &SNP_out, next_snp;   

     FILE *file; 
     file = fopen("debug_xshare.txt", "w"); 

// Remainder subjects (nextra) are included in last fold 
     ncross=floor(*nsubj/ *nfold);
     nextra=*nsubj - ncross*(*nfold);
     count1=0;
     count2=0;


     if (*verbose==1){ 
       fprintf(file,"nsub= %d\n",*nsubj);
       fprintf(file,"nobs= %d\n",*nobs);
       fprintf(file,"ncross= %i\t",ncross);
       fprintf(file,"nextra= %i\n",nextra);
     }
 
// Create indices for cv-folds :
// cvlabel identifies which observations are in which fold
// cvlabel2 identifies which subject is in which fold
// here count1 is observations, count2 is subjects

     cvlabel=long_vec(*nobs);
     cvlabel2=long_vec(*nsubj);
     test_size=long_vec(*nfold);
     nobs_train=(long *) malloc(sizeof(long));

     for (i=0;i<*nfold;i++) test_size[i]=0;

     for (i=0;i<(*nfold-1);i++) {
        for (j=0;j<ncross;j++) {
           for (k=0;k<subj_rep[count1];k++){
              cvlabel[count2]=i+1;
              count2 ++;            
           }        
           test_size[i] +=subj_rep[count1];
           cvlabel2[count1]=i+1;
           count1++;
        }
     } 
  
     for (j=0;j<(nextra+ncross);j++) {   
         for (k=0;k<subj_rep[count1];k++) {
             cvlabel[count2]=*nfold;
             count2 ++;
         }  
         test_size[*nfold-1] +=subj_rep[count1];
         cvlabel2[count1]=*nfold;
         count1++; 
     }


     if (*verbose==1) { 
       fprintf(file,"\ntest size = (obs) :\t");
       print_vector_long(test_size,*nfold, file);
       fprintf(file,"\ncvlabel2 (subs) :\n");
       print_vector_long(cvlabel2,count1,file);
       fprintf(file,"\ncvlabel (obs) :\n");
       print_vector_long(cvlabel,count2,file);
     }

     dev_mat = double_matrix(*nfold,*maxsnps*2);
     full_nhap=*nhap;

// put covariate data in matrix
	if (*nngcov>1) { 
		ngcov_mat = double_vec_to_mat(ngcov_vec, *nobs, *nngcov); 
		}

     if (*verbose==1 && *nngcov!=0) {
       fprintf(file,"\n\nngcov_vec : \n");
       print_vector_double(ngcov_vec,*nobs*2,file);
       fprintf(file,"\nngcov_mat :\n");
       print_matrix_double(ngcov_mat,*nobs,*nngcov,file); 
		}

// Start n-fold CV
	// 1. subset data on fold
	// 2. estimate new hap probabilities for training set
	// 3. compute null deviance
	// 4. build model
	// 5. prune model
	// 6. sum deviances across folds for plot 


     for (i=0;i<*nfold;i++) {        

       if (*verbose==1) fprintf(file,"\nthe iteration number in the cross valid=%i\n\n", i);    

			// Containers for observations   

			// Non-genetic covariates  
        if (*nngcov!=0)  {
		       train_ngcov= double_matrix(*nobs - test_size[i], *nngcov);        
		       test_ngcov= double_matrix(test_size[i], *nngcov);
        }

         train_hap1code=long_vec(*nobs - test_size[i]);
         train_hap2code=long_vec(*nobs - test_size[i]);
         train_csctl=long_vec(*nobs - test_size[i]);
         train_wgt= double_vec(*nobs - test_size[i]);        

         test_hap1code=long_vec(test_size[i]);
         test_csctl= long_vec(test_size[i]);
         test_hap2code=long_vec(test_size[i]);     
         test_wgt= double_vec(test_size[i]);


			// Containers for subject numbers 
			// identify train and test subjects in this fold, and no. of train observations

         if (i<*nfold-1) { 
                    train_subj_rep=long_vec(*nsubj - ncross);
                    test_subj_rep=long_vec(ncross); 
                    train_sample_size=*nsubj - ncross;
                  
                    if (*verbose==1) fprintf(file,"\nthe number of train subjects=%i\n\n", train_sample_size);  

  			// on the last fold, pick up nextra subjects                  
         } else {
                    train_subj_rep=long_vec(*nsubj - ncross-nextra);
                    test_subj_rep=long_vec(ncross+nextra); 
                  
                    train_sample_size=*nsubj - ncross - nextra;
         }
      
// Subset observations on fold 

			// haplo codes, weights and case/control status on fold and train/test 
         count1=0;
         count2=0;
         for (j=0;j<*nobs;j++)  {
                 if (cvlabel[j]!=i+1)  {
        
                          train_hap1code[count1] = hap1code[j];
                          train_hap2code[count1] = hap2code[j];
                          train_csctl[count1] = csctl[j];
                          train_wgt[count1] = wgt[j];
                          count1 ++ ;
                       
		 } else {
                          test_hap1code[count2] = hap1code[j];
                          test_hap2code[count2] = hap2code[j];                   
                          test_wgt[count2] = wgt[j];   
                          test_csctl[count2] = csctl[j];
                          count2 ++ ;
                 }         
         }


		// non-genetic covariates
       if (*nngcov!=0) {
       	for (k=0;k<*nngcov;k++)  {
       		count1=0;
       		count2=0;
       		for (j=0;j<*nobs;j++)  {
       			if (cvlabel[j]!=i+1)  {
       				train_ngcov[count1][k] = ngcov_mat[j][k]; 
       				count1 ++ ;
       			} else {
       				test_ngcov[count2][k] = ngcov_mat[j][k]; 
       				count2 ++ ;
       			}
       		}         
       	}
       }


    if (*verbose==1) {
      fprintf(file,"\ncount1(# train obs), count2(# test obs) =%i %i \t\n", count1, count2);
      fprintf(file,"\nthe last test wgt=%8.5f\n\n", test_wgt[count2-1]);
      fprintf(file,"\ntrain_hap1code initial : \n");
      print_vector_long(train_hap1code,*nobs - test_size[i],file); 
      if (*nngcov != 0) {fprintf(file,"\n test_ngcov : \n");
      print_matrix_double(test_ngcov, test_size[i], *nngcov, file); }}


      // Identify number of observations for each subject  
      count1=0;
      count2=0;

	  	for (j=0; j<*nsubj;j++) {
	   	if (cvlabel2[j]!=i+1) {
         	train_subj_rep[count1]=subj_rep[j]; 
         	count1++;
         } else {
         	test_subj_rep[count2]=subj_rep[j];                       
         	count2++;
         }
      }

       count=0;
       for (j=0;j<count1;j++) count += train_subj_rep[j];



    if (*verbose==1)    {
         fprintf(file,"\ncount1(# train subs), count2(# test subs) =%i %i \t\n", count1, count2);
         fprintf(file,"\nthe sum of train reps =%i\n", count);
         fprintf(file,"\nthe last test wgt=%8.5f\n", test_wgt[count2-1]);
         fprintf(file,"\ntest case ctl vector\n\n");
         print_vector_long(test_csctl,count2,file); }



// Estimate the new haplotype probability for training set 
  
         *nobs_train= *nobs-test_size[i];
         
         if (*verbose==1)    fprintf(file,"\nthe number of training data=%i\n\n", *nobs_train);  
 
         *nhap=full_nhap; 

         if (*phase==0) shrink_phase_infer(
				*nloci,nobs_train,nhap,train_hap1code,train_hap2code, uhap,happrob, 
				train_wgt, train_sample_size,train_subj_rep, *tol, verbose,test_size[i],
				test_hap1code, test_hap2code,file);
 
      	if (*verbose==1)   {
      	         fprintf(file,"\nAfter shrink_phase_infer :\n");
      	         fprintf(file,"\ntrain_hap1code :\n");
      	         print_vector_long(train_hap1code,*nobs_train,file);
      	         fprintf(file,"\nthe number of training data=%i\n\n", *nobs_train);
      	         fprintf(file,"\ntrain_subj_rep :\n");
      	         print_vector_long(train_subj_rep,train_sample_size,file);
      	  }

// Compute the null deviance  
         obssum=0;
         wgtsum=0;

         for (j=0; j<(*nobs - test_size[i]); j++) {
                   obssum += train_csctl[j]*train_wgt[j];
                   wgtsum += train_wgt[j];
	 		}

         meany = obssum/wgtsum;
         for (j=0;j<test_size[i];j++) {  
				dev_mat[i][0] += (test_csctl[j]==1) ? test_wgt[j]*log(meany) :  test_wgt[j]*log(1-meany);
   		}
 
         if (*verbose==1)  fprintf(file,"\nthe null deviance=%8.5f\n\n", dev_mat[i][0]);   
        
// Prepare the linked list of all snps to start the search 
         snp_set = long_vec(*maxsnps);
         SNP_out.next=NULL;

         for (j=0;j<*nloci;j++) {
                  next_snp = malloc(sizeof(Link1));
                  next_snp->snp_pos = j+1;
                  next_snp->next = NULL;
                  current_snp -> next = next_snp;
                  current_snp=next_snp;
	 		}

             
          if (*verbose==1) {
               fprintf(file,"\nSNP_out :\t\n");
               print_list_long(SNP_out,file);
             }

// Model searching 

			// given the limit of maxsnps,  
         // 1. compute mustart and min_freq 
         // 2. traverse the linked list to select a snp and subset the data on that snp
         // 3. call hap_shrink_reg to compute the deviance of the train model w/ those specific snps
         // 4. select best model for that size based on deviance
         // 5. cross-val on that model size
         // 6. repeat until have grown the model to size = maxsnp  

         // count variable tracks all snps , count1 tracks snps in snp_out
         
         count=0;
         mustart=double_vec(*nobs_train);

         for (j=0;j<*nobs_train;j++)  { mustart[j]= (train_wgt[j]*train_csctl[j]+0.5)/(train_wgt[j]+1);}

         min_freq= 2.0/(2.0* (*nsubj));
        
         if (*verbose==1)    fprintf(file,"\nthe freq cutoff = %8.9f",min_freq);

         // Build model one snp at a time up to maxsnps:
         while (count < *maxsnps) {
             current_snp = &SNP_out;
             count1=0;

             while (current_snp->next !=NULL) {
                 next_snp=current_snp -> next;
                 snp_set[count]= next_snp -> snp_pos;
                 if (*verbose==1)     fprintf(file,"\nthe current snp set :\t\n");
                 if (*verbose==1)    print_vector_long(snp_set,count+1, file);

         // call hap_shrink_reg : 
         // 1. call hap_shrink (recode haplos)
         // 2. call haplo_cluster
         // 3. find haplotype with largest frequency
         // 4. call iwls_bin (get least squares solution)
         // 5. output result_model

         temp=hap_shrink_reg(
				train_csctl, uhap,*nloci, *nhap,train_hap1code,train_hap2code,
				train_wgt,happrob,*nobs_train,train_subj_rep,snp_set, count+1,
				min_freq, mustart,
				nngcov, train_ngcov,
				maxit, *Minherit,*tol, verbose,file); 

         // Compare deviance in current model with previous best 
         
         if  (count==0) {
         	if (count1==0) {
         		best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
         		copy_hmodel(temp,best_model);
         	}
         	else if  (best_model->dev < temp -> dev)  {
         	Free_hmodel(best_model);
         	best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
         	copy_hmodel(temp,best_model);
         	}
         	              
         } else {
         	if (count1==0) {
         		best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
         		copy_hmodel(temp,best_model);
         		min_phi =  (temp->nhaps == base_model->nhaps) ? 0 : (temp -> dev - base_model->dev )/ (temp->nhaps - base_model->nhaps);  
         	} else {
         		temp_phi= (temp->nhaps == base_model->nhaps) ? 0 : (temp->dev-base_model->dev )/ (temp->nhaps - base_model->nhaps);  
         		if (temp_phi > min_phi) {
         			Free_hmodel(best_model);
         			best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
         			copy_hmodel(temp,best_model);
         			min_phi=temp_phi;    
         		}
         	}
         }
         


        if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 
	     if (*verbose==1)  {
				fprintf(file,"\n\nbest_model : \n");
				 print_hmodel(best_model,file);
			}
       
        // Increment the snp and build the next model of size count 
        current_snp = next_snp;
		  Free_hmodel(temp);                 
        count1++;
	 }                  


        // Deviance comparison of last model of size count 
        if (count==0) { 
                base_model = new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
                copy_hmodel(best_model,base_model);
        } else { 
                Free_hmodel(base_model);
                base_model = new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
                copy_hmodel(best_model,base_model);
        }            
  
        if (*verbose==1)  print_hmodel(best_model,file);
        if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 

        snp_set[count]= best_model->snp_set[count];
        if (*verbose==1)   print_vector_long(snp_set,count+1,file);
 
        if (i!=*nfold) test_nsub = ncross; else test_nsub = ncross + nextra;
       

        // Call cv function to compute deviance on model of size count w/ train data 
	     dev_mat[i][count+1]=cross_val(
				best_model,test_size[i], test_csctl,test_nsub, test_subj_rep, uhap,
				*nloci,*nhap,test_hap1code,test_hap2code, happrob,
				nngcov, test_ngcov,
				*Minherit,verbose,phase,file); 

              Free_hmodel(best_model);
           

			// Delete the selected SNPs from SNP_out (don't want same snp in model twice) 
         if (*verbose==1) {
         fprintf(file,"\n delete snps from snp_out: before = \n");
         print_list_long(SNP_out,file);}

         current_snp = &SNP_out;
         while (current_snp->next !=NULL) {                 
                  next_snp=current_snp -> next;  
                  if (next_snp -> snp_pos ==  snp_set[count])  {
                     current_snp ->next = next_snp->next;
                     free(next_snp);
                  } else {
                     current_snp = next_snp;
                  }
         }
         //if (*verbose==1)     print_list_long(SNP_out,file);
         if (*verbose==1) {
         fprintf(file,"\n delete snps from snp_out: after = \n");
         print_list_long(SNP_out,file);}

         // increment count and move to next model size
			count++;
	 }

// Pruning 
        if (*verbose==1) fprintf(file,"\nstart pruning \n") ;
        while (count>1) {

          current_snp_set = long_vec(count-1);
          count1=0;
          for (j=0;j<count;j++) {

                for (k=0;k<(count-1);k++)  {
                          if (k<j) current_snp_set[k] = snp_set[k];  
                          else current_snp_set[k] = snp_set[k+1];  
                }

				temp=hap_shrink_reg(
					train_csctl, uhap,*nloci, *nhap,train_hap1code,train_hap2code,train_wgt,happrob,
					*nobs_train,train_subj_rep,current_snp_set, count-1,min_freq, mustart,
					nngcov, train_ngcov,
					maxit,*Minherit, *tol, verbose,file); 
           

			    if (count1==0) {
             	best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
             	copy_hmodel(temp,best_model);
             	min_phi =  (temp->nhaps == base_model->nhaps) ? 0 : (temp -> dev - base_model->dev )/ (temp->nhaps - base_model->nhaps);  
             } else {
             temp_phi= (temp->nhaps == base_model->nhaps) ? 0 : (temp->dev-base_model->dev )/ (temp->nhaps - base_model->nhaps);  
            if (temp_phi < min_phi) {
            	Free_hmodel(best_model);
            	best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
            	copy_hmodel(temp,best_model);
            	min_phi=temp_phi;    
            }
            }
              count1++;
              Free_hmodel(temp);
	  }

          Free_hmodel(base_model);
          base_model = new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
          copy_hmodel(best_model,base_model);
         
     
          if (*verbose==1)  print_hmodel(best_model,file);
          if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 
          free(snp_set);
          snp_set=long_vec(count-1);
          for (k=0;k<(count-1);k++)  snp_set[k] = best_model->snp_set[k]; 
          free(current_snp_set);
          if (*verbose==1)   print_vector_long(snp_set,count-1,file);
          if (i!=*nfold) test_nsub = ncross; else test_nsub = ncross + nextra;              
	  			dev_mat[i][2* (*maxsnps)-count+1]=cross_val(
					best_model,test_size[i], test_csctl,test_nsub, test_subj_rep, uhap,*nloci,
					*nhap,test_hap1code,test_hap2code, happrob,
					nngcov, test_ngcov,
					*Minherit,verbose,phase,file); 
          Free_hmodel(best_model);     
          count--;

     }

// Clear memory 
             current_snp = &SNP_out;
             while (current_snp->next !=NULL) {                 
                         next_snp=current_snp -> next;  
                         current_snp ->next = next_snp->next;
                         free(next_snp);                     
             }
             if (*verbose==1) print_list_long(SNP_out,file);

	     free(train_hap1code);
             free(train_hap2code);
             free(train_csctl);

             free(train_wgt);

	     free(test_hap1code);   
             free(test_hap2code); 
	     free(test_csctl);
	     free(test_wgt);

             free(train_subj_rep);
             free(test_subj_rep);  
             free(snp_set);
             free(mustart);  
             Free_hmodel(base_model);

     }

// Sum deviance across folds
// The deviance matrix is nfold x maxsnps*2 
// The first half of the values are from growing and the last half are from pruning

     for (j=0;j<*maxsnps*2;j++) {    
       for (i=0; i<*nfold;i++) { 
           deviance[j] += dev_mat[i][j];
       }
     }     

     if (*verbose==1) print_matrix_double(dev_mat,*nfold, *maxsnps*2,file);
     
     free(nobs_train);
     free(cvlabel);
     free(cvlabel2);
     free(test_size);
     for (i=0;i<*nfold;i++) free(dev_mat[i]);
     free(dev_mat);     
     
	  fclose(file);

}  
     





/* this function is to find the best model with size m using all data */
void finalsubset(long *indx_subj,long *nsubj,long *nobs,long *subj_rep,long *csctl,long *nloci,long *nhap,long *hap1code,long *hap2code,long *uhap,double *happrob,double *wgt,
long *nngcov, double *ngcov_vec,
long *maxsnps, long *bestsize, long *output_snp_set, long *out_haplo_vec, double *out_hap_freq, double *tol, long *phase, double *varstore, double *coef, long *out_nhaps, long *Minherit, long *verbose){

  /*
    long *indx_subj,
    long *nsubj,
    long *nobs,
    long *subj_rep,
    long *csctl,
    long *nloci,
    long *nhap,
    long *hap1code,
    long *hap2code,
    long *uhap,
    double *happrob,
    double *wgt,
    long *nngcov, 
    double *ngcov,  
    long *maxsnps, 
    long *bestsize, 
    long *output_snp_set, 
    long *out_haplo_vec, 
    double *out_hap_freq, 
    double *tol, 
    long *phase, 
    double *varstore, 
    double *coef, 
    long *out_nhaps, 
    long *verbose
   */

 long count, count1,i,j, k, maxit=1000, *snp_set, *current_snp_set;
 double *mustart, min_phi, temp_phi, min_freq, **ngcov_mat;  
 hmodel *temp, *base_model, *best_model; 
 Node SNP_out;
 Link1 current_snp= &SNP_out, next_snp;   

 FILE *file; 
 file = fopen("debug_finalsubset.txt", "w"); 


if (*verbose==1 && *nngcov != 0) {
fprintf(file,"\n\nngcov_vec in final : \n");
print_vector_double(ngcov_vec,(*nobs * *nngcov),file);
fprintf(file,"\nngcov_mat in final :\n");
//print_matrix_double(ngcov_mat,*nobs,*nngcov,file);
 }

 
// put covariate data in matrix
if (*nngcov>1) {
ngcov_mat = double_vec_to_mat(ngcov_vec, *nobs, *nngcov);
      }


         snp_set = long_vec(*maxsnps);
         /* prepare the linked list to start the search */
         SNP_out.next=NULL;


         for (j=0;j<*nloci;j++) {
                  next_snp = malloc(sizeof(Link1));
                  next_snp->snp_pos = j+1;
                  next_snp->next = NULL;
                  current_snp -> next = next_snp;
                  current_snp=next_snp;
	 }

         if (*verbose==1)    print_list_long(SNP_out,file);
              
// Model search 
        
         mustart=double_vec(*nobs);
         for (j=0;j<*nobs;j++)  { mustart[j]= (wgt[j]*csctl[j]+0.5)/(wgt[j]+1);}

         min_freq= 2.0/(2.0* (*nsubj));
         if (*verbose==1)    fprintf(file,"\nthe freq cutoff = %8.9f",min_freq);

			//if best size is less than or equal to maxsnps
         if (*bestsize <= *maxsnps) {
            count=0;
            while (count < *bestsize) {

                current_snp = &SNP_out;
                count1=0;

             while (current_snp->next !=NULL) {
                 next_snp=current_snp -> next;
                 snp_set[count]= next_snp -> snp_pos;
                 if (*verbose==1)     fprintf(file,"\nthe current snp set\n");
                 if (*verbose==1)    print_vector_long(snp_set,count+1,file);
                 temp=hap_shrink_reg(
						csctl, uhap,*nloci, *nhap,hap1code,hap2code,wgt,happrob,*nobs,subj_rep,
						snp_set, count+1,min_freq,mustart,
						nngcov, ngcov_mat, 
						maxit,*Minherit, *tol, verbose,file); 
	         
                 if  (count==0) {
        
                    if (count1==0) {
                                     best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                     copy_hmodel(temp,best_model);
                            }
                            else if  (best_model->dev < temp -> dev)  {
                                     Free_hmodel(best_model);
                                     best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                     copy_hmodel(temp,best_model);
                            }
                                          
                 } else {
                            if (count1==0) {
                                    Free_hmodel(best_model);
                                    best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                    copy_hmodel(temp,best_model);
                                    min_phi =  (temp->nhaps == base_model->nhaps) ? 0 : (temp -> dev - base_model->dev )/ (temp->nhaps - base_model->nhaps);  
                            } else {
                                    temp_phi= (temp->nhaps == base_model->nhaps) ? 0 : (temp->dev-base_model->dev )/ (temp->nhaps - base_model->nhaps);  
                                    if (temp_phi > min_phi) {
                                          Free_hmodel(best_model);
                                          best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                          copy_hmodel(temp,best_model);
                                          min_phi=temp_phi;    
                                    }
                            }
                 }
                 if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 
	         if (*verbose==1)   print_hmodel(best_model,file);
                 current_snp = next_snp;
		 Free_hmodel(temp);                 
                 count1++;
	     }                  

             
             if (count==0) { 
                     base_model = new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
                     copy_hmodel(best_model,base_model);
             } else { 
                     Free_hmodel(base_model);
                     base_model = new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
                     copy_hmodel(best_model,base_model);
             }            
  
             if (*verbose==1)  print_hmodel(best_model,file);
             if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 

             snp_set[count]= best_model->snp_set[count];
             if (*verbose==1)   print_vector_long(snp_set,count+1,file);
 
             /*delete the selected SNPs from SNP_out */
             current_snp = &SNP_out;
             while (current_snp->next !=NULL) {                 
                      next_snp=current_snp -> next;  
                      if (next_snp -> snp_pos ==  snp_set[count])  {
                         current_snp ->next = next_snp->next;
                         free(next_snp);
                      } else {
                         current_snp = next_snp;
                      }
             }
             if (*verbose==1)     print_list_long(SNP_out,file);
             count++;
	 			}

		// if best size is not less than or equal to max snps
	 		} else {

            count=0;
            while (count < *maxsnps) {

                current_snp = &SNP_out;
                count1=0;

             while (current_snp->next !=NULL) {
                 next_snp=current_snp -> next;
                 snp_set[count]= next_snp -> snp_pos;
                 if (*verbose==1)     fprintf(file,"\nthe current snp set\n");
                 if (*verbose==1)    print_vector_long(snp_set,count+1,file);


                 temp=hap_shrink_reg(
						csctl, uhap,*nloci, *nhap,hap1code,hap2code,wgt,happrob,*nobs,subj_rep,snp_set, 
						count+1,min_freq, mustart,
						nngcov, ngcov_mat,
						maxit,*Minherit,*tol, verbose,file); 
	         
                 if  (count==0) {
                            if (count1==0) {
                                     best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                     copy_hmodel(temp,best_model);
                            }
                            else if  (best_model->dev < temp -> dev)  {
                                     Free_hmodel(best_model);
                                     best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                     copy_hmodel(temp,best_model);
                            }
                                          
                 } else {
                            if (count1==0) {
                                    Free_hmodel(best_model);
                                    best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                    copy_hmodel(temp,best_model);
                                    min_phi =  (temp->nhaps == base_model->nhaps) ? 0 : (temp -> dev - base_model->dev )/ (temp->nhaps - base_model->nhaps);  
                            } else {
                                    temp_phi= (temp->nhaps == base_model->nhaps) ? 0 : (temp->dev-base_model->dev )/ (temp->nhaps - base_model->nhaps);  
                                    if (temp_phi > min_phi) {
                                          Free_hmodel(best_model);
                                          best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                          copy_hmodel(temp,best_model);
                                          min_phi=temp_phi;    
                                    }
                            }
                 }
                 if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 
	         if (*verbose==1)   print_hmodel(best_model,file);
                 current_snp = next_snp;
		 Free_hmodel(temp);                 
                 count1++;
	     }                  

             
             if (count==0) { 
                     base_model = new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
                     copy_hmodel(best_model,base_model);
             } else { 
                     Free_hmodel(base_model);
                     base_model = new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
                     copy_hmodel(best_model,base_model);
             }            
  
             if (*verbose==1)  print_hmodel(best_model,file);
             if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 

             snp_set[count]= best_model->snp_set[count];
             if (*verbose==1)   print_vector_long(snp_set,count+1,file);
 
             /*delete the selected SNPs from SNP_out */
             current_snp = &SNP_out;
             while (current_snp->next !=NULL) {                 
                      next_snp=current_snp -> next;  
                      if (next_snp -> snp_pos ==  snp_set[count])  {
                         current_snp ->next = next_snp->next;
                         free(next_snp);
                      } else {
                         current_snp = next_snp;
                      }
             }
             if (*verbose==1)     print_list_long(SNP_out,file);
             count++;
	 }   


      
     // Pruning
         if (*verbose==1) fprintf(file,"\nstart pruning \n") ;
        while (count>(2*(*maxsnps)-*bestsize)) {

          current_snp_set = long_vec(count-1);
          count1=0;
          for (j=0;j<count;j++) {

                for (k=0;k<(count-1);k++)  {
                          if (k<j) current_snp_set[k] = snp_set[k];  
                          else current_snp_set[k] = snp_set[k+1];  
                }

		temp=hap_shrink_reg(
			csctl, uhap,*nloci, *nhap,hap1code,hap2code,wgt,happrob,*nobs,subj_rep,current_snp_set, 
			count-1,min_freq,mustart,
			nngcov, ngcov_mat,
			maxit,*Minherit,*tol, verbose,file); 
         

      if (count1==0) {
                                    Free_hmodel(best_model);
                                    best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                    copy_hmodel(temp,best_model);
                                    min_phi =  (temp->nhaps == base_model->nhaps) ? 0 : (temp -> dev - base_model->dev )/ (temp->nhaps - base_model->nhaps);  
                            } else {
                                    temp_phi= (temp->nhaps == base_model->nhaps) ? 0 : (temp->dev-base_model->dev )/ (temp->nhaps - base_model->nhaps);  
                                    if (temp_phi < min_phi) {
                                          Free_hmodel(best_model);
                                          best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
                                          copy_hmodel(temp,best_model);
                                          min_phi=temp_phi;    
                                    }
                            }
              count1++;
              Free_hmodel(temp);
	  }

          Free_hmodel(base_model);
          base_model = new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
          copy_hmodel(best_model,base_model);
          

          if (*verbose==1)  print_hmodel(best_model,file);
          if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 
          free(snp_set);
          snp_set=long_vec(count-1);
          for (k=0;k<(count-1);k++)  snp_set[k] = best_model->snp_set[k]; 
          free(current_snp_set);
          if (*verbose==1)   print_vector_long(snp_set,count-1,file);       
          count--;
     }

	 } //end loop for best size vs max snps


// snp set and variance estimation
     for (i=0;i<count;i++) output_snp_set[i]=snp_set[i]; 

     varest(best_model,nsubj,nobs,subj_rep,csctl,nloci,nhap,hap1code,hap2code,uhap,
			happrob,wgt,snp_set,count,min_freq,
			nngcov, ngcov_mat,
			*Minherit,verbose,file,phase,varstore);

// gather info for best model
     for (i=0;i<best_model->nhaps;i++) out_hap_freq[i] = best_model->hfreq[i];
     for (i=0;i<(best_model->nhaps + best_model->nngcov);i++) coef[i]=best_model->coef[i];
     for (i=0;i<(best_model->nsnps * best_model->nhaps);i++) out_haplo_vec[i]=best_model->haplo_vec[i];
     
     *out_nhaps=best_model->nhaps;
     /*     final_model=new_hmodel(best_model->nsnps,best_model->nhaps,best_model->nngcov,best_model->nrec);
	    copy_hmodel(best_model,final_model); */

// free memory
     Free_hmodel(base_model);
     Free_hmodel(best_model);
     free(snp_set);
     free(mustart);

     /*   if (*verbose==1) print_hmodel(final_model,file);*/
     current_snp = &SNP_out;
             while (current_snp->next !=NULL) {                 
                         next_snp=current_snp -> next;  
                         current_snp ->next = next_snp->next;
                         free(next_snp);                     
             }


    fclose(file);
}


/* this function is to compute the cv deviance for external data */

/*
void extern_cross_val(
      long  *nrecords,               
      long  *y,                      
      long  *nsubj,                  
      long  *subj_reps,              
      long  *haplotype_vec,          
      long  *nsnps,                  
      long  *nhaps,                  
      long  *hap1code,               
      long  *hap2code, 
      double *hfreq,
      long *verbose, 
      long *phase,
      double *dev)
{
  FILE *file;    
    file = fopen("debug2.txt", "w");
    if (*verbose==1) print_hmodel(final_model,file); 
    *dev=cross_val(final_model,*nrecords,y,*nsubj,subj_reps,haplotype_vec,*nsnps,*nhaps,hap1code,hap2code,hfreq,verbose,phase,file);
    Free_hmodel(final_model);
    fclose(file);
}

*/


/* This function computes the (sandwich) variance of final output model */

void varest(hmodel *best_model,long *nsubj,long *nobs,long *subj_rep,long *csctl,long *nloci,long *nhaps,long *hap1code,long *hap2code,long *uhap,double *happrob,double *wgt,long *snp_set, long nsnp_set, double cutoff,long *nngcov, double **ngcov, 
long Minherit, long *verbose,FILE *file, long *phase, double *varstore){


   long i,j,k,count, *out_nhaps, *out_hap1code, *out_hap2code, hapbase, **out_haplotype;
   double *out_hfreq, **amat, **bmat, *varmu, *escore; 
   double **x, temp_freq, *pred, *err;

   out_nhaps=(long *)malloc(sizeof(long));
   out_hap1code=long_vec(*nobs);
   out_hap2code=long_vec(*nobs);
   out_hfreq=double_vec(*nhaps);
   
// If a subset of snps are chosen, recode hapolotypes
   if (*nloci!=nsnp_set) out_haplotype = hap_shrink(*nobs,*nloci,*nhaps,uhap,happrob, hap1code,hap2code,snp_set,nsnp_set,out_nhaps,out_hap1code, out_hap2code, out_hfreq, verbose,file);
   else {
        *out_nhaps=*nhaps;
        for (i=0;i<*nobs;i++) {
                   out_hap1code[i]=hap1code[i];
                   out_hap2code[i]=hap2code[i];
        }
        for (i=0;i<*nhaps;i++) {
                   out_hfreq[i] = happrob[i];
        }
        out_haplotype=long_vec_to_mat(uhap,*nhaps,*nloci);
        
   }   

// Cluster low frequency haplos to high 
   count=0;
   for (i=0; i<*out_nhaps;i++) {if (out_hfreq[i]<cutoff) count ++ ; }
 
   if (count >0)  haplo_cluster(*nobs, out_hap1code, out_hap2code,out_haplotype, out_hfreq,cutoff,nsnp_set, out_nhaps,verbose, file);         
   temp_freq=out_hfreq[0];
 
// Find haplotype with the largest frequency 

   hapbase=1; 
   for (i=1;i<*out_nhaps;i++) { 
           if (temp_freq <= out_hfreq[i]) {
                  hapbase=i+1;
                  temp_freq=out_hfreq[i];
	    }
   }


// Design matrix
   //x=double_matrix(*nobs,*out_nhaps);  
   x=double_matrix(*nobs,(*out_nhaps + *nngcov));  
   pred=double_vec(*nobs);
    varmu=double_vec(*nobs);
   err=double_vec(*nobs);

// Put haplotypes in design matrix
   for (i=0; i< *nobs; i++) {       
       x[i][0]=1;
       for (j=1;j<(*out_nhaps+1);j++) {
           if (j < hapbase) {   
               x[i][j]=0;     
               if (Minherit==1) {
                  if (j==out_hap1code[i]) x[i][j] += 1;                               
                  if (j==out_hap2code[i]) x[i][j] += 1;
               } 
               if (Minherit==2) {
                  if (j==out_hap1code[i] || j==out_hap2code[i]  ) x[i][j] += 1;                               
               }
               if (Minherit==3) {
                  if (j==out_hap1code[i] && j==out_hap2code[i]  ) x[i][j] += 1;                               
               }

           }
  
           if (j > hapbase) {
               x[i][j-1]=0;
               if (Minherit==1) {
                  if (j==out_hap1code[i]) x[i][j-1] += 1;                               
                  if (j==out_hap2code[i]) x[i][j-1] += 1;
               } 
               if (Minherit==2) {
                  if (j==out_hap1code[i] || j==out_hap2code[i]  ) x[i][j-1] += 1;                               
               }
               if (Minherit==3) {
                  if (j==out_hap1code[i] && j==out_hap2code[i]  ) x[i][j-1] += 1;                               
               }

               
           } 
      
        }

   if (*verbose==1) {
     fprintf(file,"\n\n x in varest before ngcov added : \n");
     print_matrix_double(x,*nobs,(*out_nhaps + *nngcov), file);}

		 // Add non-genetic covariates to design matrix 
   	 if (*nngcov!=0)  {
        	for (j=(*out_nhaps); j<(*out_nhaps + *nngcov); j++) {
                x[i][j]=ngcov[i][(j-*out_nhaps)];
         }
    	 }
 
   if (*verbose==1) {
     fprintf(file,"\n\n out_nhaps and best_mod nhaps = %i %i \t\n", *out_nhaps, best_model->nhaps); 
     fprintf(file,"\n\n x in varest after ngcov added : \n");
     print_matrix_double(x,*nobs,(*out_nhaps + *nngcov), file);}

// Prediction
       pred[i]=0;
     //   for (j=0;j<best_model->nhaps;j++) pred[i] += x[i][j]*best_model->coef[j];
        for (j=0;j<(best_model->nhaps + best_model->nngcov);j++) pred[i] += x[i][j]*best_model->coef[j];
        pred[i]=exp(pred[i])/(1+exp(pred[i]));
        err[i]= csctl[i]-pred[i];
        varmu[i]=pred[i]*(1-pred[i]);
      
	} 
   
   if (*verbose==1) {
     fprintf(file,"\n\n after prediction : \n");
     fprintf(file,"\nx \n");
     print_matrix_double(x,*nobs,(*out_nhaps + *nngcov), file);
     fprintf(file,"\npred \n");
     print_vector_double(pred,*nobs, file);
     fprintf(file,"\nerr \n");
     print_vector_double(err,*nobs, file);
     fprintf(file,"\nvarmu nobs \n");
     print_vector_double(varmu,*nobs,file); 
     fprintf(file,"\nvarmu wgt \n");
     print_vector_double(varmu,*wgt,file); 
   }

// Create matrix a
   amat=double_matrix((best_model->nhaps + best_model->nngcov),(best_model->nhaps + best_model->nngcov));
     for (i=0;i<(best_model->nhaps+best_model->nngcov);i++) {
       for  (j=0;j<(best_model->nhaps+best_model->nngcov);j++) {
	 amat[i][j] = 0.0;
       }
     } 
   
   count=0;
   while (count<*nobs) {
     for (i=0;i<(best_model->nhaps+best_model->nngcov);i++) {
       for  (j=0;j<(best_model->nhaps+best_model->nngcov);j++) {
	 amat[i][j]  += x[count][i] * x[count][j] *varmu[count]*wgt[count];
       } 
     }
     count++;
   }

// Create matrix b 
   if (*phase==0) {
     count=0;   
     bmat=double_matrix((best_model->nhaps+best_model->nngcov),(best_model->nhaps+best_model->nngcov));
          for (i=0;i<(best_model->nhaps+best_model->nngcov);i++) {
           for  (j=0;j<(best_model->nhaps+best_model->nngcov);j++) {
	    bmat[i][j] = 0.0;
           }
          } 
     escore=double_vec(best_model->nhaps+best_model->nngcov);
     for (i=0;i<*nsubj;i++) {
       for  (j=0;j<(best_model->nhaps+best_model->nngcov);j++) escore[j]=0.0;
       for (j=0;j<subj_rep[i];j++) {
	 for (k=0;k<(best_model->nhaps+best_model->nngcov);k++) {
           escore[k]+= wgt[count]*x[count][k]*err[count];
         } 
         count++;           
       }
       for (k=0;k<(best_model->nhaps+best_model->nngcov);k++) {
         for  (j=0;j<(best_model->nhaps+best_model->nngcov);j++) {
	   bmat[k][j]  += escore[k]*escore[j];
         } 
       }
     }


   }     

// Create varstore : holds values from a and b matrices
   count=0;
   for (i=0;i<(best_model->nhaps+best_model->nngcov);i++) {
        for  (j=0;j<(best_model->nhaps+best_model->nngcov);j++) {
	    varstore[count]=amat[i][j];
            count++;
           }
   } 

   if (*phase==0) { 
        for (i=0;i<(best_model->nhaps+best_model->nngcov);i++) {
          for  (j=0;j<(best_model->nhaps+best_model->nngcov);j++) {
	    varstore[count]=bmat[i][j];
            count++;
           }
        } 
  
   }


// Clear memory
  free(out_hap1code);
  free(out_hap2code);
  free(out_hfreq);
  for (i=0;i<*nobs;i++) free(x[i]);
  free(x);
  for (i=0;i<*out_nhaps;i++) free(out_haplotype[i]);
  free(out_haplotype);  
  free(out_nhaps);
  free(pred);
  free(err);
  free(varmu);
  if (*phase==0) free(escore);
  for (i=0;i<(best_model->nhaps+best_model->nngcov);i++) {
    free(amat[i]);
    if (*phase==0) free(bmat[i]);
  }
  free(amat);
   if (*phase==0) free(bmat);

   return;
   
}

     

/* this function is to compute the estimated probability of each pair given a set of starting hap probabilities*/

void shrink_phase_infer(
     long nsnp,               /* the number of SNPs */
     long *nrecord,            /* length of hap1code */   
     long *nhap,              /* number of haplotypes */
     long *hap1code,          /* the vector of codes for haplotype 1 */
     long *hap2code,         /* the vector of codes for haplotype 2 */
     long *uhap,             /* vector of unique haplotype */
     double *happrob,        /* the vector of haplotype probabilities */
     double *wgt,            /* the posterior probabilities of each pair of haplotypes */   
     long nsubj,              /* the number of subjects */
     long *subj_rep,         /* the vector of number of hap pair, length nsubj */
     double tol, 
     long *verbose,
     long ntestdata, 
     long *test_hap1code,
     long *test_hap2code, 
     FILE *file)           

{
 
    double *uhapwgt, wgtsum, newwgt, *sub_wgt, *sub_wgtsum, maxdiff, *diff, *happrob_new;
    long maxit,count1,count2, count, count3, *hap_del,*uhap_new, i,j,temp, *hap11code, *hap22code, *newcode;   
 
    
    uhapwgt=double_vec(*nhap);
    sub_wgt=double_vec(*nrecord);
    sub_wgtsum=double_vec(nsubj);
    happrob_new=double_vec(*nhap);
    uhap_new=long_vec(*nhap * nsnp);

    hap_del=long_vec(*nhap);
    diff = double_vec(*nhap);
    hap11code=long_vec(*nrecord);
    hap22code=long_vec(*nrecord);


    // EM algorithm to get the estimated hap prob 
    // maximum iterations set at 100 

    maxdiff=1;    
    if (*verbose==1) {
     print_vector_double(happrob,*nhap,file);
     print_vector_double(wgt,*nrecord,file);
    }

    for (i=0;i<*nrecord;i++) sub_wgt[i]=wgt[i];
    maxit=0;
    while ((maxdiff > tol) & (maxit<100)) {
         maxit++;
	 		wgtsum=0;
	 		for (i=0;i<*nhap;i++)  uhapwgt[i] =0;
         for (i=0;i<*nrecord;i++) {
             uhapwgt[hap1code[i]-1] += sub_wgt[i];
             uhapwgt[hap2code[i]-1] += sub_wgt[i];
         }
  
         if (*verbose==1) print_vector_double(uhapwgt,*nhap,file);
         
         for (j=0;j<*nhap;j++)  wgtsum += uhapwgt[j] ;
         for (j=0;j<*nhap;j++)  happrob_new[j]=uhapwgt[j]/wgtsum ;    

         if (*verbose==1) print_vector_double(happrob_new,*nhap,file);

         count1=0;
         for (i=0;i<nsubj;i++)  count1 += subj_rep[i];
         
         count1 = 0;
         count2 = 0;
         for (i=0; i<nsubj;i++) {      
             sub_wgtsum[count2]=0;
             for (j=0;j<subj_rep[i];j++) {
                newwgt = happrob_new[hap1code[count1]-1] *  happrob_new[hap2code[count1]-1];
                newwgt = hap1code[count1]!=hap2code[count1] ? (2*newwgt) : newwgt; 
                sub_wgtsum[count2] += newwgt;
                sub_wgt[count1] = newwgt;
                count1 ++;
             }       
             count2 ++;
         }
      
         count1 = 0;
         count2 = 0;
         for (i=0;i<nsubj;i++) {     
             for (j=0;j<subj_rep[i];j++) {
		 				sub_wgt[count1] =  sub_wgt[count1]/sub_wgtsum[count2];
                 	count1 ++;
             }       
             count2 ++;
         }

         for (i=0;i<*nhap;i++)      {
             diff[i]= fabs(happrob_new[i] - happrob[i]);
             if (i==0) maxdiff=diff[i]; else maxdiff= maxdiff>diff[i] ? maxdiff : diff[i]; 
             happrob[i] = happrob_new[i];
         }
			if (*verbose==1) fprintf(file,"\nthe maximal difference=%8.5f\n", maxdiff);       
    }     
    if (*verbose==1) print_vector_double(happrob,*nhap,file);
    if (maxit==100) fprintf(file,"\nthe em algorithm to phase the haplotype did not converge\n");
 

	// Delete haplos with frequency 0 
    count1 = 0;      
    for (i=0;i<*nhap;i++)      {
          if (happrob[i]==0.0) {
                  count1 ++;
                  hap_del[i]=0;
	  } else hap_del[i]=1;  
    }
     
    if (*verbose==1){   
      print_vector_long(hap_del,*nhap,file);
      print_vector_long(hap1code,*nrecord,file);
      print_vector_long(hap2code,*nrecord,file);  
      print_vector_double(sub_wgt,*nrecord,file); 
    }

    if (count1>0) {            
             newcode=long_vec(*nhap);
             for (i=0;i<*nhap;i++) {
                      newcode[i]=0; 
                      for (j=0;j<(i+1);j++) newcode[i]+=hap_del[j];
             }
             print_vector_long(newcode,*nhap,file);

             count2=0;
             count3=0;
             for (i=0;i<nsubj;i++) {                
	       for (j=0;j<subj_rep[i];j++) {                                   
                  count=0;
                  if (hap_del[hap1code[count3]-1] ==1 && hap_del[hap2code[count3]-1] ==1 ) {
                         temp = hap1code[count3];
                         hap11code[count2] =newcode[temp-1];
                         temp = hap2code[count3];   
                         hap22code[count2] =newcode[temp-1];                                      
                         wgt[count2] = sub_wgt[count3];
                         count2++;                         
                  } else count++;
                  count3++;
                         
               }
               subj_rep[i] = subj_rep[i] -count;

             } 

             for (i=0;i<ntestdata;i++) {
                      test_hap1code[i]=newcode[test_hap1code[i]-1];
                      test_hap2code[i]=newcode[test_hap2code[i]-1];
             }
             if (*verbose==1) {
             print_vector_long(hap11code,count2,file);
             print_vector_long(hap22code,count2,file); 
             print_vector_double(wgt,count2,file); 
             
             print_vector_long(test_hap1code,ntestdata,file);
             print_vector_long(test_hap2code,ntestdata,file); }
             *nrecord=count2;

             for (i=0;i<count2;i++) {
                      hap1code[i]=hap11code[i];
                      hap2code[i]=hap22code[i];
             }

             count2 = 0;
             for (i=0;i<*nhap;i++) {
                  if (happrob[i]!=0) {
                         happrob_new[count2] = happrob[i]; 
                         count2 ++; 
				  		}   
             }                          
	   

             for (i=0;i<*nhap;i++) {
                  if (i<count2) happrob[i] = happrob_new[i];  else happrob[i] =0;                	     
             }
             count=0;
             count2=0;
             for (i=0;i<*nhap;i++) {
                    for (j=0;j<nsnp;j++) {  
	              			if (hap_del[i]==1) { 		    
                 				uhap_new[count2]=uhap[count];
                 				count2++;
				 	  			}
		       			count++;  
                    }
             }
	
     			if (*verbose==1){
	  			    print_vector_long(uhap,count,file);
     			    print_vector_long(uhap_new,count2,file); }
  
 
             for (i=0;i<*nhap* nsnp;i++) {
                     if (i <count2) uhap[i]=uhap_new[i] ; else uhap[i]=0;
			 	 }
            
	     		*nhap=*nhap-count1;
             free(newcode);
             
     }     

     free(uhapwgt);
     free(sub_wgt);
     free(sub_wgtsum);
     free(happrob_new);
     free(uhap_new);
     free(hap_del);
     free(diff);
     free(hap11code);
     free(hap22code);   
  
     if (*verbose==1) fprintf(file,"\n the new number of haplotypes=%i\n", *nhap);
}







 /* This function is for changing the haplotype coding, given the full haplotype coding and a set of SNPs being used */


long **hap_shrink(
		long nrecords,           /* the number of records */ 
		long nsnps,              /* the total number of SNPs */
		long nhaps,              /* the number of haplotypes */ 
		long *fhaplotype_vec,    /* the vector of haplotype coding for all SNPs  */
		double *hfreq,           /* the vector of haplotype frequencies */
		long *hap1code,          /* the hap1 coding using full set */
		long *hap2code,          /* the hap2 coding using full set */
		long *snp_set,           /* the set of SNPs considered, each element is an integer index the SNP position */
		long nsnp_set,           /* the number of SNPs in the set */
		long *out_nhaps,         /* the number of haplotypes resulted */
		long *out_hap1code,      /* the vector of output hap1code */
		long *out_hap2code,      /* the vector of output hap2code */
		double *out_hfreq,        /* the vector of output haplotype frequencies */   
		long *verbose,
		FILE *file)
    {
      long i,j, k, *index,*index1, count, *newcode, *newcode1, *order;
      long **out_haplotype, **fhaplotype, **haplotype1, *thaplotype1;

 
      fhaplotype = long_vec_to_mat(fhaplotype_vec, nhaps, nsnps);
      newcode=long_vec(nhaps); 
      order=long_vec(nhaps);
      index=long_vec(nhaps);
      for (i=0;i<nhaps;i++) order[i]=i+1;
      newcode1=long_vec(nhaps); 
      thaplotype1=long_vec(nhaps);


      if (*verbose==1) {
			fprintf(file,"\nnhaps and nsnps = %i %i \t\n", nhaps, nsnps);
		//	fprintf(file,"\nfhaplotype matrix : \n");
	//		print_matrix_long(fhaplotype,nhaps,nsnps,file);
		}

     if (*verbose==1)  fprintf(file,"\nthe number of SNPs is %i",nsnp_set);
      if (nsnp_set >1) {

                haplotype1= long_matrix(nhaps,nsnp_set);
     		for (i=0;i<nhaps;i++) index[i]=0;
	        for (i=0;i<nhaps;i++) {
	           for (j=0;j<nsnp_set;j++) {
                     haplotype1[i][j] = fhaplotype[i][snp_set[j]-1];
                     index[i] += pow(10,nsnp_set-j-1)* haplotype1[i][j];
                   }
	      	}
		if (*verbose==1) {
                 print_matrix_long(haplotype1,nhaps,nsnp_set,file);
                 print_vector_long(index,nhaps,file);
                }
                index1=insertionSort(index,nhaps,order);
                if (*verbose==1){
                print_vector_long(index1,nhaps,file);
                print_vector_long(order,nhaps,file);
                }
                newcode[0]=1;
                count=1;
                for (i=1;i<nhaps;i++) {
                    if (index1[i] != index1[i-1]) count ++;                     
                    newcode[i] = count;
                }
                for (i=0;i<nhaps;i++) newcode1[order[i]-1]=newcode[i]; 
                if (*verbose==1) {
                print_vector_long(newcode,nhaps,file);
                print_vector_long(newcode1,nhaps,file); }
                *out_nhaps=0;

                for (i=0;i<nhaps;i++)  *out_nhaps = max_long(*out_nhaps,newcode[i]);

	        out_haplotype= long_matrix(*out_nhaps,nsnp_set);

	        count=0;
                for (i=0;i<*out_nhaps;i++) {
       	           for (j=0;j<nhaps;j++) {
	             if (newcode1[j]==(i+1)) {
                       for (k=0;k<nsnp_set;k++)  out_haplotype[i][k] = haplotype1[j][k];
                       break;
		       
		      }
		   }
                }
              if (*verbose==1)   print_matrix_long(out_haplotype,*out_nhaps,nsnp_set,file);		  

                for (i=0;i<*out_nhaps;i++) { 
		  out_hfreq[i]=0;
                  for (j=0;j<nhaps;j++) {
                      if (newcode1[j]==i+1) out_hfreq[i]+=hfreq[j];
                  }
                }

                for (i=0;i<nrecords;i++) {
                    out_hap1code[i] = newcode1[hap1code[i]-1];
                    out_hap2code[i] = newcode1[hap2code[i]-1];
                }
                for (i=0;i<nhaps;i++) free(haplotype1[i]);
                free(haplotype1);  
                free(index1); 
		
         
      }
      else {
	   if (*verbose==1)  print_vector_long(snp_set,nsnp_set,file);

	    for (i=0;i<nhaps;i++)  thaplotype1[i] = fhaplotype[i][snp_set[0]-1];                              			      
            for (i=0;i<nhaps;i++) {
                  if (thaplotype1[i]==1)    newcode[i] = 2; else newcode[i] = 1;
            }

          if (*verbose==1)   print_vector_long(newcode,nhaps,file);
            *out_nhaps = 2;

           if (*verbose==1)  fprintf(file,"\nthe new number of haplotype=%i",*out_nhaps);
	    out_haplotype= long_matrix(*out_nhaps,nsnp_set);

            for (i=0;i<*out_nhaps;i++)      out_haplotype[i][0]= i;

	    if (*verbose==1)   {
                      print_matrix_long(out_haplotype,*out_nhaps,nsnp_set,file);
		      print_vector_double(hfreq,nhaps,file);
            }
            for (i=0;i<*out_nhaps;i++) { 
		out_hfreq[i]=0.0;
                for (j=0;j<nhaps;j++) {
                      if (newcode[j]==i+1) out_hfreq[i]+=hfreq[j];
                }
            }

           if (*verbose==1)  print_vector_double(out_hfreq,2,file);

            for (i=0;i<nrecords;i++) {
                    out_hap1code[i] = newcode[hap1code[i]-1];
                    out_hap2code[i] = newcode[hap2code[i]-1];
            }
 
	  }  

        free(newcode);
        free(newcode1);
        free(order);
        free(thaplotype1);
        free(index); 
       
 
        for (i=0;i<nhaps;i++) free(fhaplotype[i]);
        free(fhaplotype);
        

	if (*verbose==1){ 
             print_vector_long(out_hap1code,nrecords,file);
	     print_vector_long(hap1code,nrecords,file); 
        }
      return out_haplotype;
}	  
            




// Pair wise distance of two haplotypes 

double haplo_dist(
         long *haplotype1,     
         long *haplotype2,
         double *locus_freq,
         long nlocus)
{
    long i;
    double d;
    d=0;
    for (i=0;i<nlocus;i++){
      if ((haplotype1[i]==1) & (haplotype2[i]==1)) d = d + (1-locus_freq[i]);
      if ((haplotype1[i]==0) & (haplotype2[i]==0)) d = d + locus_freq[i];
    }
    return(1-d/nlocus);
}    




// This function clusters the rare haplotypes to common ones */

void haplo_cluster(
       long nrecords,                      /* the number of records */ 
       long *hap1code,                     /* the haplotype 1 coding */
       long *hap2code,                     /* the haplotype 2 coding */
       long **haplotype,                /* the vector of input haplotype */
       double *freq,                       /* the vector of frequencies */
       double min_freq,                    /* the freq cut-off */
       long nlocus,                        /* the number of SNPs */
       long *nhap,                         /* the number of haplotypes outputed */
       long *verbose,
       FILE *file)
{   
        long i, j, k, index1, min_hap, index2, count1, count2, count, *haplotype1, *haplotype2, *map_set;
        double locus_freq[nlocus], min_dist, *cur_dist, *out_freq;
	 	  long rare_set[*nhap], common_set[*nhap], codebook[*nhap], temp_code[*nhap], **out_haplotype;

        if (*verbose==1) print_matrix_long(haplotype,*nhap,nlocus,file);

        for (i=0; i<nlocus;i++) {
	     	locus_freq[i]=0;
        	for (j=0;j<*nhap;j++) {
        		if (haplotype[j][i]==1) locus_freq[i] += freq[j];
        	}
        }
       if (*verbose==1)  print_vector_double(locus_freq,nlocus,file);
    
        count1=0;
		  count2=0;
        for (i=0; i<*nhap;i++) {
     	  if (freq[i]<min_freq) {
                 	rare_set[count1]=i;
		 				codebook[i]=9999;
                 	count1 ++ ;
          }
          else {
             	common_set[count2]=i;
	     			codebook[i]=i;
             	count2 ++ ;
          }
        }
 
       if (*verbose==1)  print_vector_long(rare_set,count1,file);
       if (*verbose==1)  print_vector_long(common_set,count2,file);
			
        	map_set = long_vec(count1);
			out_freq= double_vec(count2);
        	cur_dist= double_vec(count2);
			haplotype1 = long_vec(nlocus);
			haplotype2 = long_vec(nlocus);
        

        for (i=0;i<count1;i++) {
                index1 = rare_set[i];
                for (j=0;j<nlocus;j++) haplotype1[j] = haplotype[index1][j];
                if (*verbose==1) print_vector_long(haplotype1,nlocus,file);
                min_dist=1;
					 min_hap=0;
                for (j=0;j<count2;j++) {
                    index2 = common_set[j]; 
                    for (k=0;k<nlocus;k++) haplotype2[k] = haplotype[index2][k];
                if (*verbose==1)     print_vector_long(haplotype2,nlocus,file);
                    cur_dist[j]=haplo_dist(haplotype1,haplotype2,locus_freq,nlocus);
		    if (min_dist>cur_dist[j]) min_hap=j;
                    min_dist=min_double(min_dist,cur_dist[j]);

                }
                if (*verbose==1) print_vector_double(cur_dist,count2,file);
                map_set[i] = common_set[min_hap];
	}

     if (*verbose==1)    print_vector_long(map_set,count1,file);
        
        for (i=0;i<*nhap;i++) {
		  temp_code[i]=i;
		  for (j=0; j< count1; j++) {
                       if (i==rare_set[j]) {
                         temp_code[i]=map_set[j];
                         break;
                       }
                  }
	}                    
       if (*verbose==1)  print_vector_long(temp_code,*nhap,file);
	for (i=0;i<count2;i++) {
            for (j=0;j<*nhap;j++) {
    			if (temp_code[j]==common_set[i]) out_freq[i] += freq[j];
			}	
        }
	
        for (i=0;i<*nhap;i++) {
			 if (i<count2) freq[i]=out_freq[i]; else freq[i]=0;
        }

       if (*verbose==1)  print_vector_double(freq,count2,file);
	for (i=0;i<*nhap;i++) {
	       count=0;
	       for (j=0;j<*nhap;j++) {
                  if (temp_code[i]>=codebook[j]) count++;
               }
	       temp_code[i]=count;
        }

        for (i=0;i<nrecords;i++) {
               hap1code[i]=temp_code[hap1code[i]-1];
	       hap2code[i]=temp_code[hap2code[i]-1];
        }
	

        out_haplotype=long_matrix(count2,nlocus);
        for (i=0;i<count2;i++) {
               for (j=0;j<nlocus;j++) out_haplotype[i][j] = haplotype[common_set[i]][j];
        }

        for (i=0;i<*nhap;i++) {
             if (i<count2) {
                     for (j=0;j<nlocus;j++)  haplotype[i][j] = out_haplotype[i][j];
             } else {
                    free(haplotype[i]);
             }
        }

      if (*verbose==1)   print_matrix_long(haplotype,count2,nlocus,file);
      if (*verbose==1)   print_vector_long(temp_code,*nhap,file);

        *nhap=count2;

        free(out_freq);
        free(map_set);
        free(cur_dist);
        free(haplotype1);
        free(haplotype2);

        for (i=0;i<count2;i++) free(out_haplotype[i]);
        free(out_haplotype);


}        









// This function computes the MLE in the logistic regression via iwls algorithm */
// iwls calls dqrls

void iwls_bin(
		long n,                      /*  the number of observations */
		double **x,                /*  the design matrix   */
		long ncov,                   /*  the number of columns in design matrix */
		long *y,                    /*  the binary outcome  */
		double *weights,            /*  the case weights assigned to each observation */
		double *mustart,            /*  the starting values of mu */
		long maxit,                  /*  the maximum number of iteration */
		double tol,                  /*  the precision aimed in iteration */
		double *coef,               /*  the output regression coefficients */
		double *deviance,            /*  the deviance at convergence */
		long *conv,                  /*  the indicator of algorithm convergence */
		long *verbose,
		FILE *file)
{ 


   long ny=1, rank=1, pivot[ncov], iter, i, j,count;
   double *eta, *mu_eta_val,lhs[n][ncov], *lhs2, *rhs, *resid, *effect, *work, *qraux, current_coef[ncov], current_dev, *mu, *varmu, *z, *w;

   eta=double_vec(n);
   mu_eta_val= double_vec(n);

   mu=double_vec(n);
   varmu=double_vec(n);
   z=double_vec(n);
   w=double_vec(n);
   lhs2=double_vec(ncov*n);
   rhs=double_vec(n);
   resid=double_vec(n);
   effect=double_vec(n);
   work=double_vec(2*ncov);
   qraux=double_vec(ncov);
   
  
   /*  print_vector_double(mustart,n);   fprintf(file,"\n");
   print_matrix_double(x,n,ncov,file);   fprintf(file,"\n");
   print_vector_long(y,n,file);   fprintf(file,"\n");
   print_vector_double(weights,n,file); */

   for (i=0;i<n;i++) eta[i] = log(mustart[i]/(1-mustart[i])); 
   
   for (iter=0;iter<maxit;iter++) {

                 if (iter==0) {
                         *deviance=0;
                         for (i=0;i<n;i++)  { 
                                mu[i]=mustart[i];                             
		                *deviance += y[i]==1?(weights[i]*log(mu[i])) : (weights[i]*log(1-mu[i]));
                         }

                 }

                 for (i=0;i<n;i++)  {
                        varmu[i] = mu[i]*(1-mu[i]); 
                        mu_eta_val[i] = mu[i]*(1-mu[i]);  
                        z[i] = eta[i] + (y[i]-mu[i])/mu_eta_val[i];
                        w[i] = sqrt(weights[i]*mu_eta_val[i]);  
                        for (j=0;j<ncov;j++) lhs[i][j] = x[i][j]*w[i];
                        rhs[i]=w[i]*z[i];
                 }


		 /*              print_vector_double(varmu,n);
                 print_vector_double(z,n);
                 print_vector_double(w,n); */
                 for (i=0; i<ncov; i++) {
                            pivot[i]=i+1; 
                            coef[i]=0;
                            qraux[i]=0;
                 }
                 for (i=0; i<n; i++) {
                            resid[i]=0;
                            effect[i]=0;
                 }
                 for (i=0;i<2*ncov;i++) work[i]=0;
 
		 /*  print_vector_double(work,2*ncov);*/
                 count=0;
                 for (i=0;i<n;i++) {
                        for (j=0;j<ncov;j++)  {
                               lhs2[count]=lhs[i][j];
                               count++;
                        }

                 }

		 /*
                 print_matrix_double(lhs2,ncov,n);                
                 print_vector_double(w,n); 
                 fprintf(file,"\n");
                 print_matrix_double(lhs,n,ncov);
                 fprintf(file,"\n");
                 for (i=0;i<n;i++) {
                        fprintf(file,"%8.5f\n",rhs[i]);
                 }
                 fprintf(file,"\n");
                 if (*verbose==1)    print_vector_double(weights,n); */
 
                 dqr_(lhs2, &n, &ncov, rhs, &ny, &tol, coef, resid, effect, &rank, pivot, qraux, work);      

              if (*verbose==1)    fprintf(file,"\npass dqrls \n");

		 /*               print_matrix_double(lhs2,ncov,n);   
                 for (i=0;i<ncov;i++) free(lhs2[i]);
                 free(lhs2); */
                 

                 /*
                 print_vector_double(coef,ncov);
                 print_vector_long(pivot,ncov); */
                 for (i=0;i<ncov;i++) {                    
                     current_coef[i]=coef[i];
                 }

	      if (*verbose==1)    print_vector_double(current_coef,ncov,file); 
                 current_dev=0.0;
                 for (i=0;i<n;i++)  { 
                   eta[i]=0;
		   		  for (j=0;j<ncov;j++)  eta[i] += x[i][j]*current_coef[j];                   
                   mu[i] = exp(eta[i])/(1+exp(eta[i]));
		   			 current_dev += y[i]==1?(weights[i]*log(mu[i])) : (weights[i]*log(1-mu[i]));
                 }
		 /*   print_vector_double(eta,n); */
                 if (fabs(*deviance-current_dev)/(0.1+fabs(current_dev)) < tol) {
                        *conv = 1;
                        *deviance = current_dev;
                      if (*verbose==1)    fprintf(file,"\nIWLS algorithm converged\n");  
                        break;
                 } else      *deviance = current_dev;
              if (*verbose==1)    fprintf(file,"\nthe current deviance is%8.5f",*deviance);
                 iter++;
   }               

   free(eta);
   free(mu_eta_val);
   free(mu);
   free(varmu);
   free(z);
   free(w);
   free(lhs2);
   free(rhs);
   free(resid);
   free(work);
   free(effect);
   free(qraux);
                
}




// This function computes the deviance of a particular model for a subset of SNPs 
hmodel *hap_shrink_reg( 
		  long *y,             	/* the vector of outcomes */ 
		  long  *haplotype_vec,	/* the vector of unique haplotype coding */
        long  nsnps,         	/* the number of SNPs to start with */ 
        long  nhaps,         	/* the number of haplotypes to start with */
        long  *hap1code,     	/* the haplotype coding for each observations */
        long  *hap2code,     	                     
        double *weight,      	/* the posterior probabilities for haplotype pair */
        double *hfreq,       	/* the vector of frequencies for initial haplotypes */ 
        long nrecords,       	/* the number of records */
        long *nreps,         	/* the vector of records for each subject */
        long *snp_set,       	/* the subset of SNPs under investigation */
        long nsnp_set,       	/* the number of snps considered in the subset */
        double cutoff,       	/* the frequency cut-off of haplotypes resulted */
        double *mustart,     	/* the vector of starting values */
        long *nngcov,     	  	/* the number of non-genetic covariates */
        double **ngcov,      	/* the matrix of non-genetic covariates */
        long maxit,          	/* the number of maximum iterations in the irls algorithm */
        long Minherit,       	/* mode of inheritance */
        double tol,
        long *verbose, 
        FILE *file)                  
{

   long ncov, i,j,c, *out_nhaps, *out_hap1code, *out_hap2code, count, hapbase,conv=0, **out_haplotype, *out_haplotype_vec;
   double *out_hfreq, *deviance; 
   double **x, temp_freq, *coef;
   hmodel *result_model;


   deviance=(double *)malloc(sizeof(double));
   out_nhaps=(long *)malloc(sizeof(long));
   out_hap1code=long_vec(nrecords);
   out_hap2code=long_vec(nrecords);
   out_hfreq=double_vec(nhaps);
  
   if (*verbose==1 && *nngcov != 0) {    
      fprintf(file,"\n nngcov inside hap_shrink_reg %i \n", *nngcov);
      fprintf(file,"\n covariates inside hap_shrink_reg \n");
      print_matrix_double(ngcov,nrecords,*nngcov, file);
      }

// Transfer haplotype data to "out" containers
   // If using a subset of snps :
   if (nsnps!=nsnp_set) {
			out_haplotype = hap_shrink(
			nrecords,nsnps,nhaps,haplotype_vec,hfreq, hap1code,hap2code,snp_set,nsnp_set,
			out_nhaps,out_hap1code, out_hap2code, out_hfreq, verbose,file); }
	// If using all snps :
   else {
        *out_nhaps=nhaps;
        for (i=0;i<nrecords;i++) {
                   out_hap1code[i]=hap1code[i];
                   out_hap2code[i]=hap2code[i];
        }
        for (i=0;i<nhaps;i++) {
                   out_hfreq[i] = hfreq[i];
        }
        out_haplotype=long_vec_to_mat(haplotype_vec,nhaps,nsnps);
   }   
   

	if (*verbose==1) {   
   fprintf(file,"\ncheckpoint2 in hapreg\n");
   fprintf(file,"\nthe freq cutoff = %8.5f",cutoff);
   print_vector_double(out_hfreq,*out_nhaps, file);
   }

	// Count number of haplos below frequency cutoff
   count=0;
   for (i=0; i<*out_nhaps;i++) {if (out_hfreq[i]<cutoff) count ++ ; }
   if (*verbose==1) fprintf(file,"\ncount=%i\n",count); 
	
	// Cluster low frequency haplos with high frequency haplos
   if (count >0)  {
		haplo_cluster(
			nrecords, out_hap1code, out_hap2code,out_haplotype, out_hfreq,
			cutoff,nsnp_set, out_nhaps,verbose, file);    }     
   temp_freq=out_hfreq[0];
 
	if (*verbose==1)   fprintf(file,"\ncheckpoint3 in hapreg\n");
   

	// Find haplotype with the largest frequency 
   hapbase=1; 
   for (i=1;i<*out_nhaps;i++) { 
           if (temp_freq <= out_hfreq[i]) {
                  hapbase=i+1;
                  temp_freq=out_hfreq[i];
	    }
   }

// Create design matrix
	// Intercept term is added here; we omit haplo with largest frequency
	// design matrix =  [ observations x (intercept, haplos, non-genetic covariates) ]
	// nngcov represents the number of non-genetic covariates
	// out_nhaps and nngcov represent the number of haplotypes and non-genetic covariates, respectively 

	if (*verbose==1) {  
			fprintf(file,"\nnrecords = %i \n",nrecords);
			fprintf(file,"\nout_nhaps = %i \n",*out_nhaps);
			fprintf(file,"\nnngcov = %i \n",*nngcov);
			fprintf(file,"\nout_nhaps + nngcov = %i \n",*out_nhaps + *nngcov);}

   //x=double_matrix(nrecords,*out_nhaps);  
   x=double_matrix(nrecords,*out_nhaps + *nngcov); 

   if(*verbose==1) fprintf(file,"\n after x allocation \n");

// Compute haplo codes based on inheritance and put in design matrix 
   for (i=0; i< nrecords; i++) {       
       x[i][0]=1;
       for (j=1;j<(*out_nhaps+1);j++) {
           if (j < hapbase) {   
               x[i][j]=0;   
               if (Minherit==1) {
                  if (j==out_hap1code[i]) x[i][j] += 1;                               
                  if (j==out_hap2code[i]) x[i][j] += 1;
               } 
               if (Minherit==2) {
                  if (j==out_hap1code[i] || j==out_hap2code[i]  ) x[i][j] += 1;                               
               }
               if (Minherit==3) {
                  if (j==out_hap1code[i] && j==out_hap2code[i]  ) x[i][j] += 1;                               
               }    
               
           }
  
           if (j > hapbase) {
               x[i][j-1]=0;  
               if (Minherit==1) {
                  if (j==out_hap1code[i]) x[i][j-1] += 1;                               
                  if (j==out_hap2code[i]) x[i][j-1] += 1;
               } 
               if (Minherit==2) {
                  if (j==out_hap1code[i] || j==out_hap2code[i]  ) x[i][j-1] += 1;                               
               }
               if (Minherit==3) {
                  if (j==out_hap1code[i] && j==out_hap2code[i]  ) x[i][j-1] += 1;                               
               }    
             
           } 
      
        }
   }    
  
   if (*verbose==1) {
   fprintf(file,"\nnngcov = %i \n", *nngcov);
   fprintf(file,"\nout_nhaps = %i \n", *out_nhaps);
   fprintf(file,"\ndesign matrix before covariates :\n");
   print_matrix_double(x,nrecords,(*out_nhaps+*nngcov),file); }


// Add non-genetic covariates to design matrix 
   if (*nngcov!=0)  {
    for (i=0; i< nrecords; i++) {       
        for (j=*out_nhaps;j<(*out_nhaps+*nngcov);j++) {
                x[i][j]=ngcov[i][j-*out_nhaps];
					//	if (*verbose==1) fprintf(file,"i and j = %i %i \t\n",i, j);
         }
    }
 	}   


   if (*verbose==1) {
   fprintf(file,"\nout_hap1code :\n");
   print_vector_long(out_hap1code,nrecords, file);
   fprintf(file,"\nout_hap2code :\n");
   print_vector_long(out_hap2code,nrecords, file);
   fprintf(file,"\ndesign matrix after covariates:\n");
   print_matrix_double(x,nrecords,(*out_nhaps+*nngcov),file); }

   // ncov is number of total covariates in the design matrix = haplos + non-genetic covariates
	//   ncov=*out_nhaps;  
   ncov=*out_nhaps + *nngcov;

   coef=double_vec(ncov);

	// call iwls for least squares estimates   	
   iwls_bin(nrecords,x,ncov,y,weight,mustart,maxit,tol,coef,deviance,&conv,verbose,file);


  out_haplotype_vec=long_mat_to_vec(out_haplotype, *out_nhaps, nsnp_set);
  

if (*verbose==1) {
  fprintf(file,"\npass iwls_bin\n");
  fprintf(file,"\nout_haplotype as matrix :\n");
  print_matrix_long(out_haplotype,*out_nhaps, nsnp_set,file);
  fprintf(file,"\nout_haplotype as vector :\n");
  print_vector_long(out_haplotype_vec,*out_nhaps*nsnp_set,file);
  }


  result_model=(hmodel *)Calloc(1, hmodel);
  
// Replace ncov with *out_nhaps, ncov now represents haplos + non-genetic covariates
 
   if (result_model) {
        //result_model -> nhaps = ncov;
        result_model -> nhaps = *out_nhaps;
        result_model -> nngcov = *nngcov;
        result_model -> nsnps = nsnp_set;
        result_model -> nrec = nrecords;
        result_model -> dev=*deviance;
        result_model -> hapbase=hapbase;
        
		result_model -> snp_set = (long *) Calloc(nsnp_set,long);
		if ( result_model -> snp_set )     {
	  	for (i=0;i<nsnp_set;i++) result_model -> snp_set[i] = snp_set[i]; }

      //result_model -> haplo_vec = (long *) Calloc(ncov*nsnp_set,long);
      result_model -> haplo_vec = (long *) Calloc(*out_nhaps*nsnp_set,long);
      if ( result_model -> haplo_vec )     {
	   //for (i=0;i<ncov*nsnp_set;i++) result_model -> haplo_vec[i] = out_haplotype_vec[i]; }
	   for (i=0;i<*out_nhaps*nsnp_set;i++) result_model -> haplo_vec[i] = out_haplotype_vec[i]; }

		if(*nngcov!=0) {
	   result_model -> ngcov_vec = (double *) Calloc(*nngcov*nrecords,double);
      	if ( result_model -> ngcov_vec )     {
		   c=0;	
			for (i=0;i<*nngcov;i++) {
				for (j=0;j<nrecords;j++) {
					result_model -> ngcov_vec[c] = ngcov[j][i];	
					c++;
				}
			}
			}
			}	

		//result_model -> hfreq = (double *) Calloc(ncov,double);
      result_model -> hfreq = (double *) Calloc(*out_nhaps,double);
		if ( result_model -> hfreq )     {
	   //for (i=0;i<ncov;i++)  result_model->hfreq[i]=out_hfreq[i]; }
	   for (i=0;i<*out_nhaps;i++)  result_model->hfreq[i]=out_hfreq[i]; }

     	result_model -> coef = (double *) Calloc(ncov,double);
     	if ( result_model -> coef )     {
	  	for (i=0;i<ncov;i++)  result_model->coef[i]=coef[i]; }             
   }
  

if (*verbose==1)  {
		fprintf(file, "\n\n Result Model in hap_shrink_reg : \n");
		print_hmodel(result_model,file);}

  free(out_hap1code);
  free(out_hap2code);
  free(out_hfreq);
  free(deviance);
  for (i=0;i<nrecords;i++) free(x[i]);
  free(x);
  free(out_haplotype_vec);
  free(coef);
  for (i=0;i<*out_nhaps;i++) free(out_haplotype[i]);
  free(out_haplotype);  
  free(out_nhaps);

   return result_model;
   
}

     
 


// This function computes the prediction deviance given the model selected and testing data 

double cross_val(
      hmodel *result_model,	/* the model to be evaluated */
      long  nrecords,      	/* length of hapcode vector for testing data after phasing using all data */
      long  *y,            	/* the response vector, length nrecords */   
      long  nsubj,         	/* the number of unique subjects in testing data */     
      long  *subj_reps,    	/* the vector of replication times for each subject */    
      long  *haplotype_vec,	/* the vector of unique haplotype coding */
      long  nsnps,         	/* the number of SNPs to start with */ 
      long  nhaps,         	/* the number of haplotypes to start with */
      long  *hap1code,     	/* the temparary haplotype coding for each observations */
      long  *hap2code,  		/* the temparary haplotype coding for each observations */ 
      double *hfreq,       	/* the haplotype frequencies */
      long *nngcov,        	/* the number of non-genetic covariates */
      double **ngcov,      	/* the matrix of non-genetic covariates */
      long Minherit,       	/* inheritance mode */
      long *verbose,       	/* 1 turns on debug file and comments */
      long *phase,         	/* 1 is phased, 0 is not */
      FILE *file)

{

      
   long *out_nhaps, *out_hap1code, *out_hap2code, *hap11code, *hap22code, *hapcomcode, nsnp_set;
   double out_hfreq[nhaps]; 
   double **x,  *pred, *outcome, *wgt, *wgtsum, dev_out, mindist, *dist_hap, *locus_freq;
   long **out_haplotype, **model_haplotype;  
   long *out_codebook, *model_codebook, *switch_code, *nreps_test;
   long i,j, count,count1,count2, nrecords_test, *temp1; 


   if (*verbose==1)   fprintf(file,"\n Inside cv function : \n");
   if (*verbose==1)   fprintf(file,"\nthe number of test data=%i\n",nrecords);
   out_hap1code = long_vec(nrecords);
   out_hap2code = long_vec(nrecords);
   hap11code = long_vec(nrecords);
   hap22code = long_vec(nrecords);
   hapcomcode= long_vec(nrecords);
   model_codebook = long_vec(result_model->nhaps);
   out_nhaps = (long *)malloc(sizeof(long));
   nsnp_set=result_model->nsnps;  
   if (*verbose==1)   fprintf(file,"\nbefore cv hap shrink\n");

// Recode haplotypes
   if (nsnps!=nsnp_set) {
		out_haplotype = hap_shrink(
			nrecords,nsnps,nhaps,haplotype_vec,hfreq, hap1code,hap2code,result_model->snp_set,nsnp_set,
			out_nhaps,out_hap1code, out_hap2code, out_hfreq,verbose, file);
	} else {
        *out_nhaps=nhaps;
        for (i=0;i<nrecords;i++) {
                   out_hap1code[i]=hap1code[i];
                   out_hap2code[i]=hap2code[i];
        }
        for (i=0;i<nhaps;i++) {
                   out_hfreq[i] = hfreq[i];
        }
        out_haplotype=long_vec_to_mat(haplotype_vec,nhaps,nsnps);
   }   
    

   model_haplotype=long_vec_to_mat(result_model->haplo_vec,result_model->nhaps,result_model->nsnps);
   out_codebook =  long_vec(*out_nhaps);

  if (*verbose==1)  {
		fprintf(file,"\npass cv hap shrink\n");
  		fprintf(file,"\nmodel_haplotype :\n");
    	print_matrix_long(model_haplotype,result_model->nhaps,result_model->nsnps, file);
  }


// Create out_codebook from all haplotypes  
   for (i=0;i<*out_nhaps;i++) {
      out_codebook[i]=0;
      for (j=0;j<nsnp_set;j++) {
             out_codebook[i] += out_haplotype[i][j] *pow(10,nsnp_set-j-1);
      }
   }

  if (*verbose==1) { 
   fprintf(file,"\nout_codebook :\n");
  	print_vector_long(out_codebook,*out_nhaps, file);
	}


// Create model_codebook from haplotypes in current model (train data) 
   for (i=0;i<result_model->nhaps;i++) {
      model_codebook[i]=0;
      for (j=0;j<result_model->nsnps;j++) {
             model_codebook[i] += model_haplotype[i][j] *pow(10,nsnp_set-j-1);
      }
   }

  if (*verbose==1) { 
	   fprintf(file,"\nmodel_codebook :\n");
  	 	print_vector_long(model_codebook,result_model->nhaps, file);
		fprintf(file,"\nresult_model :\n");
		print_hmodel(result_model,file);
	}


// Retrieve loci and haplo frequencies from results model
   switch_code = long_vec(*out_nhaps); 
   dist_hap =  double_vec(result_model->nhaps);
   locus_freq = double_vec(result_model->nsnps);
  
   for (i=0; i<result_model->nsnps;i++) {
	    locus_freq[i]=0;
            for (j=0;j<result_model->nhaps;j++) {
                if (model_haplotype[j][i]==1) locus_freq[i] += result_model->hfreq[j];
            }
   }
   if (*verbose==1) print_vector_double(locus_freq,result_model->nsnps, file);
   if (*verbose==1) print_matrix_long(model_haplotype,result_model->nhaps,result_model->nsnps, file);
 

 // Assign every haplotype in testing data to the coding used in training data 

  for (i=0;i<*out_nhaps;i++) { 
         count=0;
         for (j=0;j<result_model->nhaps;j++) {
             if (out_codebook[i]==model_codebook[j])  {
                 switch_code[i]=j+1;
                 count++;
                 break;
             }
         }
         
         if (count==0) {
               count1=0;
               mindist=0.0;
               for (j=0;j<result_model->nhaps;j++) {
                         dist_hap[j] = haplo_dist(out_haplotype[i],model_haplotype[j],locus_freq,result_model->nsnps);
                         
                         if (j>0) {
                                count1= (dist_hap[j] < mindist)? j:count1;
                                mindist= (dist_hap[j]< mindist)?dist_hap[j]:mindist;
                         } else mindist=dist_hap[j];
               
               }
               switch_code[i]=count1+1;
         }

  }
   
  if (*verbose==1) {
	fprintf(file,"\nswitch_code = \t\n");
 	print_vector_long(switch_code,*out_nhaps, file);
	fprintf(file,"\nout_hap1code = \t\n");
   print_vector_long(out_hap1code,nrecords, file);
}

// If need to phase 
  if (*phase==0) {
      for (i=0;i<nrecords;i++) {
               out_hap1code[i]=switch_code[out_hap1code[i]-1];
               out_hap2code[i]=switch_code[out_hap2code[i]-1];
               hap11code[i] = (out_hap1code[i] < out_hap2code[i])? out_hap1code[i]:out_hap2code[i];
               hap22code[i] = (out_hap1code[i] > out_hap2code[i])? out_hap1code[i]:out_hap2code[i];
               hapcomcode[i] = 1000* hap11code[i] + hap22code[i];
      }  
      if (*verbose==1) print_vector_long(hapcomcode,nrecords, file);


		// Remove redundant haplotype pairs 
  		outcome=double_vec(nrecords);
  		wgt=double_vec(nrecords);
      count=0;
      count2=0;
  
      wgtsum=double_vec(nsubj);
      nreps_test=long_vec(nsubj);
      if (*verbose==1) print_vector_long(subj_reps,nsubj, file);

      for (i=0;i<nsubj;i++) {
          temp1 = (long *) Calloc(subj_reps[i],long);
          for (j=0;j<subj_reps[i];j++) {
               temp1[j]=hapcomcode[count];
               count++;
          }
          insertSort(temp1,subj_reps[i]);

          for (j=0;j<subj_reps[i];j++) {

             if (j==0)    {
                     hap22code[count2]= temp1[0] % 1000;
                     hap11code[count2]= (temp1[0] - hap22code[count2])/1000;
                     outcome[count2]=y[count-1];
                     wgt[count2]= result_model->hfreq[hap11code[count2]-1] *  result_model->hfreq[hap22code[count2]-1];
                     wgt[count2]= ( hap11code[count2]==   hap22code[count2]) ? 2*wgt[count2]:wgt[count2];
                     wgtsum[i]=wgt[count2];
                     nreps_test[i]=1;
                     count2 ++;
             } else if (temp1[j]!=temp1[j-1]) {
                  hap22code[count2]=temp1[j] % 1000;
                  hap11code[count2]= (temp1[j] - hap22code[count2])/1000;
                  outcome[count2]=y[count-1];
                  wgt[count2]= result_model->hfreq[hap11code[count2]-1] *  result_model->hfreq[hap22code[count2]-1];
                  wgt[count2]= ( hap11code[count2]==   hap22code[count2]) ? 2*wgt[count2]:wgt[count2];
                  wgtsum[i]+=wgt[count2];
                  nreps_test[i]++;
                  count2 ++;
             }
          }      
          free(temp1);   
        }

      if (*verbose==1)  {  
			fprintf(file,"\nhap11code and hap22code : \n");
			print_vector_long(hap11code,count2, file);
			print_vector_long(hap22code,count2, file); 
		}
     
		// Recompute weights 
      count=0;
      for (i=0;i<nsubj;i++) { 
         for (j=0;j<nreps_test[i];j++) {
           wgt[count]= wgt[count]/wgtsum[i];
           count++;
         }
      }

      if (*verbose==1) print_vector_double(wgt,count, file);

      nrecords_test=count2;   
      free(wgtsum);

      free(nreps_test);

  } else {
  
// If do not need to phase   
      for (i=0;i<nrecords;i++) {
               out_hap1code[i]=switch_code[out_hap1code[i]-1];
               out_hap2code[i]=switch_code[out_hap2code[i]-1];
               hap11code[i] = (out_hap1code[i] < out_hap2code[i])? out_hap1code[i]:out_hap2code[i];
               hap22code[i] = (out_hap1code[i] > out_hap2code[i])? out_hap1code[i]:out_hap2code[i];     
               wgt[i] =1;
               outcome[i] = y[i];
      }  
      nrecords_test = nrecords;
      
  }


// Create design matrix
	// Intercept term is added here; we omit haplo with largest frequency
	// design matrix =  [ observations x (intercept, haplos, non-genetic covariates) ]
	// nngcov represents the number of non-genetic covariates
	// out_nhaps and nngcov represent the number of haplotypes and non-genetic covariates, respectively 
  
   //x=double_matrix(nrecords_test,result_model->nhaps);  
	x=double_matrix(nrecords_test, (result_model->nhaps + result_model->nngcov) );  
   pred=double_vec(nrecords_test);  
   dev_out=0;

   if (*verbose==1) {
      fprintf(file,"\nin cv, nhaps nngcov = %i %i\t\n", result_model->nhaps, result_model->nngcov);
      fprintf(file,"\nempty x in cv :\n");
      print_matrix_double(x,nrecords_test,(result_model->nhaps + result_model->nngcov), file);
   }

   for (i=0; i< nrecords_test; i++) {       
       x[i][0]=1;
		 // Compute haplo codes based on inheritance and put in design matrix 
       for (j=1;j<(result_model->nhaps+1);j++) {
           if (j <result_model->hapbase) {   
               x[i][j]=0; 
               if (Minherit==1) {
                  if (j==hap11code[i]) x[i][j] += 1;                               
                  if (j==hap22code[i]) x[i][j] += 1;
               } 
               if (Minherit==2) {
                  if (j==hap11code[i] || j==hap22code[i]  ) x[i][j] += 1;                               
               }
               if (Minherit==3) {
                  if (j==hap11code[i] && j==hap22code[i]  ) x[i][j] += 1;                               
               }    
           }
           if (j >result_model->hapbase) {
               x[i][j-1]=0;
               if (Minherit==1) {
                  if (j==hap11code[i]) x[i][j-1] += 1;                               
                  if (j==hap22code[i]) x[i][j-1] += 1;
               } 
               if (Minherit==2) {
                  if (j==hap11code[i] || j==hap22code[i]  ) x[i][j-1] += 1;                               
               }
               if (Minherit==3) {
                  if (j==hap11code[i] && j==hap22code[i]  ) x[i][j-1] += 1;                               
               }    
           } 
        }


		 // Add non-genetic covariates to design matrix 
   	 if (*nngcov!=0)  {
        	for (j=(result_model->nhaps); j<(result_model->nhaps + result_model->nngcov); j++) {
                x[i][j]=ngcov[i][(j-result_model->nhaps)];
         }
    	 }

   if (*verbose==1) {
          fprintf(file,"\n filling in x in cv :\n");
          print_matrix_double(x,nrecords_test,(result_model->nhaps + result_model->nngcov), file);
   }


	  // Prediction
	  // Multipling one row at a time
     pred[i]=0;
     //for (j=0;j<result_model->nhaps;j++) pred[i] += x[i][j]*result_model->coef[j];
     for (j=0;j<(result_model->nhaps + result_model->nngcov);j++) {
		pred[i] += x[i][j]*result_model->coef[j];
	  }
     pred[i]=exp(pred[i])/(1+exp(pred[i]));
     dev_out += (outcome[i]==1.0)? (wgt[i]*log(pred[i])):(wgt[i]*log(1-pred[i]));
   }    

   if (*verbose==1) {
                  fprintf(file,"\n\nthe number of testing data =%i\n",nrecords_test);
                  fprintf(file,"\n\nx, pred, outcome : \n");
                  print_matrix_double(x,nrecords_test,(result_model->nhaps + result_model->nngcov), file);
                  print_vector_double(pred,nrecords_test,file);
                  print_vector_double(outcome,nrecords_test,file);
   }


   free(out_hap1code);
   free(out_hap2code);
   free(hap11code);
   free(hap22code);
   free(hapcomcode);
   free(model_codebook);
   free(out_codebook);
   for (i=0;i<nrecords_test;i++)  free(x[i]);
   free(x);
   free(pred);
   free(switch_code);
   free(locus_freq);
   free(dist_hap);
   free(outcome);
   free(wgt);

   for (i=0;i<result_model->nhaps;i++)  free(model_haplotype[i]);
   free(model_haplotype);
   for (i=0;i<*out_nhaps;i++)  free(out_haplotype[i]);
   free(out_haplotype);
   free(out_nhaps);
 
   if (*verbose==1)  fprintf(file,"\nthe cv deviance is%8.5f\n",dev_out);  
   if (*verbose==1)  fprintf(file,"\n Exit cv function \n");  

   return(dev_out);
}






// This function performs the stepwise search using BIC

void stepwise_search_alpha(long *indx_subj,
			   long *nsubj,
			   long *nobs, 
			   long *subj_rep,
			   double *wgt, 
			   long *csctl,
			   long *nloci,
			   long *nhap,
			   long *hap1code,
			   long *hap2code,
			   long *uhap,
			   double *happrob,
			   long *bestsize,
			   long *Minherit, 
			   double *deviance, 
			   long  *maxsnps,
				long *nngcov,
				double *ngcov_vec,
			   double *tol,
			   double *alpha,
			   long *verbose,
			   long *phase, 
			   long *output_snp_set, 
			   long *out_haplo_vec, 
			   double *out_hap_freq, 
			   double *varstore, 
			   double *coef, 
			   long *out_nhaps) {

  /**   long *indx_subj              vector of subject ids, length nobs  **/
  /**   long *nsubj,                 total number of distinct subjects  **/
  /**   long *nobs,                  total number of observations  **/ 
  /**   long *subj_rep,             the vector of number of possible haplotype pairs  **/
  /**   double *wgt,                 the vector of probabilities of haplotype pairs    **/
  /**   long *nngcov,                the number of non-genetic covariates    **/
  /**   double **ngcov,               matrix of non-genetic covariates    **/
  /**   long *csctl,                 vector of case-control status, length nobs  **/
  /**   long *nloci,                 total number of SNPs considered  **/
  /**   long *nhap,                  number of haplotypes phased using all data , to start with **/
  /**   long *hap1code,              the vector of codes for haplotype 1, length nobs  **/
  /**   long *hap2code,              the vector of codes for haplotype 2, length nobs  **/
  /**   long *uhap,                  vector of unique haplotype  **/
  /**   double *happrob,             the vector of haplotype probabilities  **/
  /**   double *alpha,               the tuning parameter in AIC **/
  /**   double *tol                  the convergence parameter  **/
  /**   long *phase                  the indicator of whether phase is know or not **/

  long count, count1, count2, i,j, k, maxit=1000, *snp_set, *current_snp_set;
  double *mustart, min_phi, ncase, bestdev, temp_phi, mindev, min_freq,meany, **ngcov_mat;  
  hmodel *temp, *best_model, *base_model; 
  Node SNP_out;
  Link1 current_snp= &SNP_out, next_snp;   /** current_snp is the address of SNP_out **/

  FILE *file; 
  file = fopen("debug_stepwise.txt", "w"); 
     

// Put covariate data in matrix
   if (*nngcov>1) {
      ngcov_mat = double_vec_to_mat(ngcov_vec, *nobs, *nngcov);
      }
 
// Prepare the linked list of all snps to start the search 
  SNP_out.next=NULL;

// Construct SNP_out via current_snp 
  for (j=0;j<*nloci;j++) {
    /** alloc memory for next_snp **/
    next_snp = malloc(sizeof(Link1));
    next_snp->snp_pos = j+1;
    next_snp->next = NULL;

    current_snp->next = next_snp;
    current_snp=next_snp;
  }

  if (*verbose==1){
    fprintf(file,"SNP_out (start)\n");    
    print_list_long(SNP_out, file);
    fprintf(file,"SNP_out (stop)\n");    
  }
             

// Model searching

	// given the limit of maxsnps,
	// 1. compute mustart and min_freq
 	// 2. traverse the linked list to select a snp and subset the data on that snp
   // 3. call hap_shrink_reg to compute the deviance of the train model w/ those specific snps
   // 4. select best model for that size based on deviance
   // 5. cross-val on that model size
   // 6. repeat until have grown the model to size = maxsnp  

	// count variable tracks all snps, coun1 tracks snps in snp_out
	
  mustart=double_vec(*nobs);
  snp_set = long_vec(*maxsnps);

  for (j=0;j<*nobs;j++){              	  
    mustart[j]= (wgt[j]*csctl[j]+0.5)/(wgt[j]+1);
      }
  min_freq= (2.0/(2.0 * (*nsubj))); 
  if (*verbose==1) 
    fprintf(file,"\n the freq cutoff = %8.9f",min_freq);
  
  // Compute the null deviance 
  ncase =0.0;
  for (j=0; j<*nobs; j++) {
    ncase+=csctl[j]*wgt[j];
  }

  meany = ncase/(*nsubj);
  for (j=0;j<*nobs;j++)   
    deviance[0] += (csctl[j]==1) ? wgt[j]*log(meany):wgt[j]*log(1-meany);
         
  deviance[0] = (-2)*deviance[0] + *alpha;

	// Build model one snp at a time up to maxsnps :		     	 
  count=0; /** # of loop @ while (count < *maxsnps) **/	
  while (count < *maxsnps) {
    current_snp = &SNP_out;

    count1=0; /** # of loop @ while (current_snp->next !=NULL) **/
    while (current_snp->next !=NULL) {
      next_snp=current_snp -> next;
      snp_set[count]= next_snp -> snp_pos; /** QQQ: why? **/
      if (*verbose==1) {
			fprintf(file,"\nthe current snp set 'snp_set' (start)\n");
			fprintf(file,"count = %d\n", count);
			fprintf(file,"count1 = %d\n", count1);
			print_vector_long(snp_set, count+1, file);
			fprintf(file,"\nthe current snp set 'snp_set' (stop)\n");
      }
     

	// within call to hap_shrink_reg :
	// 1. calls hap_shrink (recode haplos)
	// 2. call haplo_cluster
	// 3. iwls
	// 4. output result model	
		 temp=hap_shrink_reg(
			csctl, uhap,*nloci, *nhap,hap1code,hap2code,wgt,happrob,*nobs,subj_rep,snp_set, 
			count+1,min_freq, mustart,
			nngcov, ngcov_mat,
			maxit,*Minherit, *tol, verbose,file); 
	        
	// Compare deviance in current model with previous best 
      if  (count==0) {
			if (count1==0) { /** count==0 & count1==0 **/
			  best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
			  copy_hmodel(temp,best_model);
			}
			else { /** count==0 & count1!=0 **/
			  if  (best_model->dev < temp -> dev)  {
			    Free_hmodel(best_model);
			    best_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
			    copy_hmodel(temp,best_model);
			  }
			}
   		   } else {	/** if count != 0 **/
			/** if count != 0, start to calculate phi **/
			if (count1==0) {
			  base_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
			  copy_hmodel(temp,base_model);
			  min_phi =  (-2)*(temp->dev) + *alpha * temp->nhaps;  
			} else {
			  temp_phi=  (-2)*(temp->dev) + *alpha * temp->nhaps;  
			  /** if temp_phi is less than current min_phi, the temp model is the new best model & the temp_phi is the new min_phi **/
			  if (temp_phi < min_phi) {
			    Free_hmodel(base_model);
			    base_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
			    copy_hmodel(temp,base_model);
			    min_phi=temp_phi;    
			  }
			}
      }

      if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 
      if (*verbose==1){
	fprintf(file,"\nthe current best model (start)\n");
	fprintf(file,"count = %d\n", count);
	fprintf(file,"count1 = %d\n", count1);
	print_hmodel(best_model,file);
	fprintf(file,"\nthe current best model (stop)\n");
      }

	// Increment the snp and build the next model of size count
      current_snp = next_snp;
      Free_hmodel(temp);
      count1++;
    } /** end of while (current_snp->next !=NULL) **/
                 

	// Deviance comparison of last model of size count 
    if (count==0){
      deviance[count+1]=(-2)*(best_model->dev) + *alpha * best_model->nhaps; 
      snp_set[count]= best_model->snp_set[count]; 
    } else {
      deviance[count+1]=(-2)*(base_model->dev) + *alpha * base_model->nhaps;
      snp_set[count]= base_model->snp_set[count];
    }

    // Determine the minimum deviance 
    if (count==0) {
      mindev = (-2*best_model->dev)+ *alpha * best_model->nhaps;
    } else { /** count != 0 **/
      if (deviance[count+1] < mindev) {
			Free_hmodel(best_model);
			best_model=new_hmodel(base_model->nsnps,base_model->nhaps,base_model->nngcov,base_model->nrec);
			copy_hmodel(base_model,best_model);
			mindev = deviance[count+1];
      } 
      if (*verbose==1) fprintf(file,"\nnow the base model");
      if (*verbose==1)  print_hmodel(base_model,file);
      Free_hmodel(base_model);
    }
   
	 if (*verbose==1) fprintf(file,"\nnow the best model");
    if (*verbose==1)  print_hmodel(best_model,file);
    if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of phi=%8.5f\n",min_phi); 

           
    if (*verbose==1)   print_vector_long(snp_set,count+1,file);

	// Delete the selected SNPs from SNP_out (don't want same snp in model twice) 
    current_snp = &SNP_out;
    while (current_snp->next !=NULL) {                 
      next_snp=current_snp -> next;  
      if (next_snp -> snp_pos ==  snp_set[count])  {
			current_snp ->next = next_snp->next;
			free(next_snp);
      } else {
		current_snp = next_snp;
      }
    }

    if (*verbose==1)     
      print_list_long(SNP_out,file);
     
	// Increment count and move to next model size 
    count++;
  } /** end of while (count < *maxsnps) **/


  // Pruning
  if (*verbose==1) 
    fprintf(file,"\nstart pruning \n") ;
    
  count2 =count+1;
  while (count>1) { 
    current_snp_set = long_vec(count-1);
    count1 =0;
    for (j=0;j<count;j++) {
      
		// Create current_snp_set based on 
      for (k=0;k<(count-1);k++)  {
	if (k<j) current_snp_set[k] = snp_set[k];  
	else current_snp_set[k] = snp_set[k+1];  
      }

      temp=hap_shrink_reg(
			csctl, uhap,*nloci, *nhap,hap1code,hap2code,wgt,happrob,*nobs,subj_rep,current_snp_set, 
			count-1,min_freq, mustart,
			nngcov, ngcov_mat,
			maxit, *Minherit,*tol, verbose,file); 
      
      temp_phi= (-2)*(temp->dev) + *alpha * temp->nhaps ;

      if (count1==0) { 
	/** BUG Fixed Here (2009-05-14) **/
	base_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
	copy_hmodel(temp,base_model);
	/* copy_hmodel(temp,best_model); */

	min_phi =  temp_phi;
      } else { /** count1 != 0 **/
	if (temp_phi < min_phi) {
	  Free_hmodel(base_model);
	  base_model=new_hmodel(temp->nsnps,temp->nhaps,temp->nngcov,temp->nrec);
	  copy_hmodel(temp,base_model);
	  min_phi=temp_phi;    
	}
      }
      count1++;
      Free_hmodel(temp);
    } /** end of j **/
    
    deviance[count2]=(-2)*(base_model->dev) + *alpha * base_model->nhaps;   
    if ((count>0) & (*verbose==1)) fprintf(file,"\nthe value of current dev=%8.5f\n",deviance[count2]);  
    if (deviance[count2] < mindev) {
      Free_hmodel(best_model);
      best_model=new_hmodel(base_model->nsnps,base_model->nhaps,base_model->nngcov,base_model->nrec);
      copy_hmodel(base_model,best_model);
      mindev = deviance[count2];
    }
    
	 if ((count>0) & (*verbose==1)) {
		fprintf(file,"\nthe value of minimal dev=%8.5f\n",mindev);    
      fprintf(file,"\nthe value of phi=%8.5f\n",min_phi);
	 	} 
    if (*verbose==1) {
		fprintf(file,"\nnow the base model");
      print_hmodel(base_model,file);  
      fprintf(file,"\nnow the best model");
      print_hmodel(best_model,file);
		}
    
	 free(snp_set);
    snp_set=long_vec(count-1);
    for (k=0;k<(count-1);k++)  snp_set[k] = base_model->snp_set[k]; 
    free(current_snp_set);
    Free_hmodel(base_model);
    if (*verbose==1)   print_vector_long(snp_set,count-1,file);
    count--;
    count2++;

  } /** end of while(count>1) **/

 	bestdev =  (-2)*best_model->dev + *alpha*best_model->nhaps ;
	if (deviance[0] < bestdev) *bestsize=0; else {
                       
                 for (i=0;i< best_model->nsnps;i++) output_snp_set[i]=best_model->snp_set[i]; 

                 varest(best_model,nsubj,nobs,subj_rep,csctl,nloci,nhap,hap1code,hap2code,
						uhap,happrob,wgt,output_snp_set,best_model->nsnps,min_freq,
						nngcov, ngcov_mat,
						*Minherit,verbose,file,phase,varstore);

                 for (i=0;i<(best_model->nhaps+best_model->nngcov);i++) {
                   coef[i]=best_model->coef[i];
                   out_hap_freq[i] = best_model->hfreq[i];
                 }

                 for (i=0;i<(best_model->nsnps * best_model->nhaps);i++) out_haplo_vec[i]=best_model->haplo_vec[i];
     
                 *out_nhaps=best_model->nhaps;   
                 *bestsize=best_model->nsnps;  
	}
	 
	         if (*verbose==1)   print_vector_long(snp_set,count,file);   
        
                 Free_hmodel(best_model);
		 /*if (*verbose==1)    fprintf(file,"\n BUG HERE?\n");*/

		 free(snp_set);
                 free(mustart);
     
                 current_snp = &SNP_out;
                 while (current_snp->next !=NULL) {                 
                         next_snp=current_snp -> next;  
                         current_snp ->next = next_snp->next;
                         free(next_snp);                     
                 }

       fclose(file);

}   
		    
		         
		    










long *insertionSort(long *numbers, long array_size, long *order)
{
  long i, j, index, *sorted_numbers, index1;
  
  sorted_numbers=long_vec(array_size);

  for (i=0; i < array_size; i++) sorted_numbers[i]=numbers[i];
  for (i=0; i < array_size; i++)
  {
    index = sorted_numbers[i];
    index1= order[i];
    j = i;
    while ((j > 0) && (sorted_numbers[j-1] > index))
    {
      sorted_numbers[j] = sorted_numbers[j-1];
      order[j]=order[j-1];
      j = j - 1;
    } 
    sorted_numbers[j] = index;
    order[j]=index1;
  }
  return sorted_numbers;
}





void insertSort(long *numbers, long array_size)
{
  long i, j, index;
  
  for (i=0; i < array_size; i++)
  {
    index = numbers[i];
    j = i;
    while ((j > 0) && (numbers[j-1] > index))
    {
      numbers[j] = numbers[j-1];
      j = j - 1;
    } 
    numbers[j] = index;
  }
  return;
}



double **double_vec_to_mat(double *Yvec, long nrow, long ncol){
	   long i,j,k=0;
	   double **Y;
           Y=double_matrix(nrow,ncol);
// BUG FIX 2010-07-15 : with correction will now fill by column 
//	   for (i=0;i<nrow;i++){
	   for (i=0;i<ncol;i++){
		   for (j=0;j<nrow;j++) { 
			    Y[j][i]=Yvec[k];
				k++;
		   }
	   }
	   return Y;
}

// This function was designed to work on the non-genetic covariates
// The matrix should be flattened by column
double *double_mat_to_vec(double **Ymat, long nrow, long ncol){
	   long i,j,count=0;
	   double *Y;
      Y=double_vec(nrow*ncol);
		for (j=0;j<ncol;j++) {
	   	for (i=0;i<nrow;i++){
			   Y[count]= Ymat[i][j];
				count++;
		   }
	   }
	   return Y;
}


long *long_mat_to_vec(long **Ymat, long nrow, long ncol){
	   long i,j,count=0;
	   long *Y;
           Y=long_vec(nrow*ncol);
	   for (i=0;i<nrow;i++){
		   for (j=0;j<ncol;j++) {
			   Y[count]= Ymat[i][j];
				count++;
		   }
	   }
	   return Y;
}



long **long_vec_to_mat(long *Yvec, long nrow, long ncol){
	   long i,j,k;
	   long **Y;
           Y = long_matrix(nrow,ncol);
           k=0;
	   for (i=0;i<nrow;i++){
		   for (j=0;j<ncol;j++) {
			    Y[i][j]=Yvec[k];
				k++;
		   }
	   }
	   return Y;
}





double **double_matrix(long nrow, long ncol){
	 long i;
	 double **m;
	 m=(double **) Calloc(nrow, double *);
	 if (!m) errmsg("mem alloc failure 1 in double_matrix");
	 for (i=0;i<nrow;i++) {
		  m[i]=(double *) Calloc(ncol,double);
		  if (!m[i]) errmsg("mem alloc failure 2 in double_matrix");
	 }
	 return m;
}





long **long_matrix(long nrow, long ncol){
	 long i;
	 long **m;
	 m=(long **) Calloc(nrow, long *);
	 if (!m) errmsg("mem alloc failure 1 in long_matrix");
	 for (i=0;i<nrow;i++) {
		  m[i]=(long *) Calloc(ncol,long);
		  if (!m[i]) errmsg("mem alloc failure 2 in long_matrix");
	 }
	 return m;
}




double *double_vec(long n){

	 double *v;
	 v=(double *) Calloc(n, double);
	 if (!v) errmsg("mem alloc failure in double_vec");
	 return v;
}

long *long_vec(long n){
	 long *v;
	 v=(long *) Calloc(n, long);
	 if (!v) errmsg("mem alloc failure in long_vec");
	 return v;
}


long max_long(long x, long y){
     long z;
     z=(x>=y)?x:y;
     return z;
}


long min_long(long x, long y){
     long z;
     z=(x>=y)?y:x;
     return z;
}




double max_double(double x, double y){
     double z;
     z=(x>=y)?x:y;
     return z;
}



double min_double(double x, double y){
     double z;
     z=(x>=y)?y:x;
     return z;
}


void print_matrix_long(long **m, long nrow, long ncol, FILE *file){
  long i,j;
  for (i=0;i<nrow;i++) {
    for (j=0;j<ncol;j++){
          if (j==0) fprintf(file,"\n%i",m[i][j]);
          else fprintf(file,"\t%i",m[i][j]);
    }
  }
}


void print_matrix_double(double **m, long nrow, long ncol, FILE *file){
  long i,j;
  for (i=0;i<nrow;i++) {
    for (j=0;j<ncol;j++){
          if (j==0) fprintf(file,"\n%8.5f",m[i][j]);
          else fprintf(file,"\t%8.5f",m[i][j]);
    }
  }
}



void print_vector_double(double *m, long n, FILE *file){
    long j; 
    for (j=0;j<n;j++){
          if (j==0) fprintf(file,"\n%8.5f",m[j]);
          else fprintf(file,"\t%8.5f",m[j]);
    }
  
}



void print_vector_long(long *m, long n, FILE *file){
    long j;
    for (j=0;j<n;j++){
          if (j==0) fprintf(file,"\n%i",m[j]);
          else fprintf(file,"\t%i",m[j]);
    }
    return;
  
}



void print_list_long(Node m, FILE *file){
    Link1 current, next;
    long count=0;
    current= &m;
    while (current->next !=NULL) {
          next = current->next;
          if (count==0) fprintf(file,"\n%i",next->snp_pos);
          else fprintf(file,"\t%i",next->snp_pos);
          current=next; 

          count++;
    } 
    return;
}


void print_hmodel(hmodel *m, FILE *file) {
  long **x, i;
  double **c;
    fprintf(file,"\nthe current number of haplotypes=%i\n",m->nhaps);
    fprintf(file,"\nthe current number of snps=%i\n",m->nsnps);
    fprintf(file,"\nthe current deviance=%8.5f\n",m->dev);
    fprintf(file,"\nthe current set of SNPs=\n");
    print_vector_long(m->snp_set,m->nsnps,file);
    fprintf(file,"\nthe current haplotype matrix\n");
    x = long_vec_to_mat(m->haplo_vec,m->nhaps, m->nsnps);
    print_matrix_long(x,m->nhaps,m->nsnps,file);
    fprintf(file,"\nthe current non-genetic covariates\n");
    if (m->nngcov!=0) { 
    c = double_vec_to_mat(m->ngcov_vec,m->nrec, m->nngcov);
    print_matrix_double(c,m->nrec,m->nngcov,file);
    } else { fprintf(file,"None\n"); }
    fprintf(file,"\nthe hap base is=%i\n",m->hapbase);
    fprintf(file,"\nthe haplotype frequencies=\n");
    print_vector_double(m->hfreq,m->nhaps,file);
    fprintf(file,"\nthe coefficients=\n");
    print_vector_double(m->coef,(m->nhaps+m->nngcov),file);  
    for (i=0;i<m->nhaps;i++) free(x[i]); 
    free(x);
    if (m->nngcov!=0) {
		for (i=0;i<m->nrec;i++) free(c[i]); 
    free(c);
    }
    return;
}


hmodel *new_hmodel(long nsnps,long nhaps, long nngcov, long nrec) {
    hmodel *x;
    x=(hmodel *)Calloc(1, hmodel);
    if (x) {
        x-> nhaps = nhaps;
        x-> nngcov = nngcov;
        x-> nsnps = nsnps;
        x-> nrec = nrec;
        x-> snp_set = (long *) Calloc(nsnps,long);
        x-> haplo_vec = (long *) Calloc(nhaps*nsnps,long);
   if (nngcov!=0) x-> ngcov_vec = (double *) Calloc(nrec*nngcov,double);
        x-> hfreq = (double *) Calloc(nhaps,double);
        x-> coef = (double *) Calloc(nhaps+nngcov,double);
   }
   return x;
}


void copy_hmodel(hmodel *current, hmodel *best){
     long i;
     best->nhaps=current->nhaps;
     best->nngcov=current->nngcov;
     best->nsnps=current->nsnps;
     best->dev=current->dev;
     best->hapbase=current->hapbase;
     best->nrec=current->nrec;
     for (i=0;i<current->nsnps;i++) best->snp_set[i]=current->snp_set[i];
     for (i=0;i<(current->nsnps * current->nhaps);i++) best->haplo_vec[i]=current->haplo_vec[i];
     if (current->nngcov!=0) { 
			for (i=0;i<(current->nrec * current->nngcov);i++) best->ngcov_vec[i]=current->ngcov_vec[i];
	  }
     for (i=0;i<current->nhaps;i++) best->hfreq[i]=current->hfreq[i]; 
	  for (i=0;i<(current->nhaps + current->nngcov);i++) best->coef[i]=current->coef[i]; 
}  


void Free_hmodel(hmodel *x){
  if (x!=NULL) {
     if (x->snp_set!=NULL) free(x->snp_set);
     if (x->haplo_vec!=NULL) free(x->haplo_vec);
     if (x->ngcov_vec!=NULL) free(x->ngcov_vec);
     if (x->hfreq!=NULL) free(x->hfreq);
     if (x->coef!=NULL) free(x->coef);
     free(x);
  }
  return;
}




/***********************************************************************************/

static void errmsg(char *string){

  /* Function to emulate "stop" of S+ - see page 134, S Programing, by
     Venables and Ripley */

   PROBLEM "%s", string RECOVER(NULL_ENTRY);
}

