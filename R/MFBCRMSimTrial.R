#' Provides simulation results using MFBCRM
#' @param max_samp Total number of patients recruited/will be recruited in the trial.
#' @param Cohort Number of patients within each cohort.
#' @param ptrue True toxicity probability vector.
#' @param mu_mat Prior expected toxicity probability matrix at each dose.
#' @param p_rho Prior probability that two dose-toxicity probabilities will not cluster together.
#' @param sigma Prior standard deviation for the parameter alpha.
#' @param mtd Maximum Tolerated dose toxicity probability (pre defined).
#' @param p_u Cut-off toxicity probability for first dose.
#' @param B Number of Iterations to run for MCMC.
#' @param M Number of simulations to run.
#' @return A list containing (1) Design parameters and prior hyperparameters used for simulating the trials and (2) a summary of the trial simulation results including the percent of times each dose was selected and the average number of toxicities seen in the trial.
#' @examples
#' max_samp=15   ##Change to larger size
#' Cohort=3
#' ptrue=c(0.01,0.05,0.15,0.3,0.45,0.5,0.6,0.8)
#' ##Skeletons for 8 doses
#' mu1=c(0.02,0.06,0.08,0.12,0.2,0.3,0.4,0.5)
#' mu2=c(0.01,0.05,0.09,0.14,0.18,0.22,0.26,0.30)
#' mu3=c(0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80)
#' mu4=c(0.20,0.30,0.40,0.50,0.60,0.65,0.70,0.75)
#' mu_mat=matrix(c(mu1,mu2,mu3,mu4),nrow = 4,byrow = TRUE)
#' p_rho=0.9
#' sigma = 2
#' mtd = 0.3
#' p_u=0.9
#' B=200 ##Number of iterations, change to 2k
#' M=10 ##Number of simulations, change to larger
#' Z=MFBCRMSimTrial(max_samp,Cohort,ptrue,mu_mat,p_rho,sigma,mtd,p_u,B,M)
#' Z
#'@export

MFBCRMSimTrial=function(max_samp,Cohort,ptrue,mu_mat,p_rho,sigma,mtd,p_u,B,M){

  names(ptrue)=c(paste0("Dose",1:length(ptrue)))
  colnames(mu_mat)=c(paste0("Dose",1:ncol(mu_mat)))
  rownames(mu_mat)=c(paste0("Skeleton",1:nrow(mu_mat)))

  HOLD_SAMP=0
  ##Check for errors in dimension specification
  if(max_samp%%Cohort!=0){
    HOLD_SAMP=1
  }

  HOLD_LENGTH=0

  ##Check if the dimensions of ptrue and mu are same or not
  if(length(ptrue)!=ncol(mu_mat)){
    HOLD_LENGTH=1
  }

  if(HOLD_SAMP==1){
    message("Total expected number of patients should be divisible by the cohort size.")
  }else if(HOLD_LENGTH==1){
    message("Length of vector ptrue and column length of matrix mu_mat should be same.")
    }else{


    ###Contains Design parameters
    DESIGN = as.list(rep(NA,10))

    names(DESIGN) = c("Total number of patients expected for the trial:",
                      "Cohort size for dose assignment = ",
                      "True toxicity probabilities for each dose:",
                      "Initial skeletal toxicity probability matrix for each dose level:",
                      "Probability that a dose will not cluster with another dose = ",
                      "Prior standard deviation of alpha = ",
                      "Dose Limiting Toxicity probability = ",
                      "Cut-off probability, if the 1st dose is too toxic",
                      "Number of MCMC iterations = ",
                      "Number of simulations ="
                      )

    DESIGN[[1]]=max_samp
    DESIGN[[2]]=Cohort
    DESIGN[[3]]=ptrue
    DESIGN[[4]]=mu_mat
    DESIGN[[5]]=p_rho
    DESIGN[[6]]=sigma
    DESIGN[[7]]=mtd
    DESIGN[[8]]=p_u
    DESIGN[[9]]=B
    DESIGN[[10]]=M

    TIME = paste0("Model run on: ",Sys.time())

    DESIGN=c(TIME,DESIGN)

    names(DESIGN)[[1]]="Date/Time of escalation decision"

    DESIGN = c("FBCRM Package Version: 1.0.0",DESIGN)

    ####constants to calculate the numerical integration for CRM part
    a=-3
    b=3
    n=1000

    ##This will provide a list
    RESULTS=MFBCRM_RUNTRIAL(Cohort,max_samp,ptrue,mu_mat,p_rho,sigma,mtd,p_u,B,M,a,b,n)

    dose_chosen=table(factor((RESULTS[[3]]),levels=1:length(ptrue)))

    ##Average toxicity
    Avg_Tox=sum(RESULTS[[2]])/M

    optimal_d=RESULTS[[3]]

    #############Extra for X>mtd checking###############

    #patient allocation over B simulations
    p_allocation=RESULTS[[1]]
    #Avg number of toxicity events
    p_avgtox=RESULTS[[2]]

    avg_patient=rep(NA,length(ptrue))
    avg_tox_each_dose=avg_patient

    for(j in 1:length(ptrue)){
      ##Avg. number of patient allocated in each dose
      avg_patient[j]=round(mean(p_allocation[,j]),1)
      ##Avg. number of toxicity events in each dose
      avg_tox_each_dose[j]=round(mean(p_avgtox[,j]),1)
    }

    WHICH1=which(1:length(ptrue)>(optdose(ptrue,mtd)+1))

    ##Average no. of people allotted above mtd
    avg_aboveMTD=sum(RESULTS[[1]][,WHICH1])/M


    #################End of extra ###############

    delta=rep(NA,M)
    sdelta=0

    for(l in 1:M){

      if(optimal_d[l]==0){
        delta[l]=abs(ptrue[optdose(ptrue,mtd)+1])
      }else{
        delta[l]=mean(abs(ptrue[optdose(ptrue,mtd)+1]-ptrue[optimal_d[l]]))
      }
      sdelta=sdelta+delta[l]
    }

    #####delta
    DELTA=round(sdelta/M,2)

    ########PSEL
    PSEL=as.numeric(dose_chosen[[optdose(ptrue,mtd)+1]]*100/M)

  }

  Z=as.list(rep(0,8))
OUT1=NA
  
  
    ##If the trial stopped in some of the simulated scenarios because the first dose was too toxic then
    if(sum(dose_chosen*100/M)<100){

OUT1=paste0("The trial stopped because the 1st dose was too toxic in ", (100-sum(dose_chosen*100/M)),
        "% of the times among M trials.")

    }

  Z[[1]]=PSEL

  Z[[2]]=dose_chosen*100/M
  names(Z[[2]])=c(paste0("Dose",1:length(ptrue)))

  Z[[3]]=Avg_Tox

  Z[[4]]=avg_patient

  Z[[5]]=avg_tox_each_dose

  Z[[6]]=avg_aboveMTD

  Z[[7]]=DELTA
  
  Z[[8]]=OUT1

  names(Z)=c(paste0("Probability % of selecting the correct dose as MTD over ",M," simulated trials"),
             paste0("Probability % of selecting a dose as MTD over ",M," simulated trials"),
             paste0("Average toxicity among ",max_samp," patients over ",M," simulated trials"),
             paste0("Average patients treated per dose over ",M," simulated trials"),
             paste0("Average number of DLT per dose over ",M," simulated trials"),
             "Average number of patients treated above true MTD",
             paste0("Mean deviation for true MTD and the selected MTD over ",M," simulated trials"),
             "Stopping Probability")

  ###Write the dataframe into the last item of the list
  Z1=as.list(c(1,2))
  Z1[[1]]=DESIGN
  Z1[[2]]=Z
  names(Z1)[[1]]="Design Parameters"
  names(Z1)[[2]]="Results"

return(Z1)

}


