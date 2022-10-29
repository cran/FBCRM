#' Provides the optimal dose level closest to the mtd where the next cohort of patients should be allotted based on the data
#' @param X Vector of patients allotted to each dose level .
#' @param Y Vector of toxicity events in each dose.
#' @param Cohort Number of patients within each cohort.
#' @param mu_mat Prior expected toxicity probability matrix at each dose.
#' @param p_rho Prior probability that two dose-toxicity probabilities will not cluster together.
#' @param sigma Prior standard deviation for the parameter alpha.
#' @param mtd Maximum Tolerated dose toxicity probability (pre defined).
#' @param B Number of Iterations to run for MCMC.
#' @param p_u Cut-off toxicity probability for first dose
#' @return A list containing (1) Design parameters and prior hyperparameters used for running the trials and (2) a posterior summary of the resuls, including the next dose to assign patients to.
#' @examples
#' X=c(3, 6, 3, 3, 3, 9, 15, 6)
#' Y=c(1, 0, 1, 0, 0, 2,  4, 5)
#' Cohort=3
#' ##Skeletons for 8 doses
#' mu1=c(0.02,0.06,0.08,0.12,0.2,0.3,0.4,0.5)
#' mu2=c(0.01,0.05,0.09,0.14,0.18,0.22,0.26,0.30)
#' mu3=c(0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80)
#' mu4=c(0.20,0.30,0.40,0.50,0.60,0.65,0.70,0.75)
#' mu_mat=matrix(c(mu1,mu2,mu3,mu4),nrow = 4,byrow = TRUE)
#' p_rho=0.9
#' sigma = 2
#' mtd = 0.3
#' B=2000 ##Number of iterations
#' p_u=0.9
#' Z=GetMFBCRM(X, Y, Cohort, mu_mat, p_rho, sigma, mtd, B, p_u)
#' Z
#'@export

GetMFBCRM=function(X,Y,Cohort,mu_mat,p_rho,sigma,mtd,B,p_u){

  names(X)=c(paste0("Dose",1:length(X)))
  names(Y)=c(paste0("Dose",1:length(Y)))
  colnames(mu_mat)=c(paste0("Dose",1:ncol(mu_mat)))
  rownames(mu_mat)=c(paste0("Skeleton",1:nrow(mu_mat)))

  DATAHOLD = data.frame(rbind(X,Y,mu_mat))
  colnames(DATAHOLD)=c("Subjects assigned","Toxicity events","Initial skeletal probability")

  ERRHOLD=c(length(X), length(Y), ncol(mu_mat))

  HOLD=0
  ##Check for errors in dimension specification
  for(k in 1:length(ERRHOLD)){
    for(m in 1:length(ERRHOLD)){
      if(ERRHOLD[k] != ERRHOLD[m]){
        HOLD=1
      }
    }
  }

  if(HOLD==1){
    message("Subjects assigned at each dose, toxicity events at each dose, or initial skelatl probability matrix
            has incorrect dimensions")
  }else{


    ###Contains Design parameters
    DESIGN = as.list(rep(NA,9))

    names(DESIGN) = c("Subjects assigned to each dose level:",
                      "No. of toxicity events observed in each dose level:",
                      "Cohort size for dose assignment = ",
                      "Initial skeletal toxicity probability matrix for each dose level:",
                      "Probability that a dose will not cluster with another dose = ",
                      "Prior standard deviation of alpha = ",
                      "Dose Limiting Toxicity probability = ",
                      "Number of MCMC iterations = ",
                      "Cut-off probability, if the 1st dose is too toxic")

    DESIGN[[1]]=X
    DESIGN[[2]]=Y
    DESIGN[[3]]=Cohort
    DESIGN[[4]]=mu_mat
    DESIGN[[5]]=p_rho
    DESIGN[[6]]=sigma
    DESIGN[[7]]=mtd
    DESIGN[[8]]=B
    DESIGN[[9]]=p_u

    TIME = paste0("Model run on: ",Sys.time())

    DESIGN=c(TIME,DESIGN)

    names(DESIGN)[[1]]="Date/Time of escalation decision"

    DESIGN = c("FBCRM Package Version: 1.0.0",DESIGN)

    ##This will provide a list with all the different outputs
    RESULTS=MFBCRM_MCMC(X,Y,mu_mat,p_rho,sigma,mtd,B)

    ##Posterior toxicity probabilities for each group
    Post_P=RESULTS[[2]]
    names(Post_P)=c(paste0("Dose",1:length(X)))

    ##Mean clustering
    CLUST=RESULTS[[3]]
    names(CLUST)=c(paste0("Dose",1:length(X)))

    ##Get optimal Dose
    OptDose= RESULTS[[4]]+1

    ##Average posterior toxicity probability at the 1st dose
    Avg_tox=RESULTS[[7]]

    ## % of time each skeletons were chosen in B iterations
    
    
    Skel_switch=rep(NA,nrow(mu_mat))
    for(j in 1:nrow(mu_mat)){
     Skel_switch[j]=mean(RESULTS[[16]]==(j-1))*100 
    }
    
    names(Skel_switch)=c(paste0("Skeleton",1:nrow(mu_mat)))

  }

  Z=as.list(rep(0,4))

 OUT1=NA
    ##If atleast 2 cohorts has been tried in 1st dose then only
    if(X[1]>Cohort){
      if(Avg_tox>p_u){

    OUT1=paste0("1st dose is too toxic, do not enroll patients at this time.")

      }
    }else{

OUT1=paste0("Next reccomended dose level for the trial : Dose ",OptDose)

    }

  Z[[1]]=OUT1

  Z[[2]]=round(t(Post_P),3)
  colnames(Z[[2]])=c(paste0("Dose",1:length(X)))

  Z[[3]]=Avg_tox

  Z[[4]]=Skel_switch

  names(Z)=c("Optimal Dose","Posterior Mean Toxicity Probability",
             "Average Toxicity Probability of 1st Dose",
             paste0("% of the times over ",B," iterations different skeletons are chosen"))

  ###Write the dataframe into the last item of the list
  Z1=as.list(c(1,2))
  Z1[[1]]=DESIGN
  Z1[[2]]=Z
  names(Z1)[[1]]="Design Parameters"
  names(Z1)[[2]]="Results"

return(Z1)

}


