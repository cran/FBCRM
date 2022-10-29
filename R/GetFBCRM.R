#' Provides the optimal dose level closest to the mtd where the next cohort of patients should be allotted based on the data.
#' @param X Vector of patients allotted to each dose level.
#' @param Y Vector of toxicity events in each dose.
#' @param Cohort Number of patients within each cohort.
#' @param mu Prior expected toxicity probability at each dose.
#' @param p_rho Prior probability that two dose-toxicity probabilities will not cluster together.
#' @param sigma Prior standard deviation for the parameter alpha.
#' @param mtd Maximum Tolerated dose toxicity probability (pre defined).
#' @param B Number of Iterations to run for MCMC.
#' @param p_u Cut-off toxicity probability for first dose.
#' @return A list containing (1) Design parameters and prior hyperparameters used for running the trials and (2) a posterior summary of the resuls, including the next dose to assign patients to.
#' @examples
#' X=c(3, 6, 3, 3, 3, 9, 15, 6)
#' Y=c(1, 0, 1, 0, 0, 2,  4, 5)
#' Cohort=3
#' mu=seq(0.1,0.8,0.1)
#' p_rho=0.9
#' sigma = 2
#' mtd = 0.3
#' B=2000 ##Number of iterations
#' p_u=0.9
#' Z=GetFBCRM(X, Y, Cohort, mu, p_rho, sigma, mtd, B, p_u)
#' Z
#'@export

GetFBCRM=function(X,Y,Cohort,mu,p_rho,sigma,mtd,B,p_u){

  names(X)=c(paste0("Dose",1:length(X)))
  names(Y)=c(paste0("Dose",1:length(Y)))
  names(mu)=c(paste0("Dose",1:length(mu)))

  DATAHOLD = data.frame(cbind(X,Y,mu))
  colnames(DATAHOLD)=c("Subjects assigned","Toxicity events","Initial skeletal probability")

  ERRHOLD=c(length(X), length(Y), length(mu))

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
    message("Subjects assigned at each dose, toxicity events at each dose, or initial skelatl probability has incorrect dimensions")
  }else{


    ###Contains Design parameters
    DESIGN = as.list(rep(NA,9))

    names(DESIGN) = c("Subjects assigned to each dose level:",
                      "No. of toxicity events observed in each dose level:",
                      "Cohort size for dose assignment = ",
                      "Initial skeletal toxicity probability for each dose level:",
                      "Probability that a dose will not cluster with another dose = ",
                      "Prior standard deviation of alpha = ",
                      "Dose Limiting Toxicity probability = ",
                      "Number of MCMC iterations = ",
                      "Cut-off probability, if the 1st dose is too toxic")

    DESIGN[[1]]=X
    DESIGN[[2]]=Y
    DESIGN[[3]]=Cohort
    DESIGN[[4]]=mu
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
    RESULTS=FBCRM_MCMC(X,Y,mu,p_rho,sigma,mtd,B)

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

  }

  Z=as.list(c(0,0,0))

  OUT1 = NA
  

    ##If atleast 2 cohorts has been tried in 1st dose then only
    if(X[1]>Cohort){
      if(Avg_tox>p_u){

        OUT1 = paste0("1st dose is too toxic, do not enroll patients at this time.")

      }
    }else{

       OUT1=paste0("Next reccomended dose level for the trial : Dose ",OptDose)

    }

  Z[[1]]=OUT1

  Z[[2]]=round(t(Post_P),3)
  colnames(Z[[2]])=c(paste0("Dose",1:length(X)))

  Z[[3]]=Avg_tox

  names(Z)=c("Optimal Dose","Posterior Mean Toxicity Probability",
             "Average Toxicity Probability of 1st Dose")

  ###Write the dataframe into the last item of the list
  Z1=as.list(c(1,2))
  Z1[[1]]=DESIGN
  Z1[[2]]=Z
  names(Z1)[[1]]="Design Parameters"
  names(Z1)[[2]]="Results"

return(Z1)

}


