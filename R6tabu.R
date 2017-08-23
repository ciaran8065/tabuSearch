library(R6)

tabuObj<-R6Class("tabuObj",
                 public=list(
                   size=NULL,
                   iters=NULL,
                   objFunc=NULL,
                   config=NULL,
                   neigh=NULL,
                   listSize=NULL,
                   nRestarts=NULL,
                   repeatAll=NULL,
                   verbose=NULL,
                   initialize=function(size = 10, iters = 100, objFunc = NULL,
                                       config = NULL,neigh = size, listSize = 9, nRestarts = 10,
                                       repeatAll = 1,verbose = FALSE){
                       self$size<-size
                       self$iters<-iters
                       self$objFunc<-objFunc
                       self$config<-config
                       self$neigh<-neigh
                       self$listSize<-listSize
                       self$nRestarts<-nRestarts
                       self$repeatAll<-repeatAll
                       self$verbose<-verbose
                   },
                   
                   ################################## BASE TABU SEARCH #####################################
                   tabuSearch=function(obj){
                     size<-obj$size
                     iters<-obj$iters
                     objFunc<-obj$objFunc
                     config<-obj$config
                     neigh<-obj$neigh
                     listSize<-obj$listSize
                     nRestarts<-obj$nRestarts
                     repeatAll<-obj$repeatAll
                     verbose<-obj$verbose
                     is.binary <- function(x){all(x %in% 0:1)}
                     if (is.null(size)){
                       if (is.null(config)){
                         stop("error: must supply either initial configuration or its length")
                       } else {
                         if(!is.binary(config)){
                           stop("error: config should be binary")  
                         }
                         size <- length(config)
                       }
                     } else {
                       if (size < 2){
                         stop("error: config too short!")
                       }
                       if (is.null(config)) {
                         config <- sample(0:1, size, replace = TRUE)  # a random binary configuration 
                       } else {
                         if (size != length(config)){
                           stop("error: length of the starting configuration != size")
                         }
                         if(!is.binary(config)){
                           stop("error: config should be binary")  
                         }
                       }
                     }
                     if (iters < 2){
                       stop("error: not enough iterations!")
                     }
                     if (is.null(objFunc)) {
                       stop("error: an evaluation function must be provided. See the objFunc parameter.")
                     } else{    
                       if (!is.function(objFunc)) {
                         stop("error: objFunc is not a valid function.")
                       }
                     }
                     
                     if (neigh > size){
                       stop("error: too many neighbours!")
                     } else{
                       if (neigh < 2){
                         stop("error: too few neighbours!")
                       }
                     }
                     if (is.null(listSize)){
                       listSize <- size - 1
                     } else {
                       if (listSize >= size){
                         stop("error: listSize too big!")
                       } else {
                         if(listSize < 1){
                           stop("error: listSize too small!")
                         }
                       }
                     } 
                     if (nRestarts < 0){
                       stop("error: nRestarts must be >= 0")
                     }
                     if (repeatAll < 0){
                       stop("error: repeatAll must be >= 0")
                     }
                     if (!is.logical(verbose)){
                       warning("warning: verbose should be logical")
                     }
                     iter <- 1 #not to be confused with "iters", iter is used to track the next row/position in the matrices/vectors at which we insert data.
                     configKeep <- matrix(, repeatAll * iters * (nRestarts + 3), size) #An empty matrix with rows specified and "size" columns to hold the configuration of binary string at each stage
                     eUtilityKeep <- vector(, repeatAll * iters * (nRestarts + 3)) #vector of undetermined type, to hold the values of objective function at each itteration
                     for (j in 1:repeatAll) { #ALGORITHM START POINT
                       if (j > 1) { #if it's not the first iteration through the algorithm we use a random configuration as the user can only specify *initial* config
                         config <- matrix(0, 1, size) #a single row matrix of 0's of length "size"
                         config[sample(1:size, sample(1:size, 1))] <- 1 #get a random number (from 1:size) of samples from the vector 1:size. Using these as indices set the digits of the config matrix to 1 at these indices.
                       }
                       
                       tabuList <- matrix(0, 1, size) #creating the tabuList as a matrix of 0's with 1 row and "size" columns, it denotes which moves are tabu
                       listOrder <- matrix(0, 1, listSize) #This single row matrix holds the order in which certain moves got added to the tabuList, it is used to implement the FIFO method.
                       eUtility <- objFunc(config) #objective function evaluated at initial condition defined by "config"
                       aspiration <- eUtility #initial aspiration value, if the objective function evaluated at a certain configuration is tabu but > this we can cancel it's tabu status
                       
                       
                       preliminarySearch <- function() {
                         configKeep[iter, ] <- config #put the current config into the matrix of used configurations in the current row denoted by "iter"
                         eUtilityKeep[iter] <- eUtility #put the value of the current config into the first position of the vector
                         
                         iter <- iter + 1 #increment the iteration, if it's the first access of this function it is set to 2 because we do 1 outside where we initially get eUtility etc.
                         for (i in 2:iters) { #see comment above to explanation of 2:iters
                           neighboursEUtility <- matrix(0, 1, size) #1 row matrix to hold the value of each of the neighbours' configurations evaluated at the objective function
                           
                           configTemp <- t(matrix(config, size, neigh)) #a matrix where each row is set to config, transpose used as matricies by default fill by column
                           
                           randomNeighbours <- sample(size, neigh) #if neigh<size it chooses the random neighbours as required but if neigh=size it takes all of them
                           
                           diag(configTemp[, randomNeighbours]) <- abs(diag(configTemp[, randomNeighbours]) - 1) #the diagonal of the matrix is inverted so that each row is different by 1 bit, this gives all the possible 1 bit-flip neighbour configurations
                           
                           neighboursEUtility[randomNeighbours] <- apply(configTemp, 1, objFunc)#"neighboursEUtility" will hold the value of the objective function calculated at each configuration
                           
                           maxNontaboo <- max(neighboursEUtility[tabuList == 0]) #finding the maximum non-tabu value
                           
                           maxTaboo <- max(neighboursEUtility[tabuList == 1], 0) #finding the maximum tabu value
                           
                           move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration,  #condition: if the move is tabu AND it's higher than the asp value
                                          
                                          ifelse(length(which(neighboursEUtility == maxTaboo)) == 1, #if true: if the maxTabu from neighbour values is unique
                                                 which(neighboursEUtility == maxTaboo), #take the value of the position of that neighbour (i.e the next place to move to)
                                                 sample(which(neighboursEUtility == maxTaboo), 1)), #if it is not unique take a random one (i.e ties can be broken arbitrarily)
                                          
                                          ifelse(length(which(neighboursEUtility == maxNontaboo & tabuList == 0)) == 1, #if false: (either maxTabu is < maxNonTabu OR maxTabu < aspiration) AND it's unique
                                                 which(neighboursEUtility == maxNontaboo & tabuList == 0), #take the value of the position of the maxNonTabu
                                                 sample(which(neighboursEUtility == maxNontaboo & tabuList == 0), 1)) #if it's not unique pick a random one
                                          
                           )
                           
                           if (eUtility >= neighboursEUtility[move]) { #if the initial(or previous) evaluated configuration is >= the value of the neighbour config we move to it but DO NOT update aspiration value
                             tabuList[move] <- 1 #set that move to tabu
                             if (sum(tabuList) > listSize) { #in this the default is tabuList in which 10 different moves can be tabu but at max there can only be 9 at one time so that we can always move
                               tabuList[listOrder[1]] <- 0 #set the first move which was made tabu to nonTabu (uses first in first out system)
                               listOrder[1:listSize] <- c(listOrder[2:listSize], 0) #discards the first element and moves all the other elements 1 to the left appending a 0 to the end. (like a left shift operation)
                             }
                             listOrder[min(which(listOrder == 0))] <- move #add the next move found earlier to the first available position in "listOrder"
                           }
                           else if (neighboursEUtility[move] > aspiration) #if the new move results in a configuration which when evaluated at the objective function is > aspiration value, set the new aspiration value to this value
                             aspiration <- neighboursEUtility[move]
                           eUtility <- neighboursEUtility[move] #set the max found value to the value found at this neighbour
                           
                           config[move] <- abs(config[move] - 1) #update the configuration by switching the bit at the position corresponding to the move
                           configKeep[iter, ] <- config #put the configuration used into "configKeep" which is storing the configurations
                           eUtilityKeep[iter] <- eUtility #put the value corresponding to this configuration into "eUtilityKeep"
                           iter <- iter + 1 #increment iter
                         }
                         result = list(aspiration = aspiration, configKeep = configKeep, eUtilityKeep = eUtilityKeep, iter = iter) #putting together the result of the preliminary search
                         return(result)
                       }
                       
                       
                       if (verbose) #if the user wants to see the steps as the algorithm executes
                         cat("Preliminary search stage...\n")
                       result <- preliminarySearch() #store result of prelim search, which is a list containing values visible in line 89
                       aspiration <- result$aspiration #set the new aspiration value to the value obtained from the preliminary search
                       configKeep <- result$configKeep #set the matrix of configurations to the configs gotten from the preliminary search function
                       eUtilityKeep <- result$eUtilityKeep #set the vector containing the values of the objective function evaluated at each of the configurations to the values obtained from the preliminary search function
                       iter <- result$iter #extract from the result list the number of iterations that occured
                       temp_asp <- 0 #temporary value for aspiration
                       restarts <- 0 #tracks the number of restarts of the intensification stage
                       while (temp_asp < aspiration & restarts < nRestarts) { #check if an equal or better solution has been found (when temp_asp >=aspiration) or we run out of restarts
                         if (verbose)
                           cat("Intensification stage...\n")
                         eUtility <- max(eUtilityKeep) #intensification stage begins with the best solution found thus far
                         temp_asp <- aspiration #assigning the value of the current aspiration to temp_asp
                         config <- configKeep[max(which(eUtilityKeep == max(eUtilityKeep))), ] #takes the configuration from "configKeep" where the position is the position of the most recent, best configuration.
                         result <- preliminarySearch() #re-run the preliminary search
                         aspiration <- result$aspiration #obtain results
                         configKeep <- result$configKeep
                         eUtilityKeep <- result$eUtilityKeep
                         iter <- result$iter
                         restarts <- restarts + 1 
                       }
                       if (verbose)
                         cat("Diversification stage...\n")
                       config <- matrix(0, 1, size) #reset the config matrix
                       config[sample(1:size, sample(1:size, 1))] <- 1 #take a random starting config
                       eUtility <- objFunc(config) #evaluate the objective function at this starting config
                       #Now finding the most frequent moves and sets them to tabu
                       frequent <- apply(configKeep, 2, function(x) sum(diff(x) != 0)) #For each column in configKeep, "sum(diff(x) !=0)" will count how many times this digit of the config changes.
                       tabuList <- as.numeric(rank(frequent, ties.method = "random") > (size - listSize)) #the values in frequent are ranked in order of magnitude, these are then compared to "(size-listSize)" by default this is 1. If frequent is > this 1, that move is set to tabu.
                       listOrder <- sample(which(rank(frequent, ties.method = "random") > (size - listSize)), listSize) #A random sample of size "listSize" is taken to be the order of the tabu list, we need an order to use the First in First out method.
                       result <- preliminarySearch() #re-run preliminary search in the hopes to find a better solution in a reasonably unexplored area of the solution space
                       iter <- result$iter #obtain results
                       configKeep <- result$configKeep
                       eUtilityKeep <- result$eUtilityKeep
                       
                     }#end of algorithm
                     
                     endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ], eUtilityKeep = eUtilityKeep[1:(iter - 1)], iters = iters, neigh = neigh, listSize = listSize, repeatAll = repeatAll)
                     class(endResult) = "tabu"
                     return(endResult)
                   },
                   
                   ###################### SUMMARY FUNCTION FOR BASE TABU #######################
                   
                   summ=function (object, verbose = FALSE, ...) #summary function
                   {
                     tabuObject <- object
                     nVars <- rowSums(tabuObject$configKeep)
                     nSelect <- colSums(tabuObject$configKeep)
                     uniqueConfig <- dim(unique(tabuObject$configKeep))[1]
                     output <- paste("Tabu Settings", "\n", "  Type                                       = ", 
                                     tabuObject$type, "\n", "  No of algorithm repeats                    = ", 
                                     tabuObject$repeatAll, "\n", "  No of iterations at each prelim search     = ", 
                                     tabuObject$iters, "\n", "  Total no of iterations                     = ", 
                                     length(nVars), "\n", "  No of unique best configurations           = ", 
                                     uniqueConfig, "\n", "  Tabu list size                             = ", 
                                     tabuObject$listSize, "\n", "  Configuration length                       = ", 
                                     length(nSelect), "\n", "  No of neighbours visited at each iteration = ", 
                                     tabuObject$neigh, "\n", sep = "")
                     maxObjFunction <- max(tabuObject$eUtilityKeep)
                     optimumNVars <- nVars[which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))]
                     optimumIteration <- which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))
                     optimumConfig <- tabuObject$configKeep[optimumIteration, 
                                                            ]
                     optionPart <- paste("Results:", "\n", "  Highest value of objective fn    = ", 
                                         round(maxObjFunction, 5), "\n", "  Occurs # of times                = ", 
                                         length(optimumNVars), "\n", "  Optimum number of variables      = ", 
                                         deparse(optimumNVars), "\n", sep = "")
                     cat(output)
                     cat(optionPart)
                     if (verbose) {
                       cat(paste("Optimum configuration:", "\n"))
                       print((optimumConfig))
                     }
                   },
                   
                   ############################# PROX TABU ############################
                   
                   tabuSearchProx=function(obj)
                   {
                     size<-obj$size
                     iters<-obj$iters
                     objFunc<-obj$objFunc
                     config<-obj$config
                     neigh<-obj$neigh
                     listSize<-obj$listSize
                     nRestarts<-obj$nRestarts
                     repeatAll<-obj$repeatAll
                     verbose<-obj$verbose
                     is.binary <- function(x){all(x %in% 0:1)}
                     if (is.null(size)){
                       if (is.null(config)){
                         stop("error: must supply either initial configuration or its length")
                       } else {
                         if(!is.binary(config)){
                           stop("error: config should be binary")  
                         }
                         size <- length(config)
                       }
                     } else {
                       if (size < 2){
                         stop("error: config too short!")
                       }
                       if (is.null(config)) {
                         config <- sample(0:1, size, replace = TRUE)  # a random binary configuration 
                       } else {
                         if (size != length(config)){
                           stop("error: length of the starting configuration != size")
                         }
                         if(!is.binary(config)){
                           stop("error: config should be binary")  
                         }
                       }
                     }
                     if (iters < 2){
                       stop("error: not enough iterations!")
                     }
                     if (is.null(objFunc)) {
                       stop("error: an evaluation function must be provided. See the objFunc parameter.")
                     } else{    
                       if (!is.function(objFunc)) {
                         stop("error: objFunc is not a valid function.")
                       }
                     }
                     
                     if (neigh > size){
                       stop("error: too many neighbours!")
                     } else{
                       if (neigh < 2){
                         stop("error: too few neighbours!")
                       }
                     }
                     if (is.null(listSize)){
                       listSize <- size - 1
                     } else {
                       if (listSize >= size){
                         stop("error: listSize too big!")
                       } else {
                         if(listSize < 1){
                           stop("error: listSize too small!")
                         }
                       }
                     } 
                     if (nRestarts < 0){
                       stop("error: nRestarts must be >= 0")
                     }
                     if (repeatAll < 0){
                       stop("error: repeatAll must be >= 0")
                     }
                     if (!is.logical(verbose)){
                       warning("warning: verbose should be logical")
                     }
                     
                     iter <- 1 #not to be confused with "iters", iter is used to track the next row/position in the matrices/vectors at which we insert data.
                     configKeep <- matrix(, repeatAll * iters * (nRestarts + 3), size) #An empty matrix with rows specified and "size" columns to hold the configuration of binary string at each stage
                     eUtilityKeep <- vector(, repeatAll * iters * (nRestarts + 3)) #vector of undetermined type, to hold the values of objective function at each itteration
                     for (j in 1:repeatAll) { #ALGORITHM START POINT
                       if (j > 1) { #if it's not the first iteration through the algorithm we use a random configuration as the user can only specify *initial* config
                         config <- matrix(0, 1, size) #a single row matrix of 0's of length "size"
                         config[sample(1:size, sample(1:size, 1))] <- 1 #get a random number (from 1:size) of samples from the vector 1:size. Using these as indices set the digits of the config matrix to 1 at these indices.
                       }
                       
                       tabuList <- matrix(0, 1, size) #creating the tabuList as a matrix of 0's with 1 row and "size" columns, it denotes which moves are tabu
                       listOrder <- matrix(0, 1, listSize) #This single row matrix holds the order in which certain moves got added to the tabuList, it is used to implement the FIFO method.
                       eUtility <- objFunc(config) #objective function evaluated at initial condition defined by "config"
                       aspiration <- eUtility #initial aspiration value, if the objective function evaluated at a certain configuration is tabu but > this we can cancel it's tabu status
                       
                       
                       preliminarySearch <- function() {
                         configKeep[iter, ] <- config #put the current config into the matrix of used configurations in the current row denoted by "iter"
                         eUtilityKeep[iter] <- eUtility #put the value of the current config into the first position of the vector
                         
                         iter <- iter + 1 #increment the iteration, if it's the first access of this function it is set to 2 because we do 1 outside where we initially get eUtility etc.
                         for (i in 2:iters) { #see comment above to explanation of 2:iters
                           neighboursEUtility <- matrix(0, 1, size) #1 row matrix to hold the value of each of the neighbours' configurations evaluated at the objective function
                           configTemp <- t(matrix(config, size, neigh)) #a matrix where each row is set to config, transpose used as matricies by default fill by column
                           randomNeighbours <- sample(size, neigh) #if neigh<size it chooses the random neighbours as required but if neigh=size it takes all of them
                           diag(configTemp[, randomNeighbours]) <- abs(diag(configTemp[, randomNeighbours]) - 1) #the diagonal of the matrix is inverted so that each row is different by 1 bit, this gives all the possible 1 bit-flip neighbour configurations
                           neighboursEUtility[randomNeighbours] <- apply(configTemp, 1, objFunc)#"neighboursEUtility" will hold the value of the objective function calculated at each configuration
                           maxNontaboo <- max(neighboursEUtility[tabuList == 0]) #finding the maximum non-tabu value
                           if(length(neighboursEUtility[tabuList==1])==0){ #if there are no tabu values yet set the tabu value to the largest possible negative integer "-2147483647". As 0 is the optimal answer from this algorithm setting maxTabu to 0 will cause issues.
                             maxTaboo<-(-2147483647)
                           }else{
                             maxTaboo <- max(neighboursEUtility[tabuList == 1]) #finding the maximum tabu value
                           }
                           
                           #determining where to move to next
                           move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration,  #condition: if the move is tabu AND it's higher than the asp value
                                          
                                          ifelse(length(which(neighboursEUtility == maxTaboo)) == 1, #if true: if the maxTabu from neighbour values is unique
                                                 which(neighboursEUtility == maxTaboo), #take the value of the position of that neighbour (i.e the next place to move to)
                                                 sample(which(neighboursEUtility == maxTaboo),1)), #if it is not unique take a random one (i.e ties can be broken arbitrarily)
                                          
                                          ifelse(length(which(neighboursEUtility == maxNontaboo & tabuList == 0)) == 1, #if false: (either maxTabu is < maxNonTabu OR maxTabu < aspiration) AND it's unique
                                                 which(neighboursEUtility == maxNontaboo & tabuList == 0), #take the value of the position of the maxNonTabu
                                                 sample(which(neighboursEUtility == maxNontaboo & tabuList == 0),1)) #if it's not unique pick a random one
                                          
                           )
                           
                           if (eUtility >= neighboursEUtility[move]) { #if the initial(or previous) evaluated configuration is >= the value of the neighbour config we move to it but DO NOT update aspiration value
                             tabuList[move] <- 1 #set that move to tabu
                             if (sum(tabuList) > listSize) { #in this the default is tabuList in which 10 different moves can be tabu but at max there can only be 9 at one time so that we can always move
                               tabuList[listOrder[1]] <- 0 #set the first move which was made tabu to nonTabu (uses first in first out system)
                               listOrder[1:listSize] <- c(listOrder[2:listSize], 0) #discards the first element and moves all the other elements 1 to the left appending a 0 to the end. (like a left shift operation)
                             }
                             listOrder[min(which(listOrder == 0))] <- move #add the next move found earlier to the first available position in "listOrder"
                           }
                           else if (neighboursEUtility[move] > aspiration) #if the new move results in a configuration which when evaluated at the objective function is > aspiration value, set the new aspiration value to this value
                             aspiration <- neighboursEUtility[move]
                           eUtility <- neighboursEUtility[move] #set the max found value to the value found at this neighbour
                           
                           config[move] <- abs(config[move] - 1) #update the configuration by switching the bit at the position corresponding to the move
                           configKeep[iter, ] <- config #put the configuration used into "configKeep" which is storing the configurations
                           eUtilityKeep[iter] <- eUtility #put the value corresponding to this configuration into "eUtilityKeep"
                           iter <- iter + 1 #increment iter
                         }
                         result = list(aspiration = aspiration, configKeep = configKeep, eUtilityKeep = eUtilityKeep, iter = iter) #putting together the result of the preliminary search
                         return(result)
                       }
                       
                       
                       if (verbose) #if the user wants to see the steps as the algorithm executes
                         cat("Preliminary search stage...\n")
                       result <- preliminarySearch() #store result of prelim search, which is a list containing values visible in line 89
                       aspiration <- result$aspiration #set the new aspiration value to the value obtained from the preliminary search
                       configKeep <- result$configKeep #set the matrix of configurations to the configs gotten from the preliminary search function
                       eUtilityKeep <- result$eUtilityKeep #set the vector containing the values of the objective function evaluated at each of the configurations to the values obtained from the preliminary search function
                       iter <- result$iter #extract from the result list the number of iterations that occured
                       temp_asp <- 0 #temporary value for aspiration
                       restarts <- 0 #tracks the number of restarts of the intensification stage
                       while (temp_asp < aspiration & restarts < nRestarts) { #check if an equal or better solution has been found (when temp_asp >=aspiration) or we run out of restarts
                         if (verbose)
                           cat("Intensification stage...\n")
                         eUtility <- max(eUtilityKeep) #intensification stage begins with the best solution found thus far
                         temp_asp <- aspiration #assigning the value of the current aspiration to temp_asp
                         config <- configKeep[max(which(eUtilityKeep == max(eUtilityKeep))), ] #takes the configuration from "configKeep" where the position is the position of the most recent, best configuration.
                         result <- preliminarySearch() #re-run the preliminary search
                         aspiration <- result$aspiration #obtain results
                         configKeep <- result$configKeep
                         eUtilityKeep <- result$eUtilityKeep
                         iter <- result$iter
                         restarts <- restarts + 1 
                       }
                       if (verbose)
                         cat("Diversification stage...\n")
                       config <- matrix(0, 1, size) #reset the config matrix
                       config[sample(1:size, sample(1:size, 1))] <- 1 #take a random starting config
                       eUtility <- objFunc(config) #evaluate the objective function at this starting config
                       #Now finding the most frequent moves and sets them to tabu
                       frequent <- apply(configKeep, 2, function(x) sum(diff(x) != 0)) #For each column in configKeep, "sum(diff(x) !=0)" will count how many times this digit of the config changes.
                       tabuList <- as.numeric(rank(frequent, ties.method = "random") > (size - listSize)) #the values in frequent are ranked in order of magnitude, these are then compared to "(size-listSize)" by default this is 1. If frequent is > this 1, that move is set to tabu.
                       listOrder <- sample(which(rank(frequent, ties.method = "random") > (size - listSize)), listSize) #A random sample of size "listSize" is taken to be the order of the tabu list, we need an order to use the First in First out method.
                       result <- preliminarySearch() #re-run preliminary search in the hopes to find a better solution in a reasonably unexplored area of the solution space
                       iter <- result$iter #obtain results
                       configKeep <- result$configKeep
                       eUtilityKeep <- result$eUtilityKeep
                       
                     }#end of algorithm
                     
                     
                     endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ], eUtilityKeep = eUtilityKeep[1:(iter - 1)], iters = iters, neigh = neigh, listSize = listSize, repeatAll = repeatAll)
                     class(endResult) = "tabu"
                     return(endResult)
                     
                   },
                   
                   ################ SUMMARY FUNCTION FOR PROXIMITY TABU #################
                   
                   ProxSumm=function (object, verbose = FALSE, ...) 
                   {
                     tabuObject <- object
                     nVars <- rowSums(tabuObject$configKeep)
                     nSelect <- colSums(tabuObject$configKeep)
                     uniqueConfig <- dim(unique(tabuObject$configKeep))[1]
                     output <- paste("Tabu Settings", "\n", "  Type                                       = ", 
                                     tabuObject$type, "\n", "  No of algorithm repeats                    = ", 
                                     tabuObject$repeatAll, "\n", "  No of iterations at each prelim search     = ", 
                                     tabuObject$iters, "\n", "  Total no of iterations                     = ", 
                                     length(nVars), "\n", "  No of unique best configurations           = ", 
                                     uniqueConfig, "\n", "  Tabu list size                             = ", 
                                     tabuObject$listSize, "\n", "  Configuration length                       = ", 
                                     length(nSelect), "\n", "  No of neighbours visited at each iteration = ", 
                                     tabuObject$neigh, "\n", sep = "")
                     maxObjFunction <- max(tabuObject$eUtilityKeep)
                     optimumNVars <- nVars[which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))]
                     optimumIteration <- which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))
                     if(length(optimumIteration)==1){
                       optimumConfig <- matrix(tabuObject$configKeep[optimumIteration, ],1,length(nSelect))
                     }else{
                       optimumConfig<-tabuObject$configKeep[optimumIteration, ]
                     }
                     optionPart <- paste("Results:", "\n", "  Minimum found distance to fn     = ", 
                                         round(maxObjFunction, 5), "\n", "  Occurs # of times                = ", 
                                         length(optimumNVars), "\n", "  Optimum number of variables      = ", 
                                         deparse(optimumNVars), "\n", sep = "")
                     cat(output)
                     cat(optionPart)
                     if (verbose) {
                       cat(paste("Optimum configuration:", "\n"))
                       print((optimumConfig))
                     }
                   },
                   
                   ####################### KNAPSACK TABU #########################
                   
                   tabuSearchKnap=function(obj,weights=NULL,values=NULL,limit=NULL){
                     size<-obj$size
                     iters<-obj$iters
                     objFunc<-obj$objFunc
                     config<-obj$config
                     neigh<-obj$neigh
                     listSize<-obj$listSize
                     nRestarts<-obj$nRestarts
                     repeatAll<-obj$repeatAll
                     verbose<-obj$verbose
                     is.binary <- function(x){all(x %in% 0:1)}
                     if (is.null(size)){
                       if (is.null(config)){
                         stop("error: must supply either initial configuration or its length")
                       } else {
                         if(!is.binary(config)){
                           stop("error: config should be binary")  
                         }
                         size <- length(config)
                       }
                     } else {
                       if (size < 2){
                         stop("error: config too short!")
                       }
                       if(is.null(limit)){
                         stop("A maximum weight must be provided")
                       }
                       if(length(values)!=length(weights)){
                         stop("Weights and Value vectors must be of equal length")
                       }
                       if(length(values)!=size){
                         stop("Size must equal the number of values")
                       }
                       if(length(weights)!=size){
                         stop("Size must equal the number of weights")
                       }
                       if (is.null(config)) { #if no config create a default.
                         greedy<-function(size,w,v,limit){
                           vals<-v
                           ws<-w
                           rs<-vals/ws
                           ind<-sort.int(rs,decreasing=T,index.return=T)$ix
                           reOrderedWeights<-(ws[ind])
                           curW<-0
                           Wm<-limit
                           conf<-numeric(size)
                           track<-1
                           while(track<=size){
                             if(curW+reOrderedWeights[track]>Wm){
                               break
                             }else{
                               curW<-curW+reOrderedWeights[track]
                               conf[ind[track]]<-1
                               track<-track+1
                             }
                           }
                           return(conf)
                         }
                         config<-greedy(size,weights,values,limit)
                       }else if (size != length(config)) {
                         stop("Length of the starting configuration != size")
                       }
                     }
                     if (iters < 2){
                       stop("error: not enough iterations!")
                     }
                     if (is.null(objFunc)) {
                       stop("error: an evaluation function must be provided. See the objFunc parameter.")
                     } else{    
                       if (!is.function(objFunc)) {
                         stop("error: objFunc is not a valid function.")
                       }
                     }
                     
                     if (neigh > size){
                       stop("error: too many neighbours!")
                     } else{
                       if (neigh < 2){
                         stop("error: too few neighbours!")
                       }
                     }
                     if (is.null(listSize)){
                       listSize <- size - 1
                     } else {
                       if (listSize >= size){
                         stop("error: listSize too big!")
                       } else {
                         if(listSize < 1){
                           stop("error: listSize too small!")
                         }
                       }
                     } 
                     if (nRestarts < 0){
                       stop("error: nRestarts must be >= 0")
                     }
                     if (repeatAll < 0){
                       stop("error: repeatAll must be >= 0")
                     }
                     if (!is.logical(verbose)){
                       warning("warning: verbose should be logical")
                     }
                     iter <- 1 #not to be confused with "iters", iter is used to track the next row/position in the matrices/vectors at which we insert data.
                     configKeep <- matrix(, repeatAll * iters * (nRestarts + 3), size) #An empty matrix with rows specified and "size" columns to hold the configuration of binary string at each stage
                     eUtilityKeep <- matrix(, repeatAll * iters * (nRestarts + 3),3) #vector of undetermined type, to hold the values of objective function at each itteration
                     for (j in 1:repeatAll) { #ALGORITHM START POINT
                       if (j > 1) { #if it's not the first iteration through the algorithm we use a random configuration as the user can only specify *initial* config
                         config <- matrix(0, 1, size) #a single row matrix of 0's of length "size"
                         config[sample(1:size, sample(1:size, 1))] <- 1 #get a random number (from 1:size) of samples from the vector 1:size. Using these as indices set the digits of the config matrix to 1 at these indices.
                       }
                       
                       tabuList <- matrix(0, 1, size) #creating the tabuList as a matrix of 0's with 1 row and "size" columns, it denotes which moves are tabu
                       listOrder <- matrix(0, 1, listSize) #This single row matrix holds the order in which certain moves got added to the tabuList, it is used to implement the FIFO method.
                       eUtility <- objFunc(config,weights,values,limit) #objective function evaluated at initial condition defined by "config"
                       aspiration <- eUtility #initial aspiration value, if the objective function evaluated at a certain configuration is tabu but > this we can cancel it's tabu status
                       
                       
                       preliminarySearch <- function() {
                         configKeep[iter, ] <- config #put the current config into the matrix of used configurations in the current row denoted by "iter"
                         eUtilityKeep[iter,] <- eUtility #put the value of the current config into the first position of the vector
                         
                         iter <- iter + 1 #increment the iteration, if it's the first access of this function it is set to 2 because we do 1 outside where we initially get eUtility etc.
                         for (i in 2:iters) { #see comment above to explanation of 2:iters
                           neighboursEUtility <- matrix(0, size, 3)#1 row matrix to hold the value of each of the neighbours' configurations evaluated at the objective function
                           configTemp <- t(matrix(config, size, neigh)) #a matrix where each row is set to config, transpose used as matricies by default fill by column
                           randomNeighbours <- sample(size, neigh) #if neigh<size it chooses the random neighbours as required but if neigh=size it takes all of them
                           diag(configTemp[, randomNeighbours]) <- abs(diag(configTemp[, randomNeighbours]) - 1) #the diagonal of the matrix is inverted so that each row is different by 1 bit, this gives all the possible 1 bit-flip neighbour configurations
                           neighboursEUtility[randomNeighbours,] <- t(apply(configTemp, 1, objFunc,weights=weights,values=values,limit=limit))#"neighboursEUtility" will hold the value of the objective function calculated at each configuration
                           maxNontaboo <- max(neighboursEUtility[tabuList == 0,1]) #finding the maximum non-tabu value
                           maxTaboo <- max(neighboursEUtility[tabuList == 1,1], 0) #finding the maximum tabu value
                           
                           #determining where to move to next
                           
                           move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration[1],  #condition: if the move is tabu AND it's higher than the asp value
                                          
                                          ifelse(length(which(neighboursEUtility[,1] == maxTaboo)) == 1, #if true: if the maxTabu from neighbour values is unique
                                                 which(neighboursEUtility[,1] == maxTaboo), #take the value of the position of that neighbour (i.e the next place to move to)
                                                 sample(which(neighboursEUtility[,1] == maxTaboo), 1)), #if it is not unique take a random one (i.e ties can be broken arbitrarily)
                                          
                                          ifelse(length(which(neighboursEUtility[,1] == maxNontaboo & tabuList == 0)) == 1, #if false: (either maxTabu is < maxNonTabu OR maxTabu < aspiration) AND it's unique
                                                 which(neighboursEUtility[,1] == maxNontaboo & tabuList == 0), #take the value of the position of the maxNonTabu
                                                 sample(which(neighboursEUtility[,1] == maxNontaboo & tabuList == 0), 1)) #if it's not unique pick a random one
                                          
                           ) #NEXT MOVE HAS BEEN CHOSEN
                           
                           if (eUtility[1] >= neighboursEUtility[move,1]) {#if the initial(or previous) evaluated configuration is >= the value of the neighbour config we move to it but DO NOT update aspiration value
                             
                             tabuList[move] <- 1 #set that move to tabu
                             if (sum(tabuList) > listSize) { #in this the default is tabuList in which 10 different moves can be tabu but at max there can only be 9 at one time so that we can always move
                               
                               tabuList[listOrder[1]] <- 0 #set the first move which was made tabu to nonTabu (uses first in first out system)
                               listOrder[1:listSize] <- c(listOrder[2:listSize], 0) #discards the first element and moves all the other elements 1 to the left appending a 0 to the end. (like a left shift operation)
                             }
                             listOrder[min(which(listOrder == 0))] <- move #add the next move found earlier to the first available position in "listOrder"
                           }
                           else if (neighboursEUtility[move,1] > aspiration[1]){ #if the new move results in a configuration which when evaluated at the objective function is > aspiration value, set the new aspiration value to this value
                            
                             aspiration <- neighboursEUtility[move,]
                           }
                           eUtility <- neighboursEUtility[move,] #set the max found value to the value found at this neighbour
                           
                           config[move] <- abs(config[move] - 1) #update the configuration by switching the bit at the position corresponding to the move
                           configKeep[iter, ] <- config #put the configuration used into "configKeep" which is storing the configurations
                           eUtilityKeep[iter,] <- eUtility #put the value corresponding to this configuration into "eUtilityKeep"
                           iter <- iter + 1 #increment iter
                         }
                         result = list(aspiration = aspiration, configKeep = configKeep, eUtilityKeep = eUtilityKeep, iter = iter) #putting together the result of the preliminary search
                         return(result)
                       }
                       
                       
                       if (verbose) #if the user wants to see the steps as the algorithm executes
                         cat("Preliminary search stage...\n")
                       result <- preliminarySearch() #store result of prelim search, which is a list containing values visible in line 89
                       aspiration <- result$aspiration[1] #set the new aspiration value to the value obtained from the preliminary search
                       configKeep <- result$configKeep #set the matrix of configurations to the configs gotten from the preliminary search function
                       eUtilityKeep <- result$eUtilityKeep #set the vector containing the values of the objective function evaluated at each of the configurations to the values obtained from the preliminary search function
                       iter <- result$iter #extract from the result list the number of iterations that occured
                       temp_asp <- 0 #temporary value for aspiration
                       restarts <- 0 #tracks the number of restarts of the intensification stage
                       
                       while (temp_asp < aspiration[1] && restarts < nRestarts) { #check if an equal or better solution has been found (when temp_asp >=aspiration) or we run out of restarts
                         if (verbose)
                           cat("Intensification stage...\n")
                         eUtility <- max(eUtilityKeep[,1],na.rm=TRUE) #intensification stage begins with the best solution found thus far
                         
                         temp_asp <- aspiration #assigning the value of the current aspiration to temp_asp
                         
                         na_rm<-as.matrix(eUtilityKeep[!is.na(eUtilityKeep)],iters,3)
                         
                         config <- configKeep[max(which(na_rm[,1] == max(na_rm[,1]))), ] #takes the configuration from "configKeep" where the position is the position of the most recent, best configuration.
                         
                         result <- preliminarySearch() #re-run the preliminary search
                         aspiration <- result$aspiration #obtain results
                         configKeep <- result$configKeep
                         eUtilityKeep <- result$eUtilityKeep
                         iter <- result$iter
                         restarts <- restarts + 1 
                       }
                       if (verbose)
                         cat("Diversification stage...\n")
                       config <- matrix(0, 1, size) #reset the config matrix
                       config[sample(1:size, sample(1:size, 1))] <- 1 #take a random starting config
                       eUtility <- objFunc(config,weights,values,limit) #evaluate the objective function at this starting config
                       #Now finding the most frequent moves and sets them to tabu
                       frequent <- apply(configKeep, 2, function(x) sum(diff(x) != 0)) #For each column in configKeep, "sum(diff(x) !=0)" will count how many times this digit of the config changes.
                       tabuList <- as.numeric(rank(frequent, ties.method = "random") > (size - listSize)) #the values in frequent are ranked in order of magnitude, these are then compared to "(size-listSize)" by default this is 1. If frequent is > this 1, that move is set to tabu.
                       listOrder <- sample(which(rank(frequent, ties.method = "random") > (size - listSize)), listSize) #A random sample of size "listSize" is taken to be the order of the tabu list, we need an order to use the First in First out method.
                       result <- preliminarySearch() #re-run preliminary search in the hopes to find a better solution in a reasonably unexplored area of the solution space
                       iter <- result$iter #obtain results
                       configKeep <- result$configKeep
                       eUtilityKeep <- result$eUtilityKeep
                      
                     }#end of algorithm
                     
                    
                     endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ], eUtilityKeep = eUtilityKeep[1:(iter - 1),], iters = iters, neigh = neigh, listSize = listSize, repeatAll = repeatAll,weights=weights,values=values,limit=limit)
                     class(endResult) = "tabu"
                     return(endResult)
                   },
                   KnapSumm=function (res,eval)
                   {
                     v<-res$eUtilityKeep[,3]==0
                     exc<-which(res$eUtilityKeep[,3]!=0) #where it is not 0, we can exclude these so that the which.max correctly matches up
                     fixed<-res$configKeep[-exc,]
                     finalConfig<-fixed[which.max((res$eUtilityKeep[res$eUtilityKeep[,3]==0,])[,1]),]
                     ans<-eval(finalConfig,res$weights,res$values,res$limit)
                     cat(" Value found: ",ans[1],"\n","Weight used: ",ans[2],"\n")
                     m<-(rbind(finalConfig))
                     rownames(m)<-"Final Configuration"
                     print(m)
                   },
                   
                   ###################### TABU TSP ########################
                   
                   tabuSearchTSP=function(obj,dist)
                     {
                       size<-obj$size
                       iters<-obj$iters
                       objFunc<-obj$objFunc
                       config<-obj$config
                       neigh<-obj$neigh
                       listSize<-obj$listSize
                       nRestarts<-obj$nRestarts
                       repeatAll<-obj$repeatAll
                       verbose<-obj$verbose
                       
                       tl<-function(v){
                         size<-length(v)
                         nums<-which(v==1)
                         combos<-t(combn(c(1:size),2))
                         tlist<-apply(combos,1,check,nums=nums)
                         return(tlist)
                       }
                       
                       check<-function(vec,nums){ #used in the tabuList function above
                         if(any(nums%in%vec)){
                           return(1)
                         }else{
                           return(0)
                         }
                       }
                       
                       
                       
                       evaluate<-function(v,d){
                         
                         pathLength<-0
                         for(i in 1:length(v)){
                           pathLength<-pathLength+(d[i,which(v==i)]) 
                         }
                         return(-1*pathLength)
                       }
                       retPath2<-function(v){
                         tempV<-v
                         config<-numeric(length(v)-1)
                         for(i in 1:(length(v)-1)){
                           t<-which(tempV==i)
                           
                           config[tempV[t-1]]<-i
                         }
                         return(config)
                       }
                       retPath<-function(v){
                         vec<-numeric(length(v))
                         for(i in 1:length(v)){
                           current<-which(v==i)
                           if((current+1)%%length(v)==0){
                             vec[i]<-v[length(v)]
                           }else{
                             vec[i]<-v[(current+1)%%length(v)]
                           }
                           current<-current+1
                         }
                         return(vec)
                       }
                       
                       swaps<-function(v){
                         m<-matrix(0,choose(length(v),2),length(v)+1)
                         track<-1
                         for(i in 1:(length(v)-1)){
                           for(j in (i+1):(length(v))){
                             
                             vect<-swap(getPath(length(v),v),i,j)
                             
                             m[track,]<-vect
                             track<-track+1
                           }
                         }
                         return(m)
                       }
                       
                       
                       findConfP<-function(v){ #parallel neighbour finding
                         if(!"foreach" %in% rownames(installed.packages())){
                           install.packages("foreach")
                         }
                         if("doParallel" %in% rownames(installed.packages())){
                           install.packages("doParallel")
                         }
                         nc<-detectCores()-1
                         cl<-makeCluster(nc)
                         library(doParallel) #put these outside
                         library(foreach)
                         m<-matrix(v,1,length(v))
                         avec<-c(1:length(v))
                         bvec<-c(1:length(v))
                         mt<-
                           foreach(b=bvec, .combine='rbind')%:%
                           foreach(a=avec, .combine='rbind', .export='Pswap') %dopar%{
                             Pswap(m,a,b)
                           }
                         stopCluster(cl)
                         return(mt)
                       }
                       
                       
                       swap<-function(v,i,j){ 
                         if(i==1){
                           v1<-NULL
                         }else{
                           v1<-v[1:(i-1)]
                         }
                         v2<-rev(v[i:j])
                         if(j+1>length(v)){
                           vec<-c(v1,v2)
                         }else{
                           if(i!=1){
                             v3<-v[(j+1):length(v)]
                             vec<-c(v1,v2,v3)
                           }else if(i==1 && j==(length(v)-1)){
                             vec<-c(v2,v2[1])
                           }else{
                             v3<-c(v[(j+1):(length(v)-1)],v2[1])
                             vec<-c(v1,v2,v3)
                           }
                         }
                         return(vec)
                       }
                       
                       generate<-function(size){
                         v<-numeric(size) 
                         t<-c(1:size) 
                         for(i in 1:size){
                           exc<-c(i,v,which(v==i)) 
                           
                           temp<-t[-exc] 
                           
                           if(i==size && (temp==size || length(temp)==0)){ 
                             v<-generate(size) 
                           }else{
                             if(length(temp)==1){
                               v[i]<-temp
                             }else{
                               v[i]<-sample(t[-exc],1) 
                             }
                           }
                           
                         }
                         return(v)
                       }
                       
                       getRandom<-function(size){
                         loop<-T
                         v<-numeric(size)
                         while(loop){
                           v<-generate(size)
                           if(valid(v)==T){
                             return(v)
                           }
                         }
                       }
                       
                       valid<-function(b){ 
                         if(length(b)<1){
                           return(F)
                         }else if(length(b)==1){
                           return(T)
                         }else{
                           track<-1
                           max<-length(b)
                           bool<-T 
                           visited<-numeric(length(b))
                           current<-1
                           for(i in 1:max){
                             visited[i]<-b[current] 
                             current<-b[current]
                             
                           }
                           
                         }
                         
                         return(!any(duplicated(visited)))
                       }
                       getPath<-function(size,resM){
                         visited<-numeric(size)
                         current<-1
                         for(i in 1:size){
                           visited[i]<-resM[current] 
                           current<-resM[current]
                         }
                         cycle<-c(1,visited)
                         return(rbind(cycle))
                       }
                       
                       getGreedy<-function(size,d){
                         diag(d)<-diag(d)+.Machine$integer.max
                         v<-numeric(size)
                         track<-1
                         exc<-numeric(size)
                         exc[1]<-1
                         for(i in 1:size){
                           exc<-c(track,exc)
                           if(i==size){
                             v[track]<-1
                           }else{
                             if(length(d[track,-exc])==1){ #it was being forced to a vector when there was length 1 and it was causing issues with "names()"
                               vect2<-which(d[track,]==d[track,-exc]) #PROBLEM this can produce more than 1 number when there is 2 numbers exactly the same in the same row
                               pos<-setdiff(vect2,exc) #This takes the differene between the 2 vectors i.e removes the excluded towns from vect2 (if an excluded town has exactly the same distance it can be included)
                               v[track]<-pos
                             }else{
                               v[track]<-as.integer(names(which.min(d[track,-exc])))
                             }
                             track<-v[track]
                             
                           }
                         }
                         return(v)
                       }
                       
                       
                       if (size < 2) { 
                         stop("error: config too short!")
                       }
                       if(is.null(dist)){
                         stop("error: must provide a distance matrix")
                       }
                       if(size!=length(dist[1,])){
                         stop("error: size must be equal to the number of variables")
                       }
                       if (iters < 2) { 
                         stop("error: not enough iterations!")
                       }
                       if (listSize >= size) { 
                         stop("error: listSize too big!")
                       }
                       if (neigh > size^2) { 
                         stop("error: too many neighbours!")
                       }
                       if (is.null(objFunc)) { 
                         stop("A evaluation function must be provided. See the objFunc parameter.")
                       }
                       if (is.null(config)) {
                         config<-getGreedy(size,dist)
                       }
                       else if (size != length(config)) {
                         stop("Length of the starting configuration != size")
                       }
                       if (repeatAll < 1) {
                         stop("error: repeatAll must be > 0")
                       }
                       
                       iter <- 1 
                       configKeep <- matrix(0, repeatAll * iters * (nRestarts + 3), size) 
                       eUtilityKeep <- vector(, repeatAll * iters * (nRestarts + 3)) 
                       for (j in 1:repeatAll) { 
                         if (j > 1) { 
                           config<-rbind(getRandom(size))
                         }
                         tabuList <- matrix(0, 1, choose(size,2)) 
                         listOrder <- matrix(0, 1, choose(size,2))
                         eUtility <- objFunc(config,dist) 
                         aspiration <- eUtility 
                         
                         
                         preliminarySearch <- function() {
                           configKeep[iter, ] <- config 
                           eUtilityKeep[iter] <- eUtility 
                           iter <- iter + 1 
                           for (i in 2:iters) { 
                             neighboursEUtility <- matrix(0, 1, choose(size,2)) 
                             configTemp <- t(matrix(config, size, neigh)) 
                             Neighbours <- c(1:choose(size,2))
                             if(size>=60){
                               configTemp<-findConfP(config)
                             }else{
                               configTemp<-swaps(config)
                               configTemp<-t(apply(configTemp,1,retPath2))
                             }
                             
                             neighboursEUtility[Neighbours] <- apply(configTemp, 1,objFunc,d=dist)
                             
                             
                             maxNontaboo <- max(neighboursEUtility[tabuList == 0])
                             if(length(neighboursEUtility[tabuList==1])==0){
                               maxTaboo<-(-2147483647)
                             }else{
                               maxTaboo <- max(neighboursEUtility[tabuList == 1]) 
                             }
                             move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration, 
                                            
                                            ifelse(length(which(neighboursEUtility == maxTaboo)) == 1, 
                                                   which(neighboursEUtility == maxTaboo), 
                                                   sample(which(neighboursEUtility == maxTaboo), 1)),
                                            
                                            ifelse(length(which(neighboursEUtility == maxNontaboo & tabuList == 0)) == 1,
                                                   which(neighboursEUtility == maxNontaboo & tabuList==0), 
                                                   sample(which(neighboursEUtility == maxNontaboo & tabuList==0), 1)) 
                                            
                             )
                             
                             if (eUtility >= neighboursEUtility[move]) { 
                               tabuList[move] <- 1 
                               if (sum(tabuList) > listSize) { 
                                 tabuList[listOrder[1]] <- 0 
                                 listOrder[1:listSize] <- c(listOrder[2:listSize], 0) 
                               }
                               listOrder[min(which(listOrder == 0))] <- move
                             }
                             else if (neighboursEUtility[move] > aspiration)
                               aspiration <- neighboursEUtility[move]
                             eUtility <- neighboursEUtility[move] 
                             
                             config<-configTemp[move,]
                             configKeep[iter, ] <- config 
                             eUtilityKeep[iter] <- eUtility 
                             iter <- iter + 1 
                           }
                           result = list(aspiration = aspiration, configKeep = configKeep, eUtilityKeep = eUtilityKeep, iter = iter)
                           return(result)
                         }
                         
                         
                         if (verbose) 
                           cat("Preliminary search stage...\n")
                         
                         result <- preliminarySearch() 
                         aspiration <- result$aspiration 
                         configKeep <- result$configKeep 
                         eUtilityKeep <- result$eUtilityKeep 
                         iter <- result$iter 
                         temp_asp <- -2147483647
                         restarts <- 0 
                         
                         
                         while (temp_asp < aspiration & restarts < nRestarts) { 
                           if (verbose)
                             cat("Intensification stage...\n")
                           eUtility <- max(eUtilityKeep[which(eUtilityKeep!=0)]) 
                           temp_asp <- aspiration 
                           config <- configKeep[max(which(eUtilityKeep == eUtility)), ] 
                           result <- preliminarySearch() 
                           aspiration <- result$aspiration 
                           configKeep <- result$configKeep
                           eUtilityKeep <- result$eUtilityKeep
                           iter <- result$iter
                           restarts <- restarts + 1 
                         }
                         if (verbose)
                           cat("Diversification stage...\n")
                         
                         config<-rbind(getRandom(size))
                         
                         eUtility <- objFunc(config,dist) 
                         
                         frequent <- apply(configKeep, 2, function(x) sum(diff(x) != 0)) #will return size number of variables showing how much each one changed
                         
                         tabuListTemp <- as.numeric(rank(frequent, ties.method = "random") > (size - listSize)) #the values in frequent are ranked in order of magnitude, these are then compared to "(size-listSize)" by default this is 1. If frequent is > this 1, that move is set to tabu.
                         tabuList<-tl(tabuListTemp)
                         
                         ones<-which(tabuList==1)
                         listOrder<-sample(ones,length(ones))
                         
                         result <- preliminarySearch() 
                         iter <- result$iter
                         configKeep <- result$configKeep
                         eUtilityKeep <- result$eUtilityKeep
                         
                       }
                       
                       endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ], eUtilityKeep = eUtilityKeep[1:(iter - 1)], iters = iters, neigh = neigh, listSize = listSize, repeatAll = repeatAll)
                       class(endResult) = "tabu"
                       return(endResult)
                       
                     },
                   TSPSumm=function (object, verbose = FALSE, ...)
                   {
                     getPath<-function(size,resM){
                       visited<-numeric(size)
                       current<-1
                       for(i in 1:size){
                         visited[i]<-resM[current] 
                         current<-resM[current]
                       }
                       cycle<-c(1,visited)
                       return(rbind(cycle))
                     }
                     tabuObject <- object
                     nVars <- rowSums(tabuObject$configKeep)
                     nSelect <- colSums(tabuObject$configKeep)
                     uniqueConfig <- dim(unique(tabuObject$configKeep))[1]
                     output <- paste("Tabu Settings", "\n", "  Type                                       = ",
                                     tabuObject$type, "\n", "  No of algorithm repeats                    = ",
                                     tabuObject$repeatAll, "\n", "  No of iterations at each prelim search     = ",
                                     tabuObject$iters, "\n", "  Total no of iterations                     = ",
                                     length(nVars), "\n", "  No of unique best configurations           = ",
                                     uniqueConfig, "\n", "  Tabu list size                             = ",
                                     tabuObject$listSize, "\n", "  Configuration length                       = ",
                                     length(nSelect), "\n", "  No of neighbours visited at each iteration = ",
                                     tabuObject$neigh, "\n", sep = "")
                     maxObjFunction <- max(tabuObject$eUtilityKeep)
                     optimumNVars <- nVars[which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))]
                     optimumIteration <- which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))
                     if(length(optimumIteration)==1){
                       optimumConfig <- matrix(tabuObject$configKeep[optimumIteration, ],1,length(nSelect))
                     }else{
                       optimumConfig<-tabuObject$configKeep[optimumIteration, ]
                     }
                     optionPart <- paste("Results:", "\n", "  Shortest path found              = ",
                                         -1*maxObjFunction, "\n", "  Occurs # of times                = ",
                                         length(optimumNVars), "\n", "  Optimum number of variables      = ",
                                         unique(optimumNVars), "\n", sep = "")
                     cat(output)
                     cat(optionPart)
                     if (verbose) {
                       cat(paste("Optimum configuration found:", "\n"))
                       print(unique(optimumConfig))
                       cat(paste("Optimum path found:","\n"))
                       print(getPath(length(nSelect),tabuObject$configKeep[which.max(tabuObject$eUtilityKeep),]))
                     }
                   }
                   
                   
                   
                 )
)

#TESTING BASE TABU

evaluate<-function(config){#c is the configuration
  target<-400 #maximum amount to spend, maybe put this as a function input
  vCoef<-c(150,100,200,60,40,360,20) #coefficients of x1..x4
  value<-sum(config*vCoef)
  if(sum(config)==0 || value>target){
    return(0)
  }
  return(value)
}
tabu<-tabuObj$new(size=7,iters=20,objFunc=evaluate,listSize=4,nRestarts=4)
res<-tabu$tabuSearch(tabu)
tabu$summ(res,T)

#TESTING PROX TABU

evaluate2<-function(config){#c is the configuration
  limit<-400 #maximum amount to spend, maybe put this as a function input
  vCoef<-c(350,200,210,150,160,195) #coefficients of x1..x4
  value<-sum(config*vCoef)
  if(sum(config)==0)return(-1*.Machine$integer.max)#smallest integer number possible
  differ<-abs(value-limit) #proximity of the current value to the optimal value
  return(-1*differ)
}

ptabu<-tabuObj$new(size=6,iters=20,objFunc=evaluate2,listSize=4,nRestarts=4)
pres<-ptabu$tabuSearchProx(ptabu)
ptabu$ProxSumm(pres,T)

#TESTING KNAPSACK TABU

vals<-c(50,40,30,60,100,150,120,70,40)
ws<-c(50,30,60,50,50,50,50,200,10)

evaluate3<-function(conf,weights,values,limit){
  Wm<-limit
  vals<-values
  ws<-weights
  cVal<-sum(conf*vals) #config value
  cW<-sum(conf*ws) #config weight
  lim<-Wm+mean(ws)
  res<-numeric(3)
  res[1]<-ifelse(cW>lim,0,cVal) #if the weight is over the limit return 0, else return the value
  res[2]<-cW
  res[3]<-ifelse(cW>Wm && cW<lim,1,0) #if the weight is greater than Wm AND less than the limit return 1 else return 0
  return(res)
}

ktabu<-tabuObj$new(size=9,iters=50,objFunc=evaluate3,listSize=4,nRestarts=10)
kres<-ktabu$tabuSearchKnap(ktabu,weights=ws,values=vals,limit=400)
ktabu$KnapSumm(kres,evaluate3)

#TESTING TSP TABU

Xs<-sample(1:100,21,replace=T) 
Ys<-sample(1:100,21,replace=T) 
Ps<-cbind(Xs,Ys) 
d<-as.matrix(eurodist)

row.names(d)<-c(1:21)
colnames(d)<-c(1:21)

evaluate4<-function(v,d){
  pathLength<-0
  for(i in 1:length(v)){
    pathLength<-pathLength+(d[i,which(v==i)]) 
  }
  return(-1*pathLength)
}

tsptabu<-tabuObj$new(size=21,iters=10,objFunc=evaluate4,listSize=15,nRestarts=10,repeatAll=1)
tspres<-tsptabu$tabuSearchTSP(tsptabu,d)
tsptabu$TSPSumm(tspres,T)

###################

