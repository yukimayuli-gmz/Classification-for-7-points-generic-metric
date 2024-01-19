#calculate 1-solutions for metric
get_1solution_by_ilp<-function(M_metric){
  n<-sqrt(length(M_metric))
  A0<-matrix(0,ncol = n*(n-1)/2,nrow = n+n*(n-1)/2)
  M_col_row_for_0<-matrix(0,ncol = n,nrow = n)
  for (i in 1:n) {
    M_col_row_for_i<-M_col_row_for_0
    M_col_row_for_i[i,]<-M_col_row_for_i[,i]<-1
    A0[i,]<-M_col_row_for_i[col(M_col_row_for_i)<row(M_col_row_for_i)]
  }
  A0[(n+1):(n+n*(n-1)/2),]<-diag(1,ncol = n*(n-1)/2,nrow = n*(n-1)/2)
  
  if(exists("gurobi")){
    model <- list()
    model$A <- A0
    model$obj <- M_metric[col(M_metric)<row(M_metric)]
    model$modelsense <- 'max'
    model$rhs <- c(seq(from=2,to=2,length.out=n),seq(from=0,to=0,length.out=n*(n-1)/2))
    model$sense <- c(rep("=",n),rep(">=",n*(n-1)/2))
    model$vtype <- "I"
    params <- list(OutputFlag=0)
    
    gurobi_1solution<-gurobi(model,params)
    return(gurobi_1solution$x)
  }
  if(exists("Rglpk_solve_LP")){
    lp_c<-M_metric[col(M_metric)<row(M_metric)]
    lp_matrix<-A0
    lp_rhs<-c(seq(from=2,to=2,length.out=n),seq(from=0,to=0,length.out=n*(n-1)/2))
    lp_dir<-c(rep("=",n),rep(">=",n*(n-1)/2))
    eg.rglpk<-Rglpk_solve_LP(obj = lp_c,
                             mat = lp_matrix,
                             rhs = lp_rhs,
                             dir = lp_dir,
                             max = TRUE,
                             types = "I")
    return(eg.rglpk$solution)
  }
}

#calculate sub_cycle by ilp solution
get_sub_cycle_in_1solution<-function(solution_result,taxa_set){
  m<-length(taxa_set)
  M_sub_cycle<-matrix(0,ncol = m,nrow = m)
  M_sub_cycle[col(M_sub_cycle)<row(M_sub_cycle)]<-solution_result
  M_sub_cycle<-M_sub_cycle+t(M_sub_cycle)
  #M_sub_cycle
  
  seq_set<-(1:m)
  cycle_list<-list()
  cycle_number<-0
  if(length(which(M_sub_cycle==1))==0){
    i<-seq_set[1]
  }
  if(length(which(M_sub_cycle==1))>1){
    i<-(which(M_sub_cycle==1)%/%m)[1]+1
  }
  while (length(seq_set)>0) {
    sub_cycle<-c()
    if(length(seq_set)<m){
      i<-seq_set[1]
    }
    while (length(which(seq_set==i))>0) {
      sub_cycle<-c(sub_cycle,taxa_set[i])
      seq_set<-seq_set[-which(seq_set==i)]
      j<-which(M_sub_cycle[i,]>0.0000001)[1]
      if((length(which(seq_set==j))==0)&(length(which(M_sub_cycle[i,]>0.0000001))==2)){
        j<-which(M_sub_cycle[i,]>0.0000001)[2]
      }
      i<-j
    }
    cycle_number<-cycle_number+1
    cycle_list[[cycle_number]]<-sub_cycle
  }
  return(cycle_list)
}
#remove 1 point "u" in (4,5,6,7), in the 6 points 1-solution, the choice for another point, output label which can be (0,1,2,3,4,5,6,7)
get_label_for_v_in_1solution_6points__when_remove_u<-function(M_metric,cycle_list_7points,remove_u){
  n<-sqrt(length(M_metric))
  if(length(cycle_list_7points)==2){
    if(remove_u==6){
      u<-cycle_list_7points[[2]][1]
      v<-cycle_list_7points[[2]][2]
    }
    if(remove_u==7){
      u<-cycle_list_7points[[2]][2]
      v<-cycle_list_7points[[2]][1]
    }
  }
  if(length(cycle_list_7points)==3){
    if(remove_u==4){
      u<-cycle_list_7points[[2]][1]
      v<-cycle_list_7points[[2]][2]
    }
    if(remove_u==5){
      u<-cycle_list_7points[[2]][2]
      v<-cycle_list_7points[[2]][1]
    }
    if(remove_u==6){
      u<-cycle_list_7points[[3]][1]
      v<-cycle_list_7points[[3]][2]
    }
    if(remove_u==7){
      u<-cycle_list_7points[[3]][2]
      v<-cycle_list_7points[[3]][1]
    }
  }
  
  M_metric_6points<-M_metric[-u,][,-u]
  ilp_1solution_result_6sub<-get_1solution_by_ilp(M_metric_6points)
  cycle_list_6points<-get_sub_cycle_in_1solution(ilp_1solution_result_6sub,(1:n)[-u])
  
  if(length(cycle_list_7points)==2){
    for (i in 1:3) {
      if(v %in% cycle_list_6points[[i]]){
        x<-cycle_list_6points[[i]][which(cycle_list_6points[[i]]!=v)]
        return(which(cycle_list_7points[[1]]==x))
      }
    }
  }
  
  if(length(cycle_list_7points)==3){
    if(length(cycle_list_6points)==2){
      return(0)
    }else{
      for (i in 1:3) {
        if(v %in% cycle_list_6points[[i]]){
          x<-cycle_list_6points[[i]][which(cycle_list_6points[[i]]!=v)]
          if(x %in% cycle_list_7points[[1]]){
            return(which(cycle_list_7points[[1]]==x))
            #print(which(cycle_list_7points[[1]]==x))
          }else if(x %in% cycle_list_7points[[2]]){
            return(which(cycle_list_7points[[2]]==x)+3)
            #print(which(cycle_list_7points[[2]]==x)+3)
          }else if(x %in% cycle_list_7points[[3]]){
            return(which(cycle_list_7points[[3]]==x)+5)
            #print(which(cycle_list_7points[[3]]==x)+5)
          }
        }
      }
    }
  }
}

#find original case and permutation1
#permutation[1]=2 is mean put column 2 of input matrix to column 1 of new matrix
get_casenumberin37_and_permutation<-function(label_vector,cycle_list_7points){
  permutation<-(1:7)
  if(length(label_vector)==2){
    if(label_vector[1]==label_vector[2]){
      #(5,5)
      casenumberin37<-2
      permute_cycle<-(((0:4)-label_vector[2]) %% 5)+1
      permutation[permute_cycle]<-cycle_list_7points[[1]]
      permutation[6:7]<-cycle_list_7points[[2]]
    }else if((abs(label_vector[1]-label_vector[2])==1)|(abs(label_vector[1]-label_vector[2])==4)){
      #(4,5)
      casenumberin37<-3
      if(((label_vector[1]-label_vector[2]) %% 5)==1){
        permute_cycle<-(((4:0)-label_vector[2]) %% 5)+1
      }else if(((label_vector[1]-label_vector[2]) %% 5)==4){
        permute_cycle<-(((0:4)-label_vector[2]) %% 5)+1
      }
      permutation[permute_cycle]<-cycle_list_7points[[1]]
      permutation[6:7]<-cycle_list_7points[[2]]
    }else if((abs(label_vector[1]-label_vector[2])==2)|(abs(label_vector[1]-label_vector[2])==3)){
      #(3,5)
      casenumberin37<-4
      if(((label_vector[1]-label_vector[2]) %% 5)==2){
        permute_cycle<-(((4:0)-label_vector[2]) %% 5)+1
      }else if(((label_vector[1]-label_vector[2]) %% 5)==3){
        permute_cycle<-(((0:4)-label_vector[2]) %% 5)+1
      }
      permutation[permute_cycle]<-cycle_list_7points[[1]]
      permutation[6:7]<-cycle_list_7points[[2]]
    }else{
      print("label vector error")
      return()
    }
  }else if(length(label_vector)==4){
    if(label_vector[1]==0){
      if(label_vector[2]==0){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            #(0,0,0,0)
            #(0,0,0,0)
            casenumberin37<-5
            permutation[1:3]<-cycle_list_7points[[1]]
            permutation[4:5]<-cycle_list_7points[[2]]
            permutation[6:7]<-cycle_list_7points[[3]]
          }else if(label_vector[4]<=3){
            #(0,0,0,3)
            #(3,0,0,0)
            casenumberin37<-6
            permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
            permutation[permute_cycle]<-cycle_list_7points[[1]]
            permutation[5:4]<-cycle_list_7points[[3]]
            permutation[6:7]<-cycle_list_7points[[2]]
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(0,0,3,0)
            #(3,0,0,0)
            casenumberin37<-6
            permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
            permutation[permute_cycle]<-cycle_list_7points[[1]]
            permutation[4:5]<-cycle_list_7points[[3]]
            permutation[6:7]<-cycle_list_7points[[2]]
          }else if(label_vector[4]<=3){
            #(0,0,3,3)
            if(label_vector[3]==label_vector[4]){
              #(3,3,0,0)
              casenumberin37<-7
              permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[3]]
              permutation[6:7]<-cycle_list_7points[[2]]
            }else if(label_vector[3]!=label_vector[4]){
              #(2,3,0,0)
              casenumberin37<-8
              if(((label_vector[3]-label_vector[4]) %% 3)==1){
                permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
              }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              }
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[3]]
              permutation[6:7]<-cycle_list_7points[[2]]
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            print("label vector error")
            return()
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else if(label_vector[2]<=3){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            #(0,3,0,0)
            #(3,0,0,0)
            casenumberin37<-6
            permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
            permutation[permute_cycle]<-cycle_list_7points[[1]]
            permutation[5:4]<-cycle_list_7points[[2]]
            permutation[6:7]<-cycle_list_7points[[3]]
          }else if(label_vector[4]<=3){
            #(0,3,0,3)
            if(label_vector[2]==label_vector[4]){
              #(3,0,3,0)
              casenumberin37<-9
              permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[2]]
              permutation[7:6]<-cycle_list_7points[[3]]
            }else if(label_vector[2]!=label_vector[4]){
              #(2,0,3,0)
              casenumberin37<-10
              if(((label_vector[2]-label_vector[4]) %% 3)==1){
                permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
              }else if(((label_vector[2]-label_vector[4]) %% 3)==2){
                permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              }
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[2]]
              permutation[7:6]<-cycle_list_7points[[3]]
            }
          }else if(label_vector[4]<=5){
            #(0,3,0,5)
            if(label_vector[4]==4){
              #(3,0,5,0)
              casenumberin37<-30
              permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[2]]
              permutation[7:6]<-cycle_list_7points[[3]]
            }else if(label_vector[4]==5){
              print("label vector error")
              return()
            }
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(0,3,3,0)
            if(label_vector[2]==label_vector[3]){
              #(3,0,3,0)
              casenumberin37<-9
              permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }else if(label_vector[2]!=label_vector[3]){
              #(2,0,3,0)
              casenumberin37<-10
              if(((label_vector[2]-label_vector[3]) %% 3)==1){
                permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
              }else if(((label_vector[2]-label_vector[3]) %% 3)==2){
                permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              }
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }
          }else if(label_vector[4]<=3){
            #(0,3,3,3)
            if(label_vector[3]==label_vector[4]){
              if(label_vector[2]==label_vector[3]){
                #(3,3,3,0)
                casenumberin37<-11
                permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }else if(label_vector[2]!=label_vector[3]){
                #(3,3,2,0)
                casenumberin37<-13
                if(((label_vector[2]-label_vector[4]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                }else if(((label_vector[2]-label_vector[4]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }
            }else if(label_vector[3]!=label_vector[4]){
              if(label_vector[2]==label_vector[3]){
                #(2,3,3,0)
                casenumberin37<-12
                if(((label_vector[4]-label_vector[3]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }else if(label_vector[2]!=label_vector[3]){
                if(label_vector[2]==label_vector[4]){
                  #(2,3,3,0)
                  casenumberin37<-12
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[7:6]<-cycle_list_7points[[2]]
                }else if(label_vector[2]!=label_vector[4]){
                  #(1,2,3,0)
                  casenumberin37<-14
                  if(((label_vector[4]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[7:6]<-cycle_list_7points[[2]]
                }
              }
            }
          }else if(label_vector[4]<=5){
            #(0,3,3,5)
            if(label_vector[4]==4){
              if(label_vector[2]==label_vector[3]){
                #(3,0,5,3)
                casenumberin37<-28
                permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }else if(label_vector[2]!=label_vector[3]){
                #(3,0,5,2)
                casenumberin37<-29
                if(((label_vector[3]-label_vector[2]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                }else if(((label_vector[3]-label_vector[2]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }
            }else if(label_vector[4]==5){
              print("label vector error")
              return()
            }
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            #(0,3,5,0)
            if(label_vector[3]==4){
              #(3,0,5,0)
              casenumberin37<-30
              permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }else if(label_vector[3]==5){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=3){
            #(0,3,5,3)
            if(label_vector[3]==4){
              if(label_vector[2]==label_vector[4]){
                #(3,0,5,3)
                casenumberin37<-28
                permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }else if(label_vector[2]!=label_vector[4]){
                #(3,0,5,2)
                casenumberin37<-29
                if(((label_vector[4]-label_vector[2]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                }else if(((label_vector[4]-label_vector[2]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }
            }else if(label_vector[3]==5){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            #(0,3,5,5)
            if((label_vector[3]==4)&(label_vector[4]==4)){
              #(3,0,5,5)
              casenumberin37<-33
              permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }else{
              print("label vector error")
              return()
            }
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else if(label_vector[2]>=6){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            #(0,7,0,3)
            if(label_vector[2]==6){
              #(3,0,5,0)
              casenumberin37<-30
              permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[3]]
              permutation[7:6]<-cycle_list_7points[[2]]
            }else if(label_vector[2]==7){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(0,7,3,0)
            if(label_vector[2]==7){
              #(3,0,5,0)
              casenumberin37<-30
              permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[3]]
              permutation[7:6]<-cycle_list_7points[[2]]
            }else if(label_vector[2]==6){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=3){
            #(0,7,3,3)
            if(label_vector[3]==label_vector[4]){
              #(3,3,5,0)
              casenumberin37<-23
              permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              if(label_vector[2]==6){
                permutation[5:4]<-cycle_list_7points[[3]]
              }else if(label_vector[2]==7){
                permutation[4:5]<-cycle_list_7points[[3]]
              }
              permutation[7:6]<-cycle_list_7points[[2]]
            }else if(label_vector[3]!=label_vector[4]){
              #(3,2,5,0)
              casenumberin37<-27
              if(label_vector[2]==6){
                if(((label_vector[3]-label_vector[4]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }else if(label_vector[2]==7){
                if(((label_vector[4]-label_vector[3]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            print("label vector error")
            return()
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else{
        print("label vector error")
        return()
      }
    }else if(label_vector[1]<=3){
      if(label_vector[2]==0){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            #(3,0,0,0)
            #(3,0,0,0)
            casenumberin37<-6
            permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
            permutation[permute_cycle]<-cycle_list_7points[[1]]
            permutation[4:5]<-cycle_list_7points[[2]]
            permutation[6:7]<-cycle_list_7points[[3]]
          }else if(label_vector[4]<=3){
            #(3,0,0,3)
            if(label_vector[1]==label_vector[4]){
              #(3,0,3,0)
              casenumberin37<-9
              permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[7:6]<-cycle_list_7points[[3]]
            }else if(label_vector[1]!=label_vector[4]){
              #(2,0,3,0)
              casenumberin37<-10
              if(((label_vector[1]-label_vector[4]) %% 3)==1){
                permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
              }else if(((label_vector[1]-label_vector[4]) %% 3)==2){
                permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              }
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[7:6]<-cycle_list_7points[[3]]
            }
          }else if(label_vector[4]<=5){
            #(3,0,0,5)
            if(label_vector[4]==5){
              #(3,0,5,0)
              casenumberin37<-30
              permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[7:6]<-cycle_list_7points[[3]]
            }else if(label_vector[4]==4){
              print("label vector error")
              return()
            }
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(3,0,3,0)
            if(label_vector[1]==label_vector[3]){
              #(3,0,3,0)
              casenumberin37<-9
              permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }else if(label_vector[1]!=label_vector[3]){
              #(2,0,3,0)
              casenumberin37<-10
              if(((label_vector[1]-label_vector[3]) %% 3)==1){
                permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
              }else if(((label_vector[1]-label_vector[3]) %% 3)==2){
                permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              }
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }
          }else if(label_vector[4]<=3){
            #(3,0,3,3)
            if(label_vector[3]==label_vector[4]){
              if(label_vector[1]==label_vector[3]){
                #(3,3,3,0)
                casenumberin37<-11
                permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }else if(label_vector[1]!=label_vector[3]){
                #(3,3,2,0)
                casenumberin37<-13
                if(((label_vector[1]-label_vector[4]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                }else if(((label_vector[1]-label_vector[4]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }
            }else if(label_vector[3]!=label_vector[4]){
              if(label_vector[1]==label_vector[3]){
                #(2,3,3,0)
                casenumberin37<-12
                if(((label_vector[4]-label_vector[3]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }else if(label_vector[1]!=label_vector[3]){
                if(label_vector[1]==label_vector[4]){
                  #(2,3,3,0)
                  casenumberin37<-12
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else if(label_vector[1]!=label_vector[4]){
                  #(1,2,3,0)
                  casenumberin37<-14
                  if(((label_vector[4]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }
              }
            }
          }else if(label_vector[4]<=5){
            #(3,0,3,5)
            if(label_vector[4]==5){
              if(label_vector[1]==label_vector[3]){
                #(3,0,5,3)
                casenumberin37<-28
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }else if(label_vector[1]!=label_vector[3]){
                #(3,0,5,2)
                casenumberin37<-29
                if(((label_vector[3]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[3]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }
            }else if(label_vector[4]==4){
              print("label vector error")
              return()
            }
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            #(3,0,5,0)
            if(label_vector[4]==5){
              #(3,0,5,0)
              casenumberin37<-30
              permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }else if(label_vector[4]==4){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=3){
            #(3,0,5,3)
            if(label_vector[3]==5){
              if(label_vector[1]==label_vector[4]){
                #(3,0,5,3)
                casenumberin37<-28
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }else if(label_vector[1]!=label_vector[4]){
                #(3,0,5,2)
                casenumberin37<-29
                if(((label_vector[4]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[4]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }
            }else if(label_vector[3]==4){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            #(3,0,5,5)
            if((label_vector[3]==5)&(label_vector[4]==5)){
              #(3,0,5,5)
              casenumberin37<-33
              permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }else{
              print("label vector error")
              return()
            }
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else if(label_vector[2]<=3){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            #(3,3,0,0)
            if(label_vector[1]==label_vector[2]){
              #(3,3,0,0)
              casenumberin37<-7
              permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }else if(label_vector[1]!=label_vector[2]){
              #(2,3,0,0)
              casenumberin37<-8
              if(((label_vector[1]-label_vector[2]) %% 3)==1){
                permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
              }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
              }
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[2]]
              permutation[6:7]<-cycle_list_7points[[3]]
            }
          }else if(label_vector[4]<=3){
            #(3,3,0,3)
            if(label_vector[1]==label_vector[2]){
              if(label_vector[1]==label_vector[4]){
                #(3,3,3,0)
                casenumberin37<-11
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }else if(label_vector[1]!=label_vector[4]){
                #(3,3,2,0)
                casenumberin37<-13
                if(((label_vector[4]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[4]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }
            }else if(label_vector[1]!=label_vector[2]){
              if(label_vector[2]==label_vector[4]){
                #(2,3,3,0)
                casenumberin37<-12
                if(((label_vector[1]-label_vector[2]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }else if(label_vector[2]!=label_vector[4]){
                if(label_vector[1]==label_vector[4]){
                  #(2,3,3,0)
                  casenumberin37<-12
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[7:6]<-cycle_list_7points[[3]]
                }else if(label_vector[1]!=label_vector[4]){
                  #(1,2,3,0)
                  casenumberin37<-14
                  if(((label_vector[2]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[7:6]<-cycle_list_7points[[3]]
                }
              }
            }
          }else if(label_vector[4]<=5){
            #(3,3,0,5)
            if(label_vector[1]==label_vector[2]){
              #(3,3,5,0)
              casenumberin37<-23
              permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              if(label_vector[4]==4){
                permutation[5:4]<-cycle_list_7points[[2]]
              }else if(label_vector[4]==5){
                permutation[4:5]<-cycle_list_7points[[2]]
              }
              permutation[7:6]<-cycle_list_7points[[3]]
            }else if(label_vector[1]!=label_vector[2]){
              #(3,2,5,0)
              casenumberin37<-27
              if(label_vector[4]==4){
                if(((label_vector[1]-label_vector[2]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }else if(label_vector[4]==5){
                if(((label_vector[2]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }
            }
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(3,3,3,0)
            if(label_vector[1]==label_vector[2]){
              if(label_vector[1]==label_vector[3]){
                #(3,3,3,0)
                casenumberin37<-11
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }else if(label_vector[1]!=label_vector[3]){
                #(3,3,2,0)
                casenumberin37<-13
                if(((label_vector[3]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[3]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }
            }else if(label_vector[1]!=label_vector[2]){
              if(label_vector[2]==label_vector[3]){
                #(2,3,3,0)
                casenumberin37<-12
                if(((label_vector[1]-label_vector[2]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }else if(label_vector[2]!=label_vector[3]){
                if(label_vector[1]==label_vector[3]){
                  #(2,3,3,0)
                  casenumberin37<-12
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else if(label_vector[1]!=label_vector[3]){
                  #(1,2,3,0)
                  casenumberin37<-14
                  if(((label_vector[2]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }
            }
          }else if(label_vector[4]<=3){
            #(3,3,3,3)
            if(label_vector[1]==label_vector[2]){
              if(label_vector[3]==label_vector[4]){
                if(label_vector[1]==label_vector[3]){
                  #(3,3,3,3)
                  casenumberin37<-15
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else{
                  #(2,2,3,3)
                  casenumberin37<-17
                  if(((label_vector[1]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }else if(label_vector[3]!=label_vector[4]){
                if(label_vector[3]==label_vector[1]){
                  #(2,3,3,3)
                  casenumberin37<-16
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else if(label_vector[4]==label_vector[1]){
                  #(2,3,3,3)
                  casenumberin37<-16
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else{
                  #(1,2,3,3)
                  casenumberin37<-19
                  if(((label_vector[4]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }
              }
            }else if(label_vector[1]!=label_vector[2]){
              if(label_vector[3]==label_vector[4]){
                if(label_vector[1]==label_vector[3]){
                  #(2,3,3,3)
                  casenumberin37<-16
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else if(label_vector[2]==label_vector[3]){
                  #(2,3,3,3)
                  casenumberin37<-16
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else{
                  #(1,2,3,3)
                  casenumberin37<-19
                  if(((label_vector[2]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }else if(label_vector[3]!=label_vector[4]){
                if(label_vector[1]==label_vector[3]){
                  if(label_vector[2]==label_vector[4]){
                    #(2,3,2,3)
                    casenumberin37<-18
                    if(((label_vector[2]-label_vector[1]) %% 3)==1){
                      permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                    }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                      permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                    }
                    permutation[permute_cycle]<-cycle_list_7points[[1]]
                    permutation[5:4]<-cycle_list_7points[[2]]
                    permutation[7:6]<-cycle_list_7points[[3]]
                  }else if(label_vector[2]!=label_vector[4]){
                    #(1,3,2,3)
                    casenumberin37<-20
                    if(((label_vector[4]-label_vector[3]) %% 3)==1){
                      permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                    }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                      permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                    }
                    permutation[permute_cycle]<-cycle_list_7points[[1]]
                    permutation[5:4]<-cycle_list_7points[[2]]
                    permutation[7:6]<-cycle_list_7points[[3]]
                  }
                }else if(label_vector[1]==label_vector[4]){
                  if(label_vector[2]==label_vector[3]){
                    #(2,3,2,3)
                    casenumberin37<-18
                    if(((label_vector[2]-label_vector[1]) %% 3)==1){
                      permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                    }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                      permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                    }
                    permutation[permute_cycle]<-cycle_list_7points[[1]]
                    permutation[5:4]<-cycle_list_7points[[2]]
                    permutation[6:7]<-cycle_list_7points[[3]]
                  }else if(label_vector[2]!=label_vector[3]){
                    #(1,3,2,3)
                    casenumberin37<-20
                    if(((label_vector[3]-label_vector[4]) %% 3)==1){
                      permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                    }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                      permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                    }
                    permutation[permute_cycle]<-cycle_list_7points[[1]]
                    permutation[5:4]<-cycle_list_7points[[2]]
                    permutation[6:7]<-cycle_list_7points[[3]]
                  }
                }else{
                  if(label_vector[2]==label_vector[3]){
                    #(1,3,2,3)
                    casenumberin37<-20
                    if(((label_vector[4]-label_vector[3]) %% 3)==1){
                      permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                    }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                      permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                    }
                    permutation[permute_cycle]<-cycle_list_7points[[1]]
                    permutation[4:5]<-cycle_list_7points[[2]]
                    permutation[7:6]<-cycle_list_7points[[3]]
                  }else if(label_vector[2]==label_vector[4]){
                    #(1,3,2,3)
                    casenumberin37<-20
                    if(((label_vector[3]-label_vector[4]) %% 3)==1){
                      permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                    }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                      permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                    }
                    permutation[permute_cycle]<-cycle_list_7points[[1]]
                    permutation[4:5]<-cycle_list_7points[[2]]
                    permutation[6:7]<-cycle_list_7points[[3]]
                  }else{
                    print("label vector error")
                    return()
                  }
                }
              }
            }
          }else if(label_vector[4]<=5){
            #(3,3,3,5)
            if(label_vector[1]==label_vector[2]){
              if(label_vector[1]==label_vector[3]){
                #(3,3,5,3)
                casenumberin37<-21
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                if(label_vector[4]==4){
                  permutation[5:4]<-cycle_list_7points[[2]]
                }else if(label_vector[4]==5){
                  permutation[4:5]<-cycle_list_7points[[2]]
                }
                permutation[7:6]<-cycle_list_7points[[3]]
              }else if(label_vector[1]!=label_vector[3]){
                #(3,3,5,2)
                casenumberin37<-22
                if(((label_vector[3]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[3]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                if(label_vector[4]==4){
                  permutation[5:4]<-cycle_list_7points[[2]]
                }else if(label_vector[4]==5){
                  permutation[4:5]<-cycle_list_7points[[2]]
                }
                permutation[7:6]<-cycle_list_7points[[3]]
              }
            }else if(label_vector[1]!=label_vector[2]){
              if(label_vector[4]==4){
                if(label_vector[3]==label_vector[1]){
                  #(3,2,5,2)
                  casenumberin37<-25
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[7:6]<-cycle_list_7points[[3]]
                }else if(label_vector[3]==label_vector[2]){
                  #(3,2,5,3)
                  casenumberin37<-24
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[7:6]<-cycle_list_7points[[3]]
                }else{
                  #(3,2,5,1)
                  casenumberin37<-26
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[7:6]<-cycle_list_7points[[3]]
                }
              }else if(label_vector[4]==5){
                if(label_vector[3]==label_vector[1]){
                  #(3,2,5,3)
                  casenumberin37<-24
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[7:6]<-cycle_list_7points[[3]]
                }else if(label_vector[3]==label_vector[2]){
                  #(3,2,5,2)
                  casenumberin37<-25
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[7:6]<-cycle_list_7points[[3]]
                }else{
                  #(3,2,5,1)
                  casenumberin37<-26
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[7:6]<-cycle_list_7points[[3]]
                }
              }
            }
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            #(3,3,5,0)
            if(label_vector[1]==label_vector[2]){
              #(3,3,5,0)
              casenumberin37<-23
              permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              if(label_vector[3]==4){
                permutation[5:4]<-cycle_list_7points[[2]]
              }else if(label_vector[3]==5){
                permutation[4:5]<-cycle_list_7points[[2]]
              }
              permutation[6:7]<-cycle_list_7points[[3]]
            }else if(label_vector[1]!=label_vector[2]){
              #(3,2,5,0)
              casenumberin37<-27
              if(label_vector[3]==4){
                if(((label_vector[1]-label_vector[2]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }else if(label_vector[3]==5){
                if(((label_vector[2]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }
            }
          }else if(label_vector[4]<=3){
            #(3,3,5,3)
            if(label_vector[1]==label_vector[2]){
              if(label_vector[1]==label_vector[4]){
                #(3,3,5,3)
                casenumberin37<-21
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                if(label_vector[3]==4){
                  permutation[5:4]<-cycle_list_7points[[2]]
                }else if(label_vector[3]==5){
                  permutation[4:5]<-cycle_list_7points[[2]]
                }
                permutation[6:7]<-cycle_list_7points[[3]]
              }else if(label_vector[1]!=label_vector[4]){
                #(3,3,5,2)
                casenumberin37<-22
                if(((label_vector[4]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[4]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                if(label_vector[3]==4){
                  permutation[5:4]<-cycle_list_7points[[2]]
                }else if(label_vector[3]==5){
                  permutation[4:5]<-cycle_list_7points[[2]]
                }
                permutation[6:7]<-cycle_list_7points[[3]]
              }
            }else if(label_vector[1]!=label_vector[2]){
              if(label_vector[3]==4){
                if(label_vector[4]==label_vector[1]){
                  #(3,2,5,2)
                  casenumberin37<-25
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else if(label_vector[4]==label_vector[2]){
                  #(3,2,5,3)
                  casenumberin37<-24
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else{
                  #(3,2,5,1)
                  casenumberin37<-26
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }else if(label_vector[3]==5){
                if(label_vector[4]==label_vector[1]){
                  #(3,2,5,3)
                  casenumberin37<-24
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else if(label_vector[4]==label_vector[2]){
                  #(3,2,5,2)
                  casenumberin37<-25
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else{
                  #(3,2,5,1)
                  casenumberin37<-26
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }
            }
          }else if(label_vector[4]<=5){
            #(3,3,5,5)
            if(label_vector[3]==4){
              if(label_vector[4]==4){
                if(label_vector[1]==label_vector[2]){
                  #(3,3,5,5)
                  casenumberin37<-31
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else if(label_vector[1]!=label_vector[2]){
                  #(3,2,5,5)
                  casenumberin37<-32
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }else if(label_vector[4]==5){
                if(label_vector[1]==label_vector[2]){
                  #(3,3,5,4)
                  casenumberin37<-34
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else if(label_vector[1]!=label_vector[2]){
                  #(3,2,5,4)
                  casenumberin37<-35
                  if(((label_vector[1]-label_vector[2]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                  }else if(((label_vector[1]-label_vector[2]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }
            }else if(label_vector[3]==5){
              if(label_vector[4]==4){
                if(label_vector[1]==label_vector[2]){
                  #(3,3,5,4)
                  casenumberin37<-34
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else if(label_vector[1]!=label_vector[2]){
                  #(3,2,5,4)
                  casenumberin37<-35
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }else if(label_vector[4]==5){
                if(label_vector[1]==label_vector[2]){
                  #(3,3,5,5)
                  casenumberin37<-31
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }else if(label_vector[1]!=label_vector[2]){
                  #(3,2,5,5)
                  casenumberin37<-32
                  if(((label_vector[2]-label_vector[1]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                  }else if(((label_vector[2]-label_vector[1]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[2]]
                  permutation[6:7]<-cycle_list_7points[[3]]
                }
              }
            }
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else if(label_vector[2]>=6){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            #(3,7,0,3)
            if(label_vector[2]==6){
              if(label_vector[1]==label_vector[4]){
                #(3,0,5,3)
                casenumberin37<-28
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }else if(label_vector[1]!=label_vector[4]){
                #(3,0,5,2)
                casenumberin37<-29
                if(((label_vector[1]-label_vector[4]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                }else if(((label_vector[1]-label_vector[4]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }
            }else if(label_vector[2]==7){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(3,7,3,0)
            if(label_vector[2]==7){
              if(label_vector[1]==label_vector[3]){
                #(3,0,5,3)
                casenumberin37<-28
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }else if(label_vector[1]!=label_vector[3]){
                #(3,0,5,2)
                casenumberin37<-29
                if(((label_vector[1]-label_vector[3]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                }else if(((label_vector[1]-label_vector[3]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[7:6]<-cycle_list_7points[[2]]
              }
            }else if(label_vector[2]==6){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=3){
            #(3,7,3,3)
            if(label_vector[3]==label_vector[4]){
              if(label_vector[1]==label_vector[3]){
                #(3,3,5,3)
                casenumberin37<-21
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                if(label_vector[2]==6){
                  permutation[5:4]<-cycle_list_7points[[3]]
                }else if(label_vector[2]==7){
                  permutation[4:5]<-cycle_list_7points[[3]]
                }
                permutation[7:6]<-cycle_list_7points[[2]]
              }else if(label_vector[1]!=label_vector[3]){
                #(3,3,5,2)
                casenumberin37<-22
                if(((label_vector[1]-label_vector[3]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                }else if(((label_vector[1]-label_vector[3]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                if(label_vector[2]==6){
                  permutation[5:4]<-cycle_list_7points[[3]]
                }else if(label_vector[2]==7){
                  permutation[4:5]<-cycle_list_7points[[3]]
                }
                permutation[7:6]<-cycle_list_7points[[2]]
              }
            }else if(label_vector[3]!=label_vector[4]){
              if(label_vector[2]==6){
                if(label_vector[1]==label_vector[3]){
                  #(3,2,5,2)
                  casenumberin37<-25
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[7:6]<-cycle_list_7points[[2]]
                }else if(label_vector[1]==label_vector[4]){
                  #(3,2,5,3)
                  casenumberin37<-24
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[7:6]<-cycle_list_7points[[2]]
                }else{
                  #(3,2,5,1)
                  casenumberin37<-26
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[7:6]<-cycle_list_7points[[2]]
                }
              }else if(label_vector[2]==7){
                if(label_vector[1]==label_vector[3]){
                  #(3,2,5,3)
                  casenumberin37<-24
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[7:6]<-cycle_list_7points[[2]]
                }else if(label_vector[1]==label_vector[4]){
                  #(3,2,5,2)
                  casenumberin37<-25
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[7:6]<-cycle_list_7points[[2]]
                }else{
                  #(3,2,5,1)
                  casenumberin37<-26
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[7:6]<-cycle_list_7points[[2]]
                }
              }
            }
          }else if(label_vector[4]<=5){
            #(3,7,3,5)
            if((label_vector[2]==7)&(label_vector[4]==5)){
              if(label_vector[1]==label_vector[3]){
                #(3,6,5,3)
                casenumberin37<-36
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }else if(label_vector[1]!=label_vector[3]){
                #(3,6,5,2)
                casenumberin37<-37
                if(((label_vector[3]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[3]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }
            }else{
              print("label vector error")
              return()
            }
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            #(3,7,5,3)
            if((label_vector[2]==6)&(label_vector[3]==5)){
              if(label_vector[1]==label_vector[4]){
                #(3,6,5,3)
                casenumberin37<-36
                permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }else if(label_vector[1]!=label_vector[4]){
                #(3,6,5,2)
                casenumberin37<-37
                if(((label_vector[4]-label_vector[1]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[1]) %% 3)+1
                }else if(((label_vector[4]-label_vector[1]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[1]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }
            }else{
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else{
        print("label vector error")
        return()
      }
    }else if(label_vector[1]>=6){
      if(label_vector[2]==0){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            #(7,0,0,3)
            if(label_vector[1]==6){
              #(3,0,5,0)
              casenumberin37<-30
              permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[3]]
              permutation[6:7]<-cycle_list_7points[[2]]
            }else if(label_vector[1]==7){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(7,0,3,0)
            if(label_vector[1]==7){
              #(3,0,5,0)
              casenumberin37<-30
              permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[3]]
              permutation[6:7]<-cycle_list_7points[[2]]
            }else if(label_vector[1]==6){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=3){
            #(7,0,3,3)
            if(label_vector[3]==label_vector[4]){
              #(3,3,5,0)
              casenumberin37<-23
              permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              if(label_vector[1]==6){
                permutation[5:4]<-cycle_list_7points[[3]]
              }else if(label_vector[1]==7){
                permutation[4:5]<-cycle_list_7points[[3]]
              }
              permutation[6:7]<-cycle_list_7points[[2]]
            }else if(label_vector[3]!=label_vector[4]){
              #(3,2,5,0)
              casenumberin37<-27
              if(label_vector[1]==6){
                if(((label_vector[3]-label_vector[4]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                }else if(((label_vector[4]-label_vector[4]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }else if(label_vector[1]==7){
                if(((label_vector[4]-label_vector[3]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            print("label vector error")
            return()
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else if(label_vector[2]<=3){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            #(7,3,0,3)
            if(label_vector[1]==6){
              if(label_vector[2]==label_vector[4]){
                #(3,0,5,3)
                casenumberin37<-28
                permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }else if(label_vector[2]!=label_vector[4]){
                #(3,0,5,2)
                casenumberin37<-29
                if(((label_vector[2]-label_vector[4]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                }else if(((label_vector[2]-label_vector[4]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }
            }else if(label_vector[1]==7){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(7,3,3,0)
            if(label_vector[1]==7){
              if(label_vector[2]==label_vector[3]){
                #(3,0,5,3)
                casenumberin37<-28
                permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }else if(label_vector[2]!=label_vector[3]){
                #(3,0,5,2)
                casenumberin37<-29
                if(((label_vector[2]-label_vector[3]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                }else if(((label_vector[2]-label_vector[3]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[4:5]<-cycle_list_7points[[3]]
                permutation[6:7]<-cycle_list_7points[[2]]
              }
            }else if(label_vector[1]==6){
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=3){
            #(7,3,3,3)
            if(label_vector[3]==label_vector[4]){
              if(label_vector[2]==label_vector[3]){
                #(3,3,5,3)
                casenumberin37<-21
                permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                if(label_vector[1]==6){
                  permutation[5:4]<-cycle_list_7points[[3]]
                }else if(label_vector[1]==7){
                  permutation[4:5]<-cycle_list_7points[[3]]
                }
                permutation[6:7]<-cycle_list_7points[[2]]
              }else if(label_vector[2]!=label_vector[3]){
                #(3,3,5,2)
                casenumberin37<-22
                if(((label_vector[2]-label_vector[3]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                }else if(((label_vector[2]-label_vector[3]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                if(label_vector[1]==6){
                  permutation[5:4]<-cycle_list_7points[[3]]
                }else if(label_vector[1]==7){
                  permutation[4:5]<-cycle_list_7points[[3]]
                }
                permutation[6:7]<-cycle_list_7points[[2]]
              }
            }else if(label_vector[3]!=label_vector[4]){
              if(label_vector[1]==6){
                if(label_vector[2]==label_vector[3]){
                  #(3,2,5,2)
                  casenumberin37<-25
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else if(label_vector[2]==label_vector[4]){
                  #(3,2,5,3)
                  casenumberin37<-24
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else{
                  #(3,2,5,1)
                  casenumberin37<-26
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }
              }else if(label_vector[1]==7){
                if(label_vector[2]==label_vector[3]){
                  #(3,2,5,3)
                  casenumberin37<-24
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else if(label_vector[2]==label_vector[4]){
                  #(3,2,5,2)
                  casenumberin37<-25
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else{
                  #(3,2,5,1)
                  casenumberin37<-26
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }
              }
            }
          }else if(label_vector[4]<=5){
            #(7,3,3,5)
            if((label_vector[1]==7)&(label_vector[4]==4)){
              if(label_vector[2]==label_vector[3]){
                #(3,6,5,3)
                casenumberin37<-36
                permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }else if(label_vector[2]!=label_vector[3]){
                #(3,6,5,2)
                casenumberin37<-37
                if(((label_vector[3]-label_vector[2]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                }else if(((label_vector[3]-label_vector[2]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[7:6]<-cycle_list_7points[[3]]
              }
            }else{
              print("label vector error")
              return()
            }
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            #(7,3,5,3)
            if((label_vector[1]==6)&(label_vector[3]==4)){
              if(label_vector[2]==label_vector[4]){
                #(3,6,5,3)
                casenumberin37<-36
                permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }else if(label_vector[2]!=label_vector[4]){
                #(3,6,5,2)
                casenumberin37<-37
                if(((label_vector[4]-label_vector[2]) %% 3)==1){
                  permute_cycle<-(((2:0)-label_vector[2]) %% 3)+1
                }else if(((label_vector[4]-label_vector[2]) %% 3)==2){
                  permute_cycle<-(((0:2)-label_vector[2]) %% 3)+1
                }
                permutation[permute_cycle]<-cycle_list_7points[[1]]
                permutation[5:4]<-cycle_list_7points[[2]]
                permutation[6:7]<-cycle_list_7points[[3]]
              }
            }else{
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else if(label_vector[2]>=6){
        if(label_vector[3]==0){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            #(7,7,0,3)
            if((label_vector[1]==6)&(label_vector[2]==6)){
              #(3,0,5,5)
              casenumberin37<-33
              permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[5:4]<-cycle_list_7points[[3]]
              permutation[6:7]<-cycle_list_7points[[2]]
            }else{
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=3){
          if(label_vector[4]==0){
            #(7,7,3,0)
            if((label_vector[1]==7)&(label_vector[2]==7)){
              #(3,0,5,5)
              casenumberin37<-33
              permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
              permutation[permute_cycle]<-cycle_list_7points[[1]]
              permutation[4:5]<-cycle_list_7points[[3]]
              permutation[6:7]<-cycle_list_7points[[2]]
            }else{
              print("label vector error")
              return()
            }
          }else if(label_vector[4]<=3){
            #(7,7,3,3)
            if(label_vector[1]==6){
              if(label_vector[2]==6){
                if(label_vector[3]==label_vector[4]){
                  #(3,3,5,5)
                  casenumberin37<-31
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else if(label_vector[3]!=label_vector[4]){
                  #(3,2,5,5)
                  casenumberin37<-32
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }
              }else if(label_vector[2]==7){
                if(label_vector[3]==label_vector[4]){
                  #(3,3,5,4)
                  casenumberin37<-34
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else if(label_vector[3]!=label_vector[4]){
                  #(3,2,5,4)
                  casenumberin37<-35
                  if(((label_vector[3]-label_vector[4]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[4]) %% 3)+1
                  }else if(((label_vector[3]-label_vector[4]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[4]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[5:4]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }
              }
            }else if(label_vector[1]==7){
              if(label_vector[2]==6){
                if(label_vector[3]==label_vector[4]){
                  #(3,3,5,4)
                  casenumberin37<-34
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else if(label_vector[3]!=label_vector[4]){
                  #(3,2,5,4)
                  casenumberin37<-35
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }
              }else if(label_vector[2]==7){
                if(label_vector[3]==label_vector[4]){
                  #(3,3,5,5)
                  casenumberin37<-31
                  permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }else if(label_vector[3]!=label_vector[4]){
                  #(3,2,5,5)
                  casenumberin37<-32
                  if(((label_vector[4]-label_vector[3]) %% 3)==1){
                    permute_cycle<-(((2:0)-label_vector[3]) %% 3)+1
                  }else if(((label_vector[4]-label_vector[3]) %% 3)==2){
                    permute_cycle<-(((0:2)-label_vector[3]) %% 3)+1
                  }
                  permutation[permute_cycle]<-cycle_list_7points[[1]]
                  permutation[4:5]<-cycle_list_7points[[3]]
                  permutation[6:7]<-cycle_list_7points[[2]]
                }
              }
            }
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else if(label_vector[3]<=5){
          if(label_vector[4]==0){
            print("label vector error")
            return()
          }else if(label_vector[4]<=3){
            print("label vector error")
            return()
          }else if(label_vector[4]<=5){
            print("label vector error")
            return()
          }else{
            print("label vector error")
            return()
          }
        }else{
          print("label vector error")
          return()
        }
      }else{
        print("label vector error")
        return()
      }
    }else{
      print("label vector error")
      return()
    }
  }else{
    print("label vector error")
    return()
  }
  return(list(casenumberin37,permutation))
}

#permute a matrix
permute_matrix<-function(M,permutation){
  n<-length(permutation)
  M<-M[permutation,]
  M<-M[,permutation]
  return(M)
}

#check metric and Gromov structure
check_one_gromov_structure_for_metric<-function(G_matrix,pM_metric2,fesible_taxa_list){
  for (ai in 1:7) {
    for (i in 1:(length(fesible_taxa_list[[ai]])-1)) {
      for (j in (i+1):length(fesible_taxa_list[[ai]])) {
        if((pM_metric2[ai,G_matrix[ai,2]]+pM_metric2[ai,G_matrix[ai,3]]-pM_metric2[G_matrix[ai,2],G_matrix[ai,3]])>
           (pM_metric2[ai,fesible_taxa_list[[ai]][i]]+pM_metric2[ai,fesible_taxa_list[[ai]][j]]-
            pM_metric2[fesible_taxa_list[[ai]][i],fesible_taxa_list[[ai]][j]])){
          return(0)
        }
      }
    }
  }
  return(1)
}

#find right Gromov structure, if no result, output label_i=0
search_gromov_structure<-function(G_matrix_isomorph_list,pM_metric2,fesible_taxa_list){
  label_i<-0
  for (i in 1:length(G_matrix_isomorph_list)) {
    if(check_one_gromov_structure_for_metric(G_matrix_isomorph_list[[i]],pM_metric2,fesible_taxa_list)){
      label_i<-i
      break
    }
  }
  return(label_i)
}

#check metric and extended graph
check_one_extended_graph_for_metric<-function(f_matrix,pM_metric2,fesible_taxa_list){
  for (ai in 1:7) {
    if(length(fesible_taxa_list[[ai]])==4){
      if(f_matrix[ai,4]==8){
        v_extended_graph<-pM_metric2[f_matrix[ai,2],f_matrix[ai,3]]+pM_metric2[f_matrix[ai,5],f_matrix[ai,6]]-
          pM_metric2[ai,f_matrix[ai,5]]-pM_metric2[ai,f_matrix[ai,6]]
      }else if(f_matrix[ai,4]==9){
        v_extended_graph<-pM_metric2[f_matrix[ai,2],f_matrix[ai,5]]+pM_metric2[f_matrix[ai,3],f_matrix[ai,6]]-
          pM_metric2[ai,f_matrix[ai,5]]-pM_metric2[ai,f_matrix[ai,6]]
      }
      b<-setdiff(fesible_taxa_list[[ai]],f_matrix[ai,c(2,3)])
      if(v_extended_graph<pM_metric2[f_matrix[ai,2],f_matrix[ai,3]]+pM_metric2[b[1],b[2]]-
         pM_metric2[ai,b[1]]-pM_metric2[ai,b[2]]){
        return(0)
      }
      for (i in 1:2) {
        for (j in 1:2) {
          if(v_extended_graph<pM_metric2[f_matrix[ai,2],b[i]]+pM_metric2[f_matrix[ai,3],b[j]]-
             pM_metric2[ai,b[i]]-pM_metric2[ai,b[j]]){
            return(0)
          }
        }
      }
    }
  }
  return(1)
}

#find right extended graph, if no result, output label_j=0
search_extended_graph<-function(isomorph_f_matrix_list_i,pM_metric2,fesible_taxa_list){
  label_j<-0
  for (i in 1:length(isomorph_f_matrix_list_i)) {
    if(check_one_extended_graph_for_metric(isomorph_f_matrix_list_i[[i]],pM_metric2,fesible_taxa_list)){
      label_j<-i
      break
    }
  }
  return(label_j)
}

#check metric and quartets
check_one_quartets_for_metric<-function(quartets_matrix,pM_metric2){
  for (q in 1:nrow(quartets_matrix)) {
    if(quartets_matrix[q,5]==1){
      if(((pM_metric2[quartets_matrix[q,1],quartets_matrix[q,2]]+pM_metric2[quartets_matrix[q,3],quartets_matrix[q,4]])<
          (pM_metric2[quartets_matrix[q,1],quartets_matrix[q,3]]+pM_metric2[quartets_matrix[q,2],quartets_matrix[q,4]]))|
         ((pM_metric2[quartets_matrix[q,1],quartets_matrix[q,2]]+pM_metric2[quartets_matrix[q,3],quartets_matrix[q,4]])<
          (pM_metric2[quartets_matrix[q,1],quartets_matrix[q,4]]+pM_metric2[quartets_matrix[q,2],quartets_matrix[q,3]]))){
        return(0)
      }
    }else if(quartets_matrix[q,5]==2){
      if(((pM_metric2[quartets_matrix[q,1],quartets_matrix[q,3]]+pM_metric2[quartets_matrix[q,2],quartets_matrix[q,4]])<
          (pM_metric2[quartets_matrix[q,1],quartets_matrix[q,2]]+pM_metric2[quartets_matrix[q,3],quartets_matrix[q,4]]))|
         ((pM_metric2[quartets_matrix[q,1],quartets_matrix[q,3]]+pM_metric2[quartets_matrix[q,2],quartets_matrix[q,4]])<
          (pM_metric2[quartets_matrix[q,1],quartets_matrix[q,4]]+pM_metric2[quartets_matrix[q,2],quartets_matrix[q,3]]))){
        return(0)
      }
    }else if(quartets_matrix[q,5]==3){
      if(((pM_metric2[quartets_matrix[q,1],quartets_matrix[q,4]]+pM_metric2[quartets_matrix[q,2],quartets_matrix[q,3]])<
          (pM_metric2[quartets_matrix[q,1],quartets_matrix[q,3]]+pM_metric2[quartets_matrix[q,2],quartets_matrix[q,4]]))|
         ((pM_metric2[quartets_matrix[q,1],quartets_matrix[q,4]]+pM_metric2[quartets_matrix[q,2],quartets_matrix[q,3]])<
          (pM_metric2[quartets_matrix[q,1],quartets_matrix[q,2]]+pM_metric2[quartets_matrix[q,3],quartets_matrix[q,4]]))){
        return(0)
      }
    }
  }
  return(1)
}

#find right quartets, if no result, output label_k=0
search_quartets<-function(isomorph_quartets_matrix_list_i_j,isomorph_quartets_list_i_j,pM_metric2){
  label_k<-0
  if(isomorph_quartets_list_i_j== -1){
    label_k<-1
  }else if(isomorph_quartets_list_i_j>0){
    if(length(isomorph_quartets_matrix_list_i_j)==isomorph_quartets_list_i_j){
      for (i in 1:isomorph_quartets_list_i_j) {
        if(check_one_quartets_for_metric(isomorph_quartets_matrix_list_i_j[[i]][,(1:5)],pM_metric2)){
          label_k<-i
          break
        }
      }
    }else if(ncol(isomorph_quartets_matrix_list_i_j)==(isomorph_quartets_list_i_j+4)){
      for (i in 1:isomorph_quartets_list_i_j) {
        if(check_one_quartets_for_metric(isomorph_quartets_matrix_list_i_j[,c((1:4),4+i)],pM_metric2)){
          label_k<-i
          break
        }
      }
    }
  }
  return(label_k)
}

#check metric and k23
check_one_k23_for_metric<-function(k23_matrix,pM_metric2,i){
  if(length(k23_matrix)==6){
    if(k23_matrix[6]==1){
      return(1)
    }else if(k23_matrix[6]==2){
      if(i==1){
        if((pM_metric2[k23_matrix[3],k23_matrix[5]]+pM_metric2[k23_matrix[4],k23_matrix[5]]+pM_metric2[k23_matrix[1],k23_matrix[2]])>
           (pM_metric2[k23_matrix[1],k23_matrix[5]]+pM_metric2[k23_matrix[2],k23_matrix[5]]+pM_metric2[k23_matrix[3],k23_matrix[4]])){
          return(0)
        }
      }else if(i==2){
        if((pM_metric2[k23_matrix[1],k23_matrix[5]]+pM_metric2[k23_matrix[2],k23_matrix[5]]+pM_metric2[k23_matrix[3],k23_matrix[4]])>
           (pM_metric2[k23_matrix[3],k23_matrix[5]]+pM_metric2[k23_matrix[4],k23_matrix[5]]+pM_metric2[k23_matrix[1],k23_matrix[2]])){
          return(0)
        }
      }
    }
  }else if(length(k23_matrix)>=10){
    k23_number<-which(k23_matrix[,6]==2)
    if(ncol(k23_matrix)==6){
      if(length(k23_number)==0){
        return(1)
      }else if(length(k23_number)==1){
        if(i==1){
          if((pM_metric2[k23_matrix[k23_number,3],k23_matrix[k23_number,5]]+
              pM_metric2[k23_matrix[k23_number,4],k23_matrix[k23_number,5]]+
              pM_metric2[k23_matrix[k23_number,1],k23_matrix[k23_number,2]])>
             (pM_metric2[k23_matrix[k23_number,1],k23_matrix[k23_number,5]]+
              pM_metric2[k23_matrix[k23_number,2],k23_matrix[k23_number,5]]+
              pM_metric2[k23_matrix[k23_number,3],k23_matrix[k23_number,4]])){
            return(0)
          }
        }else if(i==2){
          if((pM_metric2[k23_matrix[k23_number,1],k23_matrix[k23_number,5]]+
              pM_metric2[k23_matrix[k23_number,2],k23_matrix[k23_number,5]]+
              pM_metric2[k23_matrix[k23_number,3],k23_matrix[k23_number,4]])>
             (pM_metric2[k23_matrix[k23_number,3],k23_matrix[k23_number,5]]+
              pM_metric2[k23_matrix[k23_number,4],k23_matrix[k23_number,5]]+
              pM_metric2[k23_matrix[k23_number,1],k23_matrix[k23_number,2]])){
            return(0)
          }
        }
      }
    }else if(ncol(k23_matrix)==7){
      k23_matrix_sub<-k23_matrix[k23_number,]
      for (j in 1:length(k23_number)) {
        if(k23_matrix_sub[j,7]==1){
          if((pM_metric2[k23_matrix_sub[j,3],k23_matrix_sub[j,5]]+
              pM_metric2[k23_matrix_sub[j,4],k23_matrix_sub[j,5]]+
              pM_metric2[k23_matrix_sub[j,1],k23_matrix_sub[j,2]])>
             (pM_metric2[k23_matrix_sub[j,1],k23_matrix_sub[j,5]]+
              pM_metric2[k23_matrix_sub[j,2],k23_matrix_sub[j,5]]+
              pM_metric2[k23_matrix_sub[j,3],k23_matrix_sub[j,4]])){
            return(0)
          }
        }else if(k23_matrix_sub[j,7]==2){
          if((pM_metric2[k23_matrix_sub[j,1],k23_matrix_sub[j,5]]+
              pM_metric2[k23_matrix_sub[j,2],k23_matrix_sub[j,5]]+
              pM_metric2[k23_matrix_sub[j,3],k23_matrix_sub[j,4]])>
             (pM_metric2[k23_matrix_sub[j,3],k23_matrix_sub[j,5]]+
              pM_metric2[k23_matrix_sub[j,4],k23_matrix_sub[j,5]]+
              pM_metric2[k23_matrix_sub[j,1],k23_matrix_sub[j,2]])){
            return(0)
          }
        }
      }
    }
  }
  return(1)
}

#find right k23, if no result, output label_l=0
search_k23<-function(isomorph_k23_matrix_list_i_j_k,isomorph_k23_list_i_j_k,pM_metric2){
  label_l<-0
  if(isomorph_k23_list_i_j_k== -1){
    label_l<-1
  }else if(isomorph_k23_list_i_j_k>0){
    if(length(isomorph_k23_matrix_list_i_j_k)==isomorph_k23_list_i_j_k){
      for (i in 1:isomorph_k23_list_i_j_k) {
        k23_matrix<-isomorph_k23_matrix_list_i_j_k[[i]]
        if(length(k23_matrix)>6){
          if(ncol(k23_matrix)>6){
            k23_matrix<-k23_matrix[,(1:7)]
          }
        }
        if(check_one_k23_for_metric(k23_matrix,pM_metric2,i)){
          label_l<-i
          break
        }
      }
    }else if(ncol(isomorph_k23_matrix_list_i_j_k)!=isomorph_k23_list_i_j_k){
      for (i in 1:isomorph_k23_list_i_j_k) {
        k23_matrix<-isomorph_k23_matrix_list_i_j_k
        if(length(k23_matrix)>6){
          if(ncol(k23_matrix)>6){
            k23_matrix<-k23_matrix[,c((1:6),6+i)]
          }
        }
        if(check_one_k23_for_metric(k23_matrix,pM_metric2,i)){
          label_l<-i
          break
        }
      }
    }
  }
  return(label_l)
}

#get case number by c(label_i,label_j,label_k,label_l)
search_case_number<-function(label,isomorph_f_list,isomorph_k23_list){
  number<-0
  if(label[1]>1){
    for (i in 1:(label[1]-1)) {
      for (j in 1:isomorph_f_list[[i]]) {
        a<-sum(abs(isomorph_k23_list[[i]][[j]]))
        number<-number+a
      }
    }
  }
  if(label[2]>1){
    for (j in 1:(label[2]-1)) {
      a<-sum(abs(isomorph_k23_list[[label[1]]][[j]]))
      number<-number+a
    }
  }
  if(label[3]>1){
    number<-number+sum(abs(isomorph_k23_list[[label[1]]][[label[2]]][1:(label[3]-1)]))
  }
  number<-number+label[4]
  
  return(number)
}