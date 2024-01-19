search_case_number_by_generic_metric<-function(M_metric,data_path = getwd()){
  n<-sqrt(length(M_metric))
  #calculate 7 points 1-solution
  ilp_1solution_result<-get_1solution_by_ilp(M_metric)
  cycle_list_7points<-get_sub_cycle_in_1solution(ilp_1solution_result,(1:n))
  #3 cases for 7 points 1-solution
  if(length(cycle_list_7points[[1]])==7){
    casenumberin37<-1
    permutation<-(1:n)
    pM_metric1<-M_metric
  }
  if(length(cycle_list_7points[[1]])==5){
    #calculate 6 points 1-solution
    label7<-get_label_for_v_in_1solution_6points__when_remove_u(M_metric,cycle_list_7points,remove_u = 6)
    label6<-get_label_for_v_in_1solution_6points__when_remove_u(M_metric,cycle_list_7points,remove_u = 7)
    label_vector<-c(label6,label7)
    
    permutation_result<-get_casenumberin37_and_permutation(label_vector,cycle_list_7points)
    casenumberin37<-permutation_result[[1]]
    permutation1<-permutation_result[[2]]
    pM_metric1<-permute_matrix(M_metric,permutation1)
  }
  if(length(cycle_list_7points[[1]])==3){
    #calculate 6 points 1-solution
    label5<-get_label_for_v_in_1solution_6points__when_remove_u(M_metric,cycle_list_7points,remove_u = 4)
    label4<-get_label_for_v_in_1solution_6points__when_remove_u(M_metric,cycle_list_7points,remove_u = 5)
    label7<-get_label_for_v_in_1solution_6points__when_remove_u(M_metric,cycle_list_7points,remove_u = 6)
    label6<-get_label_for_v_in_1solution_6points__when_remove_u(M_metric,cycle_list_7points,remove_u = 7)
    label_vector<-c(label4,label5,label6,label7)
    
    permutation_result<-get_casenumberin37_and_permutation(label_vector,cycle_list_7points)
    casenumberin37<-permutation_result[[1]]
    permutation1<-permutation_result[[2]]
    pM_metric1<-permute_matrix(M_metric,permutation1)
  }
  #read data by original case number, choose permutation2_matrix and fesible_taxa_list for Gromov structure and extend graph
  #load list data
  setwd(paste0(data_path,"/rdatalist/"))
  load(paste0("G_matrix_isomorph_list_7_",casenumberin37,".rda"))
  load(paste0("isomorph_f_list_7_",casenumberin37,".rda"))
  load(paste0("isomorph_f_matrix_list_7_",casenumberin37,".rda"))
  load(paste0("isomorph_k23_list_7_",casenumberin37,".rda"))
  load(paste0("isomorph_k23_matrix_list_7_",casenumberin37,".rda"))
  load(paste0("isomorph_quartets_list_7_",casenumberin37,".rda"))
  load(paste0("isomorph_quartets_matrix_list_7_",casenumberin37,".rda"))
  #use general variable names
  G_matrix_isomorph_list<-get(paste0("G_matrix_isomorph_list_7_",casenumberin37))
  isomorph_f_list<-get(paste0("isomorph_f_list_7_",casenumberin37))
  isomorph_f_matrix_list<-get(paste0("isomorph_f_matrix_list_7_",casenumberin37))
  isomorph_k23_list<-get(paste0("isomorph_k23_list_7_",casenumberin37))
  isomorph_k23_matrix_list<-get(paste0("isomorph_k23_matrix_list_7_",casenumberin37))
  isomorph_quartets_list<-get(paste0("isomorph_quartets_list_7_",casenumberin37))
  isomorph_quartets_matrix_list<-get(paste0("isomorph_quartets_matrix_list_7_",casenumberin37))
  #remove old data
  rm(list = c(paste0("G_matrix_isomorph_list_7_",casenumberin37),paste0("isomorph_f_list_7_",casenumberin37),
              paste0("isomorph_f_matrix_list_7_",casenumberin37),paste0("isomorph_k23_list_7_",casenumberin37),
              paste0("isomorph_k23_matrix_list_7_",casenumberin37),paste0("isomorph_quartets_list_7_",casenumberin37),
              paste0("isomorph_quartets_matrix_list_7_",casenumberin37)))
  #choose permutation2_matrix and fesible_taxa_list for Gromov structure and extend graph
  if(casenumberin37==1){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    2,3,4,5,6,7,1,
                                    3,4,5,6,7,1,2,
                                    4,5,6,7,1,2,3,
                                    5,6,7,1,2,3,4,
                                    6,7,1,2,3,4,5,
                                    7,1,2,3,4,5,6,
                                    7,6,5,4,3,2,1,
                                    6,5,4,3,2,1,7,
                                    5,4,3,2,1,7,6,
                                    4,3,2,1,7,6,5,
                                    3,2,1,7,6,5,4,
                                    2,1,7,6,5,4,3,
                                    1,7,6,5,4,3,2),nrow = 7))
    a1<-c(3,4,5,6)
    a2<-c(4,5,6,7)
    a3<-c(1,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3,7)
    a6<-c(1,2,3,4)
    a7<-c(2,3,4,5)
  }else if(casenumberin37==2){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    4,3,2,1,5,6,7,
                                    4,3,2,1,5,7,6),nrow = 7))
    a1<-c(3,4,6,7)
    a2<-c(4,5,6,7)
    a3<-c(1,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(2,3,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3,4)
  }else if(casenumberin37==3){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    3,2,1,5,4,7,6),nrow = 7))
    a1<-c(3,4,6,7)
    a2<-c(4,5,6,7)
    a3<-c(1,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(2,3,6,7)
    a6<-c(1,2,3,5)
    a7<-c(1,2,3,4)
  }else if(casenumberin37==4){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    2,1,5,4,3,7,6),nrow = 7))
    a1<-c(3,4,6,7)
    a2<-c(4,5,6,7)
    a3<-c(1,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(2,3,6,7)
    a6<-c(1,2,4,5)
    a7<-c(1,2,3,4)
  }else if(casenumberin37==5){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    1,2,3,5,4,6,7,
                                    1,2,3,5,4,7,6,
                                    1,2,3,6,7,5,4,
                                    1,2,3,6,7,4,5,
                                    1,2,3,7,6,5,4,
                                    1,2,3,7,6,4,5,
                                    1,3,2,4,5,6,7,
                                    1,3,2,4,5,7,6,
                                    1,3,2,5,4,6,7,
                                    1,3,2,5,4,7,6,
                                    1,3,2,6,7,5,4,
                                    1,3,2,6,7,4,5,
                                    1,3,2,7,6,5,4,
                                    1,3,2,7,6,4,5,
                                    2,1,3,4,5,6,7,
                                    2,1,3,4,5,7,6,
                                    2,1,3,5,4,6,7,
                                    2,1,3,5,4,7,6,
                                    2,1,3,6,7,5,4,
                                    2,1,3,6,7,4,5,
                                    2,1,3,7,6,5,4,
                                    2,1,3,7,6,4,5,
                                    2,3,1,4,5,6,7,
                                    2,3,1,4,5,7,6,
                                    2,3,1,5,4,6,7,
                                    2,3,1,5,4,7,6,
                                    2,3,1,6,7,5,4,
                                    2,3,1,6,7,4,5,
                                    2,3,1,7,6,5,4,
                                    2,3,1,7,6,4,5,
                                    3,1,2,4,5,6,7,
                                    3,1,2,4,5,7,6,
                                    3,1,2,5,4,6,7,
                                    3,1,2,5,4,7,6,
                                    3,1,2,6,7,5,4,
                                    3,1,2,6,7,4,5,
                                    3,1,2,7,6,5,4,
                                    3,1,2,7,6,4,5,
                                    3,2,1,4,5,6,7,
                                    3,2,1,4,5,7,6,
                                    3,2,1,5,4,6,7,
                                    3,2,1,5,4,7,6,
                                    3,2,1,6,7,5,4,
                                    3,2,1,6,7,4,5,
                                    3,2,1,7,6,5,4,
                                    3,2,1,7,6,4,5),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,3)
    a5<-c(1,2,3)
    a6<-c(1,2,3)
    a7<-c(1,2,3)
  }else if(casenumberin37==6){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    2,1,3,4,5,6,7,
                                    2,1,3,4,5,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3)
    a6<-c(1,2,3)
    a7<-c(1,2,3)
  }else if(casenumberin37==7){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    1,2,3,5,4,6,7,
                                    1,2,3,5,4,7,6,
                                    2,1,3,4,5,6,7,
                                    2,1,3,4,5,7,6,
                                    2,1,3,5,4,6,7,
                                    2,1,3,5,4,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,3)
    a7<-c(1,2,3)
  }else if(casenumberin37==8){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    1,3,2,5,4,6,7,
                                    1,3,2,5,4,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,3,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,3)
    a7<-c(1,2,3)
  }else if(casenumberin37==9){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,6,7,4,5,
                                    2,1,3,4,5,6,7,
                                    2,1,3,6,7,4,5),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3)
    a6<-c(1,2,4,5)
    a7<-c(1,2,3)
  }else if(casenumberin37==10){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,3,2,6,7,4,5),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,3,6,7)
    a5<-c(1,2,3)
    a6<-c(1,2,4,5)
    a7<-c(1,2,3)
  }else if(casenumberin37==11){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,5,4,6,7,
                                    2,1,3,4,5,6,7,
                                    2,1,3,5,4,6,7),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,4,5)
    a7<-c(1,2,3)
  }else if(casenumberin37==12){
    permutation2_matrix<-0
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,3,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,4,5)
    a7<-c(1,2,3)
  }else if(casenumberin37==13){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,5,4,6,7),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,3,4,5)
    a7<-c(1,2,3)
  }else if(casenumberin37==14){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    2,1,3,5,4,6,7),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(2,3,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,4,5)
    a7<-c(1,2,3)
  }else if(casenumberin37==15){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    1,2,3,5,4,6,7,
                                    1,2,3,5,4,7,6,
                                    1,2,3,6,7,4,5,
                                    1,2,3,6,7,5,4,
                                    1,2,3,7,6,4,5,
                                    1,2,3,7,6,5,4,
                                    2,1,3,4,5,6,7,
                                    2,1,3,4,5,7,6,
                                    2,1,3,5,4,6,7,
                                    2,1,3,5,4,7,6,
                                    2,1,3,6,7,4,5,
                                    2,1,3,6,7,5,4,
                                    2,1,3,7,6,4,5,
                                    2,1,3,7,6,5,4),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,4,5)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==16){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,3,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,4,5)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==17){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    1,2,3,5,4,6,7,
                                    1,2,3,5,4,7,6,
                                    1,3,2,6,7,4,5,
                                    1,3,2,6,7,5,4,
                                    1,3,2,7,6,4,5,
                                    1,3,2,7,6,5,4),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,3,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,4,5)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==18){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,6,7,4,5,
                                    1,3,2,5,4,7,6,
                                    1,3,2,7,6,5,4),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,3,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,3,4,5)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==19){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    2,1,3,5,4,6,7,
                                    2,1,3,5,4,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(2,3,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,4,5)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==20){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    2,1,3,6,7,4,5),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(2,3,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,3,4,5)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==21){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    2,1,3,4,5,6,7),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==22){
    permutation2_matrix<-0
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,3,4,5)
  }else if(casenumberin37==23){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    2,1,3,4,5,6,7),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3)
  }else if(casenumberin37==24){
    permutation2_matrix<-0
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==25){
    permutation2_matrix<-0
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,3,4,5)
  }else if(casenumberin37==26){
    permutation2_matrix<-0
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,3,4)
    a7<-c(2,3,4,5)
  }else if(casenumberin37==27){
    permutation2_matrix<-0
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3)
  }else if(casenumberin37==28){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    2,1,3,4,5,6,7),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3)
    a6<-c(1,2,3,4)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==29){
    permutation2_matrix<-0
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3)
    a6<-c(1,2,3,4)
    a7<-c(1,3,4,5)
  }else if(casenumberin37==30){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    2,1,3,4,5,6,7),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3)
  }else if(casenumberin37==31){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    2,1,3,4,5,6,7,
                                    2,1,3,4,5,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3,4)
  }else if(casenumberin37==32){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3,4)
  }else if(casenumberin37==33){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,4,5,7,6,
                                    2,1,3,4,5,6,7,
                                    2,1,3,4,5,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3,4)
  }else if(casenumberin37==34){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,5,4,7,6,
                                    2,1,3,4,5,6,7,
                                    2,1,3,5,4,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3,5)
  }else if(casenumberin37==35){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,3,2,5,4,7,6),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,3,6,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,3,5)
  }else if(casenumberin37==36){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,2,3,7,6,5,4,
                                    2,1,3,4,5,6,7,
                                    2,1,3,7,6,5,4),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3,7)
    a6<-c(1,2,3,4)
    a7<-c(1,2,4,5)
  }else if(casenumberin37==37){
    permutation2_matrix<-t(matrix(c(1,2,3,4,5,6,7,
                                    1,3,2,7,6,5,4),nrow = 7))
    a1<-c(4,5,6,7)
    a2<-c(4,5,6,7)
    a3<-c(4,5,6,7)
    a4<-c(1,2,6,7)
    a5<-c(1,2,3,7)
    a6<-c(1,2,3,4)
    a7<-c(1,3,4,5)
  }
  fesible_taxa_list<-list(a1,a2,a3,a4,a5,a6,a7)
  #traverse permutation2, find label_i,j,k,l
  if(length(permutation2_matrix)==1){
    permutation2<-(1:7)
    pM_metric2<-pM_metric1
    label_i<-search_gromov_structure(G_matrix_isomorph_list,pM_metric2,fesible_taxa_list)
    label_j<-search_extended_graph(isomorph_f_matrix_list[[label_i]],pM_metric2,fesible_taxa_list)
    label_k<-search_quartets(isomorph_quartets_matrix_list[[label_i]][[label_j]],
                             isomorph_quartets_list[[label_i]][[label_j]],pM_metric2)
    label_l<-search_k23(isomorph_k23_matrix_list[[label_i]][[label_j]][[label_k]],
                        isomorph_k23_list[[label_i]][[label_j]][[label_k]],pM_metric2)
  }else{
    for (permutation2_i in 1:nrow(permutation2_matrix)) {
      permutation2<-permutation2_matrix[permutation2_i,]
      pM_metric2<-permute_matrix(pM_metric1,permutation2)
      #find right Gromov structure, record label_j
      label_i<-search_gromov_structure(G_matrix_isomorph_list,pM_metric2,fesible_taxa_list)
      if(label_i>0){
        #find right extended graph, record label_j
        label_j<-search_extended_graph(isomorph_f_matrix_list[[label_i]],pM_metric2,fesible_taxa_list)
        if(label_j>0){
          #find right quartets, record label_k
          label_k<-search_quartets(isomorph_quartets_matrix_list[[label_i]][[label_j]],
                                   isomorph_quartets_list[[label_i]][[label_j]],pM_metric2)
          if(label_k>0){
            #find right k23，record label_l
            label_l<-search_k23(isomorph_k23_matrix_list[[label_i]][[label_j]][[label_k]],
                                isomorph_k23_list[[label_i]][[label_j]][[label_k]],pM_metric2)
            if(label_l>0){
              break
            }
          }
        }
      }
    }
  }
  if((label_i*label_j*label_k*label_l)==0){
    print("label_i_j_k_l error")
  }
  case_number<-search_case_number(c(label_i,label_j,label_k,label_l),isomorph_f_list,isomorph_k23_list)
  permutation_combine<-permutation1[permutation2]
  #adjust some outputs
  #Gromov_structure_and_extended_graph<-isomorph_f_matrix_list[[label_i]][[label_j]]
  #adjust quartets output
  if(isomorph_quartets_list[[label_i]][[label_j]]== -1){
    Quartets<-0
  }else if(length(isomorph_quartets_matrix_list[[label_i]][[label_j]])==isomorph_quartets_list[[label_i]][[label_j]]){
    Quartets<-isomorph_quartets_matrix_list[[label_i]][[label_j]][[label_k]][,(1:5)]
  }else if(length(isomorph_quartets_matrix_list[[label_i]][[label_j]])!=isomorph_quartets_list[[label_i]][[label_j]]){
    Quartets<-isomorph_quartets_matrix_list[[label_i]][[label_j]][,c((1:4),4+label_k)]
  }
  rownames(Quartets)<-NULL
  colnames(Quartets)<-NULL
  if(length(Quartets)==10){
    if(!(0 %in% (Quartets[1,]==Quartets[2,]))){
      Quartets<-Quartets[1,]
    }
  }
  #adjust k23 output
  if(isomorph_k23_list[[label_i]][[label_j]][[label_k]]== -1){
    K23<-0
  }else if(length(isomorph_k23_matrix_list[[label_i]][[label_j]][[label_k]])==isomorph_k23_list[[label_i]][[label_j]][[label_k]]){
    K23<-isomorph_k23_matrix_list[[label_i]][[label_j]][[label_k]][[label_l]]
    if(length(K23)>6){
      if(ncol(K23)>6){
        K23<-K23[,(1:7)]
      }
    }
  }else if(length(isomorph_k23_matrix_list[[label_i]][[label_j]][[label_k]])!=isomorph_k23_list[[label_i]][[label_j]][[label_k]]){
    K23<-isomorph_k23_matrix_list[[label_i]][[label_j]][[label_k]]
    if(length(K23)>6){
      if(ncol(K23)>6){
        K23<-K23[,c((1:6),6+label_l)]
      }
    }
  }
  rownames(K23)<-NULL
  colnames(K23)<-NULL
  if(length(K23)==6){
    if(K23[6]==2){
      K23<-c(K23,label_l)
    }
  }else if(length(K23)>10){
    if(ncol(K23)==6){
      k23_2_sides_number<-which(K23[,6]==2)
      if(length(k23_2_sides_number)==1){
        K23<-cbind(K23,0)
        K23[k23_2_sides_number,7]<-label_l
      }
    }
  }
  #output: case_number, permutation_combine, base_case_number, f_matrix, quartets, k23
  result_list<-list(case_number = case_number,
                    permutation = permutation_combine,
                    Base_case = casenumberin37,
                    Gromov_structure_and_extended_graph = isomorph_f_matrix_list[[label_i]][[label_j]],
                    Quartets = Quartets,
                    K23 = K23)
  return(result_list)
}
