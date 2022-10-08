make_adj<-function(data, label, level, node_select ,num_cluster=7 ){
  
  margi_val<-list()
  result_list<-list()
  
  for (r in 1:num_cluster){
    if (r==0){
      x <- data[,,,,]
      x_stack <- array(0,dim=dim(data)[1:4])
      for(s in 1:length(label)){
        x_stack <- x_stack + x[,,,,s]
      }
    }else{
      x <- data[,,,,which(label==r)]
      x_stack <- array(0,dim=dim(data)[1:4])
      for(s in 1:length(which(label==r))){
        x_stack <- x_stack + x[,,,,s]
      }
    }
    
    adjacency_matrix_total<-array(0,dim = c(36,36))
    
    for (i in 1:length(level)){
      for (j in 1:length(level)){
        
        adjacency_matrix_total[(sum(level[1:i])-level[i]+1):sum(level[1:i]),(sum(level[1:j])-level[j]+1):sum(level[1:j])] <- x_stack[1:level[i],1:level[j],i,j]
        
      }
    }
    adjacency_matrix_total2<-array(0,dim=dim(adjacency_matrix_total))
    
    adjacency_matrix_total2<-t(adjacency_matrix_total)
    adjacency_matrix_total2[upper.tri(adjacency_matrix_total2)] <- adjacency_matrix_total[upper.tri(adjacency_matrix_total)]
    
    adjacency_matrix_total2[is.na(adjacency_matrix_total2)] = 0
    
    
    number<-rep(1:length(level),level)
    node_name<-rep('c',sum(level))
    
    numbering_list<-c()
    numbering_list2<-c()
    for (i in 1:length(level)){
      numbering_list<-c(numbering_list,seq(1,level[i]))
      
      if (i>=10){
        for (j in 1:level[i]){
          numbering_list2<-c(numbering_list2,',')
        }
      }else{
        for (j in 1:level[i]){
          numbering_list2<-c(numbering_list2,'')
        }
      }
    }
    
    merge<-paste0(node_name,number)
    index_name<-paste0(merge,numbering_list2,numbering_list)
    
    
    rownames(adjacency_matrix_total2)<-index_name
    colnames(adjacency_matrix_total2)<-index_name
    
    select_name<-c()
    groups <-list()
    j=1
    for (i in node_select){
      
      select_name<-c(select_name,index_name[(sum(level[1:i])-level[i]+1):sum(level[1:i])])
      
      list_name <- paste0('x',i)
      
      groups[[list_name]]<-seq(j,j+level[i]-1)
      
      j=j+level[i]
    }
    
    
    
    index_cho<-c()
    for (i in 1:length(select_name)){
      index_cho<-append(index_cho,which(rownames(adjacency_matrix_total2)==select_name[i]))
    }
    
    if (r==0){
      adjacency_matrix_strength<-(adjacency_matrix_total2/length(label))[index_cho,index_cho]
    }else{
      adjacency_matrix_strength<-(adjacency_matrix_total2/sum(label==r))[index_cho,index_cho] 
    }
    
    margi_val<- append(margi_val,list(rowSums(adjacency_matrix_strength)/(length(node_select)-1)))
    
    result_adj<-adjacency_matrix_total2[index_cho,index_cho]
    
    list_name <- paste0("select_adj",r)
    
    result_list[[list_name]] <- as.matrix(result_adj)
  }
  return(list('result_list'=result_list,'margi_val'=margi_val, 'select_name'=select_name, 'groups'=groups))
}  
