generate_plot<-function(result_list,groups, save=F, save_path = "", file_type="pdf"){
  
  if (save==T){
    setwd(save_path)
    file_type <<-file_type
  }
  
  par(mfrow=c(1,1))
  indegree_list<<-list()
  indegree_max_list<<-list()
  
  groups <<- groups
  layout_plot <<- 'circle'
  edge_color <<- 'black'
  
  for (i in 1:length(result_list)){
    code_line<- paste0("graph_",i,"<- qgraph(result_list[[",i,"]],plot=TRUE)")
    eval(parse(text=code_line), envir=.GlobalEnv)
    
    code_line2<- paste0("Centrality_",i,"<- centrality(graph_",i,")")
    eval(parse(text=code_line2), envir=.GlobalEnv)    
    
    code_line3<- paste0("indegree_list[[",i,"]]<-","Centrality_",i,"$InDegree")
    eval(parse(text=code_line3), envir=.GlobalEnv)  
    
    
    code_line4<- paste0("indegree_max_list[[",i,"]]<-","max(graph_",i,"$Edgelist$weight)")
    eval(parse(text=code_line4), envir=.GlobalEnv)  
    
  }
  
  max_indegree <<- c()
  for (i in 1:length(indegree_list)){
    max_indegree <<- append(max_indegree,max(indegree_list[[i]]))
  }
  
  for (i in 1:length(result_list)){
    title_name<<-paste("Cluster" ,as.roman(i))
    if (save==T){
      code_line5<- paste0("finalgraph_",i,"<- qgraph(result_list$select_adj",i,", title = title_name,label.prop=1,groups=groups,label.scale.equal=TRUE,layout = layout_plot,vsize = margi_val2[[",i,"]]*exp(-min(max_indegree)/max_indegree[",i,"])*10 + 5
           ,edge.color = edge_color,labels = colnames(result_list$select_adj",i,"),esize=indegree_max_list[[",i,"]]/indegree_max_list[[which.max(indegree_max_list)]]*30,legend=TRUE,filetype=file_type, filename=paste0('graph',",i,"))")
      eval(parse(text=code_line5), envir=.GlobalEnv)
    }else{
      code_line5<- paste0("finalgraph_",i,"<- qgraph(result_list$select_adj",i,", title = title_name,label.prop=1,groups=groups,label.scale.equal=TRUE,layout = layout_plot,vsize = margi_val2[[",i,"]]*exp(-min(max_indegree)/max_indegree[",i,"])*10 + 5
           ,edge.color = edge_color,labels = colnames(result_list$select_adj",i,"),esize=indegree_max_list[[",i,"]]/indegree_max_list[[which.max(indegree_max_list)]]*30,legend=TRUE)")
      eval(parse(text=code_line5), envir=.GlobalEnv)
    }
    
  }  
}
