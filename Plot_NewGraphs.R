library(tidyverse)

#load("C:/Users/Rob/Documents/!Post-doc1-Harvard!/Papers and Manuscripts/AutoBend/Example data.Rdata")

ROM.dplyrProcess<-function(ROM, Joint, grpvars,dirvars,FctOrder){
  # ROM The list of dataframes from Maya
  # Joint The name of the column containing the joint names
  # grpvars Names of the columns used to group the data for calculating averages
  # dirvars Names of the directions of ROM we want to keep
  # FctOrder A vector conssiting of the levels of the "Trial" factor in the order we want them plotted

  
  plotdat = ROM %>% 
    # Bind list of dataframes into one big dataframe
    bind_rows(.id = "df_label") %>% 
    # Group dataframe to calculate ROM mean and sd
    group_by_at((c(grpvars,Joint))) %>%
    # Select the variables (ROM directions) you want to plot
    select(all_of(dirvars)) %>% 
    # Calculate means and sds of ROM for each joint across specimens/models
    dplyr::summarise(across(where(is.numeric), list(mean = mean, sd = sd),na.rm=T)) %>% 
    # Regroup, dropping the joint variable for calculating whole-vertebral column
    # means and sds
    group_by_at(grpvars) %>% select(!(!!sym(Joint))) %>% 
    # Calculate whole-vertebral column means and sds
    dplyr::summarise(across(where(is.numeric),mean)) %>% 
    # Rearrange the dataframe to be in a more plot-friendly format
    pivot_longer(cols=where(is.numeric), names_to="Name",values_to="Value") %>% 
    separate(Name,into=c("direction","Measure")) %>% 
    pivot_wider(names_from=Measure,values_from=Value) %>% 
    # Change some column names to be more consistent with other functions
    dplyr::rename(ROM=mean,species=Species) %>% 
    # Relevel the "trial" factor so things plot in the right order
    mutate(Trial=fct_relevel(Trial,FctOrder)) %>% 
    arrange(Trial,species)
  
  plotdat
}

# ROM.ggPoints <- function(pltdf,col.scale,xlab,ylab,ylim,title){
#   # pltdf The dataframe of data to plot
#   # col.scale The color scale for our plots
#   # xlab Label for the x axis
#   # ylab Label for the y axis a vector conssiting of the levels of the "Trial" factor in the order we want them plotted
#   # ylim A vector of limits for the y-axis e.g. c(0,30)
#   # title Title of the graph
# 
#   # Plot with points
#   p1 = ggplot(pltdf, aes(x=species, y=ROM, color=direction,shape=Trial,
#                            group = interaction(Trial,direction))) + 
#     geom_point(position=position_dodge(.67),size=5) +
#     scale_color_manual("legend", values = col.scale) + 
#     xlab(xlab)+ylab(ylab)+ylim(ylim)+ggtitle(title)+
#     theme_classic()+
#     theme(axis.text.x=element_text(size=16,face="bold"),
#           axis.text.y=element_text(size=12),
#           axis.title.y=element_text(size=16,face="bold"))
# 
#   # Plot with points and error bars
#   p2 = p1 + geom_errorbar(aes(ymin=ROM-sd, ymax=ROM+sd), width=.25,size=1,
#                           position=position_dodge(.67)) 
#   p2
#   
# 
# }

ROM.ggPoints <- function(pltdf,col.scale,xlab,ylab,ylim,title){
  # pltdf The dataframe of data to plot
  # col.scale The color scale for our plots
  # xlab Label for the x axis
  # ylab Label for the y axis a vector conssiting of the levels of the "Trial" factor in the order we want them plotted
  # ylim A vector of limits for the y-axis e.g. c(0,30)
  # title Title of the graph
  # Plot with points
  p1 = ggplot(pltdf, aes(x=direction, y=ROM, color=direction,shape=Trial,
                         group = Trial)) + 
    geom_line(aes(group=Trial,linetype=Trial),color="black",position=position_dodge(.67))+ 
    geom_point(position=position_dodge(.67),size=5) +
    scale_color_manual("legend", values = col.scale) + 
    xlab(xlab)+ylab(ylab)+ylim(ylim)+ggtitle(title)+
    theme_classic()+
    theme(axis.text.x=element_text(size=16,face="bold"),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=16,face="bold"))+
    facet_wrap(~species)+
    theme(
      strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
      strip.text.x = element_text(size = 16, color = "black", face = "bold")
    )
  # Plot with points and error bars
  p2 = p1 + geom_errorbar(aes(ymin=ROM-sd, ymax=ROM+sd), width=.25,size=1,
                          position=position_dodge(.67)) 
  p2 
}

ROM.ggPointsLines <- function(pltdf,BaselineVars,col.scale,xlab,ylab,ylim,title){
  # pltdf The dataframe of data to plot
  # BaselineVars The column name and variable name e.g. c("Trial","real") we want to compare our other data to
  # col.scale The color scale for our plots
  # xlab Label for the x axis
  # ylab Label for the y axis a vector conssiting of the levels of the "Trial" factor in the order we want them plotted
  # ylim A vector of limits for the y-axis e.g. c(0,30)
  # title Title of the graph
  
  
  
  #This finds the "baseline" (real experimental data or our chosen model)
  # we want to compare each row (species,trial,direction) combination to  
  hline_dat = pltdf %>% filter(.,!!(sym(BaselineVars[1]))==BaselineVars[2])
  
  # Plot with points
  p1 = ggplot(pltdf, aes(x=species, y=ROM, color=direction,shape=Trial,
                         group = interaction(Trial,direction))) + 
    geom_point(position=position_dodge(.67),size=5) +
    scale_color_manual("legend", values = col.scale) + 
    xlab(xlab)+ylab(ylab)+ylim(ylim)+ggtitle(title)+
    theme_classic()+
    theme(axis.text.x=element_text(size=16,face="bold"),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=16,face="bold"))
  
  # Plot with points and error bars
  p2 = p1 + geom_errorbar(aes(ymin=ROM-sd, ymax=ROM+sd), width=.25,size=1,
                          position=position_dodge(.67)) 

  # Plot with points, error bars and horizontal "baselines" showing deviation from 
  # real data or our chosen the "best fit" model
  #p3 = p2 +  geom_errorbar(data=hline_dat,aes(x=species,y=NULL,ymin=ROM,ymax=ROM), width=0.67,size=1,
                           #position=position_dodge(.67),linetype="solid", color="black") 
  p3 = p2 + geom_line(linetype="dashed",position=position_dodge(.67), aes(x=species, y=ROM, color="black",
                                                                          group = interaction(direction, Trial)))
  
  p3
}

ROM.ggPointsFilledLines <- function(pltdf,col.scale,xlab,ylab,ylim,title){
  # pltdf The dataframe of data to plot
  # col.scale The color scale for our plots
  # xlab Label for the x axis
  # ylab Label for the y axis a vector conssiting of the levels of the "Trial" factor in the order we want them plotted
  # ylim A vector of limits for the y-axis e.g. c(0,30)
  # title Title of the graph
  
  # Plot with points
  p1 = ggplot(pltdf, aes(x=direction, y=ROM, color=direction,fill=direction,shape=Trial,
                         group = Trial)) + 
    geom_errorbar(aes(ymin=ROM-sd, ymax=ROM+sd), width=.25,size=1,
                  position=position_dodge(0.33)) +
    geom_line(aes(group=Trial,linetype=Trial),color="black",position=position_dodge(0.33), show.legend = FALSE)+ 
    geom_point(position=position_dodge(0.33),color="black",size=5) +
    scale_color_manual(values = col.scale, guide=FALSE) + 
    scale_fill_manual(values = col.scale, guide=guide_legend(order=1, title="Direction", 
                                                             override.aes = list(shape=21, fill=col.scale, color="Black"))) + 
    scale_shape_manual(values=c(21:25), guide=guide_legend(order=2, override.aes = list(fill="black")))+
    xlab(xlab)+ylab(ylab)+ylim(ylim)+ggtitle(title)+
    theme_classic()+
    theme(axis.text.x=element_text(size=16,face="bold"),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=16,face="bold"))+
    facet_wrap(~species)+
    theme(
      strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
      strip.text.x = element_text(size = 16, color = "black", face = "bold")
    )
  
  # Plot with points and error bars
  p1 
  
  
}

# 
# plotdat = ROM.dplyrProcess(ROM=comp1,Joint="JOINT",
#                            grpvars = c("df_label","Trial", "Species"),
#                            dirvars = c("Latero", "Sag", "Tor"),
#                            FctOrder = c("real","Error"))
# 
# ROM.ggPoints(plotdat,
#              col.scale=c("tomato","slategrey","gold"),
#                 title="Test Plot",xlab="",ylab="Range of Motion",ylim=c(0,35))
# 
# ROM.ggPointsLines(plotdat,BaselineVars = c("Trial","real"),
#                  col.scale=c("tomato","slategrey","gold"),
#                 title="Test Plot",xlab="",ylab="Range of Motion",ylim=c(0,35))
