rm(list=ls())
#Testing the Git

if(T){
  install  = T
  import   = T
  summary  = T
  plotting = F
  
  # -----------------------------------------------------------------------------
  if(install==T){
    # install required packages (if not already installed)
    list.of.packages <- c("ggplot2", "dplyr", "reshape2", "lubridate")
    new.packages     <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    
    # loading required packages
    lapply(list.of.packages, require, character.only = TRUE)
  }
  # -----------------------------------------------------------------------------
  if(import==T){
    # load data
    RLS     <- read.csv("RLS_MPA_combined_long.csv")
    classes <- read.csv("classes.csv")
    growth  <- read.csv("SeaLifeBaseGrowthData.csv")
    
    # remove useless rows
    RLS <- RLS %>% filter(RecordedSpecies != "Survey not done")
    
    # add the vert/invert column to RLS dataset
    RLS$group <- classes$Group[match(RLS$CLASS, classes$Class)]
    
    
    size.classes <- sort(unique(RLS$SizeClass))
    RLS.SC   <- c(0,2.5,5,7.5,10,12.5,15,20,25,30,35,40,50,62.5,75, 400)
    log2.SC  <- 2^c(-5:10)
    log10.SC <- 10^c(-5:10)
    
    RLS$Date <- as.Date(RLS$SurveyDate, format="%d/%m/%Y")
    
    # creating a new col that has only the major groups in it
    RLS$SizeClassNew <- as.numeric(sapply(RLS$SizeClass, function(x) RLS.SC[which(RLS.SC / x > 1)[1]-1]))
    RLS$Log2Class    <- as.numeric(sapply(RLS$SizeClass, function(x) log2.SC[which(log2.SC / x > 1)[1]-1]))
    RLS$Log10Class   <- as.numeric(sapply(RLS$SizeClass, function(x) log10.SC[which(log10.SC / x > 1)[1]-1]))
    
    # creating the sub-groups
    inverts <- RLS %>% filter(group=="Invertebrate")
    verts   <- RLS %>% filter(group=="Vertebrate")
    
  }
  # -----------------------------------------------------------------------------
  if(summary==T){
    print(head(RLS))
    print(noquote(rep("=", 15, sep="")))
    n <- data.frame(name=c("nrow(RLS)","nrow(inverts)","nrow(verts)", "n.sites", "n.divers", 
                           "n.species", "n.vert.species", "n.invert.species", "n.vert.counts",
                           "n.invert.counts", "Earliest record", "Lastest record", "Mean surveys/site"), 
                    n = c(nrow(RLS),nrow(inverts),nrow(verts), length(unique(RLS$SiteCode)),
                          length(unique(RLS$Diver)), length(unique(RLS$RecordedSpecies)),
                          length(unique(verts$RecordedSpecies)), length(unique(inverts$RecordedSpecies)),
                          sum(verts$Abundance), sum(inverts$Abundance), 
                          as.character(format.Date(min(RLS$Date), format="%d-%b-%Y")), 
                          as.character(format.Date(max(RLS$Date), format="%d-%b-%Y")),
                          round(length(unique(RLS$SurveyID))/length(unique(RLS$SiteCode)),1) ))
    print(n)
  }
  # -----------------------------------------------------------------------------
  if(plotting==T){
    # Asking user for plot type wanted
    choose <- c("Abundance for specific site",
                "Total abundance",
                "Abundance Size spectrum",
                "Abundance Size spectrum (Invertebrates)",
                "Abundance Size spectrum (Vertebrates)")
    pick <- menu(choose, graphics=TRUE, title="Choose plot")
    # -----------------------------------------------------------------------------
    # Abundance for a specific site 
    if(pick==1){
      # asking user for site wanted
      site.n <- menu(levels(RLS$Site.name), graphics=TRUE, title="Choose site")
      
      # plotting the abundance dist for the specified site
      site <- RLS %>% filter(Site.name==levels(RLS$Site.name)[site.n])
      site.group <- aggregate(list(abundance=site$Abundance), 
                              by=list(size.class=site$SizeClass), FUN=sum)
      site.df <- data.frame(size.class=size.classes, 
                            abundance=site.group$abundance[match(size.classes, site.group$size.class)])
      windows()
      barplot(height=site.df$abundance, 
              main=paste("Site ",site.n,":", unique(site$Site.name)),
              xlab= "Size class", 
              ylab= "Abundance",
              cex.lab=1.5, 
              names.arg = site.df$size.class)
    }
    # -----------------------------------------------------------------------------
    # Total abundance
    if(pick == 2){
      RLS.tab <- aggregate(list(abundance=RLS$Abundance), 
                           by=list(size.class=RLS$SizeClass), FUN=sum)
      windows()
      barplot(height=RLS.tab$abundance,
              main="Total Abundance distribution",
              xlab= "Size class", 
              ylab= "Abundance",
              cex.lab=1.5, 
              names.arg = RLS.tab$size.class)
    }
    # -----------------------------------------------------------------------------
    # Abundance size spectrum 
    if(pick == 3){
      RLS.tab <- aggregate(list(abundance=RLS$Abundance), 
                           by=list(size.class=RLS$Log2Class), FUN=sum)
      windows()
      barplot(height=log2(RLS.tab$abundance),
              main="Total abundance size spectrum",
              xlab= "Size class (Log2 scale)", 
              ylab= "Log2(Abundance)",
              cex.lab=1.5, 
              names.arg = RLS.tab$size.class)
    }
    # -----------------------------------------------------------------------------
    # Abundance size spectrum (invertebrates)
    if(pick == 4){
      inverts.tab <- aggregate(list(abundance=inverts$Abundance), 
                               by=list(size.class=inverts$Log2Class), FUN=sum)
      windows()
      barplot(height=log2(inverts.tab$abundance),
              main="Invertebrate abundance size spectrum",
              xlab= "Size class (Log2 scale)", 
              ylab= "Log(Abundance)",
              cex.lab=1.5, 
              names.arg = inverts.tab$size.class)
    }
    # -----------------------------------------------------------------------------
    # Abundance size spectrum (vertebrates)
    if(pick == 5){
      verts.tab <- aggregate(list(abundance=verts$Abundance), 
                             by=list(size.class=verts$Log2Class), FUN=sum)
      windows()
      barplot(height=log2(verts.tab$abundance),
              main="Vertebrate abundance size spectrum",
              xlab= "Size class (Log2 scale)", 
              ylab= "Log(Abundance)",
              cex.lab=1.5, 
              names.arg = verts.tab$size.class)
    }
    # -----------------------------------------------------------------------------
    # Normalised total abundance size spectrum
    if(pick == 6){
      RLS.tab <- aggregate(list(abundance=RLS$Abundance), 
                           by=list(size.class=RLS$Log2Class), FUN=sum)
      windows()
      barplot(height=(log2(RLS.tab$abundance)/RLS.tab$size.class),
              main="Total abundance size spectrum",
              xlab= "Size class (Log2 scale)", 
              ylab= "Log2(Abundance)",
              cex.lab=1.5, 
              names.arg = RLS.tab$size.class)
    }
  }
}
dev.off()