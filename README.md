## *Modeling the viability of an endangered grassland passerine on a fragmented landscape*

#### Abby E. Bratt, Gary L. Slater, Scott F. Pearson, Ilai N. Keren, Sarah J. Converse

##### Please contact the first author for questions about the code or data: Abby Bratt (abby@proteus.co.nz)
##### Secondary contact: Sarah Converse (sconver@usgs.gov)

_______________________________________________________________________________________

## Abstract

The decline of grassland birds in North America represents a conservation crisis. 
Understanding the population dynamics of these declining species is critical for assessing 
population viability and identifying management actions likely to address population declines. 
In many grassland systems, habitat loss and fragmentation have resulted in the isolation 
of occupied sites, and trends can vary substantially between these sites. In these cases, 
it is essential to consider inter-site movement, which may reveal source-sink dynamics 
that can influence site-specific extirpation risks and range-wide population viability. 
We developed a multi-site integrated population model for the Streaked Horned Lark 
(*Eremophila alpestris strigata*; SHLA), a subspecies native to lowland prairies 
in the Cascadia bioregion of North America. We modeled the population in the rapidly 
urbanizing South Puget Lowlands of Washington State, USA, where SHLA is now restricted 
to just a handful of sites. Using nest monitoring, mark-resight, and count data collected 
from nine occupied sites since 2010, we explicitly modeled inter-site movements to assess 
population trends and predict viability. Our results indicate that the regional abundance 
of SHLA has increased by approximately 2% (95% credible interval 0-3%) annually since 2010. 
However, this rate varied across sites, with some experiencing substantial growth and others
showing little or no growth. Movement dynamics revealed nuanced site-specific trends that 
inform local extirpation risk. The predicted probability of regional extirpation from the 
South Puget Lowlands over the next 20 years was 0, though local extirpation was possible 
given uncertainty in parameter values and stochasticity across years and sites. 
Results emphasize the importance of considering inter-site movement and variability 
in demographic rates to predict viability and inform conservation of species on 
fragmented landscapes. Our findings will guide continued monitoring and management efforts, 
including reintroduction of SHLA to currently unoccupied sites, ultimately enhancing the 
long-term viability of SHLA in the South Puget Lowlands.

### Details of Article 

Bratt AE, GL Slater, SF Pearson, IL Keren, and SJ Converse. *In review*. 
Modeling the viability of an endangered grassland passerine on a fragmented landscape. 

### [Scripts](./scripts)

Contains scripts to run all analyses. 

### How to Use this Repository 

This repository contains all code required to fit the integrated population model 
to processed data, conduct population viability analysis, conduct posterior predictive checks, 
and process results to create figures and tables. The scripts may be run in order of their prefix. 
The data supporting this research are restricted and not available publicly; 
contact the lead author at the address above for more information. 

### Required Packages and Versions Used 

ggpubr_0.6.3       
gridExtra_2.3.1    
cowplot_1.2.0      
wesanderson_0.3.7  
ggdist_3.3.3      
tidybayes_3.0.7    
viridis_0.6.5      
viridisLite_0.4.3  
RColorBrewer_1.1-3 
strex_2.0.1       
postpack_0.5.4     
gtools_3.9.5       
lubridate_1.9.5    
forcats_1.0.1      
stringr_1.6.0     
dplyr_1.2.0        
purrr_1.2.1        
readr_2.2.0        
tidyr_1.3.2        
tibble_3.3.1      
tidyverse_2.0.0    
ggplot2_4.0.2      
coda_0.19-4.1      
beepr_2.0          
here_1.0.2        
nimble_1.4.1 