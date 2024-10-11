################################################################################
## species weighting schemes
################################################################################

## Load data
spp_list <- read.csv(here(dirname(here()), "data", "maxent_spp_list.csv"))
spp_area_list <- read.csv(here(dirname(here()), "data", "spp_area_list.csv"))

spp_wght_fun <- function(df, scheme = 1, wv = c(1,2,4,6,8,2)){
  
  nm <- paste0("scheme_",scheme)
  
  # using {{}} := allows you to reference a column name as a string
  #https://stackoverflow.com/questions/73285003/how-do-i-pass-a-column-name-to-a-function-involving-mutate
  
  if(scheme == 1 | scheme == 3){
    
  df1 <- df %>% 
    
    mutate({{nm}} := case_when(redlistCategory == "Least Concern" ~ wv[1],
                              redlistCategory == "Near Threatened" ~ wv[2],
                              redlistCategory == "Vulnerable" ~ wv[3],
                              redlistCategory == "Endangered" ~ wv[4],
                              redlistCategory == "Critically Endangered" ~ wv[5],
                              redlistCategory == "Data Deficient" ~ wv[6]))
  
  
  }
  if(scheme == 2 | scheme == 4){
    
  df1 <- df %>% 
      
    mutate({{nm}} := case_when(redlistCategory == "Least Concern" ~ wv[1],
                                redlistCategory == "Near Threatened" ~ wv[2],
                                redlistCategory == "Vulnerable" ~ wv[3],
                                redlistCategory == "Endangered" ~ wv[4],
                                redlistCategory == "Critically Endangered" ~ wv[5],
                                redlistCategory == "Data Deficient" ~ wv[6])) %>%
    
      mutate({{nm}} := case_when(code == "8.1.1" | code == "8.1.2" ~ get(!!nm) + 1,
                            .default = get(!!nm)))
  
  }
  
  return(df1)
  
  }

spp_area_wght_fun <- function(df, scheme = 1, wv = c(1,2,4,6,8,2)){
  
  nm <- paste0("scheme_",scheme)
  
  areas <- c("Ecosystem","Community","Ramsar","Upstream")

  if(scheme == 1 | scheme == 3){
    
    df1 <- df %>% 
      
      mutate({{nm}} := case_when(threatStatus == "Least Concern" &
                                   !classGroup %in% areas ~ wv[1],
                                 threatStatus == "Near Threatened" &
                                   !classGroup %in% areas ~ wv[2],
                                 threatStatus == "Vulnerable" &
                                   !classGroup %in% areas ~ wv[3],
                                 threatStatus == "Endangered" &
                                   !classGroup %in% areas ~ wv[4],
                                 threatStatus == "Critically Endangered" &
                                   !classGroup %in% areas ~ wv[5],
                                 threatStatus == "Data Deficient" &
                                   !classGroup %in% areas ~ wv[6],
             .default = weight)) #%>%
      
      #mutate({{nm}} := case_when(classGroup %in% areas ~ get(!!weight)))
  }
  if(scheme == 2 |scheme == 4){
    
    df1 <- df %>% 
      
      mutate({{nm}} := case_when(threatStatus == "Least Concern" &
                                   !classGroup %in% areas ~ wv[1],
                                 threatStatus == "Near Threatened" &
                                   !classGroup %in% areas ~ wv[2],
                                 threatStatus == "Vulnerable" &
                                   !classGroup %in% areas ~ wv[3],
                                 threatStatus == "Endangered" &
                                   !classGroup %in% areas ~ wv[4],
                                 threatStatus == "Critically Endangered" &
                                   !classGroup %in% areas ~ wv[5],
                                 threatStatus == "Data Deficient" &
                                   !classGroup %in% areas ~ wv[6],
             .default = weight)) %>%
      
      #mutate({{nm}} := case_when(classGroup %in% areas ~ weight)) %>%
      
      mutate({{nm}} := case_when(code == "8.1.1" | code == "8.1.2" ~ get(!!nm) + 1,
                                 .default = get(!!nm)))
    
    
  }
  
  return(df1)
  
  }


## 1. Same threat ratings: no IAS
spp_list <- spp_wght_fun(spp_list)
spp_area_list <- spp_area_wght_fun(spp_area_list)

## 2. Same threat ratings: IAS = threat + 1
spp_list <- spp_wght_fun(spp_list, scheme = 2)
spp_area_list <- spp_area_wght_fun(spp_area_list, scheme = 2)

## 3. Different threat ratings (1,2,3,4,5,2): no IAS
spp_list <- spp_wght_fun(spp_list, scheme = 3, wv = c(1,2,3,4,5,2))
spp_area_list <- spp_area_wght_fun(spp_area_list, scheme = 3, wv = c(1,2,3,4,5,2))

## 4. Different threat ratings (1,2,3,4,5,2): IAS = threat + 1
spp_list <- spp_wght_fun(spp_list, scheme = 4, wv = c(1,2,3,4,5,2))
spp_area_list <- spp_area_wght_fun(spp_area_list, scheme = 4, wv = c(1,2,3,4,5,2))

## Write to disk
write.csv(spp_list, here(dirname(here()), "data", "maxent_spp_list_upd.csv"))
write.csv(spp_area_list, here(dirname(here()), "data", "spp_area_list_upd.csv"))
